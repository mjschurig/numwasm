/**
 * NumJS-WASM Broadcasting Implementation
 *
 * Adapted from NumPy's broadcasting implementation in:
 * numpy/_core/src/multiarray/nditer_constr.c lines 1480-1681
 *
 * The key insight is that broadcasting works via the "stride trick":
 * - When a dimension has size 1 but needs to broadcast to size N,
 *   we set the stride for that dimension to 0
 * - This means the pointer doesn't advance, effectively repeating the value
 */

#include "broadcast.h"
#include <stdlib.h>
#include <string.h>

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#define EXPORT EMSCRIPTEN_KEEPALIVE
#else
#define EXPORT
#endif

/**
 * Compute broadcast shape for two shapes.
 *
 * Algorithm adapted from NumPy nditer_constr.c lines 1511-1521:
 *
 * For each dimension (working from right to left):
 *   - If broadcast_shape[dim] == 1, set it to operand shape
 *   - Else if shapes don't match and operand shape != 1, error
 */
EXPORT int broadcast_shapes(const int32_t* shape1, int32_t ndim1,
                             const int32_t* shape2, int32_t ndim2,
                             int32_t* out_shape, int32_t* out_ndim)
{
    if (!out_shape || !out_ndim) return -1;

    int32_t result_ndim = (ndim1 > ndim2) ? ndim1 : ndim2;

    /* Initialize broadcast shape with 1s */
    for (int i = 0; i < result_ndim; i++) {
        out_shape[i] = 1;
    }

    /* Process first shape (align to right) */
    for (int i = 0; i < ndim1; i++) {
        int result_idx = i + result_ndim - ndim1;
        int32_t dim = shape1[i];

        if (out_shape[result_idx] == 1) {
            out_shape[result_idx] = dim;
        }
        else if (dim != 1 && dim != out_shape[result_idx]) {
            return -1; /* Incompatible shapes */
        }
    }

    /* Process second shape (align to right) */
    for (int i = 0; i < ndim2; i++) {
        int result_idx = i + result_ndim - ndim2;
        int32_t dim = shape2[i];

        if (out_shape[result_idx] == 1) {
            out_shape[result_idx] = dim;
        }
        else if (dim != 1 && dim != out_shape[result_idx]) {
            return -1; /* Incompatible shapes */
        }
    }

    *out_ndim = result_ndim;
    return 0;
}

/**
 * Compute broadcast shape for multiple shapes.
 */
EXPORT int broadcast_shapes_multi(const int32_t** shapes, const int32_t* ndims,
                                   int32_t num_arrays,
                                   int32_t* out_shape, int32_t* out_ndim)
{
    if (num_arrays == 0) {
        *out_ndim = 0;
        return 0;
    }

    if (num_arrays == 1) {
        memcpy(out_shape, shapes[0], ndims[0] * sizeof(int32_t));
        *out_ndim = ndims[0];
        return 0;
    }

    /* Start with first shape */
    int32_t current_shape[NPY_MAXDIMS];
    int32_t current_ndim;

    memcpy(current_shape, shapes[0], ndims[0] * sizeof(int32_t));
    current_ndim = ndims[0];

    /* Progressively broadcast with each additional shape */
    for (int i = 1; i < num_arrays; i++) {
        int32_t new_shape[NPY_MAXDIMS];
        int32_t new_ndim;

        if (broadcast_shapes(current_shape, current_ndim,
                             shapes[i], ndims[i],
                             new_shape, &new_ndim) != 0) {
            return -1;
        }

        memcpy(current_shape, new_shape, new_ndim * sizeof(int32_t));
        current_ndim = new_ndim;
    }

    memcpy(out_shape, current_shape, current_ndim * sizeof(int32_t));
    *out_ndim = current_ndim;
    return 0;
}

/**
 * Compute strides for broadcasting an array to a target shape.
 *
 * Adapted from NumPy nditer_constr.c lines 1594-1615:
 *
 * For each dimension in target shape:
 *   - If dimension doesn't exist in operand (prepended 1): stride = 0
 *   - If operand dimension is 1 (broadcasting): stride = 0
 *   - Otherwise: use original stride
 *
 * The stride=0 trick means the pointer doesn't advance in that dimension,
 * effectively repeating the same value across the broadcast extent.
 */
EXPORT int broadcast_strides(NDArray* arr, const int32_t* target_shape, int32_t target_ndim,
                              int32_t* out_strides)
{
    if (!arr || !target_shape || !out_strides) return -1;

    int diff = target_ndim - arr->ndim;

    /* Can't broadcast to fewer dimensions */
    if (diff < 0) return -1;

    for (int i = 0; i < target_ndim; i++) {
        int arr_idx = i - diff;

        if (arr_idx < 0) {
            /* Prepended dimension: stride = 0 (broadcast) */
            out_strides[i] = 0;
        }
        else if (arr->shape[arr_idx] == 1) {
            /* Size-1 dimension being broadcast: stride = 0 */
            out_strides[i] = 0;
        }
        else if (arr->shape[arr_idx] == target_shape[i]) {
            /* Matching dimension: use original stride */
            out_strides[i] = arr->strides[arr_idx];
        }
        else {
            /* Incompatible: operand dimension != 1 and != target */
            return -1;
        }
    }

    return 0;
}

/**
 * Create a broadcast view of an array.
 */
EXPORT NDArray* ndarray_broadcast_to(NDArray* arr, const int32_t* target_shape, int32_t target_ndim)
{
    if (!arr || !target_shape) return NULL;

    /* Validate broadcast compatibility */
    int diff = target_ndim - arr->ndim;
    if (diff < 0) return NULL;

    for (int i = 0; i < arr->ndim; i++) {
        int target_idx = i + diff;
        if (arr->shape[i] != 1 && arr->shape[i] != target_shape[target_idx]) {
            return NULL;
        }
    }

    /* Compute broadcast strides */
    int32_t new_strides[NPY_MAXDIMS];
    if (broadcast_strides(arr, target_shape, target_ndim, new_strides) != 0) {
        return NULL;
    }

    /* Create view with broadcast shape and strides */
    NDArray* view = ndarray_view(arr, target_ndim, (int32_t*)target_shape, new_strides);

    /* Mark as not writeable since writes would be repeated */
    if (view) {
        view->flags &= ~NDARRAY_WRITEABLE;
    }

    return view;
}

/**
 * Check if two shapes are broadcast-compatible.
 */
EXPORT int shapes_are_broadcastable(const int32_t* shape1, int32_t ndim1,
                                     const int32_t* shape2, int32_t ndim2)
{
    int32_t out_shape[NPY_MAXDIMS];
    int32_t out_ndim;
    return broadcast_shapes(shape1, ndim1, shape2, ndim2, out_shape, &out_ndim) == 0;
}
