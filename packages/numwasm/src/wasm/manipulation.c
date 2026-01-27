#include "manipulation.h"
#include <stdlib.h>
#include <string.h>

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#define EXPORT EMSCRIPTEN_KEEPALIVE
#else
#define EXPORT
#endif

/**
 * Concatenate arrays along an existing axis.
 *
 * Based on NumPy's concatenate implementation.
 * All arrays must have the same shape except along the concatenation axis.
 */
EXPORT NDArray* ndarray_concatenate(NDArray** arrays, int32_t n_arrays, int32_t axis)
{
    if (!arrays || n_arrays == 0) return NULL;

    NDArray* first = arrays[0];
    if (!first) return NULL;

    /* Normalize axis */
    int32_t ndim = first->ndim;
    if (axis < 0) axis += ndim;
    if (axis < 0 || axis >= ndim) return NULL;

    /* Validate all arrays have same shape except along axis */
    int32_t total_axis_size = first->shape[axis];
    DType result_dtype = first->dtype;

    for (int i = 1; i < n_arrays; i++) {
        if (!arrays[i]) return NULL;
        if (arrays[i]->ndim != ndim) return NULL;

        for (int d = 0; d < ndim; d++) {
            if (d != axis && arrays[i]->shape[d] != first->shape[d]) {
                return NULL;  /* Shape mismatch */
            }
        }
        total_axis_size += arrays[i]->shape[axis];

        /* Use first array's dtype (caller handles type promotion) */
    }

    /* Create output shape */
    int32_t* out_shape = (int32_t*)malloc((size_t)ndim * sizeof(int32_t));
    if (!out_shape) return NULL;
    memcpy(out_shape, first->shape, (size_t)ndim * sizeof(int32_t));
    out_shape[axis] = total_axis_size;

    /* Create output array */
    NDArray* result = ndarray_empty(ndim, out_shape, result_dtype);
    free(out_shape);
    if (!result) return NULL;

    /* Copy data from each input array */
    size_t axis_offset = 0;
    for (int i = 0; i < n_arrays; i++) {
        NDArray* arr = arrays[i];
        size_t arr_size = arr->size;

        /* For each element in the source array */
        for (size_t flat = 0; flat < arr_size; flat++) {
            /* Convert flat index to multi-index in source */
            int32_t multi_idx[32];
            size_t remainder = flat;
            for (int d = ndim - 1; d >= 0; d--) {
                multi_idx[d] = (int32_t)(remainder % (size_t)arr->shape[d]);
                remainder /= (size_t)arr->shape[d];
            }

            /* Adjust axis index for destination */
            multi_idx[axis] += (int32_t)axis_offset;

            /* Convert multi-index to flat index in result (C-contiguous) */
            size_t result_flat = 0;
            size_t multiplier = 1;
            for (int d = ndim - 1; d >= 0; d--) {
                result_flat += (size_t)multi_idx[d] * multiplier;
                multiplier *= (size_t)result->shape[d];
            }

            /* Copy value */
            double val = ndarray_get_flat(arr, flat);
            ndarray_set_flat(result, result_flat, val);
        }

        axis_offset += (size_t)arr->shape[axis];
    }

    return result;
}
