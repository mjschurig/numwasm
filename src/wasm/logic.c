/**
 * NumJS-WASM Logic & Comparison Functions Implementation
 *
 * Provides truth testing, element-wise predicates, and array comparison operations.
 *
 * Reference implementations:
 * - numpy/_core/fromnumeric.py (all, any)
 * - numpy/_core/numeric.py (isclose, allclose, array_equal, array_equiv)
 * - numpy/lib/_ufunclike_impl.py (isneginf, isposinf)
 */

#include "logic.h"
#include "broadcast.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#define EXPORT EMSCRIPTEN_KEEPALIVE
#else
#define EXPORT
#endif

/* ============ Element-wise Predicates ============ */

EXPORT NDArray* ndarray_isfinite(NDArray* arr)
{
    if (!arr) return NULL;

    /* Create bool result with same shape */
    NDArray* result = ndarray_empty(arr->ndim, arr->shape, DTYPE_BOOL);
    if (!result) return NULL;

    uint8_t* out = (uint8_t*)result->data;

    for (size_t i = 0; i < arr->size; i++) {
        double val = ndarray_get_flat(arr, i);
        out[i] = isfinite(val) ? 1 : 0;
    }

    return result;
}

EXPORT NDArray* ndarray_isinf(NDArray* arr)
{
    if (!arr) return NULL;

    NDArray* result = ndarray_empty(arr->ndim, arr->shape, DTYPE_BOOL);
    if (!result) return NULL;

    uint8_t* out = (uint8_t*)result->data;

    for (size_t i = 0; i < arr->size; i++) {
        double val = ndarray_get_flat(arr, i);
        out[i] = isinf(val) ? 1 : 0;
    }

    return result;
}

EXPORT NDArray* ndarray_isnan(NDArray* arr)
{
    if (!arr) return NULL;

    NDArray* result = ndarray_empty(arr->ndim, arr->shape, DTYPE_BOOL);
    if (!result) return NULL;

    uint8_t* out = (uint8_t*)result->data;

    for (size_t i = 0; i < arr->size; i++) {
        double val = ndarray_get_flat(arr, i);
        out[i] = isnan(val) ? 1 : 0;
    }

    return result;
}

EXPORT NDArray* ndarray_isneginf(NDArray* arr)
{
    if (!arr) return NULL;

    NDArray* result = ndarray_empty(arr->ndim, arr->shape, DTYPE_BOOL);
    if (!result) return NULL;

    uint8_t* out = (uint8_t*)result->data;

    for (size_t i = 0; i < arr->size; i++) {
        double val = ndarray_get_flat(arr, i);
        out[i] = (isinf(val) && val < 0) ? 1 : 0;
    }

    return result;
}

EXPORT NDArray* ndarray_isposinf(NDArray* arr)
{
    if (!arr) return NULL;

    NDArray* result = ndarray_empty(arr->ndim, arr->shape, DTYPE_BOOL);
    if (!result) return NULL;

    uint8_t* out = (uint8_t*)result->data;

    for (size_t i = 0; i < arr->size; i++) {
        double val = ndarray_get_flat(arr, i);
        out[i] = (isinf(val) && val > 0) ? 1 : 0;
    }

    return result;
}

EXPORT NDArray* ndarray_iscomplex_elem(NDArray* arr)
{
    if (!arr) return NULL;

    NDArray* result = ndarray_empty(arr->ndim, arr->shape, DTYPE_BOOL);
    if (!result) return NULL;

    uint8_t* out = (uint8_t*)result->data;

    /* Non-complex types: all False */
    if (arr->dtype != DTYPE_COMPLEX64 && arr->dtype != DTYPE_COMPLEX128) {
        memset(out, 0, result->size);
        return result;
    }

    /* Complex: check imag != 0 */
    for (size_t i = 0; i < arr->size; i++) {
        double imag = ndarray_get_complex_imag(arr, i);
        out[i] = (imag != 0.0) ? 1 : 0;
    }

    return result;
}

EXPORT NDArray* ndarray_isreal_elem(NDArray* arr)
{
    if (!arr) return NULL;

    NDArray* result = ndarray_empty(arr->ndim, arr->shape, DTYPE_BOOL);
    if (!result) return NULL;

    uint8_t* out = (uint8_t*)result->data;

    /* Non-complex types: all True */
    if (arr->dtype != DTYPE_COMPLEX64 && arr->dtype != DTYPE_COMPLEX128) {
        memset(out, 1, result->size);
        return result;
    }

    /* Complex: check imag == 0 */
    for (size_t i = 0; i < arr->size; i++) {
        double imag = ndarray_get_complex_imag(arr, i);
        out[i] = (imag == 0.0) ? 1 : 0;
    }

    return result;
}

/* ============ Reductions: all / any ============ */

EXPORT NDArray* ndarray_all(NDArray* arr)
{
    if (!arr) return NULL;

    /* Create 0-d bool result, initialized to False */
    NDArray* result = ndarray_scalar(0.0, DTYPE_BOOL);
    if (!result) return NULL;

    /* Check all elements */
    for (size_t i = 0; i < arr->size; i++) {
        double val = ndarray_get_flat(arr, i);
        if (val == 0.0) {
            /* Found a falsy value - return False */
            return result;
        }
    }

    /* All truthy - return True */
    ndarray_set_flat(result, 0, 1.0);
    return result;
}

EXPORT NDArray* ndarray_all_axis(NDArray* arr, int32_t axis, int keepdims)
{
    if (!arr) return NULL;

    /* Normalize negative axis */
    if (axis < 0) axis += arr->ndim;
    if (axis < 0 || axis >= arr->ndim) return NULL;

    /* Compute result shape */
    int32_t result_shape[NPY_MAXDIMS];
    int32_t result_ndim;

    if (keepdims) {
        result_ndim = arr->ndim;
        for (int i = 0; i < arr->ndim; i++) {
            result_shape[i] = (i == axis) ? 1 : arr->shape[i];
        }
    } else {
        result_ndim = arr->ndim - 1;
        int j = 0;
        for (int i = 0; i < arr->ndim; i++) {
            if (i != axis) {
                result_shape[j++] = arr->shape[i];
            }
        }
    }

    /* Handle 0-d result (when arr is 1-d) */
    if (result_ndim == 0) {
        return ndarray_all(arr);
    }

    /* Initialize result to True (1) */
    NDArray* result = ndarray_full(result_ndim, result_shape, DTYPE_BOOL, 1.0);
    if (!result) return NULL;

    /* Iterate and reduce along axis */
    int32_t indices[NPY_MAXDIMS] = {0};
    int32_t result_indices[NPY_MAXDIMS];

    for (size_t flat = 0; flat < arr->size; flat++) {
        /* Map arr indices to result indices (skipping axis) */
        int j = 0;
        for (int i = 0; i < arr->ndim; i++) {
            if (keepdims) {
                result_indices[i] = (i == axis) ? 0 : indices[i];
            } else if (i != axis) {
                result_indices[j++] = indices[i];
            }
        }

        /* Get current value and result value */
        double val = ndarray_get_item(arr, indices, arr->ndim);
        double cur = ndarray_get_item(result, result_indices, result_ndim);

        /* Logical AND: if val is 0, result becomes 0 */
        if (cur != 0.0 && val == 0.0) {
            ndarray_set_item(result, result_indices, result_ndim, 0.0);
        }

        /* Increment indices */
        for (int d = arr->ndim - 1; d >= 0; d--) {
            indices[d]++;
            if (indices[d] < arr->shape[d]) break;
            indices[d] = 0;
        }
    }

    return result;
}

EXPORT NDArray* ndarray_any(NDArray* arr)
{
    if (!arr) return NULL;

    /* Create 0-d bool result, initialized to False */
    NDArray* result = ndarray_scalar(0.0, DTYPE_BOOL);
    if (!result) return NULL;

    /* Check all elements */
    for (size_t i = 0; i < arr->size; i++) {
        double val = ndarray_get_flat(arr, i);
        if (val != 0.0) {
            /* Found a truthy value - return True */
            ndarray_set_flat(result, 0, 1.0);
            return result;
        }
    }

    /* All falsy - return False */
    return result;
}

EXPORT NDArray* ndarray_any_axis(NDArray* arr, int32_t axis, int keepdims)
{
    if (!arr) return NULL;

    /* Normalize negative axis */
    if (axis < 0) axis += arr->ndim;
    if (axis < 0 || axis >= arr->ndim) return NULL;

    /* Compute result shape */
    int32_t result_shape[NPY_MAXDIMS];
    int32_t result_ndim;

    if (keepdims) {
        result_ndim = arr->ndim;
        for (int i = 0; i < arr->ndim; i++) {
            result_shape[i] = (i == axis) ? 1 : arr->shape[i];
        }
    } else {
        result_ndim = arr->ndim - 1;
        int j = 0;
        for (int i = 0; i < arr->ndim; i++) {
            if (i != axis) {
                result_shape[j++] = arr->shape[i];
            }
        }
    }

    /* Handle 0-d result (when arr is 1-d) */
    if (result_ndim == 0) {
        return ndarray_any(arr);
    }

    /* Initialize result to False (0) */
    NDArray* result = ndarray_create(result_ndim, result_shape, DTYPE_BOOL);
    if (!result) return NULL;

    /* Iterate and reduce along axis */
    int32_t indices[NPY_MAXDIMS] = {0};
    int32_t result_indices[NPY_MAXDIMS];

    for (size_t flat = 0; flat < arr->size; flat++) {
        /* Map arr indices to result indices (skipping axis) */
        int j = 0;
        for (int i = 0; i < arr->ndim; i++) {
            if (keepdims) {
                result_indices[i] = (i == axis) ? 0 : indices[i];
            } else if (i != axis) {
                result_indices[j++] = indices[i];
            }
        }

        /* Get current value and result value */
        double val = ndarray_get_item(arr, indices, arr->ndim);
        double cur = ndarray_get_item(result, result_indices, result_ndim);

        /* Logical OR: if val is non-zero, result becomes 1 */
        if (cur == 0.0 && val != 0.0) {
            ndarray_set_item(result, result_indices, result_ndim, 1.0);
        }

        /* Increment indices */
        for (int d = arr->ndim - 1; d >= 0; d--) {
            indices[d]++;
            if (indices[d] < arr->shape[d]) break;
            indices[d] = 0;
        }
    }

    return result;
}

/* ============ Comparison Functions ============ */

EXPORT NDArray* ndarray_isclose(NDArray* a, NDArray* b, double rtol, double atol, int equal_nan)
{
    if (!a || !b) return NULL;

    /* Compute broadcast shape */
    int32_t bc_shape[NPY_MAXDIMS];
    int32_t bc_ndim;

    if (broadcast_shapes(a->shape, a->ndim, b->shape, b->ndim, bc_shape, &bc_ndim) != 0) {
        return NULL;
    }

    /* Create result array */
    NDArray* result = ndarray_empty(bc_ndim, bc_shape, DTYPE_BOOL);
    if (!result) return NULL;

    /* Compute broadcast strides */
    int32_t a_strides[NPY_MAXDIMS];
    int32_t b_strides[NPY_MAXDIMS];

    if (broadcast_strides(a, bc_shape, bc_ndim, a_strides) != 0 ||
        broadcast_strides(b, bc_shape, bc_ndim, b_strides) != 0) {
        ndarray_free(result);
        return NULL;
    }

    size_t a_elem_size = dtype_size(a->dtype);
    size_t b_elem_size = dtype_size(b->dtype);

    uint8_t* out = (uint8_t*)result->data;
    int32_t indices[NPY_MAXDIMS] = {0};

    for (size_t i = 0; i < result->size; i++) {
        /* Compute flat offsets using broadcast strides */
        size_t a_offset = 0;
        size_t b_offset = 0;

        for (int d = 0; d < bc_ndim; d++) {
            if (a_strides[d] != 0) {
                a_offset += (size_t)indices[d] * ((size_t)a_strides[d] / a_elem_size);
            }
            if (b_strides[d] != 0) {
                b_offset += (size_t)indices[d] * ((size_t)b_strides[d] / b_elem_size);
            }
        }

        double val_a = ndarray_get_flat(a, a_offset);
        double val_b = ndarray_get_flat(b, b_offset);

        /* Check for NaN equality */
        int a_nan = isnan(val_a);
        int b_nan = isnan(val_b);

        if (equal_nan && a_nan && b_nan) {
            out[i] = 1;
        } else if (a_nan || b_nan) {
            out[i] = 0;
        } else if (val_a == val_b) {
            /* Exact equality (handles inf == inf) */
            out[i] = 1;
        } else if (isfinite(val_b)) {
            /* Tolerance comparison: |a - b| <= atol + rtol * |b| */
            double diff = fabs(val_a - val_b);
            double threshold = atol + rtol * fabs(val_b);
            out[i] = (diff <= threshold) ? 1 : 0;
        } else {
            /* b is inf but a is not exactly equal */
            out[i] = 0;
        }

        /* Increment indices */
        for (int d = bc_ndim - 1; d >= 0; d--) {
            indices[d]++;
            if (indices[d] < bc_shape[d]) break;
            indices[d] = 0;
        }
    }

    return result;
}

EXPORT int ndarray_allclose(NDArray* a, NDArray* b, double rtol, double atol, int equal_nan)
{
    NDArray* close = ndarray_isclose(a, b, rtol, atol, equal_nan);
    if (!close) return 0;

    NDArray* all_result = ndarray_all(close);
    ndarray_free(close);

    if (!all_result) return 0;

    int result = (ndarray_get_flat(all_result, 0) != 0.0) ? 1 : 0;
    ndarray_free(all_result);

    return result;
}

EXPORT int ndarray_array_equal(NDArray* a1, NDArray* a2, int equal_nan)
{
    if (!a1 || !a2) return 0;

    /* Check shape equality */
    if (a1->ndim != a2->ndim) return 0;

    for (int i = 0; i < a1->ndim; i++) {
        if (a1->shape[i] != a2->shape[i]) return 0;
    }

    /* Check element equality */
    for (size_t i = 0; i < a1->size; i++) {
        double v1 = ndarray_get_flat(a1, i);
        double v2 = ndarray_get_flat(a2, i);

        int v1_nan = isnan(v1);
        int v2_nan = isnan(v2);

        if (v1_nan && v2_nan) {
            if (!equal_nan) return 0;
            /* else: NaN == NaN is OK when equal_nan is true */
        } else if (v1_nan || v2_nan) {
            return 0;
        } else if (v1 != v2) {
            return 0;
        }
    }

    return 1;
}

EXPORT int ndarray_array_equiv(NDArray* a1, NDArray* a2)
{
    if (!a1 || !a2) return 0;

    /* Check broadcast compatibility */
    int32_t bc_shape[NPY_MAXDIMS];
    int32_t bc_ndim;

    if (broadcast_shapes(a1->shape, a1->ndim, a2->shape, a2->ndim, bc_shape, &bc_ndim) != 0) {
        return 0;
    }

    /* Compute broadcast strides */
    int32_t a1_strides[NPY_MAXDIMS];
    int32_t a2_strides[NPY_MAXDIMS];

    if (broadcast_strides(a1, bc_shape, bc_ndim, a1_strides) != 0 ||
        broadcast_strides(a2, bc_shape, bc_ndim, a2_strides) != 0) {
        return 0;
    }

    size_t a1_elem = dtype_size(a1->dtype);
    size_t a2_elem = dtype_size(a2->dtype);

    int32_t indices[NPY_MAXDIMS] = {0};
    size_t total = 1;
    for (int i = 0; i < bc_ndim; i++) total *= (size_t)bc_shape[i];

    for (size_t i = 0; i < total; i++) {
        size_t o1 = 0, o2 = 0;

        for (int d = 0; d < bc_ndim; d++) {
            if (a1_strides[d] != 0) {
                o1 += (size_t)indices[d] * ((size_t)a1_strides[d] / a1_elem);
            }
            if (a2_strides[d] != 0) {
                o2 += (size_t)indices[d] * ((size_t)a2_strides[d] / a2_elem);
            }
        }

        double v1 = ndarray_get_flat(a1, o1);
        double v2 = ndarray_get_flat(a2, o2);

        if (v1 != v2) return 0;

        /* Increment indices */
        for (int d = bc_ndim - 1; d >= 0; d--) {
            indices[d]++;
            if (indices[d] < bc_shape[d]) break;
            indices[d] = 0;
        }
    }

    return 1;
}
