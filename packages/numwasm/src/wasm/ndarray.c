/*
 * NumJS NDArray Implementation
 *
 * A minimal n-dimensional array for WebAssembly, inspired by NumPy.
 * Uses NumPy's pairwise summation algorithm for accurate reductions.
 */

#include "ndarray.h"
#include "dtype.h"
#include "pairwise_sum.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* Mark functions for export when compiled with Emscripten */
#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#define EXPORT EMSCRIPTEN_KEEPALIVE
#else
#define EXPORT
#endif

/* ============ Utility Functions ============ */

static size_t compute_size(int32_t ndim, int32_t* shape)
{
    if (ndim == 0) return 1; /* Scalar */
    size_t size = 1;
    for (int32_t i = 0; i < ndim; i++) {
        size *= (size_t)shape[i];
    }
    return size;
}

/*
 * Compute C-contiguous (row-major) strides.
 * This matches NumPy's default memory layout.
 */
static void compute_strides(int32_t ndim, int32_t* shape, DType dtype, int32_t* strides)
{
    size_t elem_size = dtype_size(dtype);
    if (ndim == 0) {
        return;
    }
    strides[ndim - 1] = (int32_t)elem_size;
    for (int32_t i = ndim - 2; i >= 0; i--) {
        strides[i] = strides[i + 1] * shape[i + 1];
    }
}

/* ============ Contiguity Checks ============ */

EXPORT bool ndarray_is_c_contiguous(const NDArray* arr)
{
    if (!arr) return false;
    if (arr->ndim == 0) return true;
    if (arr->size == 0) return true;

    size_t elem_size = dtype_size(arr->dtype);
    size_t expected = elem_size;

    for (int32_t i = arr->ndim - 1; i >= 0; i--) {
        if (arr->shape[i] == 0) return true;
        if (arr->shape[i] != 1 && (size_t)arr->strides[i] != expected) {
            return false;
        }
        expected *= (size_t)arr->shape[i];
    }
    return true;
}

EXPORT bool ndarray_is_f_contiguous(NDArray* arr)
{
    if (!arr) return false;
    if (arr->ndim == 0) return true;
    if (arr->size == 0) return true;

    size_t elem_size = dtype_size(arr->dtype);
    size_t expected = elem_size;

    for (int32_t i = 0; i < arr->ndim; i++) {
        if (arr->shape[i] == 0) return true;
        if (arr->shape[i] != 1 && (size_t)arr->strides[i] != expected) {
            return false;
        }
        expected *= (size_t)arr->shape[i];
    }
    return true;
}

EXPORT void ndarray_update_flags(NDArray* arr)
{
    if (!arr) return;

    /* Clear contiguity flags */
    arr->flags &= ~(NDARRAY_C_CONTIGUOUS | NDARRAY_F_CONTIGUOUS);

    /* Set based on current layout */
    if (ndarray_is_c_contiguous(arr)) {
        arr->flags |= NDARRAY_C_CONTIGUOUS;
    }
    if (ndarray_is_f_contiguous(arr)) {
        arr->flags |= NDARRAY_F_CONTIGUOUS;
    }
}

/* ============ Array Creation ============ */

/*
 * Internal helper to allocate array structure without data.
 */
static NDArray* ndarray_alloc_struct(int32_t ndim, int32_t* shape, DType dtype)
{
    NDArray* arr = (NDArray*)malloc(sizeof(NDArray));
    if (!arr) return NULL;

    arr->ndim = ndim;
    arr->dtype = dtype;
    arr->size = compute_size(ndim, shape);
    arr->flags = NDARRAY_OWNDATA | NDARRAY_WRITEABLE;
    arr->base = NULL;
    arr->data = NULL;

    /* Allocate shape and strides arrays */
    if (ndim > 0) {
        arr->shape = (int32_t*)malloc((size_t)ndim * sizeof(int32_t));
        arr->strides = (int32_t*)malloc((size_t)ndim * sizeof(int32_t));
        if (!arr->shape || !arr->strides) {
            free(arr->shape);
            free(arr->strides);
            free(arr);
            return NULL;
        }
        memcpy(arr->shape, shape, (size_t)ndim * sizeof(int32_t));
        compute_strides(ndim, shape, dtype, arr->strides);
    } else {
        arr->shape = NULL;
        arr->strides = NULL;
    }

    return arr;
}

EXPORT NDArray* ndarray_create(int32_t ndim, int32_t* shape, DType dtype)
{
    NDArray* arr = ndarray_alloc_struct(ndim, shape, dtype);
    if (!arr) return NULL;

    /* Allocate and zero-initialize data buffer */
    size_t data_size = arr->size * dtype_size(dtype);
    if (data_size > 0) {
        arr->data = malloc(data_size);
        if (!arr->data) {
            free(arr->shape);
            free(arr->strides);
            free(arr);
            return NULL;
        }
        memset(arr->data, 0, data_size);
    }

    ndarray_update_flags(arr);
    return arr;
}

EXPORT NDArray* ndarray_empty(int32_t ndim, int32_t* shape, DType dtype)
{
    NDArray* arr = ndarray_alloc_struct(ndim, shape, dtype);
    if (!arr) return NULL;

    /* Allocate data buffer WITHOUT initialization (faster) */
    size_t data_size = arr->size * dtype_size(dtype);
    if (data_size > 0) {
        arr->data = malloc(data_size);
        if (!arr->data) {
            free(arr->shape);
            free(arr->strides);
            free(arr);
            return NULL;
        }
        /* NO memset - leave uninitialized */
    }

    ndarray_update_flags(arr);
    return arr;
}

EXPORT NDArray* ndarray_full(int32_t ndim, int32_t* shape, DType dtype, double value)
{
    NDArray* arr = ndarray_empty(ndim, shape, dtype);
    if (!arr) return NULL;

    ndarray_fill(arr, value);
    return arr;
}

EXPORT NDArray* ndarray_arange(double start, double stop, double step, DType dtype)
{
    if (step == 0.0) return NULL;

    /* Calculate number of elements */
    int32_t num;
    if (step > 0) {
        num = (stop > start) ? (int32_t)ceil((stop - start) / step) : 0;
    } else {
        num = (start > stop) ? (int32_t)ceil((start - stop) / (-step)) : 0;
    }

    if (num <= 0) {
        int32_t shape[] = {0};
        return ndarray_empty(1, shape, dtype);
    }

    int32_t shape[] = {num};
    NDArray* arr = ndarray_empty(1, shape, dtype);
    if (!arr) return NULL;

    for (int32_t i = 0; i < num; i++) {
        ndarray_set_flat(arr, i, start + i * step);
    }

    return arr;
}

EXPORT NDArray* ndarray_linspace(double start, double stop, int32_t num, int32_t endpoint, DType dtype)
{
    if (num < 0) return NULL;

    if (num == 0) {
        int32_t shape[] = {0};
        return ndarray_empty(1, shape, dtype);
    }

    int32_t shape[] = {num};
    NDArray* arr = ndarray_empty(1, shape, dtype);
    if (!arr) return NULL;

    if (num == 1) {
        ndarray_set_flat(arr, 0, start);
        return arr;
    }

    double div = endpoint ? (double)(num - 1) : (double)num;
    double step = (stop - start) / div;

    for (int32_t i = 0; i < num; i++) {
        ndarray_set_flat(arr, i, start + i * step);
    }

    return arr;
}

EXPORT NDArray* ndarray_eye(int32_t N, int32_t M, int32_t k, DType dtype)
{
    if (N < 0 || M < 0) return NULL;

    int32_t shape[] = {N, M};
    NDArray* arr = ndarray_create(2, shape, dtype);  /* zeros */
    if (!arr) return NULL;

    /* Fill diagonal with ones */
    for (int32_t i = 0; i < N; i++) {
        int32_t j = i + k;
        if (j >= 0 && j < M) {
            size_t flat_idx = (size_t)i * M + j;
            ndarray_set_flat(arr, flat_idx, 1.0);
        }
    }

    return arr;
}

EXPORT NDArray* ndarray_tri(int32_t N, int32_t M, int32_t k, DType dtype)
{
    if (N < 0 || M < 0) return NULL;

    int32_t shape[] = {N, M};
    NDArray* arr = ndarray_create(2, shape, dtype);  /* zeros */
    if (!arr) return NULL;

    /* Fill lower triangle with ones */
    for (int32_t i = 0; i < N; i++) {
        int32_t max_j = i + k + 1;
        if (max_j > M) max_j = M;
        for (int32_t j = 0; j < max_j; j++) {
            size_t flat_idx = (size_t)i * M + j;
            ndarray_set_flat(arr, flat_idx, 1.0);
        }
    }

    return arr;
}

EXPORT NDArray* ndarray_tril(const NDArray* arr, int32_t k)
{
    if (!arr || arr->ndim != 2) return NULL;

    NDArray* result = ndarray_copy(arr);
    if (!result) return NULL;

    int32_t rows = arr->shape[0];
    int32_t cols = arr->shape[1];

    /* Zero elements above the k-th diagonal */
    for (int32_t i = 0; i < rows; i++) {
        for (int32_t j = i + k + 1; j < cols; j++) {
            size_t flat_idx = (size_t)i * cols + j;
            ndarray_set_flat(result, flat_idx, 0.0);
        }
    }

    return result;
}

EXPORT NDArray* ndarray_triu(const NDArray* arr, int32_t k)
{
    if (!arr || arr->ndim != 2) return NULL;

    NDArray* result = ndarray_copy(arr);
    if (!result) return NULL;

    int32_t rows = arr->shape[0];
    int32_t cols = arr->shape[1];

    /* Zero elements below the k-th diagonal */
    for (int32_t i = 0; i < rows; i++) {
        int32_t max_j = i + k;
        if (max_j > cols) max_j = cols;
        for (int32_t j = 0; j < max_j; j++) {
            size_t flat_idx = (size_t)i * cols + j;
            ndarray_set_flat(result, flat_idx, 0.0);
        }
    }

    return result;
}

EXPORT NDArray* ndarray_vander(const NDArray* x, int32_t N, int32_t increasing)
{
    if (!x || x->ndim != 1) return NULL;

    int32_t rows = (int32_t)x->size;
    if (N <= 0) N = rows;

    int32_t shape[] = {rows, N};
    NDArray* result = ndarray_empty(2, shape, x->dtype);
    if (!result) return NULL;

    for (int32_t i = 0; i < rows; i++) {
        double val = ndarray_get_flat(x, i);
        for (int32_t j = 0; j < N; j++) {
            int32_t power = increasing ? j : (N - 1 - j);
            double elem = pow(val, (double)power);
            size_t flat_idx = (size_t)i * N + j;
            ndarray_set_flat(result, flat_idx, elem);
        }
    }

    return result;
}

EXPORT NDArray* ndarray_scalar(double value, DType dtype)
{
    NDArray* arr = ndarray_empty(0, NULL, dtype);
    if (!arr) return NULL;

    /* Allocate single element */
    size_t elem_size = dtype_size(dtype);
    arr->data = malloc(elem_size);
    if (!arr->data) {
        free(arr);
        return NULL;
    }

    ndarray_set_flat(arr, 0, value);
    ndarray_update_flags(arr);
    return arr;
}

EXPORT NDArray* ndarray_from_data(void* data, int32_t ndim, int32_t* shape, DType dtype)
{
    NDArray* arr = ndarray_create(ndim, shape, dtype);
    if (!arr) return NULL;

    /* Copy data into the array's buffer */
    size_t data_size = arr->size * dtype_size(dtype);
    if (data_size > 0 && data != NULL) {
        memcpy(arr->data, data, data_size);
    }

    return arr;
}

EXPORT void ndarray_free(NDArray* arr)
{
    if (!arr) return;

    if (arr->flags & NDARRAY_OWNDATA) {
        free(arr->data);
    }
    free(arr->shape);
    free(arr->strides);
    free(arr);
}

EXPORT NDArray* ndarray_copy(const NDArray* arr)
{
    if (!arr) return NULL;

    NDArray* copy = ndarray_empty(arr->ndim, arr->shape, arr->dtype);
    if (!copy) return NULL;

    /* Copy data */
    if (ndarray_is_c_contiguous(arr)) {
        /* Fast path: memcpy for contiguous arrays */
        memcpy(copy->data, arr->data, arr->size * dtype_size(arr->dtype));
    } else {
        /* Slow path: element-by-element for non-contiguous */
        for (size_t i = 0; i < arr->size; i++) {
            double val = ndarray_get_flat(arr, i);
            ndarray_set_flat(copy, i, val);
        }
    }

    return copy;
}

EXPORT int32_t ndarray_copyto(NDArray* dst, const NDArray* src, const NDArray* where)
{
    if (!dst || !src) return -1;
    if (dst->size != src->size) return -1;
    if (where && where->size != dst->size) return -1;

    /* Fast path: both contiguous, same dtype, no mask */
    if (!where && dst->dtype == src->dtype &&
        ndarray_is_c_contiguous(dst) && ndarray_is_c_contiguous(src)) {
        memcpy(dst->data, src->data, dst->size * dtype_size(dst->dtype));
        return 0;
    }

    /* General path: element-by-element */
    for (size_t i = 0; i < dst->size; i++) {
        /* Check mask if provided */
        if (where) {
            double mask_val = ndarray_get_flat(where, i);
            if (mask_val == 0.0) continue; /* Skip if mask is false */
        }

        double val = ndarray_get_flat(src, i);
        ndarray_set_flat(dst, i, val);
    }

    return 0;
}

EXPORT NDArray* ndarray_astype(NDArray* arr, DType dtype)
{
    if (!arr) return NULL;

    /* Same dtype: just copy */
    if (arr->dtype == dtype) {
        return ndarray_copy(arr);
    }

    NDArray* result = ndarray_empty(arr->ndim, arr->shape, dtype);
    if (!result) return NULL;

    /* Convert each element */
    for (size_t i = 0; i < arr->size; i++) {
        double val = ndarray_get_flat(arr, i);
        ndarray_set_flat(result, i, val);
    }

    return result;
}

/* ============ Views ============ */

EXPORT NDArray* ndarray_view(NDArray* src, int32_t ndim, int32_t* shape, int32_t* strides)
{
    return ndarray_view_with_offset(src, ndim, shape, strides, 0);
}

EXPORT NDArray* ndarray_view_with_offset(NDArray* src, int32_t ndim, int32_t* shape,
                                          int32_t* strides, size_t byte_offset)
{
    if (!src) return NULL;

    NDArray* view = (NDArray*)malloc(sizeof(NDArray));
    if (!view) return NULL;

    view->ndim = ndim;
    view->dtype = src->dtype;
    view->flags = NDARRAY_WRITEABLE;  /* No OWNDATA flag for views */
    view->size = compute_size(ndim, shape);
    view->data = (char*)src->data + byte_offset;

    /* Find ultimate base (for chained views) */
    view->base = (src->base != NULL) ? src->base : src;

    /* Allocate and copy shape/strides */
    if (ndim > 0) {
        view->shape = (int32_t*)malloc((size_t)ndim * sizeof(int32_t));
        view->strides = (int32_t*)malloc((size_t)ndim * sizeof(int32_t));
        if (!view->shape || !view->strides) {
            free(view->shape);
            free(view->strides);
            free(view);
            return NULL;
        }
        memcpy(view->shape, shape, (size_t)ndim * sizeof(int32_t));
        memcpy(view->strides, strides, (size_t)ndim * sizeof(int32_t));
    } else {
        view->shape = NULL;
        view->strides = NULL;
    }

    /* Update contiguity flags */
    ndarray_update_flags(view);

    return view;
}

/* ============ Shape Manipulation ============ */

EXPORT NDArray* ndarray_reshape(NDArray* arr, int32_t new_ndim, int32_t* new_shape)
{
    if (!arr || !new_shape) return NULL;

    /* Resolve -1 dimension */
    int32_t resolved[32];
    int unknown_idx = -1;
    size_t known_product = 1;

    for (int i = 0; i < new_ndim; i++) {
        if (new_shape[i] == -1) {
            if (unknown_idx != -1) return NULL;  /* Multiple -1 not allowed */
            unknown_idx = i;
        } else if (new_shape[i] < 0) {
            return NULL;  /* Invalid negative dimension */
        } else {
            resolved[i] = new_shape[i];
            known_product *= (size_t)new_shape[i];
        }
    }

    if (unknown_idx >= 0) {
        if (known_product == 0 || arr->size % known_product != 0) {
            return NULL;  /* Cannot determine unknown dimension */
        }
        resolved[unknown_idx] = (int32_t)(arr->size / known_product);
    }

    /* Verify total size matches */
    size_t new_size = compute_size(new_ndim, resolved);
    if (new_size != arr->size) return NULL;

    /* Only contiguous arrays can be reshaped as views */
    if (!ndarray_is_c_contiguous(arr)) return NULL;

    /* Compute new C-contiguous strides */
    int32_t new_strides[32];
    compute_strides(new_ndim, resolved, arr->dtype, new_strides);

    return ndarray_view(arr, new_ndim, resolved, new_strides);
}

EXPORT NDArray* ndarray_transpose(NDArray* arr, int32_t* axes)
{
    if (!arr) return NULL;

    int32_t new_shape[32];
    int32_t new_strides[32];
    int32_t default_axes[32];

    /* Default: reverse all axes */
    if (axes == NULL) {
        for (int i = 0; i < arr->ndim; i++) {
            default_axes[i] = arr->ndim - 1 - i;
        }
        axes = default_axes;
    }

    /* Validate axes (check for valid range and no duplicates) */
    uint32_t seen = 0;
    for (int i = 0; i < arr->ndim; i++) {
        int ax = axes[i];
        if (ax < 0) ax += arr->ndim;
        if (ax < 0 || ax >= arr->ndim) return NULL;
        if (seen & (1u << ax)) return NULL;  /* Duplicate axis */
        seen |= (1u << ax);
    }

    /* Permute shape and strides */
    for (int i = 0; i < arr->ndim; i++) {
        int ax = axes[i];
        if (ax < 0) ax += arr->ndim;
        new_shape[i] = arr->shape[ax];
        new_strides[i] = arr->strides[ax];
    }

    return ndarray_view(arr, arr->ndim, new_shape, new_strides);
}

EXPORT NDArray* ndarray_ravel(NDArray* arr)
{
    if (!arr) return NULL;

    /* If contiguous, can return view */
    if (ndarray_is_c_contiguous(arr)) {
        int32_t shape[1] = { (int32_t)arr->size };
        int32_t strides[1] = { (int32_t)dtype_size(arr->dtype) };
        return ndarray_view(arr, 1, shape, strides);
    }

    /* Non-contiguous: return NULL to signal copy needed */
    return NULL;
}

EXPORT NDArray* ndarray_flatten(NDArray* arr)
{
    if (!arr) return NULL;

    /* Always create a copy */
    int32_t shape[1] = { (int32_t)arr->size };
    NDArray* result = ndarray_create(1, shape, arr->dtype);
    if (!result) return NULL;

    /* Copy elements in row-major order */
    size_t elem_size = dtype_size(arr->dtype);

    if (ndarray_is_c_contiguous(arr)) {
        /* Fast path: memcpy */
        memcpy(result->data, arr->data, arr->size * elem_size);
    } else {
        /* Slow path: element by element using flat iteration */
        int32_t indices[32] = {0};
        for (size_t flat = 0; flat < arr->size; flat++) {
            double val = ndarray_get_item(arr, indices, arr->ndim);
            ndarray_set_flat(result, flat, val);

            /* Increment indices (row-major order) */
            for (int d = arr->ndim - 1; d >= 0; d--) {
                indices[d]++;
                if (indices[d] < arr->shape[d]) break;
                indices[d] = 0;
            }
        }
    }

    return result;
}

EXPORT NDArray* ndarray_squeeze(NDArray* arr, int32_t axis)
{
    if (!arr) return NULL;

    int32_t new_shape[32];
    int32_t new_strides[32];
    int new_ndim = 0;

    /* Use INT32_MIN as sentinel for "squeeze all" */
    if (axis == INT32_MIN) {
        /* Remove all size-1 dimensions */
        for (int i = 0; i < arr->ndim; i++) {
            if (arr->shape[i] != 1) {
                new_shape[new_ndim] = arr->shape[i];
                new_strides[new_ndim] = arr->strides[i];
                new_ndim++;
            }
        }
    } else {
        /* Remove specific axis */
        int ax = axis;
        if (ax < 0) ax += arr->ndim;
        if (ax < 0 || ax >= arr->ndim) return NULL;
        if (arr->shape[ax] != 1) return NULL;  /* Can only squeeze size-1 */

        for (int i = 0; i < arr->ndim; i++) {
            if (i != ax) {
                new_shape[new_ndim] = arr->shape[i];
                new_strides[new_ndim] = arr->strides[i];
                new_ndim++;
            }
        }
    }

    /* Handle case where all dimensions are squeezed (0-d array) */
    if (new_ndim == 0) {
        /* Create a 0-dimensional view */
        NDArray* view = (NDArray*)malloc(sizeof(NDArray));
        if (!view) return NULL;

        view->ndim = 0;
        view->dtype = arr->dtype;
        view->flags = NDARRAY_WRITEABLE;
        view->size = 1;
        view->data = arr->data;
        view->base = (arr->base != NULL) ? arr->base : arr;
        view->shape = NULL;
        view->strides = NULL;

        ndarray_update_flags(view);
        return view;
    }

    return ndarray_view(arr, new_ndim, new_shape, new_strides);
}

EXPORT NDArray* ndarray_expand_dims(NDArray* arr, int32_t axis)
{
    if (!arr) return NULL;

    /* Handle negative axis */
    int ax = axis;
    if (ax < 0) ax += arr->ndim + 1;
    if (ax < 0 || ax > arr->ndim) return NULL;

    int32_t new_shape[33];
    int32_t new_strides[33];
    int new_ndim = arr->ndim + 1;

    int j = 0;
    for (int i = 0; i < new_ndim; i++) {
        if (i == ax) {
            new_shape[i] = 1;
            /* Stride for size-1 dimension: use the next dimension's stride
             * or element size if at the end */
            if (j < arr->ndim) {
                new_strides[i] = arr->strides[j];
            } else if (arr->ndim > 0) {
                /* Use last dimension's stride * last dimension's size */
                new_strides[i] = arr->strides[arr->ndim - 1];
            } else {
                new_strides[i] = (int32_t)dtype_size(arr->dtype);
            }
        } else {
            new_shape[i] = arr->shape[j];
            new_strides[i] = arr->strides[j];
            j++;
        }
    }

    return ndarray_view(arr, new_ndim, new_shape, new_strides);
}

EXPORT NDArray* ndarray_swapaxes(NDArray* arr, int32_t axis1, int32_t axis2)
{
    if (!arr) return NULL;

    /* Handle negative axes */
    int ax1 = axis1;
    int ax2 = axis2;
    if (ax1 < 0) ax1 += arr->ndim;
    if (ax2 < 0) ax2 += arr->ndim;
    if (ax1 < 0 || ax1 >= arr->ndim) return NULL;
    if (ax2 < 0 || ax2 >= arr->ndim) return NULL;

    int32_t new_shape[32];
    int32_t new_strides[32];

    memcpy(new_shape, arr->shape, (size_t)arr->ndim * sizeof(int32_t));
    memcpy(new_strides, arr->strides, (size_t)arr->ndim * sizeof(int32_t));

    /* Swap */
    int32_t tmp = new_shape[ax1];
    new_shape[ax1] = new_shape[ax2];
    new_shape[ax2] = tmp;

    tmp = new_strides[ax1];
    new_strides[ax1] = new_strides[ax2];
    new_strides[ax2] = tmp;

    return ndarray_view(arr, arr->ndim, new_shape, new_strides);
}

/* ============ Element Access ============ */

/*
 * Compute byte offset from multi-dimensional indices using strides.
 * Returns SIZE_MAX on error.
 */
static size_t compute_byte_offset(NDArray* arr, int32_t* indices, int32_t ndim)
{
    if (!arr || !indices || ndim != arr->ndim) {
        return SIZE_MAX;
    }

    size_t byte_offset = 0;

    for (int32_t i = 0; i < ndim; i++) {
        /* Bounds check */
        if (indices[i] < 0 || indices[i] >= arr->shape[i]) {
            return SIZE_MAX;
        }
        byte_offset += (size_t)indices[i] * (size_t)arr->strides[i];
    }

    return byte_offset;
}

EXPORT size_t ndarray_flat_index(NDArray* arr, int32_t* indices, int32_t ndim)
{
    if (!arr || !indices || ndim != arr->ndim) {
        return SIZE_MAX;
    }

    /* For C-contiguous arrays, return a proper flat index */
    /* For non-contiguous, compute byte offset / elem_size (legacy behavior) */
    size_t flat = 0;
    size_t elem_size = dtype_size(arr->dtype);

    for (int32_t i = 0; i < ndim; i++) {
        /* Bounds check */
        if (indices[i] < 0 || indices[i] >= arr->shape[i]) {
            return SIZE_MAX;
        }
        flat += (size_t)indices[i] * ((size_t)arr->strides[i] / elem_size);
    }

    return flat;
}

EXPORT bool ndarray_check_bounds(NDArray* arr, int32_t* indices, int32_t ndim)
{
    if (!arr || !indices || ndim != arr->ndim) {
        return false;
    }

    for (int32_t i = 0; i < ndim; i++) {
        if (indices[i] < 0 || indices[i] >= arr->shape[i]) {
            return false;
        }
    }
    return true;
}

/*
 * Convert flat index to byte offset using strides.
 * This properly handles non-contiguous arrays (views with strides).
 */
static size_t flat_to_byte_offset(const NDArray* arr, size_t flat_idx)
{
    if (arr->ndim == 0) {
        return 0;
    }

    /* Convert flat index to multi-dimensional indices */
    size_t byte_offset = 0;
    size_t remainder = flat_idx;

    for (int d = arr->ndim - 1; d >= 0; d--) {
        int32_t dim_size = arr->shape[d];
        if (dim_size > 0) {
            size_t idx = remainder % (size_t)dim_size;
            remainder /= (size_t)dim_size;
            byte_offset += idx * (size_t)arr->strides[d];
        }
    }

    return byte_offset;
}

EXPORT double ndarray_get_flat(const NDArray* arr, size_t flat_idx)
{
    if (!arr || !arr->data || flat_idx >= arr->size) {
        return 0.0;
    }

    /* Calculate byte offset using strides */
    size_t byte_offset = flat_to_byte_offset(arr, flat_idx);
    char* ptr = (char*)arr->data + byte_offset;

    switch (arr->dtype) {
        case DTYPE_FLOAT64:
            return *(double*)ptr;
        case DTYPE_FLOAT32:
            return (double)*(float*)ptr;
        case DTYPE_INT64:
            return (double)*(int64_t*)ptr;
        case DTYPE_INT32:
            return (double)*(int32_t*)ptr;
        case DTYPE_INT16:
            return (double)*(int16_t*)ptr;
        case DTYPE_INT8:
            return (double)*(int8_t*)ptr;
        case DTYPE_UINT64:
            return (double)*(uint64_t*)ptr;
        case DTYPE_UINT32:
            return (double)*(uint32_t*)ptr;
        case DTYPE_UINT16:
            return (double)*(uint16_t*)ptr;
        case DTYPE_UINT8:
            return (double)*(uint8_t*)ptr;
        case DTYPE_BOOL:
            return (double)*(uint8_t*)ptr;
        case DTYPE_FLOAT16:
            /* Float16 not fully supported - treat as uint16 for now */
            return (double)*(uint16_t*)ptr;
        case DTYPE_COMPLEX64:
            /* Return real part */
            return (double)*(float*)ptr;
        case DTYPE_COMPLEX128:
            /* Return real part */
            return *(double*)ptr;
        default:
            return 0.0;
    }
}

EXPORT void ndarray_set_flat(NDArray* arr, size_t flat_idx, double value)
{
    if (!arr || !arr->data || flat_idx >= arr->size) {
        return;
    }
    if (!(arr->flags & NDARRAY_WRITEABLE)) {
        return;
    }

    /* Calculate byte offset using strides */
    size_t byte_offset = flat_to_byte_offset(arr, flat_idx);
    char* ptr = (char*)arr->data + byte_offset;

    switch (arr->dtype) {
        case DTYPE_FLOAT64:
            *(double*)ptr = value;
            break;
        case DTYPE_FLOAT32:
            *(float*)ptr = (float)value;
            break;
        case DTYPE_INT64:
            *(int64_t*)ptr = (int64_t)value;
            break;
        case DTYPE_INT32:
            *(int32_t*)ptr = (int32_t)value;
            break;
        case DTYPE_INT16:
            *(int16_t*)ptr = (int16_t)value;
            break;
        case DTYPE_INT8:
            *(int8_t*)ptr = (int8_t)value;
            break;
        case DTYPE_UINT64:
            *(uint64_t*)ptr = (uint64_t)value;
            break;
        case DTYPE_UINT32:
            *(uint32_t*)ptr = (uint32_t)value;
            break;
        case DTYPE_UINT16:
            *(uint16_t*)ptr = (uint16_t)value;
            break;
        case DTYPE_UINT8:
            *(uint8_t*)ptr = (uint8_t)value;
            break;
        case DTYPE_BOOL:
            *(uint8_t*)ptr = value != 0.0 ? 1 : 0;
            break;
        case DTYPE_FLOAT16:
            /* Float16 not fully supported */
            *(uint16_t*)ptr = (uint16_t)value;
            break;
        case DTYPE_COMPLEX64:
            /* Set real part, zero imag */
            *(float*)ptr = (float)value;
            *((float*)ptr + 1) = 0.0f;
            break;
        case DTYPE_COMPLEX128:
            /* Set real part, zero imag */
            *(double*)ptr = value;
            *((double*)ptr + 1) = 0.0;
            break;
        case DTYPE_COUNT:
            /* Not a valid dtype, ignore */
            break;
    }
}

/*
 * Get value at byte offset (internal helper for strided access).
 */
static double get_value_at_offset(NDArray* arr, size_t byte_offset)
{
    char* ptr = (char*)arr->data + byte_offset;

    switch (arr->dtype) {
        case DTYPE_FLOAT64:
            return *(double*)ptr;
        case DTYPE_FLOAT32:
            return (double)*(float*)ptr;
        case DTYPE_INT64:
            return (double)*(int64_t*)ptr;
        case DTYPE_INT32:
            return (double)*(int32_t*)ptr;
        case DTYPE_INT16:
            return (double)*(int16_t*)ptr;
        case DTYPE_INT8:
            return (double)*(int8_t*)ptr;
        case DTYPE_UINT64:
            return (double)*(uint64_t*)ptr;
        case DTYPE_UINT32:
            return (double)*(uint32_t*)ptr;
        case DTYPE_UINT16:
            return (double)*(uint16_t*)ptr;
        case DTYPE_UINT8:
        case DTYPE_BOOL:
            return (double)*(uint8_t*)ptr;
        case DTYPE_FLOAT16:
            return (double)*(uint16_t*)ptr;
        case DTYPE_COMPLEX64:
            return (double)*(float*)ptr;
        case DTYPE_COMPLEX128:
            return *(double*)ptr;
        default:
            return 0.0;
    }
}

/*
 * Set value at byte offset (internal helper for strided access).
 */
static void set_value_at_offset(NDArray* arr, size_t byte_offset, double value)
{
    char* ptr = (char*)arr->data + byte_offset;

    switch (arr->dtype) {
        case DTYPE_FLOAT64:
            *(double*)ptr = value;
            break;
        case DTYPE_FLOAT32:
            *(float*)ptr = (float)value;
            break;
        case DTYPE_INT64:
            *(int64_t*)ptr = (int64_t)value;
            break;
        case DTYPE_INT32:
            *(int32_t*)ptr = (int32_t)value;
            break;
        case DTYPE_INT16:
            *(int16_t*)ptr = (int16_t)value;
            break;
        case DTYPE_INT8:
            *(int8_t*)ptr = (int8_t)value;
            break;
        case DTYPE_UINT64:
            *(uint64_t*)ptr = (uint64_t)value;
            break;
        case DTYPE_UINT32:
            *(uint32_t*)ptr = (uint32_t)value;
            break;
        case DTYPE_UINT16:
            *(uint16_t*)ptr = (uint16_t)value;
            break;
        case DTYPE_UINT8:
            *(uint8_t*)ptr = (uint8_t)value;
            break;
        case DTYPE_BOOL:
            *(uint8_t*)ptr = value != 0.0 ? 1 : 0;
            break;
        case DTYPE_FLOAT16:
            *(uint16_t*)ptr = (uint16_t)value;
            break;
        case DTYPE_COMPLEX64:
            *(float*)ptr = (float)value;
            *((float*)ptr + 1) = 0.0f;
            break;
        case DTYPE_COMPLEX128:
            *(double*)ptr = value;
            *((double*)ptr + 1) = 0.0;
            break;
        default:
            break;
    }
}

EXPORT double ndarray_get_item(NDArray* arr, int32_t* indices, int32_t ndim)
{
    if (!arr || !arr->data) return 0.0;

    size_t byte_offset = compute_byte_offset(arr, indices, ndim);
    if (byte_offset == SIZE_MAX) {
        return 0.0;
    }
    return get_value_at_offset(arr, byte_offset);
}

EXPORT void ndarray_set_item(NDArray* arr, int32_t* indices, int32_t ndim, double value)
{
    if (!arr || !arr->data) return;
    if (!(arr->flags & NDARRAY_WRITEABLE)) return;

    size_t byte_offset = compute_byte_offset(arr, indices, ndim);
    if (byte_offset == SIZE_MAX) {
        return;
    }
    set_value_at_offset(arr, byte_offset, value);
}

EXPORT double ndarray_get_complex_real(NDArray* arr, size_t flat_idx)
{
    if (!arr || !arr->data || flat_idx >= arr->size) {
        return 0.0;
    }

    switch (arr->dtype) {
        case DTYPE_COMPLEX64:
            return (double)((float*)arr->data)[flat_idx * 2];
        case DTYPE_COMPLEX128:
            return ((double*)arr->data)[flat_idx * 2];
        default:
            return ndarray_get_flat(arr, flat_idx);
    }
}

EXPORT double ndarray_get_complex_imag(NDArray* arr, size_t flat_idx)
{
    if (!arr || !arr->data || flat_idx >= arr->size) {
        return 0.0;
    }

    switch (arr->dtype) {
        case DTYPE_COMPLEX64:
            return (double)((float*)arr->data)[flat_idx * 2 + 1];
        case DTYPE_COMPLEX128:
            return ((double*)arr->data)[flat_idx * 2 + 1];
        default:
            return 0.0; /* Non-complex types have zero imaginary part */
    }
}

EXPORT void ndarray_set_complex(NDArray* arr, size_t flat_idx, double real, double imag)
{
    if (!arr || !arr->data || flat_idx >= arr->size) {
        return;
    }
    if (!(arr->flags & NDARRAY_WRITEABLE)) {
        return;
    }

    switch (arr->dtype) {
        case DTYPE_COMPLEX64:
            ((float*)arr->data)[flat_idx * 2] = (float)real;
            ((float*)arr->data)[flat_idx * 2 + 1] = (float)imag;
            break;
        case DTYPE_COMPLEX128:
            ((double*)arr->data)[flat_idx * 2] = real;
            ((double*)arr->data)[flat_idx * 2 + 1] = imag;
            break;
        default:
            /* For non-complex types, just set real part */
            ndarray_set_flat(arr, flat_idx, real);
            break;
    }
}

/* ============ Array Operations ============ */

EXPORT void ndarray_fill(NDArray* arr, double value)
{
    if (!arr || !arr->data) return;

    for (size_t i = 0; i < arr->size; i++) {
        ndarray_set_flat(arr, i, value);
    }
}

EXPORT double ndarray_sum(NDArray* arr)
{
    if (!arr || !arr->data || arr->size == 0) {
        return 0.0;
    }

    /*
     * Use pairwise summation for floating point types.
     * This gives O(lg n) rounding error instead of O(n).
     * For integer types, use simple summation (exact).
     */
    switch (arr->dtype) {
        case DTYPE_FLOAT64:
            return pairwise_sum_f64((double*)arr->data, arr->size, 1);

        case DTYPE_FLOAT32:
            return (double)pairwise_sum_f32((float*)arr->data, arr->size, 1);

        case DTYPE_INT64: {
            int64_t sum = 0;
            int64_t* data = (int64_t*)arr->data;
            for (size_t i = 0; i < arr->size; i++) {
                sum += data[i];
            }
            return (double)sum;
        }

        case DTYPE_INT32: {
            int64_t sum = 0;
            int32_t* data = (int32_t*)arr->data;
            for (size_t i = 0; i < arr->size; i++) {
                sum += data[i];
            }
            return (double)sum;
        }

        case DTYPE_INT16: {
            int64_t sum = 0;
            int16_t* data = (int16_t*)arr->data;
            for (size_t i = 0; i < arr->size; i++) {
                sum += data[i];
            }
            return (double)sum;
        }

        case DTYPE_INT8: {
            int64_t sum = 0;
            int8_t* data = (int8_t*)arr->data;
            for (size_t i = 0; i < arr->size; i++) {
                sum += data[i];
            }
            return (double)sum;
        }

        case DTYPE_UINT64: {
            uint64_t sum = 0;
            uint64_t* data = (uint64_t*)arr->data;
            for (size_t i = 0; i < arr->size; i++) {
                sum += data[i];
            }
            return (double)sum;
        }

        case DTYPE_UINT32: {
            uint64_t sum = 0;
            uint32_t* data = (uint32_t*)arr->data;
            for (size_t i = 0; i < arr->size; i++) {
                sum += data[i];
            }
            return (double)sum;
        }

        case DTYPE_UINT16: {
            uint64_t sum = 0;
            uint16_t* data = (uint16_t*)arr->data;
            for (size_t i = 0; i < arr->size; i++) {
                sum += data[i];
            }
            return (double)sum;
        }

        case DTYPE_UINT8: {
            uint64_t sum = 0;
            uint8_t* data = (uint8_t*)arr->data;
            for (size_t i = 0; i < arr->size; i++) {
                sum += data[i];
            }
            return (double)sum;
        }

        case DTYPE_BOOL: {
            uint64_t sum = 0;
            uint8_t* data = (uint8_t*)arr->data;
            for (size_t i = 0; i < arr->size; i++) {
                sum += data[i] ? 1 : 0;
            }
            return (double)sum;
        }

        case DTYPE_COMPLEX64: {
            /* Sum real parts only (as NumPy does for sum()) */
            float* data = (float*)arr->data;
            double sum = 0.0;
            for (size_t i = 0; i < arr->size; i++) {
                sum += data[i * 2];
            }
            return sum;
        }

        case DTYPE_COMPLEX128: {
            double* data = (double*)arr->data;
            double sum = 0.0;
            for (size_t i = 0; i < arr->size; i++) {
                sum += data[i * 2];
            }
            return sum;
        }

        default:
            return 0.0;
    }
}

/* ============ Property Accessors ============ */

EXPORT int32_t ndarray_get_ndim(NDArray* arr)
{
    return arr ? arr->ndim : 0;
}

EXPORT int32_t* ndarray_get_shape(NDArray* arr)
{
    return arr ? arr->shape : NULL;
}

EXPORT int32_t* ndarray_get_strides(NDArray* arr)
{
    return arr ? arr->strides : NULL;
}

EXPORT void* ndarray_get_data(NDArray* arr)
{
    return arr ? arr->data : NULL;
}

EXPORT size_t ndarray_get_size(NDArray* arr)
{
    return arr ? arr->size : 0;
}

EXPORT int32_t ndarray_get_dtype(NDArray* arr)
{
    return arr ? (int32_t)arr->dtype : 0;
}

EXPORT int32_t ndarray_get_flags(NDArray* arr)
{
    return arr ? arr->flags : 0;
}

EXPORT NDArray* ndarray_get_base(NDArray* arr)
{
    return arr ? arr->base : NULL;
}

/* ============ Slicing ============ */

/*
 * Create a sliced view of an array.
 *
 * Adapted from NumPy's get_view_from_index() in mapping.c lines 865-946.
 *
 * @param arr         Source array
 * @param indices     Array of IndexSpec structs
 * @param num_indices Number of indices
 * @return            New view or NULL on error
 */
EXPORT NDArray* ndarray_slice(NDArray* arr, IndexSpec* indices, int32_t num_indices)
{
    if (!arr || !indices) return NULL;

    int32_t new_strides[32];
    int32_t new_shape[32];
    int new_dim = 0;
    int orig_dim = 0;
    size_t byte_offset = 0;

    for (int i = 0; i < num_indices; i++) {
        switch (indices[i].type) {
            case INDEX_TYPE_INTEGER: {
                /*
                 * Integer index: advances data pointer, doesn't add dimension.
                 * From NumPy mapping.c lines 879-889.
                 */
                int32_t idx = indices[i].value;

                /* Handle negative indices */
                if (idx < 0) {
                    idx += arr->shape[orig_dim];
                }

                /* Bounds check */
                if (idx < 0 || idx >= arr->shape[orig_dim]) {
                    return NULL;
                }

                byte_offset += (size_t)arr->strides[orig_dim] * (size_t)idx;
                /* new_dim unchanged - integer removes dimension */
                orig_dim++;
                break;
            }

            case INDEX_TYPE_SLICE: {
                /*
                 * Slice: advances data pointer to start, adjusts stride.
                 * From NumPy mapping.c lines 898-916.
                 */
                int32_t start = indices[i].start;
                int32_t stop = indices[i].stop;
                int32_t step = indices[i].step;
                int32_t dim_size = arr->shape[orig_dim];

                /* step defaults to 1 */
                if (step == 0) step = 1;

                /* Calculate number of steps (slice length) */
                int32_t n_steps;
                if (step > 0) {
                    n_steps = (stop > start) ? ((stop - start + step - 1) / step) : 0;
                } else {
                    n_steps = (start > stop) ? ((start - stop - step - 1) / (-step)) : 0;
                }

                if (n_steps < 0) n_steps = 0;

                /* Advance data pointer to start position */
                byte_offset += (size_t)arr->strides[orig_dim] * (size_t)start;

                /* New stride = original stride * step */
                new_strides[new_dim] = arr->strides[orig_dim] * step;
                new_shape[new_dim] = n_steps;

                new_dim++;
                orig_dim++;
                break;
            }

            case INDEX_TYPE_NEWAXIS: {
                /*
                 * Newaxis: adds dimension of size 1 with stride 0.
                 * From NumPy mapping.c lines 917-921.
                 */
                new_strides[new_dim] = 0;
                new_shape[new_dim] = 1;
                new_dim++;
                /* orig_dim unchanged - newaxis doesn't consume source dimension */
                break;
            }

            case INDEX_TYPE_ELLIPSIS: {
                /*
                 * Ellipsis: copy multiple dimensions unchanged.
                 * From NumPy mapping.c lines 890-896.
                 * indices[i].value contains number of dimensions to expand.
                 */
                int32_t expand_count = indices[i].value;
                for (int j = 0; j < expand_count; j++) {
                    new_strides[new_dim] = arr->strides[orig_dim];
                    new_shape[new_dim] = arr->shape[orig_dim];
                    new_dim++;
                    orig_dim++;
                }
                break;
            }

            default:
                return NULL;
        }
    }

    /* Handle remaining dimensions (implicit full slices) */
    while (orig_dim < arr->ndim) {
        new_strides[new_dim] = arr->strides[orig_dim];
        new_shape[new_dim] = arr->shape[orig_dim];
        new_dim++;
        orig_dim++;
    }

    /* Create view with computed offset, shape, and strides */
    return ndarray_view_with_offset(arr, new_dim, new_shape, new_strides, byte_offset);
}

/*
 * Get a sub-array by integer index along axis 0.
 * Returns a view with ndim-1 dimensions.
 *
 * @param arr   Source array
 * @param index Index along first axis (supports negative indices)
 * @return      View into arr or NULL on error
 */
EXPORT NDArray* ndarray_get_subarray(NDArray* arr, int32_t index)
{
    if (!arr) return NULL;
    if (arr->ndim == 0) return NULL;  /* Can't index scalar */

    /* Handle negative index */
    int32_t idx = index;
    if (idx < 0) {
        idx += arr->shape[0];
    }

    /* Bounds check */
    if (idx < 0 || idx >= arr->shape[0]) {
        return NULL;
    }

    /* Compute byte offset */
    size_t byte_offset = (size_t)arr->strides[0] * (size_t)idx;

    /* New shape and strides are arr's shape/strides without first element */
    int32_t new_ndim = arr->ndim - 1;

    if (new_ndim == 0) {
        /* Result is a scalar - create 0-d view */
        NDArray* view = (NDArray*)malloc(sizeof(NDArray));
        if (!view) return NULL;

        view->ndim = 0;
        view->dtype = arr->dtype;
        view->flags = NDARRAY_WRITEABLE;
        view->size = 1;
        view->data = (char*)arr->data + byte_offset;
        view->base = (arr->base != NULL) ? arr->base : arr;
        view->shape = NULL;
        view->strides = NULL;

        ndarray_update_flags(view);
        return view;
    }

    /* For ndim > 0, use rest of shape/strides */
    return ndarray_view_with_offset(arr, new_ndim, arr->shape + 1,
                                     arr->strides + 1, byte_offset);
}

/* ============ View Extensions ============ */

/*
 * Create a view with different dtype interpretation.
 * The array must be C-contiguous.
 * Last dimension is adjusted for size difference.
 *
 * @param arr   Source array
 * @param dtype New dtype
 * @return      View with new dtype or NULL on error
 */
EXPORT NDArray* ndarray_view_dtype(NDArray* arr, DType dtype)
{
    if (!arr) return NULL;

    /* Must be C-contiguous for dtype view */
    if (!ndarray_is_c_contiguous(arr)) {
        return NULL;
    }

    size_t old_itemsize = dtype_size(arr->dtype);
    size_t new_itemsize = dtype_size(dtype);

    if (old_itemsize == 0 || new_itemsize == 0) {
        return NULL;
    }

    /* Total bytes in array */
    size_t total_bytes = arr->size * old_itemsize;

    /* New total size in elements */
    if (total_bytes % new_itemsize != 0) {
        return NULL;  /* Size not divisible */
    }
    size_t new_size = total_bytes / new_itemsize;

    /* Create view struct manually (not using ndarray_view since dtype differs) */
    NDArray* view = (NDArray*)malloc(sizeof(NDArray));
    if (!view) return NULL;

    view->dtype = dtype;
    view->flags = NDARRAY_WRITEABLE;
    view->data = arr->data;
    view->base = (arr->base != NULL) ? arr->base : arr;

    if (arr->ndim == 0) {
        /* 0-d array: check if compatible */
        if (new_size != 1) {
            free(view);
            return NULL;
        }
        view->ndim = 0;
        view->size = 1;
        view->shape = NULL;
        view->strides = NULL;
    } else {
        /* Adjust last dimension */
        view->ndim = arr->ndim;
        view->shape = (int32_t*)malloc((size_t)arr->ndim * sizeof(int32_t));
        view->strides = (int32_t*)malloc((size_t)arr->ndim * sizeof(int32_t));

        if (!view->shape || !view->strides) {
            free(view->shape);
            free(view->strides);
            free(view);
            return NULL;
        }

        /* Copy all but last dimension */
        for (int i = 0; i < arr->ndim - 1; i++) {
            view->shape[i] = arr->shape[i];
        }

        /* Compute new last dimension */
        int32_t last_dim_bytes = arr->shape[arr->ndim - 1] * (int32_t)old_itemsize;
        if (last_dim_bytes % (int32_t)new_itemsize != 0) {
            free(view->shape);
            free(view->strides);
            free(view);
            return NULL;
        }
        view->shape[arr->ndim - 1] = last_dim_bytes / (int32_t)new_itemsize;

        /* Compute new strides */
        compute_strides(view->ndim, view->shape, dtype, view->strides);
        view->size = compute_size(view->ndim, view->shape);
    }

    ndarray_update_flags(view);
    return view;
}

/*
 * Return array as C-contiguous.
 * Returns view if already contiguous, otherwise copies.
 *
 * @param arr   Source array
 * @return      Contiguous array (view or copy)
 */
EXPORT NDArray* ndarray_ascontiguousarray(NDArray* arr)
{
    if (!arr) return NULL;

    /* If already C-contiguous, return a simple view */
    if (ndarray_is_c_contiguous(arr)) {
        int32_t strides[32];
        compute_strides(arr->ndim, arr->shape, arr->dtype, strides);
        return ndarray_view(arr, arr->ndim, arr->shape, strides);
    }

    /* Not contiguous: need to copy */
    return ndarray_copy(arr);
}

/*
 * Return array as Fortran-contiguous (column-major).
 * Returns view if already F-contiguous, otherwise copies.
 *
 * @param arr   Source array
 * @return      F-contiguous array (view or copy)
 */
EXPORT NDArray* ndarray_asfortranarray(NDArray* arr)
{
    if (!arr) return NULL;

    /* If already F-contiguous, return a simple view */
    if (ndarray_is_f_contiguous(arr)) {
        /* Compute F-order strides */
        int32_t strides[32];
        size_t elem_size = dtype_size(arr->dtype);

        if (arr->ndim > 0) {
            strides[0] = (int32_t)elem_size;
            for (int i = 1; i < arr->ndim; i++) {
                strides[i] = strides[i - 1] * arr->shape[i - 1];
            }
        }

        return ndarray_view(arr, arr->ndim, arr->shape, strides);
    }

    /* Not F-contiguous: need to copy with F-order layout */
    NDArray* result = ndarray_empty(arr->ndim, arr->shape, arr->dtype);
    if (!result) return NULL;

    /* Compute F-order strides for result */
    if (arr->ndim > 0) {
        size_t elem_size = dtype_size(arr->dtype);
        result->strides[0] = (int32_t)elem_size;
        for (int i = 1; i < arr->ndim; i++) {
            result->strides[i] = result->strides[i - 1] * arr->shape[i - 1];
        }
    }

    /* Copy elements in F-order */
    int32_t indices[32] = {0};
    for (size_t i = 0; i < arr->size; i++) {
        double val = ndarray_get_item(arr, indices, arr->ndim);
        ndarray_set_item(result, indices, arr->ndim, val);

        /* Increment indices (column-major order) */
        for (int d = 0; d < arr->ndim; d++) {
            indices[d]++;
            if (indices[d] < arr->shape[d]) break;
            indices[d] = 0;
        }
    }

    ndarray_update_flags(result);
    return result;
}

/* ============ WASM Memory Helpers ============ */

EXPORT void* wasm_malloc(size_t size)
{
    return malloc(size);
}

EXPORT void wasm_free(void* ptr)
{
    free(ptr);
}
