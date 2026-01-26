/**
 * NumJS-WASM Index Functions Implementation
 *
 * Adapted from NumPy's item_selection.c and multiarraymodule.c.
 */

#include "indexing.h"
#include "broadcast.h"
#include <stdlib.h>
#include <string.h>

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#define EXPORT EMSCRIPTEN_KEEPALIVE
#else
#define EXPORT
#endif

/* ============ Helper Functions ============ */

/**
 * Handle index clipping/wrapping based on mode.
 * Returns -1 if index is invalid in RAISE mode.
 */
static int32_t handle_index(int32_t idx, int32_t max_idx, int32_t clipmode)
{
    if (idx < 0) {
        /* Try wrapping negative index */
        idx += max_idx;
    }

    if (idx < 0 || idx >= max_idx) {
        switch (clipmode) {
            case CLIP_WRAP:
                /* Wrap around using modulo */
                if (max_idx == 0) return -1;
                idx = ((idx % max_idx) + max_idx) % max_idx;
                break;
            case CLIP_CLIP:
                /* Clip to valid range */
                if (idx < 0) idx = 0;
                if (idx >= max_idx) idx = max_idx - 1;
                break;
            case CLIP_RAISE:
            default:
                return -1;  /* Signal error */
        }
    }

    return idx;
}

/**
 * Get index value from array at flat position.
 */
static int32_t get_index_value(NDArray* indices, size_t flat_idx)
{
    switch (indices->dtype) {
        case DTYPE_INT32:
            return ((int32_t*)indices->data)[flat_idx];
        case DTYPE_INT64:
            return (int32_t)((int64_t*)indices->data)[flat_idx];
        case DTYPE_INT16:
            return (int32_t)((int16_t*)indices->data)[flat_idx];
        case DTYPE_INT8:
            return (int32_t)((int8_t*)indices->data)[flat_idx];
        case DTYPE_UINT32:
            return (int32_t)((uint32_t*)indices->data)[flat_idx];
        case DTYPE_UINT64:
            return (int32_t)((uint64_t*)indices->data)[flat_idx];
        case DTYPE_UINT16:
            return (int32_t)((uint16_t*)indices->data)[flat_idx];
        case DTYPE_UINT8:
            return (int32_t)((uint8_t*)indices->data)[flat_idx];
        case DTYPE_FLOAT64:
            return (int32_t)((double*)indices->data)[flat_idx];
        case DTYPE_FLOAT32:
            return (int32_t)((float*)indices->data)[flat_idx];
        default:
            return 0;
    }
}

/**
 * Check if value is nonzero.
 */
static bool is_nonzero(NDArray* arr, size_t flat_idx)
{
    switch (arr->dtype) {
        case DTYPE_FLOAT64:
            return ((double*)arr->data)[flat_idx] != 0.0;
        case DTYPE_FLOAT32:
            return ((float*)arr->data)[flat_idx] != 0.0f;
        case DTYPE_INT64:
            return ((int64_t*)arr->data)[flat_idx] != 0;
        case DTYPE_INT32:
            return ((int32_t*)arr->data)[flat_idx] != 0;
        case DTYPE_INT16:
            return ((int16_t*)arr->data)[flat_idx] != 0;
        case DTYPE_INT8:
            return ((int8_t*)arr->data)[flat_idx] != 0;
        case DTYPE_UINT64:
            return ((uint64_t*)arr->data)[flat_idx] != 0;
        case DTYPE_UINT32:
            return ((uint32_t*)arr->data)[flat_idx] != 0;
        case DTYPE_UINT16:
            return ((uint16_t*)arr->data)[flat_idx] != 0;
        case DTYPE_UINT8:
        case DTYPE_BOOL:
            return ((uint8_t*)arr->data)[flat_idx] != 0;
        case DTYPE_COMPLEX64: {
            float* data = (float*)arr->data;
            return data[flat_idx * 2] != 0.0f || data[flat_idx * 2 + 1] != 0.0f;
        }
        case DTYPE_COMPLEX128: {
            double* data = (double*)arr->data;
            return data[flat_idx * 2] != 0.0 || data[flat_idx * 2 + 1] != 0.0;
        }
        default:
            return false;
    }
}

/* ============ Take ============ */

/**
 * Take elements from flattened array.
 *
 * Adapted from NumPy's take operation.
 */
EXPORT NDArray* ndarray_take_flat(NDArray* arr, NDArray* indices, int32_t clipmode)
{
    if (!arr || !indices) return NULL;

    /* Result has same shape as indices */
    NDArray* result = ndarray_empty(indices->ndim, indices->shape, arr->dtype);
    if (!result) return NULL;

    int32_t max_idx = (int32_t)arr->size;

    for (size_t i = 0; i < indices->size; i++) {
        int32_t idx = get_index_value(indices, i);
        idx = handle_index(idx, max_idx, clipmode);

        if (idx < 0) {
            /* Out of bounds in RAISE mode */
            ndarray_free(result);
            return NULL;
        }

        double val = ndarray_get_flat(arr, (size_t)idx);
        ndarray_set_flat(result, i, val);
    }

    return result;
}

/**
 * Take elements along an axis.
 *
 * Adapted from NumPy's PyArray_TakeFrom() in item_selection.c.
 */
EXPORT NDArray* ndarray_take(NDArray* arr, NDArray* indices, int32_t axis, int32_t clipmode)
{
    if (!arr || !indices) return NULL;

    /* Handle axis */
    if (axis < 0) axis += arr->ndim;
    if (axis < 0 || axis >= arr->ndim) return NULL;

    /* For now, only support axis=0 for simplicity, or use flat version */
    if (arr->ndim == 1 || axis == 0) {
        /* Compute result shape:
         * For 1D: same as indices shape
         * For nD with axis=0: replace first dim with indices shape
         */
        int32_t result_ndim = arr->ndim - 1 + indices->ndim;
        int32_t result_shape[32];

        /* Copy indices shape */
        for (int i = 0; i < indices->ndim; i++) {
            result_shape[i] = indices->shape[i];
        }
        /* Copy remaining arr shape (after axis) */
        for (int i = 1; i < arr->ndim; i++) {
            result_shape[indices->ndim + i - 1] = arr->shape[i];
        }

        NDArray* result = ndarray_empty(result_ndim, result_shape, arr->dtype);
        if (!result) return NULL;

        /* Size of a single "chunk" (all dims after axis 0) */
        size_t chunk_size = 1;
        for (int i = 1; i < arr->ndim; i++) {
            chunk_size *= (size_t)arr->shape[i];
        }

        int32_t max_idx = arr->shape[0];
        size_t elem_size = dtype_size(arr->dtype);

        for (size_t i = 0; i < indices->size; i++) {
            int32_t idx = get_index_value(indices, i);
            idx = handle_index(idx, max_idx, clipmode);

            if (idx < 0) {
                ndarray_free(result);
                return NULL;
            }

            /* Copy chunk from arr[idx] to result[i] */
            char* src = (char*)arr->data + (size_t)idx * chunk_size * elem_size;
            char* dst = (char*)result->data + i * chunk_size * elem_size;
            memcpy(dst, src, chunk_size * elem_size);
        }

        return result;
    }

    /* For other axes, we need more complex logic */
    /* For now, return NULL (not implemented) */
    return NULL;
}

/* ============ Put ============ */

/**
 * Put values into an array at specified flat indices.
 *
 * Adapted from NumPy's PyArray_PutTo() in item_selection.c.
 */
EXPORT int ndarray_put(NDArray* arr, NDArray* indices, NDArray* values, int32_t clipmode)
{
    if (!arr || !indices || !values) return -1;
    if (!(arr->flags & NDARRAY_WRITEABLE)) return -1;

    int32_t max_idx = (int32_t)arr->size;
    size_t nv = values->size;

    if (nv == 0) return 0;  /* Nothing to do */

    for (size_t i = 0; i < indices->size; i++) {
        int32_t idx = get_index_value(indices, i);
        idx = handle_index(idx, max_idx, clipmode);

        if (idx < 0) {
            /* Out of bounds in RAISE mode */
            return -1;
        }

        /* Broadcast values if smaller than indices */
        double val = ndarray_get_flat(values, i % nv);
        ndarray_set_flat(arr, (size_t)idx, val);
    }

    return 0;
}

/* ============ Nonzero ============ */

/**
 * Count nonzero elements.
 */
EXPORT size_t ndarray_count_nonzero(NDArray* arr)
{
    if (!arr || !arr->data) return 0;

    size_t count = 0;
    for (size_t i = 0; i < arr->size; i++) {
        if (is_nonzero(arr, i)) {
            count++;
        }
    }
    return count;
}

/**
 * Find indices of nonzero elements.
 *
 * Returns 2D array of shape (num_nonzero, ndim).
 * Adapted from NumPy's PyArray_Nonzero().
 */
EXPORT NDArray* ndarray_nonzero(NDArray* arr)
{
    if (!arr || !arr->data) return NULL;
    if (arr->ndim == 0) return NULL;  /* Can't call nonzero on 0-d array */

    /* Count nonzero elements */
    size_t nonzero_count = ndarray_count_nonzero(arr);

    /* Create result array: shape (nonzero_count, ndim) */
    int32_t result_shape[2] = { (int32_t)nonzero_count, arr->ndim };
    NDArray* result = ndarray_create(2, result_shape, DTYPE_INT32);
    if (!result) return NULL;

    if (nonzero_count == 0) {
        return result;  /* Empty result */
    }

    /* Fill in indices */
    int32_t* result_data = (int32_t*)result->data;
    size_t result_row = 0;

    /* Iterate in C-order (row-major) */
    int32_t indices[32] = {0};

    for (size_t flat = 0; flat < arr->size; flat++) {
        if (is_nonzero(arr, flat)) {
            /* Store current multi-index */
            for (int d = 0; d < arr->ndim; d++) {
                result_data[result_row * arr->ndim + d] = indices[d];
            }
            result_row++;
        }

        /* Increment indices (C-order) */
        for (int d = arr->ndim - 1; d >= 0; d--) {
            indices[d]++;
            if (indices[d] < arr->shape[d]) break;
            indices[d] = 0;
        }
    }

    return result;
}

/**
 * Return flat indices of nonzero elements.
 */
EXPORT NDArray* ndarray_flatnonzero(NDArray* arr)
{
    if (!arr || !arr->data) return NULL;

    /* Count nonzero elements */
    size_t nonzero_count = ndarray_count_nonzero(arr);

    /* Create result array */
    int32_t result_shape[1] = { (int32_t)nonzero_count };
    NDArray* result = ndarray_create(1, result_shape, DTYPE_INT32);
    if (!result) return NULL;

    int32_t* result_data = (int32_t*)result->data;
    size_t result_idx = 0;

    for (size_t i = 0; i < arr->size; i++) {
        if (is_nonzero(arr, i)) {
            result_data[result_idx++] = (int32_t)i;
        }
    }

    return result;
}

/* ============ Where ============ */

/**
 * Return elements chosen from x or y depending on condition.
 *
 * Adapted from NumPy's PyArray_Where().
 */
EXPORT NDArray* ndarray_where(NDArray* condition, NDArray* x, NDArray* y)
{
    if (!condition) return NULL;

    /* If x and y are NULL, return nonzero indices */
    if (x == NULL && y == NULL) {
        return ndarray_nonzero(condition);
    }

    /* Both x and y must be provided */
    if (x == NULL || y == NULL) return NULL;

    /* Broadcast shapes */
    int32_t bc_shape[NPY_MAXDIMS];
    int32_t bc_ndim;

    const int32_t* shapes[3] = { condition->shape, x->shape, y->shape };
    int32_t ndims[3] = { condition->ndim, x->ndim, y->ndim };

    if (broadcast_shapes_multi(shapes, ndims, 3, bc_shape, &bc_ndim) != 0) {
        return NULL;  /* Shapes not broadcastable */
    }

    /* Create result array */
    NDArray* result = ndarray_empty(bc_ndim, bc_shape, x->dtype);  /* Use x's dtype */
    if (!result) return NULL;

    /* Compute broadcast strides for each input */
    int32_t cond_strides[NPY_MAXDIMS];
    int32_t x_strides[NPY_MAXDIMS];
    int32_t y_strides[NPY_MAXDIMS];

    if (broadcast_strides(condition, bc_shape, bc_ndim, cond_strides) != 0 ||
        broadcast_strides(x, bc_shape, bc_ndim, x_strides) != 0 ||
        broadcast_strides(y, bc_shape, bc_ndim, y_strides) != 0) {
        ndarray_free(result);
        return NULL;
    }

    /* Iterate over result and select from x or y based on condition */
    size_t elem_size_cond = dtype_size(condition->dtype);
    size_t elem_size_x = dtype_size(x->dtype);
    size_t elem_size_y = dtype_size(y->dtype);

    int32_t indices[32] = {0};

    for (size_t i = 0; i < result->size; i++) {
        /* Compute offsets into each array */
        size_t cond_offset = 0;
        size_t x_offset = 0;
        size_t y_offset = 0;

        for (int d = 0; d < bc_ndim; d++) {
            cond_offset += (size_t)indices[d] * ((size_t)cond_strides[d] / elem_size_cond);
            x_offset += (size_t)indices[d] * ((size_t)x_strides[d] / elem_size_x);
            y_offset += (size_t)indices[d] * ((size_t)y_strides[d] / elem_size_y);
        }

        /* Check condition */
        bool cond_val = is_nonzero(condition, cond_offset);

        /* Select from x or y */
        double val;
        if (cond_val) {
            val = ndarray_get_flat(x, x_offset);
        } else {
            val = ndarray_get_flat(y, y_offset);
        }
        ndarray_set_flat(result, i, val);

        /* Increment indices (C-order) */
        for (int d = bc_ndim - 1; d >= 0; d--) {
            indices[d]++;
            if (indices[d] < bc_shape[d]) break;
            indices[d] = 0;
        }
    }

    return result;
}

/* ============ Compress ============ */

/**
 * Return selected slices of an array along given axis.
 */
EXPORT NDArray* ndarray_compress(NDArray* condition, NDArray* arr, int32_t axis)
{
    if (!condition || !arr) return NULL;
    if (condition->ndim != 1) return NULL;  /* Condition must be 1D */

    /* Handle axis */
    if (axis < 0) axis += arr->ndim;
    if (axis < 0 || axis >= arr->ndim) return NULL;

    /* Count true values in condition */
    size_t true_count = 0;
    for (size_t i = 0; i < condition->size; i++) {
        if (is_nonzero(condition, i)) true_count++;
    }

    /* Build result shape */
    int32_t result_shape[32];
    for (int i = 0; i < arr->ndim; i++) {
        if (i == axis) {
            result_shape[i] = (int32_t)true_count;
        } else {
            result_shape[i] = arr->shape[i];
        }
    }

    NDArray* result = ndarray_empty(arr->ndim, result_shape, arr->dtype);
    if (!result) return NULL;

    /* For simplicity, only implement axis=0 */
    if (axis == 0) {
        size_t chunk_size = 1;
        for (int i = 1; i < arr->ndim; i++) {
            chunk_size *= (size_t)arr->shape[i];
        }

        size_t elem_size = dtype_size(arr->dtype);
        size_t result_idx = 0;

        for (size_t i = 0; i < condition->size && i < (size_t)arr->shape[0]; i++) {
            if (is_nonzero(condition, i)) {
                char* src = (char*)arr->data + i * chunk_size * elem_size;
                char* dst = (char*)result->data + result_idx * chunk_size * elem_size;
                memcpy(dst, src, chunk_size * elem_size);
                result_idx++;
            }
        }
    } else {
        /* For other axes, need more complex implementation */
        ndarray_free(result);
        return NULL;
    }

    return result;
}

/* ============ Extract ============ */

/**
 * Return elements of an array that satisfy some condition.
 */
EXPORT NDArray* ndarray_extract(NDArray* condition, NDArray* arr)
{
    if (!condition || !arr) return NULL;

    /* Condition and arr must have same size */
    if (condition->size != arr->size) return NULL;

    /* Count true values */
    size_t true_count = ndarray_count_nonzero(condition);

    /* Create 1D result */
    int32_t result_shape[1] = { (int32_t)true_count };
    NDArray* result = ndarray_empty(1, result_shape, arr->dtype);
    if (!result) return NULL;

    size_t result_idx = 0;
    for (size_t i = 0; i < arr->size; i++) {
        if (is_nonzero(condition, i)) {
            double val = ndarray_get_flat(arr, i);
            ndarray_set_flat(result, result_idx++, val);
        }
    }

    return result;
}

/* ============ Choose ============ */

/**
 * Construct an array from an index array and a set of arrays to choose from.
 */
EXPORT NDArray* ndarray_choose(NDArray* indices, NDArray** choices, int32_t num_choices,
                               int32_t clipmode)
{
    if (!indices || !choices || num_choices <= 0) return NULL;

    /* Verify all choices have compatible shapes */
    for (int i = 0; i < num_choices; i++) {
        if (!choices[i]) return NULL;
        /* For simplicity, require same shape as indices */
        if (choices[i]->size != indices->size) return NULL;
    }

    /* Create result with same shape as indices, dtype from first choice */
    NDArray* result = ndarray_empty(indices->ndim, indices->shape, choices[0]->dtype);
    if (!result) return NULL;

    for (size_t i = 0; i < indices->size; i++) {
        int32_t choice_idx = get_index_value(indices, i);
        choice_idx = handle_index(choice_idx, num_choices, clipmode);

        if (choice_idx < 0) {
            ndarray_free(result);
            return NULL;
        }

        double val = ndarray_get_flat(choices[choice_idx], i);
        ndarray_set_flat(result, i, val);
    }

    return result;
}

/* ============ Diagonal ============ */

/**
 * Return specified diagonals.
 */
EXPORT NDArray* ndarray_diagonal(NDArray* arr, int32_t offset, int32_t axis1, int32_t axis2)
{
    if (!arr || arr->ndim < 2) return NULL;

    /* Handle negative axes */
    if (axis1 < 0) axis1 += arr->ndim;
    if (axis2 < 0) axis2 += arr->ndim;
    if (axis1 < 0 || axis1 >= arr->ndim) return NULL;
    if (axis2 < 0 || axis2 >= arr->ndim) return NULL;
    if (axis1 == axis2) return NULL;

    /* For simplicity, only support 2D arrays with axis1=0, axis2=1 */
    if (arr->ndim != 2 || axis1 != 0 || axis2 != 1) {
        return NULL;  /* Not implemented */
    }

    int32_t n_rows = arr->shape[0];
    int32_t n_cols = arr->shape[1];

    /* Calculate diagonal length */
    int32_t diag_len;
    int32_t start_row, start_col;

    if (offset >= 0) {
        start_row = 0;
        start_col = offset;
        diag_len = (n_rows < n_cols - offset) ? n_rows : n_cols - offset;
    } else {
        start_row = -offset;
        start_col = 0;
        diag_len = (n_rows + offset < n_cols) ? n_rows + offset : n_cols;
    }

    if (diag_len <= 0) {
        /* Empty diagonal */
        int32_t shape[1] = { 0 };
        return ndarray_create(1, shape, arr->dtype);
    }

    /* Create result array */
    int32_t result_shape[1] = { diag_len };
    NDArray* result = ndarray_empty(1, result_shape, arr->dtype);
    if (!result) return NULL;

    /* Extract diagonal elements */
    for (int32_t i = 0; i < diag_len; i++) {
        int32_t indices[2] = { start_row + i, start_col + i };
        double val = ndarray_get_item(arr, indices, 2);
        ndarray_set_flat(result, (size_t)i, val);
    }

    return result;
}
