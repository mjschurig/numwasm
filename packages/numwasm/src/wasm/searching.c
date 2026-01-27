#include "searching.h"
#include <stdlib.h>
#include <math.h>

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#define EXPORT EMSCRIPTEN_KEEPALIVE
#else
#define EXPORT
#endif

static void get_reduced_shape(const int32_t* shape, int32_t ndim, int32_t axis,
                              bool keepdims, int32_t* out_shape, int32_t* out_ndim) {
    if (keepdims) {
        *out_ndim = ndim;
        for (int i = 0; i < ndim; i++)
            out_shape[i] = (i == axis) ? 1 : shape[i];
    } else {
        *out_ndim = (ndim > 1) ? ndim - 1 : 1;
        int j = 0;
        for (int i = 0; i < ndim; i++)
            if (i != axis) out_shape[j++] = shape[i];
        if (*out_ndim == 1 && j == 0) out_shape[0] = 1;
    }
}

EXPORT NDArray* ndarray_argmax(NDArray* arr, int32_t axis, bool keepdims) {
    if (!arr || !arr->data) return NULL;

    /* Flatten case */
    if (axis == -2147483648) {
        size_t max_idx = 0;
        double max_val = ndarray_get_flat(arr, 0);
        for (size_t i = 1; i < arr->size; i++) {
            double val = ndarray_get_flat(arr, i);
            if (!isnan(val) && (isnan(max_val) || val > max_val)) {
                max_val = val;
                max_idx = i;
            }
        }
        int32_t shape_one[1] = {1};
        NDArray* result;
        if (keepdims) {
            int32_t* kd_shape = malloc(arr->ndim * sizeof(int32_t));
            if (!kd_shape) return NULL;
            for (int i = 0; i < arr->ndim; i++) kd_shape[i] = 1;
            result = ndarray_create(arr->ndim, kd_shape, DTYPE_INT32);
            free(kd_shape);
        } else {
            result = ndarray_create(1, shape_one, DTYPE_INT32);
        }
        if (!result) return NULL;
        ((int32_t*)result->data)[0] = (int32_t)max_idx;
        return result;
    }

    if (axis < 0) axis += arr->ndim;
    if (axis < 0 || axis >= arr->ndim) return NULL;

    int32_t out_shape[32], out_ndim;
    get_reduced_shape(arr->shape, arr->ndim, axis, keepdims, out_shape, &out_ndim);

    NDArray* result = ndarray_create(out_ndim, out_shape, DTYPE_INT32);
    if (!result) return NULL;

    int32_t axis_size = arr->shape[axis];
    size_t outer_size = 1, inner_size = 1;
    for (int i = 0; i < axis; i++) outer_size *= arr->shape[i];
    for (int i = axis + 1; i < arr->ndim; i++) inner_size *= arr->shape[i];

    int32_t* res_data = (int32_t*)result->data;
    size_t res_idx = 0;

    for (size_t outer = 0; outer < outer_size; outer++) {
        for (size_t inner = 0; inner < inner_size; inner++) {
            int32_t max_idx = 0;
            double max_val = ndarray_get_flat(arr,
                outer * axis_size * inner_size + inner);

            for (int32_t i = 1; i < axis_size; i++) {
                size_t idx = outer * axis_size * inner_size + i * inner_size + inner;
                double val = ndarray_get_flat(arr, idx);
                if (!isnan(val) && (isnan(max_val) || val > max_val)) {
                    max_val = val;
                    max_idx = i;
                }
            }
            res_data[res_idx++] = max_idx;
        }
    }
    return result;
}

EXPORT NDArray* ndarray_argmin(NDArray* arr, int32_t axis, bool keepdims) {
    if (!arr || !arr->data) return NULL;

    if (axis == -2147483648) {
        size_t min_idx = 0;
        double min_val = ndarray_get_flat(arr, 0);
        for (size_t i = 1; i < arr->size; i++) {
            double val = ndarray_get_flat(arr, i);
            if (!isnan(val) && (isnan(min_val) || val < min_val)) {
                min_val = val;
                min_idx = i;
            }
        }
        int32_t shape_one[1] = {1};
        NDArray* result;
        if (keepdims) {
            int32_t* kd_shape = malloc(arr->ndim * sizeof(int32_t));
            if (!kd_shape) return NULL;
            for (int i = 0; i < arr->ndim; i++) kd_shape[i] = 1;
            result = ndarray_create(arr->ndim, kd_shape, DTYPE_INT32);
            free(kd_shape);
        } else {
            result = ndarray_create(1, shape_one, DTYPE_INT32);
        }
        if (!result) return NULL;
        ((int32_t*)result->data)[0] = (int32_t)min_idx;
        return result;
    }

    if (axis < 0) axis += arr->ndim;
    if (axis < 0 || axis >= arr->ndim) return NULL;

    int32_t out_shape[32], out_ndim;
    get_reduced_shape(arr->shape, arr->ndim, axis, keepdims, out_shape, &out_ndim);

    NDArray* result = ndarray_create(out_ndim, out_shape, DTYPE_INT32);
    if (!result) return NULL;

    int32_t axis_size = arr->shape[axis];
    size_t outer_size = 1, inner_size = 1;
    for (int i = 0; i < axis; i++) outer_size *= arr->shape[i];
    for (int i = axis + 1; i < arr->ndim; i++) inner_size *= arr->shape[i];

    int32_t* res_data = (int32_t*)result->data;
    size_t res_idx = 0;

    for (size_t outer = 0; outer < outer_size; outer++) {
        for (size_t inner = 0; inner < inner_size; inner++) {
            int32_t min_idx = 0;
            double min_val = ndarray_get_flat(arr,
                outer * axis_size * inner_size + inner);

            for (int32_t i = 1; i < axis_size; i++) {
                size_t idx = outer * axis_size * inner_size + i * inner_size + inner;
                double val = ndarray_get_flat(arr, idx);
                if (!isnan(val) && (isnan(min_val) || val < min_val)) {
                    min_val = val;
                    min_idx = i;
                }
            }
            res_data[res_idx++] = min_idx;
        }
    }
    return result;
}

EXPORT NDArray* ndarray_searchsorted(NDArray* sorted_arr, NDArray* values,
                                      int32_t side, NDArray* sorter) {
    if (!sorted_arr || !values || sorted_arr->ndim != 1) return NULL;

    NDArray* result = ndarray_empty(values->ndim, values->shape, DTYPE_INT32);
    if (!result) return NULL;

    int32_t* res_data = (int32_t*)result->data;

    for (size_t i = 0; i < values->size; i++) {
        double val = ndarray_get_flat(values, i);
        size_t lo = 0, hi = sorted_arr->size;

        while (lo < hi) {
            size_t mid = lo + (hi - lo) / 2;
            size_t idx = sorter ? (size_t)ndarray_get_flat(sorter, mid) : mid;
            double arr_val = ndarray_get_flat(sorted_arr, idx);

            if (side == SEARCH_RIGHT ? arr_val <= val : arr_val < val)
                lo = mid + 1;
            else
                hi = mid;
        }
        res_data[i] = (int32_t)lo;
    }
    return result;
}
