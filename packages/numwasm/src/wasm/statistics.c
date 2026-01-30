#include "statistics.h"
#include "sorting.h"
#include "pairwise_sum.h"
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
    if (axis == -2147483648) {
        if (keepdims) {
            *out_ndim = ndim;
            for (int i = 0; i < ndim; i++) out_shape[i] = 1;
        } else {
            *out_ndim = 0;
        }
        return;
    }
    if (keepdims) {
        *out_ndim = ndim;
        for (int i = 0; i < ndim; i++)
            out_shape[i] = (i == axis) ? 1 : shape[i];
    } else {
        *out_ndim = (ndim > 1) ? ndim - 1 : 0;
        int j = 0;
        for (int i = 0; i < ndim; i++)
            if (i != axis) out_shape[j++] = shape[i];
    }
}

static DType get_reduction_dtype(DType input, DType requested) {
    if (requested != (DType)-1) return requested;
    switch (input) {
        case DTYPE_INT8: case DTYPE_INT16: case DTYPE_INT32: case DTYPE_INT64:
        case DTYPE_UINT8: case DTYPE_UINT16: case DTYPE_UINT32: case DTYPE_UINT64:
        case DTYPE_BOOL:
            return DTYPE_FLOAT64;
        default:
            return input;
    }
}

EXPORT NDArray* ndarray_sum_axis(NDArray* arr, int32_t axis, bool keepdims, DType dtype) {
    if (!arr || !arr->data) return NULL;
    DType out_dtype = get_reduction_dtype(arr->dtype, dtype);

    if (axis == -2147483648) {
        double total = ndarray_sum(arr);  /* Uses pairwise sum */
        if (keepdims) {
            int32_t* shape = malloc(arr->ndim * sizeof(int32_t));
            if (!shape) return NULL;
            for (int i = 0; i < arr->ndim; i++) shape[i] = 1;
            NDArray* r = ndarray_full(arr->ndim, shape, out_dtype, total);
            free(shape);
            return r;
        }
        return ndarray_scalar(total, out_dtype);
    }

    if (axis < 0) axis += arr->ndim;
    if (axis < 0 || axis >= arr->ndim) return NULL;

    int32_t out_shape[32], out_ndim;
    get_reduced_shape(arr->shape, arr->ndim, axis, keepdims, out_shape, &out_ndim);

    NDArray* result = out_ndim == 0 ? ndarray_scalar(0, out_dtype) :
                      ndarray_create(out_ndim, out_shape, out_dtype);
    if (!result) return NULL;

    int32_t axis_size = arr->shape[axis];
    size_t outer = 1, inner = 1;
    for (int i = 0; i < axis; i++) outer *= arr->shape[i];
    for (int i = axis + 1; i < arr->ndim; i++) inner *= arr->shape[i];

    double* temp = malloc(axis_size * sizeof(double));
    if (!temp) { ndarray_free(result); return NULL; }
    size_t res_idx = 0;

    for (size_t o = 0; o < outer; o++) {
        for (size_t n = 0; n < inner; n++) {
            for (int32_t i = 0; i < axis_size; i++) {
                size_t idx = o * axis_size * inner + i * inner + n;
                temp[i] = ndarray_get_flat(arr, idx);
            }
            double sum = pairwise_sum_f64(temp, axis_size, 1);
            ndarray_set_flat(result, res_idx++, sum);
        }
    }
    free(temp);
    return result;
}

EXPORT NDArray* ndarray_mean_axis(NDArray* arr, int32_t axis, bool keepdims, DType dtype) {
    if (!arr) return NULL;

    if (axis == -2147483648) {
        double mean_val = ndarray_sum(arr) / (double)arr->size;
        DType out_dtype = get_reduction_dtype(arr->dtype, dtype);
        if (keepdims) {
            int32_t* shape = malloc(arr->ndim * sizeof(int32_t));
            if (!shape) return NULL;
            for (int i = 0; i < arr->ndim; i++) shape[i] = 1;
            NDArray* r = ndarray_full(arr->ndim, shape, out_dtype, mean_val);
            free(shape);
            return r;
        }
        return ndarray_scalar(mean_val, out_dtype);
    }

    NDArray* sum_res = ndarray_sum_axis(arr, axis, keepdims, dtype);
    if (!sum_res) return NULL;

    int32_t norm_axis = axis < 0 ? axis + arr->ndim : axis;
    int32_t count = arr->shape[norm_axis];

    for (size_t i = 0; i < sum_res->size; i++) {
        double val = ndarray_get_flat(sum_res, i);
        ndarray_set_flat(sum_res, i, val / count);
    }
    return sum_res;
}

EXPORT NDArray* ndarray_var_axis(NDArray* arr, int32_t axis, bool keepdims,
                                  DType dtype, int32_t ddof) {
    if (!arr) return NULL;
    DType out_dtype = get_reduction_dtype(arr->dtype, dtype);

    if (axis == -2147483648) {
        double mean_val = ndarray_sum(arr) / (double)arr->size;
        double sum_sq = 0;
        for (size_t i = 0; i < arr->size; i++) {
            double diff = ndarray_get_flat(arr, i) - mean_val;
            sum_sq += diff * diff;
        }
        double var_val = sum_sq / (arr->size - ddof);
        if (keepdims) {
            int32_t* shape = malloc(arr->ndim * sizeof(int32_t));
            if (!shape) return NULL;
            for (int i = 0; i < arr->ndim; i++) shape[i] = 1;
            NDArray* r = ndarray_full(arr->ndim, shape, out_dtype, var_val);
            free(shape);
            return r;
        }
        return ndarray_scalar(var_val, out_dtype);
    }

    if (axis < 0) axis += arr->ndim;
    if (axis < 0 || axis >= arr->ndim) return NULL;

    NDArray* mean_res = ndarray_mean_axis(arr, axis, true, dtype);
    if (!mean_res) return NULL;

    int32_t out_shape[32], out_ndim;
    get_reduced_shape(arr->shape, arr->ndim, axis, keepdims, out_shape, &out_ndim);

    NDArray* result = out_ndim == 0 ? ndarray_scalar(0, out_dtype) :
                      ndarray_create(out_ndim, out_shape, out_dtype);
    if (!result) { ndarray_free(mean_res); return NULL; }

    int32_t axis_size = arr->shape[axis];
    size_t outer = 1, inner = 1;
    for (int i = 0; i < axis; i++) outer *= arr->shape[i];
    for (int i = axis + 1; i < arr->ndim; i++) inner *= arr->shape[i];

    size_t res_idx = 0, mean_idx = 0;
    for (size_t o = 0; o < outer; o++) {
        for (size_t n = 0; n < inner; n++) {
            double mean_val = ndarray_get_flat(mean_res, mean_idx++);
            double sum_sq = 0;
            for (int32_t i = 0; i < axis_size; i++) {
                size_t idx = o * axis_size * inner + i * inner + n;
                double diff = ndarray_get_flat(arr, idx) - mean_val;
                sum_sq += diff * diff;
            }
            ndarray_set_flat(result, res_idx++, sum_sq / (axis_size - ddof));
        }
    }
    ndarray_free(mean_res);
    return result;
}

EXPORT NDArray* ndarray_std_axis(NDArray* arr, int32_t axis, bool keepdims,
                                  DType dtype, int32_t ddof) {
    NDArray* var_res = ndarray_var_axis(arr, axis, keepdims, dtype, ddof);
    if (!var_res) return NULL;
    for (size_t i = 0; i < var_res->size; i++) {
        ndarray_set_flat(var_res, i, sqrt(ndarray_get_flat(var_res, i)));
    }
    return var_res;
}

EXPORT NDArray* ndarray_min_axis(NDArray* arr, int32_t axis, bool keepdims) {
    if (!arr) return NULL;

    if (axis == -2147483648) {
        double min_val = ndarray_get_flat(arr, 0);
        for (size_t i = 1; i < arr->size; i++) {
            double v = ndarray_get_flat(arr, i);
            if (v < min_val) min_val = v;
        }
        if (keepdims) {
            int32_t* shape = malloc(arr->ndim * sizeof(int32_t));
            if (!shape) return NULL;
            for (int i = 0; i < arr->ndim; i++) shape[i] = 1;
            NDArray* r = ndarray_full(arr->ndim, shape, arr->dtype, min_val);
            free(shape);
            return r;
        }
        return ndarray_scalar(min_val, arr->dtype);
    }

    if (axis < 0) axis += arr->ndim;
    if (axis < 0 || axis >= arr->ndim) return NULL;

    int32_t out_shape[32], out_ndim;
    get_reduced_shape(arr->shape, arr->ndim, axis, keepdims, out_shape, &out_ndim);

    NDArray* result = out_ndim == 0 ? ndarray_scalar(0, arr->dtype) :
                      ndarray_create(out_ndim, out_shape, arr->dtype);
    if (!result) return NULL;

    int32_t axis_size = arr->shape[axis];
    size_t outer = 1, inner = 1;
    for (int i = 0; i < axis; i++) outer *= arr->shape[i];
    for (int i = axis + 1; i < arr->ndim; i++) inner *= arr->shape[i];

    size_t res_idx = 0;
    for (size_t o = 0; o < outer; o++) {
        for (size_t n = 0; n < inner; n++) {
            double min_val = ndarray_get_flat(arr, o * axis_size * inner + n);
            for (int32_t i = 1; i < axis_size; i++) {
                double v = ndarray_get_flat(arr, o * axis_size * inner + i * inner + n);
                if (v < min_val) min_val = v;
            }
            ndarray_set_flat(result, res_idx++, min_val);
        }
    }
    return result;
}

EXPORT NDArray* ndarray_max_axis(NDArray* arr, int32_t axis, bool keepdims) {
    if (!arr) return NULL;

    if (axis == -2147483648) {
        double max_val = ndarray_get_flat(arr, 0);
        for (size_t i = 1; i < arr->size; i++) {
            double v = ndarray_get_flat(arr, i);
            if (v > max_val) max_val = v;
        }
        if (keepdims) {
            int32_t* shape = malloc(arr->ndim * sizeof(int32_t));
            if (!shape) return NULL;
            for (int i = 0; i < arr->ndim; i++) shape[i] = 1;
            NDArray* r = ndarray_full(arr->ndim, shape, arr->dtype, max_val);
            free(shape);
            return r;
        }
        return ndarray_scalar(max_val, arr->dtype);
    }

    if (axis < 0) axis += arr->ndim;
    if (axis < 0 || axis >= arr->ndim) return NULL;

    int32_t out_shape[32], out_ndim;
    get_reduced_shape(arr->shape, arr->ndim, axis, keepdims, out_shape, &out_ndim);

    NDArray* result = out_ndim == 0 ? ndarray_scalar(0, arr->dtype) :
                      ndarray_create(out_ndim, out_shape, arr->dtype);
    if (!result) return NULL;

    int32_t axis_size = arr->shape[axis];
    size_t outer = 1, inner = 1;
    for (int i = 0; i < axis; i++) outer *= arr->shape[i];
    for (int i = axis + 1; i < arr->ndim; i++) inner *= arr->shape[i];

    size_t res_idx = 0;
    for (size_t o = 0; o < outer; o++) {
        for (size_t n = 0; n < inner; n++) {
            double max_val = ndarray_get_flat(arr, o * axis_size * inner + n);
            for (int32_t i = 1; i < axis_size; i++) {
                double v = ndarray_get_flat(arr, o * axis_size * inner + i * inner + n);
                if (v > max_val) max_val = v;
            }
            ndarray_set_flat(result, res_idx++, max_val);
        }
    }
    return result;
}

EXPORT NDArray* ndarray_median(NDArray* arr, int32_t axis, bool keepdims) {
    if (!arr) return NULL;

    if (axis == -2147483648) {
        NDArray* flat = ndarray_flatten(arr);
        if (!flat) return NULL;
        ndarray_sort(flat, 0, 0);
        size_t n = flat->size, mid = n / 2;
        double median_val = (n % 2 == 0) ?
            (ndarray_get_flat(flat, mid-1) + ndarray_get_flat(flat, mid)) / 2.0 :
            ndarray_get_flat(flat, mid);
        ndarray_free(flat);

        if (keepdims) {
            int32_t* shape = malloc(arr->ndim * sizeof(int32_t));
            if (!shape) return NULL;
            for (int i = 0; i < arr->ndim; i++) shape[i] = 1;
            NDArray* r = ndarray_full(arr->ndim, shape, DTYPE_FLOAT64, median_val);
            free(shape);
            return r;
        }
        return ndarray_scalar(median_val, DTYPE_FLOAT64);
    }

    /* Axis-specific median implementation */
    if (axis < 0) axis += arr->ndim;
    if (axis < 0 || axis >= arr->ndim) return NULL;

    int32_t out_shape[32], out_ndim;
    get_reduced_shape(arr->shape, arr->ndim, axis, keepdims, out_shape, &out_ndim);

    NDArray* result = out_ndim == 0 ? ndarray_scalar(0, DTYPE_FLOAT64) :
                      ndarray_create(out_ndim, out_shape, DTYPE_FLOAT64);
    if (!result) return NULL;

    int32_t axis_size = arr->shape[axis];
    size_t mid = axis_size / 2;
    size_t outer = 1, inner = 1;
    for (int i = 0; i < axis; i++) outer *= arr->shape[i];
    for (int i = axis + 1; i < arr->ndim; i++) inner *= arr->shape[i];

    double* temp = malloc(axis_size * sizeof(double));
    if (!temp) { ndarray_free(result); return NULL; }
    size_t res_idx = 0;

    for (size_t o = 0; o < outer; o++) {
        for (size_t n = 0; n < inner; n++) {
            for (int32_t i = 0; i < axis_size; i++) {
                temp[i] = ndarray_get_flat(arr, o * axis_size * inner + i * inner + n);
            }
            /* Sort temp using insertion sort */
            for (int32_t i = 1; i < axis_size; i++) {
                double key = temp[i];
                int32_t j = i - 1;
                while (j >= 0 && temp[j] > key) { temp[j+1] = temp[j]; j--; }
                temp[j+1] = key;
            }
            double median_val = (axis_size % 2 == 0) ?
                (temp[mid-1] + temp[mid]) / 2.0 : temp[mid];
            ndarray_set_flat(result, res_idx++, median_val);
        }
    }
    free(temp);
    return result;
}

EXPORT NDArray* ndarray_quantile(NDArray* arr, double q, int32_t axis, bool keepdims) {
    if (!arr || q < 0 || q > 1) return NULL;

    if (axis == -2147483648) {
        NDArray* flat = ndarray_flatten(arr);
        if (!flat) return NULL;
        ndarray_sort(flat, 0, 0);
        size_t n = flat->size;
        double idx = q * (n - 1);
        size_t lo = (size_t)floor(idx), hi = (size_t)ceil(idx);
        double frac = idx - lo;
        double qval = (lo == hi) ? ndarray_get_flat(flat, lo) :
            ndarray_get_flat(flat, lo) * (1 - frac) + ndarray_get_flat(flat, hi) * frac;
        ndarray_free(flat);

        if (keepdims) {
            int32_t* shape = malloc(arr->ndim * sizeof(int32_t));
            if (!shape) return NULL;
            for (int i = 0; i < arr->ndim; i++) shape[i] = 1;
            NDArray* r = ndarray_full(arr->ndim, shape, DTYPE_FLOAT64, qval);
            free(shape);
            return r;
        }
        return ndarray_scalar(qval, DTYPE_FLOAT64);
    }
    return NULL;  /* Axis-specific implementation TODO */
}

EXPORT NDArray* ndarray_percentile(NDArray* arr, double q, int32_t axis, bool keepdims) {
    return ndarray_quantile(arr, q / 100.0, axis, keepdims);
}

/* NaN-aware implementations */
EXPORT NDArray* ndarray_nansum(NDArray* arr, int32_t axis, bool keepdims, DType dtype) {
    if (!arr) return NULL;
    DType out_dtype = get_reduction_dtype(arr->dtype, dtype);

    if (axis == -2147483648) {
        double total = 0;
        for (size_t i = 0; i < arr->size; i++) {
            double v = ndarray_get_flat(arr, i);
            if (!isnan(v)) total += v;
        }
        if (keepdims) {
            int32_t* shape = malloc(arr->ndim * sizeof(int32_t));
            if (!shape) return NULL;
            for (int i = 0; i < arr->ndim; i++) shape[i] = 1;
            NDArray* r = ndarray_full(arr->ndim, shape, out_dtype, total);
            free(shape);
            return r;
        }
        return ndarray_scalar(total, out_dtype);
    }
    return NULL;  /* Axis-specific TODO */
}

EXPORT NDArray* ndarray_nanmean(NDArray* arr, int32_t axis, bool keepdims, DType dtype) {
    if (!arr) return NULL;
    DType out_dtype = get_reduction_dtype(arr->dtype, dtype);

    if (axis == -2147483648) {
        double total = 0;
        size_t count = 0;
        for (size_t i = 0; i < arr->size; i++) {
            double v = ndarray_get_flat(arr, i);
            if (!isnan(v)) { total += v; count++; }
        }
        double mean_val = count > 0 ? total / count : NAN;
        if (keepdims) {
            int32_t* shape = malloc(arr->ndim * sizeof(int32_t));
            if (!shape) return NULL;
            for (int i = 0; i < arr->ndim; i++) shape[i] = 1;
            NDArray* r = ndarray_full(arr->ndim, shape, out_dtype, mean_val);
            free(shape);
            return r;
        }
        return ndarray_scalar(mean_val, out_dtype);
    }
    return NULL;
}

EXPORT NDArray* ndarray_nanvar(NDArray* arr, int32_t axis, bool keepdims,
                                DType dtype, int32_t ddof) {
    (void)arr; (void)axis; (void)keepdims; (void)dtype; (void)ddof;
    return NULL;  /* TODO */
}

EXPORT NDArray* ndarray_nanstd(NDArray* arr, int32_t axis, bool keepdims,
                                DType dtype, int32_t ddof) {
    (void)arr; (void)axis; (void)keepdims; (void)dtype; (void)ddof;
    return NULL;  /* TODO */
}

/* ============ Cumulative Operations (Phase 22) ============ */

EXPORT NDArray* ndarray_cumsum_axis(NDArray* arr, int32_t axis, DType dtype) {
    if (!arr || !arr->data) return NULL;
    DType out_dtype = get_reduction_dtype(arr->dtype, dtype);

    /* axis = -2147483648 means flatten and cumsum */
    if (axis == -2147483648) {
        int32_t shape[1] = { (int32_t)arr->size };
        NDArray* result = ndarray_create(1, shape, out_dtype);
        if (!result) return NULL;

        double sum = 0.0;
        for (size_t i = 0; i < arr->size; i++) {
            sum += ndarray_get_flat(arr, i);
            ndarray_set_flat(result, i, sum);
        }
        return result;
    }

    /* Normalize negative axis */
    if (axis < 0) axis += arr->ndim;
    if (axis < 0 || axis >= arr->ndim) return NULL;

    /* Output has same shape as input */
    NDArray* result = ndarray_create(arr->ndim, arr->shape, out_dtype);
    if (!result) return NULL;

    /* Compute cumsum along axis using outer/inner pattern */
    int32_t axis_size = arr->shape[axis];
    size_t outer = 1, inner = 1;
    for (int i = 0; i < axis; i++) outer *= arr->shape[i];
    for (int i = axis + 1; i < arr->ndim; i++) inner *= arr->shape[i];

    for (size_t o = 0; o < outer; o++) {
        for (size_t n = 0; n < inner; n++) {
            double sum = 0.0;
            for (int32_t i = 0; i < axis_size; i++) {
                size_t idx = o * axis_size * inner + i * inner + n;
                sum += ndarray_get_flat(arr, idx);
                ndarray_set_flat(result, idx, sum);
            }
        }
    }
    return result;
}

EXPORT NDArray* ndarray_cumprod_axis(NDArray* arr, int32_t axis, DType dtype) {
    if (!arr || !arr->data) return NULL;
    DType out_dtype = get_reduction_dtype(arr->dtype, dtype);

    /* axis = -2147483648 means flatten and cumprod */
    if (axis == -2147483648) {
        int32_t shape[1] = { (int32_t)arr->size };
        NDArray* result = ndarray_create(1, shape, out_dtype);
        if (!result) return NULL;

        double prod = 1.0;
        for (size_t i = 0; i < arr->size; i++) {
            prod *= ndarray_get_flat(arr, i);
            ndarray_set_flat(result, i, prod);
        }
        return result;
    }

    /* Normalize negative axis */
    if (axis < 0) axis += arr->ndim;
    if (axis < 0 || axis >= arr->ndim) return NULL;

    /* Output has same shape as input */
    NDArray* result = ndarray_create(arr->ndim, arr->shape, out_dtype);
    if (!result) return NULL;

    /* Compute cumprod along axis using outer/inner pattern */
    int32_t axis_size = arr->shape[axis];
    size_t outer = 1, inner = 1;
    for (int i = 0; i < axis; i++) outer *= arr->shape[i];
    for (int i = axis + 1; i < arr->ndim; i++) inner *= arr->shape[i];

    for (size_t o = 0; o < outer; o++) {
        for (size_t n = 0; n < inner; n++) {
            double prod = 1.0;
            for (int32_t i = 0; i < axis_size; i++) {
                size_t idx = o * axis_size * inner + i * inner + n;
                prod *= ndarray_get_flat(arr, idx);
                ndarray_set_flat(result, idx, prod);
            }
        }
    }
    return result;
}

EXPORT NDArray* ndarray_nancumsum_axis(NDArray* arr, int32_t axis, DType dtype) {
    if (!arr || !arr->data) return NULL;
    DType out_dtype = get_reduction_dtype(arr->dtype, dtype);

    /* axis = -2147483648 means flatten and nancumsum */
    if (axis == -2147483648) {
        int32_t shape[1] = { (int32_t)arr->size };
        NDArray* result = ndarray_create(1, shape, out_dtype);
        if (!result) return NULL;

        double sum = 0.0;
        for (size_t i = 0; i < arr->size; i++) {
            double val = ndarray_get_flat(arr, i);
            if (!isnan(val)) {
                sum += val;
            }
            ndarray_set_flat(result, i, sum);
        }
        return result;
    }

    /* Normalize negative axis */
    if (axis < 0) axis += arr->ndim;
    if (axis < 0 || axis >= arr->ndim) return NULL;

    /* Output has same shape as input */
    NDArray* result = ndarray_create(arr->ndim, arr->shape, out_dtype);
    if (!result) return NULL;

    /* Compute nancumsum along axis using outer/inner pattern */
    int32_t axis_size = arr->shape[axis];
    size_t outer = 1, inner = 1;
    for (int i = 0; i < axis; i++) outer *= arr->shape[i];
    for (int i = axis + 1; i < arr->ndim; i++) inner *= arr->shape[i];

    for (size_t o = 0; o < outer; o++) {
        for (size_t n = 0; n < inner; n++) {
            double sum = 0.0;
            for (int32_t i = 0; i < axis_size; i++) {
                size_t idx = o * axis_size * inner + i * inner + n;
                double val = ndarray_get_flat(arr, idx);
                if (!isnan(val)) {
                    sum += val;
                }
                ndarray_set_flat(result, idx, sum);
            }
        }
    }
    return result;
}

EXPORT NDArray* ndarray_nancumprod_axis(NDArray* arr, int32_t axis, DType dtype) {
    if (!arr || !arr->data) return NULL;
    DType out_dtype = get_reduction_dtype(arr->dtype, dtype);

    /* axis = -2147483648 means flatten and nancumprod */
    if (axis == -2147483648) {
        int32_t shape[1] = { (int32_t)arr->size };
        NDArray* result = ndarray_create(1, shape, out_dtype);
        if (!result) return NULL;

        double prod = 1.0;
        for (size_t i = 0; i < arr->size; i++) {
            double val = ndarray_get_flat(arr, i);
            if (!isnan(val)) {
                prod *= val;
            }
            ndarray_set_flat(result, i, prod);
        }
        return result;
    }

    /* Normalize negative axis */
    if (axis < 0) axis += arr->ndim;
    if (axis < 0 || axis >= arr->ndim) return NULL;

    /* Output has same shape as input */
    NDArray* result = ndarray_create(arr->ndim, arr->shape, out_dtype);
    if (!result) return NULL;

    /* Compute nancumprod along axis using outer/inner pattern */
    int32_t axis_size = arr->shape[axis];
    size_t outer = 1, inner = 1;
    for (int i = 0; i < axis; i++) outer *= arr->shape[i];
    for (int i = axis + 1; i < arr->ndim; i++) inner *= arr->shape[i];

    for (size_t o = 0; o < outer; o++) {
        for (size_t n = 0; n < inner; n++) {
            double prod = 1.0;
            for (int32_t i = 0; i < axis_size; i++) {
                size_t idx = o * axis_size * inner + i * inner + n;
                double val = ndarray_get_flat(arr, idx);
                if (!isnan(val)) {
                    prod *= val;
                }
                ndarray_set_flat(result, idx, prod);
            }
        }
    }
    return result;
}

/* ============================================================================
 * Product Reductions (prod, nanprod)
 * ============================================================================ */

EXPORT NDArray* ndarray_prod_axis(NDArray* arr, int32_t axis, bool keepdims, DType dtype) {
    if (!arr || !arr->data) return NULL;
    DType out_dtype = get_reduction_dtype(arr->dtype, dtype);

    /* Global product (all elements) */
    if (axis == -2147483648) {
        double total = 1.0;
        for (size_t i = 0; i < arr->size; i++) {
            total *= ndarray_get_flat(arr, i);
        }
        if (keepdims) {
            int32_t* shape = malloc(arr->ndim * sizeof(int32_t));
            if (!shape) return NULL;
            for (int i = 0; i < arr->ndim; i++) shape[i] = 1;
            NDArray* r = ndarray_full(arr->ndim, shape, out_dtype, total);
            free(shape);
            return r;
        }
        return ndarray_scalar(total, out_dtype);
    }

    /* Normalize negative axis */
    if (axis < 0) axis += arr->ndim;
    if (axis < 0 || axis >= arr->ndim) return NULL;

    int32_t out_shape[32], out_ndim;
    get_reduced_shape(arr->shape, arr->ndim, axis, keepdims, out_shape, &out_ndim);

    NDArray* result = out_ndim == 0 ? ndarray_scalar(1.0, out_dtype) :
                      ndarray_create(out_ndim, out_shape, out_dtype);
    if (!result) return NULL;

    int32_t axis_size = arr->shape[axis];
    size_t outer = 1, inner = 1;
    for (int i = 0; i < axis; i++) outer *= arr->shape[i];
    for (int i = axis + 1; i < arr->ndim; i++) inner *= arr->shape[i];

    size_t res_idx = 0;
    for (size_t o = 0; o < outer; o++) {
        for (size_t n = 0; n < inner; n++) {
            double prod = 1.0;
            for (int32_t i = 0; i < axis_size; i++) {
                size_t idx = o * axis_size * inner + i * inner + n;
                prod *= ndarray_get_flat(arr, idx);
            }
            ndarray_set_flat(result, res_idx++, prod);
        }
    }
    return result;
}

EXPORT NDArray* ndarray_nanprod(NDArray* arr, int32_t axis, bool keepdims, DType dtype) {
    if (!arr || !arr->data) return NULL;
    DType out_dtype = get_reduction_dtype(arr->dtype, dtype);

    /* Global nanprod (all elements, skipping NaN) */
    if (axis == -2147483648) {
        double total = 1.0;
        for (size_t i = 0; i < arr->size; i++) {
            double val = ndarray_get_flat(arr, i);
            if (!isnan(val)) {
                total *= val;
            }
        }
        if (keepdims) {
            int32_t* shape = malloc(arr->ndim * sizeof(int32_t));
            if (!shape) return NULL;
            for (int i = 0; i < arr->ndim; i++) shape[i] = 1;
            NDArray* r = ndarray_full(arr->ndim, shape, out_dtype, total);
            free(shape);
            return r;
        }
        return ndarray_scalar(total, out_dtype);
    }

    /* Normalize negative axis */
    if (axis < 0) axis += arr->ndim;
    if (axis < 0 || axis >= arr->ndim) return NULL;

    int32_t out_shape[32], out_ndim;
    get_reduced_shape(arr->shape, arr->ndim, axis, keepdims, out_shape, &out_ndim);

    NDArray* result = out_ndim == 0 ? ndarray_scalar(1.0, out_dtype) :
                      ndarray_create(out_ndim, out_shape, out_dtype);
    if (!result) return NULL;

    int32_t axis_size = arr->shape[axis];
    size_t outer = 1, inner = 1;
    for (int i = 0; i < axis; i++) outer *= arr->shape[i];
    for (int i = axis + 1; i < arr->ndim; i++) inner *= arr->shape[i];

    size_t res_idx = 0;
    for (size_t o = 0; o < outer; o++) {
        for (size_t n = 0; n < inner; n++) {
            double prod = 1.0;
            for (int32_t i = 0; i < axis_size; i++) {
                size_t idx = o * axis_size * inner + i * inner + n;
                double val = ndarray_get_flat(arr, idx);
                if (!isnan(val)) {
                    prod *= val;
                }
            }
            ndarray_set_flat(result, res_idx++, prod);
        }
    }
    return result;
}

/* ============================================================================
 * Peak-to-Peak (ptp = max - min)
 * ============================================================================ */

EXPORT NDArray* ndarray_ptp_axis(NDArray* arr, int32_t axis, bool keepdims) {
    if (!arr || !arr->data) return NULL;

    /* Global ptp */
    if (axis == -2147483648) {
        double min_val = ndarray_get_flat(arr, 0);
        double max_val = min_val;
        for (size_t i = 1; i < arr->size; i++) {
            double v = ndarray_get_flat(arr, i);
            if (v < min_val) min_val = v;
            if (v > max_val) max_val = v;
        }
        double ptp_val = max_val - min_val;
        if (keepdims) {
            int32_t* shape = malloc(arr->ndim * sizeof(int32_t));
            if (!shape) return NULL;
            for (int i = 0; i < arr->ndim; i++) shape[i] = 1;
            NDArray* r = ndarray_full(arr->ndim, shape, arr->dtype, ptp_val);
            free(shape);
            return r;
        }
        return ndarray_scalar(ptp_val, arr->dtype);
    }

    /* Normalize negative axis */
    if (axis < 0) axis += arr->ndim;
    if (axis < 0 || axis >= arr->ndim) return NULL;

    int32_t out_shape[32], out_ndim;
    get_reduced_shape(arr->shape, arr->ndim, axis, keepdims, out_shape, &out_ndim);

    NDArray* result = out_ndim == 0 ? ndarray_scalar(0, arr->dtype) :
                      ndarray_create(out_ndim, out_shape, arr->dtype);
    if (!result) return NULL;

    int32_t axis_size = arr->shape[axis];
    size_t outer = 1, inner = 1;
    for (int i = 0; i < axis; i++) outer *= arr->shape[i];
    for (int i = axis + 1; i < arr->ndim; i++) inner *= arr->shape[i];

    size_t res_idx = 0;
    for (size_t o = 0; o < outer; o++) {
        for (size_t n = 0; n < inner; n++) {
            double min_val = ndarray_get_flat(arr, o * axis_size * inner + n);
            double max_val = min_val;
            for (int32_t i = 1; i < axis_size; i++) {
                double v = ndarray_get_flat(arr, o * axis_size * inner + i * inner + n);
                if (v < min_val) min_val = v;
                if (v > max_val) max_val = v;
            }
            ndarray_set_flat(result, res_idx++, max_val - min_val);
        }
    }
    return result;
}

/* ============================================================================
 * Differences (diff, ediff1d, gradient)
 * ============================================================================ */

EXPORT NDArray* ndarray_diff(NDArray* arr, int32_t n, int32_t axis) {
    if (!arr || !arr->data) return NULL;
    if (n < 0) return NULL;
    if (n == 0) return ndarray_copy(arr);
    if (arr->ndim == 0) return NULL;  /* diff requires at least 1D */

    /* Normalize negative axis */
    if (axis < 0) axis += arr->ndim;
    if (axis < 0 || axis >= arr->ndim) return NULL;

    /* Check if axis size is sufficient */
    if (arr->shape[axis] <= n) {
        /* Would result in size 0 along axis */
        int32_t out_shape[32];
        for (int i = 0; i < arr->ndim; i++) {
            out_shape[i] = (i == axis) ? 0 : arr->shape[i];
        }
        return ndarray_create(arr->ndim, out_shape, arr->dtype);
    }

    NDArray* current = arr;
    NDArray* prev = NULL;

    for (int32_t iter = 0; iter < n; iter++) {
        int32_t out_shape[32];
        for (int i = 0; i < current->ndim; i++) {
            out_shape[i] = (i == axis) ? current->shape[i] - 1 : current->shape[i];
        }

        NDArray* result = ndarray_create(current->ndim, out_shape, current->dtype);
        if (!result) {
            if (prev) ndarray_free(prev);
            return NULL;
        }

        int32_t axis_size = out_shape[axis];  /* new axis size */
        size_t outer = 1, inner = 1;
        for (int i = 0; i < axis; i++) outer *= current->shape[i];
        for (int i = axis + 1; i < current->ndim; i++) inner *= current->shape[i];

        size_t res_idx = 0;
        int32_t old_axis_size = current->shape[axis];
        for (size_t o = 0; o < outer; o++) {
            for (size_t ni = 0; ni < inner; ni++) {
                for (int32_t i = 0; i < axis_size; i++) {
                    size_t idx1 = o * old_axis_size * inner + (i + 1) * inner + ni;
                    size_t idx0 = o * old_axis_size * inner + i * inner + ni;
                    double diff_val = ndarray_get_flat(current, idx1) -
                                     ndarray_get_flat(current, idx0);
                    ndarray_set_flat(result, res_idx++, diff_val);
                }
            }
        }

        if (prev) ndarray_free(prev);
        prev = (current == arr) ? NULL : current;
        current = result;
    }

    return current;
}

EXPORT NDArray* ndarray_ediff1d(NDArray* arr, NDArray* to_begin, NDArray* to_end) {
    if (!arr || !arr->data) return NULL;

    size_t n = arr->size;
    size_t diff_size = (n > 0) ? n - 1 : 0;
    size_t begin_size = to_begin ? to_begin->size : 0;
    size_t end_size = to_end ? to_end->size : 0;
    size_t total_size = begin_size + diff_size + end_size;

    if (total_size == 0) {
        int32_t out_shape[1] = { 0 };
        return ndarray_create(1, out_shape, arr->dtype);
    }

    int32_t out_shape[1] = { (int32_t)total_size };
    NDArray* result = ndarray_create(1, out_shape, arr->dtype);
    if (!result) return NULL;

    size_t out_idx = 0;

    /* Prepend to_begin */
    if (to_begin) {
        for (size_t i = 0; i < begin_size; i++) {
            ndarray_set_flat(result, out_idx++, ndarray_get_flat(to_begin, i));
        }
    }

    /* Compute differences on flattened array */
    for (size_t i = 1; i < n; i++) {
        double diff_val = ndarray_get_flat(arr, i) - ndarray_get_flat(arr, i - 1);
        ndarray_set_flat(result, out_idx++, diff_val);
    }

    /* Append to_end */
    if (to_end) {
        for (size_t i = 0; i < end_size; i++) {
            ndarray_set_flat(result, out_idx++, ndarray_get_flat(to_end, i));
        }
    }

    return result;
}

EXPORT NDArray* ndarray_gradient(NDArray* arr, double dx, int32_t axis) {
    if (!arr || !arr->data || dx == 0.0) return NULL;
    if (arr->ndim == 0) return NULL;

    /* Normalize negative axis */
    if (axis < 0) axis += arr->ndim;
    if (axis < 0 || axis >= arr->ndim) return NULL;

    /* Output has same shape as input */
    NDArray* result = ndarray_create(arr->ndim, arr->shape, DTYPE_FLOAT64);
    if (!result) return NULL;

    int32_t axis_size = arr->shape[axis];
    if (axis_size < 2) {
        /* Cannot compute gradient with fewer than 2 points */
        ndarray_fill(result, 0.0);
        return result;
    }

    size_t outer = 1, inner = 1;
    for (int i = 0; i < axis; i++) outer *= arr->shape[i];
    for (int i = axis + 1; i < arr->ndim; i++) inner *= arr->shape[i];

    double inv_2dx = 1.0 / (2.0 * dx);
    double inv_dx = 1.0 / dx;

    for (size_t o = 0; o < outer; o++) {
        for (size_t n = 0; n < inner; n++) {
            /* Left boundary (first order forward difference) */
            {
                size_t idx0 = o * axis_size * inner + 0 * inner + n;
                size_t idx1 = o * axis_size * inner + 1 * inner + n;
                double grad = (ndarray_get_flat(arr, idx1) - ndarray_get_flat(arr, idx0)) * inv_dx;
                ndarray_set_flat(result, idx0, grad);
            }

            /* Interior points (central differences) */
            for (int32_t i = 1; i < axis_size - 1; i++) {
                size_t idx_m = o * axis_size * inner + (i - 1) * inner + n;
                size_t idx_p = o * axis_size * inner + (i + 1) * inner + n;
                size_t idx = o * axis_size * inner + i * inner + n;
                double grad = (ndarray_get_flat(arr, idx_p) - ndarray_get_flat(arr, idx_m)) * inv_2dx;
                ndarray_set_flat(result, idx, grad);
            }

            /* Right boundary (first order backward difference) */
            {
                size_t idx_m1 = o * axis_size * inner + (axis_size - 2) * inner + n;
                size_t idx = o * axis_size * inner + (axis_size - 1) * inner + n;
                double grad = (ndarray_get_flat(arr, idx) - ndarray_get_flat(arr, idx_m1)) * inv_dx;
                ndarray_set_flat(result, idx, grad);
            }
        }
    }

    return result;
}

/* ============================================================================
 * Covariance and Correlation
 * ============================================================================ */

EXPORT NDArray* ndarray_cov(NDArray* m, NDArray* y, bool rowvar, bool bias, int32_t ddof) {
    if (!m || !m->data) return NULL;

    /* Convert 1D to 2D */
    NDArray* X = NULL;
    bool free_X = false;

    if (m->ndim == 1) {
        int32_t new_shape[2];
        if (rowvar) {
            new_shape[0] = 1;
            new_shape[1] = m->shape[0];
        } else {
            new_shape[0] = m->shape[0];
            new_shape[1] = 1;
        }
        X = ndarray_reshape(m, 2, new_shape);
        if (!X) return NULL;
        free_X = true;
    } else if (m->ndim == 2) {
        X = m;
    } else {
        return NULL;  /* Only 1D or 2D supported */
    }

    /* If not rowvar, transpose */
    if (!rowvar) {
        NDArray* Xt = ndarray_transpose(X, NULL);
        if (free_X) ndarray_free(X);
        if (!Xt) return NULL;
        X = Xt;
        free_X = true;
    }

    /* Now X has shape (p, n) where p = variables, n = observations */
    int32_t p = X->shape[0];
    int32_t n = X->shape[1];

    if (n < 1) {
        if (free_X) ndarray_free(X);
        return NULL;
    }

    /* Determine ddof */
    if (ddof == -1) {
        ddof = bias ? 0 : 1;
    }

    /* Compute mean for each row (variable) */
    double* means = malloc(p * sizeof(double));
    if (!means) {
        if (free_X) ndarray_free(X);
        return NULL;
    }
    for (int32_t i = 0; i < p; i++) {
        double sum = 0.0;
        for (int32_t j = 0; j < n; j++) {
            sum += ndarray_get_flat(X, i * n + j);
        }
        means[i] = sum / n;
    }

    /* Create centered data */
    double* centered = malloc(p * n * sizeof(double));
    if (!centered) {
        free(means);
        if (free_X) ndarray_free(X);
        return NULL;
    }
    for (int32_t i = 0; i < p; i++) {
        for (int32_t j = 0; j < n; j++) {
            centered[i * n + j] = ndarray_get_flat(X, i * n + j) - means[i];
        }
    }
    free(means);

    /* Create output covariance matrix (p x p) */
    int32_t out_shape[2] = { p, p };
    NDArray* result = ndarray_create(2, out_shape, DTYPE_FLOAT64);
    if (!result) {
        free(centered);
        if (free_X) ndarray_free(X);
        return NULL;
    }

    /* Compute C = centered @ centered.T / (n - ddof) */
    double scale = (n - ddof > 0) ? 1.0 / (double)(n - ddof) : 0.0;
    double* C = (double*)result->data;

    for (int32_t i = 0; i < p; i++) {
        for (int32_t j = i; j < p; j++) {
            double dot = 0.0;
            for (int32_t k = 0; k < n; k++) {
                dot += centered[i * n + k] * centered[j * n + k];
            }
            C[i * p + j] = dot * scale;
            C[j * p + i] = C[i * p + j];  /* Symmetric */
        }
    }

    free(centered);
    if (free_X) ndarray_free(X);

    return result;
}

EXPORT NDArray* ndarray_corrcoef(NDArray* x, NDArray* y, bool rowvar) {
    /* Compute covariance matrix */
    NDArray* C = ndarray_cov(x, y, rowvar, false, 1);
    if (!C) return NULL;

    int32_t p = C->shape[0];
    double* data = (double*)C->data;

    /* Extract diagonal (variances) */
    double* d = malloc(p * sizeof(double));
    if (!d) {
        ndarray_free(C);
        return NULL;
    }
    for (int32_t i = 0; i < p; i++) {
        d[i] = data[i * p + i];
    }

    /* Normalize: R[i,j] = C[i,j] / sqrt(C[i,i] * C[j,j]) */
    for (int32_t i = 0; i < p; i++) {
        for (int32_t j = 0; j < p; j++) {
            double denom = sqrt(d[i] * d[j]);
            if (denom > 0) {
                data[i * p + j] /= denom;
            } else {
                data[i * p + j] = (i == j) ? 1.0 : 0.0;
            }
        }
    }

    /* Clip to [-1, 1] to handle floating point errors */
    for (int32_t i = 0; i < p * p; i++) {
        if (data[i] > 1.0) data[i] = 1.0;
        if (data[i] < -1.0) data[i] = -1.0;
    }

    free(d);
    return C;
}

/* ============================================================================
 * Histogram
 * ============================================================================ */

EXPORT HistogramResult* ndarray_histogram(NDArray* a, int32_t bins,
                                          double range_min, double range_max,
                                          NDArray* weights, bool density) {
    if (!a || !a->data || bins <= 0) return NULL;

    size_t n = a->size;

    /* Determine range if not specified (use NAN to indicate "auto") */
    if (isnan(range_min) || isnan(range_max)) {
        if (n == 0) return NULL;
        double data_min = ndarray_get_flat(a, 0);
        double data_max = data_min;
        for (size_t i = 1; i < n; i++) {
            double v = ndarray_get_flat(a, i);
            if (v < data_min) data_min = v;
            if (v > data_max) data_max = v;
        }
        if (isnan(range_min)) range_min = data_min;
        if (isnan(range_max)) range_max = data_max;
    }

    /* Handle edge case where all values are the same */
    if (range_max == range_min) {
        range_min -= 0.5;
        range_max += 0.5;
    }

    if (range_max < range_min) return NULL;

    double bin_width = (range_max - range_min) / bins;
    double inv_bin_width = 1.0 / bin_width;

    /* Create result structure */
    HistogramResult* result = malloc(sizeof(HistogramResult));
    if (!result) return NULL;

    int32_t hist_shape[1] = { bins };
    int32_t edge_shape[1] = { bins + 1 };

    result->hist = ndarray_create(1, hist_shape, DTYPE_FLOAT64);
    result->bin_edges = ndarray_create(1, edge_shape, DTYPE_FLOAT64);

    if (!result->hist || !result->bin_edges) {
        if (result->hist) ndarray_free(result->hist);
        if (result->bin_edges) ndarray_free(result->bin_edges);
        free(result);
        return NULL;
    }

    /* Fill bin edges */
    for (int32_t i = 0; i <= bins; i++) {
        ndarray_set_flat(result->bin_edges, i, range_min + i * bin_width);
    }

    /* Initialize histogram to zeros */
    ndarray_fill(result->hist, 0.0);

    /* Bin the data */
    for (size_t i = 0; i < n; i++) {
        double val = ndarray_get_flat(a, i);

        /* Skip values outside range */
        if (val < range_min || val > range_max) continue;

        /* Compute bin index */
        int32_t bin_idx = (int32_t)((val - range_min) * inv_bin_width);

        /* Handle right edge (last bin includes right boundary) */
        if (bin_idx >= bins) bin_idx = bins - 1;
        if (bin_idx < 0) bin_idx = 0;

        /* Add to bin */
        if (weights && weights->size > i) {
            double w = ndarray_get_flat(weights, i);
            double current = ndarray_get_flat(result->hist, bin_idx);
            ndarray_set_flat(result->hist, bin_idx, current + w);
        } else {
            double current = ndarray_get_flat(result->hist, bin_idx);
            ndarray_set_flat(result->hist, bin_idx, current + 1.0);
        }
    }

    /* Normalize if density requested */
    if (density) {
        double total = 0.0;
        for (int32_t i = 0; i < bins; i++) {
            total += ndarray_get_flat(result->hist, i);
        }
        if (total > 0) {
            double norm = total * bin_width;
            for (int32_t i = 0; i < bins; i++) {
                double v = ndarray_get_flat(result->hist, i);
                ndarray_set_flat(result->hist, i, v / norm);
            }
        }
    }

    return result;
}

EXPORT void histogram_result_free(HistogramResult* result) {
    if (!result) return;
    if (result->hist) ndarray_free(result->hist);
    if (result->bin_edges) ndarray_free(result->bin_edges);
    free(result);
}
