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
