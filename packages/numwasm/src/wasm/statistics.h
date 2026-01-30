#ifndef NUMJS_STATISTICS_H
#define NUMJS_STATISTICS_H

#include "ndarray.h"

/* Basic reductions with axis support */
NDArray* ndarray_sum_axis(NDArray* arr, int32_t axis, bool keepdims, DType dtype);
NDArray* ndarray_mean_axis(NDArray* arr, int32_t axis, bool keepdims, DType dtype);
NDArray* ndarray_var_axis(NDArray* arr, int32_t axis, bool keepdims, DType dtype, int32_t ddof);
NDArray* ndarray_std_axis(NDArray* arr, int32_t axis, bool keepdims, DType dtype, int32_t ddof);
NDArray* ndarray_min_axis(NDArray* arr, int32_t axis, bool keepdims);
NDArray* ndarray_max_axis(NDArray* arr, int32_t axis, bool keepdims);

/* Advanced statistics */
NDArray* ndarray_median(NDArray* arr, int32_t axis, bool keepdims);
NDArray* ndarray_percentile(NDArray* arr, double q, int32_t axis, bool keepdims);
NDArray* ndarray_quantile(NDArray* arr, double q, int32_t axis, bool keepdims);

/* NaN-aware versions */
NDArray* ndarray_nansum(NDArray* arr, int32_t axis, bool keepdims, DType dtype);
NDArray* ndarray_nanmean(NDArray* arr, int32_t axis, bool keepdims, DType dtype);
NDArray* ndarray_nanvar(NDArray* arr, int32_t axis, bool keepdims, DType dtype, int32_t ddof);
NDArray* ndarray_nanstd(NDArray* arr, int32_t axis, bool keepdims, DType dtype, int32_t ddof);

/* Cumulative operations (Phase 22) */
NDArray* ndarray_cumsum_axis(NDArray* arr, int32_t axis, DType dtype);
NDArray* ndarray_cumprod_axis(NDArray* arr, int32_t axis, DType dtype);
NDArray* ndarray_nancumsum_axis(NDArray* arr, int32_t axis, DType dtype);
NDArray* ndarray_nancumprod_axis(NDArray* arr, int32_t axis, DType dtype);

/* Product reductions */
NDArray* ndarray_prod_axis(NDArray* arr, int32_t axis, bool keepdims, DType dtype);
NDArray* ndarray_nanprod(NDArray* arr, int32_t axis, bool keepdims, DType dtype);

/* Peak-to-peak (max - min) */
NDArray* ndarray_ptp_axis(NDArray* arr, int32_t axis, bool keepdims);

/* Differences and gradients */
NDArray* ndarray_diff(NDArray* arr, int32_t n, int32_t axis);
NDArray* ndarray_ediff1d(NDArray* arr, NDArray* to_begin, NDArray* to_end);
NDArray* ndarray_gradient(NDArray* arr, double dx, int32_t axis);

/* Covariance and correlation */
NDArray* ndarray_cov(NDArray* m, NDArray* y, bool rowvar, bool bias, int32_t ddof);
NDArray* ndarray_corrcoef(NDArray* x, NDArray* y, bool rowvar);

/* Histogram */
typedef struct {
    NDArray* hist;
    NDArray* bin_edges;
} HistogramResult;

HistogramResult* ndarray_histogram(NDArray* a, int32_t bins,
                                   double range_min, double range_max,
                                   NDArray* weights, bool density);
void histogram_result_free(HistogramResult* result);

#endif /* NUMJS_STATISTICS_H */
