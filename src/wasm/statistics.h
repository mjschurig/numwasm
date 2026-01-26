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

#endif /* NUMJS_STATISTICS_H */
