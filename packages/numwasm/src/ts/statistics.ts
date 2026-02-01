/**
 * Statistics and searching functions for NDArray
 *
 * Provides WASM-accelerated statistical operations with full NumPy-compatible
 * axis support, keepdims, and dtype parameters.
 */

import { NDArray } from "./_core/NDArray.js";
import { DType } from "./types.js";

/**
 * Sum of array elements over a given axis.
 *
 * @param a - Input array
 * @param axis - Axis along which to sum. null sums all elements.
 * @param keepdims - If true, reduced axes are left with size 1
 * @param dtype - Output dtype (default: float64 for integers, same as input otherwise)
 * @returns Sum as number (if axis=null and keepdims=false) or NDArray
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([[1, 2], [3, 4]]);
 * const total = await sum(arr);           // 10
 * const rowSums = await sum(arr, 1);      // [3, 7]
 * const colSums = await sum(arr, 0, true); // [[4, 6]]
 * ```
 */
export async function sum(
  a: NDArray,
  axis: number | null = null,
  keepdims = false,
  dtype?: DType,
): Promise<NDArray | number> {
  const axisVal = axis === null ? -2147483648 : axis;
  const resultPtr = a._wasmModule._ndarray_sum_axis(
    a._wasmPtr,
    axisVal,
    keepdims,
    dtype ?? -1,
  );
  if (resultPtr === 0) throw new Error("sum failed");
  const result = NDArray._fromPtr(resultPtr, a._wasmModule);
  if (axis === null && !keepdims) {
    const v = result.item();
    result.dispose();
    return v;
  }
  return result;
}

/**
 * Compute the arithmetic mean along the specified axis.
 *
 * @param a - Input array
 * @param axis - Axis along which to compute mean. null computes mean of all elements.
 * @param keepdims - If true, reduced axes are left with size 1
 * @param dtype - Output dtype (default: float64 for integers, same as input otherwise)
 * @returns Mean as number (if axis=null and keepdims=false) or NDArray
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([1, 2, 3, 4, 5]);
 * const avg = await mean(arr);  // 3.0
 * ```
 */
export async function mean(
  a: NDArray,
  axis: number | null = null,
  keepdims = false,
  dtype?: DType,
): Promise<NDArray | number> {
  const axisVal = axis === null ? -2147483648 : axis;
  const resultPtr = a._wasmModule._ndarray_mean_axis(
    a._wasmPtr,
    axisVal,
    keepdims,
    dtype ?? -1,
  );
  if (resultPtr === 0) throw new Error("mean failed");
  const result = NDArray._fromPtr(resultPtr, a._wasmModule);
  if (axis === null && !keepdims) {
    const v = result.item();
    result.dispose();
    return v;
  }
  return result;
}

/**
 * Compute the variance along the specified axis.
 *
 * @param a - Input array
 * @param axis - Axis along which to compute variance. null computes variance of all elements.
 * @param ddof - Delta Degrees of Freedom. Divisor is N - ddof. Default: 0 (population variance)
 * @param keepdims - If true, reduced axes are left with size 1
 * @param dtype - Output dtype
 * @returns Variance as number (if axis=null and keepdims=false) or NDArray
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([1, 2, 3, 4, 5]);
 * const v = await variance(arr);           // 2.0 (population variance)
 * const s = await variance(arr, null, 1);  // 2.5 (sample variance)
 * ```
 */
export async function variance(
  a: NDArray,
  axis: number | null = null,
  ddof = 0,
  keepdims = false,
  dtype?: DType,
): Promise<NDArray | number> {
  const axisVal = axis === null ? -2147483648 : axis;
  const resultPtr = a._wasmModule._ndarray_var_axis(
    a._wasmPtr,
    axisVal,
    keepdims,
    dtype ?? -1,
    ddof,
  );
  if (resultPtr === 0) throw new Error("variance failed");
  const result = NDArray._fromPtr(resultPtr, a._wasmModule);
  if (axis === null && !keepdims) {
    const v = result.item();
    result.dispose();
    return v;
  }
  return result;
}

/** Alias for variance (matches NumPy's np.var) */
export const var_ = variance;

/**
 * Compute the standard deviation along the specified axis.
 *
 * @param a - Input array
 * @param axis - Axis along which to compute std. null computes std of all elements.
 * @param ddof - Delta Degrees of Freedom. Divisor is N - ddof. Default: 0
 * @param keepdims - If true, reduced axes are left with size 1
 * @param dtype - Output dtype
 * @returns Standard deviation as number (if axis=null and keepdims=false) or NDArray
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([1, 2, 3, 4, 5]);
 * const s = await std(arr);  // ~1.414
 * ```
 */
export async function std(
  a: NDArray,
  axis: number | null = null,
  ddof = 0,
  keepdims = false,
  dtype?: DType,
): Promise<NDArray | number> {
  const axisVal = axis === null ? -2147483648 : axis;
  const resultPtr = a._wasmModule._ndarray_std_axis(
    a._wasmPtr,
    axisVal,
    keepdims,
    dtype ?? -1,
    ddof,
  );
  if (resultPtr === 0) throw new Error("std failed");
  const result = NDArray._fromPtr(resultPtr, a._wasmModule);
  if (axis === null && !keepdims) {
    const v = result.item();
    result.dispose();
    return v;
  }
  return result;
}

/**
 * Return the minimum of an array or minimum along an axis.
 *
 * @param a - Input array
 * @param axis - Axis along which to find minimum. null finds minimum of all elements.
 * @param keepdims - If true, reduced axes are left with size 1
 * @returns Minimum as number (if axis=null and keepdims=false) or NDArray
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([[3, 1], [2, 4]]);
 * const minVal = await min(arr);      // 1
 * const rowMins = await min(arr, 1);  // [1, 2]
 * ```
 */
export async function min(
  a: NDArray,
  axis: number | null = null,
  keepdims = false,
): Promise<NDArray | number> {
  const axisVal = axis === null ? -2147483648 : axis;
  const resultPtr = a._wasmModule._ndarray_min_axis(
    a._wasmPtr,
    axisVal,
    keepdims,
  );
  if (resultPtr === 0) throw new Error("min failed");
  const result = NDArray._fromPtr(resultPtr, a._wasmModule);
  if (axis === null && !keepdims) {
    const v = result.item();
    result.dispose();
    return v;
  }
  return result;
}

/**
 * Return the maximum of an array or maximum along an axis.
 *
 * @param a - Input array
 * @param axis - Axis along which to find maximum. null finds maximum of all elements.
 * @param keepdims - If true, reduced axes are left with size 1
 * @returns Maximum as number (if axis=null and keepdims=false) or NDArray
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([[3, 1], [2, 4]]);
 * const maxVal = await max(arr);      // 4
 * const rowMaxs = await max(arr, 1);  // [3, 4]
 * ```
 */
export async function max(
  a: NDArray,
  axis: number | null = null,
  keepdims = false,
): Promise<NDArray | number> {
  const axisVal = axis === null ? -2147483648 : axis;
  const resultPtr = a._wasmModule._ndarray_max_axis(
    a._wasmPtr,
    axisVal,
    keepdims,
  );
  if (resultPtr === 0) throw new Error("max failed");
  const result = NDArray._fromPtr(resultPtr, a._wasmModule);
  if (axis === null && !keepdims) {
    const v = result.item();
    result.dispose();
    return v;
  }
  return result;
}

/**
 * Compute the median along the specified axis.
 *
 * @param a - Input array
 * @param axis - Axis along which to compute median. null computes median of all elements.
 * @param keepdims - If true, reduced axes are left with size 1
 * @returns Median as number (if axis=null and keepdims=false) or NDArray
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([1, 3, 2, 4, 5]);
 * const med = await median(arr);  // 3.0
 * ```
 */
export async function median(
  a: NDArray,
  axis: number | null = null,
  keepdims = false,
): Promise<NDArray | number> {
  const axisVal = axis === null ? -2147483648 : axis;
  const resultPtr = a._wasmModule._ndarray_median(
    a._wasmPtr,
    axisVal,
    keepdims,
  );
  if (resultPtr === 0) throw new Error("median failed");
  const result = NDArray._fromPtr(resultPtr, a._wasmModule);
  if (axis === null && !keepdims) {
    const v = result.item();
    result.dispose();
    return v;
  }
  return result;
}

/**
 * Returns the indices of the maximum values along an axis.
 *
 * @param a - Input array
 * @param axis - Axis along which to find argmax. null finds argmax of flattened array.
 * @param keepdims - If true, reduced axes are left with size 1
 * @returns Index as number (if axis=null and keepdims=false) or NDArray of indices
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([[1, 4], [3, 2]]);
 * const idx = await argmax(arr);       // 1 (index in flattened array)
 * const rowIdx = await argmax(arr, 1); // [1, 0]
 * ```
 */
export async function argmax(
  a: NDArray,
  axis: number | null = null,
  keepdims = false,
): Promise<NDArray | number> {
  const axisVal = axis === null ? -2147483648 : axis;
  const resultPtr = a._wasmModule._ndarray_argmax(
    a._wasmPtr,
    axisVal,
    keepdims,
  );
  if (resultPtr === 0) throw new Error("argmax failed");
  const result = NDArray._fromPtr(resultPtr, a._wasmModule);
  if (axis === null && !keepdims) {
    const v = result.item();
    result.dispose();
    return v;
  }
  return result;
}

/**
 * Returns the indices of the minimum values along an axis.
 *
 * @param a - Input array
 * @param axis - Axis along which to find argmin. null finds argmin of flattened array.
 * @param keepdims - If true, reduced axes are left with size 1
 * @returns Index as number (if axis=null and keepdims=false) or NDArray of indices
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([[3, 1], [2, 4]]);
 * const idx = await argmin(arr);       // 1 (index in flattened array)
 * const rowIdx = await argmin(arr, 1); // [1, 0]
 * ```
 */
export async function argmin(
  a: NDArray,
  axis: number | null = null,
  keepdims = false,
): Promise<NDArray | number> {
  const axisVal = axis === null ? -2147483648 : axis;
  const resultPtr = a._wasmModule._ndarray_argmin(
    a._wasmPtr,
    axisVal,
    keepdims,
  );
  if (resultPtr === 0) throw new Error("argmin failed");
  const result = NDArray._fromPtr(resultPtr, a._wasmModule);
  if (axis === null && !keepdims) {
    const v = result.item();
    result.dispose();
    return v;
  }
  return result;
}

/**
 * Find indices where elements should be inserted to maintain order.
 *
 * @param a - Sorted 1D input array
 * @param v - Values to insert
 * @param side - 'left' or 'right' side to insert
 * @returns Indices as number or NDArray
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([1, 2, 4, 5]);
 * const idx = await searchsorted(arr, 3);  // 2
 * ```
 */
export async function searchsorted(
  a: NDArray,
  v: NDArray | number | number[],
  side: "left" | "right" = "left",
): Promise<NDArray | number> {
  let vArr: NDArray;
  let disposeV = false;

  if (typeof v === "number") {
    vArr = await NDArray.fromArray([v]);
    disposeV = true;
  } else if (Array.isArray(v)) {
    vArr = await NDArray.fromArray(v);
    disposeV = true;
  } else {
    vArr = v;
  }

  const sideVal = side === "right" ? 1 : 0;
  const resultPtr = a._wasmModule._ndarray_searchsorted(
    a._wasmPtr,
    vArr._wasmPtr,
    sideVal,
    0,
  );

  if (disposeV) vArr.dispose();
  if (resultPtr === 0) throw new Error("searchsorted failed");

  const result = NDArray._fromPtr(resultPtr, a._wasmModule);
  if (typeof v === "number") {
    const val = result.item();
    result.dispose();
    return val;
  }
  return result;
}

/**
 * Return the cumulative sum of the elements along a given axis.
 *
 * @param a - Input array
 * @param axis - Axis along which to compute cumsum. null = flatten first
 * @param dtype - Type of output array (default: promotes integers to float64)
 * @returns Cumulative sum array with same shape as input (or 1D if axis=null)
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([1, 2, 3, 4]);
 * const cs = await cumsum(arr);  // [1, 3, 6, 10]
 *
 * const arr2d = await NDArray.fromArray([[1, 2], [3, 4]]);
 * const cs0 = await cumsum(arr2d, 0);  // [[1, 2], [4, 6]]
 * const cs1 = await cumsum(arr2d, 1);  // [[1, 3], [3, 7]]
 * ```
 */
export async function cumsum(
  a: NDArray,
  axis: number | null = null,
  dtype?: DType,
): Promise<NDArray> {
  const axisVal = axis === null ? -2147483648 : axis;
  const resultPtr = a._wasmModule._ndarray_cumsum_axis(
    a._wasmPtr,
    axisVal,
    dtype ?? -1,
  );
  if (resultPtr === 0) throw new Error("cumsum failed");
  return NDArray._fromPtr(resultPtr, a._wasmModule);
}

/**
 * Return the cumulative product of elements along a given axis.
 *
 * @param a - Input array
 * @param axis - Axis along which to compute cumprod. null = flatten first
 * @param dtype - Type of output array (default: promotes integers to float64)
 * @returns Cumulative product array with same shape as input (or 1D if axis=null)
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([1, 2, 3, 4]);
 * const cp = await cumprod(arr);  // [1, 2, 6, 24]
 *
 * const arr2d = await NDArray.fromArray([[1, 2], [3, 4]]);
 * const cp0 = await cumprod(arr2d, 0);  // [[1, 2], [3, 8]]
 * const cp1 = await cumprod(arr2d, 1);  // [[1, 2], [3, 12]]
 * ```
 */
export async function cumprod(
  a: NDArray,
  axis: number | null = null,
  dtype?: DType,
): Promise<NDArray> {
  const axisVal = axis === null ? -2147483648 : axis;
  const resultPtr = a._wasmModule._ndarray_cumprod_axis(
    a._wasmPtr,
    axisVal,
    dtype ?? -1,
  );
  if (resultPtr === 0) throw new Error("cumprod failed");
  return NDArray._fromPtr(resultPtr, a._wasmModule);
}

/**
 * Return the cumulative sum of array elements.
 *
 * NumPy 2.0 alias for cumsum with include_initial support.
 *
 * @param x - Input array
 * @param axis - Axis along which to compute. Default: 0
 * @param dtype - Type of output array
 * @param include_initial - If true, prepend 0 to the result (not yet implemented)
 * @returns Cumulative sum array
 */
export async function cumulative_sum(
  x: NDArray,
  axis: number = 0,
  dtype?: DType,
  _include_initial: boolean = false,
): Promise<NDArray> {
  // Note: include_initial is not yet supported
  return cumsum(x, axis, dtype);
}

/**
 * Return the cumulative product of array elements.
 *
 * NumPy 2.0 alias for cumprod with include_initial support.
 *
 * @param x - Input array
 * @param axis - Axis along which to compute. Default: 0
 * @param dtype - Type of output array
 * @param include_initial - If true, prepend 1 to the result (not yet implemented)
 * @returns Cumulative product array
 */
export async function cumulative_prod(
  x: NDArray,
  axis: number = 0,
  dtype?: DType,
  _include_initial: boolean = false,
): Promise<NDArray> {
  // Note: include_initial is not yet supported
  return cumprod(x, axis, dtype);
}

/**
 * Return the cumulative sum of array elements treating NaNs as zero.
 *
 * The cumulative sum does not change when NaNs are encountered and
 * leading NaNs are replaced by zeros.
 *
 * @param a - Input array
 * @param axis - Axis along which to compute. null = flatten first
 * @param dtype - Type of output array
 * @returns Cumulative sum with NaN treated as zero
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([1, NaN, 3, 4]);
 * const ncs = await nancumsum(arr);  // [1, 1, 4, 8]
 * ```
 */
export async function nancumsum(
  a: NDArray,
  axis: number | null = null,
  dtype?: DType,
): Promise<NDArray> {
  const axisVal = axis === null ? -2147483648 : axis;
  const resultPtr = a._wasmModule._ndarray_nancumsum_axis(
    a._wasmPtr,
    axisVal,
    dtype ?? -1,
  );
  if (resultPtr === 0) throw new Error("nancumsum failed");
  return NDArray._fromPtr(resultPtr, a._wasmModule);
}

/**
 * Return the cumulative product of array elements treating NaNs as one.
 *
 * The cumulative product does not change when NaNs are encountered and
 * leading NaNs are replaced by ones.
 *
 * @param a - Input array
 * @param axis - Axis along which to compute. null = flatten first
 * @param dtype - Type of output array
 * @returns Cumulative product with NaN treated as one
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([1, NaN, 3, 4]);
 * const ncp = await nancumprod(arr);  // [1, 1, 3, 12]
 * ```
 */
export async function nancumprod(
  a: NDArray,
  axis: number | null = null,
  dtype?: DType,
): Promise<NDArray> {
  const axisVal = axis === null ? -2147483648 : axis;
  const resultPtr = a._wasmModule._ndarray_nancumprod_axis(
    a._wasmPtr,
    axisVal,
    dtype ?? -1,
  );
  if (resultPtr === 0) throw new Error("nancumprod failed");
  return NDArray._fromPtr(resultPtr, a._wasmModule);
}

/**
 * Return the product of array elements over a given axis.
 *
 * @param a - Input array
 * @param axis - Axis along which to compute product. null computes product of all elements.
 * @param keepdims - If true, reduced axes are left with size 1
 * @param dtype - Output dtype (default: float64 for integers, same as input otherwise)
 * @returns Product as number (if axis=null and keepdims=false) or NDArray
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([1, 2, 3, 4]);
 * const total = await prod(arr);  // 24
 *
 * const arr2d = await NDArray.fromArray([[1, 2], [3, 4]]);
 * const rowProds = await prod(arr2d, 1);  // [2, 12]
 * ```
 */
export async function prod(
  a: NDArray,
  axis: number | null = null,
  keepdims = false,
  dtype?: DType,
): Promise<NDArray | number> {
  const axisVal = axis === null ? -2147483648 : axis;
  const resultPtr = a._wasmModule._ndarray_prod_axis(
    a._wasmPtr,
    axisVal,
    keepdims,
    dtype ?? -1,
  );
  if (resultPtr === 0) throw new Error("prod failed");
  const result = NDArray._fromPtr(resultPtr, a._wasmModule);
  if (axis === null && !keepdims) {
    const v = result.item();
    result.dispose();
    return v;
  }
  return result;
}

/**
 * Range of values (maximum - minimum) along an axis.
 *
 * The name of the function comes from the acronym for 'peak to peak'.
 *
 * @param a - Input array
 * @param axis - Axis along which to find the peak-to-peak value. null computes ptp of flattened array.
 * @param keepdims - If true, reduced axes are left with size 1
 * @returns Peak-to-peak value as number (if axis=null and keepdims=false) or NDArray
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([4, 9, 2, 10, 6]);
 * const range = await ptp(arr);  // 8 (10 - 2)
 *
 * const arr2d = await NDArray.fromArray([[4, 9, 2], [10, 6, 8]]);
 * const rowRanges = await ptp(arr2d, 1);  // [7, 4]
 * ```
 */
export async function ptp(
  a: NDArray,
  axis: number | null = null,
  keepdims = false,
): Promise<NDArray | number> {
  const axisVal = axis === null ? -2147483648 : axis;
  const resultPtr = a._wasmModule._ndarray_ptp_axis(
    a._wasmPtr,
    axisVal,
    keepdims,
  );
  if (resultPtr === 0) throw new Error("ptp failed");
  const result = NDArray._fromPtr(resultPtr, a._wasmModule);
  if (axis === null && !keepdims) {
    const v = result.item();
    result.dispose();
    return v;
  }
  return result;
}

/** Alias for min (matches NumPy's np.amin) */
export const amin = min;

/** Alias for max (matches NumPy's np.amax) */
export const amax = max;

/**
 * Compute the weighted average along the specified axis.
 *
 * @param a - Input array
 * @param axis - Axis along which to compute average. null computes average of all elements.
 * @param weights - Array of weights for each value. If null, all weights are equal to 1.
 * @param keepdims - If true, reduced axes are left with size 1
 * @returns Weighted average as number (if axis=null and keepdims=false) or NDArray
 *
 * @example
 * ```typescript
 * const data = await NDArray.fromArray([1, 2, 3, 4]);
 * const weights = await NDArray.fromArray([4, 3, 2, 1]);
 * const avg = await average(data, null, weights);  // 2.0
 *
 * // Without weights, same as mean
 * const mean_val = await average(data);  // 2.5
 * ```
 */
export async function average(
  a: NDArray,
  axis: number | null = null,
  weights: NDArray | null = null,
  keepdims = false,
): Promise<NDArray | number> {
  // If no weights provided, use mean
  if (weights === null) {
    return mean(a, axis, keepdims);
  }

  // Weighted average: sum(a * weights) / sum(weights)
  const module = a._wasmModule;

  // Multiply array by weights
  const weightedPtr = module._ufunc_multiply(a._wasmPtr, weights._wasmPtr);
  if (weightedPtr === 0) throw new Error("average: multiply failed");
  const weighted = NDArray._fromPtr(weightedPtr, module);

  // Sum of weighted values
  const axisVal = axis === null ? -2147483648 : axis;
  const sumWPtr = module._ndarray_sum_axis(weighted._wasmPtr, axisVal, keepdims, -1);
  if (sumWPtr === 0) {
    weighted.dispose();
    throw new Error("average: sum of weighted values failed");
  }
  const sumW = NDArray._fromPtr(sumWPtr, module);

  // Sum of weights
  const sumWeightsPtr = module._ndarray_sum_axis(weights._wasmPtr, axisVal, keepdims, -1);
  if (sumWeightsPtr === 0) {
    weighted.dispose();
    sumW.dispose();
    throw new Error("average: sum of weights failed");
  }
  const sumWeights = NDArray._fromPtr(sumWeightsPtr, module);

  // Divide to get average
  const resultPtr = module._ufunc_divide(sumW._wasmPtr, sumWeights._wasmPtr);
  weighted.dispose();
  sumW.dispose();
  sumWeights.dispose();

  if (resultPtr === 0) throw new Error("average: divide failed");
  const result = NDArray._fromPtr(resultPtr, module);

  if (axis === null && !keepdims) {
    const v = result.item();
    result.dispose();
    return v;
  }
  return result;
}

/**
 * Compute the q-th percentile of the data along the specified axis.
 *
 * @param a - Input array
 * @param q - Percentile to compute, in range [0, 100]
 * @param axis - Axis along which to compute percentile. null computes percentile of flattened array.
 * @param keepdims - If true, reduced axes are left with size 1
 * @returns Percentile value as number (if axis=null and keepdims=false) or NDArray
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([1, 2, 3, 4, 5, 6, 7, 8, 9, 10]);
 * const p50 = await percentile(arr, 50);  // 5.5 (median)
 * const p25 = await percentile(arr, 25);  // 3.25
 * const p75 = await percentile(arr, 75);  // 7.75
 * ```
 */
export async function percentile(
  a: NDArray,
  q: number,
  axis: number | null = null,
  keepdims = false,
): Promise<NDArray | number> {
  if (q < 0 || q > 100) {
    throw new Error("percentile: q must be in range [0, 100]");
  }
  const axisVal = axis === null ? -2147483648 : axis;
  const resultPtr = a._wasmModule._ndarray_percentile(
    a._wasmPtr,
    q,
    axisVal,
    keepdims,
  );
  if (resultPtr === 0) throw new Error("percentile failed");
  const result = NDArray._fromPtr(resultPtr, a._wasmModule);
  if (axis === null && !keepdims) {
    const v = result.item();
    result.dispose();
    return v;
  }
  return result;
}

/**
 * Compute the q-th quantile of the data along the specified axis.
 *
 * @param a - Input array
 * @param q - Quantile to compute, in range [0, 1]
 * @param axis - Axis along which to compute quantile. null computes quantile of flattened array.
 * @param keepdims - If true, reduced axes are left with size 1
 * @returns Quantile value as number (if axis=null and keepdims=false) or NDArray
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([1, 2, 3, 4, 5, 6, 7, 8, 9, 10]);
 * const q50 = await quantile(arr, 0.5);   // 5.5 (median)
 * const q25 = await quantile(arr, 0.25);  // 3.25
 * const q75 = await quantile(arr, 0.75);  // 7.75
 * ```
 */
export async function quantile(
  a: NDArray,
  q: number,
  axis: number | null = null,
  keepdims = false,
): Promise<NDArray | number> {
  if (q < 0 || q > 1) {
    throw new Error("quantile: q must be in range [0, 1]");
  }
  const axisVal = axis === null ? -2147483648 : axis;
  const resultPtr = a._wasmModule._ndarray_quantile(
    a._wasmPtr,
    q,
    axisVal,
    keepdims,
  );
  if (resultPtr === 0) throw new Error("quantile failed");
  const result = NDArray._fromPtr(resultPtr, a._wasmModule);
  if (axis === null && !keepdims) {
    const v = result.item();
    result.dispose();
    return v;
  }
  return result;
}

/**
 * Estimate a covariance matrix.
 *
 * Covariance indicates the level to which two variables vary together.
 *
 * @param m - A 1-D or 2-D array containing multiple variables and observations
 * @param y - An additional set of variables and observations. Has the same shape as m.
 * @param rowvar - If true, each row represents a variable, with observations in columns.
 *                 If false, the relationship is transposed.
 * @param bias - If false (default), normalization is by (N-1). If true, normalization is by N.
 * @param ddof - If not null, this overrides the bias setting. ddof=1 gives unbiased estimate,
 *               ddof=0 gives biased estimate.
 * @returns The covariance matrix of the variables
 *
 * @example
 * ```typescript
 * // Two variables with 5 observations each
 * const x = await NDArray.fromArray([[0, 2], [1, 1], [2, 0]]);
 * const c = await cov(x, null, false);  // Covariance matrix
 * ```
 */
export async function cov(
  m: NDArray,
  y: NDArray | null = null,
  rowvar = true,
  bias = false,
  ddof: number | null = null,
): Promise<NDArray> {
  const yPtr = y === null ? 0 : y._wasmPtr;
  const effectiveDdof = ddof !== null ? ddof : -1;
  const resultPtr = m._wasmModule._ndarray_cov(
    m._wasmPtr,
    yPtr,
    rowvar,
    bias,
    effectiveDdof,
  );
  if (resultPtr === 0) throw new Error("cov failed");
  return NDArray._fromPtr(resultPtr, m._wasmModule);
}

/**
 * Return Pearson product-moment correlation coefficients.
 *
 * The correlation coefficient matrix of the variables.
 *
 * @param x - A 1-D or 2-D array containing multiple variables and observations
 * @param y - An additional set of variables and observations. Has the same shape as x.
 * @param rowvar - If true, each row represents a variable, with observations in columns.
 *                 If false, the relationship is transposed.
 * @returns The correlation coefficient matrix
 *
 * @example
 * ```typescript
 * // Two variables with correlation
 * const x = await NDArray.fromArray([[0, 1, 2], [2, 1, 0]]);
 * const r = await corrcoef(x);  // Correlation matrix
 * // r[0,0] and r[1,1] are 1.0 (self-correlation)
 * // r[0,1] and r[1,0] are -1.0 (perfect negative correlation)
 * ```
 */
export async function corrcoef(
  x: NDArray,
  y: NDArray | null = null,
  rowvar = true,
): Promise<NDArray> {
  const yPtr = y === null ? 0 : y._wasmPtr;
  const resultPtr = x._wasmModule._ndarray_corrcoef(
    x._wasmPtr,
    yPtr,
    rowvar,
  );
  if (resultPtr === 0) throw new Error("corrcoef failed");
  return NDArray._fromPtr(resultPtr, x._wasmModule);
}
