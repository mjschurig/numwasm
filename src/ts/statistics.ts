/**
 * Statistics and searching functions for NDArray
 *
 * Provides WASM-accelerated statistical operations with full NumPy-compatible
 * axis support, keepdims, and dtype parameters.
 */

import { NDArray } from './NDArray.js';
import { DType } from './types.js';

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
  dtype?: DType
): Promise<NDArray | number> {
  const axisVal = axis === null ? -2147483648 : axis;
  const resultPtr = a._wasmModule._ndarray_sum_axis(
    a._wasmPtr,
    axisVal,
    keepdims,
    dtype ?? -1
  );
  if (resultPtr === 0) throw new Error('sum failed');
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
  dtype?: DType
): Promise<NDArray | number> {
  const axisVal = axis === null ? -2147483648 : axis;
  const resultPtr = a._wasmModule._ndarray_mean_axis(
    a._wasmPtr,
    axisVal,
    keepdims,
    dtype ?? -1
  );
  if (resultPtr === 0) throw new Error('mean failed');
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
  dtype?: DType
): Promise<NDArray | number> {
  const axisVal = axis === null ? -2147483648 : axis;
  const resultPtr = a._wasmModule._ndarray_var_axis(
    a._wasmPtr,
    axisVal,
    keepdims,
    dtype ?? -1,
    ddof
  );
  if (resultPtr === 0) throw new Error('variance failed');
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
  dtype?: DType
): Promise<NDArray | number> {
  const axisVal = axis === null ? -2147483648 : axis;
  const resultPtr = a._wasmModule._ndarray_std_axis(
    a._wasmPtr,
    axisVal,
    keepdims,
    dtype ?? -1,
    ddof
  );
  if (resultPtr === 0) throw new Error('std failed');
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
  keepdims = false
): Promise<NDArray | number> {
  const axisVal = axis === null ? -2147483648 : axis;
  const resultPtr = a._wasmModule._ndarray_min_axis(a._wasmPtr, axisVal, keepdims);
  if (resultPtr === 0) throw new Error('min failed');
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
  keepdims = false
): Promise<NDArray | number> {
  const axisVal = axis === null ? -2147483648 : axis;
  const resultPtr = a._wasmModule._ndarray_max_axis(a._wasmPtr, axisVal, keepdims);
  if (resultPtr === 0) throw new Error('max failed');
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
  keepdims = false
): Promise<NDArray | number> {
  const axisVal = axis === null ? -2147483648 : axis;
  const resultPtr = a._wasmModule._ndarray_median(a._wasmPtr, axisVal, keepdims);
  if (resultPtr === 0) throw new Error('median failed');
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
  keepdims = false
): Promise<NDArray | number> {
  const axisVal = axis === null ? -2147483648 : axis;
  const resultPtr = a._wasmModule._ndarray_argmax(a._wasmPtr, axisVal, keepdims);
  if (resultPtr === 0) throw new Error('argmax failed');
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
  keepdims = false
): Promise<NDArray | number> {
  const axisVal = axis === null ? -2147483648 : axis;
  const resultPtr = a._wasmModule._ndarray_argmin(a._wasmPtr, axisVal, keepdims);
  if (resultPtr === 0) throw new Error('argmin failed');
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
  side: 'left' | 'right' = 'left'
): Promise<NDArray | number> {
  let vArr: NDArray;
  let disposeV = false;

  if (typeof v === 'number') {
    vArr = await NDArray.fromArray([v]);
    disposeV = true;
  } else if (Array.isArray(v)) {
    vArr = await NDArray.fromArray(v);
    disposeV = true;
  } else {
    vArr = v;
  }

  const sideVal = side === 'right' ? 1 : 0;
  const resultPtr = a._wasmModule._ndarray_searchsorted(
    a._wasmPtr,
    vArr._wasmPtr,
    sideVal,
    0
  );

  if (disposeV) vArr.dispose();
  if (resultPtr === 0) throw new Error('searchsorted failed');

  const result = NDArray._fromPtr(resultPtr, a._wasmModule);
  if (typeof v === 'number') {
    const val = result.item();
    result.dispose();
    return val;
  }
  return result;
}
