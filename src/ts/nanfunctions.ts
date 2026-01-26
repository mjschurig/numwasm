/**
 * NumJS NaN-Handling Functions
 *
 * Provides NumPy-compatible NaN-aware statistical and aggregation functions.
 * These functions ignore NaN values when computing statistics.
 *
 * Reference: numpy/lib/_nanfunctions_impl.py
 */

import { NDArray } from './NDArray.js';
import { DType } from './types.js';
import { isnan, isposinf, isneginf, all, any } from './logic.js';
import { where, extract } from './indexing.js';
import { sum, min, max, argmin, argmax, median } from './statistics.js';
import { divide, subtract, multiply, sqrt, equal } from './ufunc.js';
import { sort } from './sorting.js';
import { applyAlongAxis } from './functional.js';
import { finfo } from './typeinfo.js';

/* ============ Internal Helpers ============ */

/**
 * Internal: Get a mask where non-NaN values are 1 and NaN values are 0.
 * Uses the property that NaN != NaN, so equal(arr, arr) returns 0 for NaN.
 * Returns the mask (caller must dispose).
 */
function _getNotNanMask(arr: NDArray): NDArray {
  // equal(arr, arr) returns 1 for non-NaN values and 0 for NaN
  // because NaN != NaN by IEEE 754
  return equal(arr, arr);
}

/**
 * Internal: Replace NaN values with a replacement value.
 * Returns the cleaned array (caller must dispose).
 */
async function _replaceNan(arr: NDArray, replacement: number): Promise<NDArray> {
  const mask = await isnan(arr);
  const replacementArr = await NDArray.fromArray([replacement]);
  const result = await where(mask, replacementArr, arr);
  mask.dispose();
  replacementArr.dispose();
  return result;
}

/**
 * Internal: Count non-NaN values along axis.
 */
async function _countNonNan(
  arr: NDArray,
  axis: number | null,
  keepdims: boolean
): Promise<NDArray | number> {
  const notNanMask = _getNotNanMask(arr);
  const result = await sum(notNanMask, axis, keepdims, DType.Float64);
  notNanMask.dispose();
  return result;
}

/**
 * Internal: Check for all-NaN slices and warn.
 */
async function _warnAllNan(arr: NDArray, axis: number | null, funcName: string): Promise<boolean> {
  const mask = await isnan(arr);
  const allNan = await all(mask, axis === null ? undefined : axis);
  mask.dispose();

  let hasAllNanSlice = false;
  if (typeof allNan === 'boolean') {
    hasAllNanSlice = allNan;
  } else {
    const anyAllNan = await any(allNan);
    hasAllNanSlice = anyAllNan as boolean;
    allNan.dispose();
  }

  if (hasAllNanSlice) {
    console.warn(`All-NaN slice encountered in ${funcName}`);
  }
  return hasAllNanSlice;
}

/* ============ 23.1 Basic NaN Aggregations ============ */

/**
 * Return the sum of array elements over a given axis treating NaNs as zero.
 *
 * @param a - Array containing numbers whose sum is desired
 * @param axis - Axis or axes along which the sum is computed (null = all elements)
 * @param keepdims - If true, reduced axes are left with size one
 * @param dtype - The type of the returned array
 * @returns Sum of array elements, with NaN treated as zero
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([1, NaN, 3]);
 * await nansum(arr);  // 4
 *
 * const mat = await NDArray.fromArray([[1, NaN], [3, 4]]);
 * await nansum(mat, 1);  // [1, 7]
 * ```
 */
export async function nansum(
  a: NDArray,
  axis: number | null = null,
  keepdims: boolean = false,
  dtype?: DType
): Promise<NDArray | number> {
  const cleaned = await _replaceNan(a, 0);
  const result = await sum(cleaned, axis, keepdims, dtype);
  cleaned.dispose();
  return result;
}

/**
 * Return the product of array elements over a given axis treating NaNs as one.
 *
 * @param a - Array containing numbers whose product is desired
 * @param axis - Axis or axes along which the product is computed (null = all elements)
 * @param keepdims - If true, reduced axes are left with size one
 * @param dtype - The type of the returned array
 * @returns Product of array elements, with NaN treated as one
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([1, NaN, 3]);
 * await nanprod(arr);  // 3
 *
 * const mat = await NDArray.fromArray([[1, NaN], [3, 4]]);
 * await nanprod(mat, 1);  // [1, 12]
 * ```
 */
export async function nanprod(
  a: NDArray,
  axis: number | null = null,
  keepdims: boolean = false,
  _dtype?: DType
): Promise<NDArray | number> {
  // Replace NaN with 1 for product
  const cleaned = await _replaceNan(a, 1);

  if (axis === null) {
    // Flatten and compute product
    const flat = cleaned.ravel();
    const data = flat.toArray() as number[];
    flat.dispose();
    cleaned.dispose();

    const product = data.reduce((acc, val) => acc * val, 1);

    if (keepdims) {
      const shape = new Array(a.ndim).fill(1);
      return NDArray.full(shape, product);
    }
    return product;
  }

  // Normalize negative axis
  const normalizedAxis = axis < 0 ? axis + a.ndim : axis;

  // Use applyAlongAxis for axis-wise product
  const result = await applyAlongAxis(
    (slice: NDArray) => {
      const data = slice.toArray() as number[];
      return data.reduce((acc, val) => acc * val, 1);
    },
    normalizedAxis,
    cleaned
  );

  cleaned.dispose();

  if (keepdims) {
    return result.expandDims(normalizedAxis);
  }
  return result;
}

/**
 * Return minimum of an array or minimum along an axis, ignoring any NaNs.
 *
 * @param a - Array containing numbers
 * @param axis - Axis or axes along which to operate (null = all elements)
 * @param keepdims - If true, reduced axes are left with size one
 * @returns Minimum values with NaN ignored
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([1, NaN, 3]);
 * await nanmin(arr);  // 1
 *
 * const mat = await NDArray.fromArray([[1, NaN], [NaN, 4]]);
 * await nanmin(mat, 1);  // [1, 4]
 * ```
 *
 * @throws Warning if all values along axis are NaN (returns NaN)
 */
export async function nanmin(
  a: NDArray,
  axis: number | null = null,
  keepdims: boolean = false
): Promise<NDArray | number> {
  await _warnAllNan(a, axis, 'nanmin');

  // Replace NaN with +Infinity for min computation
  const cleaned = await _replaceNan(a, Infinity);
  const result = await min(cleaned, axis, keepdims);
  cleaned.dispose();

  // Handle all-NaN case: Infinity -> NaN
  if (typeof result === 'number') {
    return result === Infinity ? NaN : result;
  }

  // For array result, replace Infinity with NaN
  const infArr = await NDArray.fromArray([Infinity]);
  const eqInf = equal(result, infArr);
  const nanArr = await NDArray.fromArray([NaN]);
  const finalResult = await where(eqInf, nanArr, result);

  infArr.dispose();
  eqInf.dispose();
  nanArr.dispose();
  result.dispose();

  return finalResult;
}

/**
 * Return maximum of an array or maximum along an axis, ignoring any NaNs.
 *
 * @param a - Array containing numbers
 * @param axis - Axis or axes along which to operate (null = all elements)
 * @param keepdims - If true, reduced axes are left with size one
 * @returns Maximum values with NaN ignored
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([1, NaN, 3]);
 * await nanmax(arr);  // 3
 *
 * const mat = await NDArray.fromArray([[1, NaN], [NaN, 4]]);
 * await nanmax(mat, 1);  // [1, 4]
 * ```
 *
 * @throws Warning if all values along axis are NaN (returns NaN)
 */
export async function nanmax(
  a: NDArray,
  axis: number | null = null,
  keepdims: boolean = false
): Promise<NDArray | number> {
  await _warnAllNan(a, axis, 'nanmax');

  // Replace NaN with -Infinity for max computation
  const cleaned = await _replaceNan(a, -Infinity);
  const result = await max(cleaned, axis, keepdims);
  cleaned.dispose();

  if (typeof result === 'number') {
    return result === -Infinity ? NaN : result;
  }

  const negInfArr = await NDArray.fromArray([-Infinity]);
  const eqNegInf = equal(result, negInfArr);
  const nanArr = await NDArray.fromArray([NaN]);
  const finalResult = await where(eqNegInf, nanArr, result);

  negInfArr.dispose();
  eqNegInf.dispose();
  nanArr.dispose();
  result.dispose();

  return finalResult;
}

/**
 * Return the indices of the minimum values along an axis, ignoring NaNs.
 *
 * @param a - Input array
 * @param axis - Axis along which to operate (null = flattened array)
 * @param keepdims - If true, reduced axes are left with size one
 * @returns Array of indices into the array
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([NaN, 2, 3]);
 * await nanargmin(arr);  // 1
 * ```
 */
export async function nanargmin(
  a: NDArray,
  axis: number | null = null,
  keepdims: boolean = false
): Promise<NDArray | number> {
  await _warnAllNan(a, axis, 'nanargmin');

  const cleaned = await _replaceNan(a, Infinity);
  const result = await argmin(cleaned, axis, keepdims);
  cleaned.dispose();

  return result;
}

/**
 * Return the indices of the maximum values along an axis, ignoring NaNs.
 *
 * @param a - Input array
 * @param axis - Axis along which to operate (null = flattened array)
 * @param keepdims - If true, reduced axes are left with size one
 * @returns Array of indices into the array
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([NaN, 2, 3]);
 * await nanargmax(arr);  // 2
 * ```
 */
export async function nanargmax(
  a: NDArray,
  axis: number | null = null,
  keepdims: boolean = false
): Promise<NDArray | number> {
  await _warnAllNan(a, axis, 'nanargmax');

  const cleaned = await _replaceNan(a, -Infinity);
  const result = await argmax(cleaned, axis, keepdims);
  cleaned.dispose();

  return result;
}

/* ============ 23.2 NaN Statistics ============ */

/**
 * Compute the arithmetic mean along the specified axis, ignoring NaNs.
 *
 * @param a - Array containing numbers whose mean is desired
 * @param axis - Axis or axes along which the means are computed (null = all elements)
 * @param keepdims - If true, reduced axes are left with size one
 * @param dtype - Type to use in computing the mean
 * @returns Mean of non-NaN elements
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([1, NaN, 3]);
 * await nanmean(arr);  // 2
 *
 * const mat = await NDArray.fromArray([[1, NaN], [3, 4]]);
 * await nanmean(mat, 1);  // [1, 3.5]
 * ```
 */
export async function nanmean(
  a: NDArray,
  axis: number | null = null,
  keepdims: boolean = false,
  dtype?: DType
): Promise<NDArray | number> {
  // Sum of non-NaN values
  const nanSum = await nansum(a, axis, keepdims, dtype);

  // Count of non-NaN values
  const count = await _countNonNan(a, axis, keepdims);

  // Mean = sum / count
  if (typeof nanSum === 'number' && typeof count === 'number') {
    return count === 0 ? NaN : nanSum / count;
  }

  const sumArr = typeof nanSum === 'number'
    ? await NDArray.fromArray([nanSum])
    : nanSum;
  const countArr = typeof count === 'number'
    ? await NDArray.fromArray([count])
    : count;

  const result = divide(sumArr, countArr);

  if (typeof nanSum === 'number') sumArr.dispose();
  if (typeof count === 'number') countArr.dispose();
  if (typeof nanSum !== 'number') (nanSum as NDArray).dispose();
  if (typeof count !== 'number') (count as NDArray).dispose();

  return result;
}

/**
 * Compute the variance along the specified axis, while ignoring NaNs.
 *
 * @param a - Array containing numbers
 * @param axis - Axis or axes along which the variance is computed (null = all elements)
 * @param ddof - Delta Degrees of Freedom (divisor is N - ddof)
 * @param keepdims - If true, reduced axes are left with size one
 * @param dtype - Type to use in computing the variance
 * @returns Variance of non-NaN elements
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([1, NaN, 3]);
 * await nanvar(arr);  // 1.0
 * ```
 */
export async function nanvar(
  a: NDArray,
  axis: number | null = null,
  ddof: number = 0,
  keepdims: boolean = false,
  dtype?: DType
): Promise<NDArray | number> {
  // Compute mean (with keepdims=true for broadcasting)
  const meanVal = await nanmean(a, axis, true, dtype);

  // Subtract mean from array
  const meanArr = typeof meanVal === 'number'
    ? await NDArray.full(a.shape, meanVal)
    : meanVal;

  const centered = subtract(a, meanArr);
  if (typeof meanVal === 'number') meanArr.dispose();

  // Square the centered values
  const squared = multiply(centered, centered);
  centered.dispose();

  // Replace NaN with 0 in squared values (for positions that were originally NaN)
  const mask = await isnan(a);
  const zeroArr = await NDArray.fromArray([0]);
  const squaredClean = await where(mask, zeroArr, squared);
  mask.dispose();
  zeroArr.dispose();
  squared.dispose();

  // Sum of squared deviations
  const sumSq = await sum(squaredClean, axis, keepdims, dtype);
  squaredClean.dispose();

  // Count non-NaN - ddof
  const count = await _countNonNan(a, axis, keepdims);

  // Variance = sumSq / (count - ddof)
  if (typeof sumSq === 'number' && typeof count === 'number') {
    const divisor = count - ddof;
    if (typeof meanVal !== 'number') meanVal.dispose();
    return divisor <= 0 ? NaN : sumSq / divisor;
  }

  const sumArr = typeof sumSq === 'number' ? await NDArray.fromArray([sumSq]) : sumSq;
  const countArr = typeof count === 'number' ? await NDArray.fromArray([count]) : count;
  const ddofArr = await NDArray.fromArray([ddof]);
  const divisorArr = subtract(countArr, ddofArr);
  ddofArr.dispose();

  const result = divide(sumArr, divisorArr);

  divisorArr.dispose();
  if (typeof sumSq === 'number') sumArr.dispose();
  if (typeof count === 'number') countArr.dispose();
  if (typeof sumSq !== 'number') (sumSq as NDArray).dispose();
  if (typeof count !== 'number') (count as NDArray).dispose();
  if (typeof meanVal !== 'number') meanVal.dispose();

  return result;
}

/**
 * Compute the standard deviation along the specified axis, ignoring NaNs.
 *
 * @param a - Array containing numbers
 * @param axis - Axis or axes along which the standard deviation is computed (null = all elements)
 * @param ddof - Delta Degrees of Freedom (divisor is N - ddof)
 * @param keepdims - If true, reduced axes are left with size one
 * @param dtype - Type to use in computing the standard deviation
 * @returns Standard deviation of non-NaN elements
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([1, NaN, 3]);
 * await nanstd(arr);  // 1.0
 * ```
 */
export async function nanstd(
  a: NDArray,
  axis: number | null = null,
  ddof: number = 0,
  keepdims: boolean = false,
  dtype?: DType
): Promise<NDArray | number> {
  const variance = await nanvar(a, axis, ddof, keepdims, dtype);

  if (typeof variance === 'number') {
    return Math.sqrt(variance);
  }

  const result = sqrt(variance);
  variance.dispose();
  return result;
}

/**
 * Compute the median along the specified axis, while ignoring NaNs.
 *
 * @param a - Input array
 * @param axis - Axis or axes along which the medians are computed (null = all elements)
 * @param keepdims - If true, reduced axes are left with size one
 * @returns Median of non-NaN elements
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([1, NaN, 3, 4]);
 * await nanmedian(arr);  // 3.0
 * ```
 */
export async function nanmedian(
  a: NDArray,
  axis: number | null = null,
  keepdims: boolean = false
): Promise<NDArray | number> {
  if (axis === null) {
    // Flatten, extract non-NaN, compute median
    const flat = a.ravel();
    const notNanMask = _getNotNanMask(flat);

    const validValues = await extract(notNanMask, flat);
    notNanMask.dispose();
    flat.dispose();

    if (validValues.size === 0) {
      validValues.dispose();
      if (keepdims) {
        const shape = new Array(a.ndim).fill(1);
        return NDArray.full(shape, NaN);
      }
      return NaN;
    }

    const result = await median(validValues);
    validValues.dispose();

    if (keepdims) {
      const scalar = typeof result === 'number' ? result : result.item();
      if (typeof result !== 'number') result.dispose();
      const shape = new Array(a.ndim).fill(1);
      return NDArray.full(shape, scalar);
    }
    return result;
  }

  // Normalize axis
  const normalizedAxis = axis < 0 ? axis + a.ndim : axis;

  // Axis-wise median using applyAlongAxis
  const result = await applyAlongAxis(
    async (slice: NDArray) => {
      const notNanMask = _getNotNanMask(slice);

      const validValues = await extract(notNanMask, slice);
      notNanMask.dispose();

      if (validValues.size === 0) {
        validValues.dispose();
        return NaN;
      }

      const med = await median(validValues);
      validValues.dispose();

      if (typeof med === 'number') return med;
      const val = med.item();
      med.dispose();
      return val;
    },
    normalizedAxis,
    a
  );

  if (keepdims) {
    return result.expandDims(normalizedAxis);
  }

  return result;
}

/* ============ 23.3 NaN Quantiles ============ */

/**
 * Internal: Compute quantile for a 1D sorted array.
 */
async function _computeQuantile(
  arr: NDArray,
  quantiles: number[],
  interpolation: string
): Promise<NDArray | number> {
  const sorted = await sort(arr, null);
  const n = sorted.size;
  const results: number[] = [];

  for (const q of quantiles) {
    const virtualIdx = q * (n - 1);
    const lower = Math.floor(virtualIdx);
    const upper = Math.ceil(virtualIdx);
    const frac = virtualIdx - lower;

    const lowerVal = sorted.getFlat(Math.min(lower, n - 1));
    const upperVal = sorted.getFlat(Math.min(upper, n - 1));

    let value: number;
    switch (interpolation) {
      case 'lower':
        value = lowerVal;
        break;
      case 'higher':
        value = upperVal;
        break;
      case 'nearest':
        value = frac < 0.5 ? lowerVal : upperVal;
        break;
      case 'midpoint':
        value = (lowerVal + upperVal) / 2;
        break;
      case 'linear':
      default:
        value = lowerVal * (1 - frac) + upperVal * frac;
        break;
    }
    results.push(value);
  }

  sorted.dispose();

  if (results.length === 1) {
    return results[0];
  }
  return NDArray.fromArray(results);
}

/**
 * Compute the q-th quantile of the data along the specified axis,
 * while ignoring NaN values.
 *
 * @param a - Input array
 * @param q - Quantile(s) to compute, in range [0, 1]
 * @param axis - Axis along which to compute (null = all elements)
 * @param interpolation - Interpolation method ('linear', 'lower', 'higher', 'midpoint', 'nearest')
 * @param keepdims - If true, reduced axes are left with size one
 * @returns Quantile values
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([1, 2, NaN, 4]);
 * await nanquantile(arr, 0.5);  // 2.0
 * ```
 */
export async function nanquantile(
  a: NDArray,
  q: number | number[],
  axis: number | null = null,
  interpolation: 'linear' | 'lower' | 'higher' | 'midpoint' | 'nearest' = 'linear',
  keepdims: boolean = false
): Promise<NDArray | number> {
  const qArr = Array.isArray(q) ? q : [q];

  // Validate q values
  for (const qVal of qArr) {
    if (qVal < 0 || qVal > 1) {
      throw new Error(`Quantiles must be in the range [0, 1], got ${qVal}`);
    }
  }

  if (axis === null) {
    // Flatten, extract valid, compute quantile
    const flat = a.ravel();
    const notNanMask = _getNotNanMask(flat);

    const validValues = await extract(notNanMask, flat);
    notNanMask.dispose();
    flat.dispose();

    if (validValues.size === 0) {
      validValues.dispose();
      if (qArr.length === 1) {
        if (keepdims) {
          const shape = new Array(a.ndim).fill(1);
          return NDArray.full(shape, NaN);
        }
        return NaN;
      }
      return NDArray.full([qArr.length], NaN);
    }

    const result = await _computeQuantile(validValues, qArr, interpolation);
    validValues.dispose();

    if (keepdims && qArr.length === 1 && typeof result === 'number') {
      const shape = new Array(a.ndim).fill(1);
      return NDArray.full(shape, result);
    }
    return result;
  }

  // Normalize axis
  const normalizedAxis = axis < 0 ? axis + a.ndim : axis;

  // Axis-wise quantile
  const result = await applyAlongAxis(
    async (slice: NDArray) => {
      const notNanMask = _getNotNanMask(slice);

      const validValues = await extract(notNanMask, slice);
      notNanMask.dispose();

      if (validValues.size === 0) {
        validValues.dispose();
        if (qArr.length === 1) {
          return NaN;
        }
        return NDArray.full([qArr.length], NaN);
      }

      const quantileResult = await _computeQuantile(validValues, qArr, interpolation);
      validValues.dispose();
      return quantileResult;
    },
    normalizedAxis,
    a
  );

  if (keepdims) {
    return result.expandDims(normalizedAxis);
  }

  return result;
}

/**
 * Compute the q-th percentile of the data along the specified axis,
 * while ignoring NaN values.
 *
 * @param a - Input array
 * @param q - Percentile(s) to compute, in range [0, 100]
 * @param axis - Axis along which to compute (null = all elements)
 * @param interpolation - Interpolation method ('linear', 'lower', 'higher', 'midpoint', 'nearest')
 * @param keepdims - If true, reduced axes are left with size one
 * @returns Percentile values
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([1, 2, NaN, 4]);
 * await nanpercentile(arr, 50);  // 2.0
 *
 * const mat = await NDArray.fromArray([[1, NaN], [3, 4]]);
 * await nanpercentile(mat, 50, 1);  // [1.0, 3.5]
 * ```
 */
export async function nanpercentile(
  a: NDArray,
  q: number | number[],
  axis: number | null = null,
  interpolation: 'linear' | 'lower' | 'higher' | 'midpoint' | 'nearest' = 'linear',
  keepdims: boolean = false
): Promise<NDArray | number> {
  // Convert percentile to quantile
  const qArr = Array.isArray(q) ? q : [q];

  // Validate percentile values
  for (const pVal of qArr) {
    if (pVal < 0 || pVal > 100) {
      throw new Error(`Percentiles must be in the range [0, 100], got ${pVal}`);
    }
  }

  const quantiles = qArr.map(p => p / 100);

  return nanquantile(a, quantiles.length === 1 ? quantiles[0] : quantiles, axis, interpolation, keepdims);
}

/* ============ 23.4 NaN Utilities ============ */

/**
 * Replace NaN with zero and infinity with large finite numbers.
 *
 * @param x - Input array
 * @param nan - Value to replace NaN with (default: 0.0)
 * @param posinf - Value to replace positive infinity with (default: largest finite float)
 * @param neginf - Value to replace negative infinity with (default: most negative finite float)
 * @returns Array with NaN and infinity replaced
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([NaN, Infinity, -Infinity, 1]);
 * await nan_to_num(arr);
 * // [0, 1.7976931348623157e+308, -1.7976931348623157e+308, 1]
 *
 * await nan_to_num(arr, -1, 999, -999);
 * // [-1, 999, -999, 1]
 * ```
 */
export async function nan_to_num(
  x: NDArray,
  nan: number = 0.0,
  posinf: number | null = null,
  neginf: number | null = null
): Promise<NDArray> {
  // Get dtype-appropriate infinity replacements
  let posinfVal: number;
  let neginfVal: number;

  if (posinf !== null) {
    posinfVal = posinf;
  } else {
    // Use float64 max as default
    const info = new finfo(DType.Float64);
    posinfVal = info.max;
  }

  if (neginf !== null) {
    neginfVal = neginf;
  } else {
    // Use float64 min (most negative) as default
    const info = new finfo(DType.Float64);
    neginfVal = -info.max;
  }

  // Replace NaN
  const nanMask = await isnan(x);
  const nanArr = await NDArray.fromArray([nan]);
  let result = await where(nanMask, nanArr, x);
  nanMask.dispose();
  nanArr.dispose();

  // Replace positive infinity
  const posinfMask = await isposinf(result);
  const posinfArr = await NDArray.fromArray([posinfVal]);
  const result2 = await where(posinfMask, posinfArr, result);
  posinfMask.dispose();
  posinfArr.dispose();
  result.dispose();
  result = result2;

  // Replace negative infinity
  const neginfMask = await isneginf(result);
  const neginfArr = await NDArray.fromArray([neginfVal]);
  const finalResult = await where(neginfMask, neginfArr, result);
  neginfMask.dispose();
  neginfArr.dispose();
  result.dispose();

  return finalResult;
}
