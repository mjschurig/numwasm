/**
 * NumJS Masked Arrays - Extra Functions
 *
 * Additional statistical and utility functions for masked arrays.
 * Compatible with NumPy's numpy.ma.extras module.
 */

import { NDArray } from '../NDArray.js';
import { DType } from '../types.js';
import { MaskedArray } from './core.js';
import { nomask, MaskedArrayError } from './types.js';

/**
 * Weighted average along an axis.
 *
 * @param a - Input masked array
 * @param axis - Axis along which to average (null for all elements)
 * @param weights - Weights for each element (default: equal weights)
 * @param returned - If true, return tuple of (average, sum_of_weights)
 * @returns Promise resolving to average value(s), optionally with sum of weights
 *
 * @example
 * ```typescript
 * const ma = await MaskedArray.create([1, 2, 3, 4], [false, true, false, false]);
 * await average(ma);  // (1 + 3 + 4) / 3 = 2.667
 * await average(ma, null, [1, 1, 2, 1]);  // weighted average
 * ```
 */
export async function average(
  a: MaskedArray,
  axis: number | null = null,
  weights: NDArray | number[] | null = null,
  returned: boolean = false
): Promise<number | MaskedArray | [number | MaskedArray, number | MaskedArray]> {
  if (axis !== null) {
    throw new MaskedArrayError('average with axis not yet implemented');
  }

  // Get weights array
  let w: NDArray;
  if (weights === null) {
    w = await NDArray.ones([a.size], { dtype: a.dtype });
  } else if (Array.isArray(weights)) {
    w = await NDArray.fromArray(weights, undefined, { dtype: a.dtype });
  } else {
    w = weights;
  }

  // Apply mask to weights
  let sumWeights = 0;
  let sumWx = 0;

  for (let i = 0; i < a.size; i++) {
    if (a._mask !== nomask && (a._mask as NDArray).getFlat(i)) {
      continue;
    }
    const val = a._data.getFlat(i);
    const weight = w.getFlat(i);
    sumWx += val * weight;
    sumWeights += weight;
  }

  const avg = sumWeights === 0 ? NaN : sumWx / sumWeights;

  if (returned) {
    return [avg, sumWeights];
  }

  return avg;
}

/**
 * Compute the median along an axis.
 *
 * @param a - Input masked array
 * @param axis - Axis along which to compute median (null for all elements)
 * @returns Median value(s)
 *
 * @example
 * ```typescript
 * const ma = await MaskedArray.create([1, 2, 3, 4, 5], [false, true, false, false, false]);
 * median(ma);  // median of [1, 3, 4, 5] = 3.5
 * ```
 */
export function median(
  a: MaskedArray,
  axis: number | null = null
): number | MaskedArray {
  if (axis !== null) {
    throw new MaskedArrayError('median with axis not yet implemented');
  }

  // Get non-masked values
  const values: number[] = [];
  for (let i = 0; i < a.size; i++) {
    if (a._mask !== nomask && (a._mask as NDArray).getFlat(i)) {
      continue;
    }
    values.push(a._data.getFlat(i));
  }

  if (values.length === 0) {
    return NaN;
  }

  // Sort and find median
  values.sort((x, y) => x - y);
  const mid = Math.floor(values.length / 2);

  if (values.length % 2 === 0) {
    return (values[mid - 1] + values[mid]) / 2;
  }

  return values[mid];
}

/**
 * Estimate covariance matrix.
 *
 * @param x - Input masked array (1D or 2D)
 * @param y - Optional second variable
 * @param rowvar - If true, each row is a variable (default: true)
 * @param ddof - Delta degrees of freedom (default: 1)
 * @returns Promise resolving to covariance matrix
 *
 * @example
 * ```typescript
 * const x = await MaskedArray.create([[1, 2, 3], [4, 5, 6]]);
 * await cov(x);  // 2x2 covariance matrix
 * ```
 */
export async function cov(
  x: MaskedArray,
  _y: MaskedArray | null = null,
  rowvar: boolean = true,
  ddof: number = 1
): Promise<MaskedArray> {
  // For 1D arrays, treat as single variable
  let data: MaskedArray;
  if (x.ndim === 1) {
    data = x.reshape([1, x.size]);
    rowvar = true;
  } else {
    data = rowvar ? x : x.transpose();
  }

  const nVars = data.shape[0];
  const nObs = data.shape[1];

  // Compute means
  const means: number[] = [];
  for (let i = 0; i < nVars; i++) {
    let sum = 0;
    let count = 0;
    for (let j = 0; j < nObs; j++) {
      const flatIdx = i * nObs + j;
      if (data._mask !== nomask && (data._mask as NDArray).getFlat(flatIdx)) {
        continue;
      }
      sum += data._data.getFlat(flatIdx);
      count++;
    }
    means.push(count > 0 ? sum / count : 0);
  }

  // Compute covariance matrix
  const result = await NDArray.zeros([nVars, nVars], { dtype: DType.Float64 });

  for (let i = 0; i < nVars; i++) {
    for (let j = i; j < nVars; j++) {
      let sumProd = 0;
      let count = 0;

      for (let k = 0; k < nObs; k++) {
        const flatIdxI = i * nObs + k;
        const flatIdxJ = j * nObs + k;

        // Skip if either is masked
        if (data._mask !== nomask) {
          const mask = data._mask as NDArray;
          if (mask.getFlat(flatIdxI) || mask.getFlat(flatIdxJ)) {
            continue;
          }
        }

        const xi = data._data.getFlat(flatIdxI) - means[i];
        const xj = data._data.getFlat(flatIdxJ) - means[j];
        sumProd += xi * xj;
        count++;
      }

      const covValue = count > ddof ? sumProd / (count - ddof) : NaN;
      result.setFlat(i * nVars + j, covValue);
      result.setFlat(j * nVars + i, covValue);
    }
  }

  return new MaskedArray(result);
}

/**
 * Compute correlation coefficients.
 *
 * @param x - Input masked array (1D or 2D)
 * @param y - Optional second variable
 * @param rowvar - If true, each row is a variable (default: true)
 * @returns Promise resolving to correlation matrix
 *
 * @example
 * ```typescript
 * const x = await MaskedArray.create([[1, 2, 3], [4, 5, 6]]);
 * await corrcoef(x);  // 2x2 correlation matrix
 * ```
 */
export async function corrcoef(
  x: MaskedArray,
  y: MaskedArray | null = null,
  rowvar: boolean = true
): Promise<MaskedArray> {
  // Get covariance matrix
  const c = await cov(x, y, rowvar, 1);
  const nVars = c.shape[0];

  // Normalize to correlation
  const result = await NDArray.zeros([nVars, nVars], { dtype: DType.Float64 });

  for (let i = 0; i < nVars; i++) {
    for (let j = 0; j < nVars; j++) {
      const covIJ = c._data.getFlat(i * nVars + j);
      const varI = c._data.getFlat(i * nVars + i);
      const varJ = c._data.getFlat(j * nVars + j);

      const stdI = Math.sqrt(varI);
      const stdJ = Math.sqrt(varJ);

      const corr = stdI > 0 && stdJ > 0 ? covIJ / (stdI * stdJ) : NaN;
      result.setFlat(i * nVars + j, corr);
    }
  }

  return new MaskedArray(result);
}

/**
 * Slice type for contiguous regions.
 */
export interface SliceInfo {
  start: number;
  stop: number;
}

/**
 * Find the indices of the first and last non-masked values.
 *
 * @param a - Input masked array
 * @param axis - Axis along which to search (null for flattened array)
 * @returns Tuple of [first_index, last_index] or null if all masked
 *
 * @example
 * ```typescript
 * const ma = await MaskedArray.create([1, 2, 3, 4], [true, false, false, true]);
 * notmasked_edges(ma);  // [1, 2]
 * ```
 */
export function notmasked_edges(
  a: MaskedArray,
  axis: number | null = null
): [number, number] | null {
  if (axis !== null) {
    throw new MaskedArrayError('notmasked_edges with axis not yet implemented');
  }

  return flatnotmasked_edges(a);
}

/**
 * Find contiguous non-masked regions.
 *
 * @param a - Input masked array
 * @param axis - Axis along which to search (null for flattened array)
 * @returns Array of slice info objects
 *
 * @example
 * ```typescript
 * const ma = await MaskedArray.create([1, 2, 3, 4, 5], [true, false, false, true, false]);
 * notmasked_contiguous(ma);  // [{start: 1, stop: 3}, {start: 4, stop: 5}]
 * ```
 */
export function notmasked_contiguous(
  a: MaskedArray,
  axis: number | null = null
): SliceInfo[] {
  if (axis !== null) {
    throw new MaskedArrayError(
      'notmasked_contiguous with axis not yet implemented'
    );
  }

  return flatnotmasked_contiguous(a);
}

/**
 * Find the flat indices of the first and last non-masked values.
 *
 * @param a - Input masked array
 * @returns Tuple of [first_index, last_index] or null if all masked
 *
 * @example
 * ```typescript
 * const ma = await MaskedArray.create([1, 2, 3, 4], [true, false, false, true]);
 * flatnotmasked_edges(ma);  // [1, 2]
 * ```
 */
export function flatnotmasked_edges(a: MaskedArray): [number, number] | null {
  if (a._mask === nomask) {
    if (a.size === 0) {
      return null;
    }
    return [0, a.size - 1];
  }

  const mask = a._mask as NDArray;
  let first = -1;
  let last = -1;

  for (let i = 0; i < a.size; i++) {
    if (!mask.getFlat(i)) {
      if (first === -1) {
        first = i;
      }
      last = i;
    }
  }

  if (first === -1) {
    return null;
  }

  return [first, last];
}

/**
 * Find contiguous non-masked regions in flattened array.
 *
 * @param a - Input masked array
 * @returns Array of slice info objects
 *
 * @example
 * ```typescript
 * const ma = await MaskedArray.create([1, 2, 3, 4, 5], [true, false, false, true, false]);
 * flatnotmasked_contiguous(ma);  // [{start: 1, stop: 3}, {start: 4, stop: 5}]
 * ```
 */
export function flatnotmasked_contiguous(a: MaskedArray): SliceInfo[] {
  const slices: SliceInfo[] = [];

  if (a._mask === nomask) {
    if (a.size > 0) {
      slices.push({ start: 0, stop: a.size });
    }
    return slices;
  }

  const mask = a._mask as NDArray;
  let inSlice = false;
  let start = 0;

  for (let i = 0; i < a.size; i++) {
    const isMasked = mask.getFlat(i) !== 0;

    if (!isMasked && !inSlice) {
      // Start a new slice
      start = i;
      inSlice = true;
    } else if (isMasked && inSlice) {
      // End the current slice
      slices.push({ start, stop: i });
      inSlice = false;
    }
  }

  // Handle final slice
  if (inSlice) {
    slices.push({ start, stop: a.size });
  }

  return slices;
}

/**
 * Find contiguous masked regions.
 *
 * @param a - Input masked array
 * @returns Array of slice info objects for masked regions
 *
 * @example
 * ```typescript
 * const ma = await MaskedArray.create([1, 2, 3, 4, 5], [true, false, false, true, true]);
 * clump_masked(ma);  // [{start: 0, stop: 1}, {start: 3, stop: 5}]
 * ```
 */
export function clump_masked(a: MaskedArray): SliceInfo[] {
  const slices: SliceInfo[] = [];

  if (a._mask === nomask) {
    return slices;
  }

  const mask = a._mask as NDArray;
  let inSlice = false;
  let start = 0;

  for (let i = 0; i < a.size; i++) {
    const isMasked = mask.getFlat(i) !== 0;

    if (isMasked && !inSlice) {
      // Start a new slice
      start = i;
      inSlice = true;
    } else if (!isMasked && inSlice) {
      // End the current slice
      slices.push({ start, stop: i });
      inSlice = false;
    }
  }

  // Handle final slice
  if (inSlice) {
    slices.push({ start, stop: a.size });
  }

  return slices;
}

/**
 * Find contiguous non-masked regions.
 *
 * @param a - Input masked array
 * @returns Array of slice info objects for non-masked regions
 *
 * @example
 * ```typescript
 * const ma = await MaskedArray.create([1, 2, 3, 4, 5], [true, false, false, true, true]);
 * clump_unmasked(ma);  // [{start: 1, stop: 3}]
 * ```
 */
export function clump_unmasked(a: MaskedArray): SliceInfo[] {
  return flatnotmasked_contiguous(a);
}

/**
 * Apply a function along an axis to masked array.
 *
 * @param func1d - Function to apply (takes 1D array, returns scalar)
 * @param axis - Axis along which to apply
 * @param arr - Input masked array
 * @returns Promise resolving to MaskedArray with results
 */
export async function apply_along_axis(
  func1d: (arr: MaskedArray) => number,
  axis: number,
  arr: MaskedArray
): Promise<MaskedArray> {
  const shape = arr.shape;
  const axisSize = shape[axis];

  // Calculate output shape (remove the axis)
  const outShape = shape.filter((_, i) => i !== axis);

  if (outShape.length === 0) {
    // Single result
    const result = func1d(arr.flatten());
    const resultData = await NDArray.full([1], result, { dtype: arr.dtype });
    return new MaskedArray(resultData);
  }

  // Calculate strides for iteration
  const outSize = outShape.reduce((a, b) => a * b, 1);
  const result = await NDArray.empty(outShape, { dtype: arr.dtype });

  for (let outIdx = 0; outIdx < outSize; outIdx++) {
    // Extract 1D slice along axis
    const indices = flatToMultiIndex(outIdx, outShape);

    // Build slice indices for the original array
    const sliceValues: number[] = [];
    for (let j = 0; j < axisSize; j++) {
      // Insert j at the axis position
      const fullIndices = [
        ...indices.slice(0, axis),
        j,
        ...indices.slice(axis),
      ];
      const flatIdx = multiToFlatIndex(fullIndices, shape);
      sliceValues.push(arr._data.getFlat(flatIdx));
    }

    // Create 1D masked array for this slice
    const sliceMask: number[] = [];
    if (arr._mask !== nomask) {
      const mask = arr._mask as NDArray;
      for (let j = 0; j < axisSize; j++) {
        const fullIndices = [
          ...indices.slice(0, axis),
          j,
          ...indices.slice(axis),
        ];
        const flatIdx = multiToFlatIndex(fullIndices, shape);
        sliceMask.push(mask.getFlat(flatIdx));
      }
    }

    const sliceData = await NDArray.fromArray(sliceValues, undefined, { dtype: arr.dtype });
    const sliceMaskArr = sliceMask.length > 0
      ? await NDArray.fromArray(sliceMask, undefined, { dtype: DType.Bool })
      : nomask;
    const sliceMA = new MaskedArray(sliceData, sliceMaskArr);

    // Apply function
    const value = func1d(sliceMA);
    result.setFlat(outIdx, value);
  }

  return new MaskedArray(result);
}

// ============ Helper Functions ============

function flatToMultiIndex(flatIdx: number, shape: number[]): number[] {
  const indices = new Array(shape.length);
  let remaining = flatIdx;
  for (let i = shape.length - 1; i >= 0; i--) {
    indices[i] = remaining % shape[i];
    remaining = Math.floor(remaining / shape[i]);
  }
  return indices;
}

function multiToFlatIndex(indices: number[], shape: number[]): number {
  let flatIdx = 0;
  let multiplier = 1;
  for (let i = shape.length - 1; i >= 0; i--) {
    flatIdx += indices[i] * multiplier;
    multiplier *= shape[i];
  }
  return flatIdx;
}
