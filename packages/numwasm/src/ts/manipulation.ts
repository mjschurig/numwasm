/**
 * Array manipulation functions for NumJS-WASM
 *
 * Phase 5: Joining, splitting, tiling, and rearranging arrays.
 * Based on NumPy's array manipulation routines.
 */

import { NDArray } from './NDArray.js';
import { DType } from './types.js';
import { atleast_1d, atleast_2d, atleast_3d, take } from './indexing.js';
import { Slice, slice } from './slice.js';

/* ============ Helper Functions ============ */

function arraysEqual(a: number[], b: number[]): boolean {
  if (a.length !== b.length) return false;
  for (let i = 0; i < a.length; i++) {
    if (a[i] !== b[i]) return false;
  }
  return true;
}

function flatToMulti(flatIdx: number, shape: number[]): number[] {
  const result = new Array(shape.length);
  let remainder = flatIdx;
  for (let i = shape.length - 1; i >= 0; i--) {
    result[i] = remainder % shape[i];
    remainder = Math.floor(remainder / shape[i]);
  }
  return result;
}

function multiToFlat(multiIdx: number[], shape: number[]): number {
  let flat = 0;
  let multiplier = 1;
  for (let i = shape.length - 1; i >= 0; i--) {
    flat += multiIdx[i] * multiplier;
    multiplier *= shape[i];
  }
  return flat;
}

function copyArrayData(src: NDArray, dst: NDArray): void {
  for (let i = 0; i < src.size; i++) {
    dst.setFlat(i, src.getFlat(i));
  }
}

/* ============ Joining Functions ============ */

/**
 * Join a sequence of arrays along an existing axis.
 *
 * @param arrays - Sequence of arrays to concatenate
 * @param axis - The axis along which to concatenate (default: 0)
 * @param dtype - Optional dtype for the output array
 * @returns Concatenated array
 *
 * @example
 * ```typescript
 * const a = await NDArray.fromArray([[1, 2], [3, 4]]);
 * const b = await NDArray.fromArray([[5, 6]]);
 * const result = await concatenate([a, b], 0);
 * // [[1, 2], [3, 4], [5, 6]]
 * ```
 */
export function concatenate(
  arrays: NDArray[],
  axis: number = 0,
  dtype?: DType
): NDArray {
  if (arrays.length === 0) {
    throw new Error('need at least one array to concatenate');
  }

  const module = arrays[0]._wasmModule;

  // Allocate array of pointers
  const ptrsPtr = module._malloc(arrays.length * 4);
  for (let i = 0; i < arrays.length; i++) {
    module.setValue(ptrsPtr + i * 4, arrays[i]._wasmPtr, 'i32');
  }

  const resultPtr = module._ndarray_concatenate(ptrsPtr, arrays.length, axis);
  module._free(ptrsPtr);

  if (resultPtr === 0) {
    throw new Error('concatenate failed: incompatible shapes or invalid axis');
  }

  let result = NDArray._fromPtr(resultPtr, module);

  if (dtype !== undefined && result.dtype !== dtype) {
    const converted = result.astype(dtype);
    result.dispose();
    result = converted;
  }

  return result;
}

/**
 * Join a sequence of arrays along a new axis.
 *
 * @param arrays - Sequence of arrays with the same shape
 * @param axis - The axis in the result along which to stack (default: 0)
 * @returns Stacked array with one more dimension
 *
 * @example
 * ```typescript
 * const a = await NDArray.fromArray([1, 2, 3]);
 * const b = await NDArray.fromArray([4, 5, 6]);
 * const result = await stack([a, b], 0);
 * // [[1, 2, 3], [4, 5, 6]]
 * ```
 */
export function stack(arrays: NDArray[], axis: number = 0): NDArray {
  if (arrays.length === 0) {
    throw new Error('need at least one array to stack');
  }

  // Validate all arrays have same shape
  const shape = arrays[0].shape;
  for (let i = 1; i < arrays.length; i++) {
    if (!arraysEqual(arrays[i].shape, shape)) {
      throw new Error('all input arrays must have the same shape');
    }
  }

  // Normalize axis
  const ndim = arrays[0].ndim + 1;
  axis = axis < 0 ? axis + ndim : axis;
  if (axis < 0 || axis > arrays[0].ndim) {
    throw new Error(`axis ${axis} out of bounds for array with ${ndim} dimensions`);
  }

  // Expand each array along the new axis
  const expanded = arrays.map(arr => arr.expandDims(axis));

  // Concatenate along the new axis
  const result = concatenate(expanded, axis);

  // Clean up expanded views
  for (const exp of expanded) {
    exp.dispose();
  }

  return result;
}

/**
 * Stack arrays in sequence vertically (row-wise).
 * Equivalent to concatenation along the first axis after making 1-D arrays 2-D.
 *
 * @param tup - Sequence of arrays
 * @returns Stacked array
 *
 * @example
 * ```typescript
 * const a = await NDArray.fromArray([1, 2, 3]);
 * const b = await NDArray.fromArray([4, 5, 6]);
 * const result = await vstack([a, b]);
 * // [[1, 2, 3], [4, 5, 6]]
 * ```
 */
export function vstack(tup: NDArray[]): NDArray {
  const arrs = atleast_2d(...tup);
  const arrays = Array.isArray(arrs) ? arrs : [arrs];
  return concatenate(arrays, 0);
}

// Alias for vstack
export const row_stack = vstack;

/**
 * Stack arrays in sequence horizontally (column-wise).
 * For 1-D arrays, concatenates along axis 0.
 * For 2-D+ arrays, concatenates along axis 1.
 *
 * @param tup - Sequence of arrays
 * @returns Stacked array
 *
 * @example
 * ```typescript
 * const a = await NDArray.fromArray([1, 2, 3]);
 * const b = await NDArray.fromArray([4, 5, 6]);
 * const result = await hstack([a, b]);
 * // [1, 2, 3, 4, 5, 6]
 * ```
 */
export function hstack(tup: NDArray[]): NDArray {
  const arrs = atleast_1d(...tup);
  const arrays = Array.isArray(arrs) ? arrs : [arrs];

  // For 1D arrays, concatenate along axis 0
  // For 2D+ arrays, concatenate along axis 1
  if (arrays[0].ndim === 1) {
    return concatenate(arrays, 0);
  }
  return concatenate(arrays, 1);
}

/**
 * Stack arrays in sequence along the third axis (depth-wise).
 *
 * @param tup - Sequence of arrays
 * @returns Stacked array with at least 3 dimensions
 *
 * @example
 * ```typescript
 * const a = await NDArray.fromArray([[1, 2], [3, 4]]);
 * const b = await NDArray.fromArray([[5, 6], [7, 8]]);
 * const result = await dstack([a, b]);
 * // shape: [2, 2, 2]
 * ```
 */
export function dstack(tup: NDArray[]): NDArray {
  const arrs = atleast_3d(...tup);
  const arrays = Array.isArray(arrs) ? arrs : [arrs];
  return concatenate(arrays, 2);
}

/**
 * Stack 1-D arrays as columns into a 2-D array.
 *
 * @param tup - Sequence of 1-D or 2-D arrays
 * @returns 2-D array with input arrays as columns
 *
 * @example
 * ```typescript
 * const a = await NDArray.fromArray([1, 2, 3]);
 * const b = await NDArray.fromArray([4, 5, 6]);
 * const result = await column_stack([a, b]);
 * // [[1, 4], [2, 5], [3, 6]]
 * ```
 */
export function column_stack(tup: NDArray[]): NDArray {
  const arrays: NDArray[] = [];
  const toDispose: NDArray[] = [];

  for (const arr of tup) {
    if (arr.ndim < 2) {
      // 1D to column: [N] -> [N, 1]
      const col = arr.reshape([arr.size, 1]);
      arrays.push(col);
      toDispose.push(col);
    } else {
      arrays.push(arr);
    }
  }

  const result = concatenate(arrays, 1);

  for (const arr of toDispose) {
    arr.dispose();
  }

  return result;
}

/**
 * Assemble an nd-array from nested lists of blocks.
 *
 * @param arrays - Nested list of arrays (1D list for hstack, 2D for block assembly)
 * @returns Assembled array
 *
 * @example
 * ```typescript
 * const A = await NDArray.fromArray([[1, 2], [3, 4]]);
 * const B = await NDArray.fromArray([[5], [6]]);
 * const result = await block([[A, B]]);
 * // [[1, 2, 5], [3, 4, 6]]
 * ```
 */
export function block(arrays: (NDArray | NDArray[])[]): NDArray {
  // Handle simple 1D list case
  if (arrays.every(arr => arr instanceof NDArray)) {
    return hstack(arrays as NDArray[]);
  }

  // Handle 2D block case: list of lists
  const rows: NDArray[] = [];
  const toDispose: NDArray[] = [];

  for (const row of arrays) {
    if (Array.isArray(row)) {
      const hstacked = hstack(row);
      rows.push(hstacked);
      toDispose.push(hstacked);
    } else {
      rows.push(row);
    }
  }

  const result = vstack(rows);

  for (const arr of toDispose) {
    arr.dispose();
  }

  return result;
}

/**
 * Append values to the end of an array.
 *
 * @param arr - Input array
 * @param values - Values to append
 * @param axis - Axis along which to append (default: flatten both)
 * @returns Array with appended values
 *
 * @example
 * ```typescript
 * const a = await NDArray.fromArray([1, 2, 3]);
 * const b = await NDArray.fromArray([4, 5, 6]);
 * const result = await append(a, b);
 * // [1, 2, 3, 4, 5, 6]
 * ```
 */
export function append(arr: NDArray, values: NDArray, axis?: number): NDArray {
  if (axis === undefined) {
    // Flatten both and concatenate
    const flat1 = arr.ravel();
    const flat2 = values.ravel();
    const result = concatenate([flat1, flat2], 0);
    flat1.dispose();
    flat2.dispose();
    return result;
  }

  return concatenate([arr, values], axis);
}

/* ============ Splitting Functions ============ */

/**
 * Split an array into multiple sub-arrays as views.
 * Allows unequal division (unlike split).
 *
 * @param arr - Array to split
 * @param indices_or_sections - Number of sections or split indices
 * @param axis - Axis along which to split (default: 0)
 * @returns List of sub-arrays
 */
export function array_split(
  arr: NDArray,
  indices_or_sections: number | number[],
  axis: number = 0
): NDArray[] {
  // Normalize axis
  axis = axis < 0 ? axis + arr.ndim : axis;
  if (axis < 0 || axis >= arr.ndim) {
    throw new Error(`axis ${axis} out of bounds for array with ${arr.ndim} dimensions`);
  }

  const axisSize = arr.shape[axis];
  let divPoints: number[];

  if (typeof indices_or_sections === 'number') {
    // Number of sections - compute division points
    const sections = indices_or_sections;
    const base = Math.floor(axisSize / sections);
    const extras = axisSize % sections;

    // First 'extras' sections get base+1 elements
    divPoints = [0];
    for (let i = 0; i < sections; i++) {
      divPoints.push(divPoints[i] + base + (i < extras ? 1 : 0));
    }
  } else {
    // Explicit indices
    divPoints = [0, ...indices_or_sections, axisSize];
  }

  // Create sub-arrays using slicing
  const result: NDArray[] = [];

  for (let i = 0; i < divPoints.length - 1; i++) {
    const start = divPoints[i];
    const stop = divPoints[i + 1];

    // Build slice indices
    const indices: (Slice | number)[] = [];
    for (let d = 0; d < arr.ndim; d++) {
      if (d === axis) {
        indices.push(slice(start, stop));
      } else {
        indices.push(slice(null, null)); // Full slice
      }
    }

    const subArr = arr.slice(indices);
    result.push(subArr.copy()); // Copy to own data
    subArr.dispose();
  }

  return result;
}

/**
 * Split an array into multiple sub-arrays of equal size.
 * Raises error if array cannot be split into equal sections.
 *
 * @param arr - Array to split
 * @param indices_or_sections - Number of equal sections or split indices
 * @param axis - Axis along which to split (default: 0)
 * @returns List of sub-arrays
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([1, 2, 3, 4, 5, 6]);
 * const result = await split(arr, 3);
 * // [[1, 2], [3, 4], [5, 6]]
 * ```
 */
export function split(
  arr: NDArray,
  indices_or_sections: number | number[],
  axis: number = 0
): NDArray[] {
  // Normalize axis
  const normalizedAxis = axis < 0 ? axis + arr.ndim : axis;

  if (typeof indices_or_sections === 'number') {
    const axisSize = arr.shape[normalizedAxis];
    if (axisSize % indices_or_sections !== 0) {
      throw new Error('array split does not result in an equal division');
    }
  }

  return array_split(arr, indices_or_sections, axis);
}

/**
 * Split an array into multiple sub-arrays vertically (row-wise).
 *
 * @param arr - Array with at least 2 dimensions
 * @param indices_or_sections - Number of sections or split indices
 * @returns List of sub-arrays
 */
export function vsplit(
  arr: NDArray,
  indices_or_sections: number | number[]
): NDArray[] {
  if (arr.ndim < 2) {
    throw new Error('vsplit only works on arrays of 2 or more dimensions');
  }
  return split(arr, indices_or_sections, 0);
}

/**
 * Split an array into multiple sub-arrays horizontally (column-wise).
 *
 * @param arr - Array
 * @param indices_or_sections - Number of sections or split indices
 * @returns List of sub-arrays
 */
export function hsplit(
  arr: NDArray,
  indices_or_sections: number | number[]
): NDArray[] {
  if (arr.ndim < 1) {
    throw new Error('hsplit only works on arrays of 1 or more dimensions');
  }
  if (arr.ndim > 1) {
    return split(arr, indices_or_sections, 1);
  }
  return split(arr, indices_or_sections, 0);
}

/**
 * Split array into multiple sub-arrays along the 3rd axis (depth).
 *
 * @param arr - Array with at least 3 dimensions
 * @param indices_or_sections - Number of sections or split indices
 * @returns List of sub-arrays
 */
export function dsplit(
  arr: NDArray,
  indices_or_sections: number | number[]
): NDArray[] {
  if (arr.ndim < 3) {
    throw new Error('dsplit only works on arrays of 3 or more dimensions');
  }
  return split(arr, indices_or_sections, 2);
}

/**
 * Unpack the array along a given axis.
 * The inverse of stack.
 *
 * @param arr - Array to unpack
 * @param axis - Axis along which to unpack (default: 0)
 * @returns List of arrays, each with one fewer dimension
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([[1, 2], [3, 4]]);
 * const [a, b] = await unstack(arr, 0);
 * // a: [1, 2], b: [3, 4]
 * ```
 */
export function unstack(arr: NDArray, axis: number = 0): NDArray[] {
  // Normalize axis
  axis = axis < 0 ? axis + arr.ndim : axis;

  if (arr.ndim === 0) {
    throw new Error('cannot unstack 0-d array');
  }

  const numArrays = arr.shape[axis];
  return split(arr, numArrays, axis).map(sub => {
    // Remove the singleton dimension created by split
    // squeeze() returns a view, so we need to copy it before disposing sub
    const squeezed = sub.squeeze(axis);
    const copied = squeezed.copy();
    squeezed.dispose();
    sub.dispose();
    return copied;
  });
}

/* ============ Rearranging Functions ============ */

/**
 * Reverse the order of elements in an array along the given axis.
 *
 * @param arr - Input array
 * @param axis - Axis or axes along which to flip. Default is to flip all axes.
 * @returns View with elements reversed along axis
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([[1, 2], [3, 4]]);
 * const result = flip(arr, 0);
 * // [[3, 4], [1, 2]]
 * ```
 */
export function flip(arr: NDArray, axis?: number | number[]): NDArray {
  if (axis === undefined) {
    // Flip all axes
    let result = arr;
    let isFirst = true;
    for (let d = 0; d < arr.ndim; d++) {
      const flipped = flipSingle(result, d);
      if (!isFirst) result.dispose();
      result = flipped;
      isFirst = false;
    }
    return result;
  }

  if (typeof axis === 'number') {
    return flipSingle(arr, axis);
  }

  // Multiple axes
  let result = arr;
  let isFirst = true;
  for (const ax of axis) {
    const flipped = flipSingle(result, ax);
    if (!isFirst) result.dispose();
    result = flipped;
    isFirst = false;
  }
  return result;
}

function flipSingle(arr: NDArray, axis: number): NDArray {
  axis = axis < 0 ? axis + arr.ndim : axis;
  if (axis < 0 || axis >= arr.ndim) {
    throw new Error(`axis ${axis} out of bounds for array with ${arr.ndim} dimensions`);
  }

  // Use slicing with step=-1 - creates view with negative stride
  const indices: (Slice | number)[] = [];
  for (let d = 0; d < arr.ndim; d++) {
    if (d === axis) {
      indices.push(slice(null, null, -1)); // Reverse this axis
    } else {
      indices.push(slice(null, null)); // Full slice
    }
  }

  return arr.slice(indices);
}

/**
 * Reverse the order of elements along axis 1 (left/right).
 *
 * @param arr - Input array with at least 2 dimensions
 * @returns Flipped view
 */
export function fliplr(arr: NDArray): NDArray {
  if (arr.ndim < 2) {
    throw new Error('fliplr requires array with at least 2 dimensions');
  }
  return flip(arr, 1);
}

/**
 * Reverse the order of elements along axis 0 (up/down).
 *
 * @param arr - Input array with at least 1 dimension
 * @returns Flipped view
 */
export function flipud(arr: NDArray): NDArray {
  if (arr.ndim < 1) {
    throw new Error('flipud requires array with at least 1 dimension');
  }
  return flip(arr, 0);
}

/**
 * Roll array elements along a given axis.
 *
 * @param arr - Input array
 * @param shift - Number of places to shift elements
 * @param axis - Axis along which to roll (default: flatten, roll, reshape)
 * @returns Array with elements rolled
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([1, 2, 3, 4, 5]);
 * const result = await roll(arr, 2);
 * // [4, 5, 1, 2, 3]
 * ```
 */
export function roll(
  arr: NDArray,
  shift: number | number[],
  axis?: number | number[]
): NDArray {
  if (axis === undefined) {
    // Roll flattened array
    const flat = arr.flatten();
    const rolled = rollSingle(flat, shift as number, 0);
    const result = rolled.reshape([...arr.shape]);
    flat.dispose();
    rolled.dispose();
    return result;
  }

  if (typeof axis === 'number' && typeof shift === 'number') {
    return rollSingle(arr, shift, axis);
  }

  // Multiple axes
  const shifts = Array.isArray(shift) ? shift : [shift];
  const axes = Array.isArray(axis) ? axis : [axis];

  if (shifts.length !== axes.length) {
    throw new Error('shift and axis must have same length');
  }

  let result = arr;
  let isFirst = true;
  for (let i = 0; i < axes.length; i++) {
    const rolled = rollSingle(result, shifts[i], axes[i]);
    if (!isFirst) result.dispose();
    result = rolled;
    isFirst = false;
  }

  return result;
}

function rollSingle(arr: NDArray, shift: number, axis: number): NDArray {
  axis = axis < 0 ? axis + arr.ndim : axis;
  const axisSize = arr.shape[axis];

  // Normalize shift
  shift = ((shift % axisSize) + axisSize) % axisSize;
  if (shift === 0) return arr.copy();

  // Split and concatenate
  const indices1: Slice[] = [];
  const indices2: Slice[] = [];
  for (let d = 0; d < arr.ndim; d++) {
    if (d === axis) {
      indices1.push(slice(axisSize - shift, null));
      indices2.push(slice(null, axisSize - shift));
    } else {
      indices1.push(slice(null, null));
      indices2.push(slice(null, null));
    }
  }

  const part1 = arr.slice(indices1);
  const part2 = arr.slice(indices2);
  const result = concatenate([part1, part2], axis);
  part1.dispose();
  part2.dispose();

  return result;
}

/**
 * Rotate an array by 90 degrees in the plane specified by axes.
 *
 * @param arr - Array with at least 2 dimensions
 * @param k - Number of times to rotate by 90 degrees (default: 1)
 * @param axes - The two axes that define the plane of rotation (default: [0, 1])
 * @returns Rotated array
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([[1, 2], [3, 4]]);
 * const result = rot90(arr, 1);
 * // [[2, 4], [1, 3]]
 * ```
 */
export function rot90(
  arr: NDArray,
  k: number = 1,
  axes: [number, number] = [0, 1]
): NDArray {
  if (arr.ndim < 2) {
    throw new Error('rot90 requires array with at least 2 dimensions');
  }

  k = ((k % 4) + 4) % 4; // Normalize to 0-3

  const [ax1, ax2] = axes;

  if (k === 0) return arr.copy();

  if (k === 1) {
    // Transpose then flip
    const swapped = arr.swapaxes(ax1, ax2);
    const result = flipSingle(swapped, ax1);
    swapped.dispose();
    return result;
  }

  if (k === 2) {
    // Flip both axes
    const flipped1 = flipSingle(arr, ax1);
    const result = flipSingle(flipped1, ax2);
    flipped1.dispose();
    return result;
  }

  // k === 3
  const swapped = arr.swapaxes(ax1, ax2);
  const result = flipSingle(swapped, ax2);
  swapped.dispose();
  return result;
}

/**
 * Return a new array with the specified shape.
 * If the new array is larger than the original, it is filled with repeated copies.
 *
 * @param arr - Input array
 * @param new_shape - Shape of the resized array
 * @returns Resized array
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([1, 2, 3]);
 * const result = await resize(arr, [5]);
 * // [1, 2, 3, 1, 2]
 * ```
 */
export async function resize(arr: NDArray, new_shape: number[]): Promise<NDArray> {
  const newSize = new_shape.reduce((a, b) => a * b, 1);
  const result = await NDArray.empty(new_shape, { dtype: arr.dtype });

  // Fill by cycling through original data
  for (let i = 0; i < newSize; i++) {
    result.setFlat(i, arr.getFlat(i % arr.size));
  }

  return result;
}

/**
 * Trim the leading and/or trailing zeros from a 1-D array.
 *
 * @param arr - Input 1-D array
 * @param trim - 'f' for front, 'b' for back, 'fb' for both (default: 'fb')
 * @returns Trimmed array
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([0, 0, 1, 2, 3, 0, 0]);
 * const result = trim_zeros(arr);
 * // [1, 2, 3]
 * ```
 */
export function trim_zeros(arr: NDArray, trim: string = 'fb'): NDArray {
  if (arr.ndim !== 1) {
    throw new Error('trim_zeros only works on 1-D arrays');
  }

  let start = 0;
  let end = arr.size;

  if (trim.includes('f')) {
    while (start < end && arr.getFlat(start) === 0) start++;
  }

  if (trim.includes('b')) {
    while (end > start && arr.getFlat(end - 1) === 0) end--;
  }

  return arr.slice([slice(start, end)]);
}

/* ============ Tiling Functions ============ */

/**
 * Construct an array by repeating arr the number of times given by reps.
 *
 * @param arr - Input array
 * @param reps - Number of repetitions along each axis
 * @returns Tiled array
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([1, 2]);
 * const result = await tile(arr, 3);
 * // [1, 2, 1, 2, 1, 2]
 *
 * const arr2 = await NDArray.fromArray([[1, 2], [3, 4]]);
 * const result2 = await tile(arr2, [2, 1]);
 * // [[1, 2], [3, 4], [1, 2], [3, 4]]
 * ```
 */
export async function tile(arr: NDArray, reps: number | number[]): Promise<NDArray> {
  const repsArr = Array.isArray(reps) ? reps : [reps];

  // Extend reps to match or exceed arr.ndim
  const d = Math.max(arr.ndim, repsArr.length);

  // Pad reps with 1s at the beginning if needed
  const paddedReps = new Array(d - repsArr.length).fill(1).concat(repsArr);

  // Pad array shape with 1s at the beginning if needed
  let tiled = arr;
  let needsDispose = false;
  if (arr.ndim < d) {
    const newShape = new Array(d - arr.ndim).fill(1).concat([...arr.shape]);
    tiled = arr.reshape(newShape);
    needsDispose = true;
  }

  // Compute output shape
  const outShape = tiled.shape.map((s, i) => s * paddedReps[i]);

  // Create result
  const result = await NDArray.empty(outShape, { dtype: arr.dtype });

  // Fill by iterating over output positions
  for (let outFlat = 0; outFlat < result.size; outFlat++) {
    // Convert to multi-index
    const outIdx = flatToMulti(outFlat, outShape);

    // Map to source index (modulo original shape)
    const srcIdx = outIdx.map((idx, dim) => idx % tiled.shape[dim]);
    const srcFlat = multiToFlat(srcIdx, tiled.shape);

    result.setFlat(outFlat, tiled.getFlat(srcFlat));
  }

  if (needsDispose) {
    tiled.dispose();
  }

  return result;
}

/**
 * Repeat elements of an array.
 *
 * @param arr - Input array
 * @param repeats - Number of repetitions for each element
 * @param axis - Axis along which to repeat (default: flatten first)
 * @returns Array with repeated elements
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([1, 2, 3]);
 * const result = await repeat(arr, 2);
 * // [1, 1, 2, 2, 3, 3]
 *
 * const arr2 = await NDArray.fromArray([[1, 2], [3, 4]]);
 * const result2 = await repeat(arr2, 2, 0);
 * // [[1, 2], [1, 2], [3, 4], [3, 4]]
 * ```
 */
export async function repeat(
  arr: NDArray,
  repeats: number | number[],
  axis?: number
): Promise<NDArray> {
  if (axis === undefined) {
    // Flatten and repeat
    const flat = arr.ravel();
    const result = await repeatAlongAxis(flat, repeats, 0);
    flat.dispose();
    return result;
  }

  return repeatAlongAxis(arr, repeats, axis);
}

async function repeatAlongAxis(
  arr: NDArray,
  repeats: number | number[],
  axis: number
): Promise<NDArray> {
  axis = axis < 0 ? axis + arr.ndim : axis;

  const axisSize = arr.shape[axis];
  const repsArr =
    typeof repeats === 'number'
      ? new Array(axisSize).fill(repeats)
      : repeats;

  if (repsArr.length !== axisSize) {
    throw new Error('repeats must have same length as axis');
  }

  const totalRepeats = repsArr.reduce((a, b) => a + b, 0);

  // Compute output shape
  const outShape = [...arr.shape];
  outShape[axis] = totalRepeats;

  const result = await NDArray.empty(outShape, { dtype: arr.dtype });

  // Use swapaxes trick: move axis to front
  const swapped = arr.swapaxes(axis, 0);
  const outSwapped = result.swapaxes(axis, 0);

  let outIdx = 0;
  for (let i = 0; i < axisSize; i++) {
    const srcSlice = swapped.at(i);
    const numReps = repsArr[i];

    for (let r = 0; r < numReps; r++) {
      const dstSlice = outSwapped.at(outIdx);
      copyArrayData(srcSlice, dstSlice);
      dstSlice.dispose();
      outIdx++;
    }
    srcSlice.dispose();
  }

  swapped.dispose();
  outSwapped.dispose();

  return result;
}

/**
 * Pad an array.
 *
 * @param arr - Input array
 * @param pad_width - Number of values padded to edges of each axis
 * @param mode - Padding mode (default: 'constant')
 * @param constant_values - Value to use for constant padding (default: 0)
 * @returns Padded array
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([1, 2, 3]);
 * const result = await pad(arr, 2, 'constant', 0);
 * // [0, 0, 1, 2, 3, 0, 0]
 *
 * const arr2 = await NDArray.fromArray([[1, 2], [3, 4]]);
 * const result2 = await pad(arr2, [[1, 1], [2, 2]], 'constant', 0);
 * // 0 padding around the original array
 * ```
 */
export async function pad(
  arr: NDArray,
  pad_width: number | [number, number] | [number, number][],
  mode: string = 'constant',
  constant_values: number = 0
): Promise<NDArray> {
  // Normalize pad_width to [[before, after], ...] format
  let padSpec: [number, number][];

  if (typeof pad_width === 'number') {
    padSpec = arr.shape.map(() => [pad_width, pad_width]);
  } else if (Array.isArray(pad_width) && typeof pad_width[0] === 'number') {
    const [before, after] = pad_width as [number, number];
    padSpec = arr.shape.map(() => [before, after]);
  } else {
    padSpec = pad_width as [number, number][];
  }

  if (padSpec.length !== arr.ndim) {
    throw new Error('pad_width must have same length as array dimensions');
  }

  // Compute output shape
  const outShape = arr.shape.map((s, i) => s + padSpec[i][0] + padSpec[i][1]);

  // Create padded array with constant fill
  const result = await NDArray.full(outShape, constant_values, { dtype: arr.dtype });

  // Build slice to place original data
  const indices: Slice[] = padSpec.map(
    ([before], i) => slice(before, before + arr.shape[i])
  );

  // Copy original data into the center
  const target = result.slice(indices);
  copyArrayData(arr, target);
  target.dispose();

  // Handle other modes
  if (mode !== 'constant') {
    // For edge mode, reflect, etc., we'd need additional implementation
    // For now, only constant mode is fully supported
    if (mode === 'edge') {
      padEdge(result, arr.shape, padSpec);
    }
    // Other modes could be added: 'reflect', 'symmetric', 'wrap'
  }

  return result;
}

function padEdge(
  result: NDArray,
  originalShape: number[],
  padSpec: [number, number][]
): void {
  // Edge padding - replicate edge values
  const outShape = result.shape;

  for (let i = 0; i < result.size; i++) {
    const outIdx = flatToMulti(i, outShape);

    // Map to source index, clamping to original bounds
    const srcIdx = outIdx.map((idx, d) => {
      const [before] = padSpec[d];
      const shifted = idx - before;
      return Math.max(0, Math.min(shifted, originalShape[d] - 1));
    });

    // Check if this is a padded position
    const isPadded = outIdx.some((idx, d) => {
      const [before] = padSpec[d];
      return idx < before || idx >= before + originalShape[d];
    });

    if (isPadded) {
      // Get source value and set in result
      const srcFlat = multiToFlat(
        srcIdx.map((idx, d) => idx + padSpec[d][0]),
        outShape
      );
      result.setFlat(i, result.getFlat(srcFlat));
    }
  }
}

/* ============ Insert/Delete Functions ============ */

/**
 * Insert values along the given axis before the given indices.
 *
 * @param arr - Input array
 * @param obj - Index or indices before which values are inserted
 * @param values - Values to insert
 * @param axis - Axis along which to insert (default: flatten first)
 * @returns Array with inserted values
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([1, 2, 3, 4, 5]);
 * const result = await insert(arr, 2, [10, 11]);
 * // [1, 2, 10, 11, 3, 4, 5]
 * ```
 */
export function insert(
  arr: NDArray,
  obj: number | number[],
  values: NDArray,
  axis?: number
): NDArray {
  if (axis === undefined) {
    // Insert into flattened array
    const flat = arr.ravel();
    const result = insertAlongAxis(flat, obj, values, 0);
    flat.dispose();
    return result;
  }

  return insertAlongAxis(arr, obj, values, axis);
}

function insertAlongAxis(
  arr: NDArray,
  obj: number | number[],
  values: NDArray,
  axis: number
): NDArray {
  axis = axis < 0 ? axis + arr.ndim : axis;

  const indices = typeof obj === 'number' ? [obj] : [...obj];
  const sortedIndices = indices.sort((a, b) => a - b);

  const parts: NDArray[] = [];
  const toDispose: NDArray[] = [];
  let lastIdx = 0;

  for (let i = 0; i < sortedIndices.length; i++) {
    const idx = sortedIndices[i];

    // Add part before insertion point
    if (idx > lastIdx) {
      const sliceIndices: Slice[] = [];
      for (let d = 0; d < arr.ndim; d++) {
        sliceIndices.push(d === axis ? slice(lastIdx, idx) : slice(null, null));
      }
      const part = arr.slice(sliceIndices);
      parts.push(part);
      toDispose.push(part);
    }

    // Add inserted values
    parts.push(values);
    lastIdx = idx;
  }

  // Add remaining part
  if (lastIdx < arr.shape[axis]) {
    const sliceIndices: Slice[] = [];
    for (let d = 0; d < arr.ndim; d++) {
      sliceIndices.push(
        d === axis ? slice(lastIdx, null) : slice(null, null)
      );
    }
    const part = arr.slice(sliceIndices);
    parts.push(part);
    toDispose.push(part);
  }

  const result = concatenate(parts, axis);

  for (const p of toDispose) {
    p.dispose();
  }

  return result;
}

/**
 * Return a new array with sub-arrays along an axis deleted.
 *
 * @param arr - Input array
 * @param obj - Index or indices to delete
 * @param axis - Axis along which to delete (default: flatten first)
 * @returns Array with elements deleted
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([1, 2, 3, 4, 5]);
 * const result = await deleteArr(arr, [1, 3]);
 * // [1, 3, 5]
 * ```
 */
export async function deleteArr(
  arr: NDArray,
  obj: number | number[],
  axis?: number
): Promise<NDArray> {
  if (axis === undefined) {
    const flat = arr.ravel();
    const result = await deleteAlongAxis(flat, obj, 0);
    flat.dispose();
    return result;
  }

  return deleteAlongAxis(arr, obj, axis);
}

async function deleteAlongAxis(
  arr: NDArray,
  obj: number | number[],
  axis: number
): Promise<NDArray> {
  axis = axis < 0 ? axis + arr.ndim : axis;

  const indices = typeof obj === 'number' ? [obj] : obj;
  const deleteSet = new Set(
    indices.map(i => (i < 0 ? i + arr.shape[axis] : i))
  );

  // Keep indices not in deleteSet
  const keepIndices: number[] = [];
  for (let i = 0; i < arr.shape[axis]; i++) {
    if (!deleteSet.has(i)) keepIndices.push(i);
  }

  // Use take to select kept indices
  const indicesArr = await NDArray.fromArray(keepIndices, [keepIndices.length], {
    dtype: DType.Int32,
  });
  const result = take(arr, indicesArr, axis);
  indicesArr.dispose();

  return result;
}

/* ============ Copying Functions ============ */

/**
 * Copy values from one array to another, broadcasting as necessary.
 *
 * @param dst - Destination array
 * @param src - Source array
 * @param casting - Type casting rule (default: 'same_kind')
 * @param where - Optional boolean mask array
 */
export function copyto(
  dst: NDArray,
  src: NDArray,
  _casting: string = 'same_kind',
  where?: NDArray
): void {
  if (where === undefined) {
    // Simple copy
    if (dst.size !== src.size) {
      throw new Error('dst and src must have same size');
    }
    for (let i = 0; i < dst.size; i++) {
      dst.setFlat(i, src.getFlat(i));
    }
  } else {
    // Conditional copy
    if (dst.size !== src.size || dst.size !== where.size) {
      throw new Error('dst, src, and where must have same size');
    }
    for (let i = 0; i < dst.size; i++) {
      if (where.getFlat(i) !== 0) {
        dst.setFlat(i, src.getFlat(i));
      }
    }
  }
}

/**
 * Convert the input to an array.
 *
 * @param a - Input data
 * @param dtype - Optional data type
 * @returns Array from input
 */
export async function asarray(
  a: NDArray | number[] | number,
  dtype?: DType
): Promise<NDArray> {
  if (a instanceof NDArray) {
    if (dtype !== undefined && a.dtype !== dtype) {
      return a.astype(dtype);
    }
    return a; // Return as-is (not a copy)
  }

  const data = typeof a === 'number' ? [a] : a;
  const shape = typeof a === 'number' ? [] : [data.length];
  return NDArray.fromArray(data, shape, { dtype });
}
