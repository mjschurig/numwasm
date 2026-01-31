/**
 * Array rearranging functions for NumJS-WASM
 *
 * Functions for flipping, rolling, rotating, and resizing arrays.
 */

import { NDArray } from "../_core/NDArray.js";
import { Slice, slice } from "../slice.js";
import { concatenate } from "./join.js";

/* ============ Flip Functions ============ */

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

  if (typeof axis === "number") {
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
    throw new Error(
      `axis ${axis} out of bounds for array with ${arr.ndim} dimensions`,
    );
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
    throw new Error("fliplr requires array with at least 2 dimensions");
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
    throw new Error("flipud requires array with at least 1 dimension");
  }
  return flip(arr, 0);
}

/* ============ Roll Functions ============ */

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
  axis?: number | number[],
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

  if (typeof axis === "number" && typeof shift === "number") {
    return rollSingle(arr, shift, axis);
  }

  // Multiple axes
  const shifts = Array.isArray(shift) ? shift : [shift];
  const axes = Array.isArray(axis) ? axis : [axis];

  if (shifts.length !== axes.length) {
    throw new Error("shift and axis must have same length");
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

/* ============ Rotation Functions ============ */

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
  axes: [number, number] = [0, 1],
): NDArray {
  if (arr.ndim < 2) {
    throw new Error("rot90 requires array with at least 2 dimensions");
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

/* ============ Resize Functions ============ */

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
export async function resize(
  arr: NDArray,
  new_shape: number[],
): Promise<NDArray> {
  const newSize = new_shape.reduce((a, b) => a * b, 1);

  if (newSize === 0) {
    return NDArray.empty(new_shape, { dtype: arr.dtype });
  }

  // Flatten the array
  const flat = arr.ravel();

  // Calculate how many times we need to repeat
  const repeats = Math.ceil(newSize / arr.size);

  // Use concatenate to repeat the flattened array (like NumPy)
  let expanded: NDArray;
  if (repeats === 1) {
    expanded = flat;
  } else {
    const copies: NDArray[] = [];
    for (let i = 0; i < repeats; i++) {
      copies.push(flat);
    }
    expanded = concatenate(copies, 0);
    flat.dispose();
  }

  // Take only the elements we need
  let result: NDArray;
  if (expanded.size === newSize) {
    result = expanded.reshape(new_shape);
  } else {
    // Slice to exact size, then reshape
    const sliced = expanded.slice([slice(0, newSize)]);
    result = sliced.reshape(new_shape);
    sliced.dispose();
    expanded.dispose();
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
export function trim_zeros(arr: NDArray, trim: string = "fb"): NDArray {
  if (arr.ndim !== 1) {
    throw new Error("trim_zeros only works on 1-D arrays");
  }

  let start = 0;
  let end = arr.size;

  if (trim.includes("f")) {
    while (start < end && arr.getFlat(start) === 0) start++;
  }

  if (trim.includes("b")) {
    while (end > start && arr.getFlat(end - 1) === 0) end--;
  }

  return arr.slice([slice(start, end)]);
}
