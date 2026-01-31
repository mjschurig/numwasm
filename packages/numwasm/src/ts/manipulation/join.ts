/**
 * Array joining functions for NumJS-WASM
 *
 * Functions for concatenating, stacking, and assembling arrays.
 */

import { NDArray } from "../_core/NDArray.js";
import { DType } from "../types.js";
import { atleast_1d, atleast_2d, atleast_3d } from "../indexing.js";

/* ============ Helper Functions ============ */

function arraysEqual(a: number[], b: number[]): boolean {
  if (a.length !== b.length) return false;
  for (let i = 0; i < a.length; i++) {
    if (a[i] !== b[i]) return false;
  }
  return true;
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
  dtype?: DType,
): NDArray {
  if (arrays.length === 0) {
    throw new Error("need at least one array to concatenate");
  }

  const module = arrays[0]._wasmModule;

  // Allocate array of pointers
  const ptrsPtr = module._malloc(arrays.length * 4);
  for (let i = 0; i < arrays.length; i++) {
    module.setValue(ptrsPtr + i * 4, arrays[i]._wasmPtr, "i32");
  }

  const resultPtr = module._ndarray_concatenate(ptrsPtr, arrays.length, axis);
  module._free(ptrsPtr);

  if (resultPtr === 0) {
    throw new Error("concatenate failed: incompatible shapes or invalid axis");
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
 * Join a sequence of arrays along an existing axis.
 *
 * This is an alias for `concatenate`, added for NumPy 2.0 compatibility.
 *
 * @param arrays - Sequence of arrays with compatible shapes
 * @param axis - The axis along which to join (default: 0)
 * @param dtype - Optional output dtype
 * @returns Concatenated array
 *
 * @example
 * ```typescript
 * const a = await NDArray.fromArray([1, 2, 3]);
 * const b = await NDArray.fromArray([4, 5, 6]);
 * const result = concat([a, b]);
 * // [1, 2, 3, 4, 5, 6]
 * ```
 */
export const concat = concatenate;

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
    throw new Error("need at least one array to stack");
  }

  // Validate all arrays have same shape
  const shape = arrays[0].shape;
  for (let i = 1; i < arrays.length; i++) {
    if (!arraysEqual(arrays[i].shape, shape)) {
      throw new Error("all input arrays must have the same shape");
    }
  }

  // Normalize axis
  const ndim = arrays[0].ndim + 1;
  axis = axis < 0 ? axis + ndim : axis;
  if (axis < 0 || axis > arrays[0].ndim) {
    throw new Error(
      `axis ${axis} out of bounds for array with ${ndim} dimensions`,
    );
  }

  // Expand each array along the new axis
  const expanded = arrays.map((arr) => arr.expandDims(axis));

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
  if (arrays.every((arr) => arr instanceof NDArray)) {
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
