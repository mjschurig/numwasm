/**
 * Shape manipulation functions for NumJS-WASM
 *
 * Functions for reshaping, flattening, and dimension manipulation.
 */

import { NDArray } from "../_core/NDArray.js";

/**
 * Gives a new shape to an array without changing its data.
 *
 * @param a - Array to be reshaped
 * @param newshape - New shape (can include one -1 for auto-calculation)
 * @returns Reshaped array (view if possible)
 *
 * @example
 * ```typescript
 * const a = await array([1, 2, 3, 4, 5, 6]);
 * const b = reshape(a, [2, 3]);
 * // [[1, 2, 3], [4, 5, 6]]
 * ```
 */
export function reshape(a: NDArray, newshape: number[]): NDArray {
  return a.reshape(newshape);
}

/**
 * Return a contiguous flattened array.
 *
 * @param a - Input array
 * @returns A 1-D array (view if possible, copy otherwise)
 *
 * @example
 * ```typescript
 * const a = await array([[1, 2], [3, 4]]);
 * const b = ravel(a);
 * // [1, 2, 3, 4]
 * ```
 */
export function ravel(a: NDArray): NDArray {
  return a.ravel();
}

/**
 * Return a copy of the array collapsed into one dimension.
 *
 * @param a - Input array
 * @returns A 1-D copy of the array
 *
 * @example
 * ```typescript
 * const a = await array([[1, 2], [3, 4]]);
 * const b = flatten(a);
 * // [1, 2, 3, 4]
 * ```
 */
export function flatten(a: NDArray): NDArray {
  return a.flatten();
}

/**
 * Remove axes of length one from array.
 *
 * @param a - Input array
 * @param axis - Specific axis to squeeze (optional, -1 for all size-1 axes)
 * @returns Array with size-1 dimensions removed
 *
 * @example
 * ```typescript
 * const a = await array([[[1], [2], [3]]]);  // shape: [1, 3, 1]
 * const b = squeeze(a);  // shape: [3]
 * ```
 */
export function squeeze(a: NDArray, axis?: number): NDArray {
  return a.squeeze(axis);
}

/**
 * Expand the shape of an array by inserting a new axis.
 *
 * @param a - Input array
 * @param axis - Position in the expanded axes where the new axis is placed
 * @returns Array with expanded shape
 *
 * @example
 * ```typescript
 * const a = await array([1, 2, 3]);  // shape: [3]
 * const b = expand_dims(a, 0);  // shape: [1, 3]
 * const c = expand_dims(a, 1);  // shape: [3, 1]
 * ```
 */
export function expand_dims(a: NDArray, axis: number): NDArray {
  return a.expandDims(axis);
}
