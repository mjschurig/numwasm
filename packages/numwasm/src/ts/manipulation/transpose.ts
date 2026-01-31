/**
 * Transpose and axis manipulation functions for NumJS-WASM
 *
 * Functions for transposing, swapping axes, and moving axes.
 */

import { NDArray } from "../_core/NDArray.js";

/**
 * Reverse or permute the axes of an array.
 *
 * @param a - Input array
 * @param axes - Permutation of axes (optional, reverses all axes if not provided)
 * @returns Transposed array (view)
 *
 * @example
 * ```typescript
 * const a = await array([[1, 2, 3], [4, 5, 6]]);  // shape: [2, 3]
 * const b = transpose(a);  // shape: [3, 2]
 * ```
 */
export function transpose(a: NDArray, axes?: number[]): NDArray {
  return a.transpose(axes);
}

/**
 * Permute the dimensions of an array.
 * This is an alias for transpose, added in NumPy 2.0.
 *
 * @param a - Input array
 * @param axes - Permutation of axes (optional, reverses all axes if not provided)
 * @returns Array with permuted dimensions (view)
 */
export function permute_dims(a: NDArray, axes?: number[]): NDArray {
  return a.transpose(axes);
}

/**
 * Interchange two axes of an array.
 *
 * @param a - Input array
 * @param axis1 - First axis
 * @param axis2 - Second axis
 * @returns Array with axes swapped (view)
 *
 * @example
 * ```typescript
 * const a = await array([[1, 2], [3, 4]]);
 * const b = swapaxes(a, 0, 1);
 * ```
 */
export function swapaxes(a: NDArray, axis1: number, axis2: number): NDArray {
  return a.swapaxes(axis1, axis2);
}

/**
 * Move axes of an array to new positions.
 *
 * @param a - Input array
 * @param source - Original positions of the axes to move
 * @param destination - Destination positions for each axis
 * @returns Array with moved axes (view)
 *
 * @example
 * ```typescript
 * const a = await zeros([3, 4, 5]);
 * const b = moveaxis(a, 0, -1);  // shape: [4, 5, 3]
 * const c = moveaxis(a, [0, 1], [-1, -2]);  // shape: [5, 4, 3]
 * ```
 */
export function moveaxis(
  a: NDArray,
  source: number | number[],
  destination: number | number[],
): NDArray {
  const ndim = a.ndim;

  // Normalize to arrays
  const srcArr = typeof source === "number" ? [source] : [...source];
  const dstArr =
    typeof destination === "number" ? [destination] : [...destination];

  if (srcArr.length !== dstArr.length) {
    throw new Error("source and destination must have same length");
  }

  // Normalize negative indices
  const normalizedSrc = srcArr.map((s) => (s < 0 ? s + ndim : s));
  const normalizedDst = dstArr.map((d) => (d < 0 ? d + ndim : d));

  // Validate axes
  for (const s of normalizedSrc) {
    if (s < 0 || s >= ndim) {
      throw new Error(
        `source axis ${s} out of bounds for array with ${ndim} dimensions`,
      );
    }
  }
  for (const d of normalizedDst) {
    if (d < 0 || d >= ndim) {
      throw new Error(
        `destination axis ${d} out of bounds for array with ${ndim} dimensions`,
      );
    }
  }

  // Check for duplicates
  if (new Set(normalizedSrc).size !== normalizedSrc.length) {
    throw new Error("repeated axis in source");
  }
  if (new Set(normalizedDst).size !== normalizedDst.length) {
    throw new Error("repeated axis in destination");
  }

  // Build the permutation
  // Start with axes not being moved, in order
  const order: number[] = [];
  for (let i = 0; i < ndim; i++) {
    if (!normalizedSrc.includes(i)) {
      order.push(i);
    }
  }

  // Insert moved axes at their destinations
  // Sort destinations to insert in correct order
  const moves = normalizedSrc.map((s, i) => ({
    src: s,
    dst: normalizedDst[i],
  }));
  moves.sort((a, b) => a.dst - b.dst);

  for (const move of moves) {
    order.splice(move.dst, 0, move.src);
  }

  return a.transpose(order);
}

/**
 * Roll the specified axis backwards, until it lies in a given position.
 *
 * @param a - Input array
 * @param axis - The axis to be rolled
 * @param start - Start position (default: 0)
 * @returns Array with rolled axis (view)
 *
 * @example
 * ```typescript
 * const a = await zeros([3, 4, 5]);
 * // Equivalent to: moveaxis(a, 2, 0)
 * const b = rollaxis(a, 2);  // shape: [5, 3, 4]
 * ```
 */
export function rollaxis(a: NDArray, axis: number, start: number = 0): NDArray {
  const ndim = a.ndim;

  // Normalize axis
  axis = axis < 0 ? axis + ndim : axis;
  if (axis < 0 || axis >= ndim) {
    throw new Error(
      `axis ${axis} out of bounds for array with ${ndim} dimensions`,
    );
  }

  // Normalize start
  if (start < 0) start += ndim;
  // start can be ndim (place at end)
  if (start < 0 || start > ndim) {
    throw new Error(`start ${start} out of bounds for rollaxis`);
  }

  // If axis already at start position, no-op
  if (axis === start || axis + 1 === start) {
    return a.transpose(); // Return view with same shape
  }

  // Build permutation
  const axes: number[] = [];
  for (let i = 0; i < ndim; i++) {
    if (i !== axis) {
      axes.push(i);
    }
  }

  // Insert axis at start position
  if (start > axis) {
    // If start is after axis, adjust for the removed axis
    axes.splice(start - 1, 0, axis);
  } else {
    axes.splice(start, 0, axis);
  }

  return a.transpose(axes);
}
