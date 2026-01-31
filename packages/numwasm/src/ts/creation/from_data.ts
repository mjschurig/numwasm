/**
 * Data-based Array Creation Functions
 *
 * Functions for creating arrays from data sources: fromfunction, fromiter.
 */

import { NDArray } from "../_core/NDArray.js";
import type { NDArrayOptions } from "../_core/types.js";
import { empty, array } from "./basic.js";

/**
 * Construct an array by executing a function over each coordinate.
 *
 * The resulting array has a value `fn(x, y, ...)` at coordinate `(x, y, ...)`.
 *
 * @param func - Function to execute on each coordinate (receives indices as arguments)
 * @param shape - Shape of the output array
 * @param options - Optional configuration (dtype)
 * @returns Promise resolving to the new NDArray
 *
 * @example
 * ```typescript
 * const a = await fromfunction((i, j) => i + j, [3, 3]);
 * // [[0, 1, 2],
 * //  [1, 2, 3],
 * //  [2, 3, 4]]
 *
 * const b = await fromfunction((i, j) => i * j, [3, 3]);
 * // [[0, 0, 0],
 * //  [0, 1, 2],
 * //  [0, 2, 4]]
 * ```
 */
export async function fromfunction(
  func: (...indices: number[]) => number,
  shape: number[],
  options: NDArrayOptions = {},
): Promise<NDArray> {
  const result = await empty(shape, options);
  const size = shape.reduce((a, b) => a * b, 1);

  for (let flatIdx = 0; flatIdx < size; flatIdx++) {
    // Convert flat index to multi-dimensional indices
    const indices: number[] = [];
    let rem = flatIdx;
    for (let d = shape.length - 1; d >= 0; d--) {
      indices.unshift(rem % shape[d]);
      rem = Math.floor(rem / shape[d]);
    }
    result.setFlat(flatIdx, func(...indices));
  }
  return result;
}

/**
 * Create a 1D array from an iterable object.
 *
 * @param iter - Iterable yielding numbers
 * @param options - Optional configuration (dtype)
 * @returns Promise resolving to 1D NDArray
 *
 * @example
 * ```typescript
 * function* range(n) { for (let i = 0; i < n; i++) yield i; }
 * const a = await fromiter(range(5));
 * // [0, 1, 2, 3, 4]
 * ```
 */
export async function fromiter(
  iter: Iterable<number>,
  options: NDArrayOptions = {},
): Promise<NDArray> {
  const data: number[] = [...iter];
  return array(data, undefined, options);
}
