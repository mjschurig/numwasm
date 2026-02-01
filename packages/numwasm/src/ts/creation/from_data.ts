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

/**
 * Create a 1D array from text data in a string.
 *
 * A simple way to quickly create 1D arrays from raw text. The string
 * is parsed using a separator character.
 *
 * @param string - String containing the data
 * @param sep - Separator between data items. If empty string, treats the
 *              string as binary data with each character representing a byte.
 * @param count - Number of items to read. Default is -1 (read all).
 * @param options - Optional configuration (dtype)
 * @returns Promise resolving to 1D NDArray
 *
 * @example
 * ```typescript
 * // Read whitespace-separated numbers
 * const a = await fromstring("1 2 3 4 5");
 * // [1, 2, 3, 4, 5]
 *
 * // Read comma-separated numbers
 * const b = await fromstring("1,2,3,4,5", ",");
 * // [1, 2, 3, 4, 5]
 *
 * // Read with count limit
 * const c = await fromstring("1 2 3 4 5", " ", 3);
 * // [1, 2, 3]
 *
 * // Read binary data (character codes)
 * const d = await fromstring("ABC", "");
 * // [65, 66, 67]
 * ```
 */
export async function fromstring(
  string: string,
  sep: string = " ",
  count: number = -1,
  options: NDArrayOptions = {},
): Promise<NDArray> {
  let data: number[];

  if (sep === "") {
    // Binary mode: each character is a byte
    data = [];
    for (let i = 0; i < string.length; i++) {
      data.push(string.charCodeAt(i));
    }
  } else {
    // Text mode: split by separator and parse as numbers
    const parts = string.split(sep).filter((s) => s.trim() !== "");
    data = parts.map((s) => {
      const num = parseFloat(s.trim());
      if (Number.isNaN(num)) {
        throw new Error(`Could not convert string to float: '${s}'`);
      }
      return num;
    });
  }

  // Apply count limit
  if (count >= 0 && count < data.length) {
    data = data.slice(0, count);
  }

  return array(data, undefined, options);
}
