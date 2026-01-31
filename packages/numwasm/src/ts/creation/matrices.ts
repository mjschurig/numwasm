/**
 * Matrix Construction Functions
 *
 * Functions for creating special matrices: eye, identity, diag, tri, tril, triu, vander.
 */

import { NDArray } from "../_core/NDArray.js";
import type { NDArrayOptions } from "../_core/types.js";
import { array } from "./basic.js";

/**
 * Create a 2D array with ones on the diagonal and zeros elsewhere.
 *
 * @param N - Number of rows
 * @param M - Number of columns (default: N)
 * @param k - Diagonal offset (default: 0, main diagonal)
 * @param options - Optional configuration (dtype)
 * @returns Promise resolving to the new NDArray
 *
 * @example
 * ```typescript
 * const a = await eye(3);            // 3x3 identity matrix
 * const b = await eye(3, 4);         // 3x4 matrix with ones on diagonal
 * const c = await eye(3, 3, 1);      // ones on first superdiagonal
 * ```
 */
export async function eye(
  N: number,
  M?: number,
  k: number = 0,
  options: NDArrayOptions = {},
): Promise<NDArray> {
  return NDArray.eye(N, M, k, options);
}

/**
 * Create a square identity matrix.
 *
 * @param n - Size of the matrix
 * @param options - Optional configuration (dtype)
 * @returns Promise resolving to the new NDArray
 *
 * @example
 * ```typescript
 * const I = await identity(3);  // 3x3 identity matrix
 * ```
 */
export async function identity(
  n: number,
  options: NDArrayOptions = {},
): Promise<NDArray> {
  return NDArray.identity(n, options);
}

/**
 * Extract a diagonal or construct a diagonal array.
 *
 * @param v - Input array (1D to create diagonal matrix, 2D to extract diagonal)
 * @param k - Diagonal offset (default: 0)
 * @returns Promise resolving to the new NDArray
 *
 * @example
 * ```typescript
 * const a = await fromArray([1, 2, 3]);
 * const diagMatrix = await diag(a);  // 3x3 matrix with [1,2,3] on diagonal
 * ```
 */
export async function diag(v: NDArray, k: number = 0): Promise<NDArray> {
  return NDArray.diag(v, k);
}

/**
 * Create a 2D array with ones at and below the given diagonal.
 *
 * @param N - Number of rows
 * @param M - Number of columns (default: N)
 * @param k - Diagonal offset (default: 0)
 * @param options - Optional configuration (dtype)
 * @returns Promise resolving to the new NDArray
 *
 * @example
 * ```typescript
 * const a = await tri(3);
 * // [[1, 0, 0],
 * //  [1, 1, 0],
 * //  [1, 1, 1]]
 * ```
 */
export async function tri(
  N: number,
  M?: number,
  k: number = 0,
  options: NDArrayOptions = {},
): Promise<NDArray> {
  return NDArray.tri(N, M, k, options);
}

/**
 * Lower triangle of an array.
 *
 * Return a copy of an array with elements above the k-th diagonal zeroed.
 *
 * @param arr - Input array
 * @param k - Diagonal above which to zero elements (default: 0)
 * @returns Promise resolving to the lower triangle of the array
 *
 * @example
 * ```typescript
 * const a = await array([[1, 2, 3], [4, 5, 6], [7, 8, 9]]);
 * const lower = await tril(a);
 * // [[1, 0, 0],
 * //  [4, 5, 0],
 * //  [7, 8, 9]]
 * ```
 */
export async function tril(arr: NDArray, k: number = 0): Promise<NDArray> {
  return NDArray.tril(arr, k);
}

/**
 * Upper triangle of an array.
 *
 * Return a copy of an array with elements below the k-th diagonal zeroed.
 *
 * @param arr - Input array
 * @param k - Diagonal below which to zero elements (default: 0)
 * @returns Promise resolving to the upper triangle of the array
 *
 * @example
 * ```typescript
 * const a = await array([[1, 2, 3], [4, 5, 6], [7, 8, 9]]);
 * const upper = await triu(a);
 * // [[1, 2, 3],
 * //  [0, 5, 6],
 * //  [0, 0, 9]]
 * ```
 */
export async function triu(arr: NDArray, k: number = 0): Promise<NDArray> {
  return NDArray.triu(arr, k);
}

/**
 * Create a 2D array with the flattened input as a diagonal.
 *
 * @param v - Input array, which will be flattened
 * @param k - Diagonal to set (default: 0, main diagonal)
 * @returns Promise resolving to a 2D array with the input on the k-th diagonal
 *
 * @example
 * ```typescript
 * const a = await diagflat([1, 2, 3]);
 * // [[1, 0, 0],
 * //  [0, 2, 0],
 * //  [0, 0, 3]]
 *
 * const b = await diagflat([[1, 2], [3, 4]]);
 * // [[1, 0, 0, 0],
 * //  [0, 2, 0, 0],
 * //  [0, 0, 3, 0],
 * //  [0, 0, 0, 4]]
 * ```
 */
export async function diagflat(
  v: NDArray | number[] | number[][],
  k: number = 0,
  options: NDArrayOptions = {},
): Promise<NDArray> {
  // Convert to NDArray if needed
  const arr =
    v instanceof NDArray ? v : await array(v as number[], undefined, options);
  // Flatten to 1D
  const flat = arr.flatten();
  // Create diagonal matrix from flattened array
  return diag(flat, k);
}

/**
 * Generate a Vandermonde matrix.
 *
 * The columns of the output matrix are powers of the input vector.
 * The order of the powers is determined by the `increasing` parameter.
 *
 * @param x - 1D input array
 * @param N - Number of columns (default: length of x)
 * @param increasing - Order of powers (default: false, decreasing)
 * @param options - Optional configuration (dtype)
 * @returns Promise resolving to the Vandermonde matrix
 *
 * @example
 * ```typescript
 * const v = await vander([1, 2, 3], 3);
 * // [[1, 1, 1],
 * //  [4, 2, 1],
 * //  [9, 3, 1]]
 *
 * const v2 = await vander([1, 2, 3], 3, true);
 * // [[1, 1, 1],
 * //  [1, 2, 4],
 * //  [1, 3, 9]]
 * ```
 */
export async function vander(
  x: NDArray | number[],
  N?: number,
  increasing: boolean = false,
  options: NDArrayOptions = {},
): Promise<NDArray> {
  const arr = x instanceof NDArray ? x : await array(x, undefined, options);
  return NDArray.vander(arr, N, increasing);
}
