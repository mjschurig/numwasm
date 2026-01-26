/**
 * Array creation functions.
 *
 * These are standalone wrapper functions that provide a NumPy-style API
 * for array creation operations.
 */

import { NDArray } from './NDArray.js';
import type { NDArrayOptions } from './types.js';

/**
 * Create a new NDArray filled with zeros.
 *
 * @param shape - Array dimensions
 * @param options - Optional configuration (dtype)
 * @returns Promise resolving to the new NDArray
 *
 * @example
 * ```typescript
 * const a = await zeros([3]);        // [0, 0, 0]
 * const b = await zeros([2, 3]);     // [[0, 0, 0], [0, 0, 0]]
 * ```
 */
export async function zeros(
  shape: number[],
  options: NDArrayOptions = {}
): Promise<NDArray> {
  return NDArray.zeros(shape, options);
}

/**
 * Create a new NDArray filled with ones.
 *
 * @param shape - Array dimensions
 * @param options - Optional configuration (dtype)
 * @returns Promise resolving to the new NDArray
 *
 * @example
 * ```typescript
 * const a = await ones([3]);         // [1, 1, 1]
 * const b = await ones([2, 3]);      // [[1, 1, 1], [1, 1, 1]]
 * ```
 */
export async function ones(
  shape: number[],
  options: NDArrayOptions = {}
): Promise<NDArray> {
  return NDArray.ones(shape, options);
}

/**
 * Create a new NDArray without initializing values.
 * Faster than zeros() but contains arbitrary values.
 *
 * @param shape - Array dimensions
 * @param options - Optional configuration (dtype)
 * @returns Promise resolving to the new NDArray
 */
export async function empty(
  shape: number[],
  options: NDArrayOptions = {}
): Promise<NDArray> {
  return NDArray.empty(shape, options);
}

/**
 * Create a new NDArray filled with a constant value.
 *
 * @param shape - Array dimensions
 * @param fillValue - Value to fill with
 * @param options - Optional configuration (dtype)
 * @returns Promise resolving to the new NDArray
 *
 * @example
 * ```typescript
 * const a = await full([3], 5);      // [5, 5, 5]
 * const b = await full([2, 2], 7);   // [[7, 7], [7, 7]]
 * ```
 */
export async function full(
  shape: number[],
  fillValue: number,
  options: NDArrayOptions = {}
): Promise<NDArray> {
  return NDArray.full(shape, fillValue, options);
}

/**
 * Create an NDArray with evenly spaced values within a given interval.
 *
 * @param start - Start value (or end if end is not provided)
 * @param end - End value (exclusive)
 * @param step - Step between values (default: 1)
 * @returns Promise resolving to the new NDArray
 *
 * @example
 * ```typescript
 * const a = await arange(5);         // [0, 1, 2, 3, 4]
 * const b = await arange(2, 5);      // [2, 3, 4]
 * const c = await arange(0, 10, 2);  // [0, 2, 4, 6, 8]
 * ```
 */
export async function arange(
  start: number,
  end?: number,
  step: number = 1
): Promise<NDArray> {
  return NDArray.arange(start, end, step);
}

/**
 * Create an NDArray with evenly spaced numbers over a specified interval.
 *
 * @param start - Start value
 * @param stop - End value
 * @param num - Number of samples (default: 50)
 * @param endpoint - Include stop value (default: true)
 * @param options - Optional configuration (dtype)
 * @returns Promise resolving to the new NDArray
 *
 * @example
 * ```typescript
 * const a = await linspace(0, 1, 5);  // [0, 0.25, 0.5, 0.75, 1]
 * ```
 */
export async function linspace(
  start: number,
  stop: number,
  num: number = 50,
  endpoint: boolean = true,
  options: NDArrayOptions = {}
): Promise<NDArray> {
  return NDArray.linspace(start, stop, num, endpoint, options);
}

/**
 * Create an NDArray with numbers spaced evenly on a log scale.
 *
 * @param start - Start value (power of base)
 * @param stop - End value (power of base)
 * @param num - Number of samples (default: 50)
 * @param endpoint - Include stop value (default: true)
 * @param base - Base of the log scale (default: 10)
 * @param options - Optional configuration (dtype)
 * @returns Promise resolving to the new NDArray
 *
 * @example
 * ```typescript
 * const a = await logspace(0, 2, 3);  // [1, 10, 100]
 * ```
 */
export async function logspace(
  start: number,
  stop: number,
  num: number = 50,
  endpoint: boolean = true,
  base: number = 10,
  options: NDArrayOptions = {}
): Promise<NDArray> {
  return NDArray.logspace(start, stop, num, endpoint, base, options);
}

/**
 * Create an NDArray with numbers spaced evenly on a geometric scale.
 *
 * @param start - Start value (must be non-zero)
 * @param stop - End value (must be non-zero)
 * @param num - Number of samples (default: 50)
 * @param endpoint - Include stop value (default: true)
 * @param options - Optional configuration (dtype)
 * @returns Promise resolving to the new NDArray
 *
 * @example
 * ```typescript
 * const a = await geomspace(1, 1000, 4);  // [1, 10, 100, 1000]
 * ```
 */
export async function geomspace(
  start: number,
  stop: number,
  num: number = 50,
  endpoint: boolean = true,
  options: NDArrayOptions = {}
): Promise<NDArray> {
  return NDArray.geomspace(start, stop, num, endpoint, options);
}

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
  options: NDArrayOptions = {}
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
  options: NDArrayOptions = {}
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
 * Create an NDArray from a nested array or flat array with shape.
 *
 * @param data - Nested array or flat array of numbers
 * @param shape - Optional shape (inferred from nested structure if not provided)
 * @param options - Optional configuration (dtype)
 * @returns Promise resolving to the new NDArray
 *
 * @example
 * ```typescript
 * const a = await array([1, 2, 3]);
 * const b = await array([[1, 2], [3, 4]]);
 * ```
 */
export async function array(
  data: number | number[] | number[][] | number[][][] | Float64Array | Float32Array | Int32Array,
  shape?: number[],
  options: NDArrayOptions = {}
): Promise<NDArray> {
  return NDArray.fromArray(data as Parameters<typeof NDArray.fromArray>[0], shape, options);
}

// Note: asarray is exported from manipulation.js

/**
 * Create an array of zeros with the same shape and type as a given array.
 *
 * @param a - Reference array
 * @returns Promise resolving to the new NDArray
 */
export async function zeros_like(a: NDArray): Promise<NDArray> {
  return NDArray.zeros(a.shape, { dtype: a.dtype });
}

/**
 * Create an array of ones with the same shape and type as a given array.
 *
 * @param a - Reference array
 * @returns Promise resolving to the new NDArray
 */
export async function ones_like(a: NDArray): Promise<NDArray> {
  return NDArray.ones(a.shape, { dtype: a.dtype });
}

/**
 * Create an uninitialized array with the same shape and type as a given array.
 *
 * @param a - Reference array
 * @returns Promise resolving to the new NDArray
 */
export async function empty_like(a: NDArray): Promise<NDArray> {
  return NDArray.empty(a.shape, { dtype: a.dtype });
}

/**
 * Create an array filled with a value, with the same shape and type as a given array.
 *
 * @param a - Reference array
 * @param fillValue - Value to fill with
 * @returns Promise resolving to the new NDArray
 */
export async function full_like(a: NDArray, fillValue: number): Promise<NDArray> {
  return NDArray.full(a.shape, fillValue, { dtype: a.dtype });
}
