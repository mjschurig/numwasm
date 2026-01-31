/**
 * Array creation functions.
 *
 * These are standalone wrapper functions that provide a NumPy-style API
 * for array creation operations.
 */

import { NDArray } from "../_core/NDArray.js";
import type { NDArrayOptions } from "../types.js";

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
  options: NDArrayOptions = {},
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
  options: NDArrayOptions = {},
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
  options: NDArrayOptions = {},
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
  options: NDArrayOptions = {},
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
  step: number = 1,
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
  options: NDArrayOptions = {},
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
  options: NDArrayOptions = {},
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
  options: NDArrayOptions = {},
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
  data:
    | number
    | number[]
    | number[][]
    | number[][][]
    | Float64Array
    | Float32Array
    | Int32Array,
  shape?: number[],
  options: NDArrayOptions = {},
): Promise<NDArray> {
  return NDArray.fromArray(
    data as Parameters<typeof NDArray.fromArray>[0],
    shape,
    options,
  );
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
export async function full_like(
  a: NDArray,
  fillValue: number,
): Promise<NDArray> {
  return NDArray.full(a.shape, fillValue, { dtype: a.dtype });
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
 * Convert the input to an array, passing through ndarrays unchanged.
 *
 * In NumPy, this preserves subclasses of ndarray. In numwasm, there are
 * no subclasses, so this behaves identically to asarray.
 *
 * @param a - Input data
 * @param options - Optional configuration (dtype)
 * @returns Promise resolving to NDArray
 */
export async function asanyarray(
  a: NDArray | number[] | number[][] | number,
  options: NDArrayOptions = {},
): Promise<NDArray> {
  if (a instanceof NDArray) {
    if (options.dtype !== undefined && a.dtype !== options.dtype) {
      return a.astype(options.dtype);
    }
    return a;
  }

  const data = typeof a === "number" ? [a] : a;
  const shape = typeof a === "number" ? [] : undefined;
  return array(data as number[], shape, options);
}

/**
 * Convert the input to an array, checking for NaNs or Infs.
 *
 * @param a - Input data
 * @param options - Optional configuration (dtype)
 * @returns Promise resolving to NDArray
 * @throws Error if input contains NaN or Inf values
 *
 * @example
 * ```typescript
 * const a = await asarray_chkfinite([1, 2, 3]);  // OK
 * await asarray_chkfinite([1, NaN, 3]);  // Throws error
 * await asarray_chkfinite([1, Infinity, 3]);  // Throws error
 * ```
 */
export async function asarray_chkfinite(
  a: NDArray | number[] | number[][] | number,
  options: NDArrayOptions = {},
): Promise<NDArray> {
  const arr = await asanyarray(a, options);

  // Check all values are finite
  for (let i = 0; i < arr.size; i++) {
    const val = arr.getFlat(i) as number;
    if (!Number.isFinite(val)) {
      throw new Error(
        `array must not contain infs or NaNs, found ${val} at index ${i}`,
      );
    }
  }
  return arr;
}

/**
 * Return a copy of an array.
 *
 * @param a - Input array
 * @returns Promise resolving to a copy of the array
 *
 * @example
 * ```typescript
 * const a = await array([[1, 2], [3, 4]]);
 * const b = await copy(a);
 * ```
 */
export async function copy(a: NDArray): Promise<NDArray> {
  return a.copy();
}

/**
 * Return a contiguous array (ndim >= 1) in memory (C order).
 *
 * @param a - Input data
 * @param options - Optional configuration (dtype)
 * @returns Promise resolving to a contiguous array in C order
 */
export async function ascontiguousarray(
  a: NDArray | number[] | number[][],
  options: NDArrayOptions = {},
): Promise<NDArray> {
  const arr =
    a instanceof NDArray ? a : await array(a as number[], undefined, options);
  return arr.ascontiguousarray();
}

/**
 * Return a contiguous array (ndim >= 1) in memory (Fortran order).
 *
 * @param a - Input data
 * @param options - Optional configuration (dtype)
 * @returns Promise resolving to a contiguous array in Fortran order
 */
export async function asfortranarray(
  a: NDArray | number[] | number[][],
  options: NDArrayOptions = {},
): Promise<NDArray> {
  const arr =
    a instanceof NDArray ? a : await array(a as number[], undefined, options);
  return arr.asfortranarray();
}
