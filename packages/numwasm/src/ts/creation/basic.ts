/**
 * Basic Array Creation Functions
 *
 * Functions for creating arrays with constant values: zeros, ones, empty, full, array.
 */

import { NDArray } from "../_core/NDArray.js";
import type { NDArrayOptions } from "../_core/types.js";

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
