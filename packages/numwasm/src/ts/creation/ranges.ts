/**
 * Range-based Array Creation Functions
 *
 * Functions for creating arrays with sequential values: arange, linspace, logspace, geomspace.
 */

import { NDArray } from "../_core/NDArray.js";
import type { NDArrayOptions } from "../_core/types.js";

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
