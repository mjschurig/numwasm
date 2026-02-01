/**
 * Index of Maximum Absolute Value (IDAMAX)
 *
 * Find index of element with maximum |x_i|
 */

import { getLAPACKModule } from '../loader.js';
import {
  allocateDoubles,
  allocateInts,
  freeAll,
  prepareVector,
} from '../helpers.js';
import type { Vector, IamaxResult } from './types.js';

/**
 * Find index of element with maximum absolute value.
 *
 * Uses IDAMAX from BLAS. Returns 0-based index.
 *
 * @param x - Input vector
 * @returns Index of maximum absolute value element
 *
 * @example
 * ```typescript
 * const x = [1, -5, 3, -2];
 *
 * const { index, value } = iamax(x);
 * // index = 1 (0-based)
 * // value = 5 (|−5|)
 * ```
 */
export function iamax(x: Vector): IamaxResult {
  const Module = getLAPACKModule();

  // Get vector data
  const xData = prepareVector(x);
  const n = xData.length;

  // Handle empty vector
  if (n === 0) {
    return {
      index: -1,
      value: NaN,
      n: 0,
      success: true,
      message: 'Empty vector',
    };
  }

  // Increment
  const incx = 1;

  // Allocate WASM memory
  const xPtr = allocateDoubles(Module, xData, n);

  // Parameter pointers
  const nPtr = allocateInts(Module, [n]);
  const incxPtr = allocateInts(Module, [incx]);

  try {
    // Call IDAMAX (returns 1-based index)
    const result1Based = Module._idamax_(nPtr, xPtr, incxPtr);

    // Convert to 0-based index
    const index = result1Based - 1;

    // Get the value at that index
    const value = Math.abs(xData[index]);

    return {
      index,
      value,
      n,
      success: true,
      message: 'Success',
    };
  } finally {
    freeAll(Module, [xPtr, nPtr, incxPtr]);
  }
}
