/**
 * Swap Operation (DSWAP)
 *
 * Swap vectors x and y
 */

import { getLAPACKModule } from '../loader.js';
import {
  allocateDoubles,
  allocateInts,
  readDoubles,
  freeAll,
  prepareVector,
} from '../helpers.js';
import type { Vector, SwapResult } from './types.js';

/**
 * Swap vectors x and y.
 *
 * Uses DSWAP from BLAS.
 *
 * @param x - First vector
 * @param y - Second vector (same length as x)
 * @returns Swapped vectors x and y
 *
 * @example
 * ```typescript
 * const x = [1, 2, 3];
 * const y = [4, 5, 6];
 *
 * const { x: newX, y: newY } = swap(x, y);
 * // newX = [4, 5, 6]
 * // newY = [1, 2, 3]
 * ```
 */
export function swap(x: Vector, y: Vector): SwapResult {
  const Module = getLAPACKModule();

  // Get vector data
  const xData = prepareVector(x);
  const yData = prepareVector(y);

  // Check lengths match
  if (xData.length !== yData.length) {
    throw new Error(
      `Vector lengths must match: x has ${xData.length}, y has ${yData.length}`
    );
  }

  const n = xData.length;

  // Handle empty vectors
  if (n === 0) {
    return {
      x: new Float64Array(0),
      y: new Float64Array(0),
      n: 0,
      success: true,
      message: 'Success',
    };
  }

  // Increments
  const incx = 1;
  const incy = 1;

  // Allocate WASM memory
  const xPtr = allocateDoubles(Module, xData, n);
  const yPtr = allocateDoubles(Module, yData, n);

  // Parameter pointers
  const nPtr = allocateInts(Module, [n]);
  const incxPtr = allocateInts(Module, [incx]);
  const incyPtr = allocateInts(Module, [incy]);

  try {
    // Call DSWAP
    Module._dswap_(nPtr, xPtr, incxPtr, yPtr, incyPtr);

    // Read results (both are swapped)
    const resultX = readDoubles(Module, xPtr, n);
    const resultY = readDoubles(Module, yPtr, n);

    return {
      x: resultX,
      y: resultY,
      n,
      success: true,
      message: 'Success',
    };
  } finally {
    freeAll(Module, [xPtr, yPtr, nPtr, incxPtr, incyPtr]);
  }
}
