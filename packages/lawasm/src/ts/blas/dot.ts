/**
 * Dot Product (DDOT)
 *
 * Compute x^T * y
 */

import { getLAPACKModule } from '../loader.js';
import {
  allocateDoubles,
  allocateInts,
  freeAll,
  prepareVector,
} from '../helpers.js';
import type { Vector, DotResult } from './types.js';

/**
 * Compute dot product x^T * y.
 *
 * Uses DDOT from BLAS.
 *
 * @param x - First vector
 * @param y - Second vector (same length as x)
 * @returns Dot product result
 *
 * @example
 * ```typescript
 * const x = [1, 2, 3];
 * const y = [4, 5, 6];
 *
 * const { dot } = dot(x, y);
 * // dot = 1*4 + 2*5 + 3*6 = 32
 * ```
 */
export function dot(x: Vector, y: Vector): DotResult {
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
      dot: 0,
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
    // Call DDOT
    const result = Module._ddot_(nPtr, xPtr, incxPtr, yPtr, incyPtr);

    return {
      dot: result,
      n,
      success: true,
      message: 'Success',
    };
  } finally {
    freeAll(Module, [xPtr, yPtr, nPtr, incxPtr, incyPtr]);
  }
}
