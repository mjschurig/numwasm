/**
 * AXPY Operation (DAXPY)
 *
 * Compute y = alpha*x + y
 */

import { getLAPACKModule } from '../loader.js';
import {
  allocateDoubles,
  allocateInts,
  readDoubles,
  freeAll,
  prepareVector,
} from '../helpers.js';
import type { Vector, AxpyResult } from './types.js';

/**
 * Compute y = alpha*x + y.
 *
 * Uses DAXPY from BLAS.
 *
 * @param alpha - Scalar multiplier
 * @param x - Vector to add
 * @param y - Vector to accumulate into
 * @returns Result vector y
 *
 * @example
 * ```typescript
 * const x = [1, 2, 3];
 * const y = [4, 5, 6];
 *
 * // y = 2*x + y
 * const { y: result } = axpy(2, x, y);
 * // result = [6, 9, 12]
 * ```
 */
export function axpy(alpha: number, x: Vector, y: Vector): AxpyResult {
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
  const alphaPtr = allocateDoubles(Module, [alpha], 1);
  const incxPtr = allocateInts(Module, [incx]);
  const incyPtr = allocateInts(Module, [incy]);

  try {
    // Call DAXPY
    Module._daxpy_(nPtr, alphaPtr, xPtr, incxPtr, yPtr, incyPtr);

    // Read result (y is overwritten)
    const result = readDoubles(Module, yPtr, n);

    return {
      y: result,
      n,
      success: true,
      message: 'Success',
    };
  } finally {
    freeAll(Module, [xPtr, yPtr, nPtr, alphaPtr, incxPtr, incyPtr]);
  }
}
