/**
 * Scale Operation (DSCAL)
 *
 * Compute x = alpha*x
 */

import { getLAPACKModule } from '../loader.js';
import {
  allocateDoubles,
  allocateInts,
  readDoubles,
  freeAll,
  prepareVector,
} from '../helpers.js';
import type { Vector, ScalResult } from './types.js';

/**
 * Compute x = alpha*x.
 *
 * Uses DSCAL from BLAS.
 *
 * @param alpha - Scalar multiplier
 * @param x - Vector to scale
 * @returns Scaled vector x
 *
 * @example
 * ```typescript
 * const x = [1, 2, 3];
 *
 * // x = 2*x
 * const { x: result } = scal(2, x);
 * // result = [2, 4, 6]
 * ```
 */
export function scal(alpha: number, x: Vector): ScalResult {
  const Module = getLAPACKModule();

  // Get vector data
  const xData = prepareVector(x);
  const n = xData.length;

  // Handle empty vector
  if (n === 0) {
    return {
      x: new Float64Array(0),
      n: 0,
      success: true,
      message: 'Success',
    };
  }

  // Increment
  const incx = 1;

  // Allocate WASM memory
  const xPtr = allocateDoubles(Module, xData, n);

  // Parameter pointers
  const nPtr = allocateInts(Module, [n]);
  const alphaPtr = allocateDoubles(Module, [alpha], 1);
  const incxPtr = allocateInts(Module, [incx]);

  try {
    // Call DSCAL
    Module._dscal_(nPtr, alphaPtr, xPtr, incxPtr);

    // Read result (x is overwritten)
    const result = readDoubles(Module, xPtr, n);

    return {
      x: result,
      n,
      success: true,
      message: 'Success',
    };
  } finally {
    freeAll(Module, [xPtr, nPtr, alphaPtr, incxPtr]);
  }
}
