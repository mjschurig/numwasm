/**
 * Euclidean Norm (DNRM2)
 *
 * Compute ||x||_2 = sqrt(sum(x_i^2))
 */

import { getLAPACKModule } from '../loader.js';
import {
  allocateDoubles,
  allocateInts,
  freeAll,
  prepareVector,
} from '../helpers.js';
import type { Vector, Nrm2Result } from './types.js';

/**
 * Compute Euclidean (L2) norm of a vector.
 *
 * Uses DNRM2 from BLAS.
 *
 * @param x - Input vector
 * @returns Euclidean norm ||x||_2
 *
 * @example
 * ```typescript
 * const x = [3, 4];
 *
 * const { norm } = nrm2(x);
 * // norm = 5 (sqrt(9 + 16))
 * ```
 */
export function nrm2(x: Vector): Nrm2Result {
  const Module = getLAPACKModule();

  // Get vector data
  const xData = prepareVector(x);
  const n = xData.length;

  // Handle empty vector
  if (n === 0) {
    return {
      norm: 0,
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
  const incxPtr = allocateInts(Module, [incx]);

  try {
    // Call DNRM2
    const result = Module._dnrm2_(nPtr, xPtr, incxPtr);

    return {
      norm: result,
      n,
      success: true,
      message: 'Success',
    };
  } finally {
    freeAll(Module, [xPtr, nPtr, incxPtr]);
  }
}
