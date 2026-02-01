/**
 * Sum of Absolute Values (DASUM)
 *
 * Compute sum(|x_i|)
 */

import { getLAPACKModule } from '../loader.js';
import {
  allocateDoubles,
  allocateInts,
  freeAll,
  prepareVector,
} from '../helpers.js';
import type { Vector, AsumResult } from './types.js';

/**
 * Compute sum of absolute values of vector elements.
 *
 * Uses DASUM from BLAS.
 *
 * @param x - Input vector
 * @returns Sum of absolute values
 *
 * @example
 * ```typescript
 * const x = [1, -2, 3, -4];
 *
 * const { asum } = asum(x);
 * // asum = 10 (|1| + |-2| + |3| + |-4|)
 * ```
 */
export function asum(x: Vector): AsumResult {
  const Module = getLAPACKModule();

  // Get vector data
  const xData = prepareVector(x);
  const n = xData.length;

  // Handle empty vector
  if (n === 0) {
    return {
      asum: 0,
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
    // Call DASUM
    const result = Module._dasum_(nPtr, xPtr, incxPtr);

    return {
      asum: result,
      n,
      success: true,
      message: 'Success',
    };
  } finally {
    freeAll(Module, [xPtr, nPtr, incxPtr]);
  }
}
