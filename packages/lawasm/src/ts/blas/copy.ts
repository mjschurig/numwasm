/**
 * Copy Operation (DCOPY)
 *
 * Compute y = x
 */

import { getLAPACKModule } from '../loader.js';
import {
  allocateDoubles,
  allocateInts,
  readDoubles,
  freeAll,
  prepareVector,
} from '../helpers.js';
import type { Vector, CopyResult } from './types.js';

/**
 * Copy vector x to y.
 *
 * Uses DCOPY from BLAS.
 *
 * @param x - Source vector
 * @param y - Destination vector (optional, created if not provided)
 * @returns Destination vector y
 *
 * @example
 * ```typescript
 * const x = [1, 2, 3];
 *
 * // Copy x to new vector
 * const { y } = copy(x);
 * // y = [1, 2, 3]
 *
 * // Copy x into existing vector
 * const dest = [0, 0, 0];
 * const { y: y2 } = copy(x, dest);
 * // y2 = [1, 2, 3]
 * ```
 */
export function copy(x: Vector, y?: Vector): CopyResult {
  const Module = getLAPACKModule();

  // Get vector data
  const xData = prepareVector(x);
  const n = xData.length;

  // Handle empty vector
  if (n === 0) {
    return {
      y: new Float64Array(0),
      n: 0,
      success: true,
      message: 'Success',
    };
  }

  // Prepare or create y
  let yData: Float64Array;
  if (y !== undefined) {
    yData = prepareVector(y);
    if (yData.length !== n) {
      throw new Error(
        `Destination vector length ${yData.length} must equal source length ${n}`
      );
    }
  } else {
    yData = new Float64Array(n);
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
    // Call DCOPY
    Module._dcopy_(nPtr, xPtr, incxPtr, yPtr, incyPtr);

    // Read result
    const result = readDoubles(Module, yPtr, n);

    return {
      y: result,
      n,
      success: true,
      message: 'Success',
    };
  } finally {
    freeAll(Module, [xPtr, yPtr, nPtr, incxPtr, incyPtr]);
  }
}
