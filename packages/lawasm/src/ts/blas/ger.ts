/**
 * Rank-1 Update (DGER)
 *
 * Compute A = alpha*x*y^T + A
 */

import { getLAPACKModule } from '../loader.js';
import {
  allocateDoubles,
  allocateInts,
  readDoubles,
  freeAll,
  prepareMatrix,
  prepareVector,
  getMatrixDimensions,
} from '../helpers.js';
import type { Vector, GerOptions, GerResult } from './types.js';

/**
 * Compute rank-1 update A = alpha*x*y^T + A.
 *
 * Uses DGER from BLAS.
 *
 * @param x - Column vector (length m)
 * @param y - Row vector (length n)
 * @param options - Update options
 * @returns Result matrix A (m × n)
 *
 * @example
 * ```typescript
 * const x = [1, 2, 3];
 * const y = [4, 5];
 *
 * // A = x*y^T (3×2 result)
 * const { A } = ger(x, y);
 *
 * // A = 2*x*y^T + A0
 * const A0 = [[1, 0], [0, 1], [1, 1]];
 * const { A: A2 } = ger(x, y, { alpha: 2, A: A0 });
 * ```
 */
export function ger(x: Vector, y: Vector, options: GerOptions = {}): GerResult {
  const Module = getLAPACKModule();

  const { alpha = 1.0, A } = options;

  // Get vector lengths
  const xData = prepareVector(x);
  const yData = prepareVector(y);
  const m = xData.length;
  const n = yData.length;

  // Prepare or create A
  let aData: Float64Array;
  if (A !== undefined) {
    const [rowsA, colsA] = getMatrixDimensions(A);
    if (rowsA !== m || colsA !== n) {
      throw new Error(
        `A dimensions ${rowsA}×${colsA} incompatible with result ${m}×${n}`
      );
    }
    aData = prepareMatrix(A);
  } else {
    aData = new Float64Array(m * n);
  }

  // Leading dimension and increments
  const lda = m;
  const incx = 1;
  const incy = 1;

  // Allocate WASM memory
  const xPtr = allocateDoubles(Module, xData, m);
  const yPtr = allocateDoubles(Module, yData, n);
  const aPtr = allocateDoubles(Module, aData, m * n);

  // Parameter pointers
  const mPtr = allocateInts(Module, [m]);
  const nPtr = allocateInts(Module, [n]);
  const alphaPtr = allocateDoubles(Module, [alpha], 1);
  const ldaPtr = allocateInts(Module, [lda]);
  const incxPtr = allocateInts(Module, [incx]);
  const incyPtr = allocateInts(Module, [incy]);

  try {
    // Call DGER
    Module._dger_(
      mPtr,
      nPtr,
      alphaPtr,
      xPtr,
      incxPtr,
      yPtr,
      incyPtr,
      aPtr,
      ldaPtr
    );

    // Read result
    const result = readDoubles(Module, aPtr, m * n);

    return {
      A: result,
      m,
      n,
      success: true,
      message: 'Success',
    };
  } finally {
    freeAll(Module, [xPtr, yPtr, aPtr, mPtr, nPtr, alphaPtr, ldaPtr, incxPtr, incyPtr]);
  }
}
