/**
 * Symmetric Rank-1 Update (DSYR)
 *
 * Compute A = alpha*x*x^T + A where A is symmetric
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
  CHAR,
} from '../helpers.js';
import type { Vector, SyrOptions, SyrResult } from './types.js';

/**
 * Compute symmetric rank-1 update A = alpha*x*x^T + A.
 *
 * Uses DSYR from BLAS. The result A is symmetric.
 *
 * @param x - Vector (length n)
 * @param options - Update options
 * @returns Symmetric result matrix A (n × n)
 *
 * @example
 * ```typescript
 * const x = [1, 2, 3];
 *
 * // A = x*x^T (3×3 symmetric result)
 * const { A } = syr(x);
 *
 * // A = 2*x*x^T + A0
 * const A0 = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
 * const { A: A2 } = syr(x, { alpha: 2, A: A0 });
 * ```
 */
export function syr(x: Vector, options: SyrOptions = {}): SyrResult {
  const Module = getLAPACKModule();

  const {
    uplo = 'lower',
    alpha = 1.0,
    A,
  } = options;

  // Get vector length
  const xData = prepareVector(x);
  const n = xData.length;

  // Prepare or create A
  let aData: Float64Array;
  if (A !== undefined) {
    const [rowsA, colsA] = getMatrixDimensions(A);
    if (rowsA !== n || colsA !== n) {
      throw new Error(
        `A dimensions ${rowsA}×${colsA} incompatible with result ${n}×${n}`
      );
    }
    aData = prepareMatrix(A);
  } else {
    aData = new Float64Array(n * n);
  }

  // Leading dimension and increment
  const lda = n;
  const incx = 1;

  // Allocate WASM memory
  const xPtr = allocateDoubles(Module, xData, n);
  const aPtr = allocateDoubles(Module, aData, n * n);

  // Parameter pointers
  const uploPtr = allocateInts(Module, [uplo === 'upper' ? CHAR.U : CHAR.L]);
  const nPtr = allocateInts(Module, [n]);
  const alphaPtr = allocateDoubles(Module, [alpha], 1);
  const ldaPtr = allocateInts(Module, [lda]);
  const incxPtr = allocateInts(Module, [incx]);

  try {
    // Call DSYR
    Module._dsyr_(
      uploPtr,
      nPtr,
      alphaPtr,
      xPtr,
      incxPtr,
      aPtr,
      ldaPtr
    );

    // Read result
    const result = readDoubles(Module, aPtr, n * n);

    // Symmetrize the result (DSYR only fills one triangle)
    for (let j = 0; j < n; j++) {
      for (let i = j + 1; i < n; i++) {
        if (uplo === 'lower') {
          // Lower triangle is filled, copy to upper
          result[i * n + j] = result[j * n + i];
        } else {
          // Upper triangle is filled, copy to lower
          result[j * n + i] = result[i * n + j];
        }
      }
    }

    return {
      A: result,
      n,
      success: true,
      message: 'Success',
    };
  } finally {
    freeAll(Module, [xPtr, aPtr, uploPtr, nPtr, alphaPtr, ldaPtr, incxPtr]);
  }
}
