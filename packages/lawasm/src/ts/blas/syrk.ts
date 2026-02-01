/**
 * Symmetric Rank-k Update (DSYRK)
 *
 * Compute C = alpha*A*A^T + beta*C or C = alpha*A^T*A + beta*C
 */

import { getLAPACKModule } from '../loader.js';
import {
  allocateDoubles,
  allocateInts,
  readDoubles,
  freeAll,
  prepareMatrix,
  getMatrixDimensions,
  CHAR,
} from '../helpers.js';
import type { Matrix, SyrkOptions, SyrkResult } from './types.js';

/**
 * Compute symmetric rank-k update C = alpha*A*A^T + beta*C or C = alpha*A^T*A + beta*C.
 *
 * Uses DSYRK from BLAS. The result C is symmetric.
 *
 * @param A - Input matrix
 * @param options - Update options
 * @returns Symmetric result matrix C
 *
 * @example
 * ```typescript
 * const A = [[1, 2, 3], [4, 5, 6]];
 *
 * // C = A*A^T (2×2 result since A is 2×3)
 * const { C } = syrk(A);
 *
 * // C = A^T*A (3×3 result)
 * const { C: C2 } = syrk(A, { trans: 'T' });
 *
 * // C = 2*A*A^T + 3*C0
 * const C0 = [[1, 0], [0, 1]];
 * const { C: C3 } = syrk(A, { alpha: 2, beta: 3, C: C0 });
 * ```
 */
export function syrk(A: Matrix, options: SyrkOptions = {}): SyrkResult {
  const Module = getLAPACKModule();

  const {
    uplo = 'lower',
    trans = 'N',
    alpha = 1.0,
    beta = 0.0,
    C,
  } = options;

  // Get dimensions
  const [rowsA, colsA] = getMatrixDimensions(A);

  // Determine n (dimension of result) and k (common dimension)
  const n = trans === 'N' ? rowsA : colsA;
  const k = trans === 'N' ? colsA : rowsA;

  // Prepare A matrix
  const aData = prepareMatrix(A);

  // Prepare or create C
  let cData: Float64Array;
  if (C !== undefined) {
    const [rowsC, colsC] = getMatrixDimensions(C);
    if (rowsC !== n || colsC !== n) {
      throw new Error(
        `C dimensions ${rowsC}×${colsC} incompatible with result ${n}×${n}`
      );
    }
    cData = prepareMatrix(C);
  } else {
    cData = new Float64Array(n * n);
  }

  // Leading dimensions
  const lda = rowsA;
  const ldc = n;

  // Allocate WASM memory
  const aPtr = allocateDoubles(Module, aData, rowsA * colsA);
  const cPtr = allocateDoubles(Module, cData, n * n);

  // Parameter pointers
  const uploPtr = allocateInts(Module, [uplo === 'upper' ? CHAR.U : CHAR.L]);
  const transPtr = allocateInts(Module, [trans === 'T' ? CHAR.T : CHAR.N]);
  const nPtr = allocateInts(Module, [n]);
  const kPtr = allocateInts(Module, [k]);
  const alphaPtr = allocateDoubles(Module, [alpha], 1);
  const betaPtr = allocateDoubles(Module, [beta], 1);
  const ldaPtr = allocateInts(Module, [lda]);
  const ldcPtr = allocateInts(Module, [ldc]);

  try {
    // Call DSYRK
    Module._dsyrk_(
      uploPtr,
      transPtr,
      nPtr,
      kPtr,
      alphaPtr,
      aPtr,
      ldaPtr,
      betaPtr,
      cPtr,
      ldcPtr
    );

    // Read result
    const result = readDoubles(Module, cPtr, n * n);

    // Symmetrize the result (DSYRK only fills one triangle)
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
      C: result,
      n,
      success: true,
      message: 'Success',
    };
  } finally {
    freeAll(Module, [
      aPtr,
      cPtr,
      uploPtr,
      transPtr,
      nPtr,
      kPtr,
      alphaPtr,
      betaPtr,
      ldaPtr,
      ldcPtr,
    ]);
  }
}
