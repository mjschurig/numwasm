/**
 * Symmetric Rank-2k Update (DSYR2K)
 *
 * Compute C = alpha*A*B^T + alpha*B*A^T + beta*C or
 *         C = alpha*A^T*B + alpha*B^T*A + beta*C
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
import type { Matrix, Syr2kOptions, Syr2kResult } from './types.js';

/**
 * Compute symmetric rank-2k update.
 *
 * If trans='N': C = alpha*A*B^T + alpha*B*A^T + beta*C
 * If trans='T': C = alpha*A^T*B + alpha*B^T*A + beta*C
 *
 * Uses DSYR2K from BLAS. The result C is symmetric.
 *
 * @param A - First input matrix
 * @param B - Second input matrix (same dimensions as A)
 * @param options - Update options
 * @returns Symmetric result matrix C
 *
 * @example
 * ```typescript
 * const A = [[1, 2, 3], [4, 5, 6]];
 * const B = [[1, 0, 1], [0, 1, 0]];
 *
 * // C = A*B^T + B*A^T (2×2 result)
 * const { C } = syr2k(A, B);
 *
 * // C = A^T*B + B^T*A (3×3 result)
 * const { C: C2 } = syr2k(A, B, { trans: 'T' });
 * ```
 */
export function syr2k(A: Matrix, B: Matrix, options: Syr2kOptions = {}): Syr2kResult {
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
  const [rowsB, colsB] = getMatrixDimensions(B);

  // Check A and B have same dimensions
  if (rowsA !== rowsB || colsA !== colsB) {
    throw new Error(
      `A (${rowsA}×${colsA}) and B (${rowsB}×${colsB}) must have same dimensions`
    );
  }

  // Determine n (dimension of result) and k (common dimension)
  const n = trans === 'N' ? rowsA : colsA;
  const k = trans === 'N' ? colsA : rowsA;

  // Prepare matrices
  const aData = prepareMatrix(A);
  const bData = prepareMatrix(B);

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
  const ldb = rowsB;
  const ldc = n;

  // Allocate WASM memory
  const aPtr = allocateDoubles(Module, aData, rowsA * colsA);
  const bPtr = allocateDoubles(Module, bData, rowsB * colsB);
  const cPtr = allocateDoubles(Module, cData, n * n);

  // Parameter pointers
  const uploPtr = allocateInts(Module, [uplo === 'upper' ? CHAR.U : CHAR.L]);
  const transPtr = allocateInts(Module, [trans === 'T' ? CHAR.T : CHAR.N]);
  const nPtr = allocateInts(Module, [n]);
  const kPtr = allocateInts(Module, [k]);
  const alphaPtr = allocateDoubles(Module, [alpha], 1);
  const betaPtr = allocateDoubles(Module, [beta], 1);
  const ldaPtr = allocateInts(Module, [lda]);
  const ldbPtr = allocateInts(Module, [ldb]);
  const ldcPtr = allocateInts(Module, [ldc]);

  try {
    // Call DSYR2K
    Module._dsyr2k_(
      uploPtr,
      transPtr,
      nPtr,
      kPtr,
      alphaPtr,
      aPtr,
      ldaPtr,
      bPtr,
      ldbPtr,
      betaPtr,
      cPtr,
      ldcPtr
    );

    // Read result
    const result = readDoubles(Module, cPtr, n * n);

    // Symmetrize the result
    for (let j = 0; j < n; j++) {
      for (let i = j + 1; i < n; i++) {
        if (uplo === 'lower') {
          result[i * n + j] = result[j * n + i];
        } else {
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
      bPtr,
      cPtr,
      uploPtr,
      transPtr,
      nPtr,
      kPtr,
      alphaPtr,
      betaPtr,
      ldaPtr,
      ldbPtr,
      ldcPtr,
    ]);
  }
}
