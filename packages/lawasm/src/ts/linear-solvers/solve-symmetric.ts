/**
 * Symmetric Positive Definite System Solver
 *
 * Solve Ax = b where A is symmetric positive definite using Cholesky decomposition.
 */

import { getLAPACKModule } from '../loader.js';
import {
  allocateDoubles,
  allocateInts,
  readDoubles,
  readInt,
  freeAll,
  prepareMatrix,
  is2DArray,
  getMatrixDimensions,
  getLapackErrorMessage,
  CHAR,
} from '../helpers.js';
import type { Matrix, SolveSymmetricOptions, SolveSymmetricResult } from './types.js';

/**
 * Solve a symmetric positive definite system Ax = b.
 *
 * Uses Cholesky decomposition (DPOSV) which is more efficient than LU
 * for symmetric positive definite matrices.
 *
 * @param A - Symmetric positive definite matrix (n × n). Only the upper or lower
 *            triangle is referenced based on the `upper` option.
 * @param b - Right-hand side vector or matrix.
 * @param options - Solver options including which triangle to use.
 * @returns Solution x and diagnostic information
 *
 * @throws Error if matrix is not positive definite (info > 0)
 *
 * @example
 * ```typescript
 * // Solve SPD system
 * const A = [[4, 2], [2, 5]]; // Symmetric positive definite
 * const b = [1, 2];
 * const { x } = solveSymmetric(A, b);
 * ```
 *
 * @example
 * ```typescript
 * // Use lower triangle and get Cholesky factor
 * const { x, cholesky } = solveSymmetric(A, b, {
 *   upper: false,
 *   returnFactorization: true
 * });
 * ```
 */
export function solveSymmetric(
  A: Matrix,
  b: Matrix,
  options: SolveSymmetricOptions = {}
): SolveSymmetricResult {
  const Module = getLAPACKModule();

  const {
    upper = true,
    overwriteA = false,
    overwriteB = false,
    returnFactorization = false,
  } = options;

  // Get dimensions
  const [m, n] = getMatrixDimensions(A);
  if (m !== n) {
    throw new Error(`Matrix A must be square, got ${m}×${n}`);
  }

  // Prepare matrices
  const aData = prepareMatrix(A);
  const bData = prepareMatrix(b);

  // Determine nrhs
  let nrhs: number;
  if (is2DArray(b)) {
    nrhs = b[0].length;
  } else {
    nrhs = bData.length === n ? 1 : Math.floor(bData.length / n);
  }

  // Allocate WASM memory
  const uplo = upper ? CHAR.U : CHAR.L;
  const uploPtr = allocateInts(Module, [uplo]);
  const nPtr = allocateInts(Module, [n]);
  const nrhsPtr = allocateInts(Module, [nrhs]);
  const aPtr = allocateDoubles(Module, overwriteA ? null : aData, n * n);
  const ldaPtr = allocateInts(Module, [n]);
  const bPtr = allocateDoubles(Module, overwriteB ? null : bData, n * nrhs);
  const ldbPtr = allocateInts(Module, [n]);
  const infoPtr = allocateInts(Module, [0]);

  if (!overwriteA) {
    const baseIdx = aPtr >> 3;
    for (let i = 0; i < aData.length; i++) {
      Module.HEAPF64[baseIdx + i] = aData[i];
    }
  }
  if (!overwriteB) {
    const baseIdx = bPtr >> 3;
    for (let i = 0; i < bData.length; i++) {
      Module.HEAPF64[baseIdx + i] = bData[i];
    }
  }

  try {
    // Call DPOSV
    Module._dposv_(uploPtr, nPtr, nrhsPtr, aPtr, ldaPtr, bPtr, ldbPtr, infoPtr);

    const info = readInt(Module, infoPtr);
    const success = info === 0;
    const message = getLapackErrorMessage('DPOSV', info);
    const x = readDoubles(Module, bPtr, n * nrhs);

    const result: SolveSymmetricResult = {
      x,
      nrhs,
      info,
      success,
      message,
    };

    if (returnFactorization) {
      result.cholesky = readDoubles(Module, aPtr, n * n);
    }

    return result;
  } finally {
    freeAll(Module, [uploPtr, nPtr, nrhsPtr, aPtr, ldaPtr, bPtr, ldbPtr, infoPtr]);
  }
}
