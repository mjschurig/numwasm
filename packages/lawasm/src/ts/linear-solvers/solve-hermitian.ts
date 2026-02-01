/**
 * Complex Hermitian Positive Definite System Solver
 *
 * Solve Ax = b where A is complex Hermitian positive definite.
 */

import { getLAPACKModule } from '../loader.js';
import {
  allocateDoubles,
  allocateInts,
  readDoubles,
  readInt,
  freeAll,
  getLapackErrorMessage,
  CHAR,
} from '../helpers.js';
import type { ComplexMatrix, ComplexVector, SolveHermitianOptions, SolveComplexResult } from './types.js';

/**
 * Solve a complex Hermitian positive definite system Ax = b.
 *
 * Uses complex Cholesky decomposition (ZPOSV).
 *
 * Complex values are stored as interleaved real/imaginary pairs:
 * [re0, im0, re1, im1, ...]
 *
 * @param A - Hermitian positive definite matrix in interleaved complex format.
 *            Size should be 2 * n * n where n is the matrix dimension.
 * @param b - Right-hand side in interleaved complex format.
 * @param n - Matrix dimension (number of rows/columns).
 * @param options - Solver options.
 * @returns Solution x in interleaved complex format.
 *
 * @example
 * ```typescript
 * // 2x2 Hermitian matrix: [[2, 1+i], [1-i, 3]]
 * // Stored as: [2, 0, 1, 1, 1, -1, 3, 0] (column-major, interleaved)
 * const A = [2, 0, 1, -1, 1, 1, 3, 0];
 * const b = [1, 0, 2, 0]; // [1+0i, 2+0i]
 * const { x } = solveHermitian(A, b, 2);
 * ```
 */
export function solveHermitian(
  A: ComplexMatrix,
  b: ComplexVector,
  n: number,
  options: SolveHermitianOptions = {}
): SolveComplexResult {
  const Module = getLAPACKModule();

  const { upper = true, overwriteA = false, overwriteB = false } = options;

  // Validate sizes (complex arrays have 2x the elements)
  const aLen = Array.isArray(A) ? A.length : A.length;
  const bLen = Array.isArray(b) ? b.length : b.length;

  if (aLen !== 2 * n * n) {
    throw new Error(
      `Complex matrix A should have ${2 * n * n} elements, got ${aLen}`
    );
  }

  // Determine nrhs
  const nrhs = Math.floor(bLen / (2 * n));
  if (bLen !== 2 * n * nrhs) {
    throw new Error(
      `Complex vector b length ${bLen} is not compatible with n=${n}`
    );
  }

  // Prepare data
  const aData = A instanceof Float64Array ? A : new Float64Array(A);
  const bData = b instanceof Float64Array ? b : new Float64Array(b);

  // Allocate WASM memory
  const uplo = upper ? CHAR.U : CHAR.L;
  const uploPtr = allocateInts(Module, [uplo]);
  const nPtr = allocateInts(Module, [n]);
  const nrhsPtr = allocateInts(Module, [nrhs]);
  const aPtr = allocateDoubles(Module, overwriteA ? null : aData, 2 * n * n);
  const ldaPtr = allocateInts(Module, [n]);
  const bPtr = allocateDoubles(Module, overwriteB ? null : bData, 2 * n * nrhs);
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
    // Call ZPOSV
    Module._zposv_(uploPtr, nPtr, nrhsPtr, aPtr, ldaPtr, bPtr, ldbPtr, infoPtr);

    const info = readInt(Module, infoPtr);
    const success = info === 0;
    const message = getLapackErrorMessage('ZPOSV', info);
    const x = readDoubles(Module, bPtr, 2 * n * nrhs);

    return {
      x,
      nrhs,
      info,
      success,
      message,
    };
  } finally {
    freeAll(Module, [uploPtr, nPtr, nrhsPtr, aPtr, ldaPtr, bPtr, ldbPtr, infoPtr]);
  }
}
