/**
 * Triangular System Solver
 *
 * Solve Ax = b where A is upper or lower triangular.
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
import type { Matrix, SolveTriangularOptions, SolveTriangularResult } from './types.js';

/**
 * Solve a triangular system Ax = b where A is upper or lower triangular.
 *
 * Uses DTRTRS for efficient triangular solve without factorization.
 *
 * @param A - Triangular coefficient matrix (n × n).
 * @param b - Right-hand side vector or matrix.
 * @param options - Solver options including triangle type and transpose mode.
 * @returns Solution x and diagnostic information
 *
 * @example
 * ```typescript
 * // Solve upper triangular system
 * const U = [[2, 1], [0, 3]];
 * const b = [4, 6];
 * const { x } = solveTriangular(U, b, { upper: true });
 * // x = [1, 2]
 * ```
 *
 * @example
 * ```typescript
 * // Solve lower triangular system with unit diagonal
 * const L = [[1, 0], [2, 1]];
 * const b = [1, 4];
 * const { x } = solveTriangular(L, b, { upper: false, unitDiagonal: true });
 * ```
 */
export function solveTriangular(
  A: Matrix,
  b: Matrix,
  options: SolveTriangularOptions = {}
): SolveTriangularResult {
  const Module = getLAPACKModule();

  const {
    upper = true,
    trans = 'N',
    unitDiagonal = false,
    overwriteB = false,
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

  // Character codes
  const uplo = upper ? CHAR.U : CHAR.L;
  const transCode = trans === 'T' ? CHAR.T : trans === 'C' ? CHAR.C : CHAR.N;
  const diag = unitDiagonal ? CHAR.UNIT : CHAR.N;

  // Allocate WASM memory
  const uploPtr = allocateInts(Module, [uplo]);
  const transPtr = allocateInts(Module, [transCode]);
  const diagPtr = allocateInts(Module, [diag]);
  const nPtr = allocateInts(Module, [n]);
  const nrhsPtr = allocateInts(Module, [nrhs]);
  const aPtr = allocateDoubles(Module, aData);
  const ldaPtr = allocateInts(Module, [n]);
  const bPtr = allocateDoubles(Module, overwriteB ? null : bData, n * nrhs);
  const ldbPtr = allocateInts(Module, [n]);
  const infoPtr = allocateInts(Module, [0]);

  if (!overwriteB) {
    const baseIdx = bPtr >> 3;
    for (let i = 0; i < bData.length; i++) {
      Module.HEAPF64[baseIdx + i] = bData[i];
    }
  }

  try {
    // Call DTRTRS
    Module._dtrtrs_(
      uploPtr,
      transPtr,
      diagPtr,
      nPtr,
      nrhsPtr,
      aPtr,
      ldaPtr,
      bPtr,
      ldbPtr,
      infoPtr
    );

    const info = readInt(Module, infoPtr);
    const success = info === 0;
    const message = getLapackErrorMessage('DTRTRS', info);
    const x = readDoubles(Module, bPtr, n * nrhs);

    return {
      x,
      nrhs,
      info,
      success,
      message,
    };
  } finally {
    freeAll(Module, [
      uploPtr,
      transPtr,
      diagPtr,
      nPtr,
      nrhsPtr,
      aPtr,
      ldaPtr,
      bPtr,
      ldbPtr,
      infoPtr,
    ]);
  }
}
