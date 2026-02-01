/**
 * General Linear System Solver
 *
 * Solve Ax = b for a general matrix A using LU decomposition.
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
} from '../helpers.js';
import type { Matrix, SolveOptions, SolveResult } from './types.js';

/**
 * Solve a system of linear equations Ax = b for a general matrix A.
 *
 * Uses LU decomposition with partial pivoting (DGESV).
 *
 * @param A - Coefficient matrix (n × n). Can be 2D array (row-major) or 1D array (column-major).
 * @param b - Right-hand side vector or matrix (n × nrhs).
 * @param options - Solver options
 * @returns Solution x and diagnostic information
 *
 * @example
 * ```typescript
 * // Solve 2x + y = 1, x + 3y = 2
 * const A = [[2, 1], [1, 3]];
 * const b = [1, 2];
 * const { x, success } = solve(A, b);
 * // x ≈ [0.2, 0.6]
 * ```
 *
 * @example
 * ```typescript
 * // Solve with multiple right-hand sides
 * const A = [[2, 1], [1, 3]];
 * const B = [[1, 4], [2, 5]]; // Two RHS vectors as columns
 * const { x } = solve(A, B);
 * ```
 */
export function solve(
  A: Matrix,
  b: Matrix,
  options: SolveOptions = {}
): SolveResult {
  const Module = getLAPACKModule();

  // Get dimensions
  const [m, n] = getMatrixDimensions(A);
  if (m !== n) {
    throw new Error(`Matrix A must be square, got ${m}×${n}`);
  }

  // Prepare matrices
  const aData = prepareMatrix(A);
  const bData = prepareMatrix(b);

  // Determine nrhs (number of right-hand sides)
  let nrhs: number;
  let ldb: number;
  if (is2DArray(b)) {
    nrhs = b[0].length;
    ldb = b.length;
  } else {
    // 1D vector or matrix
    if (bData.length === n) {
      nrhs = 1;
      ldb = n;
    } else {
      nrhs = Math.floor(bData.length / n);
      ldb = n;
    }
  }

  // Allocate WASM memory
  const aPtr = allocateDoubles(Module, options.overwriteA ? null : aData, n * n);
  const bPtr = allocateDoubles(Module, options.overwriteB ? null : bData, n * nrhs);
  const ipivPtr = allocateInts(Module, null, n);
  const nPtr = allocateInts(Module, [n]);
  const nrhsPtr = allocateInts(Module, [nrhs]);
  const ldaPtr = allocateInts(Module, [n]);
  const ldbPtr = allocateInts(Module, [ldb]);
  const infoPtr = allocateInts(Module, [0]);

  // Write data if not overwriting
  if (!options.overwriteA) {
    const baseIdx = aPtr >> 3;
    for (let i = 0; i < aData.length; i++) {
      Module.HEAPF64[baseIdx + i] = aData[i];
    }
  }
  if (!options.overwriteB) {
    const baseIdx = bPtr >> 3;
    for (let i = 0; i < bData.length; i++) {
      Module.HEAPF64[baseIdx + i] = bData[i];
    }
  }

  try {
    // Call DGESV
    Module._dgesv_(nPtr, nrhsPtr, aPtr, ldaPtr, ipivPtr, bPtr, ldbPtr, infoPtr);

    const info = readInt(Module, infoPtr);
    const success = info === 0;
    const message = getLapackErrorMessage('DGESV', info);

    // Read solution
    const x = readDoubles(Module, bPtr, n * nrhs);

    const result: SolveResult = {
      x,
      nrhs,
      info,
      success,
      message,
    };

    // Optionally return factorization
    if (options.returnFactorization) {
      result.lu = readDoubles(Module, aPtr, n * n);
      const ipiv = new Int32Array(n);
      const ipivBase = ipivPtr >> 2;
      for (let i = 0; i < n; i++) {
        ipiv[i] = Module.HEAP32[ipivBase + i];
      }
      result.ipiv = ipiv;
    }

    return result;
  } finally {
    freeAll(Module, [aPtr, bPtr, ipivPtr, nPtr, nrhsPtr, ldaPtr, ldbPtr, infoPtr]);
  }
}
