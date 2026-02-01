/**
 * LU Factorization
 *
 * Compute the LU factorization of a general matrix using partial pivoting.
 */

import { getLAPACKModule } from '../loader.js';
import {
  allocateDoubles,
  allocateInts,
  readDoubles,
  readInt,
  freeAll,
  prepareMatrix,
  getMatrixDimensions,
  getLapackErrorMessage,
} from '../helpers.js';
import type { Matrix, LUOptions, LUResult } from './types.js';

/**
 * Compute the LU factorization of a general m×n matrix A.
 *
 * A = P * L * U
 *
 * where P is a permutation matrix, L is lower triangular with unit diagonal,
 * and U is upper triangular.
 *
 * @param A - Input matrix (m × n). Can be 2D array (row-major) or 1D array (column-major).
 * @param options - Factorization options
 * @returns LU factorization result with L, U, and P
 *
 * @example
 * ```typescript
 * const A = [[2, 1, 1], [4, 3, 3], [8, 7, 9]];
 * const { L, U, P, success } = lu(A);
 * // L is lower triangular with unit diagonal
 * // U is upper triangular
 * // P contains pivot indices
 * ```
 */
export function lu(A: Matrix, options: LUOptions = {}): LUResult {
  const Module = getLAPACKModule();

  // Get dimensions
  const [m, n] = getMatrixDimensions(A);
  const minMN = Math.min(m, n);

  // Prepare matrix
  const aData = prepareMatrix(A);

  // Allocate WASM memory
  const aPtr = allocateDoubles(Module, options.overwriteA ? null : aData, m * n);
  const ipivPtr = allocateInts(Module, null, minMN);
  const mPtr = allocateInts(Module, [m]);
  const nPtr = allocateInts(Module, [n]);
  const ldaPtr = allocateInts(Module, [m]);
  const infoPtr = allocateInts(Module, [0]);

  // Write data if not overwriting
  if (!options.overwriteA) {
    const baseIdx = aPtr >> 3;
    for (let i = 0; i < aData.length; i++) {
      Module.HEAPF64[baseIdx + i] = aData[i];
    }
  }

  try {
    // Call DGETRF
    Module._dgetrf_(mPtr, nPtr, aPtr, ldaPtr, ipivPtr, infoPtr);

    const info = readInt(Module, infoPtr);
    const success = info >= 0; // info > 0 means singular but factorization completed
    const message = info === 0
      ? 'Success'
      : info > 0
        ? `U(${info},${info}) is exactly zero. The factorization completed, but U is singular.`
        : getLapackErrorMessage('DGETRF', info);

    // Read the combined LU factors
    const LU = readDoubles(Module, aPtr, m * n);

    // Read pivot indices
    const P = new Int32Array(minMN);
    const ipivBase = ipivPtr >> 2;
    for (let i = 0; i < minMN; i++) {
      P[i] = Module.HEAP32[ipivBase + i];
    }

    // Extract L and U from the combined LU array
    // L is m×k lower triangular with unit diagonal (k = min(m,n))
    // U is k×n upper triangular
    const L = new Float64Array(m * minMN);
    const U = new Float64Array(minMN * n);

    // Extract L (column-major)
    for (let j = 0; j < minMN; j++) {
      for (let i = 0; i < m; i++) {
        if (i > j) {
          // Below diagonal: copy from LU
          L[j * m + i] = LU[j * m + i];
        } else if (i === j) {
          // Diagonal: L has 1
          L[j * m + i] = 1.0;
        } else {
          // Above diagonal: L has 0
          L[j * m + i] = 0.0;
        }
      }
    }

    // Extract U (column-major)
    for (let j = 0; j < n; j++) {
      for (let i = 0; i < minMN; i++) {
        if (i <= j) {
          // On or above diagonal: copy from LU
          U[j * minMN + i] = LU[j * m + i];
        } else {
          // Below diagonal: U has 0
          U[j * minMN + i] = 0.0;
        }
      }
    }

    return {
      L,
      U,
      P,
      LU,
      m,
      n,
      info,
      success,
      message,
    };
  } finally {
    freeAll(Module, [aPtr, ipivPtr, mPtr, nPtr, ldaPtr, infoPtr]);
  }
}
