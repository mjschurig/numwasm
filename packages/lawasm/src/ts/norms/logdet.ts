/**
 * Log-Determinant
 *
 * Compute the natural logarithm of the absolute value of the determinant.
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
import type { Matrix, LogDetOptions, LogDetResult } from './types.js';

/**
 * Compute the natural logarithm of the absolute value of the determinant.
 *
 * This is numerically more stable than computing det() for large matrices
 * or matrices with very large/small determinants.
 *
 * log|det(A)| = sum of log|diagonal elements of U from LU factorization|
 *
 * @param A - Input square matrix (n × n)
 * @param options - Computation options
 * @returns Log of absolute determinant
 *
 * @example
 * ```typescript
 * const A = [[1, 2], [3, 4]];
 * const { logdet } = logdet(A);
 * // logdet = log(2) ≈ 0.693
 *
 * // For positive definite matrices, this equals log(det(A))
 * const S = [[4, 2], [2, 5]];
 * const { logdet: ld } = logdet(S);
 * // ld = log(16) ≈ 2.773
 * ```
 */
export function logdet(A: Matrix, options: LogDetOptions = {}): LogDetResult {
  const Module = getLAPACKModule();

  const { overwriteA: _overwriteA = false } = options;

  // Get dimensions
  const [m, n] = getMatrixDimensions(A);
  if (m !== n) {
    throw new Error(`Matrix A must be square, got ${m}×${n}`);
  }

  // Handle edge case
  if (n === 0) {
    return {
      logdet: 0, // log(1) = 0
      n: 0,
      info: 0,
      success: true,
      message: 'Success',
    };
  }

  // Prepare matrix
  const aData = prepareMatrix(A);

  // Allocate WASM memory
  const aPtr = allocateDoubles(Module, aData, n * n);
  const ipivPtr = allocateInts(Module, null, n);
  const infoPtr = allocateInts(Module, [0]);
  const mPtr = allocateInts(Module, [n]);
  const nPtr = allocateInts(Module, [n]);
  const ldaPtr = allocateInts(Module, [n]);

  try {
    // LU factorization using DGETRF
    Module._dgetrf_(mPtr, nPtr, aPtr, ldaPtr, ipivPtr, infoPtr);

    const info = readInt(Module, infoPtr);
    if (info < 0) {
      return {
        logdet: NaN,
        n,
        info,
        success: false,
        message: getLapackErrorMessage('DGETRF', info),
      };
    }

    if (info > 0) {
      // Matrix is singular (U has a zero on diagonal)
      return {
        logdet: -Infinity,
        n,
        info,
        success: true,
        message: 'Matrix is singular',
      };
    }

    // Read the LU factorization
    const luData = readDoubles(Module, aPtr, n * n);

    // Compute log|det| as sum of log|diagonal elements of U|
    let logDetValue = 0.0;

    for (let i = 0; i < n; i++) {
      const diag = luData[i * n + i]; // Column-major: element (i,i)
      if (diag === 0) {
        return {
          logdet: -Infinity,
          n,
          info: 0,
          success: true,
          message: 'Matrix is singular',
        };
      }
      logDetValue += Math.log(Math.abs(diag));
    }

    return {
      logdet: logDetValue,
      n,
      info: 0,
      success: true,
      message: 'Success',
    };
  } finally {
    freeAll(Module, [aPtr, ipivPtr, infoPtr, mPtr, nPtr, ldaPtr]);
  }
}
