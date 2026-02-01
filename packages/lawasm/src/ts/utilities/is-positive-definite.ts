/**
 * Positive Definiteness Test
 *
 * Test if a matrix is positive definite.
 */

import { getLAPACKModule } from '../loader.js';
import {
  allocateDoubles,
  allocateInts,
  readInt,
  freeAll,
  prepareMatrix,
  getMatrixDimensions,
  CHAR,
} from '../helpers.js';
import type { Matrix, IsPositiveDefiniteResult } from './types.js';

/**
 * Test if a matrix is positive definite.
 *
 * A symmetric matrix is positive definite if all its eigenvalues are positive,
 * or equivalently, if Cholesky factorization succeeds.
 *
 * @param A - Input square symmetric matrix (n × n)
 * @returns Whether the matrix is positive definite
 *
 * @example
 * ```typescript
 * // Positive definite
 * const A = [[4, 2], [2, 5]];
 * const { isPositiveDefinite } = isPositiveDefinite(A);
 * // isPositiveDefinite = true
 *
 * // Not positive definite (indefinite)
 * const B = [[1, 2], [2, 1]];
 * const { isPositiveDefinite: pd } = isPositiveDefinite(B);
 * // pd = false
 * ```
 */
export function isPositiveDefinite(A: Matrix): IsPositiveDefiniteResult {
  const Module = getLAPACKModule();

  // Get dimensions
  const [m, n] = getMatrixDimensions(A);

  // Must be square
  if (m !== n) {
    return {
      isPositiveDefinite: false,
      n: m,
      success: true,
      message: 'Matrix is not square',
    };
  }

  // Handle edge case
  if (n === 0) {
    return {
      isPositiveDefinite: true, // Empty matrix is vacuously positive definite
      n: 0,
      success: true,
      message: 'Success',
    };
  }

  // Prepare matrix
  const aData = prepareMatrix(A);

  // Allocate WASM memory
  const aPtr = allocateDoubles(Module, aData, n * n);
  const infoPtr = allocateInts(Module, [0]);

  // Parameter pointers
  const uploPtr = allocateInts(Module, [CHAR.L]); // Lower triangular
  const nPtr = allocateInts(Module, [n]);
  const ldaPtr = allocateInts(Module, [n]);

  try {
    // Attempt Cholesky factorization using DPOTRF
    // If it succeeds, the matrix is positive definite
    Module._dpotrf_(uploPtr, nPtr, aPtr, ldaPtr, infoPtr);

    const info = readInt(Module, infoPtr);

    if (info === 0) {
      // Cholesky succeeded - matrix is positive definite
      return {
        isPositiveDefinite: true,
        n,
        success: true,
        message: 'Success',
      };
    } else if (info > 0) {
      // Cholesky failed - matrix is not positive definite
      return {
        isPositiveDefinite: false,
        n,
        success: true,
        message: `Leading minor of order ${info} is not positive definite`,
      };
    } else {
      // Unexpected error
      return {
        isPositiveDefinite: false,
        n,
        success: false,
        message: `DPOTRF returned error code ${info}`,
      };
    }
  } finally {
    freeAll(Module, [aPtr, infoPtr, uploPtr, nPtr, ldaPtr]);
  }
}
