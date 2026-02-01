/**
 * Hermitian Test
 *
 * Test if a matrix is Hermitian.
 */

import {
  prepareMatrix,
  getMatrixDimensions,
} from '../helpers.js';
import type { Matrix, IsHermitianOptions, IsHermitianResult } from './types.js';

/**
 * Test if a matrix is Hermitian (A = A^H).
 *
 * A matrix is Hermitian if A[i,j] = conj(A[j,i]) for all i,j.
 * For real matrices, this is equivalent to being symmetric.
 *
 * @param A - Input square matrix (n × n)
 * @param options - Test options
 * @returns Whether the matrix is Hermitian
 *
 * @example
 * ```typescript
 * // Real symmetric matrix is Hermitian
 * const A = [[1, 2], [2, 4]];
 * const { isHermitian } = isHermitian(A);
 * // isHermitian = true
 * ```
 */
export function isHermitian(A: Matrix, options: IsHermitianOptions = {}): IsHermitianResult {
  const { tol = 1e-10 } = options;

  // Get dimensions
  const [m, n] = getMatrixDimensions(A);

  // Must be square
  if (m !== n) {
    return {
      isHermitian: false,
      maxDeviation: Infinity,
      n: m,
      success: true,
      message: 'Matrix is not square',
    };
  }

  // Handle edge case
  if (n === 0) {
    return {
      isHermitian: true,
      maxDeviation: 0,
      n: 0,
      success: true,
      message: 'Success',
    };
  }

  // Prepare matrix (column-major)
  const aData = prepareMatrix(A);

  // For real matrices, Hermitian = symmetric
  // Check A[i,j] = A[j,i]
  let maxDeviation = 0;
  let hermitian = true;

  // Also check that diagonal is real (for complex matrices, diagonal must be real)
  // For real matrices, this is always satisfied

  for (let i = 0; i < n; i++) {
    for (let j = i + 1; j < n; j++) {
      const aij = aData[i + j * n];
      const aji = aData[j + i * n];
      const deviation = Math.abs(aij - aji);

      if (deviation > maxDeviation) {
        maxDeviation = deviation;
      }

      if (deviation > tol) {
        hermitian = false;
      }
    }
  }

  return {
    isHermitian: hermitian,
    maxDeviation,
    n,
    success: true,
    message: 'Success',
  };
}
