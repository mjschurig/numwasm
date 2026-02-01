/**
 * Symmetry Test
 *
 * Test if a matrix is symmetric.
 */

import {
  prepareMatrix,
  getMatrixDimensions,
} from '../helpers.js';
import type { Matrix, IsSymmetricOptions, IsSymmetricResult } from './types.js';

/**
 * Test if a matrix is symmetric (A = A^T).
 *
 * A matrix is symmetric if A[i,j] = A[j,i] for all i,j.
 *
 * @param A - Input square matrix (n × n)
 * @param options - Test options
 * @returns Whether the matrix is symmetric
 *
 * @example
 * ```typescript
 * const A = [[1, 2, 3], [2, 4, 5], [3, 5, 6]];
 * const { isSymmetric } = isSymmetric(A);
 * // isSymmetric = true
 *
 * const B = [[1, 2], [3, 4]];
 * const { isSymmetric: sym } = isSymmetric(B);
 * // sym = false (because B[0,1]=2 but B[1,0]=3)
 * ```
 */
export function isSymmetric(A: Matrix, options: IsSymmetricOptions = {}): IsSymmetricResult {
  const { tol = 1e-10 } = options;

  // Get dimensions
  const [m, n] = getMatrixDimensions(A);

  // Must be square
  if (m !== n) {
    return {
      isSymmetric: false,
      maxDeviation: Infinity,
      n: m,
      success: true,
      message: 'Matrix is not square',
    };
  }

  // Handle edge case
  if (n === 0) {
    return {
      isSymmetric: true,
      maxDeviation: 0,
      n: 0,
      success: true,
      message: 'Success',
    };
  }

  // Prepare matrix (column-major)
  const aData = prepareMatrix(A);

  // Check symmetry
  let maxDeviation = 0;
  let symmetric = true;

  for (let i = 0; i < n; i++) {
    for (let j = i + 1; j < n; j++) {
      // A[i,j] in column-major: aData[i + j*n]
      // A[j,i] in column-major: aData[j + i*n]
      const aij = aData[i + j * n];
      const aji = aData[j + i * n];
      const deviation = Math.abs(aij - aji);

      if (deviation > maxDeviation) {
        maxDeviation = deviation;
      }

      if (deviation > tol) {
        symmetric = false;
      }
    }
  }

  return {
    isSymmetric: symmetric,
    maxDeviation,
    n,
    success: true,
    message: 'Success',
  };
}
