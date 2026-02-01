/**
 * Orthogonality Test
 *
 * Test if a matrix is orthogonal.
 */

import {
  prepareMatrix,
  getMatrixDimensions,
} from '../helpers.js';
import type { Matrix, IsOrthogonalOptions, IsOrthogonalResult } from './types.js';

/**
 * Test if a matrix is orthogonal (Q^T * Q = I).
 *
 * An orthogonal matrix has orthonormal columns and rows.
 * Q^T * Q = Q * Q^T = I, and det(Q) = ±1.
 *
 * @param A - Input square matrix (n × n)
 * @param options - Test options
 * @returns Whether the matrix is orthogonal
 *
 * @example
 * ```typescript
 * // Rotation matrix (orthogonal)
 * const theta = Math.PI / 4;
 * const Q = [[Math.cos(theta), -Math.sin(theta)],
 *            [Math.sin(theta), Math.cos(theta)]];
 * const { isOrthogonal } = isOrthogonal(Q);
 * // isOrthogonal = true
 *
 * // Not orthogonal
 * const A = [[1, 2], [3, 4]];
 * const { isOrthogonal: orth } = isOrthogonal(A);
 * // orth = false
 * ```
 */
export function isOrthogonal(A: Matrix, options: IsOrthogonalOptions = {}): IsOrthogonalResult {
  const { tol = 1e-10 } = options;

  // Get dimensions
  const [m, n] = getMatrixDimensions(A);

  // Must be square
  if (m !== n) {
    return {
      isOrthogonal: false,
      maxDeviation: Infinity,
      n: m,
      success: true,
      message: 'Matrix is not square',
    };
  }

  // Handle edge case
  if (n === 0) {
    return {
      isOrthogonal: true,
      maxDeviation: 0,
      n: 0,
      success: true,
      message: 'Success',
    };
  }

  // Prepare matrix (column-major)
  const aData = prepareMatrix(A);

  // Compute Q^T * Q and check if it equals I
  // (Q^T * Q)[i,j] = sum_k Q[k,i] * Q[k,j]
  let maxDeviation = 0;
  let orthogonal = true;

  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      // Compute (Q^T * Q)[i,j]
      let sum = 0;
      for (let k = 0; k < n; k++) {
        // Q[k,i] in column-major: aData[k + i*n]
        // Q[k,j] in column-major: aData[k + j*n]
        sum += aData[k + i * n] * aData[k + j * n];
      }

      // Expected value: 1 if i==j, 0 otherwise
      const expected = i === j ? 1 : 0;
      const deviation = Math.abs(sum - expected);

      if (deviation > maxDeviation) {
        maxDeviation = deviation;
      }

      if (deviation > tol) {
        orthogonal = false;
      }
    }
  }

  return {
    isOrthogonal: orthogonal,
    maxDeviation,
    n,
    success: true,
    message: 'Success',
  };
}
