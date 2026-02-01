/**
 * Unitarity Test
 *
 * Test if a matrix is unitary.
 */

import {
  prepareMatrix,
  getMatrixDimensions,
} from '../helpers.js';
import type { Matrix, IsUnitaryOptions, IsUnitaryResult } from './types.js';

/**
 * Test if a matrix is unitary (Q^H * Q = I).
 *
 * A unitary matrix has orthonormal columns and rows under the Hermitian inner product.
 * For real matrices, this is equivalent to being orthogonal.
 *
 * @param A - Input square matrix (n × n)
 * @param options - Test options
 * @returns Whether the matrix is unitary
 *
 * @example
 * ```typescript
 * // Rotation matrix (unitary for real matrices)
 * const theta = Math.PI / 4;
 * const Q = [[Math.cos(theta), -Math.sin(theta)],
 *            [Math.sin(theta), Math.cos(theta)]];
 * const { isUnitary } = isUnitary(Q);
 * // isUnitary = true
 * ```
 */
export function isUnitary(A: Matrix, options: IsUnitaryOptions = {}): IsUnitaryResult {
  const { tol = 1e-10 } = options;

  // Get dimensions
  const [m, n] = getMatrixDimensions(A);

  // Must be square
  if (m !== n) {
    return {
      isUnitary: false,
      maxDeviation: Infinity,
      n: m,
      success: true,
      message: 'Matrix is not square',
    };
  }

  // Handle edge case
  if (n === 0) {
    return {
      isUnitary: true,
      maxDeviation: 0,
      n: 0,
      success: true,
      message: 'Success',
    };
  }

  // Prepare matrix (column-major)
  const aData = prepareMatrix(A);

  // For real matrices, unitary = orthogonal
  // Compute Q^H * Q = Q^T * Q and check if it equals I
  let maxDeviation = 0;
  let unitary = true;

  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      // Compute (Q^H * Q)[i,j] = sum_k conj(Q[k,i]) * Q[k,j]
      // For real matrices: sum_k Q[k,i] * Q[k,j]
      let sum = 0;
      for (let k = 0; k < n; k++) {
        sum += aData[k + i * n] * aData[k + j * n];
      }

      // Expected value: 1 if i==j, 0 otherwise
      const expected = i === j ? 1 : 0;
      const deviation = Math.abs(sum - expected);

      if (deviation > maxDeviation) {
        maxDeviation = deviation;
      }

      if (deviation > tol) {
        unitary = false;
      }
    }
  }

  return {
    isUnitary: unitary,
    maxDeviation,
    n,
    success: true,
    message: 'Success',
  };
}
