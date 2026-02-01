/**
 * Conjugate Transpose (Hermitian Transpose)
 *
 * Compute the conjugate transpose A^H of a matrix.
 */

import {
  prepareMatrix,
  getMatrixDimensions,
} from '../helpers.js';
import type { Matrix, HermitianResult } from './types.js';

/**
 * Compute the conjugate transpose (Hermitian transpose) of a matrix.
 *
 * For real matrices, this is the same as the regular transpose.
 * For complex matrices, this transposes and conjugates.
 *
 * @param A - Input matrix (m × n)
 * @returns Conjugate transpose (n × m)
 *
 * @example
 * ```typescript
 * // Real matrix (same as transpose)
 * const A = [[1, 2], [3, 4]];
 * const { H } = hermitian(A);
 *
 * // For complex matrices, this computes A^H = (A^*)^T
 * ```
 */
export function hermitian(A: Matrix): HermitianResult {
  // Get dimensions
  const [m, n] = getMatrixDimensions(A);

  // Handle edge cases
  if (m === 0 || n === 0) {
    return {
      H: new Float64Array(0),
      m: n, // Transposed dimensions
      n: m,
      success: true,
      message: 'Success',
    };
  }

  // Prepare matrix (column-major)
  const aData = prepareMatrix(A);

  // Create result (n × m in column-major)
  const result = new Float64Array(m * n);

  // For real matrices, Hermitian = transpose
  // Transpose: H[j,i] = A[i,j]
  for (let i = 0; i < m; i++) {
    for (let j = 0; j < n; j++) {
      result[j + i * n] = aData[i + j * m];
    }
  }

  return {
    H: result,
    m: n, // Result has n rows
    n: m, // Result has m columns
    success: true,
    message: 'Success',
  };
}
