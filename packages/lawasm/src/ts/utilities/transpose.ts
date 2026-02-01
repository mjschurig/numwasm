/**
 * Matrix Transpose
 *
 * Compute the transpose of a matrix.
 */

import {
  prepareMatrix,
  getMatrixDimensions,
} from '../helpers.js';
import type { Matrix, TransposeResult } from './types.js';

/**
 * Compute the transpose of a matrix.
 *
 * @param A - Input matrix (m × n)
 * @returns Transposed matrix (n × m)
 *
 * @example
 * ```typescript
 * const A = [[1, 2, 3], [4, 5, 6]]; // 2×3
 * const { T } = transpose(A);
 * // T is 3×2: [[1, 4], [2, 5], [3, 6]] in column-major
 * ```
 */
export function transpose(A: Matrix): TransposeResult {
  // Get dimensions
  const [m, n] = getMatrixDimensions(A);

  // Handle edge cases
  if (m === 0 || n === 0) {
    return {
      T: new Float64Array(0),
      m: n, // Transposed dimensions
      n: m,
      success: true,
      message: 'Success',
    };
  }

  // Prepare matrix (column-major)
  const aData = prepareMatrix(A);

  // Create transposed matrix (n × m in column-major)
  const result = new Float64Array(m * n);

  // Transpose: T[j,i] = A[i,j]
  // In column-major:
  // A[i,j] is at aData[i + j*m]
  // T[j,i] is at result[j + i*n]
  for (let i = 0; i < m; i++) {
    for (let j = 0; j < n; j++) {
      result[j + i * n] = aData[i + j * m];
    }
  }

  return {
    T: result,
    m: n, // Result has n rows
    n: m, // Result has m columns
    success: true,
    message: 'Success',
  };
}
