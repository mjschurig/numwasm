/**
 * Upper Triangular Part
 *
 * Extract the upper triangular part of a matrix.
 */

import {
  prepareMatrix,
  getMatrixDimensions,
} from '../helpers.js';
import type { Matrix, TriuOptions, TriuResult } from './types.js';

/**
 * Extract the upper triangular part of a matrix.
 *
 * Returns a matrix with elements below the k-th diagonal zeroed.
 *
 * @param A - Input matrix (m × n)
 * @param options - Extraction options
 * @returns Upper triangular part
 *
 * @example
 * ```typescript
 * const A = [[1, 2, 3], [4, 5, 6], [7, 8, 9]];
 *
 * // Main diagonal and above (k=0, default)
 * const { U } = triu(A);
 * // U = [[1, 2, 3], [0, 5, 6], [0, 0, 9]]
 *
 * // First superdiagonal and above (k=1)
 * const { U: U1 } = triu(A, { k: 1 });
 * // U1 = [[0, 2, 3], [0, 0, 6], [0, 0, 0]]
 *
 * // First subdiagonal and above (k=-1)
 * const { U: Um1 } = triu(A, { k: -1 });
 * // Um1 = [[1, 2, 3], [4, 5, 6], [0, 8, 9]]
 * ```
 */
export function triu(A: Matrix, options: TriuOptions = {}): TriuResult {
  const { k = 0 } = options;

  // Get dimensions
  const [m, n] = getMatrixDimensions(A);

  // Handle edge cases
  if (m === 0 || n === 0) {
    return {
      U: new Float64Array(0),
      m,
      n,
      k,
      success: true,
      message: 'Success',
    };
  }

  // Prepare matrix (column-major)
  const aData = prepareMatrix(A);

  // Create result
  const result = new Float64Array(m * n);

  // Copy upper triangular part
  // Element (i,j) is in upper triangle if j >= i + k (equivalently i <= j - k)
  for (let j = 0; j < n; j++) {
    for (let i = 0; i < m; i++) {
      if (i <= j - k) {
        // In upper triangular part
        result[i + j * m] = aData[i + j * m];
      }
      // else: leave as 0
    }
  }

  return {
    U: result,
    m,
    n,
    k,
    success: true,
    message: 'Success',
  };
}
