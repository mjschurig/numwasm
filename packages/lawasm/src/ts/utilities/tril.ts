/**
 * Lower Triangular Part
 *
 * Extract the lower triangular part of a matrix.
 */

import {
  prepareMatrix,
  getMatrixDimensions,
} from '../helpers.js';
import type { Matrix, TrilOptions, TrilResult } from './types.js';

/**
 * Extract the lower triangular part of a matrix.
 *
 * Returns a matrix with elements above the k-th diagonal zeroed.
 *
 * @param A - Input matrix (m × n)
 * @param options - Extraction options
 * @returns Lower triangular part
 *
 * @example
 * ```typescript
 * const A = [[1, 2, 3], [4, 5, 6], [7, 8, 9]];
 *
 * // Main diagonal and below (k=0, default)
 * const { L } = tril(A);
 * // L = [[1, 0, 0], [4, 5, 0], [7, 8, 9]]
 *
 * // First subdiagonal and below (k=-1)
 * const { L: Lm1 } = tril(A, { k: -1 });
 * // Lm1 = [[0, 0, 0], [4, 0, 0], [7, 8, 0]]
 *
 * // First superdiagonal and below (k=1)
 * const { L: L1 } = tril(A, { k: 1 });
 * // L1 = [[1, 2, 0], [4, 5, 6], [7, 8, 9]]
 * ```
 */
export function tril(A: Matrix, options: TrilOptions = {}): TrilResult {
  const { k = 0 } = options;

  // Get dimensions
  const [m, n] = getMatrixDimensions(A);

  // Handle edge cases
  if (m === 0 || n === 0) {
    return {
      L: new Float64Array(0),
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

  // Copy lower triangular part
  // Element (i,j) is in lower triangle if i >= j - k (equivalently j <= i + k)
  for (let j = 0; j < n; j++) {
    for (let i = 0; i < m; i++) {
      if (i >= j - k) {
        // In lower triangular part
        result[i + j * m] = aData[i + j * m];
      }
      // else: leave as 0
    }
  }

  return {
    L: result,
    m,
    n,
    k,
    success: true,
    message: 'Success',
  };
}
