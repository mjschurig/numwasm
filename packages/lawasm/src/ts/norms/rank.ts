/**
 * Matrix Rank
 *
 * Compute the numerical rank of a matrix using SVD.
 */

import { svdvals } from '../svd/svdvals.js';
import { getMatrixDimensions } from '../helpers.js';
import type { Matrix, RankOptions, RankResult } from './types.js';

/**
 * Compute the numerical rank of a matrix.
 *
 * The rank is computed as the number of singular values greater than
 * a tolerance threshold. By default, the tolerance is:
 *   tol = max(m, n) * eps * max(singular values)
 *
 * where eps is the machine epsilon (~2.22e-16 for double precision).
 *
 * @param A - Input matrix (m × n)
 * @param options - Computation options
 * @returns The numerical rank
 *
 * @example
 * ```typescript
 * // Full rank matrix
 * const A = [[1, 2], [3, 4]];
 * const { rank: r1 } = rank(A);
 * // r1 = 2
 *
 * // Rank deficient matrix
 * const B = [[1, 2, 3], [2, 4, 6]]; // Second row is 2x first row
 * const { rank: r2 } = rank(B);
 * // r2 = 1
 *
 * // Custom tolerance
 * const C = [[1, 0], [0, 1e-10]];
 * const { rank: r3 } = rank(C, { tol: 1e-8 });
 * // r3 = 1 (second singular value below tolerance)
 * ```
 */
export function rank(A: Matrix, options: RankOptions = {}): RankResult {
  // Get dimensions
  const [m, n] = getMatrixDimensions(A);

  // Handle edge cases
  if (m === 0 || n === 0) {
    return {
      rank: 0,
      s: new Float64Array(0),
      tol: 0,
      m,
      n,
      info: 0,
      success: true,
      message: 'Success',
    };
  }

  // Compute singular values
  const svdResult = svdvals(A);
  if (!svdResult.success) {
    return {
      rank: 0,
      s: new Float64Array(0),
      tol: NaN,
      m,
      n,
      info: svdResult.info,
      success: false,
      message: svdResult.message,
    };
  }

  const s = svdResult.s;
  const minDim = Math.min(m, n);

  // Determine tolerance
  const eps = 2.220446049250313e-16; // Machine epsilon for double
  const maxDim = Math.max(m, n);
  const sMax = s[0]; // Largest singular value

  let tol: number;
  if (options.tol !== undefined) {
    tol = options.tol;
  } else {
    // Default tolerance: max(m,n) * eps * max(s)
    tol = maxDim * eps * sMax;
  }

  // Count singular values above tolerance
  let rankValue = 0;
  for (let i = 0; i < minDim; i++) {
    if (s[i] > tol) {
      rankValue++;
    }
  }

  return {
    rank: rankValue,
    s,
    tol,
    m,
    n,
    info: 0,
    success: true,
    message: 'Success',
  };
}
