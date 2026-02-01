/**
 * SVD-based Rank Computation
 *
 * Estimate the numerical rank of a matrix using SVD.
 */

import type { Matrix, SVDRankOptions, SVDRankResult } from './types.js';
import { svdvals } from './svdvals.js';

/**
 * Estimate the numerical rank of a matrix using SVD.
 *
 * The rank is estimated by counting the number of singular values above
 * a tolerance threshold. This is the most reliable method for determining
 * numerical rank.
 *
 * @param A - Input matrix (m × n). Can be 2D array (row-major) or 1D array (column-major).
 * @param options - Options including tolerance
 * @returns Rank estimate and singular values
 *
 * @example
 * ```typescript
 * // Nearly rank-deficient matrix
 * const A = [[1, 2, 3], [2, 4, 6], [1, 1, 1]];
 * const { rank, s } = svdRank(A);
 * // rank = 2 (third row is nearly a linear combination)
 * ```
 *
 * @example
 * ```typescript
 * // With custom tolerance
 * const { rank } = svdRank(A, { tol: 1e-6 });
 * ```
 */
export function svdRank(A: Matrix, options: SVDRankOptions = {}): SVDRankResult {
  const { algorithm = 'gesdd' } = options;

  // Compute singular values
  const svdResult = svdvals(A, { algorithm });

  if (!svdResult.success) {
    return {
      rank: 0,
      s: new Float64Array(0),
      tol: 0,
      m: svdResult.m,
      n: svdResult.n,
      info: svdResult.info,
      success: false,
      message: svdResult.message,
    };
  }

  const { s, m, n } = svdResult;

  // Determine tolerance
  let tol: number;
  if (options.tol !== undefined) {
    tol = options.tol;
  } else {
    // Default: max(m,n) * eps * max(s)
    const maxS = s.length > 0 ? s[0] : 0; // Singular values are in descending order
    tol = Math.max(m, n) * Number.EPSILON * maxS;
  }

  // Count singular values above tolerance
  let rank = 0;
  for (let i = 0; i < s.length; i++) {
    if (s[i] > tol) {
      rank++;
    } else {
      break; // Singular values are sorted in descending order
    }
  }

  return {
    rank,
    s,
    tol,
    m,
    n,
    info: 0,
    success: true,
    message: 'Success',
  };
}
