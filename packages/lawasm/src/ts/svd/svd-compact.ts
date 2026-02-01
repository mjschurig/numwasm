/**
 * Compact/Thin SVD
 *
 * Compute the economy-size SVD of a matrix.
 */

import type { Matrix, SVDResult } from './types.js';
import { svd } from './svd.js';

/**
 * Compute the compact (economy-size) SVD of a general m×n matrix A.
 *
 * A = U * S * V^T
 *
 * This is equivalent to `svd(A, { mode: 'reduced' })` but with a simpler interface.
 *
 * For an m×n matrix with k = min(m, n):
 * - U is m×k
 * - S is k×k diagonal (returned as vector s of length k)
 * - Vt is k×n
 *
 * @param A - Input matrix (m × n). Can be 2D array (row-major) or 1D array (column-major).
 * @param algorithm - Which algorithm to use: 'gesvd' (more accurate) or 'gesdd' (faster)
 * @returns Compact SVD result with U, s, and Vt
 *
 * @example
 * ```typescript
 * const A = [[1, 2, 3], [4, 5, 6]]; // 2×3 matrix
 * const { U, s, Vt } = svdCompact(A);
 * // U is 2×2, s has 2 elements, Vt is 2×3
 * ```
 */
export function svdCompact(
  A: Matrix,
  algorithm: 'gesvd' | 'gesdd' = 'gesdd'
): SVDResult {
  return svd(A, {
    mode: 'reduced',
    algorithm,
  });
}
