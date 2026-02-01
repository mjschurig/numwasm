/**
 * Singularity Test
 *
 * Test if a matrix is singular (non-invertible).
 */

import { svdvals } from '../svd/svdvals.js';
import { getMatrixDimensions } from '../helpers.js';
import type { Matrix, IsSingularOptions, IsSingularResult } from './types.js';

/**
 * Test if a matrix is singular (non-invertible).
 *
 * A matrix is singular if its determinant is zero, or equivalently,
 * if it has a zero singular value.
 *
 * @param A - Input square matrix (n × n)
 * @param options - Test options
 * @returns Whether the matrix is singular
 *
 * @example
 * ```typescript
 * // Singular matrix (linearly dependent columns)
 * const A = [[1, 2], [2, 4]];
 * const { isSingular } = isSingular(A);
 * // isSingular = true
 *
 * // Non-singular (invertible)
 * const B = [[1, 2], [3, 4]];
 * const { isSingular: sing } = isSingular(B);
 * // sing = false
 * ```
 */
export function isSingular(A: Matrix, options: IsSingularOptions = {}): IsSingularResult {
  // Get dimensions
  const [m, n] = getMatrixDimensions(A);

  // Must be square
  if (m !== n) {
    return {
      isSingular: true, // Non-square matrices are not invertible
      minSingularValue: 0,
      rcond: 0,
      n: m,
      success: true,
      message: 'Matrix is not square',
    };
  }

  // Handle edge case
  if (n === 0) {
    return {
      isSingular: false, // Empty matrix is vacuously invertible
      minSingularValue: Infinity,
      rcond: 1,
      n: 0,
      success: true,
      message: 'Success',
    };
  }

  // Compute singular values
  const svdResult = svdvals(A);
  if (!svdResult.success) {
    return {
      isSingular: false,
      minSingularValue: NaN,
      rcond: NaN,
      n,
      success: false,
      message: svdResult.message,
    };
  }

  const s = svdResult.s;
  const sMax = s[0];
  const sMin = s[s.length - 1];

  // Determine tolerance
  const eps = 2.220446049250313e-16; // Machine epsilon for double
  const defaultTol = n * eps * sMax;
  const tol = options.tol !== undefined ? options.tol : defaultTol;

  // Compute reciprocal condition number
  const rcond = sMax > 0 ? sMin / sMax : 0;

  // Matrix is singular if smallest singular value is below tolerance
  const singular = sMin <= tol;

  return {
    isSingular: singular,
    minSingularValue: sMin,
    rcond,
    n,
    success: true,
    message: 'Success',
  };
}
