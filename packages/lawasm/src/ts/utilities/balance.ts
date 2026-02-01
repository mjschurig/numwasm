/**
 * Matrix Balancing
 *
 * Balance a matrix to improve eigenvalue computation.
 */

import {
  prepareMatrix,
  getMatrixDimensions,
} from '../helpers.js';
import type { Matrix, BalanceResult } from './types.js';

/**
 * Balance a matrix to improve eigenvalue computation.
 *
 * Balancing uses diagonal similarity transformations to make the rows
 * and columns of the matrix as close in norm as possible. This can
 * improve the accuracy of eigenvalue computations.
 *
 * This is a pure TypeScript implementation of matrix balancing.
 *
 * @param A - Input square matrix (n × n)
 * @returns Balanced matrix and scaling information
 *
 * @example
 * ```typescript
 * const A = [[1, 1000], [0.001, 1]];
 * const { B, scale } = balance(A);
 * // B is the balanced matrix
 * // scale contains the diagonal scaling factors
 * ```
 */
export function balance(A: Matrix): BalanceResult {
  // Get dimensions
  const [m, n] = getMatrixDimensions(A);

  // Check squareness
  if (m !== n) {
    throw new Error(`Matrix must be square, got ${m}×${n}`);
  }

  // Handle edge case
  if (n === 0) {
    return {
      B: new Float64Array(0),
      scale: new Float64Array(0),
      ilo: 0,
      ihi: 0,
      n: 0,
      info: 0,
      success: true,
      message: 'Success',
    };
  }

  // Prepare matrix (column-major)
  const aData = prepareMatrix(A);

  // Copy to result
  const bData = new Float64Array(n * n);
  for (let i = 0; i < n * n; i++) {
    bData[i] = aData[i];
  }

  // Initialize scaling factors
  const scale = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    scale[i] = 1.0;
  }

  // Simple balancing algorithm:
  // Iteratively scale rows and columns to equalize their norms
  const RADIX = 2.0;
  const RADIX_SQ = RADIX * RADIX;
  const maxIter = 100;
  const eps = 1e-10;

  let converged = false;

  for (let iter = 0; iter < maxIter && !converged; iter++) {
    converged = true;

    for (let i = 0; i < n; i++) {
      // Compute row and column norms (excluding diagonal)
      let rowNorm = 0;
      let colNorm = 0;

      for (let j = 0; j < n; j++) {
        if (i !== j) {
          // B[i,j] in column-major: bData[i + j*n]
          colNorm += Math.abs(bData[j + i * n]); // Column i (elements in column i)
          rowNorm += Math.abs(bData[i + j * n]); // Row i (elements in row i)
        }
      }

      // Skip if row or column is zero
      if (rowNorm < eps || colNorm < eps) {
        continue;
      }

      // Find scaling factor that equalizes row and column norms
      // We want to find f such that rowNorm/f ≈ colNorm*f
      // This gives f² = rowNorm/colNorm, f = sqrt(rowNorm/colNorm)
      let f = 1.0;
      let s = colNorm + rowNorm;

      // Scale by powers of 2 for numerical stability (LAPACK approach)
      while (colNorm < rowNorm / RADIX_SQ) {
        f *= RADIX;
        colNorm *= RADIX_SQ;
      }
      while (colNorm > rowNorm * RADIX_SQ) {
        f /= RADIX;
        colNorm /= RADIX_SQ;
      }

      // Check if scaling improves balance
      if ((colNorm + rowNorm) < 0.95 * s * f) {
        converged = false;

        // Update scaling factor
        scale[i] *= f;

        // Scale row i by 1/f
        for (let j = 0; j < n; j++) {
          bData[i + j * n] /= f;
        }

        // Scale column i by f
        for (let j = 0; j < n; j++) {
          bData[j + i * n] *= f;
        }
      }
    }
  }

  return {
    B: bData,
    scale,
    ilo: 1, // 1-based, whole matrix is balanced
    ihi: n,
    n,
    info: 0,
    success: true,
    message: 'Success',
  };
}
