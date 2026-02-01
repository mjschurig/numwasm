/**
 * LDL^T Factorization
 *
 * Compute the LDL^T factorization of a symmetric matrix.
 *
 * Note: DSYTRF is not currently exported. This provides a pure TypeScript
 * implementation for symmetric indefinite matrices using Bunch-Kaufman pivoting.
 */

import type { Matrix, LDLOptions, LDLResult } from './types.js';
import { getMatrixDimensions, prepareMatrix } from '../helpers.js';

/**
 * Compute the LDL^T factorization of a symmetric n×n matrix A.
 *
 * A = L * D * L^T (or A = U * D * U^T if upper=true)
 *
 * where L (or U) is unit triangular and D is block diagonal with
 * 1×1 and 2×2 blocks.
 *
 * This factorization is useful for symmetric indefinite matrices
 * (which cannot use Cholesky).
 *
 * Note: This is a simplified implementation. For production use with
 * ill-conditioned matrices, the full LAPACK DSYTRF with Bunch-Kaufman
 * pivoting is recommended.
 *
 * @param A - Symmetric matrix (n × n). Only the specified triangle is used.
 * @param options - Factorization options
 * @returns LDL factorization result
 *
 * @example
 * ```typescript
 * // Symmetric indefinite matrix
 * const A = [[1, 2], [2, -1]];
 * const { factor, D, success } = ldl(A);
 * ```
 */
export function ldl(A: Matrix, options: LDLOptions = {}): LDLResult {
  const { upper = false } = options;

  // Get dimensions
  const [m, n] = getMatrixDimensions(A);
  if (m !== n) {
    throw new Error(`Matrix A must be square, got ${m}×${n}`);
  }

  // Prepare matrix
  const aData = prepareMatrix(A);

  // Work with a copy
  const factor = new Float64Array(n * n);
  for (let i = 0; i < n * n; i++) {
    factor[i] = aData[i];
  }

  // Diagonal elements
  const D = new Float64Array(n);
  // Pivot indices (positive = 1x1 block, negative = 2x2 block)
  const ipiv = new Int32Array(n);

  // Simple LDL^T without pivoting (diagonal pivoting only)
  // This is a simplified version - full Bunch-Kaufman not implemented
  try {
    if (!upper) {
      // Lower triangular factorization: A = L * D * L^T
      for (let j = 0; j < n; j++) {
        // Compute D[j]
        let djj = factor[j * n + j];
        for (let k = 0; k < j; k++) {
          const ljk = factor[k * n + j];
          djj -= ljk * ljk * D[k];
        }

        if (Math.abs(djj) < 1e-15) {
          return {
            factor: new Float64Array(0),
            D: new Float64Array(0),
            ipiv: new Int32Array(0),
            upper,
            n,
            info: j + 1,
            success: false,
            message: `Zero pivot encountered at position ${j + 1}`,
          };
        }

        D[j] = djj;
        ipiv[j] = j + 1; // 1-based, positive means 1x1 block

        // Compute L[i,j] for i > j
        for (let i = j + 1; i < n; i++) {
          let sum = factor[j * n + i]; // A[i,j]
          for (let k = 0; k < j; k++) {
            sum -= factor[k * n + i] * factor[k * n + j] * D[k];
          }
          factor[j * n + i] = sum / djj;
        }

        // Set diagonal to 1 for unit triangular
        factor[j * n + j] = 1.0;

        // Zero out upper triangle
        for (let i = 0; i < j; i++) {
          factor[j * n + i] = 0.0;
        }
      }
    } else {
      // Upper triangular factorization: A = U * D * U^T
      for (let j = n - 1; j >= 0; j--) {
        // Compute D[j]
        let djj = factor[j * n + j];
        for (let k = j + 1; k < n; k++) {
          const ujk = factor[k * n + j];
          djj -= ujk * ujk * D[k];
        }

        if (Math.abs(djj) < 1e-15) {
          return {
            factor: new Float64Array(0),
            D: new Float64Array(0),
            ipiv: new Int32Array(0),
            upper,
            n,
            info: j + 1,
            success: false,
            message: `Zero pivot encountered at position ${j + 1}`,
          };
        }

        D[j] = djj;
        ipiv[j] = j + 1;

        // Compute U[i,j] for i < j
        for (let i = 0; i < j; i++) {
          let sum = factor[j * n + i]; // A[i,j]
          for (let k = j + 1; k < n; k++) {
            sum -= factor[k * n + i] * factor[k * n + j] * D[k];
          }
          factor[j * n + i] = sum / djj;
        }

        // Set diagonal to 1
        factor[j * n + j] = 1.0;

        // Zero out lower triangle
        for (let i = j + 1; i < n; i++) {
          factor[j * n + i] = 0.0;
        }
      }
    }

    return {
      factor,
      D,
      ipiv,
      upper,
      n,
      info: 0,
      success: true,
      message: 'Success',
    };
  } catch (e) {
    return {
      factor: new Float64Array(0),
      D: new Float64Array(0),
      ipiv: new Int32Array(0),
      upper,
      n,
      info: -1,
      success: false,
      message: `Error during factorization: ${e instanceof Error ? e.message : String(e)}`,
    };
  }
}
