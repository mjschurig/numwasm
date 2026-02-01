/**
 * Symmetric Generalized Eigenvalue Problem
 *
 * Solve the symmetric generalized eigenvalue problem Ax = λBx
 * where A is symmetric and B is symmetric positive definite.
 *
 * Note: DSYGV is not currently exported from LAPACK WASM.
 * This provides a pure TypeScript implementation using Cholesky reduction.
 */

import type { Matrix, EigGeneralizedSymmetricOptions, EigGeneralizedSymmetricResult } from './types.js';
import { getMatrixDimensions, prepareMatrix } from '../helpers.js';
import { cholesky } from '../factorizations/cholesky.js';
import { eigSymmetric } from './eig-symmetric.js';

/**
 * Solve the symmetric generalized eigenvalue problem.
 *
 * Computes eigenvalues and eigenvectors of one of three forms:
 * - itype=1: A*x = λ*B*x
 * - itype=2: A*B*x = λ*x
 * - itype=3: B*A*x = λ*x
 *
 * where A is symmetric and B is symmetric positive definite.
 * All eigenvalues are guaranteed to be real.
 *
 * @param A - Symmetric matrix (n × n)
 * @param B - Symmetric positive definite matrix (n × n)
 * @param options - Computation options
 * @returns Real eigenvalues and eigenvectors
 *
 * @example
 * ```typescript
 * const A = [[4, 1], [1, 3]];
 * const B = [[2, 0], [0, 1]];
 * const { values, vectors, success } = eigGeneralizedSymmetric(A, B);
 * ```
 */
export function eigGeneralizedSymmetric(
  A: Matrix,
  B: Matrix,
  options: EigGeneralizedSymmetricOptions = {}
): EigGeneralizedSymmetricResult {
  const {
    itype = 1,
    computeVectors = true,
    uplo = 'lower',
  } = options;

  // Get dimensions
  const [mA, nA] = getMatrixDimensions(A);
  const [mB, nB] = getMatrixDimensions(B);

  if (mA !== nA) {
    throw new Error(`Matrix A must be square, got ${mA}×${nA}`);
  }
  if (mB !== nB) {
    throw new Error(`Matrix B must be square, got ${mB}×${nB}`);
  }
  if (mA !== mB) {
    throw new Error(`Matrices A and B must have the same dimensions`);
  }

  const n = mA;

  try {
    // Cholesky factorization of B: B = L * L^T
    const cholResult = cholesky(B, { upper: false });
    if (!cholResult.success) {
      return {
        values: new Float64Array(0),
        n,
        info: -1,
        success: false,
        message: 'Matrix B must be symmetric positive definite',
      };
    }

    const L = cholResult.factor;
    const aData = prepareMatrix(A);

    // Transform to standard form based on itype
    // itype=1: A*x = λ*B*x → (L^-1 * A * L^-T) * y = λ * y where x = L^-T * y
    // itype=2: A*B*x = λ*x → (L^T * A * L) * y = λ * y where x = L^-1 * y
    // itype=3: B*A*x = λ*x → (L^T * A * L) * y = λ * y where x = L * y

    let C: Float64Array;

    if (itype === 1) {
      // C = L^-1 * A * L^-T
      // First: Y = L^-1 * A (solve L * Y = A)
      const Y = new Float64Array(n * n);
      for (let j = 0; j < n; j++) {
        for (let i = 0; i < n; i++) {
          let sum = aData[j * n + i];
          for (let k = 0; k < i; k++) {
            sum -= L[k * n + i] * Y[j * n + k];
          }
          Y[j * n + i] = sum / L[i * n + i];
        }
      }

      // C = Y * L^-T (solve C * L^T = Y)
      C = new Float64Array(n * n);
      for (let i = 0; i < n; i++) {
        for (let j = 0; j < n; j++) {
          let sum = Y[j * n + i];
          for (let k = 0; k < j; k++) {
            sum -= L[k * n + j] * C[k * n + i];
          }
          C[j * n + i] = sum / L[j * n + j];
        }
      }
    } else {
      // itype = 2 or 3: C = L^T * A * L
      // First: Y = A * L
      const Y = new Float64Array(n * n);
      for (let j = 0; j < n; j++) {
        for (let i = 0; i < n; i++) {
          let sum = 0;
          for (let k = 0; k < n; k++) {
            // A[i,k] * L[k,j]
            sum += aData[k * n + i] * L[j * n + k];
          }
          Y[j * n + i] = sum;
        }
      }

      // C = L^T * Y
      C = new Float64Array(n * n);
      for (let j = 0; j < n; j++) {
        for (let i = 0; i < n; i++) {
          let sum = 0;
          for (let k = 0; k < n; k++) {
            // L^T[i,k] * Y[k,j] = L[k,i] * Y[k,j]
            sum += L[i * n + k] * Y[j * n + k];
          }
          C[j * n + i] = sum;
        }
      }
    }

    // Solve standard symmetric eigenvalue problem for C
    const eigResult = eigSymmetric(C, { computeVectors, uplo });

    if (!eigResult.success) {
      return {
        values: new Float64Array(0),
        n,
        info: eigResult.info,
        success: false,
        message: eigResult.message,
      };
    }

    // Transform eigenvectors back
    let vectors: Float64Array | undefined;
    if (computeVectors && eigResult.vectors) {
      vectors = new Float64Array(n * n);
      const eigVecs = eigResult.vectors;

      if (itype === 1) {
        // x = L^-T * y (solve L^T * x = y)
        for (let j = 0; j < n; j++) {
          for (let i = n - 1; i >= 0; i--) {
            let sum = eigVecs[j * n + i];
            for (let k = i + 1; k < n; k++) {
              sum -= L[i * n + k] * vectors[j * n + k];
            }
            vectors[j * n + i] = sum / L[i * n + i];
          }
        }
      } else if (itype === 2) {
        // x = L^-1 * y (solve L * x = y)
        for (let j = 0; j < n; j++) {
          for (let i = 0; i < n; i++) {
            let sum = eigVecs[j * n + i];
            for (let k = 0; k < i; k++) {
              sum -= L[k * n + i] * vectors[j * n + k];
            }
            vectors[j * n + i] = sum / L[i * n + i];
          }
        }
      } else {
        // itype === 3: x = L * y
        for (let j = 0; j < n; j++) {
          for (let i = 0; i < n; i++) {
            let sum = 0;
            for (let k = 0; k <= i; k++) {
              sum += L[k * n + i] * eigVecs[j * n + k];
            }
            vectors[j * n + i] = sum;
          }
        }
      }
    }

    return {
      values: eigResult.values,
      vectors,
      n,
      info: 0,
      success: true,
      message: 'Success',
    };
  } catch (e) {
    return {
      values: new Float64Array(0),
      n,
      info: -1,
      success: false,
      message: `Error: ${e instanceof Error ? e.message : String(e)}`,
    };
  }
}
