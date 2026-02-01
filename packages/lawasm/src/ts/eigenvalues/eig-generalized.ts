/**
 * Generalized Eigenvalue Problem
 *
 * Solve the generalized eigenvalue problem Ax = λBx.
 *
 * Note: DGGEV is not currently exported from LAPACK WASM.
 * This provides a pure TypeScript implementation using reduction to standard form.
 */

import type { Matrix, EigGeneralizedOptions, EigGeneralizedResult } from './types.js';
import { getMatrixDimensions, prepareMatrix } from '../helpers.js';
import { cholesky } from '../factorizations/cholesky.js';
import { eig } from './eig.js';

/**
 * Solve the generalized eigenvalue problem Ax = λBx.
 *
 * Computes eigenvalues λ and eigenvectors x such that A*x = λ*B*x.
 *
 * This implementation requires B to be symmetric positive definite.
 * It reduces the problem to standard form using Cholesky factorization:
 * If B = L*L^T, then Ax = λBx becomes (L^-1 * A * L^-T) * y = λ * y
 * where x = L^-T * y.
 *
 * @param A - First matrix (n × n)
 * @param B - Second matrix (n × n), must be symmetric positive definite
 * @param options - Computation options
 * @returns Generalized eigenvalues and eigenvectors
 *
 * @example
 * ```typescript
 * const A = [[4, 1], [1, 3]];
 * const B = [[2, 0], [0, 1]];
 * const { values, vectors, success } = eigGeneralized(A, B);
 * // Computes eigenvalues λ such that A*x = λ*B*x
 * ```
 */
export function eigGeneralized(
  A: Matrix,
  B: Matrix,
  options: EigGeneralizedOptions = {}
): EigGeneralizedResult {
  const { computeVectors = true } = options;

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
        alphar: new Float64Array(0),
        alphai: new Float64Array(0),
        beta: new Float64Array(0),
        values: [],
        n,
        info: -1,
        success: false,
        message: 'Matrix B must be symmetric positive definite for Cholesky-based reduction',
      };
    }

    const L = cholResult.factor;

    // Prepare A in column-major format
    const aData = prepareMatrix(A);

    // Compute C = L^-1 * A * L^-T
    // First solve L * Y = A to get Y = L^-1 * A
    const Y = new Float64Array(n * n);
    for (let j = 0; j < n; j++) {
      // Solve L * y_j = a_j (forward substitution)
      for (let i = 0; i < n; i++) {
        let sum = aData[j * n + i];
        for (let k = 0; k < i; k++) {
          sum -= L[k * n + i] * Y[j * n + k];
        }
        Y[j * n + i] = sum / L[i * n + i];
      }
    }

    // Now compute C = Y * L^-T
    // Solve Z * L^T = Y, or equivalently L * Z^T = Y^T
    const C = new Float64Array(n * n);
    for (let i = 0; i < n; i++) {
      // For each row of C, solve L * z = y (where y is the i-th row of Y)
      for (let j = 0; j < n; j++) {
        let sum = Y[j * n + i]; // Y[i, j] in row-major = Y[j*n + i] in column-major
        for (let k = 0; k < j; k++) {
          sum -= L[k * n + j] * C[k * n + i];
        }
        C[j * n + i] = sum / L[j * n + j];
      }
    }

    // Solve standard eigenvalue problem for C
    const eigResult = eig(C, { computeVectors });

    if (!eigResult.success) {
      return {
        alphar: new Float64Array(0),
        alphai: new Float64Array(0),
        beta: new Float64Array(0),
        values: [],
        n,
        info: eigResult.info,
        success: false,
        message: eigResult.message,
      };
    }

    // The eigenvalues are the same
    // The eigenvectors need to be transformed: x = L^-T * y
    let vectors: Float64Array | undefined;
    if (computeVectors && eigResult.vectors) {
      vectors = new Float64Array(n * n);
      const eigVecs = eigResult.vectors;

      for (let j = 0; j < n; j++) {
        // Solve L^T * x_j = y_j (back substitution)
        for (let i = n - 1; i >= 0; i--) {
          let sum = eigVecs[j * n + i];
          for (let k = i + 1; k < n; k++) {
            sum -= L[i * n + k] * vectors[j * n + k];
          }
          vectors[j * n + i] = sum / L[i * n + i];
        }
      }
    }

    // For generalized eigenvalues, we represent as α/β where β = 1
    const alphar = eigResult.wr;
    const alphai = eigResult.wi;
    const beta = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      beta[i] = 1.0;
    }

    return {
      alphar,
      alphai,
      beta,
      values: eigResult.values,
      vectors,
      n,
      info: 0,
      success: true,
      message: 'Success',
    };
  } catch (e) {
    return {
      alphar: new Float64Array(0),
      alphai: new Float64Array(0),
      beta: new Float64Array(0),
      values: [],
      n,
      info: -1,
      success: false,
      message: `Error: ${e instanceof Error ? e.message : String(e)}`,
    };
  }
}
