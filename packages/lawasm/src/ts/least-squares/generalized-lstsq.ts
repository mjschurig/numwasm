/**
 * Generalized Least Squares
 *
 * Solve: minimize ||y||_2 subject to Ax + By = d
 *
 * Note: DGGGLM is not currently exported from LAPACK WASM.
 * This provides a pure TypeScript implementation.
 */

import type { Matrix, Vector, GeneralizedLstSqOptions, GeneralizedLstSqResult } from './types.js';
import { getMatrixDimensions, prepareMatrix, prepareVector } from '../helpers.js';
import { qr } from '../factorizations/qr.js';

/**
 * Solve a generalized least squares problem.
 *
 * minimize ||y||_2  subject to Ax + By = d
 *
 * This is also known as the Gauss-Markov linear model problem.
 * It arises when fitting a model y = Ax with correlated errors,
 * where the error covariance is proportional to B*B^T.
 *
 * Requirements:
 * - [A B] must have full row rank (p = n + m)
 * - m >= p - n (enough degrees of freedom)
 *
 * @param A - Design matrix for x (p × n)
 * @param B - Design matrix for y (p × m)
 * @param d - Target vector (p)
 * @param options - Computation options
 * @returns Generalized least squares solution
 *
 * @example
 * ```typescript
 * // minimize ||y||^2 subject to x1 + y1 = 1, x1 + x2 + y2 = 2
 * const A = [[1, 0], [1, 1]];
 * const B = [[1, 0], [0, 1]];
 * const d = [1, 2];
 * const { x, y } = generalizedLstSq(A, B, d);
 * ```
 */
export function generalizedLstSq(
  A: Matrix,
  B: Matrix,
  d: Vector,
  _options: GeneralizedLstSqOptions = {}
): GeneralizedLstSqResult {
  // Get dimensions
  const [pA, n] = getMatrixDimensions(A);
  const [pB, m] = getMatrixDimensions(B);

  if (pA !== pB) {
    throw new Error(`Dimension mismatch: A has ${pA} rows but B has ${pB} rows`);
  }

  const p = pA;
  const dVec = prepareVector(d);

  if (dVec.length !== p) {
    throw new Error(`Dimension mismatch: A has ${p} rows but d has ${dVec.length} elements`);
  }

  if (n + m < p) {
    return {
      x: new Float64Array(0),
      y: new Float64Array(0),
      n,
      m,
      p,
      info: -1,
      success: false,
      message: 'System is overdetermined: n + m < p',
    };
  }

  try {
    // Method: QR factorization approach
    // Form [A B] and compute QR factorization of B
    // Then transform the problem

    const aData = prepareMatrix(A);

    // QR factorization of B: B = Q * [R; 0]
    // where Q is p×p orthogonal and R is m×m upper triangular
    const qrResult = qr(B, { mode: 'complete' });
    if (!qrResult.success || !qrResult.Q) {
      return {
        x: new Float64Array(0),
        y: new Float64Array(0),
        n,
        m,
        p,
        info: -1,
        success: false,
        message: 'QR factorization of B failed',
      };
    }

    const Q = qrResult.Q; // p×p orthogonal
    const R = qrResult.R; // k×m where k = min(p, m)

    // Transform: Q^T * d and Q^T * A
    // d_tilde = Q^T * d
    const dTilde = new Float64Array(p);
    for (let i = 0; i < p; i++) {
      let sum = 0;
      for (let j = 0; j < p; j++) {
        // Q^T[i,j] = Q[j,i] (Q is column-major: Q[i*p + j] is Q[j,i])
        sum += Q[i * p + j] * dVec[j];
      }
      dTilde[i] = sum;
    }

    // A_tilde = Q^T * A (p×n)
    const aTilde = new Float64Array(p * n);
    for (let j = 0; j < n; j++) {
      for (let i = 0; i < p; i++) {
        let sum = 0;
        for (let k = 0; k < p; k++) {
          sum += Q[i * p + k] * aData[j * p + k];
        }
        aTilde[j * p + i] = sum;
      }
    }

    // The transformed problem is:
    // Q^T * B * y = [R; 0] * y
    // So: [R; 0] * y + A_tilde * x = d_tilde
    //
    // This splits into:
    // R * y + A_tilde[0:m, :] * x = d_tilde[0:m]     (m equations)
    // A_tilde[m:p, :] * x = d_tilde[m:p]              (p-m equations)

    // If p - m >= n, the second set of equations determines x (possibly overdetermined)
    // If p - m < n, the system is underdetermined

    const numConstraints = p - m;

    if (numConstraints < 0) {
      return {
        x: new Float64Array(0),
        y: new Float64Array(0),
        n,
        m,
        p,
        info: -1,
        success: false,
        message: 'Invalid dimensions: p < m',
      };
    }

    // Extract the bottom (p-m)×n part of A_tilde and corresponding d_tilde
    // A2 = A_tilde[m:p, :]
    // d2 = d_tilde[m:p]
    let x: Float64Array;

    if (numConstraints === 0) {
      // No constraints on x, set x = 0
      x = new Float64Array(n);
    } else {
      // Solve A2 * x = d2
      const A2 = new Float64Array(numConstraints * n);
      const d2 = new Float64Array(numConstraints);

      for (let j = 0; j < n; j++) {
        for (let i = 0; i < numConstraints; i++) {
          A2[j * numConstraints + i] = aTilde[j * p + (m + i)];
        }
      }
      for (let i = 0; i < numConstraints; i++) {
        d2[i] = dTilde[m + i];
      }

      if (numConstraints >= n) {
        // Overdetermined or exact: use QR for least squares
        const qrA2 = qr(A2, { mode: 'reduced' });
        if (!qrA2.success || !qrA2.Q) {
          return {
            x: new Float64Array(0),
            y: new Float64Array(0),
            n,
            m,
            p,
            info: -2,
            success: false,
            message: 'Failed to solve for x',
          };
        }

        // Q^T * d2
        const Qd = new Float64Array(qrA2.k);
        for (let i = 0; i < qrA2.k; i++) {
          let sum = 0;
          for (let j = 0; j < numConstraints; j++) {
            sum += qrA2.Q![i * numConstraints + j] * d2[j];
          }
          Qd[i] = sum;
        }

        // Back substitution: R * x = Qd
        x = new Float64Array(n);
        for (let i = n - 1; i >= 0; i--) {
          let sum = Qd[i];
          for (let j = i + 1; j < n; j++) {
            sum -= qrA2.R[j * qrA2.k + i] * x[j];
          }
          const diag = qrA2.R[i * qrA2.k + i];
          if (Math.abs(diag) < 1e-10) {
            return {
              x: new Float64Array(0),
              y: new Float64Array(0),
              n,
              m,
              p,
              info: -2,
              success: false,
              message: 'Matrix A2 is singular',
            };
          }
          x[i] = sum / diag;
        }
      } else {
        // Underdetermined: minimum norm solution
        // For simplicity, we'll return an error for now
        return {
          x: new Float64Array(0),
          y: new Float64Array(0),
          n,
          m,
          p,
          info: -3,
          success: false,
          message: 'Underdetermined case (p - m < n) not fully implemented',
        };
      }
    }

    // Now solve for y from: R * y = d_tilde[0:m] - A_tilde[0:m, :] * x
    // First compute rhs = d_tilde[0:m] - A1 * x
    const rhs = new Float64Array(m);
    for (let i = 0; i < m; i++) {
      rhs[i] = dTilde[i];
      for (let j = 0; j < n; j++) {
        rhs[i] -= aTilde[j * p + i] * x[j];
      }
    }

    // Back substitution: R * y = rhs
    // R is k×m where k = min(p, m), stored column-major
    const k = qrResult.k;
    const y = new Float64Array(m);

    // R is upper triangular, solve by back substitution
    for (let i = Math.min(k, m) - 1; i >= 0; i--) {
      let sum = rhs[i];
      for (let j = i + 1; j < m; j++) {
        sum -= R[j * k + i] * y[j];
      }
      const diag = R[i * k + i];
      if (Math.abs(diag) < 1e-10) {
        return {
          x: new Float64Array(0),
          y: new Float64Array(0),
          n,
          m,
          p,
          info: -4,
          success: false,
          message: 'Matrix B has linearly dependent columns',
        };
      }
      y[i] = sum / diag;
    }

    return {
      x,
      y,
      n,
      m,
      p,
      info: 0,
      success: true,
      message: 'Success',
    };
  } catch (e) {
    return {
      x: new Float64Array(0),
      y: new Float64Array(0),
      n,
      m,
      p,
      info: -1,
      success: false,
      message: `Error: ${e instanceof Error ? e.message : String(e)}`,
    };
  }
}
