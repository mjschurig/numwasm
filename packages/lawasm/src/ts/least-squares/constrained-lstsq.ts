/**
 * Equality-Constrained Least Squares
 *
 * Solve: minimize ||Ax - b||_2 subject to Cx = d
 *
 * Note: DGGLSE is not currently exported from LAPACK WASM.
 * This provides a pure TypeScript implementation using the null space method.
 */

import type { Matrix, Vector, ConstrainedLstSqOptions, ConstrainedLstSqResult } from './types.js';
import { getMatrixDimensions, prepareMatrix, prepareVector } from '../helpers.js';
import { qr } from '../factorizations/qr.js';
import { lstsq } from './lstsq.js';

/**
 * Solve an equality-constrained least squares problem.
 *
 * minimize ||Ax - b||_2  subject to Cx = d
 *
 * This is useful when you have hard constraints that must be satisfied exactly,
 * along with a least squares objective.
 *
 * Requirements:
 * - C must have full row rank (p <= n)
 * - The system Cx = d must be consistent
 *
 * @param A - Objective matrix (m × n)
 * @param b - Objective vector (m)
 * @param C - Constraint matrix (p × n)
 * @param d - Constraint vector (p)
 * @param options - Computation options
 * @returns Constrained least squares solution
 *
 * @example
 * ```typescript
 * // minimize ||Ax - b||^2 subject to x1 + x2 = 1
 * const A = [[1, 0], [0, 1], [1, 1]];
 * const b = [0.5, 0.5, 1.5];
 * const C = [[1, 1]];
 * const d = [1];
 * const { x } = constrainedLstSq(A, b, C, d);
 * ```
 */
export function constrainedLstSq(
  A: Matrix,
  b: Vector,
  C: Matrix,
  d: Vector,
  _options: ConstrainedLstSqOptions = {}
): ConstrainedLstSqResult {
  // Get dimensions
  const [m, n] = getMatrixDimensions(A);
  const [p, nC] = getMatrixDimensions(C);

  if (nC !== n) {
    throw new Error(`Dimension mismatch: A has ${n} columns but C has ${nC} columns`);
  }

  const bVec = prepareVector(b);
  const dVec = prepareVector(d);

  if (bVec.length !== m) {
    throw new Error(`Dimension mismatch: A has ${m} rows but b has ${bVec.length} elements`);
  }
  if (dVec.length !== p) {
    throw new Error(`Dimension mismatch: C has ${p} rows but d has ${dVec.length} elements`);
  }

  if (p > n) {
    return {
      x: new Float64Array(0),
      residualNorm: 0,
      n,
      m,
      p,
      info: -1,
      success: false,
      message: 'Too many constraints: p > n',
    };
  }

  try {
    // Method: Null space approach
    // 1. Find a particular solution xp satisfying Cx = d
    // 2. Find the null space Z of C
    // 3. The general solution is x = xp + Z*y for any y
    // 4. Minimize ||A(xp + Zy) - b||^2 over y

    // QR factorization of C^T to get null space
    // C^T = Q * R => C = R^T * Q^T
    const cData = prepareMatrix(C);
    const cT = new Float64Array(n * p);
    for (let i = 0; i < p; i++) {
      for (let j = 0; j < n; j++) {
        cT[i * n + j] = cData[j * p + i];
      }
    }

    const qrResult = qr(cT, { mode: 'complete' });
    if (!qrResult.success || !qrResult.Q) {
      return {
        x: new Float64Array(0),
        residualNorm: 0,
        n,
        m,
        p,
        info: -1,
        success: false,
        message: 'QR factorization of C^T failed',
      };
    }

    const Q = qrResult.Q; // n×n orthogonal matrix
    const R = qrResult.R; // n×p upper triangular (only first p rows are non-zero)

    // Check rank of C (via diagonal of R)
    for (let i = 0; i < p; i++) {
      if (Math.abs(R[i * n + i]) < 1e-10) {
        return {
          x: new Float64Array(0),
          residualNorm: 0,
          n,
          m,
          p,
          info: -2,
          success: false,
          message: 'Constraint matrix C is rank deficient',
        };
      }
    }

    // Solve R^T * y1 = d for y1 (p elements)
    // R^T is p×n lower triangular (first p columns of R transposed)
    const y1 = new Float64Array(p);
    for (let i = 0; i < p; i++) {
      let sum = dVec[i];
      for (let j = 0; j < i; j++) {
        // R^T[i,j] = R[j,i] (but R is stored column-major: R[j*n + i])
        sum -= R[j * n + i] * y1[j];
      }
      y1[i] = sum / R[i * n + i];
    }

    // Particular solution: xp = Q * [y1; 0]
    const xp = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      let sum = 0;
      for (let j = 0; j < p; j++) {
        sum += Q[j * n + i] * y1[j];
      }
      xp[i] = sum;
    }

    // If n = p, we're done (no free variables)
    if (n === p) {
      // Compute residual
      const aData = prepareMatrix(A);
      let residualSq = 0;
      for (let i = 0; i < m; i++) {
        let Ax_i = 0;
        for (let j = 0; j < n; j++) {
          Ax_i += aData[j * m + i] * xp[j];
        }
        const r = Ax_i - bVec[i];
        residualSq += r * r;
      }

      return {
        x: xp,
        residualNorm: Math.sqrt(residualSq),
        n,
        m,
        p,
        info: 0,
        success: true,
        message: 'Success',
      };
    }

    // Null space: Z = last (n-p) columns of Q
    // Z is n×(n-p)
    const nullDim = n - p;
    const Z = new Float64Array(n * nullDim);
    for (let j = 0; j < nullDim; j++) {
      for (let i = 0; i < n; i++) {
        Z[j * n + i] = Q[(p + j) * n + i];
      }
    }

    // Compute AZ (m × nullDim)
    const aData = prepareMatrix(A);
    const AZ = new Float64Array(m * nullDim);
    for (let j = 0; j < nullDim; j++) {
      for (let i = 0; i < m; i++) {
        let sum = 0;
        for (let k = 0; k < n; k++) {
          sum += aData[k * m + i] * Z[j * n + k];
        }
        AZ[j * m + i] = sum;
      }
    }

    // Compute b - A*xp
    const bMinusAxp = new Float64Array(m);
    for (let i = 0; i < m; i++) {
      let Axp_i = 0;
      for (let j = 0; j < n; j++) {
        Axp_i += aData[j * m + i] * xp[j];
      }
      bMinusAxp[i] = bVec[i] - Axp_i;
    }

    // Solve least squares: minimize ||AZ*y - (b - A*xp)||^2
    const lsResult = lstsq(AZ, bMinusAxp);
    if (!lsResult.success) {
      return {
        x: new Float64Array(0),
        residualNorm: 0,
        n,
        m,
        p,
        info: -3,
        success: false,
        message: 'Least squares solve failed: ' + lsResult.message,
      };
    }

    // Final solution: x = xp + Z*y
    const y = lsResult.x;
    const x = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      x[i] = xp[i];
      for (let j = 0; j < nullDim; j++) {
        x[i] += Z[j * n + i] * y[j];
      }
    }

    // Compute residual ||Ax - b||
    let residualSq = 0;
    for (let i = 0; i < m; i++) {
      let Ax_i = 0;
      for (let j = 0; j < n; j++) {
        Ax_i += aData[j * m + i] * x[j];
      }
      const r = Ax_i - bVec[i];
      residualSq += r * r;
    }

    return {
      x,
      residualNorm: Math.sqrt(residualSq),
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
      residualNorm: 0,
      n,
      m,
      p,
      info: -1,
      success: false,
      message: `Error: ${e instanceof Error ? e.message : String(e)}`,
    };
  }
}
