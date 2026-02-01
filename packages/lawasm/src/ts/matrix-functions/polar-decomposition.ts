/**
 * Polar Decomposition
 *
 * Compute A = U*P where U is unitary and P is positive semidefinite.
 */

import { svd } from '../svd/svd.js';
import { getMatrixDimensions } from '../helpers.js';
import type { Matrix, PolarDecompositionResult } from './types.js';

/**
 * Compute the polar decomposition A = U*P.
 *
 * U is unitary (orthogonal for real matrices) and P is positive semidefinite.
 * Uses SVD: if A = W*S*Vt, then U = W*Vt and P = V*S*Vt.
 *
 * @param A - Input matrix (m × n)
 * @returns Polar decomposition A = U*P
 *
 * @example
 * ```typescript
 * const A = [[1, 2], [3, 4]];
 * const { U, P } = polarDecomposition(A);
 * // U is orthogonal, P is symmetric positive semidefinite
 * // A = U * P
 * ```
 */
export function polarDecomposition(A: Matrix): PolarDecompositionResult {
  const [m, n] = getMatrixDimensions(A);

  if (m === 0 || n === 0) {
    return {
      U: new Float64Array(0),
      P: new Float64Array(0),
      m,
      n,
      success: true,
      message: 'Success',
    };
  }

  // Compute SVD: A = W * S * Vt
  const svdResult = svd(A, { mode: 'full' });
  if (!svdResult.success) {
    return {
      U: new Float64Array(0),
      P: new Float64Array(0),
      m,
      n,
      success: false,
      message: svdResult.message,
    };
  }

  const W = svdResult.U!;    // m × m (or m × k)
  const s = svdResult.s;     // min(m,n) singular values
  const Vt = svdResult.Vt!;  // n × n (or k × n)

  const k = Math.min(m, n);

  // U = W * Vt (m × n matrix for rectangular A, or m × m for square)
  // For polar decomposition of A (m × n), U is m × n and P is n × n

  // Compute U = W[:, :k] * Vt[:k, :]
  const U = new Float64Array(m * n);
  for (let i = 0; i < m; i++) {
    for (let j = 0; j < n; j++) {
      let sum = 0;
      for (let l = 0; l < k; l++) {
        // W[i, l] * Vt[l, j]
        // W is m × m in column-major: W[i + l*m]
        // Vt is n × n in column-major: Vt[l + j*n]
        sum += W[i + l * m] * Vt[l + j * n];
      }
      U[i + j * m] = sum;
    }
  }

  // P = V * S * Vt = (Vt)^T * S * Vt
  // P is n × n symmetric positive semidefinite
  // First compute V * S where V = Vt^T
  const VS = new Float64Array(n * k);
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < k; j++) {
      // V[i, j] = Vt[j, i] = Vt[j + i*n]
      VS[i + j * n] = Vt[j + i * n] * s[j];
    }
  }

  // Then compute (V*S) * Vt
  const P = new Float64Array(n * n);
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      let sum = 0;
      for (let l = 0; l < k; l++) {
        // VS[i, l] * Vt[l, j]
        sum += VS[i + l * n] * Vt[l + j * n];
      }
      P[i + j * n] = sum;
    }
  }

  return {
    U,
    P,
    m,
    n,
    success: true,
    message: 'Success',
  };
}
