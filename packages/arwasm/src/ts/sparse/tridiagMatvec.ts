/**
 * Tridiagonal Matrix-Vector Product
 *
 * Creates a matvec function from a tridiagonal matrix.
 * Efficient O(n) implementation for this common matrix structure.
 */

import type { MatVecFunction } from '../high-level-types.js';

/**
 * Create a matrix-vector product function from a tridiagonal matrix.
 *
 * A tridiagonal matrix has non-zero entries only on the main diagonal
 * and the immediately adjacent diagonals:
 *
 * ```
 * [d0  u0   0   0  ]
 * [l0  d1  u1   0  ]
 * [ 0  l1  d2  u2  ]
 * [ 0   0  l2  d3  ]
 * ```
 *
 * @param lower - Sub-diagonal entries (length n-1)
 * @param diag - Main diagonal entries (length n)
 * @param upper - Super-diagonal entries (length n-1)
 * @returns Function computing y = A*x
 *
 * @example
 * ```ts
 * // Tridiagonal matrix (1D Laplacian-like):
 * // [[ 2, -1,  0],
 * //  [-1,  2, -1],
 * //  [ 0, -1,  2]]
 * const lower = new Float64Array([-1, -1]);
 * const diag = new Float64Array([2, 2, 2]);
 * const upper = new Float64Array([-1, -1]);
 *
 * const matvec = tridiagMatvec(lower, diag, upper);
 * const x = new Float64Array([1, 2, 1]);
 * const y = matvec(x); // [0, 0, 0] (this is the Fiedler vector for path graph)
 * ```
 */
export function tridiagMatvec(
  lower: Float64Array,
  diag: Float64Array,
  upper: Float64Array
): MatVecFunction {
  const n = diag.length;

  // Validate inputs
  if (lower.length !== n - 1) {
    throw new Error(`lower length must be n-1 (${n - 1}), got ${lower.length}`);
  }
  if (upper.length !== n - 1) {
    throw new Error(`upper length must be n-1 (${n - 1}), got ${upper.length}`);
  }

  return (x: Float64Array): Float64Array => {
    if (x.length !== n) {
      throw new Error(`Input vector length must be ${n}, got ${x.length}`);
    }

    const y = new Float64Array(n);

    // First row: y[0] = d[0]*x[0] + u[0]*x[1]
    y[0] = diag[0] * x[0];
    if (n > 1) {
      y[0] += upper[0] * x[1];
    }

    // Middle rows: y[i] = l[i-1]*x[i-1] + d[i]*x[i] + u[i]*x[i+1]
    for (let i = 1; i < n - 1; i++) {
      y[i] = lower[i - 1] * x[i - 1] + diag[i] * x[i] + upper[i] * x[i + 1];
    }

    // Last row: y[n-1] = l[n-2]*x[n-2] + d[n-1]*x[n-1]
    if (n > 1) {
      y[n - 1] = lower[n - 2] * x[n - 2] + diag[n - 1] * x[n - 1];
    }

    return y;
  };
}

/**
 * Create a symmetric tridiagonal matrix-vector product function.
 *
 * For symmetric tridiagonal matrices where lower = upper.
 * Uses only two arrays instead of three.
 *
 * @param diag - Main diagonal entries (length n)
 * @param offdiag - Off-diagonal entries (length n-1)
 * @returns Function computing y = A*x
 *
 * @example
 * ```ts
 * // Symmetric tridiagonal (1D Laplacian):
 * const diag = new Float64Array([2, 2, 2]);
 * const offdiag = new Float64Array([-1, -1]);
 *
 * const matvec = symTridiagMatvec(diag, offdiag);
 * ```
 */
export function symTridiagMatvec(
  diag: Float64Array,
  offdiag: Float64Array
): MatVecFunction {
  return tridiagMatvec(offdiag, diag, offdiag);
}

/**
 * Create a Toeplitz tridiagonal matrix-vector product function.
 *
 * For constant-coefficient tridiagonal matrices where all diagonal
 * entries are the same.
 *
 * @param lower - Sub-diagonal value
 * @param diagVal - Main diagonal value
 * @param upper - Super-diagonal value
 * @param n - Matrix dimension
 * @returns Function computing y = A*x
 *
 * @example
 * ```ts
 * // Standard 1D Laplacian: -1, 2, -1
 * const matvec = toeplitzTridiagMatvec(-1, 2, -1, 100);
 * ```
 */
export function toeplitzTridiagMatvec(
  lower: number,
  diagVal: number,
  upper: number,
  n: number
): MatVecFunction {
  return (x: Float64Array): Float64Array => {
    if (x.length !== n) {
      throw new Error(`Input vector length must be ${n}, got ${x.length}`);
    }

    const y = new Float64Array(n);

    // First row
    y[0] = diagVal * x[0];
    if (n > 1) {
      y[0] += upper * x[1];
    }

    // Middle rows
    for (let i = 1; i < n - 1; i++) {
      y[i] = lower * x[i - 1] + diagVal * x[i] + upper * x[i + 1];
    }

    // Last row
    if (n > 1) {
      y[n - 1] = lower * x[n - 2] + diagVal * x[n - 1];
    }

    return y;
  };
}
