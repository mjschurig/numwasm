/**
 * NUCLEARNORM - Approximate Nuclear Norm
 *
 * The nuclear norm (also called trace norm or Schatten 1-norm) is the
 * sum of singular values: ||A||_* = Σ σ_i
 *
 * This is useful for:
 * - Matrix completion and low-rank approximation
 * - Convex relaxation of rank minimization
 * - Regularization in machine learning
 */

import { svds } from '../svd/svds.js';
import type { MatVecFunction } from '../high-level-types.js';

/**
 * Options for nuclear norm approximation.
 */
export interface NuclearNormOptions {
  /**
   * Convergence tolerance.
   * @default 0 (machine precision)
   */
  tol?: number;

  /**
   * Maximum iterations.
   */
  maxiter?: number;

  /**
   * Number of Lanczos vectors.
   */
  ncv?: number;
}

/**
 * Result from nuclear norm approximation.
 */
export interface NuclearNormResult {
  /**
   * Approximate nuclear norm (sum of computed singular values).
   */
  norm: number;

  /**
   * The k computed singular values.
   */
  singularValues: Float64Array;

  /**
   * Whether this is a lower bound (true) or exact (k = min(m,n)).
   */
  isLowerBound: boolean;

  /**
   * Number of singular values computed.
   */
  nconv: number;

  /**
   * Whether computation was successful.
   */
  success: boolean;

  /**
   * Status message.
   */
  message: string;
}

/**
 * Approximate the nuclear norm by summing the k largest singular values.
 *
 * For matrices where most singular values are small, summing the largest k
 * gives a good approximation. This is a lower bound on the true nuclear norm.
 *
 * @param matvec - Function computing y = A*x (m×n matrix times n-vector)
 * @param matvecT - Function computing y = A^T*x (n×m matrix times m-vector)
 * @param m - Number of rows
 * @param n - Number of columns
 * @param k - Number of singular values to compute
 * @param options - Solver options
 * @returns Approximate nuclear norm
 *
 * @example
 * ```ts
 * import { nuclearNormApprox } from 'arwasm';
 *
 * const m = 1000, n = 500;
 * const matvec = (x: Float64Array): Float64Array => { ... };
 * const matvecT = (x: Float64Array): Float64Array => { ... };
 *
 * // Approximate nuclear norm using top 20 singular values
 * const result = await nuclearNormApprox(matvec, matvecT, m, n, 20);
 * console.log('Nuclear norm (lower bound):', result.norm);
 * ```
 */
export async function nuclearNormApprox(
  matvec: MatVecFunction,
  matvecT: MatVecFunction,
  m: number,
  n: number,
  k: number,
  options?: NuclearNormOptions
): Promise<NuclearNormResult> {
  const {
    tol = 0,
    maxiter,
    ncv,
  } = options ?? {};

  const minDim = Math.min(m, n);

  // Clamp k to valid range
  const actualK = Math.min(k, minDim - 1);
  if (actualK < 1) {
    throw new Error('k must be at least 1');
  }

  const result = await svds(matvec, matvecT, m, n, actualK, {
    which: 'LM',
    tol,
    maxiter,
    ncv,
    return_singular_vectors: false,
  });

  if (!result.success) {
    return {
      norm: NaN,
      singularValues: new Float64Array(0),
      isLowerBound: true,
      nconv: 0,
      success: false,
      message: result.message,
    };
  }

  // Sum singular values
  const norm = result.s.reduce((sum, s) => sum + s, 0);

  return {
    norm,
    singularValues: result.s,
    isLowerBound: result.nconv < minDim,
    nconv: result.nconv,
    success: true,
    message: result.nconv < minDim
      ? `Approximation using ${result.nconv} of ${minDim} singular values`
      : 'Exact nuclear norm computed',
  };
}
