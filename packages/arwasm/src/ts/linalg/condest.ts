/**
 * CONDEST - Condition Number Estimation
 *
 * Estimates the condition number of a matrix, which is the ratio
 * of largest to smallest singular values: κ(A) = σ_max / σ_min
 *
 * A large condition number indicates the matrix is nearly singular
 * and solutions to linear systems may be sensitive to perturbations.
 */

import { eigs } from '../core/eigs.js';
import type { MatVecFunction } from '../high-level-types.js';

/**
 * Options for condition number estimation.
 */
export interface CondestOptions {
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
 * Result from condition number estimation.
 */
export interface CondestResult {
  /**
   * Estimated condition number κ(A) = σ_max / σ_min
   */
  condition: number;

  /**
   * Largest singular value σ_max
   */
  largest: number;

  /**
   * Smallest singular value σ_min
   */
  smallest: number;

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
 * Estimate the condition number of a matrix.
 *
 * Computes both the largest and smallest singular values to estimate
 * the 2-norm condition number κ(A) = σ_max / σ_min.
 *
 * @param matvec - Function computing y = A*x (m×n matrix times n-vector)
 * @param matvecT - Function computing y = A^T*x (n×m matrix times m-vector)
 * @param n - Number of columns (for square matrices, n = m)
 * @param options - Solver options
 * @returns Condition number estimate
 *
 * @example
 * ```ts
 * import { condest } from 'arwasm';
 *
 * const n = 100;
 * const matvec = (x: Float64Array): Float64Array => { ... };
 * const matvecT = (x: Float64Array): Float64Array => { ... };
 *
 * const result = await condest(matvec, matvecT, n);
 * console.log('Condition number:', result.condition);
 *
 * if (result.condition > 1e10) {
 *   console.log('Warning: Matrix is ill-conditioned');
 * }
 * ```
 */
export async function condest(
  matvec: MatVecFunction,
  matvecT: MatVecFunction,
  n: number,
  options?: CondestOptions
): Promise<CondestResult> {
  const {
    tol = 0,
    maxiter,
    ncv,
  } = options ?? {};

  const toFloat64 = (arr: Float64Array | number[]): Float64Array => {
    return arr instanceof Float64Array ? arr : new Float64Array(arr);
  };

  // Create symmetric matvec for A^T*A
  const ataMatvec: MatVecFunction = (x: Float64Array) => {
    const Ax = toFloat64(matvec(x));
    return matvecT(Ax);
  };

  // Find largest eigenvalue of A^T*A
  const largestResult = await eigs(ataMatvec, n, 1, {
    which: 'LM',
    tol,
    maxiter,
    ncv,
    return_eigenvectors: false,
  });

  if (!largestResult.success || largestResult.nconv === 0) {
    return {
      condition: Infinity,
      largest: NaN,
      smallest: NaN,
      success: false,
      message: `Failed to compute largest singular value: ${largestResult.message}`,
    };
  }

  // Find smallest eigenvalue of A^T*A
  const smallestResult = await eigs(ataMatvec, n, 1, {
    which: 'SM',
    tol,
    maxiter,
    ncv,
    return_eigenvectors: false,
  });

  if (!smallestResult.success || smallestResult.nconv === 0) {
    return {
      condition: Infinity,
      largest: Math.sqrt(Math.max(0, largestResult.eigenvalues[0])),
      smallest: NaN,
      success: false,
      message: `Failed to compute smallest singular value: ${smallestResult.message}`,
    };
  }

  const sigmaMax = Math.sqrt(Math.max(0, largestResult.eigenvalues[0]));
  const sigmaMin = Math.sqrt(Math.max(0, smallestResult.eigenvalues[0]));

  const condition = sigmaMin > 0 ? sigmaMax / sigmaMin : Infinity;

  return {
    condition,
    largest: sigmaMax,
    smallest: sigmaMin,
    success: true,
    message: 'Converged successfully',
  };
}
