/**
 * SPECTRALRADIUS - Compute Spectral Radius
 *
 * The spectral radius of a matrix A is the largest absolute value
 * of its eigenvalues: ρ(A) = max|λ_i|
 *
 * This is important for:
 * - Convergence analysis of iterative methods
 * - Stability analysis of dynamical systems
 * - Matrix norms and condition numbers
 */

import { eign } from '../core/eign.js';
import type { MatVecFunction } from '../high-level-types.js';

/**
 * Options for spectral radius computation.
 */
export interface SpectralRadiusOptions {
  /**
   * Convergence tolerance.
   * @default 0 (machine precision)
   */
  tol?: number;

  /**
   * Maximum iterations.
   * @default n * 10
   */
  maxiter?: number;

  /**
   * Number of Arnoldi vectors.
   */
  ncv?: number;
}

/**
 * Compute the spectral radius of a matrix.
 *
 * The spectral radius is the largest magnitude eigenvalue.
 *
 * @param matvec - Function computing y = A*x
 * @param n - Dimension of the matrix
 * @param options - Solver options
 * @returns The spectral radius ρ(A)
 *
 * @example
 * ```ts
 * import { spectralRadius } from 'arwasm';
 *
 * const n = 100;
 * const matvec = (x: Float64Array): Float64Array => {
 *   // Compute y = A*x
 *   return y;
 * };
 *
 * const rho = await spectralRadius(matvec, n);
 * console.log('Spectral radius:', rho);
 *
 * // Check if iteration matrix converges
 * if (rho < 1) {
 *   console.log('Iterative method will converge');
 * }
 * ```
 */
export async function spectralRadius(
  matvec: MatVecFunction,
  n: number,
  options?: SpectralRadiusOptions
): Promise<number> {
  const {
    tol = 0,
    maxiter,
    ncv,
  } = options ?? {};

  // Find the largest magnitude eigenvalue
  const result = await eign(matvec, n, 1, {
    which: 'LM',
    tol,
    maxiter,
    ncv,
    return_eigenvectors: false,
  });

  if (!result.success || result.nconv === 0) {
    throw new Error(`Failed to compute spectral radius: ${result.message}`);
  }

  // Compute magnitude of the eigenvalue
  const re = result.eigenvaluesReal[0];
  const im = result.eigenvaluesImag[0];

  return Math.sqrt(re * re + im * im);
}
