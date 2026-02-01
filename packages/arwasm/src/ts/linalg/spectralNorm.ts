/**
 * SPECTRALNORM - Compute Spectral Norm (2-Norm)
 *
 * The spectral norm (also called 2-norm or operator norm) of a matrix A
 * is its largest singular value: ||A||_2 = σ_max(A)
 *
 * This equals sqrt(λ_max(A^T*A)) for real matrices.
 */

import { eigs } from '../core/eigs.js';
import type { MatVecFunction } from '../high-level-types.js';

/**
 * Options for spectral norm computation.
 */
export interface SpectralNormOptions {
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
 * Compute the spectral norm (2-norm) of a matrix.
 *
 * The spectral norm is the largest singular value, which equals
 * the square root of the largest eigenvalue of A^T*A.
 *
 * @param matvec - Function computing y = A*x (m×n matrix times n-vector)
 * @param matvecT - Function computing y = A^T*x (n×m matrix times m-vector)
 * @param m - Number of rows
 * @param n - Number of columns
 * @param options - Solver options
 * @returns The spectral norm ||A||_2
 *
 * @example
 * ```ts
 * import { spectralNorm } from 'arwasm';
 *
 * const m = 1000, n = 500;
 * const matvec = (x: Float64Array): Float64Array => { ... };
 * const matvecT = (x: Float64Array): Float64Array => { ... };
 *
 * const norm = await spectralNorm(matvec, matvecT, m, n);
 * console.log('Matrix 2-norm:', norm);
 * ```
 */
export async function spectralNorm(
  matvec: MatVecFunction,
  matvecT: MatVecFunction,
  m: number,
  n: number,
  options?: SpectralNormOptions
): Promise<number> {
  const {
    tol = 0,
    maxiter,
    ncv,
  } = options ?? {};

  // Helper to convert result to Float64Array
  const toFloat64 = (arr: Float64Array | number[]): Float64Array => {
    return arr instanceof Float64Array ? arr : new Float64Array(arr);
  };

  // Work with the smaller dimension for efficiency
  // If n <= m: compute λ_max(A^T*A)
  // If m < n:  compute λ_max(A*A^T)
  const useATA = n <= m;
  const dim = useATA ? n : m;

  const symmetricMatvec: MatVecFunction = useATA
    ? (x: Float64Array) => {
        const Ax = toFloat64(matvec(x));
        return matvecT(Ax);
      }
    : (x: Float64Array) => {
        const ATx = toFloat64(matvecT(x));
        return matvec(ATx);
      };

  // Find largest eigenvalue
  const result = await eigs(symmetricMatvec, dim, 1, {
    which: 'LM',
    tol,
    maxiter,
    ncv,
    return_eigenvectors: false,
  });

  if (!result.success || result.nconv === 0) {
    throw new Error(`Failed to compute spectral norm: ${result.message}`);
  }

  // Spectral norm = sqrt(largest eigenvalue of A^T*A)
  const eigval = result.eigenvalues[0];
  return eigval > 0 ? Math.sqrt(eigval) : 0;
}
