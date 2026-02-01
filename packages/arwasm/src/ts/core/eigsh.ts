/**
 * EIGSH - Unified Eigenvalue Solver API
 *
 * Provides a single entry point that dispatches to the appropriate
 * solver based on problem type.
 */

import { eigs } from './eigs.js';
import { eign } from './eign.js';
import type {
  MatVecFunction,
  EigshOptions,
  EigsResult,
  EignResult,
} from '../high-level-types.js';
import type { WhichSymmetric, WhichNonSymmetric } from '../types.js';

/**
 * Compute eigenvalues and eigenvectors of a sparse matrix.
 *
 * This is the unified API that dispatches to the appropriate solver:
 * - 'symmetric' (default): Uses Lanczos method (more efficient for symmetric matrices)
 * - 'nonsymmetric': Uses Arnoldi method (handles non-symmetric matrices, complex eigenvalues)
 *
 * @param matvec - Function computing y = A*x
 * @param n - Dimension of the matrix (A is n×n)
 * @param nev - Number of eigenvalues to compute
 * @param options - Solver options including problem type
 * @returns Eigenvalues and eigenvectors
 *
 * @example
 * ```ts
 * import { eigsh } from 'arwasm';
 *
 * // Symmetric matrix - uses faster Lanczos method
 * const symResult = await eigsh(symmetricMatvec, 1000, 10, {
 *   type: 'symmetric',
 *   which: 'SM',  // Smallest magnitude
 * });
 *
 * // Non-symmetric matrix - uses Arnoldi method
 * const nsymResult = await eigsh(nonsymmetricMatvec, 1000, 10, {
 *   type: 'nonsymmetric',
 *   which: 'LR',  // Largest real part
 * });
 * ```
 */
export async function eigsh(
  matvec: MatVecFunction,
  n: number,
  nev: number,
  options?: EigshOptions
): Promise<EigsResult | EignResult> {
  const { type = 'symmetric', ...rest } = options ?? {};

  if (type === 'symmetric') {
    return eigs(matvec, n, nev, {
      which: (rest.which as WhichSymmetric) ?? 'LM',
      tol: rest.tol,
      ncv: rest.ncv,
      maxiter: rest.maxiter,
      v0: rest.v0,
      return_eigenvectors: rest.return_eigenvectors,
      sigma: rest.sigma,
      OPinv: rest.OPinv,
      Bmatvec: rest.Bmatvec,
      mode: rest.mode as 1 | 2 | 3 | 4 | 5,
    });
  } else {
    return eign(matvec, n, nev, {
      which: (rest.which as WhichNonSymmetric) ?? 'LM',
      tol: rest.tol,
      ncv: rest.ncv,
      maxiter: rest.maxiter,
      v0: rest.v0,
      return_eigenvectors: rest.return_eigenvectors,
      sigma: rest.sigma,
      sigmai: rest.sigmai,
      OPinv: rest.OPinv,
      Bmatvec: rest.Bmatvec,
      mode: rest.mode as 1 | 2 | 3 | 4,
    });
  }
}

/**
 * Type guard to check if a result is from the symmetric solver.
 *
 * @param result - Result from eigsh
 * @returns true if result is from eigs (symmetric solver)
 *
 * @example
 * ```ts
 * const result = await eigsh(matvec, n, nev, options);
 * if (isEigsResult(result)) {
 *   // result.eigenvalues is available (real eigenvalues)
 * } else {
 *   // result.eigenvaluesReal and result.eigenvaluesImag are available
 * }
 * ```
 */
export function isEigsResult(
  result: EigsResult | EignResult
): result is EigsResult {
  return 'eigenvalues' in result;
}

/**
 * Type guard to check if a result is from the non-symmetric solver.
 *
 * @param result - Result from eigsh
 * @returns true if result is from eign (non-symmetric solver)
 */
export function isEignResult(
  result: EigsResult | EignResult
): result is EignResult {
  return 'eigenvaluesReal' in result;
}
