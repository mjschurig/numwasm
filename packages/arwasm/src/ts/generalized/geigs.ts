/**
 * GEIGS - Generalized Eigenvalue Solver
 *
 * Solves the generalized eigenvalue problem A*x = λ*B*x
 * where A is symmetric and B is symmetric positive definite.
 *
 * This is a convenience wrapper around eigs() that handles the
 * generalized problem setup more explicitly.
 */

import { eigs } from '../core/eigs.js';
import type {
  MatVecFunction,
  GeigsOptions,
  EigsResult,
} from '../high-level-types.js';

/**
 * Compute eigenvalues and eigenvectors of a generalized eigenvalue problem.
 *
 * Solves A*x = λ*B*x where:
 * - A is a symmetric matrix
 * - B is a symmetric positive definite matrix
 *
 * @param Amatvec - Function computing y = A*x
 * @param Bmatvec - Function computing y = B*x (B must be symmetric positive definite)
 * @param n - Dimension of the matrices (A, B are n×n)
 * @param nev - Number of eigenvalues to compute
 * @param options - Solver options
 * @returns Eigenvalues and eigenvectors
 *
 * @example
 * ```ts
 * import { geigs } from 'arwasm';
 *
 * // Generalized eigenvalue problem with mass matrix
 * const n = 100;
 *
 * // Stiffness matrix (tridiagonal)
 * const Amatvec = (x: Float64Array): Float64Array => {
 *   const y = new Float64Array(n);
 *   for (let i = 0; i < n; i++) {
 *     y[i] = 2 * x[i];
 *     if (i > 0) y[i] -= x[i - 1];
 *     if (i < n - 1) y[i] -= x[i + 1];
 *   }
 *   return y;
 * };
 *
 * // Mass matrix (diagonal for simplicity)
 * const Bmatvec = (x: Float64Array): Float64Array => {
 *   const y = new Float64Array(n);
 *   for (let i = 0; i < n; i++) {
 *     y[i] = (1 + 0.1 * i / n) * x[i]; // Varying mass
 *   }
 *   return y;
 * };
 *
 * const result = await geigs(Amatvec, Bmatvec, n, 6, { which: 'SM' });
 * console.log('Smallest generalized eigenvalues:', result.eigenvalues);
 * ```
 */
export async function geigs(
  Amatvec: MatVecFunction,
  Bmatvec: MatVecFunction,
  n: number,
  nev: number,
  options?: GeigsOptions
): Promise<EigsResult> {
  // Validate inputs
  if (n < 1) {
    throw new Error('Matrix dimension n must be positive');
  }
  if (nev < 1 || nev >= n) {
    throw new Error(`nev must satisfy 0 < nev < n (got nev=${nev}, n=${n})`);
  }

  const {
    which = 'LM',
    tol = 0,
    ncv,
    maxiter = n * 10,
    v0,
    return_eigenvectors = true,
    sigma,
    OPinv,
    mode: userMode,
  } = options ?? {};

  // Determine mode
  // Mode 2: Regular mode with B inner product - computes inv(B)*A*x
  // Mode 3: Shift-invert mode - computes inv(A - sigma*B)*B*x
  let mode = userMode ?? 2;

  // If sigma is provided and no mode specified, use shift-invert (mode 3)
  if (sigma !== undefined && !userMode) {
    mode = 3;
  }

  // For mode 3, OPinv is required
  if (mode === 3 && !OPinv) {
    throw new Error(
      'OPinv function required for shift-invert mode (mode 3). ' +
        'It should compute y = inv(A - sigma*B) * B * x'
    );
  }

  // Call the underlying eigs with generalized problem setup
  return eigs(Amatvec, n, nev, {
    which,
    tol,
    ncv,
    maxiter,
    v0,
    return_eigenvectors,
    sigma,
    OPinv,
    Bmatvec,
    mode,
  });
}
