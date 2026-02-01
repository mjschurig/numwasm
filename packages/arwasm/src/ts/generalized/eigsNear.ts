/**
 * EIGSNEAR - Find Eigenvalues Near a Target (Symmetric)
 *
 * Convenience function for finding eigenvalues of a symmetric matrix
 * that are closest to a specified target value (sigma).
 *
 * Uses shift-invert mode internally, which transforms the problem
 * A*x = λ*x into inv(A - σI)*x = θ*x where λ = σ + 1/θ.
 * The eigenvalues closest to σ become the largest magnitude θ values.
 */

import { eigs } from '../core/eigs.js';
import type {
  MatVecFunction,
  OperatorFunction,
  EigsNearOptions,
  EigsResult,
} from '../high-level-types.js';

/**
 * Find eigenvalues of a symmetric matrix nearest to a target value.
 *
 * This is a convenience wrapper for shift-invert mode that handles
 * the setup automatically. You provide a function to solve (A - σI)x = b
 * and this function finds eigenvalues closest to σ.
 *
 * @param matvec - Function computing y = A*x where A is symmetric
 * @param solveShifted - Function solving (A - σI)x = b, i.e., returns inv(A - σI)*b
 * @param n - Dimension of the matrix
 * @param nev - Number of eigenvalues to find
 * @param sigma - Target value - finds eigenvalues closest to this
 * @param options - Additional solver options
 * @returns Eigenvalues nearest to sigma and their eigenvectors
 *
 * @example
 * ```ts
 * import { eigsNear } from 'arwasm';
 * import { splu } from 'superluwasm'; // or another sparse solver
 *
 * const n = 1000;
 *
 * // Your symmetric matrix A
 * const matvec = (x: Float64Array): Float64Array => {
 *   // Compute y = A*x
 *   const y = new Float64Array(n);
 *   // ...
 *   return y;
 * };
 *
 * // Pre-factor A - sigma*I for efficient solves
 * const sigma = 1.5;
 * const shiftedMatrix = buildShiftedMatrix(sigma); // A - sigma*I
 * const lu = splu(shiftedMatrix);
 *
 * const solveShifted = (b: Float64Array): Float64Array => {
 *   return lu.solve(b);
 * };
 *
 * // Find 10 eigenvalues closest to 1.5
 * const result = await eigsNear(matvec, solveShifted, n, 10, sigma);
 * console.log('Eigenvalues near 1.5:', result.eigenvalues);
 * ```
 */
export async function eigsNear(
  matvec: MatVecFunction,
  solveShifted: OperatorFunction,
  n: number,
  nev: number,
  sigma: number,
  options?: EigsNearOptions
): Promise<EigsResult> {
  // Validate inputs
  if (n < 1) {
    throw new Error('Matrix dimension n must be positive');
  }
  if (nev < 1 || nev >= n) {
    throw new Error(`nev must satisfy 0 < nev < n (got nev=${nev}, n=${n})`);
  }

  const {
    tol = 0,
    ncv,
    maxiter = n * 10,
    v0,
    return_eigenvectors = true,
  } = options ?? {};

  // Use shift-invert mode (mode 3)
  // In shift-invert, we compute eigenvalues of inv(A - sigma*I)
  // The largest magnitude eigenvalues of this operator correspond
  // to eigenvalues of A closest to sigma
  return eigs(matvec, n, nev, {
    which: 'LM', // Largest magnitude of inv(A - sigma*I)
    tol,
    ncv,
    maxiter,
    v0,
    return_eigenvectors,
    sigma,
    OPinv: solveShifted,
    mode: 3,
  });
}
