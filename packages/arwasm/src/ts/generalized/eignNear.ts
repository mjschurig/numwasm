/**
 * EIGNNEAR - Find Eigenvalues Near a Target (Non-Symmetric)
 *
 * Convenience function for finding eigenvalues of a non-symmetric matrix
 * that are closest to a specified target value (sigma).
 *
 * Uses shift-invert mode internally. For non-symmetric matrices, the
 * target sigma can be complex.
 */

import { eign } from '../core/eign.js';
import type {
  MatVecFunction,
  OperatorFunction,
  EignNearOptions,
  EignResult,
  Complex,
} from '../high-level-types.js';

/**
 * Find eigenvalues of a non-symmetric matrix nearest to a target value.
 *
 * This is a convenience wrapper for shift-invert mode. You provide a
 * function to solve (A - σI)x = b and this function finds eigenvalues
 * closest to σ.
 *
 * For non-symmetric matrices, the shift σ can be complex. When σ is
 * complex, the solver uses the real arithmetic formulation with both
 * real and imaginary shift components.
 *
 * @param matvec - Function computing y = A*x
 * @param solveShifted - Function solving (A - σI)x = b for real σ,
 *                       or the appropriate shifted system for complex σ
 * @param n - Dimension of the matrix
 * @param nev - Number of eigenvalues to find
 * @param sigma - Target value (real number or complex {re, im})
 * @param options - Additional solver options
 * @returns Eigenvalues nearest to sigma and their eigenvectors
 *
 * @example
 * ```ts
 * import { eignNear } from 'arwasm';
 *
 * const n = 1000;
 *
 * // Your non-symmetric matrix A
 * const matvec = (x: Float64Array): Float64Array => {
 *   const y = new Float64Array(n);
 *   // Compute y = A*x
 *   return y;
 * };
 *
 * // Solver for (A - sigma*I)x = b
 * const sigma = 2.5;
 * const solveShifted = (b: Float64Array): Float64Array => {
 *   // Use LU decomposition or iterative solver
 *   // to solve (A - sigma*I)x = b
 *   return x;
 * };
 *
 * // Find 10 eigenvalues closest to 2.5
 * const result = await eignNear(matvec, solveShifted, n, 10, sigma);
 * console.log('Real parts:', result.eigenvaluesReal);
 * console.log('Imag parts:', result.eigenvaluesImag);
 * ```
 *
 * @example
 * ```ts
 * // With complex shift (finds eigenvalues near 1 + 0.5i)
 * const sigmaComplex = { re: 1.0, im: 0.5 };
 * const result = await eignNear(matvec, solveShifted, n, 10, sigmaComplex);
 * ```
 */
export async function eignNear(
  matvec: MatVecFunction,
  solveShifted: OperatorFunction,
  n: number,
  nev: number,
  sigma: number | Complex,
  options?: EignNearOptions
): Promise<EignResult> {
  // Validate inputs
  if (n < 1) {
    throw new Error('Matrix dimension n must be positive');
  }
  if (nev < 1 || nev >= n - 1) {
    throw new Error(`nev must satisfy 0 < nev < n-1 (got nev=${nev}, n=${n})`);
  }

  const {
    tol = 0,
    ncv,
    maxiter = n * 10,
    v0,
    return_eigenvectors = true,
  } = options ?? {};

  // Parse sigma into real and imaginary parts
  let sigmaReal: number;
  let sigmaImag: number;

  if (typeof sigma === 'number') {
    sigmaReal = sigma;
    sigmaImag = 0;
  } else {
    sigmaReal = sigma.re;
    sigmaImag = sigma.im;
  }

  // Use shift-invert mode (mode 3 for real shift, mode 4 for complex shift)
  // Mode 3: Real shift-invert
  // Mode 4: Complex shift-invert (when sigmai != 0)
  const mode = sigmaImag === 0 ? 3 : 4;

  return eign(matvec, n, nev, {
    which: 'LM', // Largest magnitude of inv(A - sigma*I)
    tol,
    ncv,
    maxiter,
    v0,
    return_eigenvectors,
    sigma: sigmaReal,
    sigmai: sigmaImag,
    OPinv: solveShifted,
    mode,
  });
}
