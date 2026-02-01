/**
 * BUCKLINGEIGS - Buckling Mode Eigenvalue Solver
 *
 * Solves the buckling eigenvalue problem:
 *   K * x = λ * Kg * x
 *
 * where K is the stiffness matrix and Kg is the geometric stiffness matrix.
 * This is ARPACK's Mode 4 for symmetric generalized problems.
 *
 * The buckling mode finds eigenvalues nearest to a shift sigma using
 * the transformation:
 *   inv(K - sigma*Kg) * K * x = nu * x
 *   lambda = sigma * nu / (nu - 1)
 */

import { eigs } from '../core/eigs.js';
import type { MatVecFunction, EigsResult } from '../high-level-types.js';

/**
 * Options for buckling eigenvalue computation.
 */
export interface BucklingEigsOptions {
  /**
   * Convergence tolerance.
   * @default 0 (machine precision)
   */
  tol?: number;

  /**
   * Number of Lanczos vectors.
   */
  ncv?: number;

  /**
   * Maximum iterations.
   */
  maxiter?: number;

  /**
   * Starting vector for the iteration.
   */
  v0?: Float64Array;

  /**
   * Whether to compute eigenvectors.
   * @default true
   */
  return_eigenvectors?: boolean;
}

/**
 * Compute buckling eigenvalues and modes.
 *
 * Solves the buckling problem K*x = λ*Kg*x where:
 * - K is the stiffness matrix (symmetric positive definite)
 * - Kg is the geometric stiffness matrix (symmetric)
 * - λ is the buckling load factor
 *
 * Uses shift-invert buckling mode (ARPACK mode 4) to find eigenvalues
 * nearest to the shift sigma.
 *
 * @param Kmatvec - Stiffness matrix operator K*x
 * @param Kgmatvec - Geometric stiffness matrix operator Kg*x
 * @param solveShifted - Solves (K - sigma*Kg)*y = x, returns y
 * @param n - Matrix dimension
 * @param nev - Number of buckling modes to compute
 * @param sigma - Shift value (find eigenvalues near sigma)
 * @param options - Solver options
 * @returns Buckling eigenvalues and modes
 *
 * @example
 * ```ts
 * // Find 6 smallest buckling load factors
 * const result = await bucklingEigs(
 *   Kmatvec,      // Stiffness
 *   Kgmatvec,     // Geometric stiffness
 *   solveKsKg,    // Solves (K - sigma*Kg)*y = x
 *   n,
 *   6,
 *   0.0           // Shift near zero for smallest eigenvalues
 * );
 * console.log('Critical buckling loads:', result.eigenvalues);
 * ```
 */
export async function bucklingEigs(
  Kmatvec: MatVecFunction,
  Kgmatvec: MatVecFunction,
  solveShifted: MatVecFunction,
  n: number,
  nev: number,
  sigma: number,
  options?: BucklingEigsOptions
): Promise<EigsResult> {
  const {
    tol = 0,
    ncv,
    maxiter,
    v0,
    return_eigenvectors = true,
  } = options ?? {};

  // Mode 4 (Buckling): OP = inv(K - sigma*Kg) * K
  // The operator for reverse communication is:
  // y = inv(K - sigma*Kg) * K * x
  const OPinv: MatVecFunction = (x: Float64Array): Float64Array => {
    // First compute K*x
    const Kx = Kmatvec(x);
    const KxArr = Kx instanceof Float64Array ? Kx : new Float64Array(Kx);
    // Then solve (K - sigma*Kg)*y = K*x
    const y = solveShifted(KxArr);
    return y instanceof Float64Array ? y : new Float64Array(y);
  };

  // Call eigs with buckling mode
  const result = await eigs(Kmatvec, n, nev, {
    which: 'LM', // Largest magnitude in transformed space
    tol,
    ncv,
    maxiter,
    v0,
    return_eigenvectors,
    sigma,
    OPinv,
    Bmatvec: Kgmatvec,
    mode: 4, // Buckling mode
  });

  return result;
}

/**
 * Result type for buckling analysis with additional information.
 */
export interface BucklingResult extends EigsResult {
  /**
   * Critical buckling load (smallest positive eigenvalue).
   */
  criticalLoad?: number;

  /**
   * Critical buckling mode (eigenvector for critical load).
   */
  criticalMode?: Float64Array;
}

/**
 * Compute critical buckling load and mode.
 *
 * Convenience function that finds the smallest positive buckling
 * eigenvalue and its corresponding mode shape.
 *
 * @param Kmatvec - Stiffness matrix operator
 * @param Kgmatvec - Geometric stiffness matrix operator
 * @param solveShifted - Solver for (K - sigma*Kg)*y = x
 * @param n - Matrix dimension
 * @param options - Solver options
 * @returns Critical buckling load and mode
 */
export async function criticalBucklingLoad(
  Kmatvec: MatVecFunction,
  Kgmatvec: MatVecFunction,
  solveShifted: MatVecFunction,
  n: number,
  options?: BucklingEigsOptions
): Promise<BucklingResult> {
  // Find several eigenvalues near zero to ensure we get the smallest positive one
  const result = await bucklingEigs(
    Kmatvec,
    Kgmatvec,
    solveShifted,
    n,
    Math.min(6, n - 2),
    0.0, // Shift at zero
    options
  );

  if (!result.success) {
    return {
      ...result,
      criticalLoad: undefined,
      criticalMode: undefined,
    };
  }

  // Find smallest positive eigenvalue
  let criticalLoad: number | undefined;
  let criticalIndex = -1;

  for (let i = 0; i < result.eigenvalues.length; i++) {
    const lambda = result.eigenvalues[i];
    if (lambda > 1e-10) { // Positive eigenvalue
      if (criticalLoad === undefined || lambda < criticalLoad) {
        criticalLoad = lambda;
        criticalIndex = i;
      }
    }
  }

  const criticalMode = criticalIndex >= 0 && result.eigenvectors
    ? result.eigenvectors[criticalIndex]
    : undefined;

  return {
    ...result,
    criticalLoad,
    criticalMode,
  };
}
