/**
 * CAYLEYEIGS - Cayley Transform Mode Eigenvalue Solver
 *
 * Solves the generalized eigenvalue problem using the Cayley transform:
 *   A * x = λ * B * x
 *
 * The Cayley transform is:
 *   OP = inv(A - sigma*B) * (A + sigma*B)
 *
 * This mode (ARPACK mode 5) is useful when eigenvalues on both sides
 * of sigma are wanted, as it maps:
 * - λ < sigma to |nu| > 1
 * - λ > sigma to |nu| < 1
 *
 * The relationship is: lambda = sigma * (1 + nu) / (1 - nu)
 */

import { eigs } from '../core/eigs.js';
import type { MatVecFunction, EigsResult } from '../high-level-types.js';
import type { WhichSymmetric } from '../types.js';

/**
 * Options for Cayley transform eigenvalue computation.
 */
export interface CayleyEigsOptions {
  /**
   * Which eigenvalues to compute in the transformed space.
   * - 'LM': Eigenvalues of original problem < sigma (|nu| > 1 in transformed)
   * - 'SM': Eigenvalues of original problem > sigma (|nu| < 1 in transformed)
   * @default 'LM'
   */
  which?: WhichSymmetric;

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
 * Compute eigenvalues using the Cayley transform.
 *
 * Solves A*x = λ*B*x using the Cayley transformation, which is useful
 * when you want eigenvalues on both sides of a shift sigma.
 *
 * The Cayley operator is: OP = inv(A - sigma*B) * (A + sigma*B)
 *
 * @param Amatvec - Matrix operator A*x
 * @param Bmatvec - Matrix operator B*x (B must be symmetric positive definite)
 * @param solveCayleyMinus - Solves (A - sigma*B)*y = x, returns y
 * @param n - Matrix dimension
 * @param nev - Number of eigenvalues to compute
 * @param sigma - Shift value
 * @param options - Solver options
 * @returns Eigenvalues and eigenvectors
 *
 * @example
 * ```ts
 * // Find eigenvalues around sigma = 5
 * const result = await cayleyEigs(
 *   Amatvec,
 *   Bmatvec,
 *   solveAmSB,    // Solves (A - 5*B)*y = x
 *   n,
 *   6,
 *   5.0
 * );
 * ```
 */
export async function cayleyEigs(
  Amatvec: MatVecFunction,
  Bmatvec: MatVecFunction,
  solveCayleyMinus: MatVecFunction,
  n: number,
  nev: number,
  sigma: number,
  options?: CayleyEigsOptions
): Promise<EigsResult> {
  const {
    which = 'LM',
    tol = 0,
    ncv,
    maxiter,
    v0,
    return_eigenvectors = true,
  } = options ?? {};

  // Mode 5 (Cayley): OP = inv(A - sigma*B) * (A + sigma*B)
  // For the reverse communication, we need:
  // y = inv(A - sigma*B) * (A + sigma*B) * x
  const OPinv: MatVecFunction = (x: Float64Array): Float64Array => {
    // Compute (A + sigma*B) * x
    const Ax = Amatvec(x);
    const Bx = Bmatvec(x);
    const AxArr = Ax instanceof Float64Array ? Ax : new Float64Array(Ax);
    const BxArr = Bx instanceof Float64Array ? Bx : new Float64Array(Bx);

    const AplusSigmaB_x = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      AplusSigmaB_x[i] = AxArr[i] + sigma * BxArr[i];
    }

    // Solve (A - sigma*B) * y = (A + sigma*B) * x
    const y = solveCayleyMinus(AplusSigmaB_x);
    return y instanceof Float64Array ? y : new Float64Array(y);
  };

  // Call eigs with Cayley mode
  const result = await eigs(Amatvec, n, nev, {
    which,
    tol,
    ncv,
    maxiter,
    v0,
    return_eigenvectors,
    sigma,
    OPinv,
    Bmatvec,
    mode: 5, // Cayley mode
  });

  return result;
}

/**
 * Compute eigenvalues in an interval using Cayley transform.
 *
 * Uses multiple Cayley transforms to find eigenvalues in a given interval
 * [a, b] by placing the shift at the midpoint.
 *
 * @param Amatvec - Matrix operator A*x
 * @param Bmatvec - Matrix operator B*x
 * @param createSolver - Function that creates a solver for (A - sigma*B)*y = x
 * @param n - Matrix dimension
 * @param nev - Number of eigenvalues to find
 * @param interval - Interval [a, b] to search for eigenvalues
 * @param options - Solver options
 * @returns Eigenvalues in the interval
 */
export async function eigsInInterval(
  Amatvec: MatVecFunction,
  Bmatvec: MatVecFunction,
  createSolver: (sigma: number) => MatVecFunction,
  n: number,
  nev: number,
  interval: [number, number],
  options?: Omit<CayleyEigsOptions, 'which'>
): Promise<EigsResult> {
  const [a, b] = interval;
  const sigma = (a + b) / 2; // Midpoint of interval

  const solver = createSolver(sigma);

  // Use 'SM' to get eigenvalues closest to sigma (inside the interval)
  const result = await cayleyEigs(
    Amatvec,
    Bmatvec,
    solver,
    n,
    nev,
    sigma,
    { ...options, which: 'SM' }
  );

  if (!result.success || !result.eigenvectors) {
    return result;
  }

  // Filter eigenvalues to only those in the interval
  const inInterval: number[] = [];
  for (let i = 0; i < result.eigenvalues.length; i++) {
    const lambda = result.eigenvalues[i];
    if (lambda >= a && lambda <= b) {
      inInterval.push(i);
    }
  }

  if (inInterval.length === 0) {
    return {
      eigenvalues: new Float64Array(0),
      eigenvectors: [],
      niter: result.niter,
      nops: result.nops,
      nconv: 0,
      info: 0,
      success: true,
      message: `No eigenvalues found in interval [${a}, ${b}]`,
    };
  }

  // Extract eigenvalues and eigenvectors in the interval
  const eigenvalues = new Float64Array(inInterval.length);
  const eigenvectors: Float64Array[] = [];

  for (let i = 0; i < inInterval.length; i++) {
    const idx = inInterval[i];
    eigenvalues[i] = result.eigenvalues[idx];
    if (result.eigenvectors) {
      eigenvectors.push(result.eigenvectors[idx]);
    }
  }

  return {
    eigenvalues,
    eigenvectors: eigenvectors.length > 0 ? eigenvectors : undefined,
    niter: result.niter,
    nops: result.nops,
    nconv: inInterval.length,
    info: 0,
    success: true,
    message: `Found ${inInterval.length} eigenvalues in interval [${a}, ${b}]`,
  };
}
