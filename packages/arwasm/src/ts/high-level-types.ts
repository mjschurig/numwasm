/**
 * ARPACK High-Level TypeScript Interfaces
 *
 * User-friendly API for computing eigenvalues and eigenvectors
 * of large sparse matrices using ARPACK's Implicitly Restarted
 * Arnoldi/Lanczos methods.
 */

import type { WhichSymmetric, WhichNonSymmetric } from './types.js';
import type { RealArray } from './helpers.js';

// Re-export for convenience
export type { RealArray } from './helpers.js';

// ============================================================
// FUNCTION TYPES
// ============================================================

/**
 * Matrix-vector product function for standard eigenvalue problems.
 *
 * Computes y = A*x where A is the matrix whose eigenvalues are sought.
 * The input x is provided as a Float64Array view into WASM memory for
 * zero-copy performance. The output can be either Float64Array or number[].
 *
 * @param x - Input vector of length n (Float64Array view)
 * @returns Output vector y = A*x of length n
 *
 * @example
 * ```ts
 * // For a tridiagonal matrix with 2 on diagonal, -1 on off-diagonals
 * const matvec: MatVecFunction = (x) => {
 *   const n = x.length;
 *   const y = new Float64Array(n);
 *   for (let i = 0; i < n; i++) {
 *     y[i] = 2 * x[i];
 *     if (i > 0) y[i] -= x[i - 1];
 *     if (i < n - 1) y[i] -= x[i + 1];
 *   }
 *   return y;
 * };
 * ```
 */
export type MatVecFunction = (x: Float64Array) => RealArray;

/**
 * Operator function for shift-invert mode.
 *
 * For shift-invert mode (finding eigenvalues near sigma):
 * - Standard problem: computes y = inv(A - sigma*I) * x
 * - Generalized problem: computes y = inv(A - sigma*B) * B * x
 *
 * This requires solving a linear system, typically using LU factorization
 * or an iterative solver.
 *
 * @param x - Input vector (Float64Array view)
 * @returns Output vector after applying the operator
 */
export type OperatorFunction = (x: Float64Array) => RealArray;

/**
 * B-matrix product function for generalized eigenvalue problems.
 *
 * For generalized problems A*x = lambda*B*x, computes y = B*x.
 * B must be symmetric positive definite.
 *
 * @param x - Input vector (Float64Array view)
 * @returns Output vector y = B*x
 */
export type BMatVecFunction = (x: Float64Array) => RealArray;

// ============================================================
// OPTIONS INTERFACES
// ============================================================

/**
 * Base options shared by all eigensolvers.
 */
export interface EigOptionsBase {
  /**
   * Convergence tolerance.
   * An eigenvalue is considered converged when the Ritz estimate
   * of the relative error is less than tol.
   * @default 0 (machine precision)
   */
  tol?: number;

  /**
   * Number of Lanczos/Arnoldi vectors to use.
   * Larger values may improve convergence but increase memory.
   * Must satisfy: nev < ncv <= n
   * @default min(n, max(2*nev+1, 20)) for non-symmetric
   * @default min(n, max(2*nev, 20)) for symmetric
   */
  ncv?: number;

  /**
   * Maximum number of Arnoldi update iterations.
   * @default n * 10
   */
  maxiter?: number;

  /**
   * Starting vector for the Arnoldi iteration.
   * If not provided, a random vector is used.
   * Can be number[] or Float64Array.
   */
  v0?: RealArray;

  /**
   * Whether to compute eigenvectors in addition to eigenvalues.
   * @default true
   */
  return_eigenvectors?: boolean;
}

/**
 * Options for symmetric eigenvalue problems (eigs).
 *
 * @example
 * ```ts
 * const options: EigsOptions = {
 *   which: 'LM',     // Largest magnitude eigenvalues
 *   tol: 1e-10,      // High precision
 *   ncv: 20,         // Use 20 Lanczos vectors
 * };
 * ```
 */
export interface EigsOptions extends EigOptionsBase {
  /**
   * Which eigenvalues to compute.
   * - 'LM': Largest Magnitude (largest |lambda|)
   * - 'SM': Smallest Magnitude (smallest |lambda|)
   * - 'LA': Largest Algebraic (most positive lambda)
   * - 'SA': Smallest Algebraic (most negative lambda)
   * - 'BE': Both Ends (half from each end)
   * @default 'LM'
   */
  which?: WhichSymmetric;

  /**
   * Shift value for shift-invert mode.
   * When provided, finds eigenvalues nearest to sigma.
   * Requires providing an OPinv function that computes
   * y = inv(A - sigma*I) * x for standard problems.
   */
  sigma?: number;

  /**
   * Operator for shift-invert mode.
   * Computes y = inv(A - sigma*I) * x (or inv(A - sigma*B) * B * x
   * for generalized problems).
   */
  OPinv?: OperatorFunction;

  /**
   * B-matrix product for generalized problems A*x = lambda*B*x.
   * If provided, solves the generalized eigenproblem.
   */
  Bmatvec?: BMatVecFunction;

  /**
   * Computational mode.
   * - 1: Standard eigenvalue problem A*x = lambda*x
   * - 2: Generalized A*x = lambda*B*x with shift-invert on B
   * - 3: Shift-invert mode: A*x = lambda*x transformed to inv(A-sigma*I)*x = theta*x
   * - 4: Buckling mode
   * - 5: Cayley mode
   * @default 1 (or 3 if sigma is provided)
   */
  mode?: 1 | 2 | 3 | 4 | 5;
}

/**
 * Options for non-symmetric eigenvalue problems (eign).
 *
 * @example
 * ```ts
 * const options: EignOptions = {
 *   which: 'LR',     // Largest real part
 *   tol: 1e-8,
 * };
 * ```
 */
export interface EignOptions extends EigOptionsBase {
  /**
   * Which eigenvalues to compute.
   * - 'LM': Largest Magnitude (largest |lambda|)
   * - 'SM': Smallest Magnitude (smallest |lambda|)
   * - 'LR': Largest Real part
   * - 'SR': Smallest Real part
   * - 'LI': Largest Imaginary part
   * - 'SI': Smallest Imaginary part
   * @default 'LM'
   */
  which?: WhichNonSymmetric;

  /**
   * Real part of shift value for shift-invert mode.
   */
  sigma?: number;

  /**
   * Imaginary part of shift (for complex shifts).
   * Only used when sigma is also provided.
   * @default 0
   */
  sigmai?: number;

  /**
   * Operator for shift-invert mode.
   */
  OPinv?: OperatorFunction;

  /**
   * B-matrix product for generalized problems.
   */
  Bmatvec?: BMatVecFunction;

  /**
   * Computational mode (1-4).
   * @default 1
   */
  mode?: 1 | 2 | 3 | 4;
}

// ============================================================
// RESULT INTERFACES
// ============================================================

/**
 * Result from symmetric eigenvalue solver (eigs).
 *
 * All eigenvalues are guaranteed to be real for symmetric matrices.
 */
export interface EigsResult {
  /**
   * Computed eigenvalues.
   * Length is equal to the number of converged eigenvalues (at most nev).
   */
  eigenvalues: Float64Array;

  /**
   * Computed eigenvectors (if return_eigenvectors=true).
   * eigenvectors[i] is the eigenvector corresponding to eigenvalues[i].
   * Each eigenvector has length n (the matrix dimension).
   */
  eigenvectors?: Float64Array[];

  /**
   * Number of Arnoldi iterations performed.
   */
  niter: number;

  /**
   * Number of matrix-vector products computed.
   */
  nops: number;

  /**
   * Number of converged eigenvalues.
   */
  nconv: number;

  /**
   * ARPACK info return code.
   * 0 = success, negative = error
   */
  info: number;

  /**
   * Whether the computation was successful.
   */
  success: boolean;

  /**
   * Human-readable status message.
   */
  message: string;
}

/**
 * Result from non-symmetric eigenvalue solver (eign).
 *
 * Eigenvalues may be complex and are returned as separate real/imaginary arrays.
 * Complex eigenvalues always appear in conjugate pairs.
 */
export interface EignResult {
  /**
   * Real parts of computed eigenvalues.
   */
  eigenvaluesReal: Float64Array;

  /**
   * Imaginary parts of computed eigenvalues.
   * If eigenvaluesImag[i] = 0, the eigenvalue is real.
   */
  eigenvaluesImag: Float64Array;

  /**
   * Computed eigenvectors (if return_eigenvectors=true).
   *
   * For real eigenvalues (eigenvaluesImag[i] = 0):
   *   eigenvectors[i] is the real eigenvector
   *
   * For complex conjugate pairs (eigenvaluesImag[i] != 0):
   *   eigenvectors[i] is the real part
   *   eigenvectors[i+1] is the imaginary part
   *   Full eigenvector for eigenvalue[i] = eigenvectors[i] + j*eigenvectors[i+1]
   *   Full eigenvector for eigenvalue[i+1] = eigenvectors[i] - j*eigenvectors[i+1]
   */
  eigenvectors?: Float64Array[];

  /**
   * Number of Arnoldi iterations performed.
   */
  niter: number;

  /**
   * Number of matrix-vector products computed.
   */
  nops: number;

  /**
   * Number of converged eigenvalues.
   */
  nconv: number;

  /**
   * ARPACK info return code.
   */
  info: number;

  /**
   * Whether the computation was successful.
   */
  success: boolean;

  /**
   * Human-readable status message.
   */
  message: string;
}

// ============================================================
// UNIFIED API TYPES
// ============================================================

/**
 * Problem type for the unified eigsh function.
 */
export type ProblemType = 'symmetric' | 'nonsymmetric';

/**
 * Unified options combining both symmetric and non-symmetric options.
 */
export interface EigshOptions extends EigOptionsBase {
  /**
   * Problem type.
   * - 'symmetric': Use symmetric solver (Lanczos) - more efficient for symmetric matrices
   * - 'nonsymmetric': Use non-symmetric solver (Arnoldi)
   * @default 'symmetric'
   */
  type?: ProblemType;

  /**
   * Which eigenvalues to compute.
   * For symmetric: 'LM', 'SM', 'LA', 'SA', 'BE'
   * For non-symmetric: 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'
   * @default 'LM'
   */
  which?: WhichSymmetric | WhichNonSymmetric;

  /**
   * Shift value (real part) for shift-invert mode.
   */
  sigma?: number;

  /**
   * Imaginary part of shift (only for non-symmetric).
   */
  sigmai?: number;

  /**
   * Operator for shift-invert mode.
   */
  OPinv?: OperatorFunction;

  /**
   * B-matrix product for generalized problems.
   */
  Bmatvec?: BMatVecFunction;

  /**
   * Computational mode.
   */
  mode?: number;
}
