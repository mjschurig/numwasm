/**
 * ARPACK High-Level TypeScript Interfaces
 *
 * User-friendly API for computing eigenvalues and eigenvectors
 * of large sparse matrices using ARPACK's Implicitly Restarted
 * Arnoldi/Lanczos methods.
 */

import type { WhichSymmetric, WhichNonSymmetric, WhichComplex } from './types.js';
import type { RealArray, ComplexArray } from './helpers.js';

// Re-export for convenience
export type { RealArray, ComplexArray } from './helpers.js';

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

// ============================================================
// COMPLEX EIGENVALUE SOLVER TYPES
// ============================================================

/**
 * Complex matrix-vector product function.
 *
 * Computes y = A*x where A is a complex matrix and x, y are complex vectors.
 * Complex arrays use the ComplexArray type with separate re and im Float64Arrays.
 *
 * @param x - Input complex vector
 * @returns Output complex vector y = A*x
 *
 * @example
 * ```ts
 * // For a complex tridiagonal matrix
 * const matvec: ComplexMatVecFunction = (x) => {
 *   const n = x.re.length;
 *   const re = new Float64Array(n);
 *   const im = new Float64Array(n);
 *   for (let i = 0; i < n; i++) {
 *     // Diagonal: 2 + 0.1i
 *     re[i] = 2 * x.re[i] - 0.1 * x.im[i];
 *     im[i] = 2 * x.im[i] + 0.1 * x.re[i];
 *     // Off-diagonals: -1
 *     if (i > 0) {
 *       re[i] -= x.re[i - 1];
 *       im[i] -= x.im[i - 1];
 *     }
 *     if (i < n - 1) {
 *       re[i] -= x.re[i + 1];
 *       im[i] -= x.im[i + 1];
 *     }
 *   }
 *   return { re, im };
 * };
 * ```
 */
export type ComplexMatVecFunction = (x: ComplexArray) => ComplexArray;

/**
 * Complex operator function for shift-invert mode.
 *
 * For shift-invert mode (finding eigenvalues near sigma):
 * - Standard problem: computes y = inv(A - sigma*I) * x
 * - Generalized problem: computes y = inv(A - sigma*B) * B * x
 *
 * @param x - Input complex vector
 * @returns Output complex vector after applying the operator
 */
export type ComplexOperatorFunction = (x: ComplexArray) => ComplexArray;

/**
 * Complex B-matrix product function for generalized eigenvalue problems.
 *
 * For generalized problems A*x = lambda*B*x, computes y = B*x.
 * B must be Hermitian positive definite for complex problems.
 *
 * @param x - Input complex vector
 * @returns Output complex vector y = B*x
 */
export type ComplexBMatVecFunction = (x: ComplexArray) => ComplexArray;

/**
 * Base options for complex eigensolvers.
 */
export interface ComplexEigOptionsBase {
  /**
   * Convergence tolerance.
   * @default 0 (machine precision)
   */
  tol?: number;

  /**
   * Number of Arnoldi vectors to use.
   * Must satisfy: nev+2 <= ncv <= n
   * @default min(n, max(2*nev+1, 20))
   */
  ncv?: number;

  /**
   * Maximum number of Arnoldi update iterations.
   * @default n * 10
   */
  maxiter?: number;

  /**
   * Whether to compute eigenvectors in addition to eigenvalues.
   * @default true
   */
  return_eigenvectors?: boolean;
}

/**
 * Options for complex eigenvalue problems (zeigs).
 *
 * @example
 * ```ts
 * const options: ZeigsOptions = {
 *   which: 'LM',     // Largest magnitude eigenvalues
 *   tol: 1e-10,      // High precision
 *   ncv: 20,         // Use 20 Arnoldi vectors
 * };
 * ```
 */
export interface ZeigsOptions extends ComplexEigOptionsBase {
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
  which?: WhichComplex;

  /**
   * Starting vector for the Arnoldi iteration (complex).
   * If not provided, a random vector is used.
   */
  v0?: ComplexArray;

  /**
   * Complex shift value for shift-invert mode.
   * When provided, finds eigenvalues nearest to sigma.
   */
  sigma?: { re: number; im: number };

  /**
   * Operator for shift-invert mode.
   * Computes y = inv(A - sigma*I) * x for standard problems.
   */
  OPinv?: ComplexOperatorFunction;

  /**
   * B-matrix product for generalized problems A*x = lambda*B*x.
   * If provided, solves the generalized eigenproblem.
   * B must be Hermitian positive definite.
   */
  Bmatvec?: ComplexBMatVecFunction;

  /**
   * Computational mode.
   * - 1: Standard eigenvalue problem A*x = lambda*x
   * - 2: Generalized A*x = lambda*B*x with shift-invert on B
   * - 3: Shift-invert mode
   * @default 1 (or 3 if sigma is provided)
   */
  mode?: 1 | 2 | 3;
}

/**
 * Options for complex Hermitian eigenvalue problems (zeigsh).
 * Hermitian matrices have real eigenvalues but complex eigenvectors.
 */
export interface ZeigshOptions extends ComplexEigOptionsBase {
  /**
   * Which eigenvalues to compute.
   * - 'LM': Largest Magnitude (largest |lambda|)
   * - 'SM': Smallest Magnitude (smallest |lambda|)
   * - 'LR': Largest Real part (same as LA for Hermitian)
   * - 'SR': Smallest Real part (same as SA for Hermitian)
   * @default 'LM'
   */
  which?: WhichComplex;

  /**
   * Starting vector for the Arnoldi iteration (complex).
   */
  v0?: ComplexArray;

  /**
   * Real shift value for shift-invert mode.
   * For Hermitian matrices, eigenvalues are real so shift is real.
   */
  sigma?: number;

  /**
   * Operator for shift-invert mode.
   */
  OPinv?: ComplexOperatorFunction;

  /**
   * B-matrix product for generalized problems.
   * B must be Hermitian positive definite.
   */
  Bmatvec?: ComplexBMatVecFunction;

  /**
   * Computational mode.
   * @default 1
   */
  mode?: 1 | 2 | 3;
}

/**
 * Result from complex eigenvalue solver (zeigs).
 *
 * All eigenvalues and eigenvectors are complex.
 */
export interface ZeigsResult {
  /**
   * Real parts of computed eigenvalues.
   */
  eigenvaluesReal: Float64Array;

  /**
   * Imaginary parts of computed eigenvalues.
   */
  eigenvaluesImag: Float64Array;

  /**
   * Computed complex eigenvectors (if return_eigenvectors=true).
   * eigenvectors[i] is the complex eigenvector for eigenvalue i.
   * Each eigenvector has re and im Float64Arrays of length n.
   */
  eigenvectors?: ComplexArray[];

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

/**
 * Result from complex Hermitian eigenvalue solver (zeigsh).
 *
 * Eigenvalues are real (since the matrix is Hermitian), but eigenvectors are complex.
 */
export interface ZeigshResult {
  /**
   * Computed eigenvalues (real for Hermitian matrices).
   */
  eigenvalues: Float64Array;

  /**
   * Computed complex eigenvectors (if return_eigenvectors=true).
   * eigenvectors[i] is the complex eigenvector for eigenvalue i.
   */
  eigenvectors?: ComplexArray[];

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
// SINGULAR VALUE DECOMPOSITION TYPES
// ============================================================

/**
 * Options for partial SVD computation (svds).
 *
 * @example
 * ```ts
 * const options: SvdsOptions = {
 *   which: 'LM',              // Largest singular values
 *   tol: 1e-10,
 *   return_singular_vectors: 'both',
 * };
 * ```
 */
export interface SvdsOptions {
  /**
   * Which singular values to compute.
   * - 'LM': Largest singular values
   * - 'SM': Smallest singular values
   * @default 'LM'
   */
  which?: 'LM' | 'SM';

  /**
   * Convergence tolerance.
   * @default 0 (machine precision)
   */
  tol?: number;

  /**
   * Number of Lanczos vectors to use.
   * @default min(min(m,n), max(2*k+1, 20))
   */
  ncv?: number;

  /**
   * Maximum number of iterations.
   * @default min(m,n) * 10
   */
  maxiter?: number;

  /**
   * Starting vector for the iteration.
   */
  v0?: RealArray;

  /**
   * Which singular vectors to compute.
   * - true or 'both': Compute both U and Vt
   * - 'u': Compute only left singular vectors U
   * - 'vh': Compute only right singular vectors Vt
   * - false: Compute only singular values
   * @default true
   */
  return_singular_vectors?: boolean | 'u' | 'vh' | 'both';
}

/**
 * Result from partial SVD computation (svds).
 */
export interface SvdsResult {
  /**
   * Left singular vectors (m x k).
   * U[i] is the i-th left singular vector of length m.
   * Only present if return_singular_vectors includes 'u'.
   */
  U?: Float64Array[];

  /**
   * Singular values (length k), sorted in descending order.
   */
  s: Float64Array;

  /**
   * Right singular vectors (k x n).
   * Vt[i] is the i-th right singular vector of length n.
   * Only present if return_singular_vectors includes 'vh'.
   */
  Vt?: Float64Array[];

  /**
   * Number of converged singular values.
   */
  nconv: number;

  /**
   * Number of iterations performed.
   */
  niter: number;

  /**
   * Number of matrix-vector products computed.
   */
  nops: number;

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
// GENERALIZED EIGENVALUE PROBLEM TYPES
// ============================================================

/**
 * Options for generalized eigenvalue problems (geigs).
 * Solves A*x = λ*B*x where B is symmetric positive definite.
 */
export interface GeigsOptions extends EigOptionsBase {
  /**
   * Which eigenvalues to compute.
   * @default 'LM'
   */
  which?: WhichSymmetric;

  /**
   * Shift value for shift-invert mode.
   * When provided, finds eigenvalues nearest to sigma.
   */
  sigma?: number;

  /**
   * Operator for shift-invert mode.
   * Must compute y = inv(A - sigma*B) * B * x
   */
  OPinv?: OperatorFunction;

  /**
   * Computational mode.
   * - 1: Regular inverse (not recommended for generalized)
   * - 2: Shift-invert on B: inv(B)*A*x = mu*x, lambda = mu
   * - 3: Shift-invert: inv(A - sigma*B)*B*x = nu*x, lambda = sigma + 1/nu
   * @default 2
   */
  mode?: 1 | 2 | 3;
}

// ============================================================
// SHIFT-INVERT CONVENIENCE TYPES
// ============================================================

/**
 * Complex number type for shift values.
 */
export interface Complex {
  re: number;
  im: number;
}

/**
 * Options for eigsNear (find eigenvalues near sigma, symmetric).
 */
export interface EigsNearOptions {
  /**
   * Convergence tolerance.
   * @default 0 (machine precision)
   */
  tol?: number;

  /**
   * Number of Lanczos vectors to use.
   */
  ncv?: number;

  /**
   * Maximum number of iterations.
   * @default n * 10
   */
  maxiter?: number;

  /**
   * Starting vector for the iteration.
   */
  v0?: RealArray;

  /**
   * Whether to compute eigenvectors.
   * @default true
   */
  return_eigenvectors?: boolean;
}

/**
 * Options for eignNear (find eigenvalues near sigma, non-symmetric).
 */
export interface EignNearOptions {
  /**
   * Convergence tolerance.
   * @default 0 (machine precision)
   */
  tol?: number;

  /**
   * Number of Arnoldi vectors to use.
   */
  ncv?: number;

  /**
   * Maximum number of iterations.
   * @default n * 10
   */
  maxiter?: number;

  /**
   * Starting vector for the iteration.
   */
  v0?: RealArray;

  /**
   * Whether to compute eigenvectors.
   * @default true
   */
  return_eigenvectors?: boolean;
}

// ============================================================
// GRAPH ANALYSIS TYPES
// ============================================================

/**
 * Options for graph Laplacian eigenvalue computation.
 */
export interface LaplacianEigsOptions {
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
   * Use normalized Laplacian (L_sym = D^{-1/2} L D^{-1/2}).
   * @default false
   */
  normalized?: boolean;
}

/**
 * Result from graph Laplacian eigenvalue computation.
 */
export interface LaplacianEigsResult {
  /**
   * Computed eigenvalues (ascending order).
   */
  eigenvalues: Float64Array;

  /**
   * Computed eigenvectors.
   * eigenvectors[i] corresponds to eigenvalues[i].
   */
  eigenvectors?: Float64Array[];

  /**
   * Fiedler value (second smallest eigenvalue).
   * Measures algebraic connectivity of the graph.
   */
  fiedlerValue: number;

  /**
   * Fiedler vector (eigenvector for second smallest eigenvalue).
   * Used for spectral partitioning.
   */
  fiedlerVector?: Float64Array;

  /**
   * Number of iterations.
   */
  niter: number;

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
 * Options for PageRank computation.
 */
export interface PagerankOptions {
  /**
   * Damping factor (probability of following a link).
   * @default 0.85
   */
  alpha?: number;

  /**
   * Convergence tolerance.
   * @default 1e-8
   */
  tol?: number;

  /**
   * Maximum iterations.
   * @default 100
   */
  maxiter?: number;

  /**
   * Personalization vector (teleport distribution).
   * If not provided, uniform distribution is used.
   */
  personalization?: Float64Array;
}

/**
 * Result from PageRank computation.
 */
export interface PagerankResult {
  /**
   * PageRank scores for each node.
   */
  ranks: Float64Array;

  /**
   * Number of iterations performed.
   */
  niter: number;

  /**
   * Whether computation converged.
   */
  converged: boolean;

  /**
   * Final residual norm.
   */
  residual: number;
}

// ============================================================
// DIMENSIONALITY REDUCTION TYPES
// ============================================================

/**
 * Options for spectral embedding (Laplacian eigenmaps).
 */
export interface SpectralEmbeddingOptions {
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
   * Use normalized Laplacian.
   * @default true
   */
  normalized?: boolean;

  /**
   * Drop the first eigenvector (constant for connected graphs).
   * @default true
   */
  drop_first?: boolean;
}

/**
 * Result from spectral embedding.
 */
export interface SpectralEmbeddingResult {
  /**
   * Embedding coordinates (nComponents x n).
   * embedding[i] is the i-th dimension of the embedding.
   */
  embedding: Float64Array[];

  /**
   * Eigenvalues corresponding to each embedding dimension.
   */
  eigenvalues: Float64Array;

  /**
   * Number of iterations.
   */
  niter: number;

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
 * Options for truncated PCA.
 */
export interface TruncatedPCAOptions {
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
}

/**
 * Result from truncated PCA.
 */
export interface TruncatedPCAResult {
  /**
   * Principal components (nComponents x n).
   * components[i] is the i-th principal component direction.
   */
  components: Float64Array[];

  /**
   * Explained variance (eigenvalues) for each component.
   */
  explainedVariance: Float64Array;

  /**
   * Ratio of variance explained by each component.
   * Sums to less than 1 for truncated PCA.
   */
  explainedVarianceRatio: Float64Array;

  /**
   * Number of iterations.
   */
  niter: number;

  /**
   * Whether computation was successful.
   */
  success: boolean;

  /**
   * Status message.
   */
  message: string;
}

// ============================================================
// LINEAR ALGEBRA UTILITY TYPES
// ============================================================

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
   * Number of Arnoldi vectors.
   */
  ncv?: number;

  /**
   * Maximum iterations.
   */
  maxiter?: number;
}

/**
 * Result from spectral radius computation.
 */
export interface SpectralRadiusResult {
  /**
   * Spectral radius (largest magnitude eigenvalue).
   */
  spectralRadius: number;

  /**
   * The dominant eigenvalue (may be complex for non-symmetric).
   */
  dominantEigenvalue: { re: number; im: number };

  /**
   * Number of iterations.
   */
  niter: number;

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
 * Options for spectral norm computation.
 */
export interface SpectralNormOptions {
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
}

/**
 * Result from spectral norm computation.
 */
export interface SpectralNormResult {
  /**
   * Spectral norm (2-norm, largest singular value).
   */
  norm: number;

  /**
   * Number of iterations.
   */
  niter: number;

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
 * Options for condition number estimation.
 */
export interface CondestOptions {
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
}

/**
 * Result from condition number estimation.
 */
export interface CondestResult {
  /**
   * Estimated 2-norm condition number (sigma_max / sigma_min).
   */
  conditionNumber: number;

  /**
   * Largest singular value.
   */
  sigmaMax: number;

  /**
   * Smallest singular value.
   */
  sigmaMin: number;

  /**
   * Number of iterations.
   */
  niter: number;

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
 * Options for nuclear norm approximation.
 */
export interface NuclearNormOptions {
  /**
   * Number of singular values to compute.
   * Higher values give better approximation.
   * @default 10
   */
  k?: number;

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
}

/**
 * Result from nuclear norm approximation.
 */
export interface NuclearNormResult {
  /**
   * Approximate nuclear norm (sum of computed singular values).
   * This is a lower bound on the true nuclear norm.
   */
  nuclearNorm: number;

  /**
   * Computed singular values used in the approximation.
   */
  singularValues: Float64Array;

  /**
   * Number of singular values computed.
   */
  k: number;

  /**
   * Number of iterations.
   */
  niter: number;

  /**
   * Whether computation was successful.
   */
  success: boolean;

  /**
   * Status message.
   */
  message: string;
}
