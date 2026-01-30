/**
 * Type definitions for sparse linear algebra module
 */

/**
 * Result returned by iterative solvers (cg, gmres, bicgstab)
 */
export interface IterativeSolverResult {
  /** Solution vector */
  x: Float64Array;
  /** Exit info: 0=converged, >0=maxiter reached, <0=breakdown/error */
  info: number;
  /** Number of iterations performed */
  iterations: number;
  /** Final residual norm */
  residualNorm: number;
}

/**
 * Options for iterative solvers
 */
export interface IterativeSolverOptions {
  /** Initial guess for solution (default: zeros) */
  x0?: Float64Array;
  /** Relative tolerance (default: 1e-5) */
  tol?: number;
  /** Absolute tolerance (default: 0) */
  atol?: number;
  /** Maximum iterations (default: n, the problem size) */
  maxiter?: number;
  /** Preconditioner M such that M approximates A^-1 */
  M?: LinearOperatorLike;
  /** Callback function called each iteration */
  callback?: (xk: Float64Array) => void;
}

/**
 * Additional options for GMRES
 */
export interface GMRESOptions extends Omit<IterativeSolverOptions, 'callback'> {
  /** Restart parameter (default: min(20, n)) */
  restart?: number;
  /** Callback receives residual norm instead of xk */
  callback?: (residualNorm: number) => void;
}

/**
 * Types that can be used as linear operators
 */
export type LinearOperatorLike = LinearOperator | import('../base.js').SparseMatrix | number[][];

/**
 * Abstract linear operator interface
 */
export interface LinearOperator {
  /** Shape [m, n] where operator maps R^n -> R^m */
  readonly shape: [number, number];
  /** Data type (always 'float64' for now) */
  readonly dtype: string;
  /** Matrix-vector product: y = A @ x */
  matvec(x: Float64Array): Float64Array;
  /** Transpose matrix-vector product: y = A.T @ x */
  rmatvec(x: Float64Array): Float64Array;
}

// ============================================================
// Eigenvalue Solver Types
// ============================================================

/**
 * Options for eigsh (symmetric eigenvalue solver)
 */
export interface EigshOptions {
  /** Number of eigenvalues to compute (default: 6) */
  k?: number;
  /** Which eigenvalues: 'LA' (largest algebraic), 'SA' (smallest algebraic),
   *  'LM' (largest magnitude), 'SM' (smallest magnitude), 'BE' (both ends) */
  which?: 'LA' | 'SA' | 'LM' | 'SM' | 'BE';
  /** Shift for shift-invert mode */
  sigma?: number;
  /** Relative tolerance (default: 0 = machine precision) */
  tol?: number;
  /** Maximum iterations (default: n*10) */
  maxiter?: number;
  /** Number of Lanczos vectors (default: min(n, max(2*k+1, 20))) */
  ncv?: number;
  /** Starting vector */
  v0?: Float64Array;
  /** Return eigenvectors (default: true) */
  return_eigenvectors?: boolean;
}

/**
 * Options for eigs (general eigenvalue solver)
 */
export interface EigsOptions {
  /** Number of eigenvalues to compute (default: 6) */
  k?: number;
  /** Which eigenvalues: 'LM', 'SM', 'LR', 'SR', 'LI', 'SI' */
  which?: 'LM' | 'SM' | 'LR' | 'SR' | 'LI' | 'SI';
  /** Shift for shift-invert mode */
  sigma?: number | { real: number; imag: number };
  /** Relative tolerance (default: 0 = machine precision) */
  tol?: number;
  /** Maximum iterations (default: n*10) */
  maxiter?: number;
  /** Number of Arnoldi vectors (default: min(n, max(2*k+1, 20))) */
  ncv?: number;
  /** Starting vector */
  v0?: Float64Array;
  /** Return eigenvectors (default: true) */
  return_eigenvectors?: boolean;
}

/**
 * Result from eigenvalue solvers
 */
export interface EigsResult {
  /** Eigenvalues (real part for non-symmetric) */
  values: Float64Array;
  /** Imaginary parts of eigenvalues (only for non-symmetric, may be undefined) */
  valuesImag?: Float64Array;
  /** Eigenvectors as columns (n x nconv), flattened column-major */
  vectors?: Float64Array;
  /** Number of converged eigenvalues */
  nconv: number;
  /** Number of iterations performed */
  iterations: number;
}

/**
 * Options for svds (truncated SVD)
 */
export interface SvdsOptions {
  /** Number of singular values to compute (default: 6) */
  k?: number;
  /** Relative tolerance */
  tol?: number;
  /** Maximum iterations */
  maxiter?: number;
  /** Number of Lanczos vectors */
  ncv?: number;
  /** Starting vector */
  v0?: Float64Array;
  /** Return singular vectors: true (both), 'u', 'vh', or false */
  return_singular_vectors?: boolean | 'u' | 'vh';
}

/**
 * Result from truncated SVD
 */
export interface SvdsResult {
  /** Left singular vectors (m x k), column-major (optional) */
  u?: Float64Array;
  /** Singular values (k) */
  s: Float64Array;
  /** Right singular vectors transposed (k x n), row-major (optional) */
  vt?: Float64Array;
}

// ============================================================
// Matrix Exponential Types
// ============================================================

/**
 * Options for matrix exponential
 */
export interface ExpmOptions {
  /** Maximum Pade approximation order (default: 13) */
  padeOrder?: number;
}

/**
 * Options for expm_multiply
 */
export interface ExpmMultiplyOptions {
  /** Krylov subspace dimension (default: min(30, n)) */
  krylovDim?: number;
  /** Tolerance for happy breakdown detection (default: 1e-12) */
  tol?: number;
  /** Time parameter t for exp(t*A)*v (default: 1.0) */
  t?: number;
}

/**
 * Result from expm_multiply
 */
export interface ExpmMultiplyResult {
  /** The result vector exp(t*A) * v */
  result: Float64Array;
  /** Number of Krylov iterations used */
  iterations: number;
  /** Estimated error bound */
  errorEstimate: number;
}
