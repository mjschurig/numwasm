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
