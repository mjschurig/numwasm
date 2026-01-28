/**
 * Type definitions for scipy.optimize.minimize
 */

/** Result of an optimization run — mirrors scipy.optimize.OptimizeResult */
export interface OptimizeResult {
  /** The solution array */
  x: number[];
  /** Whether the optimizer exited successfully */
  success: boolean;
  /** Termination status (solver-specific: 0=success, 1=maxfev, 2=maxiter) */
  status: number;
  /** Description of termination cause */
  message: string;
  /** Objective function value at x */
  fun: number;
  /** Number of function evaluations */
  nfev: number;
  /** Number of iterations */
  nit: number;
  /** Jacobian at x (if available) */
  jac?: number[];
  /** Inverse Hessian approximation (if available) */
  hess_inv?: number[][];
  /** Number of Jacobian evaluations */
  njev?: number;
  /** Number of Hessian evaluations */
  nhev?: number;
  /** Maximum constraint violation (constrained methods) */
  maxcv?: number;
  /** Final simplex (Nelder-Mead only) */
  final_simplex?: { vertices: number[][]; values: number[] };
}

/** Supported minimize methods */
export type MinimizeMethod =
  | 'Nelder-Mead'
  | 'Powell'
  | 'CG'
  | 'BFGS'
  | 'Newton-CG'
  | 'L-BFGS-B'
  | 'TNC'
  | 'SLSQP';

/** Objective function: takes array of coordinates, returns scalar */
export type ObjectiveFunction = (x: number[]) => number;

/** Jacobian function or finite-difference method */
export type JacobianOption =
  | ((x: number[]) => number[])
  | '2-point'
  | '3-point'
  | 'cs'
  | boolean;

/** Variable bounds */
export interface Bounds {
  lb: number[];
  ub: number[];
}

/** Options for minimize() */
export interface MinimizeOptions {
  method?: MinimizeMethod;
  jac?: JacobianOption;
  hess?: ((x: number[]) => number[][]) | '2-point' | '3-point' | 'cs';
  bounds?: Bounds | Array<[number | null, number | null]>;
  tol?: number;
  callback?: (intermediate: OptimizeResult) => boolean | void;
  options?: NelderMeadOptions | BFGSOptions | LBFGSBOptions | Record<string, unknown>;
}

/** Nelder-Mead specific options */
export interface NelderMeadOptions {
  maxiter?: number;
  maxfev?: number;
  xatol?: number;
  fatol?: number;
  adaptive?: boolean;
  initial_simplex?: number[][];
}

/** BFGS specific options */
export interface BFGSOptions {
  maxiter?: number;
  gtol?: number;
  norm?: number;
  eps?: number;
  xrtol?: number;
  c1?: number;
  c2?: number;
}

/** L-BFGS-B specific options */
export interface LBFGSBOptions {
  maxcor?: number;
  ftol?: number;
  gtol?: number;
  eps?: number;
  maxfun?: number;
  maxiter?: number;
  maxls?: number;
}
