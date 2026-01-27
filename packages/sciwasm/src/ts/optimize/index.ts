/**
 * Optimization and root-finding routines.
 * @module optimize
 */

import { NotImplementedError } from '../errors.js';

/** Result of an optimization routine. */
export interface OptimizeResult {
  /** Solution vector. */
  x: number[];
  /** Whether the optimizer converged. */
  success: boolean;
  /** Description of the termination. */
  message: string;
  /** Value of the objective function at x. */
  fun: number;
  /** Number of function evaluations. */
  nfev: number;
  /** Number of iterations. */
  nit: number;
}

/** Result of a root-finding routine. */
export interface RootResult {
  root: number;
  converged: boolean;
  iterations: number;
}

/**
 * Minimize a scalar function of one or more variables.
 * Mirrors scipy.optimize.minimize.
 */
export function minimize(
  _fun: (x: number[]) => number,
  _x0: number[],
  _options?: { method?: string; tol?: number; maxiter?: number }
): OptimizeResult {
  throw new NotImplementedError('sciwasm.optimize.minimize');
}

/**
 * Solve a nonlinear least-squares problem.
 * Mirrors scipy.optimize.least_squares.
 */
export function least_squares(
  _fun: (x: number[]) => number[],
  _x0: number[]
): OptimizeResult {
  throw new NotImplementedError('sciwasm.optimize.least_squares');
}

/**
 * Find a root of a scalar function.
 * Mirrors scipy.optimize.root_scalar.
 */
export function root_scalar(
  _f: (x: number) => number,
  _options?: { bracket?: [number, number]; method?: string; xtol?: number }
): RootResult {
  throw new NotImplementedError('sciwasm.optimize.root_scalar');
}

/**
 * Linear programming: minimize a linear objective subject to constraints.
 * Mirrors scipy.optimize.linprog.
 */
export function linprog(
  _c: number[],
  _options?: { A_ub?: number[][]; b_ub?: number[]; A_eq?: number[][]; b_eq?: number[] }
): OptimizeResult {
  throw new NotImplementedError('sciwasm.optimize.linprog');
}

/**
 * Fit a function to data using nonlinear least squares.
 * Mirrors scipy.optimize.curve_fit.
 */
export function curve_fit(
  _f: (x: number, ...params: number[]) => number,
  _xdata: number[],
  _ydata: number[],
  _p0?: number[]
): { popt: number[]; pcov: number[][] } {
  throw new NotImplementedError('sciwasm.optimize.curve_fit');
}
