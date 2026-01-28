/**
 * Type definitions for the integrate module.
 * Mirrors scipy.integrate types.
 */

/**
 * A scalar integrand function: f(x, ...args) → number.
 */
export type IntegrandFunction = (x: number, ...args: number[]) => number;

/**
 * Options for the quad() function.
 */
export interface QuadOptions {
  /** Extra arguments passed to the integrand function. */
  args?: number | number[];
  /** If true, return (result, error, infodict). */
  fullOutput?: boolean;
  /** Absolute error tolerance. Default: 1.49e-8. */
  epsabs?: number;
  /** Relative error tolerance. Default: 1.49e-8. */
  epsrel?: number;
  /** Upper bound on the number of subintervals. Default: 50. */
  limit?: number;
  /** If true, integrand returns complex values. Default: false. */
  complexFunc?: boolean;
}

/**
 * Information dictionary returned by quad() with fullOutput=true.
 */
export interface QuadInfoDict {
  /** Number of function evaluations. */
  neval: number;
  /** Number of subintervals produced. */
  last: number;
  /** Left endpoints of subintervals. */
  alist: Float64Array;
  /** Right endpoints of subintervals. */
  blist: Float64Array;
  /** Integral approximations on subintervals. */
  rlist: Float64Array;
  /** Absolute error estimates on subintervals. */
  elist: Float64Array;
  /** Error ordering indices. */
  iord: Int32Array;
}

/**
 * Result of quad() without fullOutput.
 * [result, abserr]
 */
export type QuadResult = [number, number];

/**
 * Result of quad() with fullOutput=true.
 * [result, abserr, infodict]
 */
export type QuadFullResult = [number, number, QuadInfoDict];

/**
 * Result of quad() with fullOutput=true when a warning was raised.
 * [result, abserr, infodict, message]
 */
export type QuadFullResultWithMessage = [number, number, QuadInfoDict, string];
