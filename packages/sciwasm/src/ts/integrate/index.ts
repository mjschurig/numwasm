/**
 * Numerical integration and ODE solvers.
 * @module integrate
 */

import { NotImplementedError } from '../errors.js';

/** Result of numerical quadrature. */
export interface QuadResult {
  /** Estimated value of the integral. */
  value: number;
  /** Estimated absolute error. */
  error: number;
}

/**
 * Compute a definite integral using adaptive quadrature.
 * Mirrors scipy.integrate.quad.
 */
export function quad(
  _func: (x: number) => number,
  _a: number,
  _b: number,
  _options?: { epsabs?: number; epsrel?: number; limit?: number }
): QuadResult {
  throw new NotImplementedError('sciwasm.integrate.quad');
}

/**
 * Compute a double integral.
 * Mirrors scipy.integrate.dblquad.
 */
export function dblquad(
  _func: (y: number, x: number) => number,
  _a: number,
  _b: number,
  _gfun: (x: number) => number,
  _hfun: (x: number) => number
): QuadResult {
  throw new NotImplementedError('sciwasm.integrate.dblquad');
}

/**
 * Compute a triple integral.
 * Mirrors scipy.integrate.tplquad.
 */
export function tplquad(
  _func: (z: number, y: number, x: number) => number,
  _a: number,
  _b: number,
  _gfun: (x: number) => number,
  _hfun: (x: number) => number,
  _qfun: (x: number, y: number) => number,
  _rfun: (x: number, y: number) => number
): QuadResult {
  throw new NotImplementedError('sciwasm.integrate.tplquad');
}

/**
 * Integrate using the composite trapezoidal rule.
 * Mirrors scipy.integrate.trapezoid.
 */
export function trapezoid(
  _y: number[],
  _x?: number[],
  _dx?: number
): number {
  throw new NotImplementedError('sciwasm.integrate.trapezoid');
}

/**
 * Integrate using the composite Simpson's rule.
 * Mirrors scipy.integrate.simpson.
 */
export function simpson(
  _y: number[],
  _x?: number[],
  _dx?: number
): number {
  throw new NotImplementedError('sciwasm.integrate.simpson');
}

/**
 * Integrate a system of ODEs.
 * Mirrors scipy.integrate.odeint.
 */
export function odeint(
  _func: (y: number[], t: number) => number[],
  _y0: number[],
  _t: number[],
  _options?: { rtol?: number; atol?: number }
): number[][] {
  throw new NotImplementedError('sciwasm.integrate.odeint');
}
