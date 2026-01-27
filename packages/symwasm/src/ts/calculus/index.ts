/**
 * Calculus operations: differentiation, integration, limits, series.
 * @module calculus
 */

import type { Expr, Symbol } from '../core/index.js';
import { NotImplementedError } from '../errors.js';

/**
 * Differentiate an expression.
 * Mirrors sympy.diff.
 */
export function diff(
  _expr: Expr,
  _symbol: Symbol,
  _n?: number
): Expr {
  throw new NotImplementedError('symwasm.calculus.diff');
}

/**
 * Compute the indefinite or definite integral.
 * Mirrors sympy.integrate.
 */
export function integrate(
  _expr: Expr,
  _symbol: Symbol | [Symbol, number | Expr, number | Expr]
): Expr {
  throw new NotImplementedError('symwasm.calculus.integrate');
}

/**
 * Compute the limit of an expression.
 * Mirrors sympy.limit.
 */
export function limit(
  _expr: Expr,
  _symbol: Symbol,
  _point: number | Expr,
  _dir?: '+' | '-'
): Expr {
  throw new NotImplementedError('symwasm.calculus.limit');
}

/**
 * Compute a power series expansion.
 * Mirrors sympy.series.
 */
export function series(
  _expr: Expr,
  _symbol: Symbol,
  _point?: number | Expr,
  _n?: number
): Expr {
  throw new NotImplementedError('symwasm.calculus.series');
}

/**
 * Compute a symbolic summation.
 * Mirrors sympy.summation.
 */
export function summation(
  _f: Expr,
  _limits: [Symbol, number | Expr, number | Expr]
): Expr {
  throw new NotImplementedError('symwasm.calculus.summation');
}
