/**
 * Calculus operations: differentiation, integration, limits, series.
 * @module calculus
 */

import type { Expr, Symbol } from '../core/index.js';
import { NotImplementedError } from '../errors.js';

// Re-export diff from core (implemented in Phase 2.2)
export { diff } from '../core/index.js';

// Re-export series from core (implemented in Phase 2.3)
export { series } from '../core/index.js';

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
 * Compute a symbolic summation.
 * Mirrors sympy.summation.
 */
export function summation(
  _f: Expr,
  _limits: [Symbol, number | Expr, number | Expr]
): Expr {
  throw new NotImplementedError('symwasm.calculus.summation');
}
