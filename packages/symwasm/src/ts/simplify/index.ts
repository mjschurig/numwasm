/**
 * Expression simplification functions.
 * @module simplify
 */

import type { Expr } from '../core/index.js';
import { NotImplementedError } from '../errors.js';

/**
 * Simplify an expression using heuristics.
 * Mirrors sympy.simplify.
 */
export function simplify(_expr: Expr): Expr {
  throw new NotImplementedError('symwasm.simplify.simplify');
}

/**
 * Expand an expression (distribute multiplication over addition).
 * Mirrors sympy.expand.
 */
export function expand(_expr: Expr): Expr {
  throw new NotImplementedError('symwasm.simplify.expand');
}

/**
 * Factor a polynomial expression.
 * Mirrors sympy.factor.
 */
export function factor(_expr: Expr): Expr {
  throw new NotImplementedError('symwasm.simplify.factor');
}

/**
 * Collect common powers of a term in an expression.
 * Mirrors sympy.collect.
 */
export function collect(_expr: Expr, _syms: Expr | Expr[]): Expr {
  throw new NotImplementedError('symwasm.simplify.collect');
}

/**
 * Cancel common factors in a rational function.
 * Mirrors sympy.cancel.
 */
export function cancel(_expr: Expr): Expr {
  throw new NotImplementedError('symwasm.simplify.cancel');
}

/**
 * Simplify trigonometric expressions.
 * Mirrors sympy.trigsimp.
 */
export function trigsimp(_expr: Expr): Expr {
  throw new NotImplementedError('symwasm.simplify.trigsimp');
}

/**
 * Simplify expressions with radicals.
 * Mirrors sympy.radsimp.
 */
export function radsimp(_expr: Expr): Expr {
  throw new NotImplementedError('symwasm.simplify.radsimp');
}

/**
 * Simplify combinatorial expressions with powers.
 * Mirrors sympy.powsimp.
 */
export function powsimp(_expr: Expr): Expr {
  throw new NotImplementedError('symwasm.simplify.powsimp');
}
