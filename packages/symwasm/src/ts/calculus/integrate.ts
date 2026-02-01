/**
 * Integration function for symbolic expressions.
 *
 * This module provides symbolic integration (antiderivative) capabilities.
 *
 * @module calculus/integrate
 */

import type { Expr } from '../core/expr.js';
import { NotImplementedError } from '../errors.js';

/**
 * Computes the indefinite or definite integral of a symbolic expression.
 *
 * For indefinite integrals, returns the antiderivative (without the constant of integration).
 * For definite integrals, evaluates the integral over the specified bounds.
 *
 * @param _expr - The expression to integrate.
 * @param _symbol - Either a symbol for indefinite integration, or a tuple
 *   `[symbol, lower, upper]` for definite integration.
 * @returns The integral result.
 *
 * @throws NotImplementedError - This function is not yet implemented.
 *
 * @remarks
 * **Note: This function is a stub and not yet implemented.**
 *
 * When implemented, will support:
 * - Indefinite integrals: `integrate(expr, x)` returns ∫f(x)dx
 * - Definite integrals: `integrate(expr, [x, a, b])` returns ∫ₐᵇf(x)dx
 *
 * Standard integration rules will apply:
 * - ∫x^n dx = x^(n+1)/(n+1) for n ≠ -1 (power rule)
 * - ∫1/x dx = log(x)
 * - ∫e^x dx = e^x
 * - ∫sin(x) dx = -cos(x)
 * - ∫cos(x) dx = sin(x)
 *
 * The Fundamental Theorem of Calculus relates differentiation and integration:
 * - d/dx [∫f(x)dx] = f(x)
 * - ∫ₐᵇf'(x)dx = f(b) - f(a)
 *
 * @example
 * ```typescript
 * import { Symbol, Integer, integrate, pow, sin } from 'symwasm';
 *
 * const x = new Symbol('x');
 *
 * // Indefinite integral (when implemented)
 * // integrate(pow(x, new Integer(2)), x);  // x³/3
 * // integrate(sin(x), x);                  // -cos(x)
 *
 * // Definite integral (when implemented)
 * // integrate(pow(x, new Integer(2)), [x, new Integer(0), new Integer(1)]);  // 1/3
 * ```
 *
 * @see {@link diff} - Symbolic differentiation (inverse operation)
 * @see {@link summation} - Symbolic summation
 */
export function integrate(
  _expr: Expr,
  _symbol: Expr | [Expr, number | Expr, number | Expr]
): Expr {
  throw new NotImplementedError('symwasm.calculus.integrate');
}
