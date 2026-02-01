/**
 * Taylor/Maclaurin series expansion function for symbolic expressions.
 *
 * This module provides symbolic series expansion capabilities, allowing
 * computation of Taylor series approximations of functions around a point.
 *
 * @module calculus/series
 */

import { getWasmModule } from '../wasm-loader.js';
import { createBasic, checkException } from '../wasm-memory.js';
import type { Expr } from '../core/expr.js';
import { exprFromWasm } from '../core/expr-factory.js';

/**
 * Computes the Taylor series expansion of an expression around a point.
 *
 * Expands a symbolic expression as a polynomial approximation using its
 * Taylor series. Currently supports expansion around x = 0 (Maclaurin series).
 *
 * @param expr - The expression to expand.
 * @param x - The expansion variable (symbol).
 * @param x0 - The expansion point (currently only 0 is supported).
 * @param n - Number of terms in the series (precision parameter, default 6).
 * @returns The series expansion as a polynomial (without the O() remainder term).
 *
 * @throws Error if x0 is not 0 (currently only Maclaurin series supported).
 *
 * @remarks
 * The Taylor series of f(x) around x = 0 is:
 * f(x) = f(0) + f'(0)*x + f''(0)*x²/2! + f'''(0)*x³/3! + ...
 *
 * Common series expansions:
 * - sin(x) = x - x³/6 + x⁵/120 - ...
 * - cos(x) = 1 - x²/2 + x⁴/24 - ...
 * - exp(x) = 1 + x + x²/2 + x³/6 + ...
 * - log(1+x) = x - x²/2 + x³/3 - ...
 * - 1/(1-x) = 1 + x + x² + x³ + ...
 *
 * The `n` parameter controls how many terms are computed. For precision `n`,
 * terms up to x^(n-1) are included.
 *
 * @example
 * ```typescript
 * import { Symbol, series, sin, cos, exp, log, add } from 'symwasm';
 *
 * const x = new Symbol('x');
 *
 * // Maclaurin series of sin(x) with 6 terms
 * series(sin(x), x);
 * // Returns: x - x³/6 + x⁵/120
 *
 * // Maclaurin series of exp(x) with 4 terms
 * series(exp(x), x, 0, 4);
 * // Returns: 1 + x + x²/2 + x³/6
 *
 * // Maclaurin series of cos(x)
 * series(cos(x), x, 0, 5);
 * // Returns: 1 - x²/2 + x⁴/24
 *
 * // Maclaurin series of 1/(1-x) = geometric series
 * series(div(new Integer(1), sub(new Integer(1), x)), x, 0, 5);
 * // Returns: 1 + x + x² + x³ + x⁴
 * ```
 *
 * @see {@link diff} - Symbolic differentiation
 * @see {@link limit} - Limit computation
 */
export function series(
  expr: Expr,
  x: Expr,
  x0: number | Expr = 0,
  n: number = 6
): Expr {
  // Validate x0 = 0 (only supported case)
  const x0Val = typeof x0 === 'number' ? x0 : null;
  if (x0Val !== 0) {
    throw new Error('series() currently only supports expansion around x=0');
  }

  const wasm = getWasmModule();
  const result = createBasic();

  try {
    const code = wasm._basic_series(
      result.getPtr(),
      expr.getWasmPtr(),
      x.getWasmPtr(),
      n
    );
    checkException(code);
    return exprFromWasm(result);
  } catch (e) {
    result.free();
    throw e;
  }
}
