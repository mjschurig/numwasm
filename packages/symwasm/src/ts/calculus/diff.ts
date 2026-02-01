/**
 * Differentiation function for symbolic expressions.
 *
 * This module provides symbolic differentiation capabilities, allowing
 * computation of derivatives with respect to one or more variables.
 *
 * @module calculus/diff
 */

import { getWasmModule } from '../wasm-loader.js';
import { createBasic, checkException } from '../wasm-memory.js';
import type { Expr } from '../core/expr.js';
import { exprFromWasm } from '../core/expr-factory.js';

/**
 * Computes the derivative of a symbolic expression.
 *
 * Performs symbolic differentiation with support for multiple variables
 * and higher-order derivatives. The function handles various calling
 * conventions to express different types of derivatives.
 *
 * @param expr - The expression to differentiate.
 * @param args - Symbols and optional derivative orders in alternating pattern.
 * @returns The derivative expression.
 *
 * @throws Error if no symbol argument is provided.
 * @throws Error if a number appears without a preceding symbol.
 * @throws Error if derivative order is negative.
 *
 * @remarks
 * Supported calling conventions:
 * - `diff(expr, x)` — First derivative with respect to x: df/dx
 * - `diff(expr, x, 2)` — Second derivative: d²f/dx²
 * - `diff(expr, x, y)` — Mixed partial: ∂²f/∂x∂y
 * - `diff(expr, x, y, z)` — Mixed partial: ∂³f/∂x∂y∂z
 * - `diff(expr, x, 2, y, 3)` — Higher mixed partial: ∂⁵f/∂x²∂y³
 *
 * Standard differentiation rules apply:
 * - d/dx [c] = 0 (constant)
 * - d/dx [x] = 1
 * - d/dx [x^n] = n*x^(n-1) (power rule)
 * - d/dx [f + g] = df/dx + dg/dx (sum rule)
 * - d/dx [f * g] = f*dg/dx + g*df/dx (product rule)
 * - d/dx [f(g(x))] = f'(g(x)) * g'(x) (chain rule)
 *
 * @example
 * ```typescript
 * import { Symbol, diff, sin, cos, pow, mul, add } from 'symwasm';
 *
 * const x = new Symbol('x');
 * const y = new Symbol('y');
 *
 * // First derivative
 * diff(pow(x, 2), x);           // 2*x
 * diff(sin(x), x);              // cos(x)
 *
 * // Second derivative
 * diff(sin(x), x, 2);           // -sin(x)
 *
 * // Mixed partial derivatives
 * const f = mul(x, y);          // f = x*y
 * diff(f, x, y);                // 1 (∂²f/∂x∂y)
 *
 * // Higher order mixed partials
 * const g = mul(pow(x, 2), pow(y, 3));  // g = x²y³
 * diff(g, x, 2, y, 2);                   // 12*y (∂⁴g/∂x²∂y²)
 * ```
 *
 * @see {@link series} - Taylor series expansion
 * @see {@link integrate} - Integration (antiderivative)
 */
export function diff(expr: Expr, ...args: (Expr | number)[]): Expr {
  if (args.length === 0) {
    throw new Error('diff() requires at least one symbol argument');
  }

  const wasm = getWasmModule();
  let current = expr;

  // Parse args: alternating [symbol, count?] pattern
  let i = 0;
  while (i < args.length) {
    const arg = args[i];

    // Must be an Expr (symbol)
    if (typeof arg === 'number') {
      throw new Error(
        `diff() argument at position ${i + 1} must be a symbol, got number without preceding symbol`
      );
    }

    const symbol = arg as Expr;
    i++;

    // Check if next argument is a number (derivative order)
    let order = 1;
    if (i < args.length && typeof args[i] === 'number') {
      order = args[i] as number;
      i++;
    }

    if (order < 0) {
      throw new Error('Derivative order must be non-negative');
    }

    // Apply differentiation 'order' times with respect to 'symbol'
    for (let j = 0; j < order; j++) {
      const result = createBasic();
      try {
        const code = wasm._basic_diff(
          result.getPtr(),
          current.getWasmPtr(),
          symbol.getWasmPtr()
        );
        checkException(code);
        current = exprFromWasm(result);
      } catch (e) {
        result.free();
        throw e;
      }
    }
  }

  return current;
}
