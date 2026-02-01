/**
 * Minimum and maximum functions for symbolic expressions.
 *
 * This module provides functions for computing the minimum and maximum
 * of multiple symbolic expressions. These functions work symbolically
 * and can handle both numeric and symbolic arguments.
 *
 * @module functions/minmax
 *
 * @example
 * ```typescript
 * import { Symbol, Integer, Max, Min } from 'symwasm';
 *
 * const x = new Symbol('x');
 * const y = new Symbol('y');
 *
 * // Numeric maximum and minimum
 * Max(new Integer(1), new Integer(5), new Integer(3));  // 5
 * Min(new Integer(1), new Integer(5), new Integer(3));  // 1
 *
 * // Symbolic expressions
 * Max(x, y);          // Max(x, y)
 * Min(x, new Integer(0));  // Min(x, 0)
 * ```
 */

import { getWasmModule } from '../../wasm-loader.js';
import { SymEngineVec, createBasic, checkException } from '../../wasm-memory.js';
import type { Expr } from '../../core/expr.js';
import { exprFromWasm } from '../../core/expr-factory.js';

/**
 * Computes the maximum of multiple symbolic expressions.
 *
 * Returns the largest value among the given arguments. For symbolic arguments
 * that cannot be compared numerically, returns a symbolic Max expression that
 * may simplify when more information is available.
 *
 * @param args - Two or more expressions to find the maximum of.
 * @returns The maximum value as a symbolic expression.
 *
 * @remarks
 * - Max(a, b) = a if a ≥ b, else b
 * - Max(a, a) = a
 * - Max(a, b) = Max(b, a) (commutative)
 * - Max(a, Max(b, c)) = Max(a, b, c) (associative)
 * - Max(a, -∞) = a
 * - Max(a, ∞) = ∞
 * - d/da Max(a, b) = 1 if a > b, 0 if a < b (undefined at a = b)
 *
 * @example
 * ```typescript
 * import { Symbol, Integer, Max, abs } from 'symwasm';
 *
 * const x = new Symbol('x');
 * const y = new Symbol('y');
 *
 * // Numeric comparison
 * Max(new Integer(3), new Integer(7));     // 7
 * Max(new Integer(-2), new Integer(5));    // 5
 *
 * // Multiple arguments
 * Max(new Integer(1), new Integer(5), new Integer(3), new Integer(2));  // 5
 *
 * // Symbolic expressions - stays as Max until values are known
 * Max(x, y);                               // Max(x, y)
 *
 * // Can simplify with known relationships
 * Max(abs(x), new Integer(0));             // abs(x) (since abs ≥ 0)
 * ```
 *
 * @see {@link Min} - Minimum function
 * @see {@link abs} - Absolute value
 */
export function Max(...args: Expr[]): Expr {
  const wasm = getWasmModule();
  const vec = new SymEngineVec();
  const result = createBasic();

  try {
    for (const a of args) {
      wasm._vecbasic_push_back(vec.getPtr(), a.getWasmPtr());
    }
    const code = wasm._basic_max(result.getPtr(), vec.getPtr());
    checkException(code);
    return exprFromWasm(result);
  } catch (e) {
    result.free();
    throw e;
  } finally {
    vec.free();
  }
}

/**
 * Computes the minimum of multiple symbolic expressions.
 *
 * Returns the smallest value among the given arguments. For symbolic arguments
 * that cannot be compared numerically, returns a symbolic Min expression that
 * may simplify when more information is available.
 *
 * @param args - Two or more expressions to find the minimum of.
 * @returns The minimum value as a symbolic expression.
 *
 * @remarks
 * - Min(a, b) = a if a ≤ b, else b
 * - Min(a, a) = a
 * - Min(a, b) = Min(b, a) (commutative)
 * - Min(a, Min(b, c)) = Min(a, b, c) (associative)
 * - Min(a, -∞) = -∞
 * - Min(a, ∞) = a
 * - Min(a, b) = -Max(-a, -b)
 * - d/da Min(a, b) = 1 if a < b, 0 if a > b (undefined at a = b)
 *
 * @example
 * ```typescript
 * import { Symbol, Integer, Min, neg } from 'symwasm';
 *
 * const x = new Symbol('x');
 * const y = new Symbol('y');
 *
 * // Numeric comparison
 * Min(new Integer(3), new Integer(7));     // 3
 * Min(new Integer(-2), new Integer(5));    // -2
 *
 * // Multiple arguments
 * Min(new Integer(1), new Integer(5), new Integer(3), new Integer(2));  // 1
 *
 * // Symbolic expressions - stays as Min until values are known
 * Min(x, y);                               // Min(x, y)
 *
 * // Can simplify with known relationships
 * Min(x, new Integer(0));                  // Min(x, 0)
 * ```
 *
 * @see {@link Max} - Maximum function
 * @see {@link abs} - Absolute value
 */
export function Min(...args: Expr[]): Expr {
  const wasm = getWasmModule();
  const vec = new SymEngineVec();
  const result = createBasic();

  try {
    for (const a of args) {
      wasm._vecbasic_push_back(vec.getPtr(), a.getWasmPtr());
    }
    const code = wasm._basic_min(result.getPtr(), vec.getPtr());
    checkException(code);
    return exprFromWasm(result);
  } catch (e) {
    result.free();
    throw e;
  } finally {
    vec.free();
  }
}
