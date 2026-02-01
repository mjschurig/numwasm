/**
 * Complex number functions for symbolic expressions.
 *
 * This module provides functions for working with complex numbers in symbolic form,
 * including extracting real and imaginary parts, computing conjugates, and calculating
 * the argument (phase angle) of complex expressions.
 *
 * @module functions/complex
 *
 * @example
 * ```typescript
 * import { Symbol, I, re, im, conjugate, arg } from 'symwasm';
 *
 * const x = new Symbol('x');
 * const y = new Symbol('y');
 * const z = add(x, mul(I, y));  // z = x + i*y
 *
 * re(z);         // Returns: x
 * im(z);         // Returns: y
 * conjugate(z);  // Returns: x - i*y
 * arg(z);        // Returns: atan2(y, x)
 * ```
 */

import { getWasmModule } from '../../wasm-loader.js';
import { createBasic, checkException } from '../../wasm-memory.js';
import type { Expr } from '../../core/expr.js';
import { exprFromWasm } from '../../core/expr-factory.js';
import { atan2 } from '../trig/index.js';

/**
 * Computes the complex conjugate of an expression.
 *
 * For a complex number z = a + bi, the conjugate is z* = a - bi.
 * This function works symbolically, so it can handle expressions
 * containing symbolic variables.
 *
 * @param x - The expression to conjugate.
 * @returns The complex conjugate of the expression.
 *
 * @remarks
 * - For real expressions, conjugate(x) = x
 * - conjugate(conjugate(x)) = x
 * - conjugate(x + y) = conjugate(x) + conjugate(y)
 * - conjugate(x * y) = conjugate(x) * conjugate(y)
 * - conjugate(I) = -I where I is the imaginary unit
 *
 * @example
 * ```typescript
 * import { Symbol, I, add, mul, conjugate } from 'symwasm';
 *
 * const x = new Symbol('x');
 * const y = new Symbol('y');
 *
 * // Conjugate of a pure imaginary: conjugate(3i) = -3i
 * conjugate(mul(new Integer(3), I));
 *
 * // Conjugate of complex expression: conjugate(x + yi) = x - yi
 * const z = add(x, mul(I, y));
 * conjugate(z);  // Returns: x - i*y
 * ```
 *
 * @see {@link re} - Extract real part
 * @see {@link im} - Extract imaginary part
 */
export function conjugate(x: Expr): Expr {
  const wasm = getWasmModule();
  const obj = createBasic();
  try {
    const code = wasm._basic_conjugate(obj.getPtr(), x.getWasmPtr());
    checkException(code);
    return exprFromWasm(obj);
  } catch (e) {
    obj.free();
    throw e;
  }
}

/**
 * Extracts the real part of a complex expression.
 *
 * For a complex number z = a + bi, re(z) returns a.
 * This function works symbolically on expressions containing
 * symbolic variables and the imaginary unit I.
 *
 * @param x - The expression to extract the real part from.
 * @returns The real part of the expression.
 *
 * @remarks
 * - re(a + bi) = a
 * - re(re(x)) = re(x)
 * - re(x + y) = re(x) + re(y)
 * - re(I) = 0 where I is the imaginary unit
 * - For purely real x, re(x) = x
 *
 * @example
 * ```typescript
 * import { Symbol, I, Integer, add, mul, re } from 'symwasm';
 *
 * // Real part of 3 + 4i
 * const z = add(new Integer(3), mul(new Integer(4), I));
 * re(z);  // Returns: 3
 *
 * // Real part of symbolic expression
 * const x = new Symbol('x');
 * const y = new Symbol('y');
 * re(add(x, mul(I, y)));  // Returns: x
 * ```
 *
 * @see {@link im} - Extract imaginary part
 * @see {@link conjugate} - Complex conjugate
 * @see {@link arg} - Argument (phase angle)
 */
export function re(x: Expr): Expr {
  const wasm = getWasmModule();
  const obj = createBasic();
  try {
    const code = wasm._complex_base_real_part(obj.getPtr(), x.getWasmPtr());
    checkException(code);
    return exprFromWasm(obj);
  } catch (e) {
    obj.free();
    throw e;
  }
}

/**
 * Extracts the imaginary part of a complex expression.
 *
 * For a complex number z = a + bi, im(z) returns b (not bi).
 * This function works symbolically on expressions containing
 * symbolic variables and the imaginary unit I.
 *
 * @param x - The expression to extract the imaginary part from.
 * @returns The imaginary part of the expression (coefficient of i).
 *
 * @remarks
 * - im(a + bi) = b (not bi)
 * - im(im(x)) = 0 (since im(x) is real)
 * - im(x + y) = im(x) + im(y)
 * - im(I) = 1 where I is the imaginary unit
 * - For purely real x, im(x) = 0
 *
 * @example
 * ```typescript
 * import { Symbol, I, Integer, add, mul, im } from 'symwasm';
 *
 * // Imaginary part of 3 + 4i
 * const z = add(new Integer(3), mul(new Integer(4), I));
 * im(z);  // Returns: 4 (not 4i)
 *
 * // Imaginary part of symbolic expression
 * const x = new Symbol('x');
 * const y = new Symbol('y');
 * im(add(x, mul(I, y)));  // Returns: y
 * ```
 *
 * @see {@link re} - Extract real part
 * @see {@link conjugate} - Complex conjugate
 * @see {@link arg} - Argument (phase angle)
 */
export function im(x: Expr): Expr {
  const wasm = getWasmModule();
  const obj = createBasic();
  try {
    const code = wasm._complex_base_imaginary_part(obj.getPtr(), x.getWasmPtr());
    checkException(code);
    return exprFromWasm(obj);
  } catch (e) {
    obj.free();
    throw e;
  }
}

/**
 * Computes the argument (phase angle) of a complex expression.
 *
 * For a complex number z = a + bi, arg(z) = atan2(b, a), which gives
 * the angle θ in the polar form z = r * e^(iθ). The result is in radians,
 * ranging from -π to π.
 *
 * @param x - The complex expression to compute the argument of.
 * @returns The argument (phase angle) in radians.
 *
 * @remarks
 * - arg(x) is computed as atan2(im(x), re(x))
 * - arg(r * e^(iθ)) = θ for r > 0
 * - arg(positive real) = 0
 * - arg(negative real) = π
 * - arg(positive imaginary) = π/2
 * - arg(negative imaginary) = -π/2
 * - arg(x * y) = arg(x) + arg(y) (mod 2π)
 *
 * @example
 * ```typescript
 * import { Symbol, I, Integer, add, mul, arg, pi } from 'symwasm';
 *
 * // Argument of 1 + i (should be π/4)
 * const z = add(new Integer(1), I);
 * arg(z);  // Returns: atan2(1, 1) = π/4
 *
 * // Argument of i (should be π/2)
 * arg(I);  // Returns: π/2
 *
 * // Argument of -1 (should be π)
 * arg(new Integer(-1));  // Returns: π
 * ```
 *
 * @see {@link re} - Extract real part
 * @see {@link im} - Extract imaginary part
 * @see {@link atan2} - Two-argument arctangent
 */
export function arg(x: Expr): Expr {
  return atan2(im(x), re(x));
}
