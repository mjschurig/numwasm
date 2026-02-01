/**
 * Hyperbolic functions for symbolic expressions.
 *
 * This module provides hyperbolic functions (sinh, cosh, tanh, etc.) and their
 * inverse functions (asinh, acosh, atanh, etc.) for symbolic computation.
 * Hyperbolic functions are analogous to trigonometric functions but are based
 * on hyperbolas rather than circles.
 *
 * @module functions/hyperbolic
 *
 * @example
 * ```typescript
 * import { Symbol, sinh, cosh, tanh, exp, diff } from 'symwasm';
 *
 * const x = new Symbol('x');
 *
 * // Hyperbolic functions in terms of exponentials
 * // sinh(x) = (e^x - e^(-x)) / 2
 * // cosh(x) = (e^x + e^(-x)) / 2
 *
 * // Derivatives
 * diff(sinh(x), x);  // cosh(x)
 * diff(cosh(x), x);  // sinh(x)
 *
 * // Identity: cosh²(x) - sinh²(x) = 1
 * ```
 */

import { getWasmModule } from '../../wasm-loader.js';
import { createBasic, checkException } from '../../wasm-memory.js';
import type { Expr } from '../../core/expr.js';
import { exprFromWasm } from '../../core/expr-factory.js';

/**
 * Computes the hyperbolic sine of a symbolic expression.
 *
 * The hyperbolic sine is defined as sinh(x) = (e^x - e^(-x)) / 2.
 * It represents the y-coordinate of a point on the unit hyperbola.
 *
 * @param x - The argument (as a symbolic expression).
 * @returns The hyperbolic sine of x as a symbolic expression.
 *
 * @remarks
 * - sinh(0) = 0
 * - sinh(-x) = -sinh(x) (odd function)
 * - sinh(x) → ∞ as x → ∞
 * - sinh(x) → -∞ as x → -∞
 * - d/dx sinh(x) = cosh(x)
 * - sinh(x) = -i * sin(ix) (relation to trig)
 * - cosh²(x) - sinh²(x) = 1
 *
 * @example
 * ```typescript
 * import { Symbol, sinh, Integer, exp, sub, div } from 'symwasm';
 *
 * const x = new Symbol('x');
 *
 * sinh(x);                 // sinh(x)
 * sinh(new Integer(0));    // 0
 *
 * // sinh(x) = (e^x - e^(-x)) / 2
 * // These are equivalent:
 * // sinh(x)
 * // div(sub(exp(x), exp(neg(x))), new Integer(2))
 * ```
 *
 * @see {@link cosh} - Hyperbolic cosine
 * @see {@link asinh} - Inverse hyperbolic sine
 * @see {@link sin} - Trigonometric sine
 */
export function sinh(x: Expr): Expr {
  const wasm = getWasmModule();
  const obj = createBasic();
  try {
    const code = wasm._basic_sinh(obj.getPtr(), x.getWasmPtr());
    checkException(code);
    return exprFromWasm(obj);
  } catch (e) {
    obj.free();
    throw e;
  }
}

/**
 * Computes the hyperbolic cosine of a symbolic expression.
 *
 * The hyperbolic cosine is defined as cosh(x) = (e^x + e^(-x)) / 2.
 * It represents the x-coordinate of a point on the unit hyperbola,
 * and also describes the shape of a hanging chain (catenary).
 *
 * @param x - The argument (as a symbolic expression).
 * @returns The hyperbolic cosine of x as a symbolic expression.
 *
 * @remarks
 * - cosh(0) = 1
 * - cosh(-x) = cosh(x) (even function)
 * - cosh(x) ≥ 1 for all real x
 * - cosh(x) → ∞ as x → ±∞
 * - d/dx cosh(x) = sinh(x)
 * - cosh(x) = cos(ix) (relation to trig)
 * - cosh²(x) - sinh²(x) = 1
 *
 * @example
 * ```typescript
 * import { Symbol, cosh, Integer } from 'symwasm';
 *
 * const x = new Symbol('x');
 *
 * cosh(x);                 // cosh(x)
 * cosh(new Integer(0));    // 1
 * ```
 *
 * @see {@link sinh} - Hyperbolic sine
 * @see {@link acosh} - Inverse hyperbolic cosine
 * @see {@link cos} - Trigonometric cosine
 */
export function cosh(x: Expr): Expr {
  const wasm = getWasmModule();
  const obj = createBasic();
  try {
    const code = wasm._basic_cosh(obj.getPtr(), x.getWasmPtr());
    checkException(code);
    return exprFromWasm(obj);
  } catch (e) {
    obj.free();
    throw e;
  }
}

/**
 * Computes the hyperbolic tangent of a symbolic expression.
 *
 * The hyperbolic tangent is defined as tanh(x) = sinh(x)/cosh(x) = (e^x - e^(-x))/(e^x + e^(-x)).
 * It is a sigmoid-like function bounded between -1 and 1.
 *
 * @param x - The argument (as a symbolic expression).
 * @returns The hyperbolic tangent of x as a symbolic expression.
 *
 * @remarks
 * - tanh(0) = 0
 * - tanh(-x) = -tanh(x) (odd function)
 * - -1 < tanh(x) < 1 for all finite x
 * - tanh(x) → 1 as x → +∞
 * - tanh(x) → -1 as x → -∞
 * - d/dx tanh(x) = sech²(x) = 1 - tanh²(x)
 * - tanh(x) = -i * tan(ix) (relation to trig)
 *
 * @example
 * ```typescript
 * import { Symbol, tanh, Integer } from 'symwasm';
 *
 * const x = new Symbol('x');
 *
 * tanh(x);                 // tanh(x)
 * tanh(new Integer(0));    // 0
 * ```
 *
 * @see {@link sinh} - Hyperbolic sine
 * @see {@link cosh} - Hyperbolic cosine
 * @see {@link atanh} - Inverse hyperbolic tangent
 * @see {@link tan} - Trigonometric tangent
 */
export function tanh(x: Expr): Expr {
  const wasm = getWasmModule();
  const obj = createBasic();
  try {
    const code = wasm._basic_tanh(obj.getPtr(), x.getWasmPtr());
    checkException(code);
    return exprFromWasm(obj);
  } catch (e) {
    obj.free();
    throw e;
  }
}

/**
 * Computes the hyperbolic cotangent of a symbolic expression.
 *
 * The hyperbolic cotangent is defined as coth(x) = cosh(x)/sinh(x) = 1/tanh(x).
 *
 * @param x - The argument (as a symbolic expression, x ≠ 0).
 * @returns The hyperbolic cotangent of x as a symbolic expression.
 *
 * @remarks
 * - coth(-x) = -coth(x) (odd function)
 * - |coth(x)| > 1 for all x ≠ 0
 * - coth(x) → 1 as x → +∞
 * - coth(x) → -1 as x → -∞
 * - coth(x) has a singularity at x = 0
 * - d/dx coth(x) = -csch²(x)
 *
 * @example
 * ```typescript
 * import { Symbol, coth, Integer } from 'symwasm';
 *
 * const x = new Symbol('x');
 *
 * coth(x);                 // coth(x)
 * coth(new Integer(1));    // coth(1) ≈ 1.313
 * ```
 *
 * @see {@link tanh} - Hyperbolic tangent (reciprocal)
 * @see {@link acoth} - Inverse hyperbolic cotangent
 * @see {@link cot} - Trigonometric cotangent
 */
export function coth(x: Expr): Expr {
  const wasm = getWasmModule();
  const obj = createBasic();
  try {
    const code = wasm._basic_coth(obj.getPtr(), x.getWasmPtr());
    checkException(code);
    return exprFromWasm(obj);
  } catch (e) {
    obj.free();
    throw e;
  }
}

/**
 * Computes the hyperbolic secant of a symbolic expression.
 *
 * The hyperbolic secant is defined as sech(x) = 1/cosh(x) = 2/(e^x + e^(-x)).
 *
 * @param x - The argument (as a symbolic expression).
 * @returns The hyperbolic secant of x as a symbolic expression.
 *
 * @remarks
 * - sech(0) = 1
 * - sech(-x) = sech(x) (even function)
 * - 0 < sech(x) ≤ 1 for all real x
 * - sech(x) → 0 as x → ±∞
 * - d/dx sech(x) = -sech(x)tanh(x)
 * - sech(x) = sec(ix) (relation to trig)
 *
 * @example
 * ```typescript
 * import { Symbol, sech, Integer } from 'symwasm';
 *
 * const x = new Symbol('x');
 *
 * sech(x);                 // sech(x)
 * sech(new Integer(0));    // 1
 * ```
 *
 * @see {@link cosh} - Hyperbolic cosine (reciprocal)
 * @see {@link asech} - Inverse hyperbolic secant
 * @see {@link sec} - Trigonometric secant
 */
export function sech(x: Expr): Expr {
  const wasm = getWasmModule();
  const obj = createBasic();
  try {
    const code = wasm._basic_sech(obj.getPtr(), x.getWasmPtr());
    checkException(code);
    return exprFromWasm(obj);
  } catch (e) {
    obj.free();
    throw e;
  }
}

/**
 * Computes the hyperbolic cosecant of a symbolic expression.
 *
 * The hyperbolic cosecant is defined as csch(x) = 1/sinh(x) = 2/(e^x - e^(-x)).
 *
 * @param x - The argument (as a symbolic expression, x ≠ 0).
 * @returns The hyperbolic cosecant of x as a symbolic expression.
 *
 * @remarks
 * - csch(-x) = -csch(x) (odd function)
 * - csch(x) has a singularity at x = 0
 * - csch(x) → 0 as x → ±∞
 * - d/dx csch(x) = -csch(x)coth(x)
 * - csch(x) = -i * csc(ix) (relation to trig)
 *
 * @example
 * ```typescript
 * import { Symbol, csch, Integer } from 'symwasm';
 *
 * const x = new Symbol('x');
 *
 * csch(x);                 // csch(x)
 * csch(new Integer(1));    // csch(1) ≈ 0.851
 * ```
 *
 * @see {@link sinh} - Hyperbolic sine (reciprocal)
 * @see {@link acsch} - Inverse hyperbolic cosecant
 * @see {@link csc} - Trigonometric cosecant
 */
export function csch(x: Expr): Expr {
  const wasm = getWasmModule();
  const obj = createBasic();
  try {
    const code = wasm._basic_csch(obj.getPtr(), x.getWasmPtr());
    checkException(code);
    return exprFromWasm(obj);
  } catch (e) {
    obj.free();
    throw e;
  }
}

/**
 * Computes the inverse hyperbolic sine of a symbolic expression.
 *
 * Returns the value y such that sinh(y) = x.
 * Also known as the area hyperbolic sine: asinh(x) = ln(x + √(x² + 1)).
 *
 * @param x - The argument (can be any real number).
 * @returns The inverse hyperbolic sine of x as a symbolic expression.
 *
 * @remarks
 * - asinh(0) = 0
 * - asinh(-x) = -asinh(x) (odd function)
 * - asinh(x) is defined for all real x
 * - sinh(asinh(x)) = x
 * - d/dx asinh(x) = 1/√(x² + 1)
 * - asinh(x) = ln(x + √(x² + 1))
 *
 * @example
 * ```typescript
 * import { Symbol, asinh, Integer } from 'symwasm';
 *
 * const x = new Symbol('x');
 *
 * asinh(x);                // asinh(x)
 * asinh(new Integer(0));   // 0
 * asinh(new Integer(1));   // ln(1 + √2)
 * ```
 *
 * @see {@link sinh} - Hyperbolic sine
 * @see {@link acosh} - Inverse hyperbolic cosine
 * @see {@link asin} - Inverse trigonometric sine
 */
export function asinh(x: Expr): Expr {
  const wasm = getWasmModule();
  const obj = createBasic();
  try {
    const code = wasm._basic_asinh(obj.getPtr(), x.getWasmPtr());
    checkException(code);
    return exprFromWasm(obj);
  } catch (e) {
    obj.free();
    throw e;
  }
}

/**
 * Computes the inverse hyperbolic cosine of a symbolic expression.
 *
 * Returns the non-negative value y such that cosh(y) = x.
 * Also known as the area hyperbolic cosine: acosh(x) = ln(x + √(x² - 1)).
 *
 * @param x - The argument (x ≥ 1 for real results).
 * @returns The inverse hyperbolic cosine of x as a symbolic expression.
 *
 * @remarks
 * - acosh(1) = 0
 * - acosh(x) ≥ 0 for x ≥ 1
 * - acosh(x) is only real for x ≥ 1
 * - cosh(acosh(x)) = x
 * - d/dx acosh(x) = 1/√(x² - 1)
 * - acosh(x) = ln(x + √(x² - 1))
 *
 * @example
 * ```typescript
 * import { Symbol, acosh, Integer } from 'symwasm';
 *
 * const x = new Symbol('x');
 *
 * acosh(x);                // acosh(x)
 * acosh(new Integer(1));   // 0
 * acosh(new Integer(2));   // ln(2 + √3)
 * ```
 *
 * @see {@link cosh} - Hyperbolic cosine
 * @see {@link asinh} - Inverse hyperbolic sine
 * @see {@link acos} - Inverse trigonometric cosine
 */
export function acosh(x: Expr): Expr {
  const wasm = getWasmModule();
  const obj = createBasic();
  try {
    const code = wasm._basic_acosh(obj.getPtr(), x.getWasmPtr());
    checkException(code);
    return exprFromWasm(obj);
  } catch (e) {
    obj.free();
    throw e;
  }
}

/**
 * Computes the inverse hyperbolic tangent of a symbolic expression.
 *
 * Returns the value y such that tanh(y) = x.
 * Also known as the area hyperbolic tangent: atanh(x) = (1/2) * ln((1+x)/(1-x)).
 *
 * @param x - The argument (-1 < x < 1 for real results).
 * @returns The inverse hyperbolic tangent of x as a symbolic expression.
 *
 * @remarks
 * - atanh(0) = 0
 * - atanh(-x) = -atanh(x) (odd function)
 * - atanh(x) is only real for |x| < 1
 * - atanh(x) → +∞ as x → 1⁻
 * - atanh(x) → -∞ as x → -1⁺
 * - tanh(atanh(x)) = x
 * - d/dx atanh(x) = 1/(1 - x²)
 *
 * @example
 * ```typescript
 * import { Symbol, atanh, Integer, Rational } from 'symwasm';
 *
 * const x = new Symbol('x');
 *
 * atanh(x);                      // atanh(x)
 * atanh(new Integer(0));         // 0
 * atanh(new Rational(1, 2));     // atanh(1/2)
 * ```
 *
 * @see {@link tanh} - Hyperbolic tangent
 * @see {@link asinh} - Inverse hyperbolic sine
 * @see {@link atan} - Inverse trigonometric tangent
 */
export function atanh(x: Expr): Expr {
  const wasm = getWasmModule();
  const obj = createBasic();
  try {
    const code = wasm._basic_atanh(obj.getPtr(), x.getWasmPtr());
    checkException(code);
    return exprFromWasm(obj);
  } catch (e) {
    obj.free();
    throw e;
  }
}

/**
 * Computes the inverse hyperbolic cotangent of a symbolic expression.
 *
 * Returns the value y such that coth(y) = x.
 * Also known as the area hyperbolic cotangent: acoth(x) = (1/2) * ln((x+1)/(x-1)).
 *
 * @param x - The argument (|x| > 1 for real results).
 * @returns The inverse hyperbolic cotangent of x as a symbolic expression.
 *
 * @remarks
 * - acoth(-x) = -acoth(x) (odd function)
 * - acoth(x) is only real for |x| > 1
 * - acoth(x) → 0 as x → ±∞
 * - coth(acoth(x)) = x
 * - d/dx acoth(x) = 1/(1 - x²) for |x| > 1
 *
 * @example
 * ```typescript
 * import { Symbol, acoth, Integer } from 'symwasm';
 *
 * const x = new Symbol('x');
 *
 * acoth(x);                // acoth(x)
 * acoth(new Integer(2));   // acoth(2)
 * ```
 *
 * @see {@link coth} - Hyperbolic cotangent
 * @see {@link atanh} - Inverse hyperbolic tangent
 * @see {@link acot} - Inverse trigonometric cotangent
 */
export function acoth(x: Expr): Expr {
  const wasm = getWasmModule();
  const obj = createBasic();
  try {
    const code = wasm._basic_acoth(obj.getPtr(), x.getWasmPtr());
    checkException(code);
    return exprFromWasm(obj);
  } catch (e) {
    obj.free();
    throw e;
  }
}

/**
 * Computes the inverse hyperbolic secant of a symbolic expression.
 *
 * Returns the non-negative value y such that sech(y) = x.
 * Also known as the area hyperbolic secant: asech(x) = ln((1 + √(1-x²))/x).
 *
 * @param x - The argument (0 < x ≤ 1 for real results).
 * @returns The inverse hyperbolic secant of x as a symbolic expression.
 *
 * @remarks
 * - asech(1) = 0
 * - asech(x) ≥ 0 for 0 < x ≤ 1
 * - asech(x) is only real for 0 < x ≤ 1
 * - asech(x) → +∞ as x → 0⁺
 * - sech(asech(x)) = x
 * - d/dx asech(x) = -1/(x√(1 - x²))
 *
 * @example
 * ```typescript
 * import { Symbol, asech, Integer, Rational } from 'symwasm';
 *
 * const x = new Symbol('x');
 *
 * asech(x);                      // asech(x)
 * asech(new Integer(1));         // 0
 * asech(new Rational(1, 2));     // asech(1/2)
 * ```
 *
 * @see {@link sech} - Hyperbolic secant
 * @see {@link acosh} - Inverse hyperbolic cosine
 * @see {@link asec} - Inverse trigonometric secant
 */
export function asech(x: Expr): Expr {
  const wasm = getWasmModule();
  const obj = createBasic();
  try {
    const code = wasm._basic_asech(obj.getPtr(), x.getWasmPtr());
    checkException(code);
    return exprFromWasm(obj);
  } catch (e) {
    obj.free();
    throw e;
  }
}

/**
 * Computes the inverse hyperbolic cosecant of a symbolic expression.
 *
 * Returns the value y such that csch(y) = x.
 * Also known as the area hyperbolic cosecant: acsch(x) = ln((1 + √(1+x²))/x) for x > 0.
 *
 * @param x - The argument (x ≠ 0 for real results).
 * @returns The inverse hyperbolic cosecant of x as a symbolic expression.
 *
 * @remarks
 * - acsch(-x) = -acsch(x) (odd function)
 * - acsch(x) is defined for all x ≠ 0
 * - acsch(x) → 0 as x → ±∞
 * - csch(acsch(x)) = x
 * - d/dx acsch(x) = -1/(|x|√(1 + x²))
 *
 * @example
 * ```typescript
 * import { Symbol, acsch, Integer } from 'symwasm';
 *
 * const x = new Symbol('x');
 *
 * acsch(x);                // acsch(x)
 * acsch(new Integer(1));   // acsch(1) = ln(1 + √2)
 * acsch(new Integer(2));   // acsch(2)
 * ```
 *
 * @see {@link csch} - Hyperbolic cosecant
 * @see {@link asinh} - Inverse hyperbolic sine
 * @see {@link acsc} - Inverse trigonometric cosecant
 */
export function acsch(x: Expr): Expr {
  const wasm = getWasmModule();
  const obj = createBasic();
  try {
    const code = wasm._basic_acsch(obj.getPtr(), x.getWasmPtr());
    checkException(code);
    return exprFromWasm(obj);
  } catch (e) {
    obj.free();
    throw e;
  }
}
