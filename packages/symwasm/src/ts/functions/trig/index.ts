/**
 * Trigonometric functions for symbolic expressions.
 *
 * This module provides the standard trigonometric functions (sine, cosine, tangent, etc.)
 * and their inverse functions (arcsine, arccosine, arctangent, etc.) for symbolic computation.
 * All functions work on symbolic expressions and return symbolic results that can be
 * simplified, differentiated, or evaluated numerically.
 *
 * @module functions/trig
 *
 * @example
 * ```typescript
 * import { Symbol, sin, cos, tan, pi, diff } from 'symwasm';
 *
 * const x = new Symbol('x');
 *
 * // Symbolic trigonometric expressions
 * sin(x);           // sin(x)
 * cos(pi);          // -1
 * tan(x);           // tan(x)
 *
 * // Derivatives work symbolically
 * diff(sin(x), x);  // cos(x)
 * diff(cos(x), x);  // -sin(x)
 * ```
 */

import { getWasmModule } from '../../wasm-loader.js';
import { createBasic, checkException } from '../../wasm-memory.js';
import type { Expr } from '../../core/expr.js';
import { exprFromWasm } from '../../core/expr-factory.js';

/**
 * Computes the sine of a symbolic expression.
 *
 * The sine function returns the ratio of the opposite side to the hypotenuse
 * in a right triangle. For symbolic expressions, it creates a symbolic sin(x)
 * that can be manipulated algebraically.
 *
 * @param x - The angle in radians (as a symbolic expression).
 * @returns The sine of x as a symbolic expression.
 *
 * @remarks
 * - sin(0) = 0
 * - sin(π/2) = 1
 * - sin(π) = 0
 * - sin(-x) = -sin(x) (odd function)
 * - sin(x + 2π) = sin(x) (periodic with period 2π)
 * - d/dx sin(x) = cos(x)
 *
 * @example
 * ```typescript
 * import { Symbol, Integer, sin, pi, div } from 'symwasm';
 *
 * const x = new Symbol('x');
 *
 * sin(x);                              // sin(x)
 * sin(div(pi, new Integer(2)));        // 1
 * sin(pi);                             // 0
 * ```
 *
 * @see {@link cos} - Cosine function
 * @see {@link asin} - Inverse sine
 * @see {@link sinh} - Hyperbolic sine
 */
export function sin(x: Expr): Expr {
  const wasm = getWasmModule();
  const obj = createBasic();
  try {
    const code = wasm._basic_sin(obj.getPtr(), x.getWasmPtr());
    checkException(code);
    return exprFromWasm(obj);
  } catch (e) {
    obj.free();
    throw e;
  }
}

/**
 * Computes the cosine of a symbolic expression.
 *
 * The cosine function returns the ratio of the adjacent side to the hypotenuse
 * in a right triangle. For symbolic expressions, it creates a symbolic cos(x)
 * that can be manipulated algebraically.
 *
 * @param x - The angle in radians (as a symbolic expression).
 * @returns The cosine of x as a symbolic expression.
 *
 * @remarks
 * - cos(0) = 1
 * - cos(π/2) = 0
 * - cos(π) = -1
 * - cos(-x) = cos(x) (even function)
 * - cos(x + 2π) = cos(x) (periodic with period 2π)
 * - d/dx cos(x) = -sin(x)
 * - sin²(x) + cos²(x) = 1
 *
 * @example
 * ```typescript
 * import { Symbol, Integer, cos, pi, div } from 'symwasm';
 *
 * const x = new Symbol('x');
 *
 * cos(x);                              // cos(x)
 * cos(new Integer(0));                 // 1
 * cos(pi);                             // -1
 * ```
 *
 * @see {@link sin} - Sine function
 * @see {@link acos} - Inverse cosine
 * @see {@link cosh} - Hyperbolic cosine
 */
export function cos(x: Expr): Expr {
  const wasm = getWasmModule();
  const obj = createBasic();
  try {
    const code = wasm._basic_cos(obj.getPtr(), x.getWasmPtr());
    checkException(code);
    return exprFromWasm(obj);
  } catch (e) {
    obj.free();
    throw e;
  }
}

/**
 * Computes the tangent of a symbolic expression.
 *
 * The tangent function is defined as tan(x) = sin(x)/cos(x), representing
 * the ratio of the opposite side to the adjacent side in a right triangle.
 *
 * @param x - The angle in radians (as a symbolic expression).
 * @returns The tangent of x as a symbolic expression.
 *
 * @remarks
 * - tan(0) = 0
 * - tan(π/4) = 1
 * - tan(-x) = -tan(x) (odd function)
 * - tan(x + π) = tan(x) (periodic with period π)
 * - tan(x) is undefined at x = π/2 + nπ
 * - d/dx tan(x) = sec²(x) = 1 + tan²(x)
 *
 * @example
 * ```typescript
 * import { Symbol, Integer, tan, pi, div } from 'symwasm';
 *
 * const x = new Symbol('x');
 *
 * tan(x);                              // tan(x)
 * tan(new Integer(0));                 // 0
 * tan(div(pi, new Integer(4)));        // 1
 * ```
 *
 * @see {@link sin} - Sine function
 * @see {@link cos} - Cosine function
 * @see {@link atan} - Inverse tangent
 * @see {@link cot} - Cotangent (reciprocal)
 */
export function tan(x: Expr): Expr {
  const wasm = getWasmModule();
  const obj = createBasic();
  try {
    const code = wasm._basic_tan(obj.getPtr(), x.getWasmPtr());
    checkException(code);
    return exprFromWasm(obj);
  } catch (e) {
    obj.free();
    throw e;
  }
}

/**
 * Computes the cotangent of a symbolic expression.
 *
 * The cotangent function is defined as cot(x) = cos(x)/sin(x) = 1/tan(x),
 * representing the ratio of the adjacent side to the opposite side.
 *
 * @param x - The angle in radians (as a symbolic expression).
 * @returns The cotangent of x as a symbolic expression.
 *
 * @remarks
 * - cot(π/4) = 1
 * - cot(π/2) = 0
 * - cot(-x) = -cot(x) (odd function)
 * - cot(x + π) = cot(x) (periodic with period π)
 * - cot(x) is undefined at x = nπ
 * - d/dx cot(x) = -csc²(x)
 *
 * @example
 * ```typescript
 * import { Symbol, cot, pi, div, Integer } from 'symwasm';
 *
 * const x = new Symbol('x');
 *
 * cot(x);                              // cot(x)
 * cot(div(pi, new Integer(4)));        // 1
 * cot(div(pi, new Integer(2)));        // 0
 * ```
 *
 * @see {@link tan} - Tangent (reciprocal)
 * @see {@link acot} - Inverse cotangent
 * @see {@link coth} - Hyperbolic cotangent
 */
export function cot(x: Expr): Expr {
  const wasm = getWasmModule();
  const obj = createBasic();
  try {
    const code = wasm._basic_cot(obj.getPtr(), x.getWasmPtr());
    checkException(code);
    return exprFromWasm(obj);
  } catch (e) {
    obj.free();
    throw e;
  }
}

/**
 * Computes the secant of a symbolic expression.
 *
 * The secant function is defined as sec(x) = 1/cos(x), representing
 * the ratio of the hypotenuse to the adjacent side.
 *
 * @param x - The angle in radians (as a symbolic expression).
 * @returns The secant of x as a symbolic expression.
 *
 * @remarks
 * - sec(0) = 1
 * - sec(-x) = sec(x) (even function)
 * - sec(x + 2π) = sec(x) (periodic with period 2π)
 * - sec(x) is undefined at x = π/2 + nπ
 * - d/dx sec(x) = sec(x)tan(x)
 * - 1 + tan²(x) = sec²(x)
 *
 * @example
 * ```typescript
 * import { Symbol, sec, Integer } from 'symwasm';
 *
 * const x = new Symbol('x');
 *
 * sec(x);                 // sec(x)
 * sec(new Integer(0));    // 1
 * ```
 *
 * @see {@link cos} - Cosine (reciprocal)
 * @see {@link asec} - Inverse secant
 * @see {@link sech} - Hyperbolic secant
 */
export function sec(x: Expr): Expr {
  const wasm = getWasmModule();
  const obj = createBasic();
  try {
    const code = wasm._basic_sec(obj.getPtr(), x.getWasmPtr());
    checkException(code);
    return exprFromWasm(obj);
  } catch (e) {
    obj.free();
    throw e;
  }
}

/**
 * Computes the cosecant of a symbolic expression.
 *
 * The cosecant function is defined as csc(x) = 1/sin(x), representing
 * the ratio of the hypotenuse to the opposite side.
 *
 * @param x - The angle in radians (as a symbolic expression).
 * @returns The cosecant of x as a symbolic expression.
 *
 * @remarks
 * - csc(π/2) = 1
 * - csc(-x) = -csc(x) (odd function)
 * - csc(x + 2π) = csc(x) (periodic with period 2π)
 * - csc(x) is undefined at x = nπ
 * - d/dx csc(x) = -csc(x)cot(x)
 * - 1 + cot²(x) = csc²(x)
 *
 * @example
 * ```typescript
 * import { Symbol, csc, pi, div, Integer } from 'symwasm';
 *
 * const x = new Symbol('x');
 *
 * csc(x);                              // csc(x)
 * csc(div(pi, new Integer(2)));        // 1
 * ```
 *
 * @see {@link sin} - Sine (reciprocal)
 * @see {@link acsc} - Inverse cosecant
 * @see {@link csch} - Hyperbolic cosecant
 */
export function csc(x: Expr): Expr {
  const wasm = getWasmModule();
  const obj = createBasic();
  try {
    const code = wasm._basic_csc(obj.getPtr(), x.getWasmPtr());
    checkException(code);
    return exprFromWasm(obj);
  } catch (e) {
    obj.free();
    throw e;
  }
}

/**
 * Computes the arcsine (inverse sine) of a symbolic expression.
 *
 * Returns the angle whose sine is x. The principal value is in [-π/2, π/2].
 *
 * @param x - The value (in the range [-1, 1] for real results).
 * @returns The arcsine of x as a symbolic expression.
 *
 * @remarks
 * - asin(0) = 0
 * - asin(1) = π/2
 * - asin(-1) = -π/2
 * - asin(-x) = -asin(x) (odd function)
 * - sin(asin(x)) = x
 * - d/dx asin(x) = 1/√(1 - x²)
 * - For |x| > 1, the result is complex
 *
 * @example
 * ```typescript
 * import { Symbol, asin, Integer } from 'symwasm';
 *
 * const x = new Symbol('x');
 *
 * asin(x);                 // asin(x)
 * asin(new Integer(0));    // 0
 * asin(new Integer(1));    // π/2
 * ```
 *
 * @see {@link sin} - Sine function
 * @see {@link acos} - Inverse cosine
 * @see {@link atan} - Inverse tangent
 */
export function asin(x: Expr): Expr {
  const wasm = getWasmModule();
  const obj = createBasic();
  try {
    const code = wasm._basic_asin(obj.getPtr(), x.getWasmPtr());
    checkException(code);
    return exprFromWasm(obj);
  } catch (e) {
    obj.free();
    throw e;
  }
}

/**
 * Computes the arccosine (inverse cosine) of a symbolic expression.
 *
 * Returns the angle whose cosine is x. The principal value is in [0, π].
 *
 * @param x - The value (in the range [-1, 1] for real results).
 * @returns The arccosine of x as a symbolic expression.
 *
 * @remarks
 * - acos(1) = 0
 * - acos(0) = π/2
 * - acos(-1) = π
 * - acos(x) + asin(x) = π/2
 * - cos(acos(x)) = x
 * - d/dx acos(x) = -1/√(1 - x²)
 * - For |x| > 1, the result is complex
 *
 * @example
 * ```typescript
 * import { Symbol, acos, Integer } from 'symwasm';
 *
 * const x = new Symbol('x');
 *
 * acos(x);                 // acos(x)
 * acos(new Integer(1));    // 0
 * acos(new Integer(0));    // π/2
 * acos(new Integer(-1));   // π
 * ```
 *
 * @see {@link cos} - Cosine function
 * @see {@link asin} - Inverse sine
 * @see {@link atan} - Inverse tangent
 */
export function acos(x: Expr): Expr {
  const wasm = getWasmModule();
  const obj = createBasic();
  try {
    const code = wasm._basic_acos(obj.getPtr(), x.getWasmPtr());
    checkException(code);
    return exprFromWasm(obj);
  } catch (e) {
    obj.free();
    throw e;
  }
}

/**
 * Computes the arctangent (inverse tangent) of a symbolic expression.
 *
 * Returns the angle whose tangent is x. The principal value is in (-π/2, π/2).
 *
 * @param x - The value (can be any real number).
 * @returns The arctangent of x as a symbolic expression.
 *
 * @remarks
 * - atan(0) = 0
 * - atan(1) = π/4
 * - atan(-x) = -atan(x) (odd function)
 * - atan(x) → π/2 as x → +∞
 * - atan(x) → -π/2 as x → -∞
 * - tan(atan(x)) = x
 * - d/dx atan(x) = 1/(1 + x²)
 *
 * @example
 * ```typescript
 * import { Symbol, atan, Integer } from 'symwasm';
 *
 * const x = new Symbol('x');
 *
 * atan(x);                 // atan(x)
 * atan(new Integer(0));    // 0
 * atan(new Integer(1));    // π/4
 * ```
 *
 * @see {@link tan} - Tangent function
 * @see {@link atan2} - Two-argument arctangent
 * @see {@link asin} - Inverse sine
 * @see {@link acos} - Inverse cosine
 */
export function atan(x: Expr): Expr {
  const wasm = getWasmModule();
  const obj = createBasic();
  try {
    const code = wasm._basic_atan(obj.getPtr(), x.getWasmPtr());
    checkException(code);
    return exprFromWasm(obj);
  } catch (e) {
    obj.free();
    throw e;
  }
}

/**
 * Computes the arc-cotangent (inverse cotangent) of a symbolic expression.
 *
 * Returns the angle whose cotangent is x. The principal value is in (0, π).
 *
 * @param x - The value (can be any real number except 0).
 * @returns The arc-cotangent of x as a symbolic expression.
 *
 * @remarks
 * - acot(1) = π/4
 * - acot(0) = π/2
 * - acot(-x) = π - acot(x)
 * - acot(x) → 0 as x → +∞
 * - acot(x) → π as x → -∞
 * - cot(acot(x)) = x
 * - d/dx acot(x) = -1/(1 + x²)
 *
 * @example
 * ```typescript
 * import { Symbol, acot, Integer } from 'symwasm';
 *
 * const x = new Symbol('x');
 *
 * acot(x);                 // acot(x)
 * acot(new Integer(1));    // π/4
 * acot(new Integer(0));    // π/2
 * ```
 *
 * @see {@link cot} - Cotangent function
 * @see {@link atan} - Inverse tangent
 */
export function acot(x: Expr): Expr {
  const wasm = getWasmModule();
  const obj = createBasic();
  try {
    const code = wasm._basic_acot(obj.getPtr(), x.getWasmPtr());
    checkException(code);
    return exprFromWasm(obj);
  } catch (e) {
    obj.free();
    throw e;
  }
}

/**
 * Computes the arc-secant (inverse secant) of a symbolic expression.
 *
 * Returns the angle whose secant is x. The principal value is in [0, π], excluding π/2.
 *
 * @param x - The value (|x| ≥ 1 for real results).
 * @returns The arc-secant of x as a symbolic expression.
 *
 * @remarks
 * - asec(1) = 0
 * - asec(-1) = π
 * - asec(2) = π/3
 * - sec(asec(x)) = x
 * - d/dx asec(x) = 1/(|x|√(x² - 1))
 * - For |x| < 1, the result is complex
 *
 * @example
 * ```typescript
 * import { Symbol, asec, Integer } from 'symwasm';
 *
 * const x = new Symbol('x');
 *
 * asec(x);                 // asec(x)
 * asec(new Integer(1));    // 0
 * asec(new Integer(2));    // π/3
 * ```
 *
 * @see {@link sec} - Secant function
 * @see {@link acos} - Inverse cosine
 */
export function asec(x: Expr): Expr {
  const wasm = getWasmModule();
  const obj = createBasic();
  try {
    const code = wasm._basic_asec(obj.getPtr(), x.getWasmPtr());
    checkException(code);
    return exprFromWasm(obj);
  } catch (e) {
    obj.free();
    throw e;
  }
}

/**
 * Computes the arc-cosecant (inverse cosecant) of a symbolic expression.
 *
 * Returns the angle whose cosecant is x. The principal value is in [-π/2, π/2], excluding 0.
 *
 * @param x - The value (|x| ≥ 1 for real results).
 * @returns The arc-cosecant of x as a symbolic expression.
 *
 * @remarks
 * - acsc(1) = π/2
 * - acsc(-1) = -π/2
 * - acsc(2) = π/6
 * - acsc(-x) = -acsc(x) (odd function)
 * - csc(acsc(x)) = x
 * - d/dx acsc(x) = -1/(|x|√(x² - 1))
 * - For |x| < 1, the result is complex
 *
 * @example
 * ```typescript
 * import { Symbol, acsc, Integer } from 'symwasm';
 *
 * const x = new Symbol('x');
 *
 * acsc(x);                 // acsc(x)
 * acsc(new Integer(1));    // π/2
 * acsc(new Integer(2));    // π/6
 * ```
 *
 * @see {@link csc} - Cosecant function
 * @see {@link asin} - Inverse sine
 */
export function acsc(x: Expr): Expr {
  const wasm = getWasmModule();
  const obj = createBasic();
  try {
    const code = wasm._basic_acsc(obj.getPtr(), x.getWasmPtr());
    checkException(code);
    return exprFromWasm(obj);
  } catch (e) {
    obj.free();
    throw e;
  }
}

/**
 * Computes the two-argument arctangent of symbolic expressions.
 *
 * Returns the angle θ in (-π, π] such that (cos(θ), sin(θ)) is proportional to (x, y).
 * This is the angle of the vector (x, y) from the positive x-axis, correctly handling
 * all four quadrants.
 *
 * @param y - The y-coordinate (opposite side / numerator in y/x).
 * @param x - The x-coordinate (adjacent side / denominator in y/x).
 * @returns The angle θ = atan2(y, x) as a symbolic expression.
 *
 * @remarks
 * - atan2(y, x) computes the angle of point (x, y) from origin
 * - atan2(0, 1) = 0 (positive x-axis)
 * - atan2(1, 0) = π/2 (positive y-axis)
 * - atan2(0, -1) = π (negative x-axis)
 * - atan2(-1, 0) = -π/2 (negative y-axis)
 * - atan2(y, x) = atan(y/x) when x > 0
 * - Unlike atan(y/x), atan2 returns the correct quadrant
 * - Used to compute the argument of complex numbers: arg(x + iy) = atan2(y, x)
 *
 * @example
 * ```typescript
 * import { Symbol, atan2, Integer } from 'symwasm';
 *
 * const y = new Symbol('y');
 * const x = new Symbol('x');
 *
 * // General symbolic case
 * atan2(y, x);
 *
 * // Point (1, 1) is at angle π/4
 * atan2(new Integer(1), new Integer(1));  // π/4
 *
 * // Point (0, 1) is at angle π/2
 * atan2(new Integer(1), new Integer(0));  // π/2
 *
 * // Point (-1, 0) is at angle π
 * atan2(new Integer(0), new Integer(-1)); // π
 * ```
 *
 * @see {@link atan} - Single-argument arctangent
 * @see {@link arg} - Complex argument (uses atan2)
 */
export function atan2(y: Expr, x: Expr): Expr {
  const wasm = getWasmModule();
  const obj = createBasic();
  try {
    const code = wasm._basic_atan2(obj.getPtr(), y.getWasmPtr(), x.getWasmPtr());
    checkException(code);
    return exprFromWasm(obj);
  } catch (e) {
    obj.free();
    throw e;
  }
}
