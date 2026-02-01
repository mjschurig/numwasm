/**
 * Elementary and special mathematical functions for symbolic expressions.
 *
 * This module provides exponential, logarithmic, root, and other fundamental
 * mathematical functions that operate on symbolic expressions. All functions
 * return symbolic results that can be further manipulated algebraically.
 *
 * @module functions
 *
 * @example
 * ```typescript
 * import { Symbol, exp, log, sqrt, abs, floor, ceiling } from 'symwasm';
 *
 * const x = new Symbol('x');
 *
 * // Exponential and logarithmic
 * exp(x);           // e^x
 * log(x);           // ln(x)
 *
 * // Roots
 * sqrt(x);          // √x
 * cbrt(x);          // ∛x
 *
 * // Other functions
 * abs(x);           // |x|
 * floor(x);         // ⌊x⌋
 * ceiling(x);       // ⌈x⌉
 * ```
 */

import { getWasmModule } from '../wasm-loader.js';
import { createBasic, checkException } from '../wasm-memory.js';
import type { Expr } from '../core/expr.js';
import { exprFromWasm } from '../core/expr-factory.js';

// Re-export trig functions
export * from './trig/index.js';

// Re-export hyperbolic functions
export * from './hyperbolic/index.js';

// Re-export special functions
export * from './special/index.js';

// Re-export complex functions
export * from './complex/index.js';

// Re-export min/max functions
export * from './minmax/index.js';

// ============================================================================
// Exponential & Logarithmic Functions
// ============================================================================

/**
 * Computes the exponential function of a symbolic expression.
 *
 * Returns e raised to the power of x, where e is Euler's number
 * (approximately 2.71828). For symbolic expressions, creates a
 * symbolic exp(x) that can be manipulated algebraically.
 *
 * @param x - The exponent (as a symbolic expression).
 * @returns e^x as a symbolic expression.
 *
 * @remarks
 * - exp(0) = 1
 * - exp(1) = e
 * - exp(-x) = 1/exp(x)
 * - exp(a + b) = exp(a) * exp(b)
 * - exp(log(x)) = x
 * - d/dx exp(x) = exp(x)
 * - exp(iπ) = -1 (Euler's identity)
 *
 * @example
 * ```typescript
 * import { Symbol, Integer, exp, log, I, pi, mul } from 'symwasm';
 *
 * const x = new Symbol('x');
 *
 * exp(x);                       // e^x
 * exp(new Integer(0));          // 1
 * exp(new Integer(1));          // e
 * exp(log(x));                  // x (simplifies)
 *
 * // Euler's identity: e^(i*π) = -1
 * exp(mul(I, pi));              // -1
 * ```
 *
 * @see {@link log} - Natural logarithm (inverse function)
 * @see {@link E} - Euler's number constant
 */
export function exp(x: Expr): Expr {
  const wasm = getWasmModule();
  const obj = createBasic();
  try {
    const code = wasm._basic_exp(obj.getPtr(), x.getWasmPtr());
    checkException(code);
    return exprFromWasm(obj);
  } catch (e) {
    obj.free();
    throw e;
  }
}

/**
 * Computes the natural logarithm of a symbolic expression.
 *
 * Returns the natural (base e) logarithm of x. For symbolic expressions,
 * creates a symbolic log(x) that can be manipulated algebraically.
 *
 * @param x - The argument (positive real or complex expression).
 * @returns ln(x) as a symbolic expression.
 *
 * @remarks
 * - log(1) = 0
 * - log(e) = 1
 * - log(x*y) = log(x) + log(y)
 * - log(x/y) = log(x) - log(y)
 * - log(x^n) = n*log(x)
 * - log(exp(x)) = x
 * - d/dx log(x) = 1/x
 * - log(-1) = iπ (principal value)
 *
 * @example
 * ```typescript
 * import { Symbol, Integer, log, exp, E } from 'symwasm';
 *
 * const x = new Symbol('x');
 *
 * log(x);                       // ln(x)
 * log(new Integer(1));          // 0
 * log(E);                       // 1
 * log(exp(x));                  // x (simplifies)
 *
 * // For other bases, use: log_b(x) = log(x) / log(b)
 * ```
 *
 * @see {@link exp} - Exponential function (inverse function)
 * @see {@link E} - Euler's number constant
 */
export function log(x: Expr): Expr {
  const wasm = getWasmModule();
  const obj = createBasic();
  try {
    const code = wasm._basic_log(obj.getPtr(), x.getWasmPtr());
    checkException(code);
    return exprFromWasm(obj);
  } catch (e) {
    obj.free();
    throw e;
  }
}

/**
 * Computes the square root of a symbolic expression.
 *
 * Returns the principal square root of x. Equivalent to pow(x, 1/2) or x^(1/2).
 * For negative real numbers, returns a complex result involving I.
 *
 * @param x - The radicand (as a symbolic expression).
 * @returns √x as a symbolic expression.
 *
 * @remarks
 * - sqrt(0) = 0
 * - sqrt(1) = 1
 * - sqrt(4) = 2
 * - sqrt(x)^2 = x
 * - sqrt(x*y) = sqrt(x)*sqrt(y) for non-negative x, y
 * - sqrt(-1) = I (imaginary unit)
 * - d/dx sqrt(x) = 1/(2*sqrt(x))
 *
 * @example
 * ```typescript
 * import { Symbol, Integer, sqrt, pow, Rational } from 'symwasm';
 *
 * const x = new Symbol('x');
 *
 * sqrt(x);                        // √x
 * sqrt(new Integer(4));           // 2
 * sqrt(new Integer(2));           // √2 (irrational, left symbolic)
 *
 * // sqrt(x) is equivalent to x^(1/2)
 * // sqrt(x) equals pow(x, new Rational(1, 2))
 * ```
 *
 * @see {@link cbrt} - Cube root
 * @see {@link pow} - General power function
 */
export function sqrt(x: Expr): Expr {
  const wasm = getWasmModule();
  const obj = createBasic();
  try {
    const code = wasm._basic_sqrt(obj.getPtr(), x.getWasmPtr());
    checkException(code);
    return exprFromWasm(obj);
  } catch (e) {
    obj.free();
    throw e;
  }
}

/**
 * Computes the cube root of a symbolic expression.
 *
 * Returns the real cube root of x. Unlike square root, the cube root
 * of a negative number is a negative real number. Equivalent to pow(x, 1/3).
 *
 * @param x - The radicand (as a symbolic expression).
 * @returns ∛x as a symbolic expression.
 *
 * @remarks
 * - cbrt(0) = 0
 * - cbrt(1) = 1
 * - cbrt(8) = 2
 * - cbrt(-8) = -2
 * - cbrt(x)^3 = x
 * - cbrt(x*y) = cbrt(x)*cbrt(y)
 * - d/dx cbrt(x) = 1/(3*cbrt(x)^2)
 *
 * @example
 * ```typescript
 * import { Symbol, Integer, cbrt } from 'symwasm';
 *
 * const x = new Symbol('x');
 *
 * cbrt(x);                        // ∛x
 * cbrt(new Integer(8));           // 2
 * cbrt(new Integer(-8));          // -2
 * cbrt(new Integer(27));          // 3
 * ```
 *
 * @see {@link sqrt} - Square root
 * @see {@link pow} - General power function
 */
export function cbrt(x: Expr): Expr {
  const wasm = getWasmModule();
  const obj = createBasic();
  try {
    const code = wasm._basic_cbrt(obj.getPtr(), x.getWasmPtr());
    checkException(code);
    return exprFromWasm(obj);
  } catch (e) {
    obj.free();
    throw e;
  }
}

/**
 * Computes the Lambert W function of a symbolic expression.
 *
 * The Lambert W function W(x) is defined as the inverse of f(w) = w*e^w.
 * That is, if y = W(x), then y*e^y = x. This function returns the
 * principal branch W₀(x).
 *
 * @param x - The argument (as a symbolic expression).
 * @returns W(x) as a symbolic expression.
 *
 * @remarks
 * - W(0) = 0
 * - W(e) = 1 (since 1*e^1 = e)
 * - W(x)*exp(W(x)) = x (definition)
 * - W(-1/e) = -1 (branch point)
 * - W(x) is real for x ≥ -1/e
 * - d/dx W(x) = W(x)/(x*(1 + W(x))) for x ≠ 0
 * - Useful for solving equations like x*e^x = a
 *
 * @example
 * ```typescript
 * import { Symbol, Integer, lambertw, E } from 'symwasm';
 *
 * const x = new Symbol('x');
 *
 * lambertw(x);                    // W(x)
 * lambertw(new Integer(0));       // 0
 * lambertw(E);                    // 1
 *
 * // To solve x*e^x = 5, the answer is x = W(5)
 * // lambertw(new Integer(5)) gives the solution
 * ```
 *
 * @see {@link exp} - Exponential function
 * @see {@link log} - Natural logarithm
 */
export function lambertw(x: Expr): Expr {
  const wasm = getWasmModule();
  const obj = createBasic();
  try {
    const code = wasm._basic_lambertw(obj.getPtr(), x.getWasmPtr());
    checkException(code);
    return exprFromWasm(obj);
  } catch (e) {
    obj.free();
    throw e;
  }
}

// ============================================================================
// Other Mathematical Functions
// ============================================================================

/**
 * Computes the absolute value of a symbolic expression.
 *
 * Returns |x|, the magnitude of x. For real numbers, this is the non-negative
 * value. For complex numbers, this is the modulus sqrt(re(x)^2 + im(x)^2).
 *
 * @param x - The expression to take the absolute value of.
 * @returns |x| as a symbolic expression.
 *
 * @remarks
 * - abs(x) ≥ 0 for all x
 * - abs(x) = x for x ≥ 0
 * - abs(x) = -x for x < 0
 * - abs(-x) = abs(x)
 * - abs(x*y) = abs(x)*abs(y)
 * - abs(x/y) = abs(x)/abs(y)
 * - For complex z = a + bi: abs(z) = sqrt(a² + b²)
 * - d/dx abs(x) = sign(x) for real x ≠ 0
 *
 * @example
 * ```typescript
 * import { Symbol, Integer, abs, I, add, mul } from 'symwasm';
 *
 * const x = new Symbol('x');
 *
 * abs(x);                         // |x|
 * abs(new Integer(-5));           // 5
 * abs(new Integer(5));            // 5
 *
 * // Complex absolute value: |3 + 4i| = 5
 * abs(add(new Integer(3), mul(new Integer(4), I)));  // 5
 * ```
 *
 * @see {@link sign} - Sign function
 * @see {@link re} - Real part of complex number
 * @see {@link im} - Imaginary part of complex number
 */
export function abs(x: Expr): Expr {
  const wasm = getWasmModule();
  const obj = createBasic();
  try {
    const code = wasm._basic_abs(obj.getPtr(), x.getWasmPtr());
    checkException(code);
    return exprFromWasm(obj);
  } catch (e) {
    obj.free();
    throw e;
  }
}

/**
 * Computes the sign (signum) function of a symbolic expression.
 *
 * Returns -1 if x < 0, 0 if x = 0, and +1 if x > 0. For complex numbers,
 * returns x/|x| (the unit complex number in the direction of x).
 *
 * @param x - The expression to compute the sign of.
 * @returns sign(x) as a symbolic expression.
 *
 * @remarks
 * - sign(x) = -1 for x < 0
 * - sign(0) = 0
 * - sign(x) = 1 for x > 0
 * - sign(-x) = -sign(x)
 * - abs(x) = x * sign(x)
 * - x = abs(x) * sign(x)
 * - For complex z: sign(z) = z/|z| (unit phasor)
 *
 * @example
 * ```typescript
 * import { Symbol, Integer, sign } from 'symwasm';
 *
 * const x = new Symbol('x');
 *
 * sign(x);                        // sign(x)
 * sign(new Integer(-5));          // -1
 * sign(new Integer(0));           // 0
 * sign(new Integer(5));           // 1
 * ```
 *
 * @see {@link abs} - Absolute value
 */
export function sign(x: Expr): Expr {
  const wasm = getWasmModule();
  const obj = createBasic();
  try {
    const code = wasm._basic_sign(obj.getPtr(), x.getWasmPtr());
    checkException(code);
    return exprFromWasm(obj);
  } catch (e) {
    obj.free();
    throw e;
  }
}

/**
 * Computes the floor function of a symbolic expression.
 *
 * Returns the greatest integer less than or equal to x, denoted ⌊x⌋.
 * Also known as the integer part or greatest integer function.
 *
 * @param x - The expression to compute the floor of.
 * @returns ⌊x⌋ as a symbolic expression.
 *
 * @remarks
 * - floor(x) ≤ x < floor(x) + 1
 * - floor(n) = n for integers n
 * - floor(2.7) = 2
 * - floor(-2.7) = -3
 * - floor(x + n) = floor(x) + n for integer n
 * - x - 1 < floor(x) ≤ x
 * - ceiling(x) = -floor(-x)
 *
 * @example
 * ```typescript
 * import { Symbol, Integer, Float, floor } from 'symwasm';
 *
 * const x = new Symbol('x');
 *
 * floor(x);                       // ⌊x⌋
 * floor(new Float(2.7));          // 2
 * floor(new Float(-2.7));         // -3
 * floor(new Integer(5));          // 5
 * ```
 *
 * @see {@link ceiling} - Ceiling function
 */
export function floor(x: Expr): Expr {
  const wasm = getWasmModule();
  const obj = createBasic();
  try {
    const code = wasm._basic_floor(obj.getPtr(), x.getWasmPtr());
    checkException(code);
    return exprFromWasm(obj);
  } catch (e) {
    obj.free();
    throw e;
  }
}

/**
 * Computes the ceiling function of a symbolic expression.
 *
 * Returns the smallest integer greater than or equal to x, denoted ⌈x⌉.
 * Also known as the least integer function.
 *
 * @param x - The expression to compute the ceiling of.
 * @returns ⌈x⌉ as a symbolic expression.
 *
 * @remarks
 * - ceiling(x) - 1 < x ≤ ceiling(x)
 * - ceiling(n) = n for integers n
 * - ceiling(2.3) = 3
 * - ceiling(-2.3) = -2
 * - ceiling(x + n) = ceiling(x) + n for integer n
 * - x ≤ ceiling(x) < x + 1
 * - floor(x) = -ceiling(-x)
 *
 * @example
 * ```typescript
 * import { Symbol, Integer, Float, ceiling } from 'symwasm';
 *
 * const x = new Symbol('x');
 *
 * ceiling(x);                     // ⌈x⌉
 * ceiling(new Float(2.3));        // 3
 * ceiling(new Float(-2.3));       // -2
 * ceiling(new Integer(5));        // 5
 * ```
 *
 * @see {@link floor} - Floor function
 */
export function ceiling(x: Expr): Expr {
  const wasm = getWasmModule();
  const obj = createBasic();
  try {
    const code = wasm._basic_ceiling(obj.getPtr(), x.getWasmPtr());
    checkException(code);
    return exprFromWasm(obj);
  } catch (e) {
    obj.free();
    throw e;
  }
}
