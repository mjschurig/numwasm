/**
 * Expression simplification and manipulation functions.
 *
 * This module provides functions for simplifying, expanding, and transforming
 * symbolic expressions. These operations help reduce expressions to simpler
 * forms or rewrite them in different equivalent representations.
 *
 * @module simplify
 *
 * @example
 * ```typescript
 * import { Symbol, expand, simplify, numer, denom, add, mul, pow, div } from 'symwasm';
 *
 * const x = new Symbol('x');
 * const y = new Symbol('y');
 *
 * // Expand (x + y)²
 * expand(pow(add(x, y), new Integer(2)));  // x² + 2xy + y²
 *
 * // Simplify sin²(x) + cos²(x)
 * simplify(add(pow(sin(x), 2), pow(cos(x), 2)));  // 1
 *
 * // Extract numerator and denominator
 * const frac = div(add(x, y), mul(x, y));
 * numer(frac);  // x + y
 * denom(frac);  // x*y
 * ```
 */

import { getWasmModule } from '../wasm-loader.js';
import { createBasic, checkException } from '../wasm-memory.js';
import type { Expr } from '../core/expr.js';
import { exprFromWasm } from '../core/expr-factory.js';

/**
 * Expands an expression by distributing multiplication over addition.
 *
 * Performs algebraic expansion of products and powers, converting expressions
 * like (a + b)² into a² + 2ab + b².
 *
 * @param expr - The expression to expand.
 * @returns The expanded expression.
 *
 * @remarks
 * Expansion rules applied:
 * - (a + b)^n expands using the binomial theorem
 * - a * (b + c) = a*b + a*c (distributive law)
 * - (a + b) * (c + d) = a*c + a*d + b*c + b*d
 *
 * Note that expand() does not simplify; use simplify() to combine like terms
 * after expansion if needed.
 *
 * @example
 * ```typescript
 * import { Symbol, Integer, expand, add, mul, pow } from 'symwasm';
 *
 * const x = new Symbol('x');
 * const y = new Symbol('y');
 *
 * // Expand (x + 1)²
 * expand(pow(add(x, new Integer(1)), new Integer(2)));
 * // Returns: x² + 2*x + 1
 *
 * // Expand (x + y)(x - y)
 * expand(mul(add(x, y), sub(x, y)));
 * // Returns: x² - y²
 *
 * // Expand x(x + y + z)
 * const z = new Symbol('z');
 * expand(mul(x, add(add(x, y), z)));
 * // Returns: x² + x*y + x*z
 * ```
 *
 * @see {@link simplify} - Simplify expressions using heuristics
 */
export function expand(expr: Expr): Expr {
  const wasm = getWasmModule();
  const result = createBasic();
  try {
    const code = wasm._basic_expand(result.getPtr(), expr.getWasmPtr());
    checkException(code);
    return exprFromWasm(result);
  } catch (e) {
    result.free();
    throw e;
  }
}

/**
 * Simplifies an expression using heuristic algorithms.
 *
 * Applies various simplification rules to reduce an expression to a simpler
 * equivalent form. This includes combining like terms, canceling common factors,
 * and applying algebraic identities.
 *
 * @param expr - The expression to simplify.
 * @returns The simplified expression.
 *
 * @remarks
 * Simplification is heuristic and may not always find the "simplest" form.
 * Different simplifications may be appropriate for different contexts.
 *
 * Operations that may be performed:
 * - Combining like terms: 2x + 3x → 5x
 * - Canceling common factors: (x² - 1)/(x - 1) → x + 1
 * - Applying identities: sin²(x) + cos²(x) → 1
 * - Reducing fractions: 6/4 → 3/2
 *
 * @example
 * ```typescript
 * import { Symbol, Integer, simplify, add, mul, pow, sin, cos } from 'symwasm';
 *
 * const x = new Symbol('x');
 *
 * // Simplify x + x + x
 * simplify(add(add(x, x), x));  // 3*x
 *
 * // Simplify sin²(x) + cos²(x)
 * simplify(add(pow(sin(x), new Integer(2)), pow(cos(x), new Integer(2))));
 * // Returns: 1
 *
 * // Simplify x/x
 * simplify(div(x, x));  // 1
 * ```
 *
 * @see {@link expand} - Expand products and powers
 * @see {@link trigsimp} - Trigonometric simplification
 * @see {@link powsimp} - Power simplification
 */
export function simplify(expr: Expr): Expr {
  const wasm = getWasmModule();
  const result = createBasic();
  try {
    const code = wasm._basic_simplify(result.getPtr(), expr.getWasmPtr());
    checkException(code);
    return exprFromWasm(result);
  } catch (e) {
    result.free();
    throw e;
  }
}

/**
 * Extracts the numerator of a rational expression.
 *
 * For an expression representing a fraction p/q, returns p.
 * For non-fractional expressions, returns the expression itself
 * (since any expression can be viewed as itself over 1).
 *
 * @param expr - The expression to extract the numerator from.
 * @returns The numerator of the expression.
 *
 * @remarks
 * The expression is converted to the form numerator/denominator internally,
 * and the numerator part is returned. For expressions that are not explicitly
 * fractions, the "fraction form" is computed first.
 *
 * @example
 * ```typescript
 * import { Symbol, Integer, numer, div, add, mul } from 'symwasm';
 *
 * const x = new Symbol('x');
 * const y = new Symbol('y');
 *
 * // Numerator of (x + y)/(x - y)
 * numer(div(add(x, y), sub(x, y)));  // x + y
 *
 * // Numerator of x (which is x/1)
 * numer(x);  // x
 *
 * // Numerator of 2/3
 * numer(new Rational(2, 3));  // 2
 * ```
 *
 * @see {@link denom} - Extract the denominator
 */
export function numer(expr: Expr): Expr {
  const wasm = getWasmModule();
  const numerResult = createBasic();
  const denomResult = createBasic();
  try {
    const code = wasm._basic_as_numer_denom(
      numerResult.getPtr(),
      denomResult.getPtr(),
      expr.getWasmPtr()
    );
    checkException(code);
    denomResult.free();
    return exprFromWasm(numerResult);
  } catch (e) {
    numerResult.free();
    denomResult.free();
    throw e;
  }
}

/**
 * Extracts the denominator of a rational expression.
 *
 * For an expression representing a fraction p/q, returns q.
 * For non-fractional expressions, returns 1 (since any expression
 * can be viewed as itself over 1).
 *
 * @param expr - The expression to extract the denominator from.
 * @returns The denominator of the expression.
 *
 * @remarks
 * The expression is converted to the form numerator/denominator internally,
 * and the denominator part is returned. For non-fractional expressions,
 * the denominator is 1.
 *
 * @example
 * ```typescript
 * import { Symbol, Integer, denom, div, add, sub } from 'symwasm';
 *
 * const x = new Symbol('x');
 * const y = new Symbol('y');
 *
 * // Denominator of (x + y)/(x - y)
 * denom(div(add(x, y), sub(x, y)));  // x - y
 *
 * // Denominator of x (which is x/1)
 * denom(x);  // 1
 *
 * // Denominator of 2/3
 * denom(new Rational(2, 3));  // 3
 * ```
 *
 * @see {@link numer} - Extract the numerator
 */
export function denom(expr: Expr): Expr {
  const wasm = getWasmModule();
  const numerResult = createBasic();
  const denomResult = createBasic();
  try {
    const code = wasm._basic_as_numer_denom(
      numerResult.getPtr(),
      denomResult.getPtr(),
      expr.getWasmPtr()
    );
    checkException(code);
    numerResult.free();
    return exprFromWasm(denomResult);
  } catch (e) {
    numerResult.free();
    denomResult.free();
    throw e;
  }
}

/**
 * Simplifies trigonometric expressions using trigonometric identities.
 *
 * Applies identities such as sin²(x) + cos²(x) = 1 to simplify
 * expressions involving trigonometric functions.
 *
 * @param expr - The expression to simplify.
 * @returns The simplified expression.
 *
 * @remarks
 * Currently delegates to the general simplify() function.
 * Common identities that may be applied:
 * - sin²(x) + cos²(x) = 1
 * - 1 + tan²(x) = sec²(x)
 * - 1 + cot²(x) = csc²(x)
 * - sin(2x) = 2sin(x)cos(x)
 * - cos(2x) = cos²(x) - sin²(x)
 *
 * @example
 * ```typescript
 * import { Symbol, trigsimp, add, pow, sin, cos, Integer } from 'symwasm';
 *
 * const x = new Symbol('x');
 *
 * // Simplify sin²(x) + cos²(x)
 * trigsimp(add(pow(sin(x), new Integer(2)), pow(cos(x), new Integer(2))));
 * // Returns: 1
 * ```
 *
 * @see {@link simplify} - General simplification
 * @see {@link rewrite_as_exp} - Rewrite trig as exponentials
 */
export function trigsimp(expr: Expr): Expr {
  return simplify(expr);
}

/**
 * Simplifies expressions involving radicals (roots).
 *
 * Rationalizes denominators and simplifies expressions containing
 * square roots, cube roots, and other radicals.
 *
 * @param expr - The expression to simplify.
 * @returns The simplified expression.
 *
 * @remarks
 * Currently delegates to the general simplify() function.
 * Operations that may be performed:
 * - Rationalize denominators: 1/√2 → √2/2
 * - Combine radicals: √2 * √3 → √6
 * - Simplify nested radicals
 * - Denest radicals when possible
 *
 * @example
 * ```typescript
 * import { Symbol, radsimp, sqrt, div, Integer, mul } from 'symwasm';
 *
 * // Simplify 1/√2
 * radsimp(div(new Integer(1), sqrt(new Integer(2))));
 * // Returns: √2/2
 *
 * // Simplify √8
 * radsimp(sqrt(new Integer(8)));
 * // Returns: 2√2
 * ```
 *
 * @see {@link simplify} - General simplification
 * @see {@link sqrt} - Square root function
 */
export function radsimp(expr: Expr): Expr {
  return simplify(expr);
}

/**
 * Simplifies expressions with powers by combining bases and exponents.
 *
 * Applies power rules to simplify expressions involving exponents,
 * such as combining x^a * x^b = x^(a+b).
 *
 * @param expr - The expression to simplify.
 * @returns The simplified expression.
 *
 * @remarks
 * Currently delegates to the general simplify() function.
 * Power rules that may be applied:
 * - x^a * x^b = x^(a+b)
 * - (x^a)^b = x^(a*b)
 * - (xy)^a = x^a * y^a
 * - x^a / x^b = x^(a-b)
 * - x^0 = 1 (for x ≠ 0)
 * - x^1 = x
 *
 * @example
 * ```typescript
 * import { Symbol, Integer, powsimp, mul, pow } from 'symwasm';
 *
 * const x = new Symbol('x');
 *
 * // Simplify x² * x³
 * powsimp(mul(pow(x, new Integer(2)), pow(x, new Integer(3))));
 * // Returns: x^5
 *
 * // Simplify (x²)³
 * powsimp(pow(pow(x, new Integer(2)), new Integer(3)));
 * // Returns: x^6
 * ```
 *
 * @see {@link simplify} - General simplification
 * @see {@link pow} - Power function
 */
export function powsimp(expr: Expr): Expr {
  return simplify(expr);
}

/**
 * Rewrites trigonometric functions in terms of complex exponentials.
 *
 * Uses Euler's formula to convert sine and cosine to exponential form:
 * - sin(x) = (e^(ix) - e^(-ix)) / (2i)
 * - cos(x) = (e^(ix) + e^(-ix)) / 2
 *
 * @param expr - The expression to rewrite.
 * @returns The expression with trig functions replaced by exponentials.
 *
 * @remarks
 * This transformation is useful for:
 * - Integration of trigonometric functions
 * - Proving trigonometric identities
 * - Simplifying products of sines and cosines
 * - Converting to a form suitable for Fourier analysis
 *
 * Euler's formula: e^(ix) = cos(x) + i*sin(x)
 *
 * @example
 * ```typescript
 * import { Symbol, rewrite_as_exp, sin, cos, add } from 'symwasm';
 *
 * const x = new Symbol('x');
 *
 * // Rewrite sin(x) as exponentials
 * rewrite_as_exp(sin(x));
 * // Returns: (e^(i*x) - e^(-i*x))/(2*i)
 *
 * // Rewrite cos(x) as exponentials
 * rewrite_as_exp(cos(x));
 * // Returns: (e^(i*x) + e^(-i*x))/2
 * ```
 *
 * @see {@link rewrite_as_sin} - Rewrite in terms of sine
 * @see {@link rewrite_as_cos} - Rewrite in terms of cosine
 */
export function rewrite_as_exp(expr: Expr): Expr {
  const wasm = getWasmModule();
  const result = createBasic();
  checkException(wasm._basic_rewrite_as_exp(result.getPtr(), expr.getWasmPtr()));
  return exprFromWasm(result);
}

/**
 * Rewrites trigonometric functions in terms of sine only.
 *
 * Converts all trigonometric functions in an expression to use only
 * the sine function, using identities like cos(x) = sin(π/2 - x).
 *
 * @param expr - The expression to rewrite.
 * @returns The expression with all trig functions written as sines.
 *
 * @remarks
 * Conversion identities used:
 * - cos(x) = sin(π/2 - x)
 * - tan(x) = sin(x)/sin(π/2 - x)
 * - cot(x) = sin(π/2 - x)/sin(x)
 * - sec(x) = 1/sin(π/2 - x)
 * - csc(x) = 1/sin(x)
 *
 * @example
 * ```typescript
 * import { Symbol, rewrite_as_sin, cos, tan } from 'symwasm';
 *
 * const x = new Symbol('x');
 *
 * // Rewrite cos(x) in terms of sin
 * rewrite_as_sin(cos(x));
 * // Returns expression equivalent to sin(π/2 - x)
 *
 * // Rewrite tan(x) in terms of sin
 * rewrite_as_sin(tan(x));
 * // Returns sin(x)/sin(π/2 - x)
 * ```
 *
 * @see {@link rewrite_as_cos} - Rewrite in terms of cosine
 * @see {@link rewrite_as_exp} - Rewrite as exponentials
 */
export function rewrite_as_sin(expr: Expr): Expr {
  const wasm = getWasmModule();
  const result = createBasic();
  checkException(wasm._basic_rewrite_as_sin(result.getPtr(), expr.getWasmPtr()));
  return exprFromWasm(result);
}

/**
 * Rewrites trigonometric functions in terms of cosine only.
 *
 * Converts all trigonometric functions in an expression to use only
 * the cosine function, using identities like sin(x) = cos(π/2 - x).
 *
 * @param expr - The expression to rewrite.
 * @returns The expression with all trig functions written as cosines.
 *
 * @remarks
 * Conversion identities used:
 * - sin(x) = cos(π/2 - x)
 * - tan(x) = cos(π/2 - x)/cos(x)
 * - cot(x) = cos(x)/cos(π/2 - x)
 * - sec(x) = 1/cos(x)
 * - csc(x) = 1/cos(π/2 - x)
 *
 * @example
 * ```typescript
 * import { Symbol, rewrite_as_cos, sin, tan } from 'symwasm';
 *
 * const x = new Symbol('x');
 *
 * // Rewrite sin(x) in terms of cos
 * rewrite_as_cos(sin(x));
 * // Returns expression equivalent to cos(π/2 - x)
 *
 * // Rewrite tan(x) in terms of cos
 * rewrite_as_cos(tan(x));
 * // Returns cos(π/2 - x)/cos(x)
 * ```
 *
 * @see {@link rewrite_as_sin} - Rewrite in terms of sine
 * @see {@link rewrite_as_exp} - Rewrite as exponentials
 */
export function rewrite_as_cos(expr: Expr): Expr {
  const wasm = getWasmModule();
  const result = createBasic();
  checkException(wasm._basic_rewrite_as_cos(result.getPtr(), expr.getWasmPtr()));
  return exprFromWasm(result);
}

/**
 * Extracts the real and imaginary parts of a complex expression.
 *
 * Separates a complex expression into its real and imaginary components,
 * returning both as separate expressions.
 *
 * @param expr - The complex expression to decompose.
 * @returns An object with `real` and `imag` properties containing the respective parts.
 *
 * @remarks
 * For a complex expression z = a + bi:
 * - real part = a
 * - imaginary part = b (not bi)
 *
 * This function is useful for:
 * - Analyzing complex-valued functions
 * - Computing absolute value |z| = √(real² + imag²)
 * - Computing argument arg(z) = atan2(imag, real)
 *
 * @example
 * ```typescript
 * import { Symbol, Integer, as_real_imag, add, mul, I, exp } from 'symwasm';
 *
 * const x = new Symbol('x');
 *
 * // Decompose 3 + 4i
 * const z1 = add(new Integer(3), mul(new Integer(4), I));
 * const parts1 = as_real_imag(z1);
 * // parts1.real = 3, parts1.imag = 4
 *
 * // Decompose e^(ix)
 * const z2 = exp(mul(I, x));
 * const parts2 = as_real_imag(z2);
 * // parts2.real = cos(x), parts2.imag = sin(x)
 * ```
 *
 * @see {@link re} - Extract just the real part
 * @see {@link im} - Extract just the imaginary part
 * @see {@link expand_complex} - Alias for this function
 */
export function as_real_imag(expr: Expr): { real: Expr; imag: Expr } {
  const wasm = getWasmModule();
  const realResult = createBasic();
  const imagResult = createBasic();
  checkException(wasm._basic_as_real_imag(realResult.getPtr(), imagResult.getPtr(), expr.getWasmPtr()));
  return {
    real: exprFromWasm(realResult),
    imag: exprFromWasm(imagResult),
  };
}

/**
 * Expands trigonometric functions using exponential identities.
 *
 * Converts trigonometric functions to their exponential form and then
 * expands the result. Useful for simplifying products of trig functions.
 *
 * @param expr - The expression to expand.
 * @returns The expanded expression.
 *
 * @remarks
 * This is a two-step process:
 * 1. Rewrite trig functions as exponentials using Euler's formula
 * 2. Expand the resulting expression
 *
 * Useful for:
 * - Computing products like sin(a)cos(b)
 * - Deriving multiple-angle formulas
 * - Converting between trig forms
 *
 * @example
 * ```typescript
 * import { Symbol, expand_trig, mul, sin, cos } from 'symwasm';
 *
 * const a = new Symbol('a');
 * const b = new Symbol('b');
 *
 * // Expand sin(a)cos(b)
 * expand_trig(mul(sin(a), cos(b)));
 * // Useful for product-to-sum conversion
 * ```
 *
 * @see {@link rewrite_as_exp} - Convert trig to exponentials
 * @see {@link expand} - General expansion
 * @see {@link trigsimp} - Simplify trigonometric expressions
 */
export function expand_trig(expr: Expr): Expr {
  return expand(rewrite_as_exp(expr));
}

/**
 * Expands a complex expression into its real and imaginary parts.
 *
 * Alias for {@link as_real_imag}. Separates a complex expression into
 * its real and imaginary components.
 *
 * @param expr - The complex expression to expand.
 * @returns An object with `real` and `imag` properties containing the respective parts.
 *
 * @example
 * ```typescript
 * import { Symbol, expand_complex, exp, mul, I } from 'symwasm';
 *
 * const x = new Symbol('x');
 *
 * // Expand e^(ix) into real and imaginary parts
 * const parts = expand_complex(exp(mul(I, x)));
 * // parts.real = cos(x)
 * // parts.imag = sin(x)
 * ```
 *
 * @see {@link as_real_imag} - Same function
 * @see {@link re} - Extract just the real part
 * @see {@link im} - Extract just the imaginary part
 */
export function expand_complex(expr: Expr): { real: Expr; imag: Expr } {
  return as_real_imag(expr);
}

/**
 * Result of common subexpression elimination.
 *
 * Contains the replacement symbols, their corresponding expressions,
 * and the reduced expressions that use those symbols.
 */
export interface CSEResult {
  /**
   * The generated symbols used as replacements (e.g., x0, x1, x2, ...).
   */
  replacementSymbols: Expr[];
  /**
   * The subexpressions that each replacement symbol stands for.
   * replacementSymbols[i] = replacementExprs[i]
   */
  replacementExprs: Expr[];
  /**
   * The original expressions rewritten using the replacement symbols.
   */
  reducedExprs: Expr[];
}

/**
 * Performs common subexpression elimination on a list of expressions.
 *
 * Identifies repeated subexpressions and replaces them with generated symbols
 * to reduce redundant computation. This is useful for optimizing numerical
 * evaluation of complex expressions.
 *
 * @param exprs - The expressions to analyze for common subexpressions.
 * @returns An object containing:
 *   - `replacementSymbols`: Generated symbols (x0, x1, ...) for common subexpressions
 *   - `replacementExprs`: The subexpressions each symbol replaces
 *   - `reducedExprs`: The original expressions rewritten using the symbols
 *
 * @remarks
 * Common subexpression elimination (CSE) is a compiler optimization technique
 * that identifies expressions that are computed multiple times and replaces
 * them with a single computed value.
 *
 * For example, if you have:
 * - expr1 = sin(x) + cos(x)
 * - expr2 = sin(x) * cos(x)
 *
 * CSE might produce:
 * - x0 = sin(x)
 * - x1 = cos(x)
 * - expr1 = x0 + x1
 * - expr2 = x0 * x1
 *
 * This is particularly useful when:
 * - Generating code for numerical evaluation
 * - Optimizing expressions for repeated evaluation
 * - Reducing the computational cost of complex expressions
 *
 * @example
 * ```typescript
 * import { Symbol, Integer, cse, sin, cos, add, mul, pow } from 'symwasm';
 *
 * const x = new Symbol('x');
 *
 * // Two expressions with common subexpressions
 * const expr1 = add(sin(x), cos(x));
 * const expr2 = mul(sin(x), cos(x));
 *
 * const result = cse([expr1, expr2]);
 * // result.replacementSymbols = [x0, x1]  (generated symbols)
 * // result.replacementExprs = [sin(x), cos(x)]  (what they represent)
 * // result.reducedExprs = [x0 + x1, x0 * x1]  (rewritten expressions)
 * ```
 *
 * @example
 * ```typescript
 * import { Symbol, Integer, cse, pow, add, mul } from 'symwasm';
 *
 * const x = new Symbol('x');
 * const y = new Symbol('y');
 *
 * // Expression with repeated subexpression (x + y)
 * const expr = add(pow(add(x, y), new Integer(2)), mul(add(x, y), new Integer(3)));
 *
 * const result = cse([expr]);
 * // x0 = x + y
 * // reduced = x0^2 + 3*x0
 * ```
 *
 * @see {@link simplify} - General simplification
 * @see {@link expand} - Expand products and powers
 */
export function cse(exprs: Expr[]): CSEResult {
  const wasm = getWasmModule();

  // Create CVecBasic containers for inputs and outputs
  const inputVec = wasm._vecbasic_new();
  const replacementSymsVec = wasm._vecbasic_new();
  const replacementExprsVec = wasm._vecbasic_new();
  const reducedExprsVec = wasm._vecbasic_new();

  try {
    // Add input expressions to the input vector
    for (const expr of exprs) {
      wasm._vecbasic_push_back(inputVec, expr.getWasmPtr());
    }

    // Call the CSE function
    const code = wasm._basic_cse(
      replacementSymsVec,
      replacementExprsVec,
      reducedExprsVec,
      inputVec
    );
    checkException(code);

    // Extract replacement symbols
    const replacementSymbols: Expr[] = [];
    const numReplacements = wasm._vecbasic_size(replacementSymsVec);
    for (let i = 0; i < numReplacements; i++) {
      const temp = createBasic();
      wasm._vecbasic_get(replacementSymsVec, i, temp.getPtr());
      replacementSymbols.push(exprFromWasm(temp));
    }

    // Extract replacement expressions
    const replacementExprs: Expr[] = [];
    for (let i = 0; i < numReplacements; i++) {
      const temp = createBasic();
      wasm._vecbasic_get(replacementExprsVec, i, temp.getPtr());
      replacementExprs.push(exprFromWasm(temp));
    }

    // Extract reduced expressions
    const reducedExprs: Expr[] = [];
    const numReduced = wasm._vecbasic_size(reducedExprsVec);
    for (let i = 0; i < numReduced; i++) {
      const temp = createBasic();
      wasm._vecbasic_get(reducedExprsVec, i, temp.getPtr());
      reducedExprs.push(exprFromWasm(temp));
    }

    return {
      replacementSymbols,
      replacementExprs,
      reducedExprs,
    };
  } finally {
    // Free all CVecBasic containers
    wasm._vecbasic_free(inputVec);
    wasm._vecbasic_free(replacementSymsVec);
    wasm._vecbasic_free(replacementExprsVec);
    wasm._vecbasic_free(reducedExprsVec);
  }
}
