/**
 * Special mathematical functions for symbolic expressions.
 *
 * This module provides special functions commonly used in mathematical physics,
 * statistics, and number theory, including gamma functions, error functions,
 * and the Riemann zeta function. All functions work symbolically and can be
 * differentiated, simplified, and evaluated numerically.
 *
 * @module functions/special
 *
 * @example
 * ```typescript
 * import { Symbol, Integer, gamma, erf, zeta, diff } from 'symwasm';
 *
 * const x = new Symbol('x');
 * const n = new Integer(5);
 *
 * // Gamma function: Γ(5) = 4! = 24
 * gamma(n);
 *
 * // Error function for probability
 * erf(x);
 *
 * // Riemann zeta function
 * zeta(new Integer(2));  // π²/6
 * ```
 */

import { getWasmModule } from '../../wasm-loader.js';
import { createBasic, checkException } from '../../wasm-memory.js';
import type { Expr } from '../../core/expr.js';
import { exprFromWasm } from '../../core/expr-factory.js';

/**
 * Computes the gamma function of a symbolic expression.
 *
 * The gamma function Γ(x) extends the factorial to complex numbers.
 * For positive integers, Γ(n) = (n-1)!. It satisfies the recurrence
 * relation Γ(x+1) = x·Γ(x).
 *
 * @param x - The argument (as a symbolic expression).
 * @returns Γ(x) as a symbolic expression.
 *
 * @remarks
 * - Γ(1) = 1
 * - Γ(n) = (n-1)! for positive integers n
 * - Γ(1/2) = √π
 * - Γ(x+1) = x·Γ(x) (recurrence relation)
 * - Γ(x)·Γ(1-x) = π/sin(πx) (reflection formula)
 * - Γ(x) has poles at x = 0, -1, -2, -3, ...
 * - d/dx Γ(x) = Γ(x)·ψ(x) where ψ is the digamma function
 *
 * @example
 * ```typescript
 * import { Symbol, Integer, Rational, gamma } from 'symwasm';
 *
 * // Γ(5) = 4! = 24
 * gamma(new Integer(5));  // 24
 *
 * // Γ(1/2) = √π
 * gamma(new Rational(1, 2));  // √π
 *
 * // Symbolic gamma
 * const x = new Symbol('x');
 * gamma(x);  // Γ(x)
 * ```
 *
 * @see {@link loggamma} - Logarithm of gamma function
 * @see {@link digamma} - Derivative of log-gamma
 * @see {@link beta} - Beta function (related to gamma)
 */
export function gamma(x: Expr): Expr {
  const wasm = getWasmModule();
  const obj = createBasic();
  try {
    const code = wasm._basic_gamma(obj.getPtr(), x.getWasmPtr());
    checkException(code);
    return exprFromWasm(obj);
  } catch (e) {
    obj.free();
    throw e;
  }
}

/**
 * Computes the log-gamma function of a symbolic expression.
 *
 * The log-gamma function is defined as log(Γ(x)), which is useful for
 * numerical computation when Γ(x) would overflow, and for computing
 * ratios of gamma functions.
 *
 * @param x - The argument (as a symbolic expression).
 * @returns log(Γ(x)) as a symbolic expression.
 *
 * @remarks
 * - loggamma(1) = 0
 * - loggamma(x) is defined for x > 0
 * - loggamma(n) = log((n-1)!) for positive integers
 * - loggamma(x+1) = loggamma(x) + log(x)
 * - d/dx loggamma(x) = ψ(x) (digamma function)
 * - More numerically stable than log(gamma(x)) for large x
 *
 * @example
 * ```typescript
 * import { Symbol, Integer, loggamma } from 'symwasm';
 *
 * // loggamma(5) = log(24)
 * loggamma(new Integer(5));
 *
 * // Symbolic usage
 * const x = new Symbol('x');
 * loggamma(x);  // log(Γ(x))
 * ```
 *
 * @see {@link gamma} - Gamma function
 * @see {@link digamma} - Derivative of log-gamma
 */
export function loggamma(x: Expr): Expr {
  const wasm = getWasmModule();
  const obj = createBasic();
  try {
    const code = wasm._basic_loggamma(obj.getPtr(), x.getWasmPtr());
    checkException(code);
    return exprFromWasm(obj);
  } catch (e) {
    obj.free();
    throw e;
  }
}

/**
 * Computes the digamma function of a symbolic expression.
 *
 * The digamma function ψ(x) is the logarithmic derivative of the gamma function:
 * ψ(x) = d/dx ln(Γ(x)) = Γ'(x)/Γ(x). It is also known as the psi function.
 *
 * @param x - The argument (as a symbolic expression).
 * @returns ψ(x) as a symbolic expression.
 *
 * @remarks
 * - ψ(1) = -γ (negative Euler-Mascheroni constant)
 * - ψ(n) = -γ + Σ(k=1 to n-1) 1/k for positive integers n
 * - ψ(x+1) = ψ(x) + 1/x (recurrence relation)
 * - ψ(1-x) - ψ(x) = π·cot(πx) (reflection formula)
 * - d/dx ψ(x) = ψ¹(x) (trigamma function)
 *
 * @example
 * ```typescript
 * import { Symbol, Integer, digamma } from 'symwasm';
 *
 * // ψ(1) = -γ ≈ -0.5772
 * digamma(new Integer(1));
 *
 * // Symbolic usage
 * const x = new Symbol('x');
 * digamma(x);  // ψ(x)
 * ```
 *
 * @see {@link gamma} - Gamma function
 * @see {@link loggamma} - Log-gamma function
 * @see {@link polygamma} - Higher derivatives of digamma
 */
export function digamma(x: Expr): Expr {
  const wasm = getWasmModule();
  const obj = createBasic();
  try {
    const code = wasm._basic_digamma(obj.getPtr(), x.getWasmPtr());
    checkException(code);
    return exprFromWasm(obj);
  } catch (e) {
    obj.free();
    throw e;
  }
}

/**
 * Computes the polygamma function of symbolic expressions.
 *
 * The polygamma function ψ^(n)(x) is the (n+1)th derivative of the log-gamma function,
 * or equivalently the nth derivative of the digamma function:
 * ψ^(n)(x) = d^n/dx^n ψ(x) = d^(n+1)/dx^(n+1) ln(Γ(x))
 *
 * @param n - The order of the derivative (n ≥ 0).
 * @param x - The argument at which to evaluate.
 * @returns ψ^(n)(x) as a symbolic expression.
 *
 * @remarks
 * - polygamma(0, x) = digamma(x) = ψ(x)
 * - polygamma(1, x) = trigamma(x) = ψ¹(x)
 * - polygamma(2, x) = tetragamma(x) = ψ²(x)
 * - ψ^(n)(x) = (-1)^(n+1) n! Σ(k=0 to ∞) 1/(x+k)^(n+1) for n ≥ 1
 * - ψ^(n)(1) = (-1)^(n+1) n! ζ(n+1) (relates to zeta function)
 *
 * @example
 * ```typescript
 * import { Symbol, Integer, polygamma, digamma } from 'symwasm';
 *
 * const x = new Symbol('x');
 *
 * // Digamma is polygamma of order 0
 * polygamma(new Integer(0), x);  // Same as digamma(x)
 *
 * // Trigamma function (second derivative of log-gamma)
 * polygamma(new Integer(1), x);  // ψ¹(x)
 *
 * // At x=1: ψ¹(1) = π²/6 = ζ(2)
 * polygamma(new Integer(1), new Integer(1));
 * ```
 *
 * @see {@link digamma} - Digamma function (polygamma of order 0)
 * @see {@link gamma} - Gamma function
 * @see {@link zeta} - Riemann zeta function
 */
export function polygamma(n: Expr, x: Expr): Expr {
  const wasm = getWasmModule();
  const obj = createBasic();
  try {
    const code = wasm._basic_polygamma(obj.getPtr(), n.getWasmPtr(), x.getWasmPtr());
    checkException(code);
    return exprFromWasm(obj);
  } catch (e) {
    obj.free();
    throw e;
  }
}

/**
 * Computes the beta function of symbolic expressions.
 *
 * The beta function B(a, b) is defined as the integral of t^(a-1)(1-t)^(b-1) from 0 to 1,
 * and can be expressed in terms of the gamma function as B(a, b) = Γ(a)Γ(b)/Γ(a+b).
 *
 * @param a - The first parameter (as a symbolic expression).
 * @param b - The second parameter (as a symbolic expression).
 * @returns B(a, b) = Γ(a)Γ(b)/Γ(a+b) as a symbolic expression.
 *
 * @remarks
 * - B(a, b) = B(b, a) (symmetric)
 * - B(1, 1) = 1
 * - B(a, 1) = 1/a
 * - B(m, n) = (m-1)!(n-1)!/(m+n-1)! for positive integers
 * - B(a, b) = ∫₀¹ t^(a-1)(1-t)^(b-1) dt
 * - Important in probability theory (beta distribution)
 *
 * @example
 * ```typescript
 * import { Symbol, Integer, beta } from 'symwasm';
 *
 * // B(2, 3) = Γ(2)Γ(3)/Γ(5) = 1!2!/4! = 1/12
 * beta(new Integer(2), new Integer(3));
 *
 * // Symbolic usage
 * const a = new Symbol('a');
 * const b = new Symbol('b');
 * beta(a, b);  // B(a, b)
 * ```
 *
 * @see {@link gamma} - Gamma function
 * @see {@link lowergamma} - Incomplete gamma function
 */
export function beta(a: Expr, b: Expr): Expr {
  const wasm = getWasmModule();
  const obj = createBasic();
  try {
    const code = wasm._basic_beta(obj.getPtr(), a.getWasmPtr(), b.getWasmPtr());
    checkException(code);
    return exprFromWasm(obj);
  } catch (e) {
    obj.free();
    throw e;
  }
}

/**
 * Computes the lower incomplete gamma function of symbolic expressions.
 *
 * The lower incomplete gamma function γ(s, x) is defined as the integral
 * of t^(s-1)e^(-t) from 0 to x. It represents the cumulative distribution
 * function of the gamma distribution (up to normalization).
 *
 * @param s - The shape parameter (as a symbolic expression).
 * @param x - The upper limit of integration (as a symbolic expression).
 * @returns γ(s, x) = ∫₀ˣ t^(s-1)e^(-t) dt as a symbolic expression.
 *
 * @remarks
 * - γ(s, 0) = 0
 * - γ(s, x) → Γ(s) as x → ∞
 * - γ(s, x) + Γ(s, x) = Γ(s) (relation to upper incomplete gamma)
 * - γ(1, x) = 1 - e^(-x)
 * - γ(n, x) = (n-1)! × (1 - e^(-x) × Σ(k=0 to n-1) x^k/k!) for positive integers n
 * - Used in the regularized incomplete gamma function P(s, x) = γ(s, x)/Γ(s)
 *
 * @example
 * ```typescript
 * import { Symbol, Integer, lowergamma } from 'symwasm';
 *
 * const s = new Symbol('s');
 * const x = new Symbol('x');
 *
 * // General symbolic form
 * lowergamma(s, x);  // γ(s, x)
 *
 * // γ(1, x) = 1 - e^(-x)
 * lowergamma(new Integer(1), x);
 * ```
 *
 * @see {@link uppergamma} - Upper incomplete gamma function
 * @see {@link gamma} - Complete gamma function
 */
export function lowergamma(s: Expr, x: Expr): Expr {
  const wasm = getWasmModule();
  const obj = createBasic();
  try {
    const code = wasm._basic_lowergamma(obj.getPtr(), s.getWasmPtr(), x.getWasmPtr());
    checkException(code);
    return exprFromWasm(obj);
  } catch (e) {
    obj.free();
    throw e;
  }
}

/**
 * Computes the upper incomplete gamma function of symbolic expressions.
 *
 * The upper incomplete gamma function Γ(s, x) is defined as the integral
 * of t^(s-1)e^(-t) from x to infinity. It's the complement of the lower
 * incomplete gamma function.
 *
 * @param s - The shape parameter (as a symbolic expression).
 * @param x - The lower limit of integration (as a symbolic expression).
 * @returns Γ(s, x) = ∫ₓ^∞ t^(s-1)e^(-t) dt as a symbolic expression.
 *
 * @remarks
 * - Γ(s, 0) = Γ(s)
 * - Γ(s, x) → 0 as x → ∞
 * - γ(s, x) + Γ(s, x) = Γ(s) (relation to lower incomplete gamma)
 * - Γ(1, x) = e^(-x)
 * - Γ(s+1, x) = s·Γ(s, x) + x^s·e^(-x) (recurrence)
 * - Used in the complementary regularized incomplete gamma function Q(s, x) = Γ(s, x)/Γ(s)
 *
 * @example
 * ```typescript
 * import { Symbol, Integer, uppergamma } from 'symwasm';
 *
 * const s = new Symbol('s');
 * const x = new Symbol('x');
 *
 * // General symbolic form
 * uppergamma(s, x);  // Γ(s, x)
 *
 * // Γ(1, x) = e^(-x)
 * uppergamma(new Integer(1), x);
 * ```
 *
 * @see {@link lowergamma} - Lower incomplete gamma function
 * @see {@link gamma} - Complete gamma function
 */
export function uppergamma(s: Expr, x: Expr): Expr {
  const wasm = getWasmModule();
  const obj = createBasic();
  try {
    const code = wasm._basic_uppergamma(obj.getPtr(), s.getWasmPtr(), x.getWasmPtr());
    checkException(code);
    return exprFromWasm(obj);
  } catch (e) {
    obj.free();
    throw e;
  }
}

/**
 * Computes the error function of a symbolic expression.
 *
 * The error function erf(x) is defined as (2/√π) × ∫₀ˣ e^(-t²) dt.
 * It is the integral of the Gaussian distribution and is widely used
 * in probability, statistics, and partial differential equations.
 *
 * @param x - The argument (as a symbolic expression).
 * @returns erf(x) = (2/√π) × ∫₀ˣ e^(-t²) dt as a symbolic expression.
 *
 * @remarks
 * - erf(0) = 0
 * - erf(-x) = -erf(x) (odd function)
 * - erf(x) → 1 as x → +∞
 * - erf(x) → -1 as x → -∞
 * - erf(x) + erfc(x) = 1
 * - d/dx erf(x) = (2/√π)e^(-x²)
 * - The probability P(|X| < x) for standard normal X is erf(x/√2)
 *
 * @example
 * ```typescript
 * import { Symbol, Integer, erf } from 'symwasm';
 *
 * const x = new Symbol('x');
 *
 * erf(x);                 // erf(x)
 * erf(new Integer(0));    // 0
 *
 * // Probability that |X| < 1 for standard normal ≈ erf(1/√2) ≈ 0.683
 * ```
 *
 * @see {@link erfc} - Complementary error function
 */
export function erf(x: Expr): Expr {
  const wasm = getWasmModule();
  const obj = createBasic();
  try {
    const code = wasm._basic_erf(obj.getPtr(), x.getWasmPtr());
    checkException(code);
    return exprFromWasm(obj);
  } catch (e) {
    obj.free();
    throw e;
  }
}

/**
 * Computes the complementary error function of a symbolic expression.
 *
 * The complementary error function erfc(x) is defined as 1 - erf(x),
 * or equivalently (2/√π) × ∫ₓ^∞ e^(-t²) dt. It is more accurate than
 * computing 1 - erf(x) directly for large x.
 *
 * @param x - The argument (as a symbolic expression).
 * @returns erfc(x) = 1 - erf(x) as a symbolic expression.
 *
 * @remarks
 * - erfc(0) = 1
 * - erfc(-x) = 2 - erfc(x)
 * - erfc(x) → 0 as x → +∞
 * - erfc(x) → 2 as x → -∞
 * - erfc(x) + erf(x) = 1
 * - d/dx erfc(x) = -(2/√π)e^(-x²)
 * - More numerically stable than 1 - erf(x) for large positive x
 *
 * @example
 * ```typescript
 * import { Symbol, Integer, erfc } from 'symwasm';
 *
 * const x = new Symbol('x');
 *
 * erfc(x);                // erfc(x)
 * erfc(new Integer(0));   // 1
 *
 * // For large x, erfc(x) ≈ e^(-x²)/(x√π)
 * ```
 *
 * @see {@link erf} - Error function
 */
export function erfc(x: Expr): Expr {
  const wasm = getWasmModule();
  const obj = createBasic();
  try {
    const code = wasm._basic_erfc(obj.getPtr(), x.getWasmPtr());
    checkException(code);
    return exprFromWasm(obj);
  } catch (e) {
    obj.free();
    throw e;
  }
}

/**
 * Computes the Riemann zeta function of a symbolic expression.
 *
 * The Riemann zeta function ζ(s) is defined as Σ(n=1 to ∞) 1/n^s for Re(s) > 1,
 * and extended to other values by analytic continuation. It is fundamental
 * in number theory and has deep connections to the distribution of prime numbers.
 *
 * @param s - The argument (as a symbolic expression).
 * @returns ζ(s) as a symbolic expression.
 *
 * @remarks
 * - ζ(2) = π²/6 (Basel problem)
 * - ζ(4) = π⁴/90
 * - ζ(2n) = (-1)^(n+1) × B_{2n} × (2π)^{2n} / (2 × (2n)!) for positive integers n
 * - ζ(-n) = -B_{n+1}/(n+1) for non-negative integers n (Bernoulli numbers)
 * - ζ(s) has a simple pole at s = 1
 * - ζ(s) = 0 at s = -2, -4, -6, ... (trivial zeros)
 * - The Riemann Hypothesis concerns non-trivial zeros with Re(s) = 1/2
 *
 * @example
 * ```typescript
 * import { Symbol, Integer, zeta } from 'symwasm';
 *
 * // ζ(2) = π²/6
 * zeta(new Integer(2));
 *
 * // ζ(4) = π⁴/90
 * zeta(new Integer(4));
 *
 * // Symbolic usage
 * const s = new Symbol('s');
 * zeta(s);
 * ```
 *
 * @see {@link dirichlet_eta} - Dirichlet eta function
 */
export function zeta(s: Expr): Expr {
  const wasm = getWasmModule();
  const obj = createBasic();
  try {
    const code = wasm._basic_zeta(obj.getPtr(), s.getWasmPtr());
    checkException(code);
    return exprFromWasm(obj);
  } catch (e) {
    obj.free();
    throw e;
  }
}

/**
 * Computes the Dirichlet eta function of a symbolic expression.
 *
 * The Dirichlet eta function η(s) is defined as Σ(n=1 to ∞) (-1)^(n-1)/n^s,
 * which is the alternating version of the Riemann zeta function.
 * It is related to the zeta function by η(s) = (1 - 2^(1-s)) × ζ(s).
 *
 * @param s - The argument (as a symbolic expression).
 * @returns η(s) as a symbolic expression.
 *
 * @remarks
 * - η(s) = (1 - 2^(1-s)) × ζ(s)
 * - η(1) = ln(2) (alternating harmonic series)
 * - η(2) = π²/12
 * - η(s) converges for Re(s) > 0, unlike ζ(s) which requires Re(s) > 1
 * - η(s) is entire (no poles), unlike ζ(s)
 * - Also known as the alternating zeta function
 *
 * @example
 * ```typescript
 * import { Symbol, Integer, dirichlet_eta } from 'symwasm';
 *
 * // η(1) = ln(2) (alternating harmonic series)
 * dirichlet_eta(new Integer(1));
 *
 * // η(2) = π²/12
 * dirichlet_eta(new Integer(2));
 *
 * // Symbolic usage
 * const s = new Symbol('s');
 * dirichlet_eta(s);
 * ```
 *
 * @see {@link zeta} - Riemann zeta function
 */
export function dirichlet_eta(s: Expr): Expr {
  const wasm = getWasmModule();
  const obj = createBasic();
  try {
    const code = wasm._basic_dirichlet_eta(obj.getPtr(), s.getWasmPtr());
    checkException(code);
    return exprFromWasm(obj);
  } catch (e) {
    obj.free();
    throw e;
  }
}

/**
 * Computes the Kronecker delta of two symbolic expressions.
 *
 * The Kronecker delta δ(i, j) equals 1 if i = j and 0 otherwise.
 * It is the discrete analog of the Dirac delta function and is
 * fundamental in linear algebra and tensor analysis.
 *
 * @param i - The first index (as a symbolic expression).
 * @param j - The second index (as a symbolic expression).
 * @returns δ(i, j) = 1 if i = j, else 0.
 *
 * @remarks
 * - δ(i, j) = δ(j, i) (symmetric)
 * - δ(i, i) = 1
 * - Σⱼ δ(i, j) × a_j = a_i (index substitution)
 * - δ(i, j) × δ(j, k) = δ(i, k) × δ(j, k)
 * - The identity matrix I_{ij} = δ(i, j)
 * - Used in Einstein summation notation
 *
 * @example
 * ```typescript
 * import { Symbol, Integer, kronecker_delta } from 'symwasm';
 *
 * const i = new Symbol('i');
 * const j = new Symbol('j');
 *
 * // Symbolic Kronecker delta
 * kronecker_delta(i, j);  // δ(i, j)
 *
 * // Numeric values
 * kronecker_delta(new Integer(1), new Integer(1));  // 1
 * kronecker_delta(new Integer(1), new Integer(2));  // 0
 * ```
 *
 * @see {@link dirichlet_eta} - Dirichlet eta function
 */
export function kronecker_delta(i: Expr, j: Expr): Expr {
  const wasm = getWasmModule();
  const obj = createBasic();
  try {
    const code = wasm._basic_kronecker_delta(obj.getPtr(), i.getWasmPtr(), j.getWasmPtr());
    checkException(code);
    return exprFromWasm(obj);
  } catch (e) {
    obj.free();
    throw e;
  }
}
