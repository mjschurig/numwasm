/**
 * Gamma function family
 *
 * Includes gamma, log-gamma, reciprocal gamma, beta, log-beta, and digamma functions.
 */

import { ensureLoaded, toFloat64Array, type ArrayInput } from '../core/utils.js';

// =============================================================================
// GAMMA FUNCTIONS
// =============================================================================

/**
 * Gamma function Γ(x).
 *
 * The gamma function extends the factorial to real and complex numbers:
 * - For positive integers: gamma(n) = (n-1)!
 * - gamma(1) = 1
 * - gamma(0.5) = √π ≈ 1.7724538509
 * - gamma(x+1) = x * gamma(x)
 *
 * @param x - Input value (must not be a non-positive integer)
 * @returns Γ(x)
 *
 * @example
 * ```ts
 * gamma(5);           // 24 (= 4!)
 * gamma(0.5);         // 1.7724... (= √π)
 * gamma([1, 2, 3]);   // Float64Array([1, 1, 2])
 * ```
 */
export function gamma(x: number): number;
export function gamma(x: ArrayInput): Float64Array;
export function gamma(x: number | ArrayInput): number | Float64Array {
  const xsf = ensureLoaded();

  if (typeof x === 'number') {
    return xsf._wasm_gamma(x);
  }

  const arr = toFloat64Array(x);
  const result = new Float64Array(arr.length);
  for (let i = 0; i < arr.length; i++) {
    result[i] = xsf._wasm_gamma(arr[i]);
  }
  return result;
}

/**
 * Natural logarithm of the absolute value of the gamma function.
 *
 * gammaln(x) = ln(|Γ(x)|)
 *
 * More numerically stable than log(gamma(x)) for large x.
 * Essential for computing with large factorials in statistical distributions.
 *
 * @param x - Input value (must not be a non-positive integer)
 * @returns ln(|Γ(x)|)
 *
 * @example
 * ```ts
 * gammaln(100);   // 359.13... (log of 99!)
 * gammaln(1000);  // 5905.22... (would overflow with gamma())
 * ```
 */
export function gammaln(x: number): number;
export function gammaln(x: ArrayInput): Float64Array;
export function gammaln(x: number | ArrayInput): number | Float64Array {
  const xsf = ensureLoaded();

  if (typeof x === 'number') {
    return xsf._wasm_gammaln(x);
  }

  const arr = toFloat64Array(x);
  const result = new Float64Array(arr.length);
  for (let i = 0; i < arr.length; i++) {
    result[i] = xsf._wasm_gammaln(arr[i]);
  }
  return result;
}

/**
 * Reciprocal of the gamma function: 1/Γ(x).
 *
 * Unlike gamma(x), rgamma(x) is an entire function (analytic everywhere).
 * It equals zero at non-positive integers where gamma() has poles.
 *
 * @param x - Input value
 * @returns 1/Γ(x)
 *
 * @example
 * ```ts
 * rgamma(5);    // 0.04166... (= 1/24)
 * rgamma(0);    // 0 (gamma(0) is infinite)
 * rgamma(-1);   // 0 (gamma(-1) is infinite)
 * ```
 */
export function rgamma(x: number): number;
export function rgamma(x: ArrayInput): Float64Array;
export function rgamma(x: number | ArrayInput): number | Float64Array {
  const xsf = ensureLoaded();

  if (typeof x === 'number') {
    return xsf._wasm_rgamma(x);
  }

  const arr = toFloat64Array(x);
  const result = new Float64Array(arr.length);
  for (let i = 0; i < arr.length; i++) {
    result[i] = xsf._wasm_rgamma(arr[i]);
  }
  return result;
}

// =============================================================================
// BETA FUNCTIONS
// =============================================================================

/**
 * Beta function B(a, b).
 *
 * The beta function is defined as:
 *   B(a, b) = ∫₀¹ t^(a-1) * (1-t)^(b-1) dt = Γ(a) * Γ(b) / Γ(a+b)
 *
 * Properties:
 * - B(a, b) = B(b, a) (symmetric)
 * - B(1, 1) = 1
 * - B(a, 1) = 1/a
 *
 * @param a - First parameter (must be positive)
 * @param b - Second parameter (must be positive)
 * @returns B(a, b)
 *
 * @example
 * ```ts
 * beta(2, 3);      // 0.08333... (= 1/12)
 * beta(0.5, 0.5);  // π ≈ 3.14159...
 *
 * // Vectorized with broadcasting
 * beta([1, 2, 3], 2);  // Float64Array([0.5, 0.1666..., 0.0833...])
 * ```
 */
export function beta(a: number, b: number): number;
export function beta(a: ArrayInput, b: number): Float64Array;
export function beta(a: number, b: ArrayInput): Float64Array;
export function beta(a: ArrayInput, b: ArrayInput): Float64Array;
export function beta(
  a: number | ArrayInput,
  b: number | ArrayInput
): number | Float64Array {
  const xsf = ensureLoaded();

  const aIsScalar = typeof a === 'number';
  const bIsScalar = typeof b === 'number';

  if (aIsScalar && bIsScalar) {
    return xsf._wasm_beta(a as number, b as number);
  }

  const aArr = aIsScalar ? null : toFloat64Array(a as ArrayInput);
  const bArr = bIsScalar ? null : toFloat64Array(b as ArrayInput);

  const length = aArr?.length ?? bArr!.length;
  if (aArr && bArr && aArr.length !== bArr.length) {
    throw new Error(
      `Array length mismatch: got arrays of length ${aArr.length} and ${bArr.length}`
    );
  }

  const result = new Float64Array(length);
  for (let i = 0; i < length; i++) {
    const aVal = aIsScalar ? (a as number) : aArr![i];
    const bVal = bIsScalar ? (b as number) : bArr![i];
    result[i] = xsf._wasm_beta(aVal, bVal);
  }
  return result;
}

/**
 * Natural logarithm of the beta function.
 *
 * betaln(a, b) = ln(B(a, b)) = ln(Γ(a)) + ln(Γ(b)) - ln(Γ(a+b))
 *
 * More numerically stable than log(beta(a, b)) for large arguments.
 * Essential for computing binomial coefficients with large n.
 *
 * @param a - First parameter (must be positive)
 * @param b - Second parameter (must be positive)
 * @returns ln(B(a, b))
 *
 * @example
 * ```ts
 * betaln(100, 200);  // -202.58... (beta would underflow)
 * betaln(2, 3);      // -2.4849... (= ln(1/12))
 * ```
 */
export function betaln(a: number, b: number): number;
export function betaln(a: ArrayInput, b: number): Float64Array;
export function betaln(a: number, b: ArrayInput): Float64Array;
export function betaln(a: ArrayInput, b: ArrayInput): Float64Array;
export function betaln(
  a: number | ArrayInput,
  b: number | ArrayInput
): number | Float64Array {
  const xsf = ensureLoaded();

  const aIsScalar = typeof a === 'number';
  const bIsScalar = typeof b === 'number';

  if (aIsScalar && bIsScalar) {
    return xsf._wasm_betaln(a as number, b as number);
  }

  const aArr = aIsScalar ? null : toFloat64Array(a as ArrayInput);
  const bArr = bIsScalar ? null : toFloat64Array(b as ArrayInput);

  const length = aArr?.length ?? bArr!.length;
  if (aArr && bArr && aArr.length !== bArr.length) {
    throw new Error(
      `Array length mismatch: got arrays of length ${aArr.length} and ${bArr.length}`
    );
  }

  const result = new Float64Array(length);
  for (let i = 0; i < length; i++) {
    const aVal = aIsScalar ? (a as number) : aArr![i];
    const bVal = bIsScalar ? (b as number) : bArr![i];
    result[i] = xsf._wasm_betaln(aVal, bVal);
  }
  return result;
}

// =============================================================================
// DIGAMMA FUNCTION
// =============================================================================

/**
 * Computes the digamma function ψ(x) = d/dx ln(Γ(x)).
 *
 * The digamma function is the logarithmic derivative of the gamma function.
 *
 * Properties:
 * - ψ(1) = -γ (Euler-Mascheroni constant ≈ -0.5772)
 * - ψ(x+1) = ψ(x) + 1/x
 * - Has poles at non-positive integers
 *
 * @param x - Input value (must not be a non-positive integer)
 * @returns ψ(x)
 *
 * @example
 * ```ts
 * digamma(1);    // -0.5772... (Euler-Mascheroni constant)
 * digamma(2);    // 0.4227...
 * ```
 */
export function digamma(x: number): number;
export function digamma(x: ArrayInput): Float64Array;
export function digamma(x: number | ArrayInput): number | Float64Array {
  const xsf = ensureLoaded();

  if (typeof x === 'number') {
    return xsf._wasm_digamma(x);
  }

  const arr = toFloat64Array(x);
  const result = new Float64Array(arr.length);
  for (let i = 0; i < arr.length; i++) {
    result[i] = xsf._wasm_digamma(arr[i]);
  }
  return result;
}
