/**
 * Gamma function family - high-level API
 *
 * Provides user-friendly wrappers for gamma-related special functions
 * with support for both scalar and array inputs.
 */

import { getXSFModule, isXSFLoaded } from './loader.js';
import type { XSFModule } from './types.js';

type ArrayInput = number[] | Float64Array;

function ensureLoaded(): XSFModule {
  if (!isXSFLoaded()) {
    throw new Error(
      'XSF WASM module not loaded. Call "await loadXSFModule()" before using special functions.'
    );
  }
  return getXSFModule();
}

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
 * import { loadXSFModule, gamma } from 'xsfwasm';
 *
 * await loadXSFModule();
 *
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

  const arr = x instanceof Float64Array ? x : new Float64Array(x);
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

  const arr = x instanceof Float64Array ? x : new Float64Array(x);
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

  const arr = x instanceof Float64Array ? x : new Float64Array(x);
  const result = new Float64Array(arr.length);
  for (let i = 0; i < arr.length; i++) {
    result[i] = xsf._wasm_rgamma(arr[i]);
  }
  return result;
}
