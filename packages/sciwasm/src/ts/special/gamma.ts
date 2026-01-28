/**
 * Gamma function and related special functions via WebAssembly.
 * Mirrors scipy.special.gamma, gammaln functions.
 */

import { loadWasmModule } from '../wasm-loader.js';

/**
 * Gamma function: Γ(x)
 *
 * Returns the gamma function of the argument. The gamma function is defined as:
 *
 *     Γ(n) = (n-1)!  for positive integers n
 *     Γ(x) = ∫₀^∞ t^(x-1) e^(-t) dt  for x > 0
 *
 * The gamma function is extended to the entire complex plane except for
 * negative integers and zero.
 *
 * @param x - Input value(s). Can be a number or array of numbers.
 * @returns Gamma function of x. Returns number for scalar input, array for array input.
 *
 * @example
 * ```ts
 * import { special } from 'sciwasm';
 *
 * // Gamma(5) = 4! = 24
 * const result = await special.gamma(5);
 * console.log(result); // 24
 *
 * // Array input
 * const results = await special.gamma([1, 2, 3, 4, 5]);
 * console.log(results); // [1, 1, 2, 6, 24]
 *
 * // Gamma of half-integer
 * const halfInt = await special.gamma(0.5);
 * console.log(halfInt); // ~1.772 (√π)
 * ```
 */
export async function gamma(x: number | number[]): Promise<number | number[]> {
  const Module = await loadWasmModule();

  // Array input
  if (Array.isArray(x)) {
    return x.map(val => Module._wasm_gamma(val));
  }

  // Scalar input
  return Module._wasm_gamma(x);
}

/**
 * Natural logarithm of the absolute value of the gamma function.
 *
 * For large arguments, gammaln(x) is more accurate than log(gamma(x))
 * because it avoids overflow of the gamma function.
 *
 * Returns: ln(|Γ(x)|)
 *
 * @param x - Input value(s). Can be a number or array of numbers.
 * @returns Natural log of absolute value of gamma function.
 *
 * @example
 * ```ts
 * import { special } from 'sciwasm';
 *
 * // For large x, gamma(x) would overflow
 * const largeX = 100;
 * const logGamma = await special.gammaln(largeX);
 * console.log(logGamma); // ~359.13 (ln(99!))
 *
 * // Compare with log(gamma(x)) for small x
 * const x = 5;
 * const method1 = await special.gammaln(x);
 * const method2 = Math.log(await special.gamma(x));
 * console.log(method1, method2); // Both ~3.178
 * ```
 */
export async function gammaln(x: number | number[]): Promise<number | number[]> {
  const Module = await loadWasmModule();

  // Array input
  if (Array.isArray(x)) {
    return x.map(val => Module._wasm_gammaln(val));
  }

  // Scalar input
  return Module._wasm_gammaln(x);
}

/**
 * Reciprocal of the gamma function: 1/Γ(x)
 *
 * This function is more numerically stable than computing 1/gamma(x),
 * especially for large |x| where gamma(x) might overflow or underflow.
 *
 * @param x - Input value(s). Can be a number or array of numbers.
 * @returns Reciprocal of the gamma function.
 *
 * @example
 * ```ts
 * import { special } from 'sciwasm';
 *
 * // At negative integers, gamma has poles but rgamma is well-defined (returns 0)
 * const atPole = await special.rgamma(-2);
 * console.log(atPole); // 0
 *
 * // For positive values
 * const x = 3;
 * const reciprocal = await special.rgamma(x);
 * console.log(reciprocal); // 0.5 (since Γ(3) = 2)
 * ```
 */
export async function rgamma(x: number | number[]): Promise<number | number[]> {
  const Module = await loadWasmModule();

  // Array input
  if (Array.isArray(x)) {
    return x.map(val => Module._wasm_rgamma(val));
  }

  // Scalar input
  return Module._wasm_rgamma(x);
}
