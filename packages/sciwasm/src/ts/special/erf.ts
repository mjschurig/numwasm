/**
 * Error functions via WebAssembly.
 * Mirrors scipy.special.erf, erfc, erfcx, erfi functions.
 */

import { loadXSFModule } from 'xsfwasm';

/**
 * Error function: erf(x) = 2/√π ∫₀ˣ e^(-t²) dt
 *
 * The error function is the integral of the Gaussian distribution
 * and is widely used in probability, statistics, and partial differential equations.
 *
 * Properties:
 *   - erf(0) = 0
 *   - erf(∞) = 1
 *   - erf(-x) = -erf(x)
 *
 * @param x - Input value(s). Can be a number or array of numbers.
 * @returns Error function value(s).
 *
 * @example
 * ```ts
 * import { special } from 'sciwasm';
 *
 * // erf(1) ≈ 0.8427
 * const result = await special.erf(1);
 * console.log(result);
 *
 * // Array input
 * const results = await special.erf([0, 0.5, 1, 2]);
 * console.log(results); // [0, 0.52..., 0.84..., 0.99...]
 * ```
 */
export async function erf(x: number | number[]): Promise<number | number[]> {
  const Module = await loadXSFModule();

  if (Array.isArray(x)) {
    return x.map((val) => Module._wasm_erf(val));
  }

  return Module._wasm_erf(x);
}

/**
 * Complementary error function: erfc(x) = 1 - erf(x)
 *
 * For large x, erfc(x) is more accurate than computing 1 - erf(x)
 * directly because erf(x) → 1 causes precision loss.
 *
 * @param x - Input value(s). Can be a number or array of numbers.
 * @returns Complementary error function value(s).
 *
 * @example
 * ```ts
 * import { special } from 'sciwasm';
 *
 * // erfc(1) = 1 - erf(1) ≈ 0.1573
 * const result = await special.erfc(1);
 * console.log(result);
 *
 * // For large x, erfc is more accurate
 * const largeX = await special.erfc(5);
 * console.log(largeX); // ~1.54e-12
 * ```
 */
export async function erfc(x: number | number[]): Promise<number | number[]> {
  const Module = await loadXSFModule();

  if (Array.isArray(x)) {
    return x.map((val) => Module._wasm_erfc(val));
  }

  return Module._wasm_erfc(x);
}

/**
 * Scaled complementary error function: erfcx(x) = e^(x²) * erfc(x)
 *
 * This function avoids the underflow/overflow issues that occur
 * when computing erfc(x) for large |x|.
 *
 * Properties:
 *   - erfcx(0) = 1
 *   - erfcx(x) → 1/(√π * x) as x → ∞
 *
 * @param x - Input value(s). Can be a number or array of numbers.
 * @returns Scaled complementary error function value(s).
 *
 * @example
 * ```ts
 * import { special } from 'sciwasm';
 *
 * // erfcx(0) = 1
 * const result = await special.erfcx(0);
 * console.log(result); // 1
 *
 * // For large x, erfcx is well-behaved
 * const largeX = await special.erfcx(100);
 * console.log(largeX); // ~0.00564 (would overflow with exp(10000) * erfc(100))
 * ```
 */
export async function erfcx(x: number | number[]): Promise<number | number[]> {
  const Module = await loadXSFModule();

  if (Array.isArray(x)) {
    return x.map((val) => Module._wasm_erfcx(val));
  }

  return Module._wasm_erfcx(x);
}

/**
 * Imaginary error function: erfi(x) = -i * erf(ix)
 *
 * The imaginary error function is the error function evaluated
 * at an imaginary argument, made real.
 *
 * @param x - Input value(s). Can be a number or array of numbers.
 * @returns Imaginary error function value(s).
 *
 * @example
 * ```ts
 * import { special } from 'sciwasm';
 *
 * // erfi(1) ≈ 1.6505
 * const result = await special.erfi(1);
 * console.log(result);
 *
 * // erfi(0) = 0
 * const zero = await special.erfi(0);
 * console.log(zero); // 0
 * ```
 */
export async function erfi(x: number | number[]): Promise<number | number[]> {
  const Module = await loadXSFModule();

  if (Array.isArray(x)) {
    return x.map((val) => Module._wasm_erfi(val));
  }

  return Module._wasm_erfi(x);
}
