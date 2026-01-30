/**
 * Beta function and related special functions via WebAssembly.
 * Mirrors scipy.special.beta, betaln functions.
 */

import { loadXSFModule } from 'xsfwasm';

/**
 * Beta function: B(a, b) = Γ(a)Γ(b) / Γ(a+b)
 *
 * The beta function is related to the gamma function by:
 *   B(a, b) = Γ(a) * Γ(b) / Γ(a + b)
 *
 * @param a - First parameter. Can be a number or array of numbers.
 * @param b - Second parameter. Can be a number or array of numbers.
 * @returns Beta function value. Returns number for scalar inputs, array for array inputs.
 *
 * @example
 * ```ts
 * import { special } from 'sciwasm';
 *
 * // B(2, 3) = 1/12 ≈ 0.0833
 * const result = await special.beta(2, 3);
 * console.log(result); // 0.08333...
 *
 * // Array input
 * const results = await special.beta([1, 2, 3], [1, 2, 3]);
 * console.log(results); // [1, 0.166..., 0.0333...]
 * ```
 */
export async function beta(
  a: number | number[],
  b: number | number[]
): Promise<number | number[]> {
  const Module = await loadXSFModule();

  // Both arrays
  if (Array.isArray(a) && Array.isArray(b)) {
    if (a.length !== b.length) {
      throw new Error('Arrays must have the same length');
    }
    return a.map((ai, i) => Module._wasm_beta(ai, b[i]));
  }

  // a is array, b is scalar
  if (Array.isArray(a)) {
    return a.map((ai) => Module._wasm_beta(ai, b as number));
  }

  // b is array, a is scalar
  if (Array.isArray(b)) {
    return b.map((bi) => Module._wasm_beta(a as number, bi));
  }

  // Both scalars
  return Module._wasm_beta(a, b);
}

/**
 * Natural logarithm of the beta function: ln(B(a, b))
 *
 * For large arguments, betaln(a, b) is more accurate than log(beta(a, b))
 * because it avoids overflow of the beta function.
 *
 * @param a - First parameter. Can be a number or array of numbers.
 * @param b - Second parameter. Can be a number or array of numbers.
 * @returns Natural log of beta function. Returns number for scalar inputs, array for array inputs.
 *
 * @example
 * ```ts
 * import { special } from 'sciwasm';
 *
 * // For large values, beta would overflow but betaln is stable
 * const result = await special.betaln(100, 200);
 * console.log(result); // -205.026...
 *
 * // Compare with log(beta) for small values
 * const a = 2, b = 3;
 * const method1 = await special.betaln(a, b);
 * const method2 = Math.log(await special.beta(a, b));
 * // Both ≈ -2.485
 * ```
 */
export async function betaln(
  a: number | number[],
  b: number | number[]
): Promise<number | number[]> {
  const Module = await loadXSFModule();

  // Both arrays
  if (Array.isArray(a) && Array.isArray(b)) {
    if (a.length !== b.length) {
      throw new Error('Arrays must have the same length');
    }
    return a.map((ai, i) => Module._wasm_betaln(ai, b[i]));
  }

  // a is array, b is scalar
  if (Array.isArray(a)) {
    return a.map((ai) => Module._wasm_betaln(ai, b as number));
  }

  // b is array, a is scalar
  if (Array.isArray(b)) {
    return b.map((bi) => Module._wasm_betaln(a as number, bi));
  }

  // Both scalars
  return Module._wasm_betaln(a, b);
}
