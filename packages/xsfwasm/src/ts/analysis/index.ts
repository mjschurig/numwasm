/**
 * Analysis functions
 *
 * Includes Riemann zeta function and Lambert W function.
 */

import { ensureLoaded, toFloat64Array, type ArrayInput } from '../core/utils.js';

// =============================================================================
// RIEMANN ZETA FUNCTION
// =============================================================================

/**
 * Riemann zeta function ζ(x).
 *
 * ζ(x) = Σₙ₌₁^∞ 1/n^x  for x > 1
 *
 * Extended to other values by analytic continuation.
 *
 * @param x - Input value
 * @returns ζ(x)
 *
 * @example
 * ```ts
 * zeta(2);    // π²/6 ≈ 1.6449
 * zeta(4);    // π⁴/90 ≈ 1.0823
 * zeta(-1);   // -1/12 (by analytic continuation)
 * ```
 */
export function zeta(x: number): number;
export function zeta(x: ArrayInput): Float64Array;
export function zeta(x: number | ArrayInput): number | Float64Array {
  const xsf = ensureLoaded();

  if (typeof x === 'number') {
    return xsf._wasm_zeta(x);
  }

  const arr = toFloat64Array(x);
  const result = new Float64Array(arr.length);
  for (let i = 0; i < arr.length; i++) {
    result[i] = xsf._wasm_zeta(arr[i]);
  }
  return result;
}

/**
 * Riemann zeta complement ζ(x) - 1.
 *
 * More accurate than zeta(x) - 1 for large x where ζ(x) ≈ 1.
 *
 * @param x - Input value
 * @returns ζ(x) - 1
 *
 * @example
 * ```ts
 * zetac(2);    // π²/6 - 1 ≈ 0.6449
 * zetac(10);   // ≈ 0.0009945
 * ```
 */
export function zetac(x: number): number;
export function zetac(x: ArrayInput): Float64Array;
export function zetac(x: number | ArrayInput): number | Float64Array {
  const xsf = ensureLoaded();

  if (typeof x === 'number') {
    return xsf._wasm_zetac(x);
  }

  const arr = toFloat64Array(x);
  const result = new Float64Array(arr.length);
  for (let i = 0; i < arr.length; i++) {
    result[i] = xsf._wasm_zetac(arr[i]);
  }
  return result;
}

// =============================================================================
// LAMBERT W FUNCTION
// =============================================================================

/**
 * Lambert W function W(x, k).
 *
 * The Lambert W function is the inverse of f(w) = w·e^w.
 * It satisfies: W(x)·e^(W(x)) = x
 *
 * @param x - Input value
 * @param k - Branch (0 for principal branch, -1 for secondary real branch)
 * @returns W(x) on branch k
 *
 * @example
 * ```ts
 * lambertw(1);       // 0.5671... (W(1) = Ω, the Omega constant)
 * lambertw(0);       // 0
 * lambertw(-0.3, -1); // -1.781... (secondary branch)
 * ```
 */
export function lambertw(x: number, k?: number): number;
export function lambertw(x: ArrayInput, k?: number): Float64Array;
export function lambertw(
  x: number | ArrayInput,
  k: number = 0
): number | Float64Array {
  const xsf = ensureLoaded();

  if (typeof x === 'number') {
    return xsf._wasm_lambertw(x, k);
  }

  const arr = toFloat64Array(x);
  const result = new Float64Array(arr.length);
  for (let i = 0; i < arr.length; i++) {
    result[i] = xsf._wasm_lambertw(arr[i], k);
  }
  return result;
}
