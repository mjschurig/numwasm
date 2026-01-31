/**
 * Riemann zeta function - high-level API
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

  const arr = x instanceof Float64Array ? x : new Float64Array(x);
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

  const arr = x instanceof Float64Array ? x : new Float64Array(x);
  const result = new Float64Array(arr.length);
  for (let i = 0; i < arr.length; i++) {
    result[i] = xsf._wasm_zetac(arr[i]);
  }
  return result;
}
