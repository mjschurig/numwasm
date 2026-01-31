/**
 * Legendre polynomials - high-level API
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
 * Legendre polynomial P_n(x).
 *
 * The Legendre polynomials are orthogonal polynomials on [-1, 1].
 *
 * P_0(x) = 1
 * P_1(x) = x
 * P_2(x) = (3x² - 1)/2
 * ...
 *
 * @param n - Degree (non-negative integer)
 * @param x - Argument (-1 ≤ x ≤ 1)
 * @returns P_n(x)
 *
 * @example
 * ```ts
 * legendre_p(0, 0.5);  // 1
 * legendre_p(1, 0.5);  // 0.5
 * legendre_p(2, 0.5);  // -0.125
 * ```
 */
export function legendre_p(n: number, x: number): number;
export function legendre_p(n: number, x: ArrayInput): Float64Array;
export function legendre_p(
  n: number,
  x: number | ArrayInput
): number | Float64Array {
  const xsf = ensureLoaded();

  if (typeof x === 'number') {
    return xsf._wasm_legendre_p(n, x);
  }

  const arr = x instanceof Float64Array ? x : new Float64Array(x);
  const result = new Float64Array(arr.length);
  for (let i = 0; i < arr.length; i++) {
    result[i] = xsf._wasm_legendre_p(n, arr[i]);
  }
  return result;
}
