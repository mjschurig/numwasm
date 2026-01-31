/**
 * Exponential integrals - high-level API
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
 * Exponential integral E₁(x).
 *
 * E₁(x) = ∫ₓ^∞ e^(-t)/t dt  for x > 0
 *
 * @param x - Input value (x > 0)
 * @returns E₁(x)
 *
 * @example
 * ```ts
 * exp1(1);    // 0.2193...
 * exp1(0);    // Infinity
 * ```
 */
export function exp1(x: number): number;
export function exp1(x: ArrayInput): Float64Array;
export function exp1(x: number | ArrayInput): number | Float64Array {
  const xsf = ensureLoaded();

  if (typeof x === 'number') {
    return xsf._wasm_exp1(x);
  }

  const arr = x instanceof Float64Array ? x : new Float64Array(x);
  const result = new Float64Array(arr.length);
  for (let i = 0; i < arr.length; i++) {
    result[i] = xsf._wasm_exp1(arr[i]);
  }
  return result;
}

/**
 * Exponential integral Ei(x).
 *
 * Ei(x) = ∫₋∞^x e^t/t dt  (principal value)
 *
 * @param x - Input value
 * @returns Ei(x)
 *
 * @example
 * ```ts
 * expi(1);    // 1.8951...
 * ```
 */
export function expi(x: number): number;
export function expi(x: ArrayInput): Float64Array;
export function expi(x: number | ArrayInput): number | Float64Array {
  const xsf = ensureLoaded();

  if (typeof x === 'number') {
    return xsf._wasm_expi(x);
  }

  const arr = x instanceof Float64Array ? x : new Float64Array(x);
  const result = new Float64Array(arr.length);
  for (let i = 0; i < arr.length; i++) {
    result[i] = xsf._wasm_expi(arr[i]);
  }
  return result;
}
