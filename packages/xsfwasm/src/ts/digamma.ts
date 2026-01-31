/**
 * Digamma (psi) function - high-level API
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

  const arr = x instanceof Float64Array ? x : new Float64Array(x);
  const result = new Float64Array(arr.length);
  for (let i = 0; i < arr.length; i++) {
    result[i] = xsf._wasm_digamma(arr[i]);
  }
  return result;
}
