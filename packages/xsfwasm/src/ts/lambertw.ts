/**
 * Lambert W function - high-level API
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

  const arr = x instanceof Float64Array ? x : new Float64Array(x);
  const result = new Float64Array(arr.length);
  for (let i = 0; i < arr.length; i++) {
    result[i] = xsf._wasm_lambertw(arr[i], k);
  }
  return result;
}
