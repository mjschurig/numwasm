/**
 * Fresnel integrals - high-level API
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

export interface FresnelResult {
  s: number;
  c: number;
}

export interface FresnelArrayResult {
  s: Float64Array;
  c: Float64Array;
}

/**
 * Fresnel integrals S(x) and C(x).
 *
 * S(x) = ∫₀^x sin(πt²/2) dt
 * C(x) = ∫₀^x cos(πt²/2) dt
 *
 * @param x - Input value
 * @returns Object with s (sine integral) and c (cosine integral)
 *
 * @example
 * ```ts
 * const { s, c } = fresnel(1);  // s ≈ 0.438, c ≈ 0.780
 * ```
 */
export function fresnel(x: number): FresnelResult;
export function fresnel(x: ArrayInput): FresnelArrayResult;
export function fresnel(x: number | ArrayInput): FresnelResult | FresnelArrayResult {
  const xsf = ensureLoaded();

  if (typeof x === 'number') {
    const sPtr = xsf._malloc(8);
    const cPtr = xsf._malloc(8);
    try {
      xsf._wasm_fresnel(x, sPtr, cPtr);
      return {
        s: xsf.HEAPF64[sPtr >> 3],
        c: xsf.HEAPF64[cPtr >> 3],
      };
    } finally {
      xsf._free(sPtr);
      xsf._free(cPtr);
    }
  }

  const arr = x instanceof Float64Array ? x : new Float64Array(x);
  const s = new Float64Array(arr.length);
  const c = new Float64Array(arr.length);

  const sPtr = xsf._malloc(8);
  const cPtr = xsf._malloc(8);
  try {
    for (let i = 0; i < arr.length; i++) {
      xsf._wasm_fresnel(arr[i], sPtr, cPtr);
      s[i] = xsf.HEAPF64[sPtr >> 3];
      c[i] = xsf.HEAPF64[cPtr >> 3];
    }
  } finally {
    xsf._free(sPtr);
    xsf._free(cPtr);
  }

  return { s, c };
}
