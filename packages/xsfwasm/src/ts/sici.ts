/**
 * Sine and cosine integrals - high-level API
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

export interface SiCiResult {
  si: number;
  ci: number;
}

export interface SiCiArrayResult {
  si: Float64Array;
  ci: Float64Array;
}

/**
 * Sine and cosine integrals Si(x) and Ci(x).
 *
 * Si(x) = ∫₀^x sin(t)/t dt
 * Ci(x) = γ + ln(x) + ∫₀^x (cos(t) - 1)/t dt
 *
 * @param x - Input value
 * @returns Object with si (sine integral) and ci (cosine integral)
 *
 * @example
 * ```ts
 * const { si, ci } = sici(1);  // si ≈ 0.946, ci ≈ 0.337
 * ```
 */
export function sici(x: number): SiCiResult;
export function sici(x: ArrayInput): SiCiArrayResult;
export function sici(x: number | ArrayInput): SiCiResult | SiCiArrayResult {
  const xsf = ensureLoaded();

  if (typeof x === 'number') {
    const siPtr = xsf._malloc(8);
    const ciPtr = xsf._malloc(8);
    try {
      xsf._wasm_sici(x, siPtr, ciPtr);
      return {
        si: xsf.HEAPF64[siPtr >> 3],
        ci: xsf.HEAPF64[ciPtr >> 3],
      };
    } finally {
      xsf._free(siPtr);
      xsf._free(ciPtr);
    }
  }

  const arr = x instanceof Float64Array ? x : new Float64Array(x);
  const si = new Float64Array(arr.length);
  const ci = new Float64Array(arr.length);

  const siPtr = xsf._malloc(8);
  const ciPtr = xsf._malloc(8);
  try {
    for (let i = 0; i < arr.length; i++) {
      xsf._wasm_sici(arr[i], siPtr, ciPtr);
      si[i] = xsf.HEAPF64[siPtr >> 3];
      ci[i] = xsf.HEAPF64[ciPtr >> 3];
    }
  } finally {
    xsf._free(siPtr);
    xsf._free(ciPtr);
  }

  return { si, ci };
}

export interface ShiChiResult {
  shi: number;
  chi: number;
}

export interface ShiChiArrayResult {
  shi: Float64Array;
  chi: Float64Array;
}

/**
 * Hyperbolic sine and cosine integrals Shi(x) and Chi(x).
 *
 * Shi(x) = ∫₀^x sinh(t)/t dt
 * Chi(x) = γ + ln(x) + ∫₀^x (cosh(t) - 1)/t dt
 *
 * @param x - Input value
 * @returns Object with shi and chi
 *
 * @example
 * ```ts
 * const { shi, chi } = shichi(1);  // shi ≈ 1.057, chi ≈ 0.837
 * ```
 */
export function shichi(x: number): ShiChiResult;
export function shichi(x: ArrayInput): ShiChiArrayResult;
export function shichi(x: number | ArrayInput): ShiChiResult | ShiChiArrayResult {
  const xsf = ensureLoaded();

  if (typeof x === 'number') {
    const shiPtr = xsf._malloc(8);
    const chiPtr = xsf._malloc(8);
    try {
      xsf._wasm_shichi(x, shiPtr, chiPtr);
      return {
        shi: xsf.HEAPF64[shiPtr >> 3],
        chi: xsf.HEAPF64[chiPtr >> 3],
      };
    } finally {
      xsf._free(shiPtr);
      xsf._free(chiPtr);
    }
  }

  const arr = x instanceof Float64Array ? x : new Float64Array(x);
  const shi = new Float64Array(arr.length);
  const chi = new Float64Array(arr.length);

  const shiPtr = xsf._malloc(8);
  const chiPtr = xsf._malloc(8);
  try {
    for (let i = 0; i < arr.length; i++) {
      xsf._wasm_shichi(arr[i], shiPtr, chiPtr);
      shi[i] = xsf.HEAPF64[shiPtr >> 3];
      chi[i] = xsf.HEAPF64[chiPtr >> 3];
    }
  } finally {
    xsf._free(shiPtr);
    xsf._free(chiPtr);
  }

  return { shi, chi };
}
