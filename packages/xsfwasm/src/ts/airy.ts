/**
 * Airy functions - high-level API
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

export interface AiryResult {
  ai: number;
  aip: number;
  bi: number;
  bip: number;
}

export interface AiryArrayResult {
  ai: Float64Array;
  aip: Float64Array;
  bi: Float64Array;
  bip: Float64Array;
}

/**
 * Computes Airy functions Ai(x), Ai'(x), Bi(x), Bi'(x).
 *
 * The Airy functions are solutions to y'' - xy = 0.
 * - Ai(x): decays exponentially for x > 0
 * - Bi(x): grows exponentially for x > 0
 *
 * @param x - Input value
 * @returns Object with ai, aip (derivative), bi, bip (derivative)
 *
 * @example
 * ```ts
 * const { ai, bi } = airy(0);  // ai ≈ 0.355, bi ≈ 0.615
 * ```
 */
export function airy(x: number): AiryResult;
export function airy(x: ArrayInput): AiryArrayResult;
export function airy(x: number | ArrayInput): AiryResult | AiryArrayResult {
  const xsf = ensureLoaded();

  if (typeof x === 'number') {
    const aiPtr = xsf._malloc(8);
    const aipPtr = xsf._malloc(8);
    const biPtr = xsf._malloc(8);
    const bipPtr = xsf._malloc(8);
    try {
      xsf._wasm_airy(x, aiPtr, aipPtr, biPtr, bipPtr);
      return {
        ai: xsf.HEAPF64[aiPtr >> 3],
        aip: xsf.HEAPF64[aipPtr >> 3],
        bi: xsf.HEAPF64[biPtr >> 3],
        bip: xsf.HEAPF64[bipPtr >> 3],
      };
    } finally {
      xsf._free(aiPtr);
      xsf._free(aipPtr);
      xsf._free(biPtr);
      xsf._free(bipPtr);
    }
  }

  const arr = x instanceof Float64Array ? x : new Float64Array(x);
  const ai = new Float64Array(arr.length);
  const aip = new Float64Array(arr.length);
  const bi = new Float64Array(arr.length);
  const bip = new Float64Array(arr.length);

  const aiPtr = xsf._malloc(8);
  const aipPtr = xsf._malloc(8);
  const biPtr = xsf._malloc(8);
  const bipPtr = xsf._malloc(8);
  try {
    for (let i = 0; i < arr.length; i++) {
      xsf._wasm_airy(arr[i], aiPtr, aipPtr, biPtr, bipPtr);
      ai[i] = xsf.HEAPF64[aiPtr >> 3];
      aip[i] = xsf.HEAPF64[aipPtr >> 3];
      bi[i] = xsf.HEAPF64[biPtr >> 3];
      bip[i] = xsf.HEAPF64[bipPtr >> 3];
    }
  } finally {
    xsf._free(aiPtr);
    xsf._free(aipPtr);
    xsf._free(biPtr);
    xsf._free(bipPtr);
  }

  return { ai, aip, bi, bip };
}
