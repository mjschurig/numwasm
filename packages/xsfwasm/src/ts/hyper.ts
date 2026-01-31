/**
 * Hypergeometric function - high-level API
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
 * Gauss hypergeometric function ₂F₁(a, b; c; x).
 *
 * The hypergeometric function is defined by the series:
 *   ₂F₁(a, b; c; x) = Σₙ (a)ₙ(b)ₙ / ((c)ₙ n!) xⁿ
 *
 * Many special functions are special cases of ₂F₁.
 *
 * @param a - First numerator parameter
 * @param b - Second numerator parameter
 * @param c - Denominator parameter (c ≠ 0, -1, -2, ...)
 * @param x - Argument (|x| < 1 for convergence, can be extended)
 * @returns ₂F₁(a, b; c; x)
 *
 * @example
 * ```ts
 * hyp2f1(1, 1, 2, 0.5);  // ln(2) ≈ 0.693
 * hyp2f1(0.5, 0.5, 1.5, 0.5);  // Elliptic K related
 * ```
 */
export function hyp2f1(a: number, b: number, c: number, x: number): number;
export function hyp2f1(
  a: number,
  b: number,
  c: number,
  x: ArrayInput
): Float64Array;
export function hyp2f1(
  a: number,
  b: number,
  c: number,
  x: number | ArrayInput
): number | Float64Array {
  const xsf = ensureLoaded();

  if (typeof x === 'number') {
    return xsf._wasm_hyp2f1(a, b, c, x);
  }

  const arr = x instanceof Float64Array ? x : new Float64Array(x);
  const result = new Float64Array(arr.length);
  for (let i = 0; i < arr.length; i++) {
    result[i] = xsf._wasm_hyp2f1(a, b, c, arr[i]);
  }
  return result;
}
