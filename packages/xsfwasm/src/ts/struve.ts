/**
 * Struve functions - high-level API
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
 * Struve function H_v(x).
 *
 * The Struve function is a solution to the inhomogeneous Bessel equation.
 *
 * @param v - Order
 * @param x - Argument (x ≥ 0)
 * @returns H_v(x)
 *
 * @example
 * ```ts
 * struve_h(0, 1);   // ≈ 0.568
 * struve_h(1, 1);   // ≈ 0.198
 * ```
 */
export function struve_h(v: number, x: number): number;
export function struve_h(v: number, x: ArrayInput): Float64Array;
export function struve_h(v: ArrayInput, x: ArrayInput): Float64Array;
export function struve_h(
  v: number | ArrayInput,
  x: number | ArrayInput
): number | Float64Array {
  const xsf = ensureLoaded();

  const vIsScalar = typeof v === 'number';
  const xIsScalar = typeof x === 'number';

  if (vIsScalar && xIsScalar) {
    return xsf._wasm_struve_h(v as number, x as number);
  }

  const vArr = vIsScalar
    ? null
    : v instanceof Float64Array
      ? v
      : new Float64Array(v as number[]);
  const xArr = xIsScalar
    ? null
    : x instanceof Float64Array
      ? x
      : new Float64Array(x as number[]);

  const length = vArr?.length ?? xArr!.length;
  if (vArr && xArr && vArr.length !== xArr.length) {
    throw new Error(`Array length mismatch`);
  }

  const result = new Float64Array(length);
  for (let i = 0; i < length; i++) {
    const vVal = vIsScalar ? (v as number) : vArr![i];
    const xVal = xIsScalar ? (x as number) : xArr![i];
    result[i] = xsf._wasm_struve_h(vVal, xVal);
  }
  return result;
}

/**
 * Modified Struve function L_v(x).
 *
 * The modified Struve function is related to the Struve H function.
 *
 * @param v - Order
 * @param x - Argument (x ≥ 0)
 * @returns L_v(x)
 *
 * @example
 * ```ts
 * struve_l(0, 1);   // ≈ 0.710
 * struve_l(1, 1);   // ≈ 0.288
 * ```
 */
export function struve_l(v: number, x: number): number;
export function struve_l(v: number, x: ArrayInput): Float64Array;
export function struve_l(v: ArrayInput, x: ArrayInput): Float64Array;
export function struve_l(
  v: number | ArrayInput,
  x: number | ArrayInput
): number | Float64Array {
  const xsf = ensureLoaded();

  const vIsScalar = typeof v === 'number';
  const xIsScalar = typeof x === 'number';

  if (vIsScalar && xIsScalar) {
    return xsf._wasm_struve_l(v as number, x as number);
  }

  const vArr = vIsScalar
    ? null
    : v instanceof Float64Array
      ? v
      : new Float64Array(v as number[]);
  const xArr = xIsScalar
    ? null
    : x instanceof Float64Array
      ? x
      : new Float64Array(x as number[]);

  const length = vArr?.length ?? xArr!.length;
  if (vArr && xArr && vArr.length !== xArr.length) {
    throw new Error(`Array length mismatch`);
  }

  const result = new Float64Array(length);
  for (let i = 0; i < length; i++) {
    const vVal = vIsScalar ? (v as number) : vArr![i];
    const xVal = xIsScalar ? (x as number) : xArr![i];
    result[i] = xsf._wasm_struve_l(vVal, xVal);
  }
  return result;
}
