/**
 * Spherical Bessel functions - high-level API
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
 * Spherical Bessel function of the first kind j_n(x).
 *
 * j_n(x) = √(π/(2x)) J_{n+1/2}(x)
 *
 * @param n - Order (non-negative integer)
 * @param x - Argument
 * @returns j_n(x)
 *
 * @example
 * ```ts
 * spherical_jn(0, 1);  // sin(1)/1 ≈ 0.841
 * spherical_jn(1, 1);  // sin(1) - cos(1) ≈ 0.301
 * ```
 */
export function spherical_jn(n: number, x: number): number;
export function spherical_jn(n: number, x: ArrayInput): Float64Array;
export function spherical_jn(
  n: number,
  x: number | ArrayInput
): number | Float64Array {
  const xsf = ensureLoaded();

  if (typeof x === 'number') {
    return xsf._wasm_spherical_jn(n, x);
  }

  const arr = x instanceof Float64Array ? x : new Float64Array(x);
  const result = new Float64Array(arr.length);
  for (let i = 0; i < arr.length; i++) {
    result[i] = xsf._wasm_spherical_jn(n, arr[i]);
  }
  return result;
}

/**
 * Spherical Bessel function of the second kind y_n(x).
 *
 * y_n(x) = √(π/(2x)) Y_{n+1/2}(x)
 *
 * @param n - Order (non-negative integer)
 * @param x - Argument (x > 0)
 * @returns y_n(x)
 *
 * @example
 * ```ts
 * spherical_yn(0, 1);  // -cos(1)/1 ≈ -0.540
 * ```
 */
export function spherical_yn(n: number, x: number): number;
export function spherical_yn(n: number, x: ArrayInput): Float64Array;
export function spherical_yn(
  n: number,
  x: number | ArrayInput
): number | Float64Array {
  const xsf = ensureLoaded();

  if (typeof x === 'number') {
    return xsf._wasm_spherical_yn(n, x);
  }

  const arr = x instanceof Float64Array ? x : new Float64Array(x);
  const result = new Float64Array(arr.length);
  for (let i = 0; i < arr.length; i++) {
    result[i] = xsf._wasm_spherical_yn(n, arr[i]);
  }
  return result;
}

/**
 * Modified spherical Bessel function of the first kind i_n(x).
 *
 * i_n(x) = √(π/(2x)) I_{n+1/2}(x)
 *
 * @param n - Order (non-negative integer)
 * @param x - Argument
 * @returns i_n(x)
 *
 * @example
 * ```ts
 * spherical_in(0, 1);  // sinh(1)/1 ≈ 1.175
 * ```
 */
export function spherical_in(n: number, x: number): number;
export function spherical_in(n: number, x: ArrayInput): Float64Array;
export function spherical_in(
  n: number,
  x: number | ArrayInput
): number | Float64Array {
  const xsf = ensureLoaded();

  if (typeof x === 'number') {
    return xsf._wasm_spherical_in(n, x);
  }

  const arr = x instanceof Float64Array ? x : new Float64Array(x);
  const result = new Float64Array(arr.length);
  for (let i = 0; i < arr.length; i++) {
    result[i] = xsf._wasm_spherical_in(n, arr[i]);
  }
  return result;
}

/**
 * Modified spherical Bessel function of the second kind k_n(x).
 *
 * k_n(x) = √(π/(2x)) K_{n+1/2}(x)
 *
 * @param n - Order (non-negative integer)
 * @param x - Argument (x > 0)
 * @returns k_n(x)
 *
 * @example
 * ```ts
 * spherical_kn(0, 1);  // π/(2e) ≈ 0.578
 * ```
 */
export function spherical_kn(n: number, x: number): number;
export function spherical_kn(n: number, x: ArrayInput): Float64Array;
export function spherical_kn(
  n: number,
  x: number | ArrayInput
): number | Float64Array {
  const xsf = ensureLoaded();

  if (typeof x === 'number') {
    return xsf._wasm_spherical_kn(n, x);
  }

  const arr = x instanceof Float64Array ? x : new Float64Array(x);
  const result = new Float64Array(arr.length);
  for (let i = 0; i < arr.length; i++) {
    result[i] = xsf._wasm_spherical_kn(n, arr[i]);
  }
  return result;
}
