/**
 * Kelvin functions - high-level API
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
 * Kelvin function ber(x).
 *
 * ber(x) = Re[J₀(x·e^(3πi/4))]
 *
 * Related to Bessel functions of complex argument.
 *
 * @param x - Input value (x ≥ 0)
 * @returns ber(x)
 */
export function ber(x: number): number;
export function ber(x: ArrayInput): Float64Array;
export function ber(x: number | ArrayInput): number | Float64Array {
  const xsf = ensureLoaded();

  if (typeof x === 'number') {
    return xsf._wasm_ber(x);
  }

  const arr = x instanceof Float64Array ? x : new Float64Array(x);
  const result = new Float64Array(arr.length);
  for (let i = 0; i < arr.length; i++) {
    result[i] = xsf._wasm_ber(arr[i]);
  }
  return result;
}

/**
 * Kelvin function bei(x).
 *
 * bei(x) = Im[J₀(x·e^(3πi/4))]
 *
 * @param x - Input value (x ≥ 0)
 * @returns bei(x)
 */
export function bei(x: number): number;
export function bei(x: ArrayInput): Float64Array;
export function bei(x: number | ArrayInput): number | Float64Array {
  const xsf = ensureLoaded();

  if (typeof x === 'number') {
    return xsf._wasm_bei(x);
  }

  const arr = x instanceof Float64Array ? x : new Float64Array(x);
  const result = new Float64Array(arr.length);
  for (let i = 0; i < arr.length; i++) {
    result[i] = xsf._wasm_bei(arr[i]);
  }
  return result;
}

/**
 * Kelvin function ker(x).
 *
 * ker(x) = Re[K₀(x·e^(πi/4))]
 *
 * @param x - Input value (x > 0)
 * @returns ker(x)
 */
export function ker(x: number): number;
export function ker(x: ArrayInput): Float64Array;
export function ker(x: number | ArrayInput): number | Float64Array {
  const xsf = ensureLoaded();

  if (typeof x === 'number') {
    return xsf._wasm_ker(x);
  }

  const arr = x instanceof Float64Array ? x : new Float64Array(x);
  const result = new Float64Array(arr.length);
  for (let i = 0; i < arr.length; i++) {
    result[i] = xsf._wasm_ker(arr[i]);
  }
  return result;
}

/**
 * Kelvin function kei(x).
 *
 * kei(x) = Im[K₀(x·e^(πi/4))]
 *
 * @param x - Input value (x > 0)
 * @returns kei(x)
 */
export function kei(x: number): number;
export function kei(x: ArrayInput): Float64Array;
export function kei(x: number | ArrayInput): number | Float64Array {
  const xsf = ensureLoaded();

  if (typeof x === 'number') {
    return xsf._wasm_kei(x);
  }

  const arr = x instanceof Float64Array ? x : new Float64Array(x);
  const result = new Float64Array(arr.length);
  for (let i = 0; i < arr.length; i++) {
    result[i] = xsf._wasm_kei(arr[i]);
  }
  return result;
}

/**
 * Derivative of Kelvin function ber'(x).
 */
export function berp(x: number): number;
export function berp(x: ArrayInput): Float64Array;
export function berp(x: number | ArrayInput): number | Float64Array {
  const xsf = ensureLoaded();

  if (typeof x === 'number') {
    return xsf._wasm_berp(x);
  }

  const arr = x instanceof Float64Array ? x : new Float64Array(x);
  const result = new Float64Array(arr.length);
  for (let i = 0; i < arr.length; i++) {
    result[i] = xsf._wasm_berp(arr[i]);
  }
  return result;
}

/**
 * Derivative of Kelvin function bei'(x).
 */
export function beip(x: number): number;
export function beip(x: ArrayInput): Float64Array;
export function beip(x: number | ArrayInput): number | Float64Array {
  const xsf = ensureLoaded();

  if (typeof x === 'number') {
    return xsf._wasm_beip(x);
  }

  const arr = x instanceof Float64Array ? x : new Float64Array(x);
  const result = new Float64Array(arr.length);
  for (let i = 0; i < arr.length; i++) {
    result[i] = xsf._wasm_beip(arr[i]);
  }
  return result;
}

/**
 * Derivative of Kelvin function ker'(x).
 */
export function kerp(x: number): number;
export function kerp(x: ArrayInput): Float64Array;
export function kerp(x: number | ArrayInput): number | Float64Array {
  const xsf = ensureLoaded();

  if (typeof x === 'number') {
    return xsf._wasm_kerp(x);
  }

  const arr = x instanceof Float64Array ? x : new Float64Array(x);
  const result = new Float64Array(arr.length);
  for (let i = 0; i < arr.length; i++) {
    result[i] = xsf._wasm_kerp(arr[i]);
  }
  return result;
}

/**
 * Derivative of Kelvin function kei'(x).
 */
export function keip(x: number): number;
export function keip(x: ArrayInput): Float64Array;
export function keip(x: number | ArrayInput): number | Float64Array {
  const xsf = ensureLoaded();

  if (typeof x === 'number') {
    return xsf._wasm_keip(x);
  }

  const arr = x instanceof Float64Array ? x : new Float64Array(x);
  const result = new Float64Array(arr.length);
  for (let i = 0; i < arr.length; i++) {
    result[i] = xsf._wasm_keip(arr[i]);
  }
  return result;
}
