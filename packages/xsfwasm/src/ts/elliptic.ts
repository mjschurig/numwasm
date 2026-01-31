/**
 * Elliptic integrals and Jacobi elliptic functions - high-level API
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
 * Complete elliptic integral of the first kind K(m).
 *
 * K(m) = ∫₀^(π/2) 1/√(1 - m·sin²θ) dθ
 *
 * @param m - Parameter (0 ≤ m < 1)
 * @returns K(m)
 *
 * @example
 * ```ts
 * ellipk(0);     // π/2 ≈ 1.5708
 * ellipk(0.5);   // 1.8541...
 * ```
 */
export function ellipk(m: number): number;
export function ellipk(m: ArrayInput): Float64Array;
export function ellipk(m: number | ArrayInput): number | Float64Array {
  const xsf = ensureLoaded();

  if (typeof m === 'number') {
    return xsf._wasm_ellipk(m);
  }

  const arr = m instanceof Float64Array ? m : new Float64Array(m);
  const result = new Float64Array(arr.length);
  for (let i = 0; i < arr.length; i++) {
    result[i] = xsf._wasm_ellipk(arr[i]);
  }
  return result;
}

/**
 * Complete elliptic integral of the second kind E(m).
 *
 * E(m) = ∫₀^(π/2) √(1 - m·sin²θ) dθ
 *
 * @param m - Parameter (0 ≤ m ≤ 1)
 * @returns E(m)
 *
 * @example
 * ```ts
 * ellipe(0);     // π/2 ≈ 1.5708
 * ellipe(1);     // 1
 * ```
 */
export function ellipe(m: number): number;
export function ellipe(m: ArrayInput): Float64Array;
export function ellipe(m: number | ArrayInput): number | Float64Array {
  const xsf = ensureLoaded();

  if (typeof m === 'number') {
    return xsf._wasm_ellipe(m);
  }

  const arr = m instanceof Float64Array ? m : new Float64Array(m);
  const result = new Float64Array(arr.length);
  for (let i = 0; i < arr.length; i++) {
    result[i] = xsf._wasm_ellipe(arr[i]);
  }
  return result;
}

/**
 * Incomplete elliptic integral of the first kind F(φ|m).
 *
 * F(φ|m) = ∫₀^φ 1/√(1 - m·sin²θ) dθ
 *
 * @param phi - Amplitude (radians)
 * @param m - Parameter (0 ≤ m < 1)
 * @returns F(φ|m)
 */
export function ellipkinc(phi: number, m: number): number;
export function ellipkinc(phi: ArrayInput, m: number): Float64Array;
export function ellipkinc(phi: number, m: ArrayInput): Float64Array;
export function ellipkinc(phi: ArrayInput, m: ArrayInput): Float64Array;
export function ellipkinc(
  phi: number | ArrayInput,
  m: number | ArrayInput
): number | Float64Array {
  const xsf = ensureLoaded();

  const phiIsScalar = typeof phi === 'number';
  const mIsScalar = typeof m === 'number';

  if (phiIsScalar && mIsScalar) {
    return xsf._wasm_ellipkinc(phi as number, m as number);
  }

  const phiArr = phiIsScalar
    ? null
    : phi instanceof Float64Array
      ? phi
      : new Float64Array(phi as number[]);
  const mArr = mIsScalar
    ? null
    : m instanceof Float64Array
      ? m
      : new Float64Array(m as number[]);

  const length = phiArr?.length ?? mArr!.length;
  if (phiArr && mArr && phiArr.length !== mArr.length) {
    throw new Error(`Array length mismatch`);
  }

  const result = new Float64Array(length);
  for (let i = 0; i < length; i++) {
    const phiVal = phiIsScalar ? (phi as number) : phiArr![i];
    const mVal = mIsScalar ? (m as number) : mArr![i];
    result[i] = xsf._wasm_ellipkinc(phiVal, mVal);
  }
  return result;
}

/**
 * Incomplete elliptic integral of the second kind E(φ|m).
 *
 * E(φ|m) = ∫₀^φ √(1 - m·sin²θ) dθ
 *
 * @param phi - Amplitude (radians)
 * @param m - Parameter (0 ≤ m ≤ 1)
 * @returns E(φ|m)
 */
export function ellipeinc(phi: number, m: number): number;
export function ellipeinc(phi: ArrayInput, m: number): Float64Array;
export function ellipeinc(phi: number, m: ArrayInput): Float64Array;
export function ellipeinc(phi: ArrayInput, m: ArrayInput): Float64Array;
export function ellipeinc(
  phi: number | ArrayInput,
  m: number | ArrayInput
): number | Float64Array {
  const xsf = ensureLoaded();

  const phiIsScalar = typeof phi === 'number';
  const mIsScalar = typeof m === 'number';

  if (phiIsScalar && mIsScalar) {
    return xsf._wasm_ellipeinc(phi as number, m as number);
  }

  const phiArr = phiIsScalar
    ? null
    : phi instanceof Float64Array
      ? phi
      : new Float64Array(phi as number[]);
  const mArr = mIsScalar
    ? null
    : m instanceof Float64Array
      ? m
      : new Float64Array(m as number[]);

  const length = phiArr?.length ?? mArr!.length;
  if (phiArr && mArr && phiArr.length !== mArr.length) {
    throw new Error(`Array length mismatch`);
  }

  const result = new Float64Array(length);
  for (let i = 0; i < length; i++) {
    const phiVal = phiIsScalar ? (phi as number) : phiArr![i];
    const mVal = mIsScalar ? (m as number) : mArr![i];
    result[i] = xsf._wasm_ellipeinc(phiVal, mVal);
  }
  return result;
}

export interface EllipjResult {
  sn: number;
  cn: number;
  dn: number;
  ph: number;
}

export interface EllipjArrayResult {
  sn: Float64Array;
  cn: Float64Array;
  dn: Float64Array;
  ph: Float64Array;
}

/**
 * Jacobi elliptic functions sn, cn, dn.
 *
 * @param u - Argument
 * @param m - Parameter (0 ≤ m ≤ 1)
 * @returns Object with sn, cn, dn, ph
 *
 * @example
 * ```ts
 * const { sn, cn, dn } = ellipj(0, 0.5);  // sn=0, cn=1, dn=1
 * ```
 */
export function ellipj(u: number, m: number): EllipjResult;
export function ellipj(u: ArrayInput, m: number): EllipjArrayResult;
export function ellipj(u: number, m: ArrayInput): EllipjArrayResult;
export function ellipj(u: ArrayInput, m: ArrayInput): EllipjArrayResult;
export function ellipj(
  u: number | ArrayInput,
  m: number | ArrayInput
): EllipjResult | EllipjArrayResult {
  const xsf = ensureLoaded();

  const uIsScalar = typeof u === 'number';
  const mIsScalar = typeof m === 'number';

  if (uIsScalar && mIsScalar) {
    const snPtr = xsf._malloc(8);
    const cnPtr = xsf._malloc(8);
    const dnPtr = xsf._malloc(8);
    const phPtr = xsf._malloc(8);
    try {
      xsf._wasm_ellipj(u as number, m as number, snPtr, cnPtr, dnPtr, phPtr);
      return {
        sn: xsf.HEAPF64[snPtr >> 3],
        cn: xsf.HEAPF64[cnPtr >> 3],
        dn: xsf.HEAPF64[dnPtr >> 3],
        ph: xsf.HEAPF64[phPtr >> 3],
      };
    } finally {
      xsf._free(snPtr);
      xsf._free(cnPtr);
      xsf._free(dnPtr);
      xsf._free(phPtr);
    }
  }

  const uArr = uIsScalar
    ? null
    : u instanceof Float64Array
      ? u
      : new Float64Array(u as number[]);
  const mArr = mIsScalar
    ? null
    : m instanceof Float64Array
      ? m
      : new Float64Array(m as number[]);

  const length = uArr?.length ?? mArr!.length;
  if (uArr && mArr && uArr.length !== mArr.length) {
    throw new Error(`Array length mismatch`);
  }

  const sn = new Float64Array(length);
  const cn = new Float64Array(length);
  const dn = new Float64Array(length);
  const ph = new Float64Array(length);

  const snPtr = xsf._malloc(8);
  const cnPtr = xsf._malloc(8);
  const dnPtr = xsf._malloc(8);
  const phPtr = xsf._malloc(8);
  try {
    for (let i = 0; i < length; i++) {
      const uVal = uIsScalar ? (u as number) : uArr![i];
      const mVal = mIsScalar ? (m as number) : mArr![i];
      xsf._wasm_ellipj(uVal, mVal, snPtr, cnPtr, dnPtr, phPtr);
      sn[i] = xsf.HEAPF64[snPtr >> 3];
      cn[i] = xsf.HEAPF64[cnPtr >> 3];
      dn[i] = xsf.HEAPF64[dnPtr >> 3];
      ph[i] = xsf.HEAPF64[phPtr >> 3];
    }
  } finally {
    xsf._free(snPtr);
    xsf._free(cnPtr);
    xsf._free(dnPtr);
    xsf._free(phPtr);
  }

  return { sn, cn, dn, ph };
}
