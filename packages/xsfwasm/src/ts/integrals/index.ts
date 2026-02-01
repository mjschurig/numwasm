/**
 * Integral functions
 *
 * Includes elliptic integrals, exponential integrals, Fresnel integrals, and sine/cosine integrals.
 */

import { ensureLoaded, toFloat64Array, type ArrayInput } from '../core/utils.js';

// =============================================================================
// ELLIPTIC INTEGRALS
// =============================================================================

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

  const arr = toFloat64Array(m);
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

  const arr = toFloat64Array(m);
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

  const phiArr = phiIsScalar ? null : toFloat64Array(phi as ArrayInput);
  const mArr = mIsScalar ? null : toFloat64Array(m as ArrayInput);

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

  const phiArr = phiIsScalar ? null : toFloat64Array(phi as ArrayInput);
  const mArr = mIsScalar ? null : toFloat64Array(m as ArrayInput);

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

  const uArr = uIsScalar ? null : toFloat64Array(u as ArrayInput);
  const mArr = mIsScalar ? null : toFloat64Array(m as ArrayInput);

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

// =============================================================================
// EXPONENTIAL INTEGRALS
// =============================================================================

/**
 * Exponential integral E₁(x).
 *
 * E₁(x) = ∫ₓ^∞ e^(-t)/t dt  for x > 0
 *
 * @param x - Input value (x > 0)
 * @returns E₁(x)
 *
 * @example
 * ```ts
 * exp1(1);    // 0.2193...
 * exp1(0);    // Infinity
 * ```
 */
export function exp1(x: number): number;
export function exp1(x: ArrayInput): Float64Array;
export function exp1(x: number | ArrayInput): number | Float64Array {
  const xsf = ensureLoaded();

  if (typeof x === 'number') {
    return xsf._wasm_exp1(x);
  }

  const arr = toFloat64Array(x);
  const result = new Float64Array(arr.length);
  for (let i = 0; i < arr.length; i++) {
    result[i] = xsf._wasm_exp1(arr[i]);
  }
  return result;
}

/**
 * Exponential integral Ei(x).
 *
 * Ei(x) = ∫₋∞^x e^t/t dt  (principal value)
 *
 * @param x - Input value
 * @returns Ei(x)
 *
 * @example
 * ```ts
 * expi(1);    // 1.8951...
 * ```
 */
export function expi(x: number): number;
export function expi(x: ArrayInput): Float64Array;
export function expi(x: number | ArrayInput): number | Float64Array {
  const xsf = ensureLoaded();

  if (typeof x === 'number') {
    return xsf._wasm_expi(x);
  }

  const arr = toFloat64Array(x);
  const result = new Float64Array(arr.length);
  for (let i = 0; i < arr.length; i++) {
    result[i] = xsf._wasm_expi(arr[i]);
  }
  return result;
}

// =============================================================================
// FRESNEL INTEGRALS
// =============================================================================

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

  const arr = toFloat64Array(x);
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

// =============================================================================
// SINE AND COSINE INTEGRALS
// =============================================================================

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

  const arr = toFloat64Array(x);
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

  const arr = toFloat64Array(x);
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
