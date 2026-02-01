/**
 * Special functions
 *
 * Includes Airy functions, Struve functions, Legendre polynomials, and hypergeometric functions.
 */

import { ensureLoaded, toFloat64Array, type ArrayInput } from '../core/utils.js';

// =============================================================================
// AIRY FUNCTIONS
// =============================================================================

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

  const arr = toFloat64Array(x);
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

// =============================================================================
// STRUVE FUNCTIONS
// =============================================================================

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
 * struveH(0, 1);   // ≈ 0.568
 * struveH(1, 1);   // ≈ 0.198
 * ```
 */
export function struveH(v: number, x: number): number;
export function struveH(v: number, x: ArrayInput): Float64Array;
export function struveH(v: ArrayInput, x: ArrayInput): Float64Array;
export function struveH(
  v: number | ArrayInput,
  x: number | ArrayInput
): number | Float64Array {
  const xsf = ensureLoaded();

  const vIsScalar = typeof v === 'number';
  const xIsScalar = typeof x === 'number';

  if (vIsScalar && xIsScalar) {
    return xsf._wasm_struve_h(v as number, x as number);
  }

  const vArr = vIsScalar ? null : toFloat64Array(v as ArrayInput);
  const xArr = xIsScalar ? null : toFloat64Array(x as ArrayInput);

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
 * struveL(0, 1);   // ≈ 0.710
 * struveL(1, 1);   // ≈ 0.288
 * ```
 */
export function struveL(v: number, x: number): number;
export function struveL(v: number, x: ArrayInput): Float64Array;
export function struveL(v: ArrayInput, x: ArrayInput): Float64Array;
export function struveL(
  v: number | ArrayInput,
  x: number | ArrayInput
): number | Float64Array {
  const xsf = ensureLoaded();

  const vIsScalar = typeof v === 'number';
  const xIsScalar = typeof x === 'number';

  if (vIsScalar && xIsScalar) {
    return xsf._wasm_struve_l(v as number, x as number);
  }

  const vArr = vIsScalar ? null : toFloat64Array(v as ArrayInput);
  const xArr = xIsScalar ? null : toFloat64Array(x as ArrayInput);

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

// =============================================================================
// LEGENDRE POLYNOMIALS
// =============================================================================

/**
 * Legendre polynomial P_n(x).
 *
 * The Legendre polynomials are orthogonal polynomials on [-1, 1].
 *
 * P_0(x) = 1
 * P_1(x) = x
 * P_2(x) = (3x² - 1)/2
 * ...
 *
 * @param n - Degree (non-negative integer)
 * @param x - Argument (-1 ≤ x ≤ 1)
 * @returns P_n(x)
 *
 * @example
 * ```ts
 * legendreP(0, 0.5);  // 1
 * legendreP(1, 0.5);  // 0.5
 * legendreP(2, 0.5);  // -0.125
 * ```
 */
export function legendreP(n: number, x: number): number;
export function legendreP(n: number, x: ArrayInput): Float64Array;
export function legendreP(
  n: number,
  x: number | ArrayInput
): number | Float64Array {
  const xsf = ensureLoaded();

  if (typeof x === 'number') {
    return xsf._wasm_legendre_p(n, x);
  }

  const arr = toFloat64Array(x);
  const result = new Float64Array(arr.length);
  for (let i = 0; i < arr.length; i++) {
    result[i] = xsf._wasm_legendre_p(n, arr[i]);
  }
  return result;
}

// =============================================================================
// HYPERGEOMETRIC FUNCTION
// =============================================================================

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

  const arr = toFloat64Array(x);
  const result = new Float64Array(arr.length);
  for (let i = 0; i < arr.length; i++) {
    result[i] = xsf._wasm_hyp2f1(a, b, c, arr[i]);
  }
  return result;
}
