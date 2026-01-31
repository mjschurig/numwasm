/**
 * Bessel functions - high-level API
 *
 * Provides user-friendly wrappers for Bessel functions of various kinds:
 * - J (first kind)
 * - Y (second kind / Neumann)
 * - I (modified first kind)
 * - K (modified second kind / MacDonald)
 *
 * All functions support both scalar and array inputs.
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

// =============================================================================
// BESSEL FUNCTIONS OF THE FIRST KIND (J)
// =============================================================================

/**
 * Bessel function of the first kind of order 0.
 *
 * J₀(x) is the solution to Bessel's equation x²y'' + xy' + x²y = 0
 * that is finite at the origin.
 *
 * Properties:
 * - J₀(0) = 1
 * - J₀(x) is an even function
 * - Has infinitely many zeros: ~2.405, 5.520, 8.654, ...
 * - |J₀(x)| ≤ 1 for all real x
 *
 * @param x - Input value (any real number)
 * @returns J₀(x)
 *
 * @example
 * ```ts
 * j0(0);      // 1
 * j0(2.405);  // ~0 (first zero)
 * j0([0, 1, 2]);  // Float64Array([1, 0.7651..., 0.2238...])
 * ```
 */
export function j0(x: number): number;
export function j0(x: ArrayInput): Float64Array;
export function j0(x: number | ArrayInput): number | Float64Array {
  const xsf = ensureLoaded();

  if (typeof x === 'number') {
    return xsf._wasm_j0(x);
  }

  const arr = x instanceof Float64Array ? x : new Float64Array(x);
  const result = new Float64Array(arr.length);
  for (let i = 0; i < arr.length; i++) {
    result[i] = xsf._wasm_j0(arr[i]);
  }
  return result;
}

/**
 * Bessel function of the first kind of order 1.
 *
 * J₁(x) is the solution to Bessel's equation for n=1.
 *
 * Properties:
 * - J₁(0) = 0
 * - J₁(x) is an odd function
 * - J₁(x) = -J₀'(x) (derivative relation)
 *
 * @param x - Input value (any real number)
 * @returns J₁(x)
 *
 * @example
 * ```ts
 * j1(0);      // 0
 * j1(1);      // 0.4400...
 * j1(3.832);  // ~0 (first zero)
 * ```
 */
export function j1(x: number): number;
export function j1(x: ArrayInput): Float64Array;
export function j1(x: number | ArrayInput): number | Float64Array {
  const xsf = ensureLoaded();

  if (typeof x === 'number') {
    return xsf._wasm_j1(x);
  }

  const arr = x instanceof Float64Array ? x : new Float64Array(x);
  const result = new Float64Array(arr.length);
  for (let i = 0; i < arr.length; i++) {
    result[i] = xsf._wasm_j1(arr[i]);
  }
  return result;
}

/**
 * Bessel function of the first kind of order v.
 *
 * Jᵥ(x) generalizes J₀ and J₁ to arbitrary real order v.
 *
 * Properties:
 * - Jᵥ(0) = 0 for v > 0, J₀(0) = 1
 * - Recurrence: Jᵥ₊₁(x) = (2v/x)Jᵥ(x) - Jᵥ₋₁(x)
 *
 * @param v - Order of the Bessel function (any real number)
 * @param x - Argument (typically x ≥ 0 for non-integer v)
 * @returns Jᵥ(x)
 *
 * @example
 * ```ts
 * jv(0, 1);      // Same as j0(1)
 * jv(0.5, 1);    // 0.6713...
 * jv(2, [1, 2, 3]);  // Float64Array for J₂ at multiple points
 * ```
 */
export function jv(v: number, x: number): number;
export function jv(v: number, x: ArrayInput): Float64Array;
export function jv(v: ArrayInput, x: ArrayInput): Float64Array;
export function jv(
  v: number | ArrayInput,
  x: number | ArrayInput
): number | Float64Array {
  const xsf = ensureLoaded();

  const vIsScalar = typeof v === 'number';
  const xIsScalar = typeof x === 'number';

  if (vIsScalar && xIsScalar) {
    return xsf._wasm_jv(v as number, x as number);
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
    throw new Error(
      `Array length mismatch: got arrays of length ${vArr.length} and ${xArr.length}`
    );
  }

  const result = new Float64Array(length);
  for (let i = 0; i < length; i++) {
    const vVal = vIsScalar ? (v as number) : vArr![i];
    const xVal = xIsScalar ? (x as number) : xArr![i];
    result[i] = xsf._wasm_jv(vVal, xVal);
  }
  return result;
}

// =============================================================================
// BESSEL FUNCTIONS OF THE SECOND KIND (Y) - NEUMANN FUNCTIONS
// =============================================================================

/**
 * Bessel function of the second kind of order 0 (Neumann function).
 *
 * Y₀(x) is the second linearly independent solution to Bessel's equation for n=0.
 *
 * Properties:
 * - Y₀(x) → -∞ as x → 0⁺ (logarithmic singularity)
 * - Y₀(x) oscillates with decreasing amplitude as x → ∞
 *
 * @param x - Input value (must be positive, x > 0)
 * @returns Y₀(x)
 *
 * @example
 * ```ts
 * y0(1);      // 0.0882...
 * y0(0.894);  // ~0 (first zero)
 * ```
 */
export function y0(x: number): number;
export function y0(x: ArrayInput): Float64Array;
export function y0(x: number | ArrayInput): number | Float64Array {
  const xsf = ensureLoaded();

  if (typeof x === 'number') {
    return xsf._wasm_y0(x);
  }

  const arr = x instanceof Float64Array ? x : new Float64Array(x);
  const result = new Float64Array(arr.length);
  for (let i = 0; i < arr.length; i++) {
    result[i] = xsf._wasm_y0(arr[i]);
  }
  return result;
}

/**
 * Bessel function of the second kind of order 1 (Neumann function).
 *
 * Y₁(x) is the Neumann function of order 1.
 *
 * Properties:
 * - Y₁(x) → -∞ as x → 0⁺
 * - Y₁(x) = -Y₀'(x) (derivative relation)
 *
 * @param x - Input value (must be positive, x > 0)
 * @returns Y₁(x)
 *
 * @example
 * ```ts
 * y1(1);      // -0.7812...
 * y1(2.197);  // ~0 (first zero)
 * ```
 */
export function y1(x: number): number;
export function y1(x: ArrayInput): Float64Array;
export function y1(x: number | ArrayInput): number | Float64Array {
  const xsf = ensureLoaded();

  if (typeof x === 'number') {
    return xsf._wasm_y1(x);
  }

  const arr = x instanceof Float64Array ? x : new Float64Array(x);
  const result = new Float64Array(arr.length);
  for (let i = 0; i < arr.length; i++) {
    result[i] = xsf._wasm_y1(arr[i]);
  }
  return result;
}

/**
 * Bessel function of the second kind of order v (Neumann function).
 *
 * Yᵥ(x) generalizes Y₀ and Y₁ to arbitrary real order v.
 *
 * Properties:
 * - Singular at x = 0 for all v
 * - Recurrence: Yᵥ₊₁(x) = (2v/x)Yᵥ(x) - Yᵥ₋₁(x)
 *
 * @param v - Order of the Bessel function (any real number)
 * @param x - Argument (must be positive, x > 0)
 * @returns Yᵥ(x)
 *
 * @example
 * ```ts
 * yv(0, 1);      // Same as y0(1)
 * yv(0.5, 1);    // -0.4310...
 * yv(2, [1, 2, 3]);  // Float64Array for Y₂ at multiple points
 * ```
 */
export function yv(v: number, x: number): number;
export function yv(v: number, x: ArrayInput): Float64Array;
export function yv(v: ArrayInput, x: ArrayInput): Float64Array;
export function yv(
  v: number | ArrayInput,
  x: number | ArrayInput
): number | Float64Array {
  const xsf = ensureLoaded();

  const vIsScalar = typeof v === 'number';
  const xIsScalar = typeof x === 'number';

  if (vIsScalar && xIsScalar) {
    return xsf._wasm_yv(v as number, x as number);
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
    throw new Error(
      `Array length mismatch: got arrays of length ${vArr.length} and ${xArr.length}`
    );
  }

  const result = new Float64Array(length);
  for (let i = 0; i < length; i++) {
    const vVal = vIsScalar ? (v as number) : vArr![i];
    const xVal = xIsScalar ? (x as number) : xArr![i];
    result[i] = xsf._wasm_yv(vVal, xVal);
  }
  return result;
}

// =============================================================================
// MODIFIED BESSEL FUNCTIONS OF THE FIRST KIND (I)
// =============================================================================

/**
 * Modified Bessel function of the first kind of order 0.
 *
 * I₀(x) is the solution to the modified Bessel equation that is finite at the origin.
 *
 * Properties:
 * - I₀(0) = 1
 * - I₀(x) is an even function
 * - I₀(x) > 0 for all real x
 * - I₀(x) ~ e^x / √(2πx) as x → ∞
 *
 * @param x - Input value (any real number)
 * @returns I₀(x)
 *
 * @example
 * ```ts
 * i0(0);      // 1
 * i0(1);      // 1.2660...
 * i0(10);     // 2815.71...
 * ```
 */
export function i0(x: number): number;
export function i0(x: ArrayInput): Float64Array;
export function i0(x: number | ArrayInput): number | Float64Array {
  const xsf = ensureLoaded();

  if (typeof x === 'number') {
    return xsf._wasm_i0(x);
  }

  const arr = x instanceof Float64Array ? x : new Float64Array(x);
  const result = new Float64Array(arr.length);
  for (let i = 0; i < arr.length; i++) {
    result[i] = xsf._wasm_i0(arr[i]);
  }
  return result;
}

/**
 * Modified Bessel function of the first kind of order 1.
 *
 * I₁(x) is related to I₀ by: I₁(x) = I₀'(x)
 *
 * Properties:
 * - I₁(0) = 0
 * - I₁(x) is an odd function
 * - I₁(x) ~ e^x / √(2πx) as x → ∞
 *
 * @param x - Input value (any real number)
 * @returns I₁(x)
 *
 * @example
 * ```ts
 * i1(0);      // 0
 * i1(1);      // 0.5651...
 * i1(10);     // 2670.98...
 * ```
 */
export function i1(x: number): number;
export function i1(x: ArrayInput): Float64Array;
export function i1(x: number | ArrayInput): number | Float64Array {
  const xsf = ensureLoaded();

  if (typeof x === 'number') {
    return xsf._wasm_i1(x);
  }

  const arr = x instanceof Float64Array ? x : new Float64Array(x);
  const result = new Float64Array(arr.length);
  for (let i = 0; i < arr.length; i++) {
    result[i] = xsf._wasm_i1(arr[i]);
  }
  return result;
}

/**
 * Modified Bessel function of the first kind of order v.
 *
 * Iᵥ(x) generalizes I₀ and I₁ to arbitrary real order v.
 *
 * Properties:
 * - Iᵥ(0) = 0 for v > 0, I₀(0) = 1
 * - Monotonically increasing for v ≥ 0 and x > 0
 *
 * @param v - Order of the Bessel function (any real number)
 * @param x - Argument (any real number for integer v; x ≥ 0 for non-integer v)
 * @returns Iᵥ(x)
 *
 * @example
 * ```ts
 * iv(0, 1);      // Same as i0(1)
 * iv(0.5, 1);    // 0.9376...
 * iv(2, [1, 2, 3]);  // Float64Array for I₂ at multiple points
 * ```
 */
export function iv(v: number, x: number): number;
export function iv(v: number, x: ArrayInput): Float64Array;
export function iv(v: ArrayInput, x: ArrayInput): Float64Array;
export function iv(
  v: number | ArrayInput,
  x: number | ArrayInput
): number | Float64Array {
  const xsf = ensureLoaded();

  const vIsScalar = typeof v === 'number';
  const xIsScalar = typeof x === 'number';

  if (vIsScalar && xIsScalar) {
    return xsf._wasm_iv(v as number, x as number);
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
    throw new Error(
      `Array length mismatch: got arrays of length ${vArr.length} and ${xArr.length}`
    );
  }

  const result = new Float64Array(length);
  for (let i = 0; i < length; i++) {
    const vVal = vIsScalar ? (v as number) : vArr![i];
    const xVal = xIsScalar ? (x as number) : xArr![i];
    result[i] = xsf._wasm_iv(vVal, xVal);
  }
  return result;
}

// =============================================================================
// MODIFIED BESSEL FUNCTIONS OF THE SECOND KIND (K) - MACDONALD FUNCTIONS
// =============================================================================

/**
 * Modified Bessel function of the second kind of order 0 (MacDonald function).
 *
 * K₀(x) is the second linearly independent solution to the modified Bessel equation.
 *
 * Properties:
 * - K₀(x) → -ln(x/2) - γ as x → 0⁺ (logarithmic singularity)
 * - K₀(x) ~ √(π/(2x)) * e^(-x) as x → ∞
 * - K₀(x) > 0 for all x > 0
 * - K₀(x) is monotonically decreasing
 *
 * @param x - Input value (must be positive, x > 0)
 * @returns K₀(x)
 *
 * @example
 * ```ts
 * k0(1);      // 0.4210...
 * k0(0.1);    // 2.4270...
 * k0(10);     // 1.778e-5
 * ```
 */
export function k0(x: number): number;
export function k0(x: ArrayInput): Float64Array;
export function k0(x: number | ArrayInput): number | Float64Array {
  const xsf = ensureLoaded();

  if (typeof x === 'number') {
    return xsf._wasm_k0(x);
  }

  const arr = x instanceof Float64Array ? x : new Float64Array(x);
  const result = new Float64Array(arr.length);
  for (let i = 0; i < arr.length; i++) {
    result[i] = xsf._wasm_k0(arr[i]);
  }
  return result;
}

/**
 * Modified Bessel function of the second kind of order 1 (MacDonald function).
 *
 * K₁(x) is the MacDonald function of order 1.
 *
 * Properties:
 * - K₁(x) → 1/x as x → 0⁺
 * - K₁(x) ~ √(π/(2x)) * e^(-x) as x → ∞
 * - K₁(x) = -K₀'(x) (derivative relation)
 *
 * @param x - Input value (must be positive, x > 0)
 * @returns K₁(x)
 *
 * @example
 * ```ts
 * k1(1);      // 0.6019...
 * k1(0.1);    // 9.8538...
 * k1(10);     // 1.864e-5
 * ```
 */
export function k1(x: number): number;
export function k1(x: ArrayInput): Float64Array;
export function k1(x: number | ArrayInput): number | Float64Array {
  const xsf = ensureLoaded();

  if (typeof x === 'number') {
    return xsf._wasm_k1(x);
  }

  const arr = x instanceof Float64Array ? x : new Float64Array(x);
  const result = new Float64Array(arr.length);
  for (let i = 0; i < arr.length; i++) {
    result[i] = xsf._wasm_k1(arr[i]);
  }
  return result;
}
