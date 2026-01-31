/**
 * Beta function family - high-level API
 *
 * Provides user-friendly wrappers for beta-related special functions
 * with support for both scalar and array inputs.
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
 * Beta function B(a, b).
 *
 * The beta function is defined as:
 *   B(a, b) = ∫₀¹ t^(a-1) * (1-t)^(b-1) dt = Γ(a) * Γ(b) / Γ(a+b)
 *
 * Properties:
 * - B(a, b) = B(b, a) (symmetric)
 * - B(1, 1) = 1
 * - B(a, 1) = 1/a
 *
 * @param a - First parameter (must be positive)
 * @param b - Second parameter (must be positive)
 * @returns B(a, b)
 *
 * @example
 * ```ts
 * beta(2, 3);      // 0.08333... (= 1/12)
 * beta(0.5, 0.5);  // π ≈ 3.14159...
 *
 * // Vectorized with broadcasting
 * beta([1, 2, 3], 2);  // Float64Array([0.5, 0.1666..., 0.0833...])
 * ```
 */
export function beta(a: number, b: number): number;
export function beta(a: ArrayInput, b: number): Float64Array;
export function beta(a: number, b: ArrayInput): Float64Array;
export function beta(a: ArrayInput, b: ArrayInput): Float64Array;
export function beta(
  a: number | ArrayInput,
  b: number | ArrayInput
): number | Float64Array {
  const xsf = ensureLoaded();

  const aIsScalar = typeof a === 'number';
  const bIsScalar = typeof b === 'number';

  if (aIsScalar && bIsScalar) {
    return xsf._wasm_beta(a as number, b as number);
  }

  const aArr = aIsScalar
    ? null
    : a instanceof Float64Array
      ? a
      : new Float64Array(a as number[]);
  const bArr = bIsScalar
    ? null
    : b instanceof Float64Array
      ? b
      : new Float64Array(b as number[]);

  const length = aArr?.length ?? bArr!.length;
  if (aArr && bArr && aArr.length !== bArr.length) {
    throw new Error(
      `Array length mismatch: got arrays of length ${aArr.length} and ${bArr.length}`
    );
  }

  const result = new Float64Array(length);
  for (let i = 0; i < length; i++) {
    const aVal = aIsScalar ? (a as number) : aArr![i];
    const bVal = bIsScalar ? (b as number) : bArr![i];
    result[i] = xsf._wasm_beta(aVal, bVal);
  }
  return result;
}

/**
 * Natural logarithm of the beta function.
 *
 * betaln(a, b) = ln(B(a, b)) = ln(Γ(a)) + ln(Γ(b)) - ln(Γ(a+b))
 *
 * More numerically stable than log(beta(a, b)) for large arguments.
 * Essential for computing binomial coefficients with large n.
 *
 * @param a - First parameter (must be positive)
 * @param b - Second parameter (must be positive)
 * @returns ln(B(a, b))
 *
 * @example
 * ```ts
 * betaln(100, 200);  // -202.58... (beta would underflow)
 * betaln(2, 3);      // -2.4849... (= ln(1/12))
 * ```
 */
export function betaln(a: number, b: number): number;
export function betaln(a: ArrayInput, b: number): Float64Array;
export function betaln(a: number, b: ArrayInput): Float64Array;
export function betaln(a: ArrayInput, b: ArrayInput): Float64Array;
export function betaln(
  a: number | ArrayInput,
  b: number | ArrayInput
): number | Float64Array {
  const xsf = ensureLoaded();

  const aIsScalar = typeof a === 'number';
  const bIsScalar = typeof b === 'number';

  if (aIsScalar && bIsScalar) {
    return xsf._wasm_betaln(a as number, b as number);
  }

  const aArr = aIsScalar
    ? null
    : a instanceof Float64Array
      ? a
      : new Float64Array(a as number[]);
  const bArr = bIsScalar
    ? null
    : b instanceof Float64Array
      ? b
      : new Float64Array(b as number[]);

  const length = aArr?.length ?? bArr!.length;
  if (aArr && bArr && aArr.length !== bArr.length) {
    throw new Error(
      `Array length mismatch: got arrays of length ${aArr.length} and ${bArr.length}`
    );
  }

  const result = new Float64Array(length);
  for (let i = 0; i < length; i++) {
    const aVal = aIsScalar ? (a as number) : aArr![i];
    const bVal = bIsScalar ? (b as number) : bArr![i];
    result[i] = xsf._wasm_betaln(aVal, bVal);
  }
  return result;
}
