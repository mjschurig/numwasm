/**
 * Error functions
 *
 * Includes erf, erfc, erfcx, and erfi functions.
 */

import { ensureLoaded, toFloat64Array, type ArrayInput } from '../core/utils.js';

/**
 * Error function erf(x).
 *
 * The error function is defined as:
 *   erf(x) = (2/√π) * ∫₀ˣ e^(-t²) dt
 *
 * Properties:
 * - erf(0) = 0
 * - erf(∞) = 1
 * - erf(-x) = -erf(x) (odd function)
 * - Related to normal CDF: Φ(x) = (1 + erf(x/√2))/2
 *
 * @param x - Input value
 * @returns erf(x), ranging from -1 to 1
 *
 * @example
 * ```ts
 * erf(0);     // 0
 * erf(1);     // 0.8427...
 * erf(-1);    // -0.8427...
 * erf([0, 1, 2]);  // Float64Array([0, 0.8427..., 0.9953...])
 * ```
 */
export function erf(x: number): number;
export function erf(x: ArrayInput): Float64Array;
export function erf(x: number | ArrayInput): number | Float64Array {
  const xsf = ensureLoaded();

  if (typeof x === 'number') {
    return xsf._wasm_erf(x);
  }

  const arr = toFloat64Array(x);
  const result = new Float64Array(arr.length);
  for (let i = 0; i < arr.length; i++) {
    result[i] = xsf._wasm_erf(arr[i]);
  }
  return result;
}

/**
 * Complementary error function erfc(x).
 *
 * erfc(x) = 1 - erf(x) = (2/√π) * ∫ₓ^∞ e^(-t²) dt
 *
 * More accurate than 1 - erf(x) for large x where erf(x) ≈ 1.
 * Essential for computing tail probabilities.
 *
 * Properties:
 * - erfc(0) = 1
 * - erfc(∞) = 0
 * - erfc(-∞) = 2
 *
 * @param x - Input value
 * @returns erfc(x) = 1 - erf(x)
 *
 * @example
 * ```ts
 * erfc(0);    // 1
 * erfc(1);    // 0.1572...
 * erfc(5);    // 1.5374e-12 (more accurate than 1 - erf(5))
 * ```
 */
export function erfc(x: number): number;
export function erfc(x: ArrayInput): Float64Array;
export function erfc(x: number | ArrayInput): number | Float64Array {
  const xsf = ensureLoaded();

  if (typeof x === 'number') {
    return xsf._wasm_erfc(x);
  }

  const arr = toFloat64Array(x);
  const result = new Float64Array(arr.length);
  for (let i = 0; i < arr.length; i++) {
    result[i] = xsf._wasm_erfc(arr[i]);
  }
  return result;
}

/**
 * Scaled complementary error function erfcx(x).
 *
 * erfcx(x) = e^(x²) * erfc(x)
 *
 * This scaled version avoids underflow for large positive x where
 * erfc(x) becomes extremely small.
 *
 * Properties:
 * - erfcx(0) = 1
 * - erfcx(x) → 1/(x√π) as x → ∞
 * - erfcx(x) is always positive
 *
 * @param x - Input value
 * @returns erfcx(x) = exp(x²) * erfc(x)
 *
 * @example
 * ```ts
 * erfcx(0);    // 1
 * erfcx(10);   // 0.0561... (erfc(10) ≈ 2e-45 would underflow)
 * erfcx(100);  // 0.00564... (stable computation)
 * ```
 */
export function erfcx(x: number): number;
export function erfcx(x: ArrayInput): Float64Array;
export function erfcx(x: number | ArrayInput): number | Float64Array {
  const xsf = ensureLoaded();

  if (typeof x === 'number') {
    return xsf._wasm_erfcx(x);
  }

  const arr = toFloat64Array(x);
  const result = new Float64Array(arr.length);
  for (let i = 0; i < arr.length; i++) {
    result[i] = xsf._wasm_erfcx(arr[i]);
  }
  return result;
}

/**
 * Imaginary error function erfi(x).
 *
 * erfi(x) = -i * erf(ix) = (2/√π) * ∫₀ˣ e^(t²) dt
 *
 * Unlike erf(x), erfi(x) grows without bound for large |x|.
 *
 * Properties:
 * - erfi(0) = 0
 * - erfi(-x) = -erfi(x) (odd function)
 * - erfi(x) ~ e^(x²)/(x√π) for large x
 *
 * @param x - Input value
 * @returns erfi(x)
 *
 * @example
 * ```ts
 * erfi(0);    // 0
 * erfi(1);    // 1.6504...
 * erfi(2);    // 18.564...
 * ```
 */
export function erfi(x: number): number;
export function erfi(x: ArrayInput): Float64Array;
export function erfi(x: number | ArrayInput): number | Float64Array {
  const xsf = ensureLoaded();

  if (typeof x === 'number') {
    return xsf._wasm_erfi(x);
  }

  const arr = toFloat64Array(x);
  const result = new Float64Array(arr.length);
  for (let i = 0; i < arr.length; i++) {
    result[i] = xsf._wasm_erfi(arr[i]);
  }
  return result;
}
