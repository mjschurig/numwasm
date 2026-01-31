/**
 * Special mathematical ufuncs for NumJS-WASM
 *
 * Miscellaneous special functions: sign operations, complex numbers, etc.
 */

import { NDArray } from "../_core/NDArray.js";
import { DType } from "../types.js";
import { isIntegerDType } from "../dtype.js";
import { getWasmModule } from "../wasm-loader.js";
import { applyUnary, applyBinary } from "./helpers.js";

/* ============ Sign and Comparison ============ */

/** Change the sign of x1 to that of x2, element-wise. */
export function copysign(x1: NDArray, x2: NDArray): NDArray {
  return applyBinary(x1, x2, (p1, p2) =>
    getWasmModule()._ufunc_copysign(p1, p2),
  );
}

/** Return element-wise True where signbit is set (less than zero). */
export function signbit(x: NDArray): NDArray {
  return applyUnary(x, (ptr) => getWasmModule()._ufunc_signbit(ptr));
}

/* ============ Floating Point Operations ============ */

/**
 * Decompose the elements of x into mantissa and twos exponent.
 * Returns (mantissa, exponent), where x = mantissa * 2^exponent.
 * The mantissa lies in the open-closed interval (-1, 1) for float64,
 * and the twos exponent is a signed integer.
 *
 * @param x - Input array
 * @returns Tuple of [mantissa, exponent] arrays
 *
 * @example
 * const [m, e] = frexp(arr);
 * // For [1.0, 2.0, 4.0, 8.0]:
 * // m: [0.5, 0.5, 0.5, 0.5]
 * // e: [1, 2, 3, 4]
 */
export function frexp(x: NDArray): [NDArray, NDArray] {
  const module = getWasmModule();
  const resultPtr = module._ufunc_frexp(x._wasmPtr);
  if (resultPtr === 0) {
    throw new Error("frexp operation failed");
  }

  const mantissaPtr = module._ufunc_tuple_get_first(resultPtr);
  const exponentPtr = module._ufunc_tuple_get_second(resultPtr);
  module._ufunc_tuple_result_free(resultPtr);

  return [
    NDArray._fromPtr(mantissaPtr, module),
    NDArray._fromPtr(exponentPtr, module),
  ];
}

/**
 * Returns x1 * 2^x2, element-wise.
 *
 * The mantissas x1 and twos exponents x2 are used to construct
 * floating point numbers x = x1 * 2^x2.
 *
 * @param x1 - Array of multipliers (mantissas)
 * @param x2 - Array of twos exponents
 * @returns x1 * 2^x2
 *
 * @example
 * ldexp([0.5, 0.5, 0.5], [1, 2, 3])  // [1.0, 2.0, 4.0]
 */
export function ldexp(x1: NDArray, x2: NDArray): NDArray {
  return applyBinary(x1, x2, (p1, p2) => getWasmModule()._ufunc_ldexp(p1, p2));
}

/**
 * Return the next floating-point value after x1 towards x2, element-wise.
 *
 * @param x1 - Starting values
 * @param x2 - Direction values
 * @returns Next representable floating-point values
 *
 * @example
 * nextafter(1.0, 2.0)  // Slightly larger than 1.0
 * nextafter(1.0, 0.0)  // Slightly smaller than 1.0
 */
export function nextafter(x1: NDArray, x2: NDArray): NDArray {
  return applyBinary(x1, x2, (p1, p2) =>
    getWasmModule()._ufunc_nextafter(p1, p2),
  );
}

/**
 * Return the distance between x and the nearest adjacent number.
 * This is effectively the ulp (unit in the last place) of x.
 *
 * @param x - Input values
 * @returns Distance to next representable value
 *
 * @example
 * spacing(1.0)  // ~2.22e-16 (machine epsilon at 1.0)
 */
export function spacing(x: NDArray): NDArray {
  return applyUnary(x, (ptr) => getWasmModule()._ufunc_spacing(ptr));
}

/**
 * Return the fractional and integral parts of an array, element-wise.
 * The fractional and integral parts are negative if the given number is negative.
 *
 * @param x - Input array
 * @returns Tuple of [fractional part, integral part]
 *
 * @example
 * const [frac, intg] = modf([3.5, -2.7]);
 * // frac: [0.5, -0.7]
 * // intg: [3.0, -2.0]
 */
export function modf(x: NDArray): [NDArray, NDArray] {
  const module = getWasmModule();
  const resultPtr = module._ufunc_modf(x._wasmPtr);
  if (resultPtr === 0) {
    throw new Error("modf operation failed");
  }

  const fracPtr = module._ufunc_tuple_get_first(resultPtr);
  const intPtr = module._ufunc_tuple_get_second(resultPtr);
  module._ufunc_tuple_result_free(resultPtr);

  return [NDArray._fromPtr(fracPtr, module), NDArray._fromPtr(intPtr, module)];
}

/* ============ Integer Operations ============ */

/**
 * Returns the greatest common divisor of |x1| and |x2|.
 * Input arrays must be integers.
 *
 * @param x1 - First array of integers
 * @param x2 - Second array of integers
 * @returns Greatest common divisor
 *
 * @example
 * gcd(12, 8)  // 4
 * gcd([12, 15, 20], [8, 10, 15])  // [4, 5, 5]
 * gcd(0, 5)  // 5
 */
export function gcd(x1: NDArray, x2: NDArray): NDArray {
  if (!isIntegerDType(x1.dtype) && x1.dtype !== DType.Bool) {
    throw new TypeError("gcd requires integer inputs");
  }
  if (!isIntegerDType(x2.dtype) && x2.dtype !== DType.Bool) {
    throw new TypeError("gcd requires integer inputs");
  }
  return applyBinary(x1, x2, (p1, p2) => getWasmModule()._ufunc_gcd(p1, p2));
}

/**
 * Returns the lowest common multiple of |x1| and |x2|.
 * Input arrays must be integers.
 *
 * @param x1 - First array of integers
 * @param x2 - Second array of integers
 * @returns Lowest common multiple
 *
 * @example
 * lcm(12, 8)  // 24
 * lcm([4, 6], [8, 9])  // [8, 18]
 * lcm(0, 5)  // 0
 */
export function lcm(x1: NDArray, x2: NDArray): NDArray {
  if (!isIntegerDType(x1.dtype) && x1.dtype !== DType.Bool) {
    throw new TypeError("lcm requires integer inputs");
  }
  if (!isIntegerDType(x2.dtype) && x2.dtype !== DType.Bool) {
    throw new TypeError("lcm requires integer inputs");
  }
  return applyBinary(x1, x2, (p1, p2) => getWasmModule()._ufunc_lcm(p1, p2));
}

/* ============ Special Functions ============ */

/**
 * Return the normalized sinc function: sin(pi*x) / (pi*x).
 * The sinc function is used in various signal processing applications.
 *
 * @param x - Input array
 * @returns sinc(x)
 *
 * @example
 * sinc(0)    // 1.0 (limit as x -> 0)
 * sinc(1)    // ~0 (sin(pi)/pi)
 * sinc(0.5)  // ~0.637 (2/pi)
 */
export function sinc(x: NDArray): NDArray {
  return applyUnary(x, (ptr) => getWasmModule()._ufunc_sinc(ptr));
}

/**
 * Compute the Heaviside step function.
 *
 * H(x) = 0 if x < 0
 * H(x) = h0 if x == 0
 * H(x) = 1 if x > 0
 *
 * @param x1 - Input values
 * @param x2 - Value of the function at x1 == 0 (h0)
 * @returns Heaviside step function values
 *
 * @example
 * heaviside([-1, 0, 1], 0.5)  // [0, 0.5, 1]
 */
export function heaviside(x1: NDArray, x2: NDArray): NDArray {
  return applyBinary(x1, x2, (p1, p2) =>
    getWasmModule()._ufunc_heaviside(p1, p2),
  );
}

/**
 * Return element-wise quotient and remainder simultaneously.
 * Equivalent to (floor_divide(x1, x2), remainder(x1, x2)), but
 * returns both in one call.
 *
 * @param x1 - Dividend array
 * @param x2 - Divisor array
 * @returns Tuple of [quotient, remainder]
 *
 * @example
 * const [q, r] = divmod(10, 3);
 * // q: 3, r: 1
 *
 * const [q, r] = divmod([10, 11, 12], 3);
 * // q: [3, 3, 4], r: [1, 2, 0]
 */
export function divmod(x1: NDArray, x2: NDArray): [NDArray, NDArray] {
  const module = getWasmModule();
  const resultPtr = module._ufunc_divmod(x1._wasmPtr, x2._wasmPtr);
  if (resultPtr === 0) {
    throw new Error("divmod operation failed");
  }

  const quotPtr = module._ufunc_tuple_get_first(resultPtr);
  const remPtr = module._ufunc_tuple_get_second(resultPtr);
  module._ufunc_tuple_result_free(resultPtr);

  return [NDArray._fromPtr(quotPtr, module), NDArray._fromPtr(remPtr, module)];
}

/* ============ Complex Number Operations ============ */

/**
 * Return the complex conjugate, element-wise.
 *
 * The complex conjugate of a complex number is obtained by changing
 * the sign of its imaginary part.
 *
 * For non-complex arrays, returns a copy of the input.
 *
 * @param x - Input array
 * @returns The complex conjugate of x, with same dtype as input
 *
 * @example
 * ```typescript
 * const c = await NDArray.fromArray([1+2j, 3-4j], { dtype: DType.Complex128 });
 * const conj = conjugate(c);  // [1-2j, 3+4j]
 * ```
 */
export function conjugate(x: NDArray): NDArray {
  // For non-complex types, just return a copy
  if (x.dtype !== DType.Complex64 && x.dtype !== DType.Complex128) {
    return x.copy();
  }

  // For complex types, negate the imaginary part
  const result = x.copy();
  const size = x.size;
  const module = x._wasmModule;

  for (let flatIdx = 0; flatIdx < size; flatIdx++) {
    const realVal = module._ndarray_get_complex_real(x._wasmPtr, flatIdx);
    const imagVal = module._ndarray_get_complex_imag(x._wasmPtr, flatIdx);
    module._ndarray_set_complex(result._wasmPtr, flatIdx, realVal, -imagVal);
  }

  return result;
}

/** Alias for conjugate. */
export const conj = conjugate;
