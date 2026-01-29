/**
 * NumJS Universal Functions (Ufuncs)
 *
 * Element-wise mathematical operations on NDArrays with broadcasting support.
 * These functions call optimized WASM implementations for performance.
 */

import { NDArray } from './NDArray.js';
import { getWasmModule } from './wasm-loader.js';
import { isIntegerDType } from './dtype.js';
import { DType } from './types.js';

/* ============ Helper Functions ============ */

/**
 * Helper to apply a unary ufunc.
 * @internal
 */
function applyUnary(
  input: NDArray,
  wasmFunc: (ptr: number) => number
): NDArray {
  const module = getWasmModule();
  const resultPtr = wasmFunc(input._wasmPtr);
  if (resultPtr === 0) {
    throw new Error('Ufunc operation failed');
  }
  return NDArray._fromPtr(resultPtr, module);
}

/**
 * Helper to apply a binary ufunc.
 * @internal
 */
function applyBinary(
  x1: NDArray,
  x2: NDArray,
  wasmFunc: (ptr1: number, ptr2: number) => number
): NDArray {
  const module = getWasmModule();
  const resultPtr = wasmFunc(x1._wasmPtr, x2._wasmPtr);
  if (resultPtr === 0) {
    throw new Error('Ufunc operation failed (possibly due to incompatible shapes for broadcasting)');
  }
  return NDArray._fromPtr(resultPtr, module);
}

/* ============ Arithmetic Operations ============ */

/** Numerical negative, element-wise. */
export function negative(x: NDArray): NDArray {
  return applyUnary(x, (ptr) => getWasmModule()._ufunc_negative(ptr));
}

/** Numerical positive (returns a copy). */
export function positive(x: NDArray): NDArray {
  return applyUnary(x, (ptr) => getWasmModule()._ufunc_positive(ptr));
}

/** Calculate the absolute value element-wise. */
export function absolute(x: NDArray): NDArray {
  return applyUnary(x, (ptr) => getWasmModule()._ufunc_absolute(ptr));
}

/** Alias for absolute. */
export const abs = absolute;

/** Returns element-wise sign: -1, 0, or +1. */
export function sign(x: NDArray): NDArray {
  return applyUnary(x, (ptr) => getWasmModule()._ufunc_sign(ptr));
}

/* ============ Powers and Roots ============ */

/** Return the non-negative square-root of an array, element-wise. */
export function sqrt(x: NDArray): NDArray {
  return applyUnary(x, (ptr) => getWasmModule()._ufunc_sqrt(ptr));
}

/** Return the element-wise square of the input. */
export function square(x: NDArray): NDArray {
  return applyUnary(x, (ptr) => getWasmModule()._ufunc_square(ptr));
}

/** Return the cube-root of an array, element-wise. */
export function cbrt(x: NDArray): NDArray {
  return applyUnary(x, (ptr) => getWasmModule()._ufunc_cbrt(ptr));
}

/** Return the reciprocal of the argument, element-wise. */
export function reciprocal(x: NDArray): NDArray {
  return applyUnary(x, (ptr) => getWasmModule()._ufunc_reciprocal(ptr));
}

/* ============ Exponential Functions ============ */

/** Calculate the exponential of all elements in the input array. */
export function exp(x: NDArray): NDArray {
  return applyUnary(x, (ptr) => getWasmModule()._ufunc_exp(ptr));
}

/** Calculate 2**x for all elements in the array. */
export function exp2(x: NDArray): NDArray {
  return applyUnary(x, (ptr) => getWasmModule()._ufunc_exp2(ptr));
}

/** Calculate exp(x) - 1 for all elements in the array. */
export function expm1(x: NDArray): NDArray {
  return applyUnary(x, (ptr) => getWasmModule()._ufunc_expm1(ptr));
}

/* ============ Logarithmic Functions ============ */

/** Natural logarithm, element-wise. */
export function log(x: NDArray): NDArray {
  return applyUnary(x, (ptr) => getWasmModule()._ufunc_log(ptr));
}

/** Base-2 logarithm of x. */
export function log2(x: NDArray): NDArray {
  return applyUnary(x, (ptr) => getWasmModule()._ufunc_log2(ptr));
}

/** Base-10 logarithm, element-wise. */
export function log10(x: NDArray): NDArray {
  return applyUnary(x, (ptr) => getWasmModule()._ufunc_log10(ptr));
}

/** Return the natural logarithm of one plus the input array, element-wise. */
export function log1p(x: NDArray): NDArray {
  return applyUnary(x, (ptr) => getWasmModule()._ufunc_log1p(ptr));
}

/* ============ Trigonometric Functions ============ */

/** Trigonometric sine, element-wise. */
export function sin(x: NDArray): NDArray {
  return applyUnary(x, (ptr) => getWasmModule()._ufunc_sin(ptr));
}

/** Cosine element-wise. */
export function cos(x: NDArray): NDArray {
  return applyUnary(x, (ptr) => getWasmModule()._ufunc_cos(ptr));
}

/** Compute tangent element-wise. */
export function tan(x: NDArray): NDArray {
  return applyUnary(x, (ptr) => getWasmModule()._ufunc_tan(ptr));
}

/** Inverse sine, element-wise. */
export function arcsin(x: NDArray): NDArray {
  return applyUnary(x, (ptr) => getWasmModule()._ufunc_arcsin(ptr));
}

/** Trigonometric inverse cosine, element-wise. */
export function arccos(x: NDArray): NDArray {
  return applyUnary(x, (ptr) => getWasmModule()._ufunc_arccos(ptr));
}

/** Trigonometric inverse tangent, element-wise. */
export function arctan(x: NDArray): NDArray {
  return applyUnary(x, (ptr) => getWasmModule()._ufunc_arctan(ptr));
}

/* ============ Hyperbolic Functions ============ */

/** Hyperbolic sine, element-wise. */
export function sinh(x: NDArray): NDArray {
  return applyUnary(x, (ptr) => getWasmModule()._ufunc_sinh(ptr));
}

/** Hyperbolic cosine, element-wise. */
export function cosh(x: NDArray): NDArray {
  return applyUnary(x, (ptr) => getWasmModule()._ufunc_cosh(ptr));
}

/** Compute hyperbolic tangent element-wise. */
export function tanh(x: NDArray): NDArray {
  return applyUnary(x, (ptr) => getWasmModule()._ufunc_tanh(ptr));
}

/** Inverse hyperbolic sine element-wise. */
export function arcsinh(x: NDArray): NDArray {
  return applyUnary(x, (ptr) => getWasmModule()._ufunc_arcsinh(ptr));
}

/** Inverse hyperbolic cosine, element-wise. */
export function arccosh(x: NDArray): NDArray {
  return applyUnary(x, (ptr) => getWasmModule()._ufunc_arccosh(ptr));
}

/** Inverse hyperbolic tangent element-wise. */
export function arctanh(x: NDArray): NDArray {
  return applyUnary(x, (ptr) => getWasmModule()._ufunc_arctanh(ptr));
}

/* ============ Rounding Functions ============ */

/** Return the floor of the input, element-wise. */
export function floor(x: NDArray): NDArray {
  return applyUnary(x, (ptr) => getWasmModule()._ufunc_floor(ptr));
}

/** Return the ceiling of the input, element-wise. */
export function ceil(x: NDArray): NDArray {
  return applyUnary(x, (ptr) => getWasmModule()._ufunc_ceil(ptr));
}

/** Return the truncated value of the input, element-wise. */
export function trunc(x: NDArray): NDArray {
  return applyUnary(x, (ptr) => getWasmModule()._ufunc_trunc(ptr));
}

/** Round elements of the array to the nearest integer. */
export function rint(x: NDArray): NDArray {
  return applyUnary(x, (ptr) => getWasmModule()._ufunc_rint(ptr));
}

/** Round to nearest even value (banker's rounding). */
export function round(x: NDArray): NDArray {
  return applyUnary(x, (ptr) => getWasmModule()._ufunc_round(ptr));
}

/* ============ Angle Conversion ============ */

/** Convert angles from radians to degrees. */
export function degrees(x: NDArray): NDArray {
  return applyUnary(x, (ptr) => getWasmModule()._ufunc_degrees(ptr));
}

/** Alias for degrees. */
export const rad2deg = degrees;

/** Convert angles from degrees to radians. */
export function radians(x: NDArray): NDArray {
  return applyUnary(x, (ptr) => getWasmModule()._ufunc_radians(ptr));
}

/** Alias for radians. */
export const deg2rad = radians;

/* ============ Logical Operations (Unary) ============ */

/** Compute the truth value of NOT x element-wise. */
export function logical_not(x: NDArray): NDArray {
  return applyUnary(x, (ptr) => getWasmModule()._ufunc_logical_not(ptr));
}

/* ============ Bitwise Operations (Unary) ============ */

/** Compute bit-wise inversion, element-wise. */
export function invert(x: NDArray): NDArray {
  return applyUnary(x, (ptr) => getWasmModule()._ufunc_invert(ptr));
}

/** Alias for invert. */
export const bitwise_not = invert;

/* ============ Binary Arithmetic Operations ============ */

/** Add arguments element-wise. */
export function add(x1: NDArray, x2: NDArray): NDArray {
  return applyBinary(x1, x2, (p1, p2) => getWasmModule()._ufunc_add(p1, p2));
}

/** Subtract arguments, element-wise. */
export function subtract(x1: NDArray, x2: NDArray): NDArray {
  return applyBinary(x1, x2, (p1, p2) => getWasmModule()._ufunc_subtract(p1, p2));
}

/** Multiply arguments element-wise. */
export function multiply(x1: NDArray, x2: NDArray): NDArray {
  return applyBinary(x1, x2, (p1, p2) => getWasmModule()._ufunc_multiply(p1, p2));
}

/** Divide arguments element-wise (true division). */
export function divide(x1: NDArray, x2: NDArray): NDArray {
  return applyBinary(x1, x2, (p1, p2) => getWasmModule()._ufunc_divide(p1, p2));
}

/** Alias for divide. */
export const true_divide = divide;

/** Return the largest integer smaller or equal to the division of the inputs. */
export function floor_divide(x1: NDArray, x2: NDArray): NDArray {
  return applyBinary(x1, x2, (p1, p2) => getWasmModule()._ufunc_floor_divide(p1, p2));
}

/** Return element-wise remainder of division (Python-style modulo). */
export function remainder(x1: NDArray, x2: NDArray): NDArray {
  return applyBinary(x1, x2, (p1, p2) => getWasmModule()._ufunc_remainder(p1, p2));
}

/** Alias for remainder. */
export const mod = remainder;

/** Returns element-wise remainder of floor_divide (C-style modulo). */
export function fmod(x1: NDArray, x2: NDArray): NDArray {
  return applyBinary(x1, x2, (p1, p2) => getWasmModule()._ufunc_fmod(p1, p2));
}

/** First array elements raised to powers from second array, element-wise. */
export function power(x1: NDArray, x2: NDArray): NDArray {
  return applyBinary(x1, x2, (p1, p2) => getWasmModule()._ufunc_power(p1, p2));
}

/* ============ Comparison Operations ============ */

/** Return (x1 == x2) element-wise. */
export function equal(x1: NDArray, x2: NDArray): NDArray {
  return applyBinary(x1, x2, (p1, p2) => getWasmModule()._ufunc_equal(p1, p2));
}

/** Return (x1 != x2) element-wise. */
export function not_equal(x1: NDArray, x2: NDArray): NDArray {
  return applyBinary(x1, x2, (p1, p2) => getWasmModule()._ufunc_not_equal(p1, p2));
}

/** Return (x1 < x2) element-wise. */
export function less(x1: NDArray, x2: NDArray): NDArray {
  return applyBinary(x1, x2, (p1, p2) => getWasmModule()._ufunc_less(p1, p2));
}

/** Return (x1 <= x2) element-wise. */
export function less_equal(x1: NDArray, x2: NDArray): NDArray {
  return applyBinary(x1, x2, (p1, p2) => getWasmModule()._ufunc_less_equal(p1, p2));
}

/** Return (x1 > x2) element-wise. */
export function greater(x1: NDArray, x2: NDArray): NDArray {
  return applyBinary(x1, x2, (p1, p2) => getWasmModule()._ufunc_greater(p1, p2));
}

/** Return (x1 >= x2) element-wise. */
export function greater_equal(x1: NDArray, x2: NDArray): NDArray {
  return applyBinary(x1, x2, (p1, p2) => getWasmModule()._ufunc_greater_equal(p1, p2));
}

/* ============ Extrema Operations ============ */

/** Element-wise maximum of array elements (propagates NaN). */
export function maximum(x1: NDArray, x2: NDArray): NDArray {
  return applyBinary(x1, x2, (p1, p2) => getWasmModule()._ufunc_maximum(p1, p2));
}

/** Element-wise minimum of array elements (propagates NaN). */
export function minimum(x1: NDArray, x2: NDArray): NDArray {
  return applyBinary(x1, x2, (p1, p2) => getWasmModule()._ufunc_minimum(p1, p2));
}

/** Element-wise maximum (ignores NaN). */
export function fmax(x1: NDArray, x2: NDArray): NDArray {
  return applyBinary(x1, x2, (p1, p2) => getWasmModule()._ufunc_fmax(p1, p2));
}

/** Element-wise minimum (ignores NaN). */
export function fmin(x1: NDArray, x2: NDArray): NDArray {
  return applyBinary(x1, x2, (p1, p2) => getWasmModule()._ufunc_fmin(p1, p2));
}

/* ============ Logical Operations (Binary) ============ */

/** Compute the truth value of x1 AND x2 element-wise. */
export function logical_and(x1: NDArray, x2: NDArray): NDArray {
  return applyBinary(x1, x2, (p1, p2) => getWasmModule()._ufunc_logical_and(p1, p2));
}

/** Compute the truth value of x1 OR x2 element-wise. */
export function logical_or(x1: NDArray, x2: NDArray): NDArray {
  return applyBinary(x1, x2, (p1, p2) => getWasmModule()._ufunc_logical_or(p1, p2));
}

/** Compute the truth value of x1 XOR x2, element-wise. */
export function logical_xor(x1: NDArray, x2: NDArray): NDArray {
  return applyBinary(x1, x2, (p1, p2) => getWasmModule()._ufunc_logical_xor(p1, p2));
}

/* ============ Bitwise Operations (Binary) ============ */

/** Compute the bit-wise AND of two arrays element-wise. */
export function bitwise_and(x1: NDArray, x2: NDArray): NDArray {
  return applyBinary(x1, x2, (p1, p2) => getWasmModule()._ufunc_bitwise_and(p1, p2));
}

/** Compute the bit-wise OR of two arrays element-wise. */
export function bitwise_or(x1: NDArray, x2: NDArray): NDArray {
  return applyBinary(x1, x2, (p1, p2) => getWasmModule()._ufunc_bitwise_or(p1, p2));
}

/** Compute the bit-wise XOR of two arrays element-wise. */
export function bitwise_xor(x1: NDArray, x2: NDArray): NDArray {
  return applyBinary(x1, x2, (p1, p2) => getWasmModule()._ufunc_bitwise_xor(p1, p2));
}

/** Shift the bits of an integer to the left. */
export function left_shift(x1: NDArray, x2: NDArray): NDArray {
  return applyBinary(x1, x2, (p1, p2) => getWasmModule()._ufunc_left_shift(p1, p2));
}

/** Shift the bits of an integer to the right. */
export function right_shift(x1: NDArray, x2: NDArray): NDArray {
  return applyBinary(x1, x2, (p1, p2) => getWasmModule()._ufunc_right_shift(p1, p2));
}

/* ============ Special Mathematical Functions ============ */

/** Element-wise arc tangent of x1/x2 choosing the quadrant correctly. */
export function arctan2(x1: NDArray, x2: NDArray): NDArray {
  return applyBinary(x1, x2, (p1, p2) => getWasmModule()._ufunc_arctan2(p1, p2));
}

/** Given the "legs" of a right triangle, return its hypotenuse. */
export function hypot(x1: NDArray, x2: NDArray): NDArray {
  return applyBinary(x1, x2, (p1, p2) => getWasmModule()._ufunc_hypot(p1, p2));
}

/** Change the sign of x1 to that of x2, element-wise. */
export function copysign(x1: NDArray, x2: NDArray): NDArray {
  return applyBinary(x1, x2, (p1, p2) => getWasmModule()._ufunc_copysign(p1, p2));
}

/** Return element-wise True where signbit is set (less than zero). */
export function signbit(x: NDArray): NDArray {
  return applyUnary(x, (ptr) => getWasmModule()._ufunc_signbit(ptr));
}

/** Compute the logarithm of the sum of exponentials of the inputs. */
export function logaddexp(x1: NDArray, x2: NDArray): NDArray {
  return applyBinary(x1, x2, (p1, p2) => getWasmModule()._ufunc_logaddexp(p1, p2));
}

/** Compute log(2**x1 + 2**x2). */
export function logaddexp2(x1: NDArray, x2: NDArray): NDArray {
  return applyBinary(x1, x2, (p1, p2) => getWasmModule()._ufunc_logaddexp2(p1, p2));
}

/* ============ Phase 26: Miscellaneous Ufuncs ============ */

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
    throw new Error('frexp operation failed');
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
  return applyBinary(x1, x2, (p1, p2) => getWasmModule()._ufunc_nextafter(p1, p2));
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
    throw new Error('modf operation failed');
  }

  const fracPtr = module._ufunc_tuple_get_first(resultPtr);
  const intPtr = module._ufunc_tuple_get_second(resultPtr);
  module._ufunc_tuple_result_free(resultPtr);

  return [NDArray._fromPtr(fracPtr, module), NDArray._fromPtr(intPtr, module)];
}

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
    throw new TypeError('gcd requires integer inputs');
  }
  if (!isIntegerDType(x2.dtype) && x2.dtype !== DType.Bool) {
    throw new TypeError('gcd requires integer inputs');
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
    throw new TypeError('lcm requires integer inputs');
  }
  if (!isIntegerDType(x2.dtype) && x2.dtype !== DType.Bool) {
    throw new TypeError('lcm requires integer inputs');
  }
  return applyBinary(x1, x2, (p1, p2) => getWasmModule()._ufunc_lcm(p1, p2));
}

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
  return applyBinary(x1, x2, (p1, p2) => getWasmModule()._ufunc_heaviside(p1, p2));
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
    throw new Error('divmod operation failed');
  }

  const quotPtr = module._ufunc_tuple_get_first(resultPtr);
  const remPtr = module._ufunc_tuple_get_second(resultPtr);
  module._ufunc_tuple_result_free(resultPtr);

  return [NDArray._fromPtr(quotPtr, module), NDArray._fromPtr(remPtr, module)];
}

/**
 * Computes the number of 1-bits in the absolute value of x.
 * Also known as popcount or population count.
 * Input array must be integers.
 *
 * @param x - Input array of integers
 * @returns Number of 1-bits in each element (uint8)
 *
 * @example
 * bitwise_count(7)  // 3 (binary: 111)
 * bitwise_count([0, 1, 2, 3, 4, 5, 6, 7])
 * // [0, 1, 1, 2, 1, 2, 2, 3]
 */
export function bitwise_count(x: NDArray): NDArray {
  if (!isIntegerDType(x.dtype) && x.dtype !== DType.Bool) {
    throw new TypeError('bitwise_count requires integer inputs');
  }
  return applyUnary(x, (ptr) => getWasmModule()._ufunc_bitwise_count(ptr));
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
