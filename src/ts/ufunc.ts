/**
 * NumJS Universal Functions (Ufuncs)
 *
 * Element-wise mathematical operations on NDArrays with broadcasting support.
 * These functions call optimized WASM implementations for performance.
 */

import { NDArray } from './NDArray.js';
import { getWasmModule } from './wasm-loader.js';

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
