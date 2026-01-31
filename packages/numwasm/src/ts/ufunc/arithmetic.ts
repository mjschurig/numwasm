/**
 * Arithmetic ufuncs for NumJS-WASM
 *
 * Basic arithmetic operations: add, subtract, multiply, divide, etc.
 */

import { NDArray } from "../_core/NDArray.js";
import { getWasmModule } from "../wasm-loader.js";
import { applyUnary, applyBinary } from "./helpers.js";

/* ============ Unary Arithmetic Operations ============ */

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

/** Compute the absolute values element-wise (always returns float). Same as absolute in numwasm. */
export const fabs = absolute;

/** Returns element-wise sign: -1, 0, or +1. */
export function sign(x: NDArray): NDArray {
  return applyUnary(x, (ptr) => getWasmModule()._ufunc_sign(ptr));
}

/* ============ Binary Arithmetic Operations ============ */

/** Add arguments element-wise. */
export function add(x1: NDArray, x2: NDArray): NDArray {
  return applyBinary(x1, x2, (p1, p2) => getWasmModule()._ufunc_add(p1, p2));
}

/** Subtract arguments, element-wise. */
export function subtract(x1: NDArray, x2: NDArray): NDArray {
  return applyBinary(x1, x2, (p1, p2) =>
    getWasmModule()._ufunc_subtract(p1, p2),
  );
}

/** Multiply arguments element-wise. */
export function multiply(x1: NDArray, x2: NDArray): NDArray {
  return applyBinary(x1, x2, (p1, p2) =>
    getWasmModule()._ufunc_multiply(p1, p2),
  );
}

/** Divide arguments element-wise (true division). */
export function divide(x1: NDArray, x2: NDArray): NDArray {
  return applyBinary(x1, x2, (p1, p2) => getWasmModule()._ufunc_divide(p1, p2));
}

/** Alias for divide. */
export const true_divide = divide;

/** Return the largest integer smaller or equal to the division of the inputs. */
export function floor_divide(x1: NDArray, x2: NDArray): NDArray {
  return applyBinary(x1, x2, (p1, p2) =>
    getWasmModule()._ufunc_floor_divide(p1, p2),
  );
}

/** Return element-wise remainder of division (Python-style modulo). */
export function remainder(x1: NDArray, x2: NDArray): NDArray {
  return applyBinary(x1, x2, (p1, p2) =>
    getWasmModule()._ufunc_remainder(p1, p2),
  );
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

/** Alias for power. */
export const pow = power;

/**
 * First array elements raised to powers from second array, element-wise.
 * Unlike power, float_power always returns float64 output.
 * In numwasm, this is equivalent to power since we use float64 by default.
 */
export function float_power(x1: NDArray, x2: NDArray): NDArray {
  // In NumPy, float_power promotes to float64. Our power already uses float64.
  return applyBinary(x1, x2, (p1, p2) => getWasmModule()._ufunc_power(p1, p2));
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
