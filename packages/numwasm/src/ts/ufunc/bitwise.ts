/**
 * Bitwise ufuncs for NumJS-WASM
 *
 * Element-wise bitwise operations.
 */

import { NDArray } from "../_core/NDArray.js";
import { DType } from "../types.js";
import { isIntegerDType } from "../dtype.js";
import { getWasmModule } from "../wasm-loader.js";
import { applyUnary, applyBinary } from "./helpers.js";

/* ============ Bitwise Operations (Unary) ============ */

/** Compute bit-wise inversion, element-wise. */
export function invert(x: NDArray): NDArray {
  return applyUnary(x, (ptr) => getWasmModule()._ufunc_invert(ptr));
}

/** Alias for invert. */
export const bitwise_not = invert;

/* ============ Bitwise Operations (Binary) ============ */

/** Compute the bit-wise AND of two arrays element-wise. */
export function bitwise_and(x1: NDArray, x2: NDArray): NDArray {
  return applyBinary(x1, x2, (p1, p2) =>
    getWasmModule()._ufunc_bitwise_and(p1, p2),
  );
}

/** Compute the bit-wise OR of two arrays element-wise. */
export function bitwise_or(x1: NDArray, x2: NDArray): NDArray {
  return applyBinary(x1, x2, (p1, p2) =>
    getWasmModule()._ufunc_bitwise_or(p1, p2),
  );
}

/** Compute the bit-wise XOR of two arrays element-wise. */
export function bitwise_xor(x1: NDArray, x2: NDArray): NDArray {
  return applyBinary(x1, x2, (p1, p2) =>
    getWasmModule()._ufunc_bitwise_xor(p1, p2),
  );
}

/** Shift the bits of an integer to the left. */
export function left_shift(x1: NDArray, x2: NDArray): NDArray {
  return applyBinary(x1, x2, (p1, p2) =>
    getWasmModule()._ufunc_left_shift(p1, p2),
  );
}

/** Shift the bits of an integer to the right. */
export function right_shift(x1: NDArray, x2: NDArray): NDArray {
  return applyBinary(x1, x2, (p1, p2) =>
    getWasmModule()._ufunc_right_shift(p1, p2),
  );
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
    throw new TypeError("bitwise_count requires integer inputs");
  }
  return applyUnary(x, (ptr) => getWasmModule()._ufunc_bitwise_count(ptr));
}
