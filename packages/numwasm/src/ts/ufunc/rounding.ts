/**
 * Rounding ufuncs for NumJS-WASM
 *
 * Element-wise rounding operations.
 */

import { NDArray } from "../_core/NDArray.js";
import { getWasmModule } from "../wasm-loader.js";
import { applyUnary } from "./helpers.js";

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

/** Alias for round. */
export const around = round;

/** Round to nearest integer towards zero (same as trunc). */
export const fix = trunc;
