/**
 * Exponential and logarithmic ufuncs for NumJS-WASM
 *
 * Exponential, logarithmic, and related functions.
 */

import { NDArray } from "../_core/NDArray.js";
import { getWasmModule } from "../wasm-loader.js";
import { applyUnary, applyBinary } from "./helpers.js";

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

/** Compute the logarithm of the sum of exponentials of the inputs. */
export function logaddexp(x1: NDArray, x2: NDArray): NDArray {
  return applyBinary(x1, x2, (p1, p2) =>
    getWasmModule()._ufunc_logaddexp(p1, p2),
  );
}

/** Compute log(2**x1 + 2**x2). */
export function logaddexp2(x1: NDArray, x2: NDArray): NDArray {
  return applyBinary(x1, x2, (p1, p2) =>
    getWasmModule()._ufunc_logaddexp2(p1, p2),
  );
}
