/**
 * Logical ufuncs for NumJS-WASM
 *
 * Element-wise logical operations.
 */

import { NDArray } from "../_core/NDArray.js";
import { getWasmModule } from "../wasm-loader.js";
import { applyUnary, applyBinary } from "./helpers.js";

/* ============ Logical Operations (Unary) ============ */

/** Compute the truth value of NOT x element-wise. */
export function logical_not(x: NDArray): NDArray {
  return applyUnary(x, (ptr) => getWasmModule()._ufunc_logical_not(ptr));
}

/* ============ Logical Operations (Binary) ============ */

/** Compute the truth value of x1 AND x2 element-wise. */
export function logical_and(x1: NDArray, x2: NDArray): NDArray {
  return applyBinary(x1, x2, (p1, p2) =>
    getWasmModule()._ufunc_logical_and(p1, p2),
  );
}

/** Compute the truth value of x1 OR x2 element-wise. */
export function logical_or(x1: NDArray, x2: NDArray): NDArray {
  return applyBinary(x1, x2, (p1, p2) =>
    getWasmModule()._ufunc_logical_or(p1, p2),
  );
}

/** Compute the truth value of x1 XOR x2, element-wise. */
export function logical_xor(x1: NDArray, x2: NDArray): NDArray {
  return applyBinary(x1, x2, (p1, p2) =>
    getWasmModule()._ufunc_logical_xor(p1, p2),
  );
}
