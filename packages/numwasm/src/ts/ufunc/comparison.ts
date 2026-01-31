/**
 * Comparison ufuncs for NumJS-WASM
 *
 * Element-wise comparison operations.
 */

import { NDArray } from "../_core/NDArray.js";
import { getWasmModule } from "../wasm-loader.js";
import { applyBinary } from "./helpers.js";

/* ============ Comparison Operations ============ */

/** Return (x1 == x2) element-wise. */
export function equal(x1: NDArray, x2: NDArray): NDArray {
  return applyBinary(x1, x2, (p1, p2) => getWasmModule()._ufunc_equal(p1, p2));
}

/** Return (x1 != x2) element-wise. */
export function not_equal(x1: NDArray, x2: NDArray): NDArray {
  return applyBinary(x1, x2, (p1, p2) =>
    getWasmModule()._ufunc_not_equal(p1, p2),
  );
}

/** Return (x1 < x2) element-wise. */
export function less(x1: NDArray, x2: NDArray): NDArray {
  return applyBinary(x1, x2, (p1, p2) => getWasmModule()._ufunc_less(p1, p2));
}

/** Return (x1 <= x2) element-wise. */
export function less_equal(x1: NDArray, x2: NDArray): NDArray {
  return applyBinary(x1, x2, (p1, p2) =>
    getWasmModule()._ufunc_less_equal(p1, p2),
  );
}

/** Return (x1 > x2) element-wise. */
export function greater(x1: NDArray, x2: NDArray): NDArray {
  return applyBinary(x1, x2, (p1, p2) =>
    getWasmModule()._ufunc_greater(p1, p2),
  );
}

/** Return (x1 >= x2) element-wise. */
export function greater_equal(x1: NDArray, x2: NDArray): NDArray {
  return applyBinary(x1, x2, (p1, p2) =>
    getWasmModule()._ufunc_greater_equal(p1, p2),
  );
}

/* ============ Extrema Operations ============ */

/** Element-wise maximum of array elements (propagates NaN). */
export function maximum(x1: NDArray, x2: NDArray): NDArray {
  return applyBinary(x1, x2, (p1, p2) =>
    getWasmModule()._ufunc_maximum(p1, p2),
  );
}

/** Element-wise minimum of array elements (propagates NaN). */
export function minimum(x1: NDArray, x2: NDArray): NDArray {
  return applyBinary(x1, x2, (p1, p2) =>
    getWasmModule()._ufunc_minimum(p1, p2),
  );
}

/** Element-wise maximum (ignores NaN). */
export function fmax(x1: NDArray, x2: NDArray): NDArray {
  return applyBinary(x1, x2, (p1, p2) => getWasmModule()._ufunc_fmax(p1, p2));
}

/** Element-wise minimum (ignores NaN). */
export function fmin(x1: NDArray, x2: NDArray): NDArray {
  return applyBinary(x1, x2, (p1, p2) => getWasmModule()._ufunc_fmin(p1, p2));
}
