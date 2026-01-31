/**
 * Trigonometric ufuncs for NumJS-WASM
 *
 * Trigonometric and hyperbolic functions.
 */

import { NDArray } from "../_core/NDArray.js";
import { getWasmModule } from "../wasm-loader.js";
import { applyUnary, applyBinary } from "./helpers.js";

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

/** NumPy 2.0 alias for arcsin. */
export const asin = arcsin;

/** NumPy 2.0 alias for arccos. */
export const acos = arccos;

/** NumPy 2.0 alias for arctan. */
export const atan = arctan;

/** Element-wise arc tangent of x1/x2 choosing the quadrant correctly. */
export function arctan2(x1: NDArray, x2: NDArray): NDArray {
  return applyBinary(x1, x2, (p1, p2) =>
    getWasmModule()._ufunc_arctan2(p1, p2),
  );
}

/** NumPy 2.0 alias for arctan2. */
export const atan2 = arctan2;

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

/** NumPy 2.0 alias for arcsinh. */
export const asinh = arcsinh;

/** NumPy 2.0 alias for arccosh. */
export const acosh = arccosh;

/** NumPy 2.0 alias for arctanh. */
export const atanh = arctanh;

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

/** Given the "legs" of a right triangle, return its hypotenuse. */
export function hypot(x1: NDArray, x2: NDArray): NDArray {
  return applyBinary(x1, x2, (p1, p2) => getWasmModule()._ufunc_hypot(p1, p2));
}
