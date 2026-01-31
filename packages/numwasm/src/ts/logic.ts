/**
 * NumJS Logic & Comparison Functions
 *
 * TypeScript wrappers for WASM logic and comparison operations.
 * Follows NumPy's logic and comparison API.
 *
 * Reference: numpy/_core/fromnumeric.py, numpy/_core/numeric.py
 */

import { NDArray } from "./_core/NDArray.js";
import { DType } from "./types.js";
import { absolute } from "./ufunc.js";

/* ============ Element-wise Predicates ============ */

/**
 * Test element-wise for finiteness (not infinity and not NaN).
 *
 * @param x - Input array
 * @returns Boolean array where True indicates finite value
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([1, Infinity, -Infinity, NaN]);
 * const result = await isfinite(arr);  // [true, false, false, false]
 * ```
 */
export async function isfinite(x: NDArray): Promise<NDArray> {
  const module = x._wasmModule;
  const resultPtr = module._ndarray_isfinite(x._wasmPtr);
  if (resultPtr === 0) {
    throw new Error("isfinite failed");
  }
  return NDArray._fromPtr(resultPtr, module);
}

/**
 * Test element-wise for positive or negative infinity.
 *
 * @param x - Input array
 * @returns Boolean array where True indicates infinity
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([1, Infinity, -Infinity, NaN]);
 * const result = await isinf(arr);  // [false, true, true, false]
 * ```
 */
export async function isinf(x: NDArray): Promise<NDArray> {
  const module = x._wasmModule;
  const resultPtr = module._ndarray_isinf(x._wasmPtr);
  if (resultPtr === 0) {
    throw new Error("isinf failed");
  }
  return NDArray._fromPtr(resultPtr, module);
}

/**
 * Test element-wise for NaN.
 *
 * @param x - Input array
 * @returns Boolean array where True indicates NaN
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([1, Infinity, NaN]);
 * const result = await isnan(arr);  // [false, false, true]
 * ```
 */
export async function isnan(x: NDArray): Promise<NDArray> {
  const module = x._wasmModule;
  const resultPtr = module._ndarray_isnan(x._wasmPtr);
  if (resultPtr === 0) {
    throw new Error("isnan failed");
  }
  return NDArray._fromPtr(resultPtr, module);
}

/**
 * Test element-wise for negative infinity.
 *
 * @param x - Input array
 * @returns Boolean array where True indicates negative infinity
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([1, Infinity, -Infinity, NaN]);
 * const result = await isneginf(arr);  // [false, false, true, false]
 * ```
 */
export async function isneginf(x: NDArray): Promise<NDArray> {
  const module = x._wasmModule;
  const resultPtr = module._ndarray_isneginf(x._wasmPtr);
  if (resultPtr === 0) {
    throw new Error("isneginf failed");
  }
  return NDArray._fromPtr(resultPtr, module);
}

/**
 * Test element-wise for positive infinity.
 *
 * @param x - Input array
 * @returns Boolean array where True indicates positive infinity
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([1, Infinity, -Infinity, NaN]);
 * const result = await isposinf(arr);  // [false, true, false, false]
 * ```
 */
export async function isposinf(x: NDArray): Promise<NDArray> {
  const module = x._wasmModule;
  const resultPtr = module._ndarray_isposinf(x._wasmPtr);
  if (resultPtr === 0) {
    throw new Error("isposinf failed");
  }
  return NDArray._fromPtr(resultPtr, module);
}

/**
 * Returns a bool array where True if input element is complex (imaginary part != 0).
 * For non-complex dtype arrays, returns all False.
 *
 * @param x - Input array
 * @returns Boolean array
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([1+2j, 3+0j], { dtype: DType.Complex128 });
 * const result = await iscomplex(arr);  // [true, false]
 * ```
 */
export async function iscomplex(x: NDArray): Promise<NDArray> {
  const module = x._wasmModule;
  const resultPtr = module._ndarray_iscomplex_elem(x._wasmPtr);
  if (resultPtr === 0) {
    throw new Error("iscomplex failed");
  }
  return NDArray._fromPtr(resultPtr, module);
}

/**
 * Returns a bool array where True if input element is real (imaginary part == 0).
 * For non-complex dtype arrays, returns all True.
 *
 * @param x - Input array
 * @returns Boolean array
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([1+2j, 3+0j], { dtype: DType.Complex128 });
 * const result = await isreal(arr);  // [false, true]
 * ```
 */
export async function isreal(x: NDArray): Promise<NDArray> {
  const module = x._wasmModule;
  const resultPtr = module._ndarray_isreal_elem(x._wasmPtr);
  if (resultPtr === 0) {
    throw new Error("isreal failed");
  }
  return NDArray._fromPtr(resultPtr, module);
}

/* ============ Type Checking (TypeScript only) ============ */

/**
 * Check if the array has a complex dtype.
 *
 * @param x - Input array
 * @returns True if dtype is Complex64 or Complex128
 *
 * @example
 * ```typescript
 * const arr = await NDArray.zeros([2, 2], { dtype: DType.Complex128 });
 * iscomplexobj(arr);  // true
 * ```
 */
export function iscomplexobj(x: NDArray): boolean {
  return x.dtype === DType.Complex64 || x.dtype === DType.Complex128;
}

/**
 * Check if the array has a non-complex dtype.
 *
 * @param x - Input array
 * @returns True if dtype is not complex
 *
 * @example
 * ```typescript
 * const arr = await NDArray.zeros([2, 2]);
 * isrealobj(arr);  // true
 * ```
 */
export function isrealobj(x: NDArray): boolean {
  return !iscomplexobj(x);
}

/**
 * Check if the array is Fortran contiguous but not C contiguous.
 *
 * @param a - Input array
 * @returns True if F-contiguous and not C-contiguous
 *
 * @example
 * ```typescript
 * const arr = await NDArray.zeros([3, 4]);
 * const fortran = arr.asfortranarray();
 * isfortran(fortran);  // true
 * ```
 */
export function isfortran(a: NDArray): boolean {
  const flags = a.flags;
  return flags.f_contiguous && !flags.c_contiguous;
}

/**
 * Check if an element is a scalar type.
 *
 * @param element - Value to check
 * @returns True if element is number, boolean, or bigint
 *
 * @example
 * ```typescript
 * isscalar(5);        // true
 * isscalar(3.14);     // true
 * isscalar(true);     // true
 * isscalar([1, 2]);   // false
 * ```
 */
export function isscalar(element: unknown): boolean {
  return (
    typeof element === "number" ||
    typeof element === "boolean" ||
    typeof element === "bigint"
  );
}

/* ============ Truth Testing (Reductions) ============ */

/**
 * Test whether all array elements along a given axis evaluate to True.
 *
 * @param a - Input array
 * @param axis - Axis along which to test (undefined for entire array)
 * @param keepdims - If true, keep reduced axis as size-1 dimension
 * @returns Boolean result (scalar or array)
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([[true, true], [true, false]]);
 * await all(arr);           // false
 * await all(arr, 0);        // [true, false]
 * await all(arr, 1);        // [true, false]
 * ```
 */
export async function all(
  a: NDArray,
  axis?: number,
  keepdims: boolean = false,
): Promise<NDArray | boolean> {
  const module = a._wasmModule;

  let resultPtr: number;

  if (axis === undefined) {
    resultPtr = module._ndarray_all(a._wasmPtr);
  } else {
    resultPtr = module._ndarray_all_axis(a._wasmPtr, axis, keepdims ? 1 : 0);
  }

  if (resultPtr === 0) {
    throw new Error("all failed");
  }

  const result = NDArray._fromPtr(resultPtr, module);

  // For scalar result without keepdims, return boolean
  if (axis === undefined || result.ndim === 0) {
    const val = result.item() !== 0;
    result.dispose();
    return val;
  }

  return result;
}

/**
 * Test whether any array element along a given axis evaluates to True.
 *
 * @param a - Input array
 * @param axis - Axis along which to test (undefined for entire array)
 * @param keepdims - If true, keep reduced axis as size-1 dimension
 * @returns Boolean result (scalar or array)
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([[false, false], [false, true]]);
 * await any(arr);           // true
 * await any(arr, 0);        // [false, true]
 * await any(arr, 1);        // [false, true]
 * ```
 */
export async function any(
  a: NDArray,
  axis?: number,
  keepdims: boolean = false,
): Promise<NDArray | boolean> {
  const module = a._wasmModule;

  let resultPtr: number;

  if (axis === undefined) {
    resultPtr = module._ndarray_any(a._wasmPtr);
  } else {
    resultPtr = module._ndarray_any_axis(a._wasmPtr, axis, keepdims ? 1 : 0);
  }

  if (resultPtr === 0) {
    throw new Error("any failed");
  }

  const result = NDArray._fromPtr(resultPtr, module);

  // For scalar result without keepdims, return boolean
  if (axis === undefined || result.ndim === 0) {
    const val = result.item() !== 0;
    result.dispose();
    return val;
  }

  return result;
}

/* ============ Comparison Functions ============ */

/**
 * Returns a boolean array where two arrays are element-wise equal within tolerance.
 *
 * Formula: |a - b| <= (atol + rtol * |b|)
 *
 * @param a - First array
 * @param b - Second array (reference)
 * @param rtol - Relative tolerance (default 1e-5)
 * @param atol - Absolute tolerance (default 1e-8)
 * @param equal_nan - Whether to compare NaN's as equal (default false)
 * @returns Boolean array of element-wise close status
 *
 * @example
 * ```typescript
 * const a = await NDArray.fromArray([1e10, 1e-7]);
 * const b = await NDArray.fromArray([1.00001e10, 1e-8]);
 * const result = await isclose(a, b);  // [true, false]
 * ```
 */
export async function isclose(
  a: NDArray,
  b: NDArray,
  rtol: number = 1e-5,
  atol: number = 1e-8,
  equal_nan: boolean = false,
): Promise<NDArray> {
  const module = a._wasmModule;

  const resultPtr = module._ndarray_isclose(
    a._wasmPtr,
    b._wasmPtr,
    rtol,
    atol,
    equal_nan ? 1 : 0,
  );

  if (resultPtr === 0) {
    throw new Error("isclose failed: shapes may not be broadcast-compatible");
  }

  return NDArray._fromPtr(resultPtr, module);
}

/**
 * Returns True if two arrays are element-wise equal within a tolerance.
 *
 * @param a - First array
 * @param b - Second array
 * @param rtol - Relative tolerance (default 1e-5)
 * @param atol - Absolute tolerance (default 1e-8)
 * @param equal_nan - Whether to compare NaN's as equal (default false)
 * @returns True if all elements are close
 *
 * @example
 * ```typescript
 * const a = await NDArray.fromArray([1e10, 1e-8]);
 * const b = await NDArray.fromArray([1.00001e10, 1e-9]);
 * await allclose(a, b);  // true
 * ```
 */
export async function allclose(
  a: NDArray,
  b: NDArray,
  rtol: number = 1e-5,
  atol: number = 1e-8,
  equal_nan: boolean = false,
): Promise<boolean> {
  const module = a._wasmModule;

  const result = module._ndarray_allclose(
    a._wasmPtr,
    b._wasmPtr,
    rtol,
    atol,
    equal_nan ? 1 : 0,
  );

  return result !== 0;
}

/**
 * True if two arrays have the same shape and elements, False otherwise.
 *
 * @param a1 - First array
 * @param a2 - Second array
 * @param equal_nan - Whether to compare NaN's as equal (default false)
 * @returns True if arrays are equal
 *
 * @example
 * ```typescript
 * const a = await NDArray.fromArray([1, 2, 3]);
 * const b = await NDArray.fromArray([1, 2, 3]);
 * await array_equal(a, b);  // true
 * ```
 */
export async function array_equal(
  a1: NDArray,
  a2: NDArray,
  equal_nan: boolean = false,
): Promise<boolean> {
  const module = a1._wasmModule;

  const result = module._ndarray_array_equal(
    a1._wasmPtr,
    a2._wasmPtr,
    equal_nan ? 1 : 0,
  );

  return result !== 0;
}

/**
 * Returns True if input arrays are shape consistent and all elements equal.
 *
 * Shape consistent means they are either the same shape, or one input array
 * can be broadcast to create the same shape as the other one.
 *
 * @param a1 - First array
 * @param a2 - Second array
 * @returns True if arrays are equivalent
 *
 * @example
 * ```typescript
 * const a = await NDArray.fromArray([1, 2]);
 * const b = await NDArray.fromArray([[1, 2], [1, 2]]);
 * await array_equiv(a, b);  // true
 * ```
 */
export async function array_equiv(a1: NDArray, a2: NDArray): Promise<boolean> {
  const module = a1._wasmModule;

  const result = module._ndarray_array_equiv(a1._wasmPtr, a2._wasmPtr);

  return result !== 0;
}

/**
 * If input is complex with all imaginary parts close to zero, return real parts.
 *
 * "Close to zero" is defined as `tol` times the machine epsilon of the type.
 *
 * @param a - Input array
 * @param tol - Tolerance in machine epsilons for complex128, or absolute
 *              tolerance if tol < 1. Default is 100.
 * @returns Real array if imaginary parts are negligible, otherwise the input
 *
 * @example
 * ```typescript
 * // Imaginary parts within tolerance - returns real
 * const c1 = await NDArray.fromArray([1+1e-16j, 2+1e-16j], { dtype: DType.Complex128 });
 * const r1 = await real_if_close(c1);  // [1, 2] (Float64)
 *
 * // Imaginary parts too large - returns complex unchanged
 * const c2 = await NDArray.fromArray([1+1j, 2+2j], { dtype: DType.Complex128 });
 * const r2 = await real_if_close(c2);  // [1+1j, 2+2j] (Complex128)
 * ```
 */
export async function real_if_close(
  a: NDArray,
  tol: number = 100,
): Promise<NDArray> {
  // For non-complex arrays, return as-is
  if (!iscomplexobj(a)) {
    return a;
  }

  // Get machine epsilon for this type
  // Float32 epsilon: ~1.19e-7, Float64 epsilon: ~2.22e-16
  const eps =
    a.dtype === DType.Complex64 ? 1.1920929e-7 : 2.220446049250313e-16;
  const tolerance = tol > 1 ? eps * tol : tol;

  // Check if all imaginary parts are below tolerance
  const imagPart = a.imag;
  const absImag = absolute(imagPart);

  // Find maximum absolute imaginary value
  const module = a._wasmModule;
  const size = absImag.size;
  let maxImag = 0;

  for (let i = 0; i < size; i++) {
    const val = module._ndarray_get_flat(absImag._wasmPtr, i);
    if (val > maxImag) {
      maxImag = val;
    }
  }

  if (maxImag < tolerance) {
    // All imaginary parts are negligible, return real parts
    return a.real;
  }

  // Imaginary parts are significant, return as-is
  return a;
}
