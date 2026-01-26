/**
 * String comparison functions for NumJS.
 *
 * Provides vectorized string comparison operations on arrays of strings,
 * compatible with NumPy's numpy.strings module.
 *
 * Note: These functions require the WASM module to be loaded first.
 * Call loadWasmModule() before using these functions.
 */

import { NDArray } from '../NDArray.js';
import { DType } from '../types.js';
import { getWasmModule } from '../wasm-loader.js';

/**
 * Return (x1 == x2) element-wise for string arrays.
 *
 * @param x1 - First string array or JS string array
 * @param x2 - Second string array or JS string array
 * @returns Boolean array of comparison results
 *
 * @example
 * ```typescript
 * equal(['hello', 'world'], ['hello', 'test'])
 * // returns NDArray with [true, false] (as 1, 0)
 * ```
 */
export function equal(
  x1: NDArray | string[],
  x2: NDArray | string[]
): NDArray {
  const [a1, a2] = _prepareStringArrays(x1, x2);
  const result = _createBoolArray(a1.shape);

  for (let i = 0; i < a1.size; i++) {
    result.setFlat(i, a1.getStringFlat(i) === a2.getStringFlat(i) ? 1 : 0);
  }

  return result;
}

/**
 * Return (x1 != x2) element-wise for string arrays.
 *
 * @param x1 - First string array
 * @param x2 - Second string array
 * @returns Boolean array of comparison results
 */
export function not_equal(
  x1: NDArray | string[],
  x2: NDArray | string[]
): NDArray {
  const [a1, a2] = _prepareStringArrays(x1, x2);
  const result = _createBoolArray(a1.shape);

  for (let i = 0; i < a1.size; i++) {
    result.setFlat(i, a1.getStringFlat(i) !== a2.getStringFlat(i) ? 1 : 0);
  }

  return result;
}

/**
 * Return (x1 < x2) element-wise for string arrays.
 * Comparison is lexicographic (dictionary order).
 *
 * @param x1 - First string array
 * @param x2 - Second string array
 * @returns Boolean array of comparison results
 */
export function less(
  x1: NDArray | string[],
  x2: NDArray | string[]
): NDArray {
  const [a1, a2] = _prepareStringArrays(x1, x2);
  const result = _createBoolArray(a1.shape);

  for (let i = 0; i < a1.size; i++) {
    result.setFlat(i, a1.getStringFlat(i) < a2.getStringFlat(i) ? 1 : 0);
  }

  return result;
}

/**
 * Return (x1 <= x2) element-wise for string arrays.
 * Comparison is lexicographic.
 *
 * @param x1 - First string array
 * @param x2 - Second string array
 * @returns Boolean array of comparison results
 */
export function less_equal(
  x1: NDArray | string[],
  x2: NDArray | string[]
): NDArray {
  const [a1, a2] = _prepareStringArrays(x1, x2);
  const result = _createBoolArray(a1.shape);

  for (let i = 0; i < a1.size; i++) {
    result.setFlat(i, a1.getStringFlat(i) <= a2.getStringFlat(i) ? 1 : 0);
  }

  return result;
}

/**
 * Return (x1 > x2) element-wise for string arrays.
 * Comparison is lexicographic.
 *
 * @param x1 - First string array
 * @param x2 - Second string array
 * @returns Boolean array of comparison results
 */
export function greater(
  x1: NDArray | string[],
  x2: NDArray | string[]
): NDArray {
  const [a1, a2] = _prepareStringArrays(x1, x2);
  const result = _createBoolArray(a1.shape);

  for (let i = 0; i < a1.size; i++) {
    result.setFlat(i, a1.getStringFlat(i) > a2.getStringFlat(i) ? 1 : 0);
  }

  return result;
}

/**
 * Return (x1 >= x2) element-wise for string arrays.
 * Comparison is lexicographic.
 *
 * @param x1 - First string array
 * @param x2 - Second string array
 * @returns Boolean array of comparison results
 */
export function greater_equal(
  x1: NDArray | string[],
  x2: NDArray | string[]
): NDArray {
  const [a1, a2] = _prepareStringArrays(x1, x2);
  const result = _createBoolArray(a1.shape);

  for (let i = 0; i < a1.size; i++) {
    result.setFlat(i, a1.getStringFlat(i) >= a2.getStringFlat(i) ? 1 : 0);
  }

  return result;
}

/**
 * Compare two string arrays element-wise.
 *
 * @param x1 - First string array
 * @param x2 - Second string array
 * @returns Integer array: -1 if x1 < x2, 0 if equal, 1 if x1 > x2
 */
export function compare_chararrays(
  x1: NDArray | string[],
  x2: NDArray | string[]
): NDArray {
  const [a1, a2] = _prepareStringArrays(x1, x2);
  const result = _createInt32Array(a1.shape);

  for (let i = 0; i < a1.size; i++) {
    const s1 = a1.getStringFlat(i);
    const s2 = a2.getStringFlat(i);
    const cmp = s1 < s2 ? -1 : s1 > s2 ? 1 : 0;
    result.setFlat(i, cmp);
  }

  return result;
}

/* ============ Helper Functions ============ */

/**
 * Prepare string arrays for element-wise operations.
 * Converts JS arrays to NDArrays and validates shapes match.
 */
function _prepareStringArrays(
  x1: NDArray | string[],
  x2: NDArray | string[]
): [NDArray, NDArray] {
  const a1 = Array.isArray(x1) ? NDArray.fromStringArray(x1) : x1;
  const a2 = Array.isArray(x2) ? NDArray.fromStringArray(x2) : x2;

  // Validate both are string arrays
  if (!a1.isStringArray) {
    throw new Error('First argument must be a string array');
  }
  if (!a2.isStringArray) {
    throw new Error('Second argument must be a string array');
  }

  // Check shapes match (no broadcasting for now)
  if (a1.shape.length !== a2.shape.length) {
    throw new Error(
      `Shape mismatch: ${JSON.stringify(a1.shape)} vs ${JSON.stringify(a2.shape)}`
    );
  }
  for (let i = 0; i < a1.shape.length; i++) {
    if (a1.shape[i] !== a2.shape[i]) {
      throw new Error(
        `Shape mismatch: ${JSON.stringify(a1.shape)} vs ${JSON.stringify(a2.shape)}`
      );
    }
  }

  return [a1, a2];
}

/**
 * Create a boolean array synchronously using the pre-loaded WASM module.
 */
function _createBoolArray(shape: number[]): NDArray {
  const module = getWasmModule();

  // Allocate shape in WASM
  const shapePtr = module._malloc(shape.length * 4);
  for (let i = 0; i < shape.length; i++) {
    module.setValue(shapePtr + i * 4, shape[i], 'i32');
  }

  const ptr = module._ndarray_create(shape.length, shapePtr, DType.Bool);
  module._free(shapePtr);

  if (ptr === 0) {
    throw new Error('Failed to create boolean array');
  }

  return NDArray._fromPtr(ptr, module);
}

/**
 * Create an int32 array synchronously using the pre-loaded WASM module.
 */
function _createInt32Array(shape: number[]): NDArray {
  const module = getWasmModule();

  // Allocate shape in WASM
  const shapePtr = module._malloc(shape.length * 4);
  for (let i = 0; i < shape.length; i++) {
    module.setValue(shapePtr + i * 4, shape[i], 'i32');
  }

  const ptr = module._ndarray_create(shape.length, shapePtr, DType.Int32);
  module._free(shapePtr);

  if (ptr === 0) {
    throw new Error('Failed to create int32 array');
  }

  return NDArray._fromPtr(ptr, module);
}
