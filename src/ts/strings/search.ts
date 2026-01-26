/**
 * String search and indexing functions for NumJS.
 *
 * Provides vectorized string search operations on arrays of strings,
 * compatible with NumPy's numpy.strings module.
 *
 * Note: These functions require the WASM module to be loaded first.
 * Call loadWasmModule() before using these functions.
 */

import { NDArray } from '../NDArray.js';
import { DType } from '../types.js';
import { getWasmModule } from '../wasm-loader.js';
import { ValueError } from './errors.js';

/**
 * For each element, return the lowest index where substring is found.
 *
 * @param a - Input string array
 * @param sub - Substring to search for
 * @param start - Start position for search (default 0)
 * @param end - End position for search (default string length)
 * @returns Integer array with indices, -1 where not found
 *
 * @example
 * ```typescript
 * find(['hello world', 'test'], 'o')
 * // returns [4, -1]
 * ```
 */
export function find(
  a: NDArray | string[],
  sub: string,
  start: number = 0,
  end: number | null = null
): NDArray {
  const arr = Array.isArray(a) ? NDArray.fromStringArray(a) : a;

  if (!arr.isStringArray) {
    throw new Error('Input must be a string array');
  }

  const result = _createInt32Array(arr.shape);

  for (let i = 0; i < arr.size; i++) {
    const s = arr.getStringFlat(i);
    const searchEnd = end ?? s.length;
    const searchStr = s.slice(start, searchEnd);
    const idx = searchStr.indexOf(sub);
    result.setFlat(i, idx === -1 ? -1 : idx + start);
  }

  return result;
}

/**
 * For each element, return the highest index where substring is found.
 *
 * @param a - Input string array
 * @param sub - Substring to search for
 * @param start - Start position for search (default 0)
 * @param end - End position for search (default string length)
 * @returns Integer array with indices, -1 where not found
 */
export function rfind(
  a: NDArray | string[],
  sub: string,
  start: number = 0,
  end: number | null = null
): NDArray {
  const arr = Array.isArray(a) ? NDArray.fromStringArray(a) : a;

  if (!arr.isStringArray) {
    throw new Error('Input must be a string array');
  }

  const result = _createInt32Array(arr.shape);

  for (let i = 0; i < arr.size; i++) {
    const s = arr.getStringFlat(i);
    const searchEnd = end ?? s.length;
    const searchStr = s.slice(start, searchEnd);
    const idx = searchStr.lastIndexOf(sub);
    result.setFlat(i, idx === -1 ? -1 : idx + start);
  }

  return result;
}

/**
 * Like find, but raises ValueError if substring is not found.
 *
 * @param a - Input string array
 * @param sub - Substring to search for
 * @param start - Start position for search (default 0)
 * @param end - End position for search (default string length)
 * @returns Integer array with indices
 * @throws {ValueError} If substring is not found in any element
 */
export function index(
  a: NDArray | string[],
  sub: string,
  start: number = 0,
  end: number | null = null
): NDArray {
  const result = find(a, sub, start, end);

  for (let i = 0; i < result.size; i++) {
    if (result.getFlat(i) === -1) {
      throw new ValueError('substring not found');
    }
  }

  return result;
}

/**
 * Like rfind, but raises ValueError if substring is not found.
 *
 * @param a - Input string array
 * @param sub - Substring to search for
 * @param start - Start position for search (default 0)
 * @param end - End position for search (default string length)
 * @returns Integer array with indices
 * @throws {ValueError} If substring is not found in any element
 */
export function rindex(
  a: NDArray | string[],
  sub: string,
  start: number = 0,
  end: number | null = null
): NDArray {
  const result = rfind(a, sub, start, end);

  for (let i = 0; i < result.size; i++) {
    if (result.getFlat(i) === -1) {
      throw new ValueError('substring not found');
    }
  }

  return result;
}

/**
 * Return the number of non-overlapping occurrences of substring.
 *
 * @param a - Input string array
 * @param sub - Substring to count
 * @param start - Start position for search (default 0)
 * @param end - End position for search (default string length)
 * @returns Integer array with occurrence counts
 *
 * @example
 * ```typescript
 * count(['aabababa', 'foo'], 'aba')
 * // returns [2, 0]
 * ```
 */
export function count(
  a: NDArray | string[],
  sub: string,
  start: number = 0,
  end: number | null = null
): NDArray {
  const arr = Array.isArray(a) ? NDArray.fromStringArray(a) : a;

  if (!arr.isStringArray) {
    throw new Error('Input must be a string array');
  }

  const result = _createInt32Array(arr.shape);

  for (let i = 0; i < arr.size; i++) {
    const s = arr.getStringFlat(i);
    const searchEnd = end ?? s.length;
    const searchStr = s.slice(start, searchEnd);

    if (sub.length === 0) {
      // Empty string matches between each character
      result.setFlat(i, searchStr.length + 1);
    } else {
      let cnt = 0;
      let pos = 0;
      while ((pos = searchStr.indexOf(sub, pos)) !== -1) {
        cnt++;
        pos += sub.length; // Non-overlapping
      }
      result.setFlat(i, cnt);
    }
  }

  return result;
}

/**
 * Return true for each element if the string starts with prefix.
 *
 * @param a - Input string array
 * @param prefix - Prefix to check
 * @param start - Start position (default 0)
 * @param end - End position (default string length)
 * @returns Boolean array
 *
 * @example
 * ```typescript
 * startswith(['hello', 'world', 'help'], 'hel')
 * // returns [true, false, true]
 * ```
 */
export function startswith(
  a: NDArray | string[],
  prefix: string,
  start: number = 0,
  end: number | null = null
): NDArray {
  const arr = Array.isArray(a) ? NDArray.fromStringArray(a) : a;

  if (!arr.isStringArray) {
    throw new Error('Input must be a string array');
  }

  const result = _createBoolArray(arr.shape);

  for (let i = 0; i < arr.size; i++) {
    const s = arr.getStringFlat(i);
    const searchEnd = end ?? s.length;
    const searchStr = s.slice(start, searchEnd);
    result.setFlat(i, searchStr.startsWith(prefix) ? 1 : 0);
  }

  return result;
}

/**
 * Return true for each element if the string ends with suffix.
 *
 * @param a - Input string array
 * @param suffix - Suffix to check
 * @param start - Start position (default 0)
 * @param end - End position (default string length)
 * @returns Boolean array
 *
 * @example
 * ```typescript
 * endswith(['hello', 'world', 'jello'], 'llo')
 * // returns [true, false, true]
 * ```
 */
export function endswith(
  a: NDArray | string[],
  suffix: string,
  start: number = 0,
  end: number | null = null
): NDArray {
  const arr = Array.isArray(a) ? NDArray.fromStringArray(a) : a;

  if (!arr.isStringArray) {
    throw new Error('Input must be a string array');
  }

  const result = _createBoolArray(arr.shape);

  for (let i = 0; i < arr.size; i++) {
    const s = arr.getStringFlat(i);
    const searchEnd = end ?? s.length;
    const searchStr = s.slice(start, searchEnd);
    result.setFlat(i, searchStr.endsWith(suffix) ? 1 : 0);
  }

  return result;
}

/* ============ Helper Functions ============ */

/**
 * Create a boolean array synchronously using the pre-loaded WASM module.
 */
function _createBoolArray(shape: number[]): NDArray {
  const module = getWasmModule();

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
