/**
 * String property testing functions for NumJS.
 *
 * Provides vectorized string property tests on arrays of strings,
 * compatible with NumPy's numpy.strings module.
 *
 * Note: These functions require the WASM module to be loaded first.
 * Call loadWasmModule() before using these functions.
 */

import { NDArray } from "../_core/NDArray.js";
import { DType } from "../types.js";
import { getWasmModule } from "../wasm-loader.js";

/**
 * Return true for each element if all characters are alphabetic.
 *
 * @param a - Input string array
 * @returns Boolean array
 *
 * @example
 * ```typescript
 * isalpha(['hello', '123', 'abc123'])
 * // returns [true, false, false]
 * ```
 */
export function isalpha(a: NDArray | string[]): NDArray {
  return _applyStringTest(a, (s) => s.length > 0 && /^[a-zA-Z]+$/.test(s));
}

/**
 * Return true for each element if all characters are digits (0-9).
 *
 * @param a - Input string array
 * @returns Boolean array
 */
export function isdigit(a: NDArray | string[]): NDArray {
  return _applyStringTest(a, (s) => s.length > 0 && /^[0-9]+$/.test(s));
}

/**
 * Return true for each element if all characters are alphanumeric.
 *
 * @param a - Input string array
 * @returns Boolean array
 */
export function isalnum(a: NDArray | string[]): NDArray {
  return _applyStringTest(a, (s) => s.length > 0 && /^[a-zA-Z0-9]+$/.test(s));
}

/**
 * Return true for each element if all characters are whitespace.
 *
 * @param a - Input string array
 * @returns Boolean array
 */
export function isspace(a: NDArray | string[]): NDArray {
  return _applyStringTest(a, (s) => s.length > 0 && /^\s+$/.test(s));
}

/**
 * Return true for each element if all cased characters are lowercase
 * and there is at least one cased character.
 *
 * @param a - Input string array
 * @returns Boolean array
 */
export function islower(a: NDArray | string[]): NDArray {
  return _applyStringTest(a, (s) => {
    if (!/[a-zA-Z]/.test(s)) return false; // No cased chars
    return s === s.toLowerCase();
  });
}

/**
 * Return true for each element if all cased characters are uppercase
 * and there is at least one cased character.
 *
 * @param a - Input string array
 * @returns Boolean array
 */
export function isupper(a: NDArray | string[]): NDArray {
  return _applyStringTest(a, (s) => {
    if (!/[a-zA-Z]/.test(s)) return false; // No cased chars
    return s === s.toUpperCase();
  });
}

/**
 * Return true for each element if the string is titlecased.
 * A string is titlecased if uppercase characters follow uncased
 * characters and lowercase characters follow cased characters.
 *
 * @param a - Input string array
 * @returns Boolean array
 */
export function istitle(a: NDArray | string[]): NDArray {
  return _applyStringTest(a, (s) => {
    if (s.length === 0) return false;

    // Check each word: first letter upper, rest lower
    const words = s.split(/(\s+)/);
    let hasCase = false;

    for (const word of words) {
      if (/^\s+$/.test(word)) continue;
      if (word.length === 0) continue;

      // Check if word follows title pattern
      let prevCased = false;
      for (const char of word) {
        const isUpper =
          char === char.toUpperCase() && char !== char.toLowerCase();
        const isLower =
          char === char.toLowerCase() && char !== char.toUpperCase();

        if (isUpper || isLower) hasCase = true;
        if (isUpper && prevCased) return false;
        if (isLower && !prevCased) return false;

        prevCased = isUpper || isLower;
      }
    }
    return hasCase;
  });
}

/**
 * Return true for each element if all characters are decimal.
 * Decimal characters are those in Unicode category "Nd".
 *
 * @param a - Input string array
 * @returns Boolean array
 */
export function isdecimal(a: NDArray | string[]): NDArray {
  return _applyStringTest(a, (s) => s.length > 0 && /^[0-9]+$/.test(s));
}

/**
 * Return true for each element if all characters are numeric.
 * Numeric characters include digit characters and characters
 * that have the Unicode numeric value property (subscripts, fractions, etc.).
 *
 * @param a - Input string array
 * @returns Boolean array
 */
export function isnumeric(a: NDArray | string[]): NDArray {
  return _applyStringTest(a, (s) => {
    if (s.length === 0) return false;

    // Extended numeric check: digits, superscripts, subscripts, fractions
    for (const char of s) {
      const code = char.charCodeAt(0);

      // Basic digits 0-9
      if (code >= 0x30 && code <= 0x39) continue;
      // Superscript digits
      if (code === 0xb2 || code === 0xb3 || code === 0xb9) continue;
      // Subscript digits
      if (code >= 0x2080 && code <= 0x2089) continue;
      // Vulgar fractions
      if (code >= 0xbc && code <= 0xbe) continue;
      // Fullwidth digits
      if (code >= 0xff10 && code <= 0xff19) continue;

      return false;
    }
    return true;
  });
}

/**
 * Return the length of each element.
 *
 * @param a - Input string array
 * @returns Integer array with string lengths
 *
 * @example
 * ```typescript
 * str_len(['hello', 'hi', ''])
 * // returns [5, 2, 0]
 * ```
 */
export function str_len(a: NDArray | string[]): NDArray {
  const arr = Array.isArray(a) ? NDArray.fromStringArray(a) : a;

  if (!arr.isStringArray) {
    throw new Error("Input must be a string array");
  }

  const result = _createInt32Array(arr.shape);

  for (let i = 0; i < arr.size; i++) {
    result.setFlat(i, arr.getStringFlat(i).length);
  }

  return result;
}

/* ============ Helper Functions ============ */

/**
 * Apply a boolean test to each string element.
 */
function _applyStringTest(
  a: NDArray | string[],
  test: (s: string) => boolean,
): NDArray {
  const arr = Array.isArray(a) ? NDArray.fromStringArray(a) : a;

  if (!arr.isStringArray) {
    throw new Error("Input must be a string array");
  }

  const result = _createBoolArray(arr.shape);

  for (let i = 0; i < arr.size; i++) {
    result.setFlat(i, test(arr.getStringFlat(i)) ? 1 : 0);
  }

  return result;
}

/**
 * Create a boolean array synchronously using the pre-loaded WASM module.
 */
function _createBoolArray(shape: number[]): NDArray {
  const module = getWasmModule();

  const shapePtr = module._malloc(shape.length * 4);
  for (let i = 0; i < shape.length; i++) {
    module.setValue(shapePtr + i * 4, shape[i], "i32");
  }

  const ptr = module._ndarray_create(shape.length, shapePtr, DType.Bool);
  module._free(shapePtr);

  if (ptr === 0) {
    throw new Error("Failed to create boolean array");
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
    module.setValue(shapePtr + i * 4, shape[i], "i32");
  }

  const ptr = module._ndarray_create(shape.length, shapePtr, DType.Int32);
  module._free(shapePtr);

  if (ptr === 0) {
    throw new Error("Failed to create int32 array");
  }

  return NDArray._fromPtr(ptr, module);
}
