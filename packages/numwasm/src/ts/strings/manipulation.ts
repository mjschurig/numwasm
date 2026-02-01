/**
 * String manipulation functions for NumJS.
 *
 * Provides vectorized string manipulation operations on arrays of strings,
 * compatible with NumPy's numpy.strings module.
 *
 * Note: These functions require the WASM module to be loaded first.
 * Call loadWasmModule() before using these functions.
 */

import { NDArray } from "../_core/NDArray.js";

/* ============ Case Conversion ============ */

/**
 * Return element-wise copy with lowercase characters converted to uppercase.
 *
 * @param a - Input string array
 * @returns String array with uppercase strings
 */
export function upper(a: NDArray | string[]): NDArray {
  return _applyTransform(a, (s) => s.toUpperCase());
}

/**
 * Return element-wise copy with uppercase characters converted to lowercase.
 *
 * @param a - Input string array
 * @returns String array with lowercase strings
 */
export function lower(a: NDArray | string[]): NDArray {
  return _applyTransform(a, (s) => s.toLowerCase());
}

/**
 * Return element-wise copy with uppercase and lowercase swapped.
 *
 * @param a - Input string array
 * @returns String array with swapped case
 */
export function swapcase(a: NDArray | string[]): NDArray {
  return _applyTransform(a, (s) => {
    return s
      .split("")
      .map((c) => {
        if (c === c.toLowerCase() && c !== c.toUpperCase()) {
          return c.toUpperCase();
        }
        if (c === c.toUpperCase() && c !== c.toLowerCase()) {
          return c.toLowerCase();
        }
        return c;
      })
      .join("");
  });
}

/**
 * Return element-wise copy with first character capitalized and rest lowercase.
 *
 * @param a - Input string array
 * @returns String array with capitalized strings
 */
export function capitalize(a: NDArray | string[]): NDArray {
  return _applyTransform(a, (s) => {
    if (s.length === 0) return s;
    return s[0].toUpperCase() + s.slice(1).toLowerCase();
  });
}

/**
 * Return element-wise titlecased version.
 * Words start with uppercase, remaining characters are lowercase.
 *
 * @param a - Input string array
 * @returns String array with titlecased strings
 */
export function title(a: NDArray | string[]): NDArray {
  return _applyTransform(a, (s) => {
    return s
      .replace(/\b\w/g, (c) => c.toUpperCase())
      .replace(/\B\w/g, (c) => c.toLowerCase());
  });
}

/* ============ Concatenation ============ */

/**
 * Return element-wise string concatenation.
 *
 * @param x1 - First string array
 * @param x2 - Second string array
 * @returns String array with concatenated strings
 *
 * @example
 * ```typescript
 * add(['hello', 'good'], [' world', 'bye'])
 * // returns ['hello world', 'goodbye']
 * ```
 */
export function add(x1: NDArray | string[], x2: NDArray | string[]): NDArray {
  const a1 = Array.isArray(x1) ? NDArray.fromStringArray(x1) : x1;
  const a2 = Array.isArray(x2) ? NDArray.fromStringArray(x2) : x2;

  if (!a1.isStringArray || !a2.isStringArray) {
    throw new Error("Both arguments must be string arrays");
  }

  // Check shapes match
  _validateShapes(a1, a2);

  const result = NDArray.emptyString(a1.shape);

  for (let i = 0; i < a1.size; i++) {
    result.setStringFlat(i, a1.getStringFlat(i) + a2.getStringFlat(i));
  }

  return result;
}

/**
 * Return element-wise (a * i), string repeated i times.
 *
 * @param a - Input string array
 * @param i - Number of repetitions (or array of repetition counts)
 * @returns String array with repeated strings
 *
 * @example
 * ```typescript
 * multiply(['ab', 'cd'], 3)
 * // returns ['ababab', 'cdcdcd']
 * ```
 */
export function multiply(a: NDArray | string[], i: number | number[]): NDArray {
  const arr = Array.isArray(a) ? NDArray.fromStringArray(a) : a;

  if (!arr.isStringArray) {
    throw new Error("First argument must be a string array");
  }

  const result = NDArray.emptyString(arr.shape);
  const counts = Array.isArray(i) ? i : null;

  for (let j = 0; j < arr.size; j++) {
    const s = arr.getStringFlat(j);
    const count = counts ? counts[j % counts.length] : (i as number);
    result.setStringFlat(j, count > 0 ? s.repeat(count) : "");
  }

  return result;
}

/* ============ Whitespace ============ */

/**
 * Return element-wise copy with leading and trailing characters removed.
 *
 * @param a - Input string array
 * @param chars - Characters to remove (default: whitespace)
 * @returns String array with stripped strings
 */
export function strip(
  a: NDArray | string[],
  chars: string | null = null,
): NDArray {
  return _applyTransform(a, (s) => {
    if (chars === null) return s.trim();
    const regex = new RegExp(
      `^[${_escapeRegex(chars)}]+|[${_escapeRegex(chars)}]+$`,
      "g",
    );
    return s.replace(regex, "");
  });
}

/**
 * Return element-wise copy with leading characters removed.
 *
 * @param a - Input string array
 * @param chars - Characters to remove (default: whitespace)
 * @returns String array with left-stripped strings
 */
export function lstrip(
  a: NDArray | string[],
  chars: string | null = null,
): NDArray {
  return _applyTransform(a, (s) => {
    if (chars === null) return s.trimStart();
    const regex = new RegExp(`^[${_escapeRegex(chars)}]+`);
    return s.replace(regex, "");
  });
}

/**
 * Return element-wise copy with trailing characters removed.
 *
 * @param a - Input string array
 * @param chars - Characters to remove (default: whitespace)
 * @returns String array with right-stripped strings
 */
export function rstrip(
  a: NDArray | string[],
  chars: string | null = null,
): NDArray {
  return _applyTransform(a, (s) => {
    if (chars === null) return s.trimEnd();
    const regex = new RegExp(`[${_escapeRegex(chars)}]+$`);
    return s.replace(regex, "");
  });
}

/**
 * Return element-wise copy with tabs expanded to spaces.
 *
 * @param a - Input string array
 * @param tabsize - Tab stop spacing (default 8)
 * @returns String array with expanded tabs
 */
export function expandtabs(
  a: NDArray | string[],
  tabsize: number = 8,
): NDArray {
  return _applyTransform(a, (s) => {
    let result = "";
    let col = 0;

    for (const char of s) {
      if (char === "\t") {
        const spaces = tabsize - (col % tabsize);
        result += " ".repeat(spaces);
        col += spaces;
      } else if (char === "\n" || char === "\r") {
        result += char;
        col = 0;
      } else {
        result += char;
        col++;
      }
    }
    return result;
  });
}

/* ============ Replacement ============ */

/**
 * Return element-wise copy with occurrences of old replaced by new.
 *
 * @param a - Input string array
 * @param old - Substring to replace
 * @param new_ - Replacement string
 * @param count - Maximum replacements per element (-1 for all)
 * @returns String array with replaced strings
 */
export function replace(
  a: NDArray | string[],
  old: string,
  new_: string,
  count: number = -1,
): NDArray {
  return _applyTransform(a, (s) => {
    if (count === 0) return s;
    if (count < 0) return s.split(old).join(new_);

    let result = s;
    let remaining = count;
    while (remaining > 0 && result.includes(old)) {
      result = result.replace(old, new_);
      remaining--;
    }
    return result;
  });
}

/* ============ Alignment ============ */

/**
 * Return element-wise centered string of given width.
 *
 * @param a - Input string array
 * @param width - Total width of output string
 * @param fillchar - Padding character (default space)
 * @returns String array with centered strings
 */
export function center(
  a: NDArray | string[],
  width: number,
  fillchar: string = " ",
): NDArray {
  if (fillchar.length !== 1) {
    throw new TypeError("fillchar must be exactly one character");
  }

  return _applyTransform(a, (s) => {
    if (s.length >= width) return s;
    const totalPad = width - s.length;
    const leftPad = Math.floor(totalPad / 2);
    const rightPad = totalPad - leftPad;
    return fillchar.repeat(leftPad) + s + fillchar.repeat(rightPad);
  });
}

/**
 * Return element-wise left-justified string of given width.
 *
 * @param a - Input string array
 * @param width - Total width of output string
 * @param fillchar - Padding character (default space)
 * @returns String array with left-justified strings
 */
export function ljust(
  a: NDArray | string[],
  width: number,
  fillchar: string = " ",
): NDArray {
  if (fillchar.length !== 1) {
    throw new TypeError("fillchar must be exactly one character");
  }

  return _applyTransform(a, (s) => {
    if (s.length >= width) return s;
    return s + fillchar.repeat(width - s.length);
  });
}

/**
 * Return element-wise right-justified string of given width.
 *
 * @param a - Input string array
 * @param width - Total width of output string
 * @param fillchar - Padding character (default space)
 * @returns String array with right-justified strings
 */
export function rjust(
  a: NDArray | string[],
  width: number,
  fillchar: string = " ",
): NDArray {
  if (fillchar.length !== 1) {
    throw new TypeError("fillchar must be exactly one character");
  }

  return _applyTransform(a, (s) => {
    if (s.length >= width) return s;
    return fillchar.repeat(width - s.length) + s;
  });
}

/**
 * Return element-wise numeric string left-padded with zeros.
 * Sign prefix is handled correctly.
 *
 * @param a - Input string array
 * @param width - Total width of output string
 * @returns String array with zero-filled strings
 *
 * @example
 * ```typescript
 * zfill(['-42', '42'], 5)
 * // returns ['-0042', '00042']
 * ```
 */
export function zfill(a: NDArray | string[], width: number): NDArray {
  return _applyTransform(a, (s) => {
    if (s.length >= width) return s;

    const sign = s[0] === "+" || s[0] === "-" ? s[0] : "";
    const rest = sign ? s.slice(1) : s;
    const padLength = width - s.length;

    return sign + "0".repeat(padLength) + rest;
  });
}

/* ============ Partitioning ============ */

/**
 * Partition each element around the first occurrence of sep.
 *
 * @param a - Input string array
 * @param sep - Separator string
 * @returns Array with shape (..., 3) containing (before, sep, after)
 *
 * @example
 * ```typescript
 * partition(['hello-world', 'foo'], '-')
 * // returns [['hello', '-', 'world'], ['foo', '', '']]
 * ```
 */
export function partition(a: NDArray | string[], sep: string): NDArray {
  const arr = Array.isArray(a) ? NDArray.fromStringArray(a) : a;

  if (!arr.isStringArray) {
    throw new Error("Input must be a string array");
  }

  const newShape = [...arr.shape, 3];
  const result = NDArray.emptyString(newShape);

  for (let i = 0; i < arr.size; i++) {
    const s = arr.getStringFlat(i);
    const idx = s.indexOf(sep);
    const base = i * 3;

    if (idx === -1) {
      result.setStringFlat(base, s);
      result.setStringFlat(base + 1, "");
      result.setStringFlat(base + 2, "");
    } else {
      result.setStringFlat(base, s.slice(0, idx));
      result.setStringFlat(base + 1, sep);
      result.setStringFlat(base + 2, s.slice(idx + sep.length));
    }
  }

  return result;
}

/**
 * Partition each element around the last occurrence of sep.
 *
 * @param a - Input string array
 * @param sep - Separator string
 * @returns Array with shape (..., 3) containing (before, sep, after)
 */
export function rpartition(a: NDArray | string[], sep: string): NDArray {
  const arr = Array.isArray(a) ? NDArray.fromStringArray(a) : a;

  if (!arr.isStringArray) {
    throw new Error("Input must be a string array");
  }

  const newShape = [...arr.shape, 3];
  const result = NDArray.emptyString(newShape);

  for (let i = 0; i < arr.size; i++) {
    const s = arr.getStringFlat(i);
    const idx = s.lastIndexOf(sep);
    const base = i * 3;

    if (idx === -1) {
      result.setStringFlat(base, "");
      result.setStringFlat(base + 1, "");
      result.setStringFlat(base + 2, s);
    } else {
      result.setStringFlat(base, s.slice(0, idx));
      result.setStringFlat(base + 1, sep);
      result.setStringFlat(base + 2, s.slice(idx + sep.length));
    }
  }

  return result;
}

/* ============ Encoding ============ */

/**
 * Encode strings to bytes using specified encoding.
 * Returns an array of Uint8Arrays (one per string element).
 *
 * @param a - Input string array
 * @param encoding - Character encoding (default 'utf-8')
 * @returns Array of Uint8Arrays
 */
export function encode(
  a: NDArray | string[],
  _encoding: string = "utf-8",
): Uint8Array[] {
  const arr = Array.isArray(a) ? NDArray.fromStringArray(a) : a;

  if (!arr.isStringArray) {
    throw new Error("Input must be a string array");
  }

  const encoder = new TextEncoder(); // Always UTF-8 in browsers
  const result: Uint8Array[] = [];

  for (let i = 0; i < arr.size; i++) {
    const s = arr.getStringFlat(i);
    result.push(encoder.encode(s));
  }

  return result;
}

/**
 * Decode bytes to strings using specified encoding.
 *
 * @param a - Array of byte arrays (Uint8Arrays)
 * @param encoding - Character encoding (default 'utf-8')
 * @returns String NDArray
 */
export function decode(a: Uint8Array[], encoding: string = "utf-8"): NDArray {
  const decoder = new TextDecoder(encoding);
  const strings: string[] = [];

  for (const bytes of a) {
    strings.push(decoder.decode(bytes));
  }

  return NDArray.fromStringArray(strings);
}

/* ============ Formatting ============ */

/**
 * Return element-wise (a % i), i.e. string formatting.
 *
 * This is the printf-style string formatting operation.
 * Each string in `a` is interpreted as a format string,
 * and the corresponding values in `values` are substituted.
 *
 * @param a - Input format string array
 * @param values - Values to substitute (array or single value)
 * @returns String array with formatted strings
 *
 * @example
 * ```typescript
 * mod(['%d items', '%.2f dollars'], [42, 19.99])
 * // returns ['42 items', '19.99 dollars']
 *
 * mod(['Hello, %s!'], 'World')
 * // returns ['Hello, World!']
 * ```
 */
export function mod(
  a: NDArray | string[],
  values: unknown | unknown[],
): NDArray {
  const arr = Array.isArray(a) ? NDArray.fromStringArray(a) : a;

  if (!arr.isStringArray) {
    throw new Error("First argument must be a string array");
  }

  const valArray = Array.isArray(values) ? values : [values];
  const result = NDArray.emptyString(arr.shape);

  for (let i = 0; i < arr.size; i++) {
    const formatStr = arr.getStringFlat(i);
    const val = valArray[i % valArray.length];
    result.setStringFlat(i, _sprintf(formatStr, val));
  }

  return result;
}

/**
 * Simple sprintf-like formatting for a single value.
 */
function _sprintf(format: string, value: unknown): string {
  // Handle the most common format specifiers
  return format.replace(/%([+\-#0 ]*)(\d*)(\.(\d+))?([diouxXeEfFgGcrsb%])/g,
    (match, flags, width, _precisionGroup, precision, specifier) => {
      if (specifier === '%') return '%';

      const leftAlign = flags.includes('-');
      const padZero = flags.includes('0') && !leftAlign;
      const showSign = flags.includes('+');
      const spaceSign = flags.includes(' ');
      const minWidth = parseInt(width) || 0;
      const prec = precision !== undefined ? parseInt(precision) : undefined;

      let result: string;

      switch (specifier) {
        case 'd':
        case 'i': {
          const num = Math.trunc(Number(value));
          const sign = num < 0 ? '-' : (showSign ? '+' : (spaceSign ? ' ' : ''));
          result = sign + Math.abs(num).toString();
          break;
        }
        case 'o': {
          const num = Math.trunc(Number(value));
          result = Math.abs(num).toString(8);
          if (flags.includes('#') && num !== 0) result = '0' + result;
          break;
        }
        case 'x': {
          const num = Math.trunc(Number(value));
          result = Math.abs(num).toString(16);
          if (flags.includes('#') && num !== 0) result = '0x' + result;
          break;
        }
        case 'X': {
          const num = Math.trunc(Number(value));
          result = Math.abs(num).toString(16).toUpperCase();
          if (flags.includes('#') && num !== 0) result = '0X' + result;
          break;
        }
        case 'e': {
          const num = Number(value);
          result = num.toExponential(prec ?? 6);
          break;
        }
        case 'E': {
          const num = Number(value);
          result = num.toExponential(prec ?? 6).toUpperCase();
          break;
        }
        case 'f':
        case 'F': {
          const num = Number(value);
          const sign = num < 0 ? '-' : (showSign ? '+' : (spaceSign ? ' ' : ''));
          result = sign + Math.abs(num).toFixed(prec ?? 6);
          break;
        }
        case 'g': {
          const num = Number(value);
          const p = prec ?? 6;
          result = num.toPrecision(p > 0 ? p : 1);
          // Remove trailing zeros
          if (!flags.includes('#')) {
            result = result.replace(/\.?0+$/, '');
          }
          break;
        }
        case 'G': {
          const num = Number(value);
          const p = prec ?? 6;
          result = num.toPrecision(p > 0 ? p : 1).toUpperCase();
          if (!flags.includes('#')) {
            result = result.replace(/\.?0+$/, '');
          }
          break;
        }
        case 'c':
          result = String.fromCharCode(Number(value));
          break;
        case 's':
          result = String(value);
          if (prec !== undefined) result = result.slice(0, prec);
          break;
        case 'r':
          result = JSON.stringify(value);
          break;
        case 'b': {
          const num = Math.trunc(Number(value));
          result = Math.abs(num).toString(2);
          if (flags.includes('#') && num !== 0) result = '0b' + result;
          break;
        }
        default:
          result = match;
      }

      // Apply width padding
      if (minWidth > result.length) {
        const padChar = padZero ? '0' : ' ';
        const padding = padChar.repeat(minWidth - result.length);
        result = leftAlign ? result + padding : padding + result;
      }

      return result;
    }
  );
}

/**
 * Return element-wise translation of string using a translation table.
 *
 * Each character in the string is mapped through the translation table.
 * Characters not in the table are left unchanged.
 *
 * @param a - Input string array
 * @param table - Translation table as a Map, object, or string pairs
 * @param deletechars - Characters to delete (optional)
 * @returns String array with translated strings
 *
 * @example
 * ```typescript
 * // Using an object as translation table
 * translate(['hello'], {'h': 'H', 'e': '3'})
 * // returns ['H3llo']
 *
 * // Deleting characters
 * translate(['hello world'], {}, 'eo')
 * // returns ['hll wrld']
 * ```
 */
export function translate(
  a: NDArray | string[],
  table: Map<string, string> | Record<string, string> | null,
  deletechars: string = "",
): NDArray {
  // Convert table to Map if needed
  let tableMap: Map<string, string>;
  if (table === null) {
    tableMap = new Map();
  } else if (table instanceof Map) {
    tableMap = table;
  } else {
    tableMap = new Map(Object.entries(table));
  }

  // Create set of characters to delete
  const deleteSet = new Set(deletechars);

  return _applyTransform(a, (s) => {
    let result = "";
    for (const char of s) {
      if (deleteSet.has(char)) {
        continue; // Skip deleted characters
      }
      const mapped = tableMap.get(char);
      result += mapped !== undefined ? mapped : char;
    }
    return result;
  });
}

/**
 * Return element-wise sliced strings.
 *
 * This is the NumPy 2.0 string slicing function, providing
 * Python-style string slicing for each element.
 *
 * @param a - Input string array
 * @param start - Start index (default 0)
 * @param stop - Stop index (default: end of string)
 * @param step - Step size (default 1)
 * @returns String array with sliced strings
 *
 * @example
 * ```typescript
 * slice(['hello', 'world'], 1, 4)
 * // returns ['ell', 'orl']
 *
 * slice(['hello', 'world'], null, null, -1)
 * // returns ['olleh', 'dlrow'] (reversed)
 *
 * slice(['hello'], 0, null, 2)
 * // returns ['hlo'] (every 2nd char)
 * ```
 */
export function slice(
  a: NDArray | string[],
  start: number | null = null,
  stop: number | null = null,
  step: number | null = null,
): NDArray {
  const actualStep = step ?? 1;

  if (actualStep === 0) {
    throw new Error("slice step cannot be zero");
  }

  return _applyTransform(a, (s) => {
    const len = s.length;

    // Handle negative indices and defaults
    let startIdx: number;
    let stopIdx: number;

    if (actualStep > 0) {
      startIdx = start === null ? 0 : _normalizeIndex(start, len);
      stopIdx = stop === null ? len : _normalizeIndex(stop, len);
    } else {
      startIdx = start === null ? len - 1 : _normalizeIndex(start, len, true);
      stopIdx = stop === null ? -1 : _normalizeIndex(stop, len, true);
    }

    // Build result string
    let result = "";

    if (actualStep > 0) {
      for (let i = startIdx; i < stopIdx && i < len; i += actualStep) {
        if (i >= 0) {
          result += s[i];
        }
      }
    } else {
      for (let i = startIdx; i > stopIdx && i >= 0; i += actualStep) {
        if (i < len) {
          result += s[i];
        }
      }
    }

    return result;
  });
}

/**
 * Normalize a slice index to be within bounds.
 */
function _normalizeIndex(idx: number, len: number, forNegativeStep: boolean = false): number {
  if (idx < 0) {
    idx = len + idx;
  }

  if (forNegativeStep) {
    // For negative step, we want different clamping behavior
    return Math.max(-1, Math.min(idx, len - 1));
  }

  return Math.max(0, Math.min(idx, len));
}

/* ============ Helper Functions ============ */

/**
 * Apply a string transform to each element.
 */
function _applyTransform(
  a: NDArray | string[],
  transform: (s: string) => string,
): NDArray {
  const arr = Array.isArray(a) ? NDArray.fromStringArray(a) : a;

  if (!arr.isStringArray) {
    throw new Error("Input must be a string array");
  }

  const result = NDArray.emptyString(arr.shape);

  for (let i = 0; i < arr.size; i++) {
    result.setStringFlat(i, transform(arr.getStringFlat(i)));
  }

  return result;
}

/**
 * Escape special regex characters in a string.
 */
function _escapeRegex(s: string): string {
  return s.replace(/[.*+?^${}()|[\]\\]/g, "\\$&");
}

/**
 * Validate that two arrays have matching shapes.
 */
function _validateShapes(a1: NDArray, a2: NDArray): void {
  if (a1.shape.length !== a2.shape.length) {
    throw new Error(
      `Shape mismatch: ${JSON.stringify(a1.shape)} vs ${JSON.stringify(a2.shape)}`,
    );
  }
  for (let i = 0; i < a1.shape.length; i++) {
    if (a1.shape[i] !== a2.shape[i]) {
      throw new Error(
        `Shape mismatch: ${JSON.stringify(a1.shape)} vs ${JSON.stringify(a2.shape)}`,
      );
    }
  }
}
