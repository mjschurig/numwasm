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
