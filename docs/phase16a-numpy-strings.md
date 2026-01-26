# Phase 16a: numpy.strings Implementation Plan

Complete implementation roadmap for the NumJS-WASM string operations module, providing NumPy-compatible vectorized string functions.

---

## ⚠️ Implementation Guideline

**The original NumPy source code is available in `/numpy`.**

Key reference files:
- `numpy/_core/strings.py` - Main string operations (1,813 lines)
- `numpy/strings/__init__.py` - Public exports
- `numpy/_core/src/umath/string_*.cpp` - C implementations

Implementation should follow NumPy's algorithms and error handling for consistency.

---

## Current State (Pre-Phase 16a)

```
src/ts/
├── NDArray.ts         # Core array class
├── types.ts           # Type definitions
├── dtype.ts           # Type utilities
├── broadcast.ts       # Broadcasting functions
└── index.ts           # Public exports
```

**Required Infrastructure:**
- NDArray with string dtype support (needs to be added)
- Broadcasting for shape compatibility
- Element-wise operations pattern

---

## Phase 16a Dependency Tree

```
PHASE 16a: NUMPY.STRINGS
│
├── 16a.1 String Array Infrastructure
│   ├── String dtype in type system
│   ├── NDArray string storage (TypedArray of string refs)
│   ├── getStringFlat, setStringFlat methods
│   └── fromStringArray factory method
│
├── 16a.2 Comparison Functions
│   ├── equal(x1, x2) → boolean array
│   ├── not_equal(x1, x2) → boolean array
│   ├── less(x1, x2) → boolean array
│   ├── less_equal(x1, x2) → boolean array
│   ├── greater(x1, x2) → boolean array
│   └── greater_equal(x1, x2) → boolean array
│
├── 16a.3 String Property Testing
│   ├── isalpha(x) → boolean array
│   ├── isdigit(x) → boolean array
│   ├── isalnum(x) → boolean array
│   ├── isspace(x) → boolean array
│   ├── islower(x) → boolean array
│   ├── isupper(x) → boolean array
│   ├── istitle(x) → boolean array
│   ├── isdecimal(x) → boolean array
│   ├── isnumeric(x) → boolean array
│   └── str_len(x) → integer array
│
├── 16a.4 Search and Indexing
│   ├── find(a, sub, start, end) → index or -1
│   ├── rfind(a, sub, start, end) → last index or -1
│   ├── index(a, sub, start, end) → index (raises on not found)
│   ├── rindex(a, sub, start, end) → index (raises on not found)
│   ├── count(a, sub, start, end) → occurrence count
│   ├── startswith(a, prefix, start, end) → boolean array
│   └── endswith(a, suffix, start, end) → boolean array
│
└── 16a.5 String Manipulation
    ├── Case Conversion
    │   ├── upper(a) → uppercase
    │   ├── lower(a) → lowercase
    │   ├── swapcase(a) → swap case
    │   ├── capitalize(a) → capitalize first
    │   └── title(a) → titlecase
    │
    ├── Concatenation
    │   ├── add(x1, x2) → concatenate strings
    │   └── multiply(a, i) → repeat strings
    │
    ├── Whitespace
    │   ├── strip(a, chars) → remove leading/trailing
    │   ├── lstrip(a, chars) → remove leading
    │   ├── rstrip(a, chars) → remove trailing
    │   └── expandtabs(a, tabsize) → expand tabs
    │
    ├── Replacement
    │   └── replace(a, old, new, count) → replace substrings
    │
    ├── Alignment
    │   ├── center(a, width, fillchar) → center in width
    │   ├── ljust(a, width, fillchar) → left justify
    │   ├── rjust(a, width, fillchar) → right justify
    │   └── zfill(a, width) → pad with zeros
    │
    ├── Partitioning
    │   ├── partition(a, sep) → (before, sep, after)
    │   └── rpartition(a, sep) → (before, sep, after)
    │
    └── Encoding
        ├── encode(a, encoding) → bytes
        └── decode(a, encoding) → strings

Dependencies: NDArray core, DType system, Broadcasting
```

---

## Detailed Implementation Specifications

### 16a.1 String Array Infrastructure

**File:** `src/ts/NDArray.ts` (additions)

```typescript
// Add to NDArray class

/**
 * Internal string storage for string arrays.
 * Maps flat index to string value.
 */
private _stringData: Map<number, string> | null = null;

/**
 * Check if this array has string dtype.
 */
get isStringArray(): boolean {
  return this.dtype === DType.String;
}

/**
 * Get string value at flat index.
 */
getStringFlat(index: number): string {
  if (!this._stringData) {
    throw new Error('Not a string array');
  }
  return this._stringData.get(index) ?? '';
}

/**
 * Set string value at flat index.
 */
setStringFlat(index: number, value: string): void {
  if (!this._stringData) {
    this._stringData = new Map();
  }
  this._stringData.set(index, value);
}

/**
 * Create string array from JavaScript string array.
 */
static fromStringArray(
  data: string[] | string[][],
  shape?: number[]
): NDArray {
  const flat = Array.isArray(data[0])
    ? (data as string[][]).flat()
    : data as string[];

  const inferredShape = shape ?? [flat.length];
  const arr = new NDArray();

  arr._shape = inferredShape;
  arr._dtype = DType.String;
  arr._size = flat.reduce((a, b) => a * b, 1);
  arr._stringData = new Map();

  for (let i = 0; i < flat.length; i++) {
    arr._stringData.set(i, flat[i]);
  }

  return arr;
}

/**
 * Create empty string array.
 */
static emptyString(shape: number[]): NDArray {
  const arr = new NDArray();
  arr._shape = shape;
  arr._dtype = DType.String;
  arr._size = shape.reduce((a, b) => a * b, 1);
  arr._stringData = new Map();
  return arr;
}
```

**File:** `src/ts/types.ts` (additions)

```typescript
// Add to DType enum
export enum DType {
  // ... existing types ...
  String = 'string',
}
```

---

### 16a.2 Comparison Functions

**File:** `src/ts/strings/compare.ts`

```typescript
import { NDArray } from '../NDArray.js';
import { DType } from '../types.js';
import { broadcastArrays } from '../broadcast.js';

/**
 * Return (x1 == x2) element-wise for string arrays.
 *
 * @param x1 - First string array
 * @param x2 - Second string array
 * @returns Boolean array of comparison results
 *
 * @example
 * equal(['hello', 'world'], ['hello', 'test'])
 * // returns [true, false]
 */
export function equal(x1: NDArray | string[], x2: NDArray | string[]): NDArray {
  const [a1, a2] = _prepareStringArrays(x1, x2);
  const result = NDArray.empty(a1.shape, DType.Bool);

  for (let i = 0; i < a1.size; i++) {
    result.setFlat(i, a1.getStringFlat(i) === a2.getStringFlat(i) ? 1 : 0);
  }

  return result;
}

/**
 * Return (x1 != x2) element-wise for string arrays.
 */
export function not_equal(x1: NDArray | string[], x2: NDArray | string[]): NDArray {
  const [a1, a2] = _prepareStringArrays(x1, x2);
  const result = NDArray.empty(a1.shape, DType.Bool);

  for (let i = 0; i < a1.size; i++) {
    result.setFlat(i, a1.getStringFlat(i) !== a2.getStringFlat(i) ? 1 : 0);
  }

  return result;
}

/**
 * Return (x1 < x2) element-wise for string arrays.
 * Comparison is lexicographic (dictionary order).
 */
export function less(x1: NDArray | string[], x2: NDArray | string[]): NDArray {
  const [a1, a2] = _prepareStringArrays(x1, x2);
  const result = NDArray.empty(a1.shape, DType.Bool);

  for (let i = 0; i < a1.size; i++) {
    result.setFlat(i, a1.getStringFlat(i) < a2.getStringFlat(i) ? 1 : 0);
  }

  return result;
}

/**
 * Return (x1 <= x2) element-wise for string arrays.
 */
export function less_equal(x1: NDArray | string[], x2: NDArray | string[]): NDArray {
  const [a1, a2] = _prepareStringArrays(x1, x2);
  const result = NDArray.empty(a1.shape, DType.Bool);

  for (let i = 0; i < a1.size; i++) {
    result.setFlat(i, a1.getStringFlat(i) <= a2.getStringFlat(i) ? 1 : 0);
  }

  return result;
}

/**
 * Return (x1 > x2) element-wise for string arrays.
 */
export function greater(x1: NDArray | string[], x2: NDArray | string[]): NDArray {
  const [a1, a2] = _prepareStringArrays(x1, x2);
  const result = NDArray.empty(a1.shape, DType.Bool);

  for (let i = 0; i < a1.size; i++) {
    result.setFlat(i, a1.getStringFlat(i) > a2.getStringFlat(i) ? 1 : 0);
  }

  return result;
}

/**
 * Return (x1 >= x2) element-wise for string arrays.
 */
export function greater_equal(x1: NDArray | string[], x2: NDArray | string[]): NDArray {
  const [a1, a2] = _prepareStringArrays(x1, x2);
  const result = NDArray.empty(a1.shape, DType.Bool);

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
  const result = NDArray.empty(a1.shape, DType.Int32);

  for (let i = 0; i < a1.size; i++) {
    const s1 = a1.getStringFlat(i);
    const s2 = a2.getStringFlat(i);
    result.setFlat(i, s1 < s2 ? -1 : (s1 > s2 ? 1 : 0));
  }

  return result;
}

/* ============ Helper Functions ============ */

function _prepareStringArrays(
  x1: NDArray | string[],
  x2: NDArray | string[]
): [NDArray, NDArray] {
  const a1 = Array.isArray(x1) ? NDArray.fromStringArray(x1) : x1;
  const a2 = Array.isArray(x2) ? NDArray.fromStringArray(x2) : x2;

  // Broadcast to common shape
  return broadcastArrays(a1, a2) as [NDArray, NDArray];
}
```

---

### 16a.3 String Property Testing

**File:** `src/ts/strings/properties.ts`

```typescript
import { NDArray } from '../NDArray.js';
import { DType } from '../types.js';

/**
 * Return true for each element if all characters are alphabetic.
 *
 * @example
 * isalpha(['hello', '123', 'abc123'])
 * // returns [true, false, false]
 */
export function isalpha(a: NDArray | string[]): NDArray {
  return _applyStringTest(a, (s) => s.length > 0 && /^[a-zA-Z]+$/.test(s));
}

/**
 * Return true for each element if all characters are digits.
 */
export function isdigit(a: NDArray | string[]): NDArray {
  return _applyStringTest(a, (s) => s.length > 0 && /^[0-9]+$/.test(s));
}

/**
 * Return true for each element if all characters are alphanumeric.
 */
export function isalnum(a: NDArray | string[]): NDArray {
  return _applyStringTest(a, (s) => s.length > 0 && /^[a-zA-Z0-9]+$/.test(s));
}

/**
 * Return true for each element if all characters are whitespace.
 */
export function isspace(a: NDArray | string[]): NDArray {
  return _applyStringTest(a, (s) => s.length > 0 && /^\s+$/.test(s));
}

/**
 * Return true for each element if all cased characters are lowercase
 * and there is at least one cased character.
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
        const isUpper = char === char.toUpperCase() && char !== char.toLowerCase();
        const isLower = char === char.toLowerCase() && char !== char.toUpperCase();
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
 */
export function isdecimal(a: NDArray | string[]): NDArray {
  return _applyStringTest(a, (s) => s.length > 0 && /^[0-9]+$/.test(s));
}

/**
 * Return true for each element if all characters are numeric.
 * Numeric characters include digit characters and characters
 * that have the Unicode numeric value property (subscripts, fractions, etc.).
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
      if (code === 0xB2 || code === 0xB3 || code === 0xB9) continue;
      // Subscript digits
      if (code >= 0x2080 && code <= 0x2089) continue;
      // Vulgar fractions
      if (code >= 0xBC && code <= 0xBE) continue;
      // Fullwidth digits
      if (code >= 0xFF10 && code <= 0xFF19) continue;
      return false;
    }
    return true;
  });
}

/**
 * Return the length of each element.
 *
 * @example
 * str_len(['hello', 'hi', ''])
 * // returns [5, 2, 0]
 */
export function str_len(a: NDArray | string[]): NDArray {
  const arr = Array.isArray(a) ? NDArray.fromStringArray(a) : a;
  const result = NDArray.empty(arr.shape, DType.Int32);

  for (let i = 0; i < arr.size; i++) {
    result.setFlat(i, arr.getStringFlat(i).length);
  }

  return result;
}

/* ============ Helper Functions ============ */

function _applyStringTest(
  a: NDArray | string[],
  test: (s: string) => boolean
): NDArray {
  const arr = Array.isArray(a) ? NDArray.fromStringArray(a) : a;
  const result = NDArray.empty(arr.shape, DType.Bool);

  for (let i = 0; i < arr.size; i++) {
    result.setFlat(i, test(arr.getStringFlat(i)) ? 1 : 0);
  }

  return result;
}
```

---

### 16a.4 Search and Indexing

**File:** `src/ts/strings/search.ts`

```typescript
import { NDArray } from '../NDArray.js';
import { DType } from '../types.js';

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
 * find(['hello world', 'test'], 'o')
 * // returns [4, -1]
 */
export function find(
  a: NDArray | string[],
  sub: string,
  start: number = 0,
  end: number | null = null
): NDArray {
  const arr = Array.isArray(a) ? NDArray.fromStringArray(a) : a;
  const result = NDArray.empty(arr.shape, DType.Int32);

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
  const result = NDArray.empty(arr.shape, DType.Int32);

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
 * @example
 * count(['aabababa', 'foo'], 'aba')
 * // returns [2, 0]
 */
export function count(
  a: NDArray | string[],
  sub: string,
  start: number = 0,
  end: number | null = null
): NDArray {
  const arr = Array.isArray(a) ? NDArray.fromStringArray(a) : a;
  const result = NDArray.empty(arr.shape, DType.Int32);

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
 * @example
 * startswith(['hello', 'world', 'help'], 'hel')
 * // returns [true, false, true]
 */
export function startswith(
  a: NDArray | string[],
  prefix: string,
  start: number = 0,
  end: number | null = null
): NDArray {
  const arr = Array.isArray(a) ? NDArray.fromStringArray(a) : a;
  const result = NDArray.empty(arr.shape, DType.Bool);

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
 * @example
 * endswith(['hello', 'world', 'jello'], 'llo')
 * // returns [true, false, true]
 */
export function endswith(
  a: NDArray | string[],
  suffix: string,
  start: number = 0,
  end: number | null = null
): NDArray {
  const arr = Array.isArray(a) ? NDArray.fromStringArray(a) : a;
  const result = NDArray.empty(arr.shape, DType.Bool);

  for (let i = 0; i < arr.size; i++) {
    const s = arr.getStringFlat(i);
    const searchEnd = end ?? s.length;
    const searchStr = s.slice(start, searchEnd);
    result.setFlat(i, searchStr.endsWith(suffix) ? 1 : 0);
  }

  return result;
}

/* ============ Error Classes ============ */

export class ValueError extends Error {
  constructor(message: string) {
    super(message);
    this.name = 'ValueError';
  }
}
```

---

### 16a.5 String Manipulation

**File:** `src/ts/strings/manipulation.ts`

```typescript
import { NDArray } from '../NDArray.js';
import { DType } from '../types.js';
import { broadcastArrays } from '../broadcast.js';

/* ============ Case Conversion ============ */

/**
 * Return element-wise copy with uppercase characters converted to lowercase.
 */
export function lower(a: NDArray | string[]): NDArray {
  return _applyTransform(a, (s) => s.toLowerCase());
}

/**
 * Return element-wise copy with lowercase characters converted to uppercase.
 */
export function upper(a: NDArray | string[]): NDArray {
  return _applyTransform(a, (s) => s.toUpperCase());
}

/**
 * Return element-wise copy with uppercase and lowercase swapped.
 */
export function swapcase(a: NDArray | string[]): NDArray {
  return _applyTransform(a, (s) => {
    return s.split('').map(c => {
      if (c === c.toLowerCase() && c !== c.toUpperCase()) {
        return c.toUpperCase();
      }
      if (c === c.toUpperCase() && c !== c.toLowerCase()) {
        return c.toLowerCase();
      }
      return c;
    }).join('');
  });
}

/**
 * Return element-wise copy with first character capitalized and rest lowercase.
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
 */
export function title(a: NDArray | string[]): NDArray {
  return _applyTransform(a, (s) => {
    return s.replace(/\b\w/g, (c) => c.toUpperCase())
            .replace(/\B\w/g, (c) => c.toLowerCase());
  });
}

/* ============ Concatenation ============ */

/**
 * Return element-wise string concatenation.
 *
 * @example
 * add(['hello', 'good'], [' world', 'bye'])
 * // returns ['hello world', 'goodbye']
 */
export function add(x1: NDArray | string[], x2: NDArray | string[]): NDArray {
  const a1 = Array.isArray(x1) ? NDArray.fromStringArray(x1) : x1;
  const a2 = Array.isArray(x2) ? NDArray.fromStringArray(x2) : x2;
  const [b1, b2] = broadcastArrays(a1, a2);

  const result = NDArray.emptyString(b1.shape);

  for (let i = 0; i < b1.size; i++) {
    result.setStringFlat(i, b1.getStringFlat(i) + b2.getStringFlat(i));
  }

  return result;
}

/**
 * Return element-wise (a * i), string repeated i times.
 *
 * @example
 * multiply(['ab', 'cd'], 3)
 * // returns ['ababab', 'cdcdcd']
 */
export function multiply(a: NDArray | string[], i: number | NDArray): NDArray {
  const arr = Array.isArray(a) ? NDArray.fromStringArray(a) : a;
  const result = NDArray.emptyString(arr.shape);

  if (typeof i === 'number') {
    for (let j = 0; j < arr.size; j++) {
      const s = arr.getStringFlat(j);
      result.setStringFlat(j, i > 0 ? s.repeat(i) : '');
    }
  } else {
    for (let j = 0; j < arr.size; j++) {
      const s = arr.getStringFlat(j);
      const count = i.getFlat(j % i.size);
      result.setStringFlat(j, count > 0 ? s.repeat(count) : '');
    }
  }

  return result;
}

/* ============ Whitespace ============ */

/**
 * Return element-wise copy with leading and trailing characters removed.
 *
 * @param a - Input string array
 * @param chars - Characters to remove (default: whitespace)
 */
export function strip(a: NDArray | string[], chars: string | null = null): NDArray {
  return _applyTransform(a, (s) => {
    if (chars === null) return s.trim();
    const regex = new RegExp(
      `^[${_escapeRegex(chars)}]+|[${_escapeRegex(chars)}]+$`,
      'g'
    );
    return s.replace(regex, '');
  });
}

/**
 * Return element-wise copy with leading characters removed.
 */
export function lstrip(a: NDArray | string[], chars: string | null = null): NDArray {
  return _applyTransform(a, (s) => {
    if (chars === null) return s.trimStart();
    const regex = new RegExp(`^[${_escapeRegex(chars)}]+`);
    return s.replace(regex, '');
  });
}

/**
 * Return element-wise copy with trailing characters removed.
 */
export function rstrip(a: NDArray | string[], chars: string | null = null): NDArray {
  return _applyTransform(a, (s) => {
    if (chars === null) return s.trimEnd();
    const regex = new RegExp(`[${_escapeRegex(chars)}]+$`);
    return s.replace(regex, '');
  });
}

/**
 * Return element-wise copy with tabs expanded to spaces.
 *
 * @param a - Input string array
 * @param tabsize - Tab stop spacing (default 8)
 */
export function expandtabs(a: NDArray | string[], tabsize: number = 8): NDArray {
  return _applyTransform(a, (s) => {
    let result = '';
    let col = 0;
    for (const char of s) {
      if (char === '\t') {
        const spaces = tabsize - (col % tabsize);
        result += ' '.repeat(spaces);
        col += spaces;
      } else if (char === '\n' || char === '\r') {
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
 */
export function replace(
  a: NDArray | string[],
  old: string,
  new_: string,
  count: number = -1
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
 */
export function center(
  a: NDArray | string[],
  width: number,
  fillchar: string = ' '
): NDArray {
  if (fillchar.length !== 1) {
    throw new TypeError('fillchar must be exactly one character');
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
 */
export function ljust(
  a: NDArray | string[],
  width: number,
  fillchar: string = ' '
): NDArray {
  if (fillchar.length !== 1) {
    throw new TypeError('fillchar must be exactly one character');
  }

  return _applyTransform(a, (s) => {
    if (s.length >= width) return s;
    return s + fillchar.repeat(width - s.length);
  });
}

/**
 * Return element-wise right-justified string of given width.
 */
export function rjust(
  a: NDArray | string[],
  width: number,
  fillchar: string = ' '
): NDArray {
  if (fillchar.length !== 1) {
    throw new TypeError('fillchar must be exactly one character');
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
 * @example
 * zfill(['-42', '42'], 5)
 * // returns ['-0042', '00042']
 */
export function zfill(a: NDArray | string[], width: number): NDArray {
  return _applyTransform(a, (s) => {
    if (s.length >= width) return s;

    const sign = (s[0] === '+' || s[0] === '-') ? s[0] : '';
    const rest = sign ? s.slice(1) : s;
    const padLength = width - s.length;

    return sign + '0'.repeat(padLength) + rest;
  });
}

/* ============ Partitioning ============ */

/**
 * Partition each element around the first occurrence of sep.
 *
 * @returns Array with shape (..., 3) containing (before, sep, after)
 *
 * @example
 * partition(['hello-world', 'foo'], '-')
 * // returns [['hello', '-', 'world'], ['foo', '', '']]
 */
export function partition(a: NDArray | string[], sep: string): NDArray {
  const arr = Array.isArray(a) ? NDArray.fromStringArray(a) : a;
  const newShape = [...arr.shape, 3];
  const result = NDArray.emptyString(newShape);

  for (let i = 0; i < arr.size; i++) {
    const s = arr.getStringFlat(i);
    const idx = s.indexOf(sep);
    const base = i * 3;

    if (idx === -1) {
      result.setStringFlat(base, s);
      result.setStringFlat(base + 1, '');
      result.setStringFlat(base + 2, '');
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
 * @returns Array with shape (..., 3) containing (before, sep, after)
 */
export function rpartition(a: NDArray | string[], sep: string): NDArray {
  const arr = Array.isArray(a) ? NDArray.fromStringArray(a) : a;
  const newShape = [...arr.shape, 3];
  const result = NDArray.emptyString(newShape);

  for (let i = 0; i < arr.size; i++) {
    const s = arr.getStringFlat(i);
    const idx = s.lastIndexOf(sep);
    const base = i * 3;

    if (idx === -1) {
      result.setStringFlat(base, '');
      result.setStringFlat(base + 1, '');
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
 *
 * @param a - Input string array
 * @param encoding - Character encoding (default 'utf-8')
 */
export function encode(
  a: NDArray | string[],
  encoding: string = 'utf-8'
): NDArray {
  const arr = Array.isArray(a) ? NDArray.fromStringArray(a) : a;
  const result = NDArray.emptyBytes(arr.shape);

  const encoder = new TextEncoder(); // Always UTF-8 in browsers

  for (let i = 0; i < arr.size; i++) {
    const s = arr.getStringFlat(i);
    const bytes = encoder.encode(s);
    result.setBytesFlat(i, bytes);
  }

  return result;
}

/**
 * Decode bytes to strings using specified encoding.
 *
 * @param a - Input byte array
 * @param encoding - Character encoding (default 'utf-8')
 */
export function decode(a: NDArray, encoding: string = 'utf-8'): NDArray {
  const result = NDArray.emptyString(a.shape);
  const decoder = new TextDecoder(encoding);

  for (let i = 0; i < a.size; i++) {
    const bytes = a.getBytesFlat(i);
    result.setStringFlat(i, decoder.decode(bytes));
  }

  return result;
}

/* ============ Helper Functions ============ */

function _applyTransform(
  a: NDArray | string[],
  transform: (s: string) => string
): NDArray {
  const arr = Array.isArray(a) ? NDArray.fromStringArray(a) : a;
  const result = NDArray.emptyString(arr.shape);

  for (let i = 0; i < arr.size; i++) {
    result.setStringFlat(i, transform(arr.getStringFlat(i)));
  }

  return result;
}

function _escapeRegex(s: string): string {
  return s.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
}
```

---

### Module Index

**File:** `src/ts/strings/index.ts`

```typescript
/**
 * NumJS String Operations Module
 *
 * Provides vectorized string operations on arrays of strings,
 * compatible with NumPy's numpy.strings module.
 */

// Comparison functions
export {
  equal,
  not_equal,
  less,
  less_equal,
  greater,
  greater_equal,
  compare_chararrays,
} from './compare.js';

// Property testing
export {
  isalpha,
  isdigit,
  isalnum,
  isspace,
  islower,
  isupper,
  istitle,
  isdecimal,
  isnumeric,
  str_len,
} from './properties.js';

// Search and indexing
export {
  find,
  rfind,
  index,
  rindex,
  count,
  startswith,
  endswith,
  ValueError,
} from './search.js';

// String manipulation
export {
  // Case conversion
  lower,
  upper,
  swapcase,
  capitalize,
  title,
  // Concatenation
  add,
  multiply,
  // Whitespace
  strip,
  lstrip,
  rstrip,
  expandtabs,
  // Replacement
  replace,
  // Alignment
  center,
  ljust,
  rjust,
  zfill,
  // Partitioning
  partition,
  rpartition,
  // Encoding
  encode,
  decode,
} from './manipulation.js';
```

---

## File Changes Summary

### New Files to Create

```
src/ts/strings/
├── index.ts           # Public exports (~50 lines)
├── compare.ts         # Comparison functions (~120 lines)
├── properties.ts      # Property testing (~150 lines)
├── search.ts          # Search functions (~180 lines)
└── manipulation.ts    # Manipulation functions (~350 lines)
```

### Files to Modify

```
src/ts/types.ts
├── Add DType.String enum value
└── Add string-related type definitions

src/ts/NDArray.ts
├── Add _stringData private field
├── Add isStringArray property
├── Add getStringFlat, setStringFlat methods
├── Add fromStringArray static method
├── Add emptyString static method
└── Add emptyBytes, getBytesFlat, setBytesFlat for encoding

src/ts/index.ts
└── Export strings module
```

---

## Implementation Order

```
Week 1: Foundation & Core Functions
├── Day 1: String array infrastructure in NDArray
│   ├── Add DType.String to types.ts
│   ├── Add _stringData storage
│   ├── Add getStringFlat, setStringFlat
│   └── Add fromStringArray factory
│
├── Day 2: Comparison functions
│   ├── equal, not_equal
│   ├── less, less_equal
│   ├── greater, greater_equal
│   └── compare_chararrays
│
├── Day 3: Property testing functions
│   ├── isalpha, isdigit, isalnum
│   ├── isspace, islower, isupper
│   ├── istitle, isdecimal, isnumeric
│   └── str_len
│
├── Day 4: Search functions
│   ├── find, rfind
│   ├── index, rindex
│   ├── count
│   └── startswith, endswith
│
└── Day 5: Basic manipulation
    ├── upper, lower, swapcase
    ├── capitalize, title
    └── strip, lstrip, rstrip

Week 2: Advanced Functions & Testing
├── Day 1: Alignment functions
│   ├── center, ljust, rjust
│   └── zfill
│
├── Day 2: Concatenation & replacement
│   ├── add, multiply
│   └── replace
│
├── Day 3: Partitioning & tabs
│   ├── partition, rpartition
│   └── expandtabs
│
├── Day 4: Encoding
│   ├── encode
│   └── decode
│
└── Day 5: Tests & documentation
    ├── Unit tests for all functions
    ├── Integration tests
    └── API documentation
```

---

## Verification Plan

```bash
# Build
npm run build

# Run tests
npm test -- --grep "strings"

# Test cases to verify:

# Comparison
✓ equal(['a', 'b'], ['a', 'c']) → [true, false]
✓ less(['a', 'b'], ['b', 'a']) → [true, false]
✓ Broadcasting works with scalar string

# Properties
✓ isalpha(['abc', '123', 'a1']) → [true, false, false]
✓ isupper(['ABC', 'Abc', '123']) → [true, false, false]
✓ str_len(['', 'a', 'hello']) → [0, 1, 5]

# Search
✓ find(['hello', 'world'], 'o') → [4, 1]
✓ count(['ababa', 'foo'], 'aba') → [1, 0] (non-overlapping)
✓ startswith(['hello', 'world'], 'hel') → [true, false]

# Manipulation
✓ upper(['hello', 'World']) → ['HELLO', 'WORLD']
✓ strip(['  hi  ', '\tfoo\n']) → ['hi', 'foo']
✓ center(['ab', 'x'], 5) → [' ab  ', '  x  ']
✓ zfill(['-42', '7'], 5) → ['-0042', '00007']
✓ partition(['a-b-c', 'foo'], '-') → [['a', '-', 'b-c'], ['foo', '', '']]
```

---

## API Compatibility Notes

### NumPy API Mapping

```python
# NumPy
np.strings.equal(a, b)
np.strings.upper(a)
np.strings.find(a, 'sub')

# NumJS
strings.equal(a, b)
strings.upper(a)
strings.find(a, 'sub')
```

### Differences from NumPy

1. **String storage**: NumJS uses Map<number, string> instead of fixed-width character arrays
2. **Unicode handling**: NumJS relies on JavaScript's native Unicode support
3. **Encoding**: Limited to encodings supported by TextEncoder/TextDecoder

---

## Performance Considerations

### TypeScript Implementation
- Leverage JavaScript's optimized string methods
- Avoid regex compilation in loops (pre-compile patterns)
- Use Map for sparse string arrays

### Future WASM Optimization
If performance becomes critical, consider WASM for:
- Large-scale string comparisons
- Repeated search operations
- Batch transformations

However, for most use cases, TypeScript implementation should be sufficient since JavaScript engines are highly optimized for string operations.

---

## Estimated Lines of Code

| File | Lines |
|------|-------|
| compare.ts | ~120 |
| properties.ts | ~150 |
| search.ts | ~180 |
| manipulation.ts | ~350 |
| index.ts | ~50 |
| **Total** | **~850** |

Plus modifications to NDArray.ts (~100 lines) and types.ts (~10 lines).
