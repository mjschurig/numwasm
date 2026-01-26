# Phase 16: Additional NumPy Modules Implementation Plan

Complete implementation roadmap for the NumJS-WASM additional modules, providing NumPy-compatible string operations, polynomial mathematics, masked arrays, record arrays, and testing utilities.

---

## ⚠️ Implementation Guideline

**The original NumPy source code is available in `/numpy`.**

Key reference files:
- `numpy/_core/strings.py` - String operations (1,813 lines)
- `numpy/polynomial/` - Polynomial module (12,417 lines total)
- `numpy/ma/` - Masked arrays (11,502 lines total)
- `numpy/_core/records.py` - Record arrays (1,086 lines)
- `numpy/testing/` - Testing utilities (13,227 lines total)

Implementation should follow NumPy's algorithms and error handling for consistency.

---

## Current State (Pre-Phase 16)

```
src/wasm/
├── ndarray.h/c        # Core NDArray with views, slicing
├── dtype.h/c          # DType system
├── broadcast.h/c      # Broadcasting
├── indexing.h/c       # Index operations
├── pairwise_sum.h/c   # Accurate summation
└── logic.c            # Logical operations

src/ts/
├── NDArray.ts         # Core array class
├── types.ts           # Type definitions
├── dtype.ts           # Type utilities
├── broadcast.ts       # Broadcasting functions
├── indexing.ts        # Index operations
├── slice.ts           # Slicing utilities
├── iterators.ts       # Iterators
└── index.ts           # Public exports
```

**Existing Infrastructure Used:**
- NDArray with contiguity flags (C_CONTIGUOUS, F_CONTIGUOUS)
- DType system with all numeric types
- Broadcasting for shape compatibility
- View system for efficient memory sharing

---

## Phase 16 Dependency Tree

```
PHASE 16: ADDITIONAL MODULES
│
├── 16.1 numpy.strings (TypeScript + Optional WASM)
│   ├── 16.1.1 Comparison Functions
│   │   ├── equal(x1, x2) → boolean array
│   │   ├── not_equal(x1, x2) → boolean array
│   │   ├── less(x1, x2), less_equal(x1, x2)
│   │   └── greater(x1, x2), greater_equal(x1, x2)
│   │
│   ├── 16.1.2 String Property Testing
│   │   ├── isalpha(x), isdigit(x), isalnum(x)
│   │   ├── isspace(x), islower(x), isupper(x)
│   │   ├── istitle(x), isdecimal(x), isnumeric(x)
│   │   └── str_len(x) → integer array
│   │
│   ├── 16.1.3 Search and Indexing
│   │   ├── find(a, sub, start, end) → index or -1
│   │   ├── rfind(a, sub, start, end) → last index or -1
│   │   ├── index(a, sub, start, end) → index (raises on not found)
│   │   ├── rindex(a, sub, start, end)
│   │   ├── count(a, sub, start, end) → occurrence count
│   │   └── startswith(a, prefix), endswith(a, suffix)
│   │
│   └── 16.1.4 String Manipulation
│       ├── upper(a), lower(a), swapcase(a)
│       ├── capitalize(a), title(a)
│       ├── add(x1, x2), multiply(a, i)
│       ├── strip(a, chars), lstrip(a), rstrip(a)
│       ├── replace(a, old, new, count)
│       ├── center(a, width), ljust(a, width), rjust(a, width)
│       ├── zfill(a, width), expandtabs(a, tabsize)
│       ├── partition(a, sep), rpartition(a, sep)
│       └── encode(a, encoding), decode(a, encoding)
│
│   Dependencies: NDArray core, DType system, Broadcasting
│
├── 16.2 numpy.rec (TypeScript)
│   ├── 16.2.1 Format Parser
│   │   └── format_parser(formats, names, titles, aligned, byteorder)
│   │
│   ├── 16.2.2 Record Class
│   │   ├── record (single record access)
│   │   └── Field access via attributes
│   │
│   ├── 16.2.3 RecArray Class
│   │   ├── recarray (full record array)
│   │   ├── Attribute-based field access
│   │   └── field(name) method
│   │
│   └── 16.2.4 Creation Functions
│       ├── array(data, dtype, formats, names)
│       ├── fromarrays(arrayList, names)
│       ├── fromrecords(recList, names)
│       ├── fromstring(datastring, dtype)
│       └── fromfile(fd, dtype)
│
│   Dependencies: NDArray core, Structured DTypes
│
├── 16.3 numpy.polynomial (TypeScript + WASM)
│   ├── 16.3.1 ABCPolyBase (Abstract Base Class)
│   │   ├── Constructor, copy, basis, identity
│   │   ├── Arithmetic: __add__, __sub__, __mul__, __truediv__, __pow__
│   │   ├── Comparison: __eq__, __ne__
│   │   ├── Conversion: cast, convert, mapparms
│   │   ├── Calculus: deriv, integ
│   │   ├── Analysis: roots, trim, cutdeg, truncate
│   │   └── Fitting: fit (least squares)
│   │
│   ├── 16.3.2 Polynomial (Power Series)
│   │   ├── polyadd, polysub, polymul, polydiv, polypow
│   │   ├── polyval, polyval2d, polyval3d
│   │   ├── polyder, polyint
│   │   ├── polyfit, polyvander
│   │   ├── polyroots, polycompanion
│   │   └── polyfromroots
│   │
│   ├── 16.3.3 Chebyshev Polynomials
│   │   ├── chebadd, chebsub, chebmul, chebdiv, chebpow
│   │   ├── chebval, chebval2d, chebval3d
│   │   ├── chebder, chebint
│   │   ├── chebfit, chebvander
│   │   ├── chebroots, chebcompanion
│   │   ├── chebinterpolate
│   │   └── Conversions: poly2cheb, cheb2poly
│   │
│   ├── 16.3.4 Legendre Polynomials
│   │   └── (Same pattern as Chebyshev)
│   │
│   ├── 16.3.5 Hermite Polynomials (Physicist's)
│   │   └── (Same pattern as Chebyshev)
│   │
│   ├── 16.3.6 HermiteE Polynomials (Probabilist's)
│   │   └── (Same pattern as Chebyshev)
│   │
│   ├── 16.3.7 Laguerre Polynomials
│   │   └── (Same pattern as Chebyshev)
│   │
│   └── 16.3.8 Utility Functions (polyutils)
│       ├── as_series, trimseq, trimcoef
│       ├── getdomain, mapdomain, mapparms
│       └── PolyError, PolyDomainWarning
│
│   Dependencies: NDArray, Broadcasting, Linear Algebra (Phase 13)
│
├── 16.4 numpy.ma (Masked Arrays) (TypeScript)
│   ├── 16.4.1 MaskedArray Class
│   │   ├── _data, _mask, _fill_value, _hardmask attributes
│   │   ├── Constructor with mask handling
│   │   ├── View and copy semantics
│   │   └── Mask propagation rules
│   │
│   ├── 16.4.2 Array Creation
│   │   ├── masked_array(data, mask)
│   │   ├── masked_equal, masked_greater, masked_less
│   │   ├── masked_inside, masked_outside
│   │   ├── masked_where, masked_invalid
│   │   ├── masked_values, masked_object
│   │   └── zeros, ones, empty (masked versions)
│   │
│   ├── 16.4.3 Mask Operations
│   │   ├── make_mask, make_mask_none, make_mask_descr
│   │   ├── getmask, getmaskarray, getdata
│   │   ├── is_mask, is_masked
│   │   ├── mask_or, nomask
│   │   └── harden_mask, soften_mask
│   │
│   ├── 16.4.4 Arithmetic Operations
│   │   ├── add, subtract, multiply, divide (with mask propagation)
│   │   ├── power, mod, remainder
│   │   ├── Comparison operators (masked)
│   │   └── Logical and bitwise operators
│   │
│   ├── 16.4.5 Mathematical Functions
│   │   ├── Trigonometric: sin, cos, tan, arcsin, arccos, arctan
│   │   ├── Hyperbolic: sinh, cosh, tanh, arcsinh, arccosh, arctanh
│   │   ├── Exponential/Log: exp, log, log10, log2
│   │   └── Other: sqrt, ceil, floor, round, fabs
│   │
│   ├── 16.4.6 Reductions
│   │   ├── sum, prod, mean, std, var (respecting masks)
│   │   ├── min, max, ptp, argmin, argmax
│   │   ├── count (non-masked elements)
│   │   ├── cumsum, cumprod
│   │   └── all, any
│   │
│   ├── 16.4.7 Shape Manipulation
│   │   ├── reshape, ravel, flatten, squeeze
│   │   ├── transpose, swapaxes
│   │   ├── compressed() → 1D non-masked
│   │   └── filled(fill_value) → regular array
│   │
│   ├── 16.4.8 Extras Module
│   │   ├── average, median, cov, corrcoef
│   │   ├── polyfit, vander
│   │   ├── unique, intersect1d, union1d
│   │   └── notmasked_edges, notmasked_contiguous
│   │
│   └── 16.4.9 Masked Records (mrecords)
│       ├── mrecarray class
│       └── MaskedRecord class
│
│   Dependencies: NDArray, Broadcasting, All math operations
│
└── 16.5 numpy.testing (TypeScript)
    ├── 16.5.1 Equality Assertions
    │   ├── assert_equal(actual, desired)
    │   ├── assert_array_equal(x, y)
    │   ├── assert_string_equal(actual, desired)
    │   └── assert_allclose(actual, desired, rtol, atol)
    │
    ├── 16.5.2 Approximate Assertions
    │   ├── assert_almost_equal(actual, desired, decimal)
    │   ├── assert_approx_equal(actual, desired, significant)
    │   ├── assert_array_almost_equal(x, y, decimal)
    │   └── assert_array_almost_equal_nulp(x, y, nulp)
    │
    ├── 16.5.3 Ordering Assertions
    │   ├── assert_array_less(x, y)
    │   ├── assert_array_max_ulp(a, b, maxulp)
    │   └── assert_(val, msg)
    │
    ├── 16.5.4 Exception/Warning Assertions
    │   ├── assert_raises(exception_class, callable)
    │   ├── assert_raises_regex(exception_class, regex, callable)
    │   ├── assert_warns(warning_class, func)
    │   └── assert_no_warnings(func)
    │
    ├── 16.5.5 Context Managers
    │   ├── suppress_warnings(message, category)
    │   └── clear_and_catch_warnings()
    │
    └── 16.5.6 Utility Functions
        ├── measure(code_snippet, times)
        ├── tempdir(), temppath()
        └── Environment flags (IS_WASM, etc.)

    Dependencies: Core types only
```

---

## Detailed Implementation Specifications

### 16.1 numpy.strings Module

#### 16.1.1 String Comparison Functions

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
 */
export function equal(x1: NDArray | string[], x2: NDArray | string[]): NDArray {
  const [a1, a2] = broadcastArrays(asStringArray(x1), asStringArray(x2));
  const result = NDArray.empty(a1.shape, DType.Bool);

  for (let i = 0; i < a1.size; i++) {
    result.setFlat(i, getString(a1, i) === getString(a2, i) ? 1 : 0);
  }

  return result;
}

/**
 * Return (x1 != x2) element-wise for string arrays.
 */
export function not_equal(x1: NDArray | string[], x2: NDArray | string[]): NDArray {
  const [a1, a2] = broadcastArrays(asStringArray(x1), asStringArray(x2));
  const result = NDArray.empty(a1.shape, DType.Bool);

  for (let i = 0; i < a1.size; i++) {
    result.setFlat(i, getString(a1, i) !== getString(a2, i) ? 1 : 0);
  }

  return result;
}

/**
 * Return (x1 < x2) element-wise for string arrays.
 * Comparison is lexicographic.
 */
export function less(x1: NDArray | string[], x2: NDArray | string[]): NDArray {
  const [a1, a2] = broadcastArrays(asStringArray(x1), asStringArray(x2));
  const result = NDArray.empty(a1.shape, DType.Bool);

  for (let i = 0; i < a1.size; i++) {
    result.setFlat(i, getString(a1, i) < getString(a2, i) ? 1 : 0);
  }

  return result;
}

/**
 * Return (x1 <= x2) element-wise for string arrays.
 */
export function less_equal(x1: NDArray | string[], x2: NDArray | string[]): NDArray {
  const [a1, a2] = broadcastArrays(asStringArray(x1), asStringArray(x2));
  const result = NDArray.empty(a1.shape, DType.Bool);

  for (let i = 0; i < a1.size; i++) {
    result.setFlat(i, getString(a1, i) <= getString(a2, i) ? 1 : 0);
  }

  return result;
}

/**
 * Return (x1 > x2) element-wise for string arrays.
 */
export function greater(x1: NDArray | string[], x2: NDArray | string[]): NDArray {
  const [a1, a2] = broadcastArrays(asStringArray(x1), asStringArray(x2));
  const result = NDArray.empty(a1.shape, DType.Bool);

  for (let i = 0; i < a1.size; i++) {
    result.setFlat(i, getString(a1, i) > getString(a2, i) ? 1 : 0);
  }

  return result;
}

/**
 * Return (x1 >= x2) element-wise for string arrays.
 */
export function greater_equal(x1: NDArray | string[], x2: NDArray | string[]): NDArray {
  const [a1, a2] = broadcastArrays(asStringArray(x1), asStringArray(x2));
  const result = NDArray.empty(a1.shape, DType.Bool);

  for (let i = 0; i < a1.size; i++) {
    result.setFlat(i, getString(a1, i) >= getString(a2, i) ? 1 : 0);
  }

  return result;
}

/* ============ Helper Functions ============ */

function asStringArray(x: NDArray | string[]): NDArray {
  if (Array.isArray(x)) {
    // Convert JS array to NDArray with string dtype
    return NDArray.fromStringArray(x);
  }
  return x;
}

function getString(arr: NDArray, index: number): string {
  // Get string value at flat index
  return arr.getStringFlat(index);
}
```

---

#### 16.1.2 String Property Testing

**File:** `src/ts/strings/properties.ts`

```typescript
import { NDArray } from '../NDArray.js';
import { DType } from '../types.js';

/**
 * Return true for each element if all characters are alphabetic.
 */
export function isalpha(a: NDArray): NDArray {
  return _applyStringTest(a, (s) => /^[a-zA-Z]+$/.test(s) && s.length > 0);
}

/**
 * Return true for each element if all characters are digits.
 */
export function isdigit(a: NDArray): NDArray {
  return _applyStringTest(a, (s) => /^[0-9]+$/.test(s) && s.length > 0);
}

/**
 * Return true for each element if all characters are alphanumeric.
 */
export function isalnum(a: NDArray): NDArray {
  return _applyStringTest(a, (s) => /^[a-zA-Z0-9]+$/.test(s) && s.length > 0);
}

/**
 * Return true for each element if all characters are whitespace.
 */
export function isspace(a: NDArray): NDArray {
  return _applyStringTest(a, (s) => /^\s+$/.test(s) && s.length > 0);
}

/**
 * Return true for each element if all cased characters are lowercase.
 */
export function islower(a: NDArray): NDArray {
  return _applyStringTest(a, (s) => {
    if (!/[a-zA-Z]/.test(s)) return false;
    return s === s.toLowerCase();
  });
}

/**
 * Return true for each element if all cased characters are uppercase.
 */
export function isupper(a: NDArray): NDArray {
  return _applyStringTest(a, (s) => {
    if (!/[a-zA-Z]/.test(s)) return false;
    return s === s.toUpperCase();
  });
}

/**
 * Return true for each element if the string is titlecased.
 */
export function istitle(a: NDArray): NDArray {
  return _applyStringTest(a, (s) => {
    if (s.length === 0) return false;
    // Check if first letter of each word is uppercase, rest lowercase
    const words = s.split(/\s+/);
    return words.every(word => {
      if (word.length === 0) return true;
      const first = word[0];
      const rest = word.slice(1);
      return first === first.toUpperCase() && rest === rest.toLowerCase();
    });
  });
}

/**
 * Return true for each element if all characters are decimal.
 */
export function isdecimal(a: NDArray): NDArray {
  return _applyStringTest(a, (s) => /^[0-9]+$/.test(s) && s.length > 0);
}

/**
 * Return true for each element if all characters are numeric.
 * Includes digits, superscripts, fractions, etc.
 */
export function isnumeric(a: NDArray): NDArray {
  return _applyStringTest(a, (s) => {
    if (s.length === 0) return false;
    // Unicode numeric category - simplified check
    return /^[\d\u00B2\u00B3\u00B9\u00BC-\u00BE\u2070-\u2079\u2080-\u2089]+$/.test(s);
  });
}

/**
 * Return the length of each element.
 */
export function str_len(a: NDArray): NDArray {
  const result = NDArray.empty(a.shape, DType.Int32);

  for (let i = 0; i < a.size; i++) {
    const s = a.getStringFlat(i);
    result.setFlat(i, s.length);
  }

  return result;
}

/* ============ Helper Functions ============ */

function _applyStringTest(a: NDArray, test: (s: string) => boolean): NDArray {
  const result = NDArray.empty(a.shape, DType.Bool);

  for (let i = 0; i < a.size; i++) {
    const s = a.getStringFlat(i);
    result.setFlat(i, test(s) ? 1 : 0);
  }

  return result;
}
```

---

#### 16.1.3 String Search and Indexing

**File:** `src/ts/strings/search.ts`

```typescript
import { NDArray } from '../NDArray.js';
import { DType } from '../types.js';
import { broadcastArrays } from '../broadcast.js';

/**
 * For each element, return the lowest index where substring is found.
 * Returns -1 if substring is not found.
 *
 * @param a - Input string array
 * @param sub - Substring to search for
 * @param start - Start index for search (default 0)
 * @param end - End index for search (default string length)
 */
export function find(
  a: NDArray,
  sub: string | NDArray,
  start: number = 0,
  end: number | null = null
): NDArray {
  const result = NDArray.empty(a.shape, DType.Int32);
  const subStr = typeof sub === 'string' ? sub : null;

  for (let i = 0; i < a.size; i++) {
    const s = a.getStringFlat(i);
    const searchSub = subStr ?? (sub as NDArray).getStringFlat(i % (sub as NDArray).size);
    const searchEnd = end ?? s.length;
    const searchStr = s.slice(start, searchEnd);
    const idx = searchStr.indexOf(searchSub);
    result.setFlat(i, idx === -1 ? -1 : idx + start);
  }

  return result;
}

/**
 * For each element, return the highest index where substring is found.
 * Returns -1 if substring is not found.
 */
export function rfind(
  a: NDArray,
  sub: string | NDArray,
  start: number = 0,
  end: number | null = null
): NDArray {
  const result = NDArray.empty(a.shape, DType.Int32);
  const subStr = typeof sub === 'string' ? sub : null;

  for (let i = 0; i < a.size; i++) {
    const s = a.getStringFlat(i);
    const searchSub = subStr ?? (sub as NDArray).getStringFlat(i % (sub as NDArray).size);
    const searchEnd = end ?? s.length;
    const searchStr = s.slice(start, searchEnd);
    const idx = searchStr.lastIndexOf(searchSub);
    result.setFlat(i, idx === -1 ? -1 : idx + start);
  }

  return result;
}

/**
 * Like find, but raises ValueError if substring is not found.
 */
export function index(
  a: NDArray,
  sub: string | NDArray,
  start: number = 0,
  end: number | null = null
): NDArray {
  const result = find(a, sub, start, end);

  // Check for -1 values
  for (let i = 0; i < result.size; i++) {
    if (result.getFlat(i) === -1) {
      throw new ValueError('substring not found');
    }
  }

  return result;
}

/**
 * Like rfind, but raises ValueError if substring is not found.
 */
export function rindex(
  a: NDArray,
  sub: string | NDArray,
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
 */
export function count(
  a: NDArray,
  sub: string | NDArray,
  start: number = 0,
  end: number | null = null
): NDArray {
  const result = NDArray.empty(a.shape, DType.Int32);
  const subStr = typeof sub === 'string' ? sub : null;

  for (let i = 0; i < a.size; i++) {
    const s = a.getStringFlat(i);
    const searchSub = subStr ?? (sub as NDArray).getStringFlat(i % (sub as NDArray).size);
    const searchEnd = end ?? s.length;
    const searchStr = s.slice(start, searchEnd);

    if (searchSub.length === 0) {
      result.setFlat(i, searchStr.length + 1);
    } else {
      let cnt = 0;
      let pos = 0;
      while ((pos = searchStr.indexOf(searchSub, pos)) !== -1) {
        cnt++;
        pos += searchSub.length;
      }
      result.setFlat(i, cnt);
    }
  }

  return result;
}

/**
 * Return true for each element if the string starts with prefix.
 */
export function startswith(
  a: NDArray,
  prefix: string | NDArray,
  start: number = 0,
  end: number | null = null
): NDArray {
  const result = NDArray.empty(a.shape, DType.Bool);
  const prefixStr = typeof prefix === 'string' ? prefix : null;

  for (let i = 0; i < a.size; i++) {
    const s = a.getStringFlat(i);
    const searchPrefix = prefixStr ?? (prefix as NDArray).getStringFlat(i % (prefix as NDArray).size);
    const searchEnd = end ?? s.length;
    const searchStr = s.slice(start, searchEnd);
    result.setFlat(i, searchStr.startsWith(searchPrefix) ? 1 : 0);
  }

  return result;
}

/**
 * Return true for each element if the string ends with suffix.
 */
export function endswith(
  a: NDArray,
  suffix: string | NDArray,
  start: number = 0,
  end: number | null = null
): NDArray {
  const result = NDArray.empty(a.shape, DType.Bool);
  const suffixStr = typeof suffix === 'string' ? suffix : null;

  for (let i = 0; i < a.size; i++) {
    const s = a.getStringFlat(i);
    const searchSuffix = suffixStr ?? (suffix as NDArray).getStringFlat(i % (suffix as NDArray).size);
    const searchEnd = end ?? s.length;
    const searchStr = s.slice(start, searchEnd);
    result.setFlat(i, searchStr.endsWith(searchSuffix) ? 1 : 0);
  }

  return result;
}

class ValueError extends Error {
  constructor(message: string) {
    super(message);
    this.name = 'ValueError';
  }
}
```

---

#### 16.1.4 String Manipulation

**File:** `src/ts/strings/manipulation.ts`

```typescript
import { NDArray } from '../NDArray.js';
import { DType } from '../types.js';
import { broadcastArrays } from '../broadcast.js';

/**
 * Return element-wise string concatenation.
 */
export function add(x1: NDArray, x2: NDArray): NDArray {
  const [a1, a2] = broadcastArrays(x1, x2);
  const result = NDArray.emptyString(a1.shape);

  for (let i = 0; i < a1.size; i++) {
    const s1 = a1.getStringFlat(i);
    const s2 = a2.getStringFlat(i);
    result.setStringFlat(i, s1 + s2);
  }

  return result;
}

/**
 * Return string multiplied by integer element-wise.
 */
export function multiply(a: NDArray, i: number | NDArray): NDArray {
  const result = NDArray.emptyString(a.shape);
  const isScalar = typeof i === 'number';

  for (let j = 0; j < a.size; j++) {
    const s = a.getStringFlat(j);
    const count = isScalar ? i : (i as NDArray).getFlat(j % (i as NDArray).size);
    result.setStringFlat(j, count > 0 ? s.repeat(count as number) : '');
  }

  return result;
}

/**
 * Return a copy with uppercase characters converted to lowercase.
 */
export function lower(a: NDArray): NDArray {
  return _applyStringTransform(a, (s) => s.toLowerCase());
}

/**
 * Return a copy with lowercase characters converted to uppercase.
 */
export function upper(a: NDArray): NDArray {
  return _applyStringTransform(a, (s) => s.toUpperCase());
}

/**
 * Return a copy with uppercase and lowercase swapped.
 */
export function swapcase(a: NDArray): NDArray {
  return _applyStringTransform(a, (s) => {
    return s.split('').map(c => {
      if (c === c.toLowerCase()) return c.toUpperCase();
      return c.toLowerCase();
    }).join('');
  });
}

/**
 * Return a copy with first character capitalized.
 */
export function capitalize(a: NDArray): NDArray {
  return _applyStringTransform(a, (s) => {
    if (s.length === 0) return s;
    return s[0].toUpperCase() + s.slice(1).toLowerCase();
  });
}

/**
 * Return a titlecased version (words start with uppercase).
 */
export function title(a: NDArray): NDArray {
  return _applyStringTransform(a, (s) => {
    return s.replace(/\b\w/g, (c) => c.toUpperCase());
  });
}

/**
 * Return a copy with leading and trailing characters removed.
 */
export function strip(a: NDArray, chars: string | null = null): NDArray {
  return _applyStringTransform(a, (s) => {
    if (chars === null) return s.trim();
    const regex = new RegExp(`^[${escapeRegex(chars)}]+|[${escapeRegex(chars)}]+$`, 'g');
    return s.replace(regex, '');
  });
}

/**
 * Return a copy with leading characters removed.
 */
export function lstrip(a: NDArray, chars: string | null = null): NDArray {
  return _applyStringTransform(a, (s) => {
    if (chars === null) return s.trimStart();
    const regex = new RegExp(`^[${escapeRegex(chars)}]+`);
    return s.replace(regex, '');
  });
}

/**
 * Return a copy with trailing characters removed.
 */
export function rstrip(a: NDArray, chars: string | null = null): NDArray {
  return _applyStringTransform(a, (s) => {
    if (chars === null) return s.trimEnd();
    const regex = new RegExp(`[${escapeRegex(chars)}]+$`);
    return s.replace(regex, '');
  });
}

/**
 * Return a copy with all occurrences of substring replaced.
 */
export function replace(
  a: NDArray,
  old: string,
  new_: string,
  count: number = -1
): NDArray {
  return _applyStringTransform(a, (s) => {
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

/**
 * Return centered string of given width.
 */
export function center(a: NDArray, width: number, fillchar: string = ' '): NDArray {
  return _applyStringTransform(a, (s) => {
    if (s.length >= width) return s;
    const totalPad = width - s.length;
    const leftPad = Math.floor(totalPad / 2);
    const rightPad = totalPad - leftPad;
    return fillchar.repeat(leftPad) + s + fillchar.repeat(rightPad);
  });
}

/**
 * Return left-justified string of given width.
 */
export function ljust(a: NDArray, width: number, fillchar: string = ' '): NDArray {
  return _applyStringTransform(a, (s) => {
    if (s.length >= width) return s;
    return s + fillchar.repeat(width - s.length);
  });
}

/**
 * Return right-justified string of given width.
 */
export function rjust(a: NDArray, width: number, fillchar: string = ' '): NDArray {
  return _applyStringTransform(a, (s) => {
    if (s.length >= width) return s;
    return fillchar.repeat(width - s.length) + s;
  });
}

/**
 * Return numeric string left-padded with zeros.
 */
export function zfill(a: NDArray, width: number): NDArray {
  return _applyStringTransform(a, (s) => {
    if (s.length >= width) return s;
    const sign = (s[0] === '+' || s[0] === '-') ? s[0] : '';
    const rest = sign ? s.slice(1) : s;
    const padLength = width - s.length;
    return sign + '0'.repeat(padLength) + rest;
  });
}

/**
 * Return a copy with tabs expanded to spaces.
 */
export function expandtabs(a: NDArray, tabsize: number = 8): NDArray {
  return _applyStringTransform(a, (s) => {
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

/**
 * Partition each element around the first occurrence of sep.
 * Returns array of 3-tuples: (before, sep, after)
 */
export function partition(a: NDArray, sep: string): NDArray {
  // Returns shape (..., 3) with string dtype
  const newShape = [...a.shape, 3];
  const result = NDArray.emptyString(newShape);

  for (let i = 0; i < a.size; i++) {
    const s = a.getStringFlat(i);
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
 */
export function rpartition(a: NDArray, sep: string): NDArray {
  const newShape = [...a.shape, 3];
  const result = NDArray.emptyString(newShape);

  for (let i = 0; i < a.size; i++) {
    const s = a.getStringFlat(i);
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

/**
 * Encode to bytes using specified encoding.
 */
export function encode(a: NDArray, encoding: string = 'utf-8'): NDArray {
  // Returns byte array
  return _applyStringTransform(a, (s) => {
    const encoder = new TextEncoder();
    return encoder.encode(s);
  }, true);
}

/**
 * Decode from bytes using specified encoding.
 */
export function decode(a: NDArray, encoding: string = 'utf-8'): NDArray {
  return _applyStringTransform(a, (s) => {
    const decoder = new TextDecoder(encoding);
    // Assume s is byte data
    return decoder.decode(new Uint8Array(s as unknown as ArrayBuffer));
  });
}

/* ============ Helper Functions ============ */

function _applyStringTransform(
  a: NDArray,
  transform: (s: string) => string | Uint8Array,
  returnBytes: boolean = false
): NDArray {
  const result = returnBytes ? NDArray.emptyBytes(a.shape) : NDArray.emptyString(a.shape);

  for (let i = 0; i < a.size; i++) {
    const s = a.getStringFlat(i);
    const transformed = transform(s);
    if (returnBytes) {
      result.setBytesFlat(i, transformed as Uint8Array);
    } else {
      result.setStringFlat(i, transformed as string);
    }
  }

  return result;
}

function escapeRegex(s: string): string {
  return s.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
}
```

---

### 16.2 numpy.rec Module

#### 16.2.1 Format Parser

**File:** `src/ts/rec/format_parser.ts`

```typescript
import { DType, StructuredDType, FieldDescriptor } from '../types.js';

/**
 * Class to convert format, names, titles to a structured dtype.
 *
 * Parses format strings like 'f8, i4, S10' into proper field descriptors.
 */
export class format_parser {
  private _dtype: StructuredDType;

  /**
   * Create a format parser.
   *
   * @param formats - Format specification (string or list)
   * @param names - Field names (string or list)
   * @param titles - Optional field titles
   * @param aligned - If true, pad fields for C alignment
   * @param byteorder - Byte order ('>', '<', '=')
   */
  constructor(
    formats: string | string[],
    names: string | string[] | null = null,
    titles: string | string[] | null = null,
    aligned: boolean = false,
    byteorder: string | null = null
  ) {
    const parsedFormats = this._parseFormats(formats, aligned);
    const parsedNames = this._parseNames(names, parsedFormats.length);
    const parsedTitles = this._parseTitles(titles, parsedFormats.length);

    this._dtype = this._createDtype(parsedFormats, parsedNames, parsedTitles, byteorder);
  }

  /**
   * The resulting structured dtype.
   */
  get dtype(): StructuredDType {
    return this._dtype;
  }

  /**
   * Parse format specification into list of dtype strings.
   */
  private _parseFormats(formats: string | string[], aligned: boolean): string[] {
    let formatList: string[];

    if (typeof formats === 'string') {
      // Split by comma, handling whitespace
      formatList = formats.split(',').map(f => f.trim());
    } else {
      formatList = [...formats];
    }

    // Validate each format
    for (const fmt of formatList) {
      if (!this._isValidFormat(fmt)) {
        throw new Error(`Invalid format string: ${fmt}`);
      }
    }

    return formatList;
  }

  /**
   * Parse names into list of field names.
   */
  private _parseNames(names: string | string[] | null, count: number): string[] {
    if (names === null) {
      // Generate default names: f0, f1, f2, ...
      return Array.from({ length: count }, (_, i) => `f${i}`);
    }

    let nameList: string[];
    if (typeof names === 'string') {
      nameList = names.split(',').map(n => n.trim());
    } else {
      nameList = [...names];
    }

    if (nameList.length !== count) {
      throw new Error(
        `Number of names (${nameList.length}) does not match number of formats (${count})`
      );
    }

    // Check for duplicates
    const duplicates = this._findDuplicates(nameList);
    if (duplicates.length > 0) {
      throw new Error(`Duplicate field names: ${duplicates.join(', ')}`);
    }

    return nameList;
  }

  /**
   * Parse titles into list of field titles.
   */
  private _parseTitles(titles: string | string[] | null, count: number): (string | null)[] {
    if (titles === null) {
      return Array(count).fill(null);
    }

    let titleList: (string | null)[];
    if (typeof titles === 'string') {
      titleList = titles.split(',').map(t => t.trim() || null);
    } else {
      titleList = [...titles];
    }

    if (titleList.length !== count) {
      throw new Error(
        `Number of titles (${titleList.length}) does not match number of formats (${count})`
      );
    }

    return titleList;
  }

  /**
   * Create the structured dtype from parsed components.
   */
  private _createDtype(
    formats: string[],
    names: string[],
    titles: (string | null)[],
    byteorder: string | null
  ): StructuredDType {
    const fields: FieldDescriptor[] = [];
    let offset = 0;

    for (let i = 0; i < formats.length; i++) {
      const dtype = this._formatToDtype(formats[i], byteorder);
      const size = this._dtypeSize(dtype);

      fields.push({
        name: names[i],
        dtype: dtype,
        offset: offset,
        title: titles[i],
      });

      offset += size;
    }

    return {
      names: names,
      fields: fields,
      itemsize: offset,
    };
  }

  /**
   * Validate format string.
   */
  private _isValidFormat(fmt: string): boolean {
    // Valid formats: i4, f8, S10, U20, b1, etc.
    return /^[<>=|]?[biufcmMOSUV]\d*$/.test(fmt);
  }

  /**
   * Convert format string to DType.
   */
  private _formatToDtype(fmt: string, byteorder: string | null): DType {
    // Parse format: optional byteorder + type char + optional size
    const match = fmt.match(/^([<>=|])?([biufcmMOSUV])(\d*)$/);
    if (!match) {
      throw new Error(`Invalid format: ${fmt}`);
    }

    const [, order, typeChar, sizeStr] = match;
    const effectiveOrder = order || byteorder || '=';
    const size = sizeStr ? parseInt(sizeStr, 10) : this._defaultSize(typeChar);

    return this._charToDtype(typeChar, size);
  }

  /**
   * Get default size for type character.
   */
  private _defaultSize(typeChar: string): number {
    const defaults: Record<string, number> = {
      'b': 1, 'i': 4, 'u': 4, 'f': 8, 'c': 16,
      'S': 1, 'U': 4, 'V': 1, 'O': 8, 'm': 8, 'M': 8,
    };
    return defaults[typeChar] || 1;
  }

  /**
   * Convert type character to DType enum.
   */
  private _charToDtype(typeChar: string, size: number): DType {
    switch (typeChar) {
      case 'b': return DType.Int8;
      case 'i':
        if (size === 1) return DType.Int8;
        if (size === 2) return DType.Int16;
        if (size === 4) return DType.Int32;
        return DType.Int64;
      case 'u':
        if (size === 1) return DType.UInt8;
        if (size === 2) return DType.UInt16;
        if (size === 4) return DType.UInt32;
        return DType.UInt64;
      case 'f':
        if (size === 2) return DType.Float16;
        if (size === 4) return DType.Float32;
        return DType.Float64;
      case 'c':
        if (size === 8) return DType.Complex64;
        return DType.Complex128;
      case 'S':
      case 'U':
        return DType.String;
      default:
        throw new Error(`Unsupported type character: ${typeChar}`);
    }
  }

  /**
   * Get size of dtype in bytes.
   */
  private _dtypeSize(dtype: DType): number {
    const sizes: Record<DType, number> = {
      [DType.Bool]: 1,
      [DType.Int8]: 1,
      [DType.Int16]: 2,
      [DType.Int32]: 4,
      [DType.Int64]: 8,
      [DType.UInt8]: 1,
      [DType.UInt16]: 2,
      [DType.UInt32]: 4,
      [DType.UInt64]: 8,
      [DType.Float16]: 2,
      [DType.Float32]: 4,
      [DType.Float64]: 8,
      [DType.Complex64]: 8,
      [DType.Complex128]: 16,
      [DType.String]: 0, // Variable
    };
    return sizes[dtype] || 0;
  }

  /**
   * Find duplicate elements in array.
   */
  private _findDuplicates(arr: string[]): string[] {
    const seen = new Set<string>();
    const duplicates: string[] = [];

    for (const item of arr) {
      if (seen.has(item) && !duplicates.includes(item)) {
        duplicates.push(item);
      }
      seen.add(item);
    }

    return duplicates;
  }
}
```

---

#### 16.2.2-3 Record and RecArray Classes

**File:** `src/ts/rec/recarray.ts`

```typescript
import { NDArray } from '../NDArray.js';
import { StructuredDType, FieldDescriptor } from '../types.js';
import { format_parser } from './format_parser.js';

/**
 * A single record from a record array.
 * Allows field access via attributes.
 */
export class record {
  private _data: NDArray;
  private _index: number;
  private _dtype: StructuredDType;

  constructor(data: NDArray, index: number) {
    this._data = data;
    this._index = index;
    this._dtype = data.dtype as unknown as StructuredDType;

    // Create property getters for each field
    for (const field of this._dtype.fields) {
      Object.defineProperty(this, field.name, {
        get: () => this._getField(field.name),
        set: (value) => this._setField(field.name, value),
        enumerable: true,
      });
    }
  }

  /**
   * Get field value by name.
   */
  private _getField(name: string): any {
    return this._data.getField(name).getFlat(this._index);
  }

  /**
   * Set field value by name.
   */
  private _setField(name: string, value: any): void {
    this._data.getField(name).setFlat(this._index, value);
  }

  /**
   * Number of fields.
   */
  get length(): number {
    return this._dtype.fields.length;
  }

  /**
   * String representation.
   */
  toString(): string {
    const values = this._dtype.fields.map(f => `${f.name}=${this._getField(f.name)}`);
    return `record(${values.join(', ')})`;
  }

  /**
   * Pretty print the record.
   */
  pprint(): void {
    console.log(this.toString());
  }
}

/**
 * Record array with named fields.
 * Allows field access as attributes.
 */
export class recarray extends NDArray {
  private _fields: Map<string, NDArray>;

  /**
   * Create a new record array.
   */
  static create(
    shape: number | number[],
    dtype: StructuredDType | null = null,
    formats: string | string[] | null = null,
    names: string | string[] | null = null,
    titles: string | string[] | null = null,
    aligned: boolean = false,
    byteorder: string | null = null
  ): recarray {
    // Determine dtype
    let structDtype: StructuredDType;
    if (dtype !== null) {
      structDtype = dtype;
    } else if (formats !== null) {
      const parser = new format_parser(formats, names, titles, aligned, byteorder);
      structDtype = parser.dtype;
    } else {
      throw new Error('Must specify either dtype or formats');
    }

    // Create base array
    const shapeArr = typeof shape === 'number' ? [shape] : shape;
    const arr = NDArray.emptyStructured(shapeArr, structDtype);

    // Convert to recarray
    return new recarray(arr, structDtype);
  }

  constructor(baseArray: NDArray, dtype: StructuredDType) {
    super();
    // Copy properties from base array
    Object.assign(this, baseArray);

    this._fields = new Map();

    // Create field accessors
    for (const field of dtype.fields) {
      Object.defineProperty(this, field.name, {
        get: () => this.field(field.name),
        set: (value) => this._setFieldArray(field.name, value),
        enumerable: true,
      });
    }
  }

  /**
   * Get a single field as an array.
   */
  field(name: string): NDArray {
    if (!this._fields.has(name)) {
      this._fields.set(name, this.getField(name));
    }
    return this._fields.get(name)!;
  }

  /**
   * Set a field array.
   */
  private _setFieldArray(name: string, value: NDArray | number[]): void {
    const fieldArr = this.field(name);
    if (Array.isArray(value)) {
      for (let i = 0; i < value.length; i++) {
        fieldArr.setFlat(i, value[i]);
      }
    } else {
      // Copy from NDArray
      for (let i = 0; i < value.size; i++) {
        fieldArr.setFlat(i, value.getFlat(i));
      }
    }
  }

  /**
   * Get a single record.
   */
  getRecord(index: number): record {
    return new record(this, index);
  }

  /**
   * Convert to list of tuples.
   */
  tolist(): any[][] {
    const result: any[][] = [];
    const dtype = this.dtype as unknown as StructuredDType;

    for (let i = 0; i < this.shape[0]; i++) {
      const row: any[] = [];
      for (const field of dtype.fields) {
        row.push(this.field(field.name).getFlat(i));
      }
      result.push(row);
    }

    return result;
  }

  /**
   * List field names for tab completion.
   */
  __dir__(): string[] {
    const dtype = this.dtype as unknown as StructuredDType;
    return dtype.fields.map(f => f.name);
  }
}
```

---

#### 16.2.4 Creation Functions

**File:** `src/ts/rec/index.ts`

```typescript
import { NDArray } from '../NDArray.js';
import { StructuredDType } from '../types.js';
import { format_parser } from './format_parser.js';
import { recarray, record } from './recarray.js';

export { format_parser, recarray, record };

/**
 * Create a record array from data.
 *
 * @param data - List of tuples or existing array
 * @param dtype - Explicit dtype
 * @param formats - Format specification
 * @param names - Field names
 * @param titles - Field titles
 */
export function array(
  data: any[][] | NDArray,
  dtype: StructuredDType | null = null,
  shape: number[] | null = null,
  formats: string | string[] | null = null,
  names: string | string[] | null = null,
  titles: string | string[] | null = null,
  aligned: boolean = false,
  byteorder: string | null = null
): recarray {
  // Determine dtype
  let structDtype: StructuredDType;
  if (dtype !== null) {
    structDtype = dtype;
  } else if (formats !== null) {
    const parser = new format_parser(formats, names, titles, aligned, byteorder);
    structDtype = parser.dtype;
  } else {
    throw new Error('Must specify either dtype or formats');
  }

  // Determine shape
  let dataShape: number[];
  if (shape !== null) {
    dataShape = shape;
  } else if (Array.isArray(data)) {
    dataShape = [data.length];
  } else {
    dataShape = data.shape;
  }

  // Create record array
  const result = recarray.create(dataShape, structDtype);

  // Fill with data
  if (Array.isArray(data)) {
    for (let i = 0; i < data.length; i++) {
      const row = data[i];
      for (let j = 0; j < structDtype.fields.length; j++) {
        const field = structDtype.fields[j];
        result.field(field.name).setFlat(i, row[j]);
      }
    }
  }

  return result;
}

/**
 * Create a record array from a list of arrays (columns).
 */
export function fromarrays(
  arrayList: NDArray[],
  dtype: StructuredDType | null = null,
  shape: number[] | null = null,
  formats: string | string[] | null = null,
  names: string | string[] | null = null,
  titles: string | string[] | null = null,
  aligned: boolean = false,
  byteorder: string | null = null
): recarray {
  // Infer formats from arrays if not provided
  if (formats === null && dtype === null) {
    formats = arrayList.map(arr => dtypeToFormat(arr.dtype));
  }

  // Determine shape from first array
  const dataShape = shape ?? arrayList[0].shape;

  // Determine dtype
  let structDtype: StructuredDType;
  if (dtype !== null) {
    structDtype = dtype;
  } else {
    const parser = new format_parser(formats!, names, titles, aligned, byteorder);
    structDtype = parser.dtype;
  }

  // Create and fill
  const result = recarray.create(dataShape, structDtype);

  for (let j = 0; j < arrayList.length; j++) {
    const field = structDtype.fields[j];
    const srcArr = arrayList[j];
    const dstArr = result.field(field.name);

    for (let i = 0; i < srcArr.size; i++) {
      dstArr.setFlat(i, srcArr.getFlat(i));
    }
  }

  return result;
}

/**
 * Create a record array from a list of records (rows).
 */
export function fromrecords(
  recList: any[][],
  dtype: StructuredDType | null = null,
  shape: number[] | null = null,
  formats: string | string[] | null = null,
  names: string | string[] | null = null,
  titles: string | string[] | null = null,
  aligned: boolean = false,
  byteorder: string | null = null
): recarray {
  return array(recList, dtype, shape, formats, names, titles, aligned, byteorder);
}

/**
 * Create a record array from a binary string.
 */
export function fromstring(
  datastring: ArrayBuffer | string,
  dtype: StructuredDType,
  shape: number[] = [],
  offset: number = 0
): recarray {
  // Convert string to ArrayBuffer if needed
  let buffer: ArrayBuffer;
  if (typeof datastring === 'string') {
    const encoder = new TextEncoder();
    buffer = encoder.encode(datastring).buffer;
  } else {
    buffer = datastring;
  }

  // Calculate shape from buffer size if not provided
  const itemsize = dtype.itemsize;
  const numRecords = shape.length > 0
    ? shape.reduce((a, b) => a * b, 1)
    : Math.floor((buffer.byteLength - offset) / itemsize);

  const dataShape = shape.length > 0 ? shape : [numRecords];

  // Create record array
  const result = recarray.create(dataShape, dtype);

  // Parse binary data
  const view = new DataView(buffer, offset);
  let byteOffset = 0;

  for (let i = 0; i < numRecords; i++) {
    for (const field of dtype.fields) {
      const value = readFieldValue(view, byteOffset + field.offset, field.dtype);
      result.field(field.name).setFlat(i, value);
    }
    byteOffset += itemsize;
  }

  return result;
}

/**
 * Create a record array from a file.
 */
export async function fromfile(
  file: File | string,
  dtype: StructuredDType,
  shape: number[] = [],
  offset: number = 0
): Promise<recarray> {
  let buffer: ArrayBuffer;

  if (typeof file === 'string') {
    // Node.js path - use fs
    if (typeof process !== 'undefined') {
      const fs = await import('fs/promises');
      const data = await fs.readFile(file);
      buffer = data.buffer;
    } else {
      throw new Error('File paths not supported in browser');
    }
  } else {
    // Browser File object
    buffer = await file.arrayBuffer();
  }

  return fromstring(buffer, dtype, shape, offset);
}

/* ============ Helper Functions ============ */

import { DType } from '../types.js';

function dtypeToFormat(dtype: DType): string {
  const formatMap: Record<DType, string> = {
    [DType.Bool]: 'b1',
    [DType.Int8]: 'i1',
    [DType.Int16]: 'i2',
    [DType.Int32]: 'i4',
    [DType.Int64]: 'i8',
    [DType.UInt8]: 'u1',
    [DType.UInt16]: 'u2',
    [DType.UInt32]: 'u4',
    [DType.UInt64]: 'u8',
    [DType.Float16]: 'f2',
    [DType.Float32]: 'f4',
    [DType.Float64]: 'f8',
    [DType.Complex64]: 'c8',
    [DType.Complex128]: 'c16',
    [DType.String]: 'U',
  };
  return formatMap[dtype] || 'V';
}

function readFieldValue(view: DataView, offset: number, dtype: DType): any {
  switch (dtype) {
    case DType.Int8: return view.getInt8(offset);
    case DType.Int16: return view.getInt16(offset, true);
    case DType.Int32: return view.getInt32(offset, true);
    case DType.UInt8: return view.getUint8(offset);
    case DType.UInt16: return view.getUint16(offset, true);
    case DType.UInt32: return view.getUint32(offset, true);
    case DType.Float32: return view.getFloat32(offset, true);
    case DType.Float64: return view.getFloat64(offset, true);
    default: return 0;
  }
}
```

---

### 16.3 numpy.polynomial Module

Due to length constraints, I'll provide the key structures. The full implementation follows similar patterns.

#### 16.3.1 ABCPolyBase Abstract Base Class

**File:** `src/ts/polynomial/_polybase.ts`

```typescript
import { NDArray } from '../NDArray.js';

/**
 * Abstract base class for polynomial series.
 *
 * All polynomial classes (Polynomial, Chebyshev, Legendre, etc.)
 * inherit from this class and implement the abstract methods.
 */
export abstract class ABCPolyBase {
  /** Maximum allowed polynomial degree */
  static maxpower: number = 100;

  /** Coefficient array */
  protected _coef: number[];

  /** Domain for the polynomial */
  protected _domain: [number, number];

  /** Window for the polynomial */
  protected _window: [number, number];

  /** Symbol used in string representation */
  protected _symbol: string;

  /* ============ Abstract Properties ============ */

  /** Default domain for this polynomial type */
  abstract get defaultDomain(): [number, number];

  /** Default window for this polynomial type */
  abstract get defaultWindow(): [number, number];

  /** Name of the basis for string representation */
  abstract get basisName(): string;

  /* ============ Abstract Static Methods ============ */

  /** Add two coefficient arrays */
  protected abstract _add(c1: number[], c2: number[]): number[];

  /** Subtract two coefficient arrays */
  protected abstract _sub(c1: number[], c2: number[]): number[];

  /** Multiply two coefficient arrays */
  protected abstract _mul(c1: number[], c2: number[]): number[];

  /** Divide two coefficient arrays */
  protected abstract _div(c1: number[], c2: number[]): [number[], number[]];

  /** Raise coefficient array to power */
  protected abstract _pow(c: number[], n: number): number[];

  /** Evaluate polynomial at x */
  protected abstract _val(x: number | number[], c: number[]): number | number[];

  /** Compute derivative of coefficient array */
  protected abstract _der(c: number[], m: number, scl: number): number[];

  /** Compute integral of coefficient array */
  protected abstract _int(c: number[], m: number, k: number[], lbnd: number, scl: number): number[];

  /** Find roots of polynomial */
  protected abstract _roots(c: number[]): number[];

  /** Create Vandermonde matrix */
  protected abstract _vander(x: number[], deg: number): number[][];

  /* ============ Constructor ============ */

  constructor(
    coef: number[],
    domain: [number, number] | null = null,
    window: [number, number] | null = null,
    symbol: string = 'x'
  ) {
    this._coef = trimCoef(coef);
    this._domain = domain ?? this.defaultDomain;
    this._window = window ?? this.defaultWindow;
    this._symbol = symbol;

    if (this._coef.length > ABCPolyBase.maxpower + 1) {
      throw new Error(`Polynomial degree exceeds maxpower (${ABCPolyBase.maxpower})`);
    }
  }

  /* ============ Properties ============ */

  get coef(): number[] {
    return [...this._coef];
  }

  get domain(): [number, number] {
    return [...this._domain] as [number, number];
  }

  get window(): [number, number] {
    return [...this._window] as [number, number];
  }

  get degree(): number {
    return this._coef.length - 1;
  }

  get symbol(): string {
    return this._symbol;
  }

  /* ============ Arithmetic Operators ============ */

  add(other: ABCPolyBase | number): this {
    if (typeof other === 'number') {
      const newCoef = [...this._coef];
      newCoef[0] = (newCoef[0] || 0) + other;
      return this._createInstance(newCoef);
    }

    this._checkSameDomain(other);
    const newCoef = this._add(this._coef, other._coef);
    return this._createInstance(newCoef);
  }

  sub(other: ABCPolyBase | number): this {
    if (typeof other === 'number') {
      const newCoef = [...this._coef];
      newCoef[0] = (newCoef[0] || 0) - other;
      return this._createInstance(newCoef);
    }

    this._checkSameDomain(other);
    const newCoef = this._sub(this._coef, other._coef);
    return this._createInstance(newCoef);
  }

  mul(other: ABCPolyBase | number): this {
    if (typeof other === 'number') {
      const newCoef = this._coef.map(c => c * other);
      return this._createInstance(newCoef);
    }

    this._checkSameDomain(other);
    const newCoef = this._mul(this._coef, other._coef);
    return this._createInstance(newCoef);
  }

  div(other: ABCPolyBase | number): [this, this] {
    if (typeof other === 'number') {
      const newCoef = this._coef.map(c => c / other);
      return [this._createInstance(newCoef), this._createInstance([0])];
    }

    this._checkSameDomain(other);
    const [quo, rem] = this._div(this._coef, other._coef);
    return [this._createInstance(quo), this._createInstance(rem)];
  }

  pow(n: number): this {
    if (!Number.isInteger(n) || n < 0) {
      throw new Error('Power must be a non-negative integer');
    }

    if (n === 0) {
      return this._createInstance([1]);
    }

    const newCoef = this._pow(this._coef, n);
    return this._createInstance(newCoef);
  }

  neg(): this {
    return this._createInstance(this._coef.map(c => -c));
  }

  /* ============ Evaluation ============ */

  call(x: number | number[]): number | number[] {
    // Map x from domain to window
    const [off, scl] = this.mapparms();
    const xMapped = Array.isArray(x)
      ? x.map(v => off + scl * v)
      : off + scl * x;

    return this._val(xMapped, this._coef);
  }

  /* ============ Calculus ============ */

  deriv(m: number = 1): this {
    if (m < 0) {
      throw new Error('Derivative order must be non-negative');
    }

    const [off, scl] = this.mapparms();
    const newCoef = this._der(this._coef, m, scl);
    return this._createInstance(newCoef);
  }

  integ(m: number = 1, k: number[] = [], lbnd: number | null = null): this {
    if (m < 0) {
      throw new Error('Integral order must be non-negative');
    }

    const [off, scl] = this.mapparms();
    const effectiveLbnd = lbnd ?? this._domain[0];

    // Pad k with zeros if needed
    const kPadded = [...k];
    while (kPadded.length < m) {
      kPadded.push(0);
    }

    const newCoef = this._int(this._coef, m, kPadded, effectiveLbnd, 1 / scl);
    return this._createInstance(newCoef);
  }

  /* ============ Analysis ============ */

  roots(): number[] {
    return this._roots(this._coef);
  }

  trim(tol: number = 0): this {
    const newCoef = trimCoef(this._coef, tol);
    return this._createInstance(newCoef);
  }

  truncate(size: number): this {
    if (size < 1) {
      throw new Error('Size must be at least 1');
    }
    const newCoef = this._coef.slice(0, size);
    return this._createInstance(newCoef);
  }

  cutdeg(deg: number): this {
    return this.truncate(deg + 1);
  }

  /* ============ Conversion ============ */

  mapparms(): [number, number] {
    // Returns [offset, scale] to map from domain to window
    const [d0, d1] = this._domain;
    const [w0, w1] = this._window;

    const scl = (w1 - w0) / (d1 - d0);
    const off = w0 - d0 * scl;

    return [off, scl];
  }

  convert<T extends ABCPolyBase>(
    TargetClass: new (coef: number[], domain?: [number, number], window?: [number, number]) => T,
    domain: [number, number] | null = null,
    window: [number, number] | null = null
  ): T {
    // Convert to target polynomial type
    // This requires implementing basis conversion (complex)
    throw new Error('Not implemented');
  }

  /* ============ Fitting ============ */

  static fit<T extends ABCPolyBase>(
    this: new (coef: number[], domain?: [number, number], window?: [number, number]) => T,
    x: number[],
    y: number[],
    deg: number,
    domain: [number, number] | null = null,
    rcond: number | null = null,
    full: boolean = false,
    w: number[] | null = null
  ): T | [T, any[]] {
    // Least squares fit using Vandermonde matrix
    // Delegates to specific implementation
    throw new Error('Implemented in subclasses');
  }

  /* ============ Creation Methods ============ */

  static basis<T extends ABCPolyBase>(
    this: new (coef: number[], domain?: [number, number], window?: [number, number]) => T,
    deg: number,
    domain: [number, number] | null = null,
    window: [number, number] | null = null
  ): T {
    const coef = Array(deg + 1).fill(0);
    coef[deg] = 1;
    return new this(coef, domain ?? undefined, window ?? undefined);
  }

  static identity<T extends ABCPolyBase>(
    this: new (coef: number[], domain?: [number, number], window?: [number, number]) => T,
    domain: [number, number] | null = null,
    window: [number, number] | null = null
  ): T {
    return new this([0, 1], domain ?? undefined, window ?? undefined);
  }

  static fromroots<T extends ABCPolyBase>(
    this: new (coef: number[], domain?: [number, number], window?: [number, number]) => T,
    roots: number[],
    domain: [number, number] | null = null,
    window: [number, number] | null = null
  ): T {
    // Construct polynomial from roots
    // Implementation depends on basis type
    throw new Error('Implemented in subclasses');
  }

  copy(): this {
    return this._createInstance([...this._coef]);
  }

  /* ============ Comparison ============ */

  hassamedomain(other: ABCPolyBase): boolean {
    return this._domain[0] === other._domain[0] && this._domain[1] === other._domain[1];
  }

  hassamewindow(other: ABCPolyBase): boolean {
    return this._window[0] === other._window[0] && this._window[1] === other._window[1];
  }

  hassamecoef(other: ABCPolyBase): boolean {
    if (this._coef.length !== other._coef.length) return false;
    return this._coef.every((c, i) => c === other._coef[i]);
  }

  /* ============ String Representation ============ */

  toString(): string {
    const terms: string[] = [];

    for (let i = 0; i < this._coef.length; i++) {
      const c = this._coef[i];
      if (c === 0) continue;

      let term = '';
      if (i === 0) {
        term = c.toString();
      } else {
        const coefStr = c === 1 ? '' : (c === -1 ? '-' : `${c}*`);
        const basisStr = i === 1 ? this._symbol : `${this.basisName}(${i}, ${this._symbol})`;
        term = `${coefStr}${basisStr}`;
      }
      terms.push(term);
    }

    return terms.length > 0 ? terms.join(' + ').replace(/\+ -/g, '- ') : '0';
  }

  /* ============ Helper Methods ============ */

  protected abstract _createInstance(coef: number[]): this;

  protected _checkSameDomain(other: ABCPolyBase): void {
    if (!this.hassamedomain(other) || !this.hassamewindow(other)) {
      throw new Error('Polynomials must have same domain and window');
    }
  }
}

/* ============ Utility Functions ============ */

function trimCoef(coef: number[], tol: number = 0): number[] {
  let end = coef.length;
  while (end > 1 && Math.abs(coef[end - 1]) <= tol) {
    end--;
  }
  return coef.slice(0, end);
}
```

---

### 16.4 numpy.ma Module (Masked Arrays)

Due to the complexity (~9,000 lines in NumPy), I'll provide the core structure.

#### 16.4.1 MaskedArray Class

**File:** `src/ts/ma/core.ts`

```typescript
import { NDArray } from '../NDArray.js';
import { DType } from '../types.js';
import { broadcastArrays, broadcastShapes } from '../broadcast.js';

/** Sentinel for no mask */
export const nomask = Symbol('nomask');

/** Type for mask: boolean array or nomask sentinel */
export type MaskType = NDArray | typeof nomask;

/**
 * MaskedArray: Array with masked (invalid) elements.
 *
 * Masked elements are excluded from computations and
 * can represent missing or invalid data.
 */
export class MaskedArray {
  /** Underlying data array */
  private _data: NDArray;

  /** Boolean mask array (true = masked) */
  private _mask: MaskType;

  /** Fill value for masked elements */
  private _fill_value: number;

  /** If true, mask cannot be changed */
  private _hardmask: boolean;

  /** If true, mask is shared with parent */
  private _sharedmask: boolean;

  /**
   * Create a masked array.
   *
   * @param data - Input data
   * @param mask - Boolean mask (true = masked)
   * @param dtype - Data type
   * @param copy - If true, copy data
   * @param fill_value - Value to use for masked elements
   * @param hard_mask - If true, mask cannot be unset
   * @param shrink - If true, shrink mask to nomask if all false
   */
  constructor(
    data: NDArray | number[],
    mask: NDArray | boolean[] | boolean | typeof nomask = nomask,
    dtype: DType | null = null,
    copy: boolean = false,
    fill_value: number | null = null,
    hard_mask: boolean = false,
    shrink: boolean = true
  ) {
    // Convert data to NDArray
    if (Array.isArray(data)) {
      this._data = NDArray.fromArray(data, dtype ?? DType.Float64);
    } else if (copy) {
      this._data = data.copy();
    } else {
      this._data = data;
    }

    // Apply dtype conversion if specified
    if (dtype !== null && this._data.dtype !== dtype) {
      this._data = this._data.astype(dtype);
    }

    // Process mask
    this._mask = this._processMask(mask, shrink);

    // Set fill value
    this._fill_value = fill_value ?? this._defaultFillValue(this._data.dtype);

    this._hardmask = hard_mask;
    this._sharedmask = false;
  }

  /* ============ Properties ============ */

  get data(): NDArray {
    return this._data;
  }

  get mask(): MaskType {
    return this._mask;
  }

  set mask(value: MaskType) {
    if (this._hardmask && this._mask !== nomask) {
      throw new Error('Cannot modify hard mask');
    }
    this._mask = this._processMask(value, true);
  }

  get fill_value(): number {
    return this._fill_value;
  }

  set fill_value(value: number) {
    this._fill_value = value;
  }

  get shape(): number[] {
    return this._data.shape;
  }

  get ndim(): number {
    return this._data.ndim;
  }

  get size(): number {
    return this._data.size;
  }

  get dtype(): DType {
    return this._data.dtype;
  }

  /* ============ Mask Operations ============ */

  /**
   * Return the number of non-masked elements.
   */
  count(axis: number | null = null): number | NDArray {
    if (this._mask === nomask) {
      if (axis === null) {
        return this._data.size;
      }
      return NDArray.full(
        this._data.shape.filter((_, i) => i !== axis),
        this._data.shape[axis],
        DType.Int32
      );
    }

    // Count non-masked (mask == false)
    const notMask = this._logicalNot(this._mask as NDArray);
    if (axis === null) {
      return notMask.sum() as number;
    }
    return notMask.sum(axis);
  }

  /**
   * Return data with masked values replaced by fill_value.
   */
  filled(fill_value: number | null = null): NDArray {
    const fv = fill_value ?? this._fill_value;

    if (this._mask === nomask) {
      return this._data.copy();
    }

    const result = this._data.copy();
    const mask = this._mask as NDArray;

    for (let i = 0; i < result.size; i++) {
      if (mask.getFlat(i)) {
        result.setFlat(i, fv);
      }
    }

    return result;
  }

  /**
   * Return 1D array of non-masked data.
   */
  compressed(): NDArray {
    if (this._mask === nomask) {
      return this._data.ravel();
    }

    const mask = this._mask as NDArray;
    const values: number[] = [];

    for (let i = 0; i < this._data.size; i++) {
      if (!mask.getFlat(i)) {
        values.push(this._data.getFlat(i));
      }
    }

    return NDArray.fromArray(values, this._data.dtype);
  }

  /**
   * Harden the mask (cannot unmask).
   */
  harden_mask(): this {
    this._hardmask = true;
    return this;
  }

  /**
   * Soften the mask (can unmask).
   */
  soften_mask(): this {
    this._hardmask = false;
    return this;
  }

  /* ============ Arithmetic Operations ============ */

  add(other: MaskedArray | NDArray | number): MaskedArray {
    return this._binaryOp(other, (a, b) => a + b);
  }

  subtract(other: MaskedArray | NDArray | number): MaskedArray {
    return this._binaryOp(other, (a, b) => a - b);
  }

  multiply(other: MaskedArray | NDArray | number): MaskedArray {
    return this._binaryOp(other, (a, b) => a * b);
  }

  divide(other: MaskedArray | NDArray | number): MaskedArray {
    return this._binaryOp(other, (a, b) => a / b);
  }

  power(other: MaskedArray | NDArray | number): MaskedArray {
    return this._binaryOp(other, (a, b) => Math.pow(a, b));
  }

  /* ============ Reductions ============ */

  sum(axis: number | null = null, keepdims: boolean = false): number | MaskedArray {
    return this._reduction((data) => data.sum(axis, keepdims), axis, keepdims);
  }

  mean(axis: number | null = null, keepdims: boolean = false): number | MaskedArray {
    const sumResult = this.sum(axis, keepdims);
    const countResult = this.count(axis);

    if (typeof sumResult === 'number' && typeof countResult === 'number') {
      return sumResult / countResult;
    }

    // Element-wise division
    return this._divideArrays(sumResult as MaskedArray, countResult as NDArray);
  }

  std(axis: number | null = null, ddof: number = 0, keepdims: boolean = false): number | MaskedArray {
    const variance = this.var(axis, ddof, keepdims);
    if (typeof variance === 'number') {
      return Math.sqrt(variance);
    }
    return this._sqrtArray(variance);
  }

  var(axis: number | null = null, ddof: number = 0, keepdims: boolean = false): number | MaskedArray {
    const mean = this.mean(axis, true);
    const centered = this.subtract(mean as MaskedArray);
    const squared = centered.multiply(centered);
    const sumSquared = squared.sum(axis, keepdims);
    const count = this.count(axis);

    if (typeof sumSquared === 'number' && typeof count === 'number') {
      return sumSquared / (count - ddof);
    }

    // Element-wise division with ddof adjustment
    return this._divideArraysWithDdof(sumSquared as MaskedArray, count as NDArray, ddof);
  }

  min(axis: number | null = null, keepdims: boolean = false): number | MaskedArray {
    return this._reduction((data) => data.min(axis, keepdims), axis, keepdims);
  }

  max(axis: number | null = null, keepdims: boolean = false): number | MaskedArray {
    return this._reduction((data) => data.max(axis, keepdims), axis, keepdims);
  }

  /* ============ Shape Manipulation ============ */

  reshape(newshape: number[]): MaskedArray {
    const newData = this._data.reshape(newshape);
    const newMask = this._mask === nomask
      ? nomask
      : (this._mask as NDArray).reshape(newshape);

    return new MaskedArray(newData, newMask, null, false, this._fill_value, this._hardmask);
  }

  ravel(): MaskedArray {
    return this.reshape([this._data.size]);
  }

  flatten(): MaskedArray {
    return new MaskedArray(
      this._data.flatten(),
      this._mask === nomask ? nomask : (this._mask as NDArray).flatten(),
      null, false, this._fill_value, this._hardmask
    );
  }

  transpose(axes: number[] | null = null): MaskedArray {
    const newData = this._data.transpose(axes);
    const newMask = this._mask === nomask
      ? nomask
      : (this._mask as NDArray).transpose(axes);

    return new MaskedArray(newData, newMask, null, false, this._fill_value, this._hardmask);
  }

  /* ============ Helper Methods ============ */

  private _processMask(mask: NDArray | boolean[] | boolean | typeof nomask, shrink: boolean): MaskType {
    if (mask === nomask) {
      return nomask;
    }

    if (typeof mask === 'boolean') {
      if (!mask) {
        return nomask;
      }
      // All masked
      return NDArray.full(this._data.shape, 1, DType.Bool);
    }

    let maskArr: NDArray;
    if (Array.isArray(mask)) {
      maskArr = NDArray.fromArray(mask, DType.Bool);
    } else {
      maskArr = mask.astype(DType.Bool);
    }

    // Broadcast mask to data shape if needed
    if (!arraysEqual(maskArr.shape, this._data.shape)) {
      const [broadcastedMask] = broadcastArrays(maskArr, this._data);
      maskArr = broadcastedMask;
    }

    // Shrink to nomask if all false
    if (shrink && !this._anyTrue(maskArr)) {
      return nomask;
    }

    return maskArr;
  }

  private _defaultFillValue(dtype: DType): number {
    const defaults: Partial<Record<DType, number>> = {
      [DType.Bool]: true as unknown as number,
      [DType.Int8]: 999999,
      [DType.Int16]: 999999,
      [DType.Int32]: 999999,
      [DType.Int64]: 999999,
      [DType.UInt8]: 999999,
      [DType.UInt16]: 999999,
      [DType.UInt32]: 999999,
      [DType.UInt64]: 999999,
      [DType.Float16]: 1e20,
      [DType.Float32]: 1e20,
      [DType.Float64]: 1e20,
    };
    return defaults[dtype] ?? 1e20;
  }

  private _binaryOp(
    other: MaskedArray | NDArray | number,
    op: (a: number, b: number) => number
  ): MaskedArray {
    let otherData: NDArray;
    let otherMask: MaskType;

    if (typeof other === 'number') {
      otherData = NDArray.full(this._data.shape, other, this._data.dtype);
      otherMask = nomask;
    } else if (other instanceof MaskedArray) {
      otherData = other._data;
      otherMask = other._mask;
    } else {
      otherData = other;
      otherMask = nomask;
    }

    // Broadcast arrays
    const [bData1, bData2] = broadcastArrays(this._data, otherData);

    // Compute result
    const resultData = NDArray.empty(bData1.shape, this._data.dtype);
    for (let i = 0; i < resultData.size; i++) {
      resultData.setFlat(i, op(bData1.getFlat(i), bData2.getFlat(i)));
    }

    // Combine masks (OR)
    const resultMask = this._combineMasks(this._mask, otherMask, bData1.shape);

    return new MaskedArray(resultData, resultMask, null, false, this._fill_value);
  }

  private _combineMasks(mask1: MaskType, mask2: MaskType, shape: number[]): MaskType {
    if (mask1 === nomask && mask2 === nomask) {
      return nomask;
    }

    const m1 = mask1 === nomask ? NDArray.zeros(shape, DType.Bool) : mask1 as NDArray;
    const m2 = mask2 === nomask ? NDArray.zeros(shape, DType.Bool) : mask2 as NDArray;

    // Broadcast masks
    const [bm1, bm2] = broadcastArrays(m1, m2);

    // OR the masks
    const result = NDArray.empty(bm1.shape, DType.Bool);
    for (let i = 0; i < result.size; i++) {
      result.setFlat(i, bm1.getFlat(i) || bm2.getFlat(i) ? 1 : 0);
    }

    return result;
  }

  private _reduction(
    reductionFn: (data: NDArray) => number | NDArray,
    axis: number | null,
    keepdims: boolean
  ): number | MaskedArray {
    // Replace masked values with appropriate neutral value before reduction
    const filledData = this.filled(this._neutralValue());
    const result = reductionFn(filledData);

    if (typeof result === 'number') {
      return result;
    }

    // For axis reductions, need to track which results are fully masked
    return new MaskedArray(result as NDArray);
  }

  private _neutralValue(): number {
    // Return neutral value for the most common reduction (sum)
    return 0;
  }

  private _anyTrue(arr: NDArray): boolean {
    for (let i = 0; i < arr.size; i++) {
      if (arr.getFlat(i)) return true;
    }
    return false;
  }

  private _logicalNot(arr: NDArray): NDArray {
    const result = NDArray.empty(arr.shape, DType.Bool);
    for (let i = 0; i < arr.size; i++) {
      result.setFlat(i, arr.getFlat(i) ? 0 : 1);
    }
    return result;
  }

  private _divideArrays(a: MaskedArray, b: NDArray): MaskedArray {
    // Element-wise division
    return a.divide(b);
  }

  private _divideArraysWithDdof(a: MaskedArray, b: NDArray, ddof: number): MaskedArray {
    // Subtract ddof from count before division
    const adjusted = NDArray.empty(b.shape, b.dtype);
    for (let i = 0; i < b.size; i++) {
      adjusted.setFlat(i, b.getFlat(i) - ddof);
    }
    return a.divide(adjusted);
  }

  private _sqrtArray(a: MaskedArray): MaskedArray {
    const result = NDArray.empty(a.shape, a.dtype);
    for (let i = 0; i < a.size; i++) {
      result.setFlat(i, Math.sqrt(a._data.getFlat(i)));
    }
    return new MaskedArray(result, a._mask);
  }
}

/* ============ Utility Functions ============ */

function arraysEqual(a: number[], b: number[]): boolean {
  if (a.length !== b.length) return false;
  return a.every((v, i) => v === b[i]);
}
```

---

### 16.5 numpy.testing Module

**File:** `src/ts/testing/assertions.ts`

```typescript
import { NDArray } from '../NDArray.js';
import { DType } from '../types.js';

/**
 * Assert that two objects are equal.
 *
 * @param actual - Actual value
 * @param desired - Expected value
 * @param err_msg - Error message on failure
 * @param verbose - If true, show detailed diff
 */
export function assert_equal(
  actual: any,
  desired: any,
  err_msg: string = '',
  verbose: boolean = true
): void {
  if (actual instanceof NDArray && desired instanceof NDArray) {
    assert_array_equal(actual, desired, err_msg, verbose);
    return;
  }

  if (Array.isArray(actual) && Array.isArray(desired)) {
    if (actual.length !== desired.length) {
      throw new AssertionError(
        buildErrorMessage(actual, desired, 'Arrays have different lengths', verbose)
      );
    }
    for (let i = 0; i < actual.length; i++) {
      assert_equal(actual[i], desired[i], err_msg, verbose);
    }
    return;
  }

  if (!Object.is(actual, desired)) {
    // Handle NaN comparison
    if (typeof actual === 'number' && typeof desired === 'number' &&
        Number.isNaN(actual) && Number.isNaN(desired)) {
      return; // NaN == NaN for testing
    }

    throw new AssertionError(
      buildErrorMessage(actual, desired, err_msg || 'Values are not equal', verbose)
    );
  }
}

/**
 * Assert that two arrays are element-wise equal.
 */
export function assert_array_equal(
  x: NDArray | number[],
  y: NDArray | number[],
  err_msg: string = '',
  verbose: boolean = true
): void {
  const arr1 = Array.isArray(x) ? NDArray.fromArray(x) : x;
  const arr2 = Array.isArray(y) ? NDArray.fromArray(y) : y;

  // Check shapes
  if (!arraysEqual(arr1.shape, arr2.shape)) {
    throw new AssertionError(
      `Arrays have different shapes: ${arr1.shape} vs ${arr2.shape}`
    );
  }

  // Check values
  for (let i = 0; i < arr1.size; i++) {
    const v1 = arr1.getFlat(i);
    const v2 = arr2.getFlat(i);

    // Handle NaN
    if (Number.isNaN(v1) && Number.isNaN(v2)) {
      continue;
    }

    if (v1 !== v2) {
      throw new AssertionError(
        buildErrorMessage(arr1, arr2, err_msg || `Arrays differ at index ${i}`, verbose)
      );
    }
  }
}

/**
 * Assert that two arrays are element-wise equal within tolerance.
 *
 * |actual - desired| <= atol + rtol * |desired|
 *
 * @param actual - Actual array
 * @param desired - Expected array
 * @param rtol - Relative tolerance (default 1e-7)
 * @param atol - Absolute tolerance (default 0)
 * @param equal_nan - If true, NaN == NaN
 */
export function assert_allclose(
  actual: NDArray | number[] | number,
  desired: NDArray | number[] | number,
  rtol: number = 1e-7,
  atol: number = 0,
  equal_nan: boolean = true,
  err_msg: string = '',
  verbose: boolean = true
): void {
  const arr1 = toNDArray(actual);
  const arr2 = toNDArray(desired);

  // Check shapes
  if (!arraysEqual(arr1.shape, arr2.shape)) {
    throw new AssertionError(
      `Arrays have different shapes: ${arr1.shape} vs ${arr2.shape}`
    );
  }

  const mismatches: number[] = [];

  for (let i = 0; i < arr1.size; i++) {
    const v1 = arr1.getFlat(i);
    const v2 = arr2.getFlat(i);

    // Handle NaN
    if (Number.isNaN(v1) && Number.isNaN(v2)) {
      if (equal_nan) continue;
      mismatches.push(i);
      continue;
    }

    // Handle Inf
    if (!Number.isFinite(v1) || !Number.isFinite(v2)) {
      if (v1 === v2) continue;
      mismatches.push(i);
      continue;
    }

    // Check tolerance
    const diff = Math.abs(v1 - v2);
    const tolerance = atol + rtol * Math.abs(v2);

    if (diff > tolerance) {
      mismatches.push(i);
    }
  }

  if (mismatches.length > 0) {
    const maxDiff = mismatches.reduce((max, i) => {
      const diff = Math.abs(arr1.getFlat(i) - arr2.getFlat(i));
      return Math.max(max, diff);
    }, 0);

    throw new AssertionError(
      `Arrays are not close within tolerance.\n` +
      `Max absolute difference: ${maxDiff}\n` +
      `Number of mismatches: ${mismatches.length} / ${arr1.size}\n` +
      (err_msg ? `\n${err_msg}` : '')
    );
  }
}

/**
 * Assert that two values are approximately equal.
 *
 * @param actual - Actual value
 * @param desired - Expected value
 * @param decimal - Number of decimal places to compare
 */
export function assert_almost_equal(
  actual: number | NDArray,
  desired: number | NDArray,
  decimal: number = 7,
  err_msg: string = '',
  verbose: boolean = true
): void {
  const atol = Math.pow(10, -decimal);

  if (typeof actual === 'number' && typeof desired === 'number') {
    if (Math.abs(actual - desired) > atol) {
      throw new AssertionError(
        `Values differ by more than ${decimal} decimal places: ` +
        `${actual} vs ${desired}`
      );
    }
    return;
  }

  assert_allclose(actual, desired, 0, atol, true, err_msg, verbose);
}

/**
 * Assert that two values are approximately equal to significant figures.
 */
export function assert_approx_equal(
  actual: number,
  desired: number,
  significant: number = 7,
  err_msg: string = '',
  verbose: boolean = true
): void {
  if (desired === 0) {
    if (actual !== 0) {
      throw new AssertionError(`Expected 0, got ${actual}`);
    }
    return;
  }

  const scale = Math.pow(10, Math.floor(Math.log10(Math.abs(desired))));
  const rtol = Math.pow(10, -significant);

  if (Math.abs(actual - desired) / scale > rtol) {
    throw new AssertionError(
      `Values differ by more than ${significant} significant figures: ` +
      `${actual} vs ${desired}`
    );
  }
}

/**
 * Assert that x < y element-wise.
 */
export function assert_array_less(
  x: NDArray | number[] | number,
  y: NDArray | number[] | number,
  err_msg: string = '',
  verbose: boolean = true
): void {
  const arr1 = toNDArray(x);
  const arr2 = toNDArray(y);

  for (let i = 0; i < arr1.size; i++) {
    const v1 = arr1.getFlat(i);
    const v2 = arr2.getFlat(i);

    if (!(v1 < v2)) {
      throw new AssertionError(
        `Arrays are not less: x[${i}] = ${v1} >= y[${i}] = ${v2}` +
        (err_msg ? `\n${err_msg}` : '')
      );
    }
  }
}

/**
 * Assert that a callable raises an exception.
 */
export function assert_raises<T extends Error>(
  exception_class: new (...args: any[]) => T,
  callable: (...args: any[]) => any,
  ...args: any[]
): void {
  try {
    callable(...args);
  } catch (e) {
    if (e instanceof exception_class) {
      return;
    }
    throw new AssertionError(
      `Expected ${exception_class.name}, but got ${(e as Error).constructor.name}: ${(e as Error).message}`
    );
  }

  throw new AssertionError(`Expected ${exception_class.name} to be raised, but no exception was thrown`);
}

/**
 * Assert that a callable raises an exception with matching message.
 */
export function assert_raises_regex<T extends Error>(
  exception_class: new (...args: any[]) => T,
  expected_regexp: string | RegExp,
  callable: (...args: any[]) => any,
  ...args: any[]
): void {
  const regex = typeof expected_regexp === 'string'
    ? new RegExp(expected_regexp)
    : expected_regexp;

  try {
    callable(...args);
  } catch (e) {
    if (e instanceof exception_class) {
      if (regex.test((e as Error).message)) {
        return;
      }
      throw new AssertionError(
        `Exception message "${(e as Error).message}" does not match "${regex}"`
      );
    }
    throw new AssertionError(
      `Expected ${exception_class.name}, but got ${(e as Error).constructor.name}`
    );
  }

  throw new AssertionError(`Expected ${exception_class.name} to be raised`);
}

/**
 * Assert that no warnings are raised.
 */
export function assert_no_warnings(
  func: (...args: any[]) => any,
  ...args: any[]
): any {
  // In JavaScript, we don't have a built-in warning system like Python
  // This would need to integrate with a custom warning system
  return func(...args);
}

/**
 * Basic assertion that works in optimized mode.
 */
export function assert_(val: any, msg: string = ''): void {
  if (!val) {
    throw new AssertionError(msg || 'Assertion failed');
  }
}

/* ============ Error Classes ============ */

export class AssertionError extends Error {
  constructor(message: string) {
    super(message);
    this.name = 'AssertionError';
  }
}

export class SkipTest extends Error {
  constructor(message: string) {
    super(message);
    this.name = 'SkipTest';
  }
}

export class KnownFailureException extends Error {
  constructor(message: string) {
    super(message);
    this.name = 'KnownFailureException';
  }
}

/* ============ Helper Functions ============ */

function toNDArray(x: NDArray | number[] | number): NDArray {
  if (x instanceof NDArray) return x;
  if (Array.isArray(x)) return NDArray.fromArray(x);
  return NDArray.fromArray([x]);
}

function arraysEqual(a: number[], b: number[]): boolean {
  if (a.length !== b.length) return false;
  return a.every((v, i) => v === b[i]);
}

function buildErrorMessage(
  actual: any,
  desired: any,
  header: string,
  verbose: boolean
): string {
  if (!verbose) {
    return header;
  }

  return `${header}\n` +
    `Actual:\n${formatValue(actual)}\n` +
    `Desired:\n${formatValue(desired)}`;
}

function formatValue(val: any): string {
  if (val instanceof NDArray) {
    return val.toString();
  }
  if (Array.isArray(val)) {
    return JSON.stringify(val);
  }
  return String(val);
}
```

---

## File Changes Summary

### New Files to Create

```
src/ts/
├── strings/
│   ├── index.ts           # Public exports
│   ├── compare.ts         # Comparison functions (6)
│   ├── properties.ts      # Property testing (10)
│   ├── search.ts          # Search functions (7)
│   └── manipulation.ts    # Manipulation functions (20+)
│
├── rec/
│   ├── index.ts           # Public exports and creation functions
│   ├── format_parser.ts   # Format string parsing
│   └── recarray.ts        # record and recarray classes
│
├── polynomial/
│   ├── index.ts           # Public exports
│   ├── _polybase.ts       # Abstract base class (~400 lines)
│   ├── polyutils.ts       # Utility functions (~200 lines)
│   ├── polynomial.ts      # Polynomial class (~600 lines)
│   ├── chebyshev.ts       # Chebyshev class (~800 lines)
│   ├── legendre.ts        # Legendre class (~700 lines)
│   ├── hermite.ts         # Hermite class (~700 lines)
│   ├── hermite_e.ts       # HermiteE class (~700 lines)
│   └── laguerre.ts        # Laguerre class (~700 lines)
│
├── ma/
│   ├── index.ts           # Public exports
│   ├── core.ts            # MaskedArray class (~2000 lines)
│   ├── extras.ts          # Extra utility functions (~1200 lines)
│   └── mrecords.ts        # Masked record arrays (~300 lines)
│
└── testing/
    ├── index.ts           # Public exports
    ├── assertions.ts      # Assertion functions (~800 lines)
    ├── exceptions.ts      # Exception classes (~100 lines)
    └── utils.ts           # Utility functions (~400 lines)

src/wasm/ (optional)
├── strings.h/c            # WASM string operations
└── polynomial.h/c         # WASM polynomial evaluation
```

### Files to Modify

```
src/ts/types.ts
├── Add StructuredDType interface
├── Add FieldDescriptor interface
├── Add string dtype support
└── Add mask-related types

src/ts/NDArray.ts
├── Add getStringFlat, setStringFlat methods
├── Add getField method for structured arrays
├── Add emptyString, emptyStructured factory methods
└── Add string array support

src/ts/index.ts
├── Export strings module
├── Export rec module
├── Export polynomial module
├── Export ma module
└── Export testing module
```

---

## Implementation Order

```
Phase 16.1: numpy.strings (Week 1-2)
├── Week 1: Core Functions
│   ├── Day 1: String array infrastructure in NDArray
│   ├── Day 2: Comparison functions (equal, not_equal, etc.)
│   ├── Day 3: Property testing (isalpha, isdigit, etc.)
│   ├── Day 4: Search functions (find, rfind, index, count)
│   └── Day 5: Basic manipulation (upper, lower, strip)
│
└── Week 2: Advanced Functions
    ├── Day 1: replace, center, ljust, rjust, zfill
    ├── Day 2: partition, rpartition, expandtabs
    ├── Day 3: add, multiply for strings
    ├── Day 4: encode, decode
    └── Day 5: Tests and documentation

Phase 16.2: numpy.rec (Week 3)
├── Day 1: StructuredDType and FieldDescriptor types
├── Day 2: format_parser class
├── Day 3: record class
├── Day 4: recarray class
└── Day 5: Creation functions + tests

Phase 16.3: numpy.polynomial (Weeks 4-6)
├── Week 4: Foundation
│   ├── Day 1: polyutils (trimcoef, as_series, etc.)
│   ├── Day 2: ABCPolyBase structure
│   ├── Day 3: ABCPolyBase arithmetic
│   ├── Day 4: ABCPolyBase calculus
│   └── Day 5: ABCPolyBase fitting
│
├── Week 5: Polynomial Class
│   ├── Day 1: polyadd, polysub, polymul
│   ├── Day 2: polydiv, polypow
│   ├── Day 3: polyval, polyval2d, polyval3d
│   ├── Day 4: polyder, polyint
│   └── Day 5: polyfit, polyvander, polyroots
│
└── Week 6: Other Bases
    ├── Day 1: Chebyshev (most complex)
    ├── Day 2: Legendre
    ├── Day 3: Hermite, HermiteE
    ├── Day 4: Laguerre
    └── Day 5: Integration tests

Phase 16.4: numpy.ma (Weeks 7-9)
├── Week 7: MaskedArray Core
│   ├── Day 1: MaskedArray class structure
│   ├── Day 2: Mask handling (nomask, make_mask, etc.)
│   ├── Day 3: filled, compressed, count
│   ├── Day 4: Basic arithmetic with mask propagation
│   └── Day 5: Comparison operators
│
├── Week 8: Operations
│   ├── Day 1: Math functions (sin, cos, exp, log)
│   ├── Day 2: Reductions (sum, mean, std, var)
│   ├── Day 3: Shape manipulation
│   ├── Day 4: Indexing with masks
│   └── Day 5: masked_equal, masked_where, etc.
│
└── Week 9: Extras
    ├── Day 1: average, median
    ├── Day 2: cov, corrcoef
    ├── Day 3: Set operations with masks
    ├── Day 4: mrecords module
    └── Day 5: Integration tests

Phase 16.5: numpy.testing (Week 10)
├── Day 1: assert_equal, assert_array_equal
├── Day 2: assert_allclose, assert_almost_equal
├── Day 3: assert_raises, assert_warns
├── Day 4: Context managers, utilities
└── Day 5: Integration with test framework
```

---

## Verification Plan

After Phase 16 completion, verify:

```bash
# Build
npm run build

# Run all tests
npm test

# Phase 16 specific tests:

# Strings
✓ strings.equal works with string arrays
✓ strings.upper/lower transform correctly
✓ strings.find returns correct indices
✓ strings.replace handles count parameter
✓ strings.isalpha returns correct boolean array

# Record Arrays
✓ format_parser creates correct dtype
✓ recarray allows field access as attributes
✓ rec.fromarrays creates from column arrays
✓ rec.fromrecords creates from row tuples
✓ Binary I/O works correctly

# Polynomial
✓ Polynomial arithmetic is correct
✓ Polynomial.fit returns least-squares fit
✓ Chebyshev conversions are accurate
✓ All bases have correct derivatives/integrals
✓ Root finding returns correct roots

# Masked Arrays
✓ MaskedArray propagates masks through operations
✓ Reductions respect masks (exclude masked values)
✓ compressed() returns only non-masked values
✓ filled() replaces masked with fill_value
✓ Hard masks cannot be unset

# Testing
✓ assert_equal works for arrays and scalars
✓ assert_allclose respects tolerances
✓ assert_raises catches correct exceptions
✓ NaN handling is correct in comparisons
```

Generate NumPy comparison vectors:

```python
# tests/python/generate_phase16_tests.py
import numpy as np
import json

tests = {
    "strings": {
        "equal": [
            {"a": ["hello", "world"], "b": ["hello", "test"],
             "result": [True, False]},
        ],
        "upper": [
            {"a": ["hello", "World"], "result": ["HELLO", "WORLD"]},
        ],
        "find": [
            {"a": ["hello world", "test"], "sub": "o",
             "result": [4, -1]},
        ],
    },
    "polynomial": {
        "polyval": [
            {"c": [1, 2, 3], "x": [0, 1, 2],
             "result": [1, 6, 17]},  # 1 + 2x + 3x^2
        ],
        "polyder": [
            {"c": [1, 2, 3], "result": [2, 6]},  # d/dx(1 + 2x + 3x^2) = 2 + 6x
        ],
    },
    "ma": {
        "mean": [
            {"data": [1, 2, 3, 4], "mask": [False, False, True, False],
             "result": 2.333333},  # (1+2+4)/3
        ],
        "sum": [
            {"data": [1, 2, 3, 4], "mask": [False, True, False, False],
             "result": 8},  # 1+3+4
        ],
    },
}

with open("tests/fixtures/phase16_vectors.json", "w") as f:
    json.dump(tests, f, indent=2)
```

---

## Complexity Summary

| Module | Estimated Lines | Difficulty | Key Challenges |
|--------|----------------|------------|----------------|
| **strings** | ~1,200 | Low-Medium | Unicode handling, variable-length strings |
| **rec** | ~600 | Low | Structured dtype integration |
| **polynomial** | ~4,500 | High | Basis conversions, numerical stability |
| **ma** | ~3,500 | Very High | Mask propagation, ndarray semantics |
| **testing** | ~1,300 | Medium | Clear error messages, tolerance handling |

**Total Estimated: ~11,100 lines TypeScript**

---

## API Compatibility Notes

### String Arrays
```typescript
// NumPy:
// np.char.equal(a, b)

// NumJS:
strings.equal(a, b);
```

### Polynomial
```typescript
// NumPy:
// p = np.polynomial.Polynomial([1, 2, 3])
// p(x)

// NumJS:
const p = new Polynomial([1, 2, 3]);
p.call(x);  // or p.evaluate(x)
```

### Masked Arrays
```typescript
// NumPy:
// ma = np.ma.array([1, 2, 3], mask=[False, True, False])
// ma.mean()

// NumJS:
const ma = new MaskedArray([1, 2, 3], [false, true, false]);
ma.mean();
```

### Testing
```typescript
// NumPy:
// np.testing.assert_allclose(a, b, rtol=1e-5)

// NumJS:
testing.assert_allclose(a, b, 1e-5);
```
