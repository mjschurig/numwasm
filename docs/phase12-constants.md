# Phase 12: Constants Implementation Plan

Phase 12 implements NumPy's mathematical and special constants for NumJS-WASM. These constants provide convenient access to commonly used mathematical values and special floating-point values.

---

## Current State (Phases 1-11 Complete)

```
src/ts/
├── index.ts           # Main exports (300+ exports)
├── NDArray.ts         # Core array class
├── types.ts           # DType enum, type definitions
├── slice.ts           # Slice, ellipsis, newaxis (for indexing)
├── broadcast.ts       # Broadcasting functions
├── indexing.ts        # Index functions
├── manipulation.ts    # Array manipulation
├── statistics.ts      # Statistics functions
├── ufunc.ts           # Universal functions
├── logic.ts           # Logical operations
├── sorting.ts         # Sorting functions
├── functional.ts      # Functional programming utilities
└── ...
```

**Note:** The `newaxis` symbol already exists in `slice.ts` for array indexing purposes. Phase 12 will re-export it as part of the constants module for API consistency with NumPy.

---

## Phase 12 Implementation Tree

```
PHASE 12: CONSTANTS
│
├── 12.1 Mathematical Constants (TypeScript)
│   ├── 12.1.1 e → Euler's number (2.718281828459045)
│   ├── 12.1.2 pi → Pi (3.141592653589793)
│   └── 12.1.3 euler_gamma → Euler-Mascheroni constant (0.5772156649015329)
│
├── 12.2 Special Floating-Point Values (TypeScript)
│   ├── 12.2.1 inf → Positive infinity
│   ├── 12.2.2 PINF → Alias for positive infinity
│   ├── 12.2.3 NINF → Negative infinity
│   ├── 12.2.4 nan → Not a Number (quiet NaN)
│   └── 12.2.5 NAN → Alias for nan
│
├── 12.3 Indexing Constants (Re-export)
│   └── 12.3.1 newaxis → Re-export from slice.ts (None/null equivalent)
│
└── 12.4 Additional NumPy Constants (TypeScript)
    ├── 12.4.1 PZERO → Positive zero (+0.0)
    ├── 12.4.2 NZERO → Negative zero (-0.0)
    └── 12.4.3 finfo/iinfo types → Float/integer type information
```

---

## Detailed Implementation Specifications

### 12.1 Mathematical Constants

**File:** `src/ts/constants.ts` (new file)

```typescript
/**
 * NumJS Constants
 *
 * Mathematical and special constants matching NumPy's constant exports.
 * Reference: numpy/_core/numeric.py, numpy/_core/umath.py
 */

// =============================================================================
// Mathematical Constants
// =============================================================================

/**
 * Euler's number, the base of natural logarithms.
 *
 * `e = 2.718281828459045`
 *
 * @example
 * import { e, exp } from 'numjs';
 * console.log(e);  // 2.718281828459045
 * console.log(exp(fromArray([1]))); // e^1 ≈ 2.718
 *
 * @see https://numpy.org/doc/stable/reference/constants.html#numpy.e
 */
export const e: number = 2.718281828459045;

/**
 * Pi, the ratio of a circle's circumference to its diameter.
 *
 * `π = 3.141592653589793`
 *
 * @example
 * import { pi, sin } from 'numjs';
 * console.log(pi);  // 3.141592653589793
 * console.log(sin(fromArray([pi / 2]))); // sin(π/2) = 1
 *
 * @see https://numpy.org/doc/stable/reference/constants.html#numpy.pi
 */
export const pi: number = 3.141592653589793;

/**
 * Euler-Mascheroni constant (γ).
 *
 * `γ = lim(n→∞) [1 + 1/2 + 1/3 + ... + 1/n - ln(n)]`
 * `γ ≈ 0.5772156649015329`
 *
 * @example
 * import { euler_gamma } from 'numjs';
 * console.log(euler_gamma);  // 0.5772156649015329
 *
 * @see https://numpy.org/doc/stable/reference/constants.html#numpy.euler_gamma
 */
export const euler_gamma: number = 0.5772156649015329;
```

### 12.2 Special Floating-Point Values

**File:** `src/ts/constants.ts` (continued)

```typescript
// =============================================================================
// Special Floating-Point Values
// =============================================================================

/**
 * IEEE 754 positive infinity.
 *
 * @example
 * import { inf, isinf, isposinf } from 'numjs';
 * console.log(inf);           // Infinity
 * console.log(1 / 0 === inf); // true
 * console.log(inf + 1 === inf); // true
 *
 * // Use with arrays
 * const arr = fromArray([1, inf, -inf, 0]);
 * console.log(isinf(arr).toArray());  // [false, true, true, false]
 * console.log(isposinf(arr).toArray()); // [false, true, false, false]
 *
 * @see https://numpy.org/doc/stable/reference/constants.html#numpy.inf
 */
export const inf: number = Infinity;

/**
 * IEEE 754 positive infinity (alias).
 *
 * @see inf
 */
export const PINF: number = Infinity;

/**
 * IEEE 754 negative infinity.
 *
 * @example
 * import { NINF, isneginf } from 'numjs';
 * console.log(NINF);            // -Infinity
 * console.log(-1 / 0 === NINF); // true
 *
 * @see https://numpy.org/doc/stable/reference/constants.html#numpy.NINF
 */
export const NINF: number = -Infinity;

/**
 * IEEE 754 quiet Not-a-Number.
 *
 * NaN is used to represent undefined or unrepresentable results
 * in floating-point calculations (e.g., 0/0, sqrt(-1)).
 *
 * Note: NaN is not equal to anything, including itself.
 * Use `isnan()` to check for NaN values.
 *
 * @example
 * import { nan, isnan } from 'numjs';
 * console.log(nan);              // NaN
 * console.log(nan === nan);      // false (NaN never equals NaN)
 * console.log(Number.isNaN(nan)); // true
 *
 * // Use with arrays
 * const arr = fromArray([1, nan, 3]);
 * console.log(isnan(arr).toArray()); // [false, true, false]
 *
 * @see https://numpy.org/doc/stable/reference/constants.html#numpy.nan
 */
export const nan: number = NaN;

/**
 * IEEE 754 quiet Not-a-Number (alias).
 *
 * @see nan
 */
export const NAN: number = NaN;

/**
 * IEEE 754 positive zero.
 *
 * JavaScript/IEEE 754 distinguishes between +0 and -0.
 * They compare equal, but have different behavior in edge cases.
 *
 * @example
 * import { PZERO, NZERO } from 'numjs';
 * console.log(PZERO);          // 0
 * console.log(NZERO);          // -0
 * console.log(PZERO === NZERO); // true (compare equal)
 * console.log(1 / PZERO);       // Infinity
 * console.log(1 / NZERO);       // -Infinity
 *
 * @see https://numpy.org/doc/stable/reference/constants.html#numpy.PZERO
 */
export const PZERO: number = 0.0;

/**
 * IEEE 754 negative zero.
 *
 * @see PZERO for examples
 * @see https://numpy.org/doc/stable/reference/constants.html#numpy.NZERO
 */
export const NZERO: number = -0.0;
```

### 12.3 Indexing Constants

**File:** `src/ts/constants.ts` (continued)

```typescript
// =============================================================================
// Indexing Constants
// =============================================================================

// Re-export newaxis from slice.ts for API consistency
// newaxis is already defined as a unique symbol in slice.ts
// It's used in array indexing to add a new dimension

/**
 * Alias for newaxis, used to expand array dimensions.
 *
 * When used in an index expression, adds a new axis of length 1
 * at that position.
 *
 * @example
 * import { newaxis, fromArray, slice } from 'numjs';
 *
 * const arr = fromArray([1, 2, 3]); // shape: [3]
 * const row = arr.slice([newaxis, slice()]); // shape: [1, 3]
 * const col = arr.slice([slice(), newaxis]); // shape: [3, 1]
 *
 * // Useful for broadcasting:
 * // If a has shape (3,) and b has shape (4,):
 * // a[:, newaxis] has shape (3, 1)
 * // a[:, newaxis] + b produces shape (3, 4) via broadcasting
 *
 * @see https://numpy.org/doc/stable/reference/constants.html#numpy.newaxis
 */
export { newaxis } from './slice.js';
```

### 12.4 Type Information Classes

**File:** `src/ts/typeinfo.ts` (new file)

```typescript
/**
 * Type information utilities for NumJS
 *
 * Provides finfo and iinfo classes for inspecting floating-point
 * and integer type properties, similar to NumPy's numpy.finfo and numpy.iinfo.
 */

import { DType, DTYPE_SIZES } from './types.js';

// =============================================================================
// Floating-Point Type Information
// =============================================================================

/**
 * Machine limits for floating-point types.
 *
 * Provides properties describing the numerical limits and precision
 * of floating-point data types.
 *
 * @example
 * import { finfo, DType } from 'numjs';
 *
 * const f32 = finfo(DType.Float32);
 * console.log(f32.eps);     // Machine epsilon
 * console.log(f32.min);     // Smallest positive normal
 * console.log(f32.max);     // Largest representable number
 * console.log(f32.bits);    // Number of bits
 *
 * @see https://numpy.org/doc/stable/reference/generated/numpy.finfo.html
 */
export class finfo {
  /** Number of bits in the type */
  readonly bits: number;

  /** Machine epsilon: smallest x such that 1.0 + x != 1.0 */
  readonly eps: number;

  /** The largest usable exponent for this type */
  readonly maxexp: number;

  /** The smallest usable exponent for this type */
  readonly minexp: number;

  /** The largest representable number */
  readonly max: number;

  /** The smallest positive usable number (normal) */
  readonly min: number;

  /** The smallest positive number (subnormal) */
  readonly tiny: number;

  /** Number of decimal digits of precision */
  readonly precision: number;

  /** Number of significant digits in base 10 */
  readonly iexp: number;

  /** Machine parameters for the dtype */
  readonly dtype: DType;

  constructor(dtype: DType | 'float16' | 'float32' | 'float64') {
    // Resolve string dtype to enum
    let resolvedDtype: DType;
    if (typeof dtype === 'string') {
      switch (dtype) {
        case 'float16':
          resolvedDtype = DType.Float16;
          break;
        case 'float32':
          resolvedDtype = DType.Float32;
          break;
        case 'float64':
          resolvedDtype = DType.Float64;
          break;
        default:
          throw new Error(`Unknown float dtype: ${dtype}`);
      }
    } else {
      resolvedDtype = dtype;
    }

    this.dtype = resolvedDtype;

    switch (resolvedDtype) {
      case DType.Float16:
        this.bits = 16;
        this.eps = 0.0009765625; // 2^-10
        this.maxexp = 16;
        this.minexp = -14;
        this.max = 65504.0;
        this.min = 6.103515625e-5; // Smallest positive normal
        this.tiny = 6.103515625e-5;
        this.precision = 3;
        this.iexp = 5;
        break;

      case DType.Float32:
        this.bits = 32;
        this.eps = 1.1920929e-7; // 2^-23
        this.maxexp = 128;
        this.minexp = -126;
        this.max = 3.4028235e38;
        this.min = 1.1754944e-38;
        this.tiny = 1.1754944e-38;
        this.precision = 6;
        this.iexp = 8;
        break;

      case DType.Float64:
        this.bits = 64;
        this.eps = 2.220446049250313e-16; // 2^-52
        this.maxexp = 1024;
        this.minexp = -1022;
        this.max = 1.7976931348623157e308;
        this.min = 2.2250738585072014e-308;
        this.tiny = 2.2250738585072014e-308;
        this.precision = 15;
        this.iexp = 11;
        break;

      default:
        throw new Error(`finfo not available for dtype: ${resolvedDtype}. Use float16, float32, or float64.`);
    }
  }

  /**
   * String representation
   */
  toString(): string {
    return (
      `finfo(dtype=${this.dtype}, ` +
      `bits=${this.bits}, ` +
      `eps=${this.eps}, ` +
      `max=${this.max}, ` +
      `min=${this.min})`
    );
  }
}

// =============================================================================
// Integer Type Information
// =============================================================================

/**
 * Machine limits for integer types.
 *
 * Provides properties describing the limits of integer data types.
 *
 * @example
 * import { iinfo, DType } from 'numjs';
 *
 * const i32 = iinfo(DType.Int32);
 * console.log(i32.min);  // -2147483648
 * console.log(i32.max);  // 2147483647
 * console.log(i32.bits); // 32
 *
 * @see https://numpy.org/doc/stable/reference/generated/numpy.iinfo.html
 */
export class iinfo {
  /** Number of bits in the type */
  readonly bits: number;

  /** Minimum value of the type */
  readonly min: number;

  /** Maximum value of the type */
  readonly max: number;

  /** Integer type being queried */
  readonly dtype: DType;

  constructor(dtype: DType | 'int8' | 'int16' | 'int32' | 'int64' | 'uint8' | 'uint16' | 'uint32' | 'uint64') {
    // Resolve string dtype to enum
    let resolvedDtype: DType;
    if (typeof dtype === 'string') {
      switch (dtype) {
        case 'int8':
          resolvedDtype = DType.Int8;
          break;
        case 'int16':
          resolvedDtype = DType.Int16;
          break;
        case 'int32':
          resolvedDtype = DType.Int32;
          break;
        case 'int64':
          resolvedDtype = DType.Int64;
          break;
        case 'uint8':
          resolvedDtype = DType.UInt8;
          break;
        case 'uint16':
          resolvedDtype = DType.UInt16;
          break;
        case 'uint32':
          resolvedDtype = DType.UInt32;
          break;
        case 'uint64':
          resolvedDtype = DType.UInt64;
          break;
        default:
          throw new Error(`Unknown integer dtype: ${dtype}`);
      }
    } else {
      resolvedDtype = dtype;
    }

    this.dtype = resolvedDtype;

    switch (resolvedDtype) {
      case DType.Bool:
        this.bits = 8;
        this.min = 0;
        this.max = 1;
        break;

      case DType.Int8:
        this.bits = 8;
        this.min = -128;
        this.max = 127;
        break;

      case DType.UInt8:
        this.bits = 8;
        this.min = 0;
        this.max = 255;
        break;

      case DType.Int16:
        this.bits = 16;
        this.min = -32768;
        this.max = 32767;
        break;

      case DType.UInt16:
        this.bits = 16;
        this.min = 0;
        this.max = 65535;
        break;

      case DType.Int32:
        this.bits = 32;
        this.min = -2147483648;
        this.max = 2147483647;
        break;

      case DType.UInt32:
        this.bits = 32;
        this.min = 0;
        this.max = 4294967295;
        break;

      case DType.Int64:
        // Note: JavaScript can't represent full int64 range precisely
        // Using BigInt internally may be needed for full precision
        this.bits = 64;
        this.min = -9223372036854775808;
        this.max = 9223372036854775807;
        break;

      case DType.UInt64:
        this.bits = 64;
        this.min = 0;
        this.max = 18446744073709551615;
        break;

      default:
        throw new Error(`iinfo not available for dtype: ${resolvedDtype}. Use an integer dtype.`);
    }
  }

  /**
   * String representation
   */
  toString(): string {
    return `iinfo(dtype=${this.dtype}, bits=${this.bits}, min=${this.min}, max=${this.max})`;
  }
}
```

---

## NumPy Reference

**Key NumPy source files for constants:**

| Constant | NumPy Source Location |
|----------|----------------------|
| `e`, `pi`, `euler_gamma` | `numpy/_core/umath.py` (from `_multiarray_umath` C extension) |
| `inf`, `nan` | `numpy/_core/numeric.py` (aliases PINF, NAN from umath) |
| `PINF`, `NINF`, `NAN` | `numpy/_core/umath.py` |
| `PZERO`, `NZERO` | `numpy/_core/umath.py` |
| `newaxis` | `numpy/_core/numeric.py` (defined as `None`) |
| `finfo` | `numpy/_core/getlimits.py` |
| `iinfo` | `numpy/_core/getlimits.py` |

**NumPy constant values (for validation):**

```python
>>> import numpy as np
>>> np.e
2.718281828459045
>>> np.pi
3.141592653589793
>>> np.euler_gamma
0.5772156649015329
>>> np.inf
inf
>>> np.PINF
inf
>>> np.NINF
-inf
>>> np.nan
nan
>>> np.NAN
nan
>>> np.PZERO
0.0
>>> np.NZERO
-0.0
>>> np.newaxis is None
True
>>> 1 / np.PZERO
inf
>>> 1 / np.NZERO
-inf
```

---

## File Changes Summary

### New Files to Create

```
src/ts/
├── constants.ts       # Mathematical constants, special values
└── typeinfo.ts        # finfo and iinfo classes
```

### Files to Modify

```
src/ts/index.ts
├── Import and export all constants from constants.ts
├── Import and export finfo, iinfo from typeinfo.ts
└── Ensure newaxis is exported (already exported from slice.ts)
```

---

## Detailed File Specifications

### `src/ts/constants.ts`

```typescript
/**
 * NumJS Constants
 *
 * Mathematical and special constants matching NumPy's constant exports.
 * Reference: numpy/_core/numeric.py, numpy/_core/umath.py
 */

// =============================================================================
// Mathematical Constants
// =============================================================================

/**
 * Euler's number, the base of natural logarithms.
 * e ≈ 2.718281828459045
 */
export const e: number = 2.718281828459045;

/**
 * Pi, the ratio of a circle's circumference to its diameter.
 * π ≈ 3.141592653589793
 */
export const pi: number = 3.141592653589793;

/**
 * Euler-Mascheroni constant.
 * γ ≈ 0.5772156649015329
 */
export const euler_gamma: number = 0.5772156649015329;

// =============================================================================
// Special Floating-Point Values
// =============================================================================

/** IEEE 754 positive infinity */
export const inf: number = Infinity;

/** IEEE 754 positive infinity (alias) */
export const PINF: number = Infinity;

/** IEEE 754 negative infinity */
export const NINF: number = -Infinity;

/** IEEE 754 quiet NaN */
export const nan: number = NaN;

/** IEEE 754 quiet NaN (alias) */
export const NAN: number = NaN;

/** IEEE 754 positive zero */
export const PZERO: number = 0.0;

/** IEEE 754 negative zero */
export const NZERO: number = -0.0;

// =============================================================================
// Indexing Constants (re-exported)
// =============================================================================

// newaxis is exported from slice.ts
export { newaxis } from './slice.js';
```

### `src/ts/index.ts` (additions)

```typescript
// Constants (Phase 12)
export {
  // Mathematical constants
  e,
  pi,
  euler_gamma,
  // Special floating-point values
  inf,
  PINF,
  NINF,
  nan,
  NAN,
  PZERO,
  NZERO,
} from './constants.js';

// Type information classes
export { finfo, iinfo } from './typeinfo.js';
```

---

## Dependencies

Phase 12 has minimal dependencies:

- **12.1-12.2**: No dependencies (pure JavaScript constants)
- **12.3**: Depends on `slice.ts` (already implemented in Level 2)
- **12.4**: Depends on `types.ts` for DType enum (already implemented)

This phase can be implemented independently of other incomplete phases.

---

## Implementation Order

```
Day 1: Core Constants
├── Create src/ts/constants.ts
├── Add mathematical constants (e, pi, euler_gamma)
├── Add special floating-point values (inf, nan, etc.)
├── Re-export newaxis from slice.ts
└── Add unit tests for constant values

Day 2: Type Information Classes
├── Create src/ts/typeinfo.ts
├── Implement finfo class for float16, float32, float64
├── Implement iinfo class for all integer types
└── Add unit tests for type limits

Day 3: Integration & Polish
├── Update src/ts/index.ts with exports
├── Add comprehensive test suite
├── Verify values match NumPy exactly
└── Add documentation/examples
```

---

## Verification Plan

### Unit Tests

```typescript
// tests/constants.test.ts
import { describe, it, expect } from 'vitest';
import {
  e, pi, euler_gamma,
  inf, PINF, NINF, nan, NAN,
  PZERO, NZERO, newaxis,
  finfo, iinfo, DType
} from '../src';

describe('Mathematical Constants', () => {
  it('e equals Euler\'s number', () => {
    expect(e).toBeCloseTo(2.718281828459045, 15);
    expect(e).toBe(Math.E);
  });

  it('pi equals π', () => {
    expect(pi).toBeCloseTo(3.141592653589793, 15);
    expect(pi).toBe(Math.PI);
  });

  it('euler_gamma equals Euler-Mascheroni constant', () => {
    expect(euler_gamma).toBeCloseTo(0.5772156649015329, 15);
  });
});

describe('Special Floating-Point Values', () => {
  it('inf is positive infinity', () => {
    expect(inf).toBe(Infinity);
    expect(PINF).toBe(Infinity);
    expect(inf).toBe(PINF);
    expect(Number.isFinite(inf)).toBe(false);
    expect(inf > Number.MAX_VALUE).toBe(true);
  });

  it('NINF is negative infinity', () => {
    expect(NINF).toBe(-Infinity);
    expect(Number.isFinite(NINF)).toBe(false);
    expect(NINF < -Number.MAX_VALUE).toBe(true);
  });

  it('nan is NaN', () => {
    expect(Number.isNaN(nan)).toBe(true);
    expect(Number.isNaN(NAN)).toBe(true);
    // NaN is not equal to itself
    expect(nan === nan).toBe(false);
    expect(NAN === NAN).toBe(false);
  });

  it('PZERO and NZERO are signed zeros', () => {
    expect(PZERO).toBe(0);
    expect(NZERO).toBe(-0);
    // They compare equal
    expect(PZERO === NZERO).toBe(true);
    // But behave differently with division
    expect(1 / PZERO).toBe(Infinity);
    expect(1 / NZERO).toBe(-Infinity);
    // Object.is can distinguish them
    expect(Object.is(PZERO, 0)).toBe(true);
    expect(Object.is(NZERO, -0)).toBe(true);
    expect(Object.is(PZERO, NZERO)).toBe(false);
  });
});

describe('newaxis', () => {
  it('newaxis is a unique symbol', () => {
    expect(typeof newaxis).toBe('symbol');
  });
});

describe('finfo', () => {
  it('provides Float32 limits', () => {
    const f32 = new finfo(DType.Float32);
    expect(f32.bits).toBe(32);
    expect(f32.eps).toBeCloseTo(1.1920929e-7, 12);
    expect(f32.max).toBeCloseTo(3.4028235e38, 30);
    expect(f32.min).toBeCloseTo(1.1754944e-38, 45);
  });

  it('provides Float64 limits', () => {
    const f64 = new finfo(DType.Float64);
    expect(f64.bits).toBe(64);
    expect(f64.eps).toBeCloseTo(2.220446049250313e-16, 25);
    expect(f64.max).toBe(Number.MAX_VALUE);
    expect(f64.min).toBe(Number.MIN_VALUE);
  });

  it('accepts string dtype', () => {
    const f32 = new finfo('float32');
    expect(f32.bits).toBe(32);
  });
});

describe('iinfo', () => {
  it('provides Int8 limits', () => {
    const i8 = new iinfo(DType.Int8);
    expect(i8.bits).toBe(8);
    expect(i8.min).toBe(-128);
    expect(i8.max).toBe(127);
  });

  it('provides UInt8 limits', () => {
    const u8 = new iinfo(DType.UInt8);
    expect(u8.bits).toBe(8);
    expect(u8.min).toBe(0);
    expect(u8.max).toBe(255);
  });

  it('provides Int32 limits', () => {
    const i32 = new iinfo(DType.Int32);
    expect(i32.bits).toBe(32);
    expect(i32.min).toBe(-2147483648);
    expect(i32.max).toBe(2147483647);
  });

  it('accepts string dtype', () => {
    const i32 = new iinfo('int32');
    expect(i32.bits).toBe(32);
  });
});
```

### NumPy Comparison Validation

```python
# tests/python/generate_constants_test.py
import numpy as np
import json

constants = {
    "mathematical": {
        "e": float(np.e),
        "pi": float(np.pi),
        "euler_gamma": float(np.euler_gamma),
    },
    "special_values": {
        "inf": float('inf'),  # np.inf
        "PINF": float('inf'),
        "NINF": float('-inf'),
        # NaN cannot be compared with JSON, verify with isnan
    },
    "finfo_float32": {
        "bits": 32,
        "eps": float(np.finfo(np.float32).eps),
        "max": float(np.finfo(np.float32).max),
        "min": float(np.finfo(np.float32).tiny),
    },
    "finfo_float64": {
        "bits": 64,
        "eps": float(np.finfo(np.float64).eps),
        "max": float(np.finfo(np.float64).max),
        "min": float(np.finfo(np.float64).tiny),
    },
    "iinfo_int32": {
        "bits": 32,
        "min": int(np.iinfo(np.int32).min),
        "max": int(np.iinfo(np.int32).max),
    },
    "iinfo_uint8": {
        "bits": 8,
        "min": int(np.iinfo(np.uint8).min),
        "max": int(np.iinfo(np.uint8).max),
    },
}

with open("tests/fixtures/constants_vectors.json", "w") as f:
    json.dump(constants, f, indent=2)
```

---

## API Surface After Phase 12

```typescript
// Mathematical Constants
export const e: number;              // 2.718281828459045
export const pi: number;             // 3.141592653589793
export const euler_gamma: number;    // 0.5772156649015329

// Special Floating-Point Values
export const inf: number;            // Infinity
export const PINF: number;           // Infinity (alias)
export const NINF: number;           // -Infinity
export const nan: number;            // NaN
export const NAN: number;            // NaN (alias)
export const PZERO: number;          // 0.0
export const NZERO: number;          // -0.0

// Indexing (re-exported from slice.ts)
export const newaxis: unique symbol; // For expanding dimensions

// Type Information Classes
export class finfo {
  constructor(dtype: DType | string);
  readonly bits: number;
  readonly eps: number;
  readonly max: number;
  readonly min: number;
  readonly tiny: number;
  readonly precision: number;
  readonly dtype: DType;
}

export class iinfo {
  constructor(dtype: DType | string);
  readonly bits: number;
  readonly min: number;
  readonly max: number;
  readonly dtype: DType;
}
```

---

## Notes

### JavaScript/IEEE 754 Alignment

The constants in this phase leverage JavaScript's built-in IEEE 754 floating-point representation:

- `Math.E` and `Math.PI` provide the same values as NumPy's `e` and `pi`
- `Infinity` and `-Infinity` are IEEE 754 infinities
- `NaN` is the IEEE 754 quiet NaN
- Signed zeros (`+0` and `-0`) are supported natively

### Precision Considerations

- Float64 values match JavaScript's native `number` type exactly
- Float32 and Float16 values are stored with reduced precision when used in typed arrays
- The `finfo` class reports the limits for stored values, not JavaScript's native precision

### newaxis Behavior

Unlike NumPy where `newaxis` is literally `None`, in NumJS it's a unique symbol. This provides:
- Type safety (can't accidentally use `null` or `undefined`)
- Clear intent in code
- Runtime distinction from other values

---

## Critical Dependencies for Later Phases

Phase 12 constants are used by:

- **Phase 13 (numpy.linalg)**: Uses `pi` for angle calculations
- **Phase 14 (numpy.fft)**: Uses `pi` for frequency calculations
- **Phase 15 (numpy.random)**: Uses `e` for exponential distribution
- **All phases**: Use `inf`, `nan` for edge case handling and validation

Phase 12 should be completed before Phases 13-15.
