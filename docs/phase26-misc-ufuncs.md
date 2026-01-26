# Phase 26: Miscellaneous Ufuncs Implementation Plan

Complete implementation roadmap for remaining universal functions including floating-point utilities, rational functions, and special functions.

---

## ⚠️ Implementation Guideline

**The original NumPy source code is available in `/numpy`.**

Key reference files:
- `numpy/_core/umath.py` - Ufunc definitions
- `numpy/_core/src/umath/loops.c.src` - Ufunc implementations
- `numpy/_core/src/npymath/npy_math.c` - Math functions

Implementation should follow NumPy's algorithms and error handling for consistency.

---

## Current State (Pre-Phase 26)

```
Already Implemented:
├── Arithmetic (add, subtract, multiply, divide, power, mod, etc.)
├── Trigonometric (sin, cos, tan, arcsin, etc.)
├── Hyperbolic (sinh, cosh, tanh, etc.)
├── Exponential/Logarithmic (exp, log, log2, log10, etc.)
├── Comparison (greater, less, equal, etc.)
├── Logical (and, or, xor, not)
├── Bitwise (and, or, xor, shift)
├── Rounding (floor, ceil, round, trunc)
└── Predicates (isnan, isinf, isfinite)

Missing:
├── Floating Point
│   ├── frexp(x) → (mantissa, exponent)
│   ├── ldexp(x, i) → x * 2^i
│   ├── nextafter(x1, x2)
│   ├── spacing(x)
│   └── modf(x) → (fractional, integer)
├── Rational
│   ├── gcd(x1, x2)
│   └── lcm(x1, x2)
├── Special
│   ├── sinc(x)
│   ├── heaviside(x1, x2)
│   └── divmod(x1, x2)
└── Bitwise
    └── bitwise_count(x)
```

---

## Phase 26 Dependency Tree

```
PHASE 26: MISCELLANEOUS UFUNCS
│
├── 26.1 Floating Point Utilities (TypeScript + WASM)
│   ├── 26.1.1 frexp(x) → (mantissa, exponent)
│   ├── 26.1.2 ldexp(x1, x2) → x1 * 2^x2
│   ├── 26.1.3 nextafter(x1, x2)
│   ├── 26.1.4 spacing(x)
│   └── 26.1.5 modf(x) → (fractional, integer)
│
│   Dependencies: NDArray core, dtype system
│
├── 26.2 Rational Functions (TypeScript + WASM)
│   ├── 26.2.1 gcd(x1, x2)
│   └── 26.2.2 lcm(x1, x2)
│
│   Dependencies: Integer ufunc infrastructure
│
├── 26.3 Special Functions (TypeScript)
│   ├── 26.3.1 sinc(x)
│   ├── 26.3.2 heaviside(x1, x2)
│   └── 26.3.3 divmod(x1, x2)
│
│   Dependencies: sin, floor_divide, mod
│
└── 26.4 Bitwise Functions (TypeScript + WASM)
    └── 26.4.1 bitwise_count(x)

    Dependencies: Integer operations
```

---

## Detailed Implementation Specifications

### 26.1 Floating Point Utilities

#### 26.1.1 frexp

**File:** `src/ts/ufunc.ts` (additions)

```typescript
/**
 * Decompose the elements of x into mantissa and twos exponent.
 *
 * Returns (mantissa, exponent), where x = mantissa * 2^exponent.
 * The mantissa lies in the open interval(-1, 1), while the twos
 * exponent is a signed integer.
 *
 * @param x - Input array
 * @returns Tuple of (mantissa, exponent) arrays
 *
 * @example
 * const [m, e] = frexp([1.0, 2.0, 4.0, 8.0]);
 * // m: [0.5, 0.5, 0.5, 0.5]
 * // e: [1, 2, 3, 4]
 *
 * // Verify: x = m * 2^e
 * ldexp(m, e)  // [1.0, 2.0, 4.0, 8.0]
 */
export function frexp(x: NDArray | ArrayLike<number>): [NDArray, NDArray] {
  const arr = asarray(x);
  const dtype = _floatType(arr.dtype);

  const mantissa = empty(arr.shape, dtype);
  const exponent = empty(arr.shape, DType.Int32);

  if (arr.size > 1000) {
    // Use WASM for large arrays
    return _frexpWASM(arr, mantissa, exponent);
  }

  // TypeScript implementation
  const data = arr.toArray();
  const mData: number[] = [];
  const eData: number[] = [];

  for (let i = 0; i < data.length; i++) {
    const [m, e] = _frexpScalar(data[i]);
    mData.push(m);
    eData.push(e);
  }

  return [
    fromArray(mData, arr.shape, dtype),
    fromArray(eData, arr.shape, DType.Int32)
  ];
}

/**
 * Decompose a single float into mantissa and exponent.
 */
function _frexpScalar(x: number): [number, number] {
  if (x === 0 || !isFinite(x)) {
    return [x, 0];
  }

  const absX = Math.abs(x);
  let exponent = Math.floor(Math.log2(absX)) + 1;
  let mantissa = x / Math.pow(2, exponent);

  // Ensure mantissa is in [0.5, 1)
  while (Math.abs(mantissa) >= 1) {
    mantissa /= 2;
    exponent++;
  }
  while (Math.abs(mantissa) < 0.5 && mantissa !== 0) {
    mantissa *= 2;
    exponent--;
  }

  return [mantissa, exponent];
}
```

#### 26.1.2 ldexp

**File:** `src/ts/ufunc.ts` (additions)

```typescript
/**
 * Returns x1 * 2^x2, element-wise.
 *
 * The mantissas x1 and twos exponents x2 are used to construct
 * floating point numbers x = x1 * 2^x2.
 *
 * @param x1 - Array of multipliers (mantissas)
 * @param x2 - Array of twos exponents (integers)
 * @returns x1 * 2^x2
 *
 * @example
 * ldexp([0.5, 0.5, 0.5], [1, 2, 3])  // [1.0, 2.0, 4.0]
 *
 * // Inverse of frexp
 * const [m, e] = frexp([1.0, 2.0, 4.0]);
 * ldexp(m, e)  // [1.0, 2.0, 4.0]
 */
export function ldexp(
  x1: NDArray | ArrayLike<number>,
  x2: NDArray | ArrayLike<number>
): NDArray {
  const arr1 = asarray(x1);
  const arr2 = asarray(x2);

  // Broadcast
  const [a1, a2] = broadcastArrays(arr1, arr2);
  const dtype = _floatType(a1.dtype);

  if (a1.size > 1000) {
    return _ldexpWASM(a1, a2);
  }

  // TypeScript implementation
  const data1 = a1.toArray();
  const data2 = a2.toArray();
  const result: number[] = [];

  for (let i = 0; i < data1.length; i++) {
    result.push(data1[i] * Math.pow(2, data2[i]));
  }

  return fromArray(result, a1.shape, dtype);
}
```

#### 26.1.3 nextafter

**File:** `src/ts/ufunc.ts` (additions)

```typescript
/**
 * Return the next floating-point value after x1 towards x2, element-wise.
 *
 * @param x1 - Starting values
 * @param x2 - Direction values
 * @returns Next representable floating-point values
 *
 * @example
 * nextafter(1.0, 2.0)  // Slightly larger than 1.0
 * nextafter(1.0, 0.0)  // Slightly smaller than 1.0
 * nextafter(1.0, 1.0)  // 1.0 (no change)
 */
export function nextafter(
  x1: NDArray | ArrayLike<number>,
  x2: NDArray | ArrayLike<number>
): NDArray {
  const arr1 = asarray(x1);
  const arr2 = asarray(x2);

  const [a1, a2] = broadcastArrays(arr1, arr2);
  const dtype = _floatType(a1.dtype);

  const data1 = a1.toArray();
  const data2 = a2.toArray();
  const result: number[] = [];

  for (let i = 0; i < data1.length; i++) {
    result.push(_nextafterScalar(data1[i], data2[i]));
  }

  return fromArray(result, a1.shape, dtype);
}

/**
 * Compute nextafter for scalar values.
 */
function _nextafterScalar(x1: number, x2: number): number {
  if (isNaN(x1) || isNaN(x2)) {
    return NaN;
  }

  if (x1 === x2) {
    return x1;
  }

  if (x1 === 0) {
    // Smallest subnormal in direction of x2
    return x2 > 0 ? Number.MIN_VALUE : -Number.MIN_VALUE;
  }

  // Use typed array to manipulate bits
  const buffer = new ArrayBuffer(8);
  const f64 = new Float64Array(buffer);
  const u64 = new BigUint64Array(buffer);

  f64[0] = x1;
  const bits = u64[0];

  // Increment or decrement based on direction
  if ((x1 > 0 && x2 > x1) || (x1 < 0 && x2 > x1)) {
    u64[0] = bits + 1n;
  } else {
    u64[0] = bits - 1n;
  }

  return f64[0];
}
```

#### 26.1.4 spacing

**File:** `src/ts/ufunc.ts` (additions)

```typescript
/**
 * Return the distance between x and the nearest adjacent number.
 *
 * This is effectively the ulp (unit in the last place) of x.
 *
 * @param x - Input values
 * @returns Distance to next representable value
 *
 * @example
 * spacing(1.0)  // ~2.22e-16 (machine epsilon at 1.0)
 * spacing(1e10)  // Larger, because floats are sparser at larger values
 */
export function spacing(x: NDArray | ArrayLike<number>): NDArray {
  const arr = asarray(x);
  const dtype = _floatType(arr.dtype);

  const data = arr.toArray();
  const result: number[] = [];

  for (let i = 0; i < data.length; i++) {
    result.push(_spacingScalar(data[i]));
  }

  return fromArray(result, arr.shape, dtype);
}

/**
 * Compute spacing for scalar value.
 */
function _spacingScalar(x: number): number {
  if (!isFinite(x)) {
    return NaN;
  }

  if (x === 0) {
    return Number.MIN_VALUE;
  }

  const next = _nextafterScalar(x, Infinity);
  return Math.abs(next - x);
}
```

#### 26.1.5 modf

**File:** `src/ts/ufunc.ts` (additions)

```typescript
/**
 * Return the fractional and integral parts of an array, element-wise.
 *
 * The fractional and integral parts are negative if the given number is negative.
 *
 * @param x - Input array
 * @returns Tuple of (fractional part, integral part)
 *
 * @example
 * const [frac, intg] = modf([3.5, -2.7]);
 * // frac: [0.5, -0.7]
 * // intg: [3.0, -2.0]
 */
export function modf(x: NDArray | ArrayLike<number>): [NDArray, NDArray] {
  const arr = asarray(x);
  const dtype = _floatType(arr.dtype);

  const data = arr.toArray();
  const fracData: number[] = [];
  const intgData: number[] = [];

  for (let i = 0; i < data.length; i++) {
    const val = data[i];

    if (!isFinite(val)) {
      if (isNaN(val)) {
        fracData.push(NaN);
        intgData.push(NaN);
      } else {
        // Infinity
        fracData.push(val > 0 ? 0 : -0);
        intgData.push(val);
      }
    } else {
      const intPart = Math.trunc(val);
      const fracPart = val - intPart;

      intgData.push(intPart);
      fracData.push(fracPart);
    }
  }

  return [
    fromArray(fracData, arr.shape, dtype),
    fromArray(intgData, arr.shape, dtype)
  ];
}
```

---

### 26.2 Rational Functions

#### 26.2.1 gcd

**File:** `src/ts/ufunc.ts` (additions)

```typescript
/**
 * Returns the greatest common divisor of |x1| and |x2|.
 *
 * @param x1 - First array of integers
 * @param x2 - Second array of integers
 * @returns Greatest common divisor
 *
 * @example
 * gcd(12, 8)  // 4
 * gcd([12, 15, 20], [8, 10, 15])  // [4, 5, 5]
 * gcd(0, 5)  // 5
 * gcd(0, 0)  // 0
 */
export function gcd(
  x1: NDArray | ArrayLike<number>,
  x2: NDArray | ArrayLike<number>
): NDArray {
  const arr1 = asarray(x1);
  const arr2 = asarray(x2);

  // Validate integer types
  if (!_isIntegerDtype(arr1.dtype) || !_isIntegerDtype(arr2.dtype)) {
    throw new TypeError('gcd requires integer inputs');
  }

  const [a1, a2] = broadcastArrays(arr1, arr2);
  const dtype = _promoteIntegerType(a1.dtype, a2.dtype);

  if (a1.size > 1000) {
    return _gcdWASM(a1, a2, dtype);
  }

  const data1 = a1.toArray();
  const data2 = a2.toArray();
  const result: number[] = [];

  for (let i = 0; i < data1.length; i++) {
    result.push(_gcdScalar(Math.abs(data1[i]), Math.abs(data2[i])));
  }

  return fromArray(result, a1.shape, dtype);
}

/**
 * Euclidean algorithm for GCD.
 */
function _gcdScalar(a: number, b: number): number {
  a = Math.abs(Math.floor(a));
  b = Math.abs(Math.floor(b));

  while (b !== 0) {
    const t = b;
    b = a % b;
    a = t;
  }

  return a;
}
```

#### 26.2.2 lcm

**File:** `src/ts/ufunc.ts` (additions)

```typescript
/**
 * Returns the lowest common multiple of |x1| and |x2|.
 *
 * @param x1 - First array of integers
 * @param x2 - Second array of integers
 * @returns Lowest common multiple
 *
 * @example
 * lcm(12, 8)  // 24
 * lcm([4, 6], [8, 9])  // [8, 18]
 * lcm(0, 5)  // 0
 */
export function lcm(
  x1: NDArray | ArrayLike<number>,
  x2: NDArray | ArrayLike<number>
): NDArray {
  const arr1 = asarray(x1);
  const arr2 = asarray(x2);

  if (!_isIntegerDtype(arr1.dtype) || !_isIntegerDtype(arr2.dtype)) {
    throw new TypeError('lcm requires integer inputs');
  }

  const [a1, a2] = broadcastArrays(arr1, arr2);
  const dtype = _promoteIntegerType(a1.dtype, a2.dtype);

  const data1 = a1.toArray();
  const data2 = a2.toArray();
  const result: number[] = [];

  for (let i = 0; i < data1.length; i++) {
    result.push(_lcmScalar(Math.abs(data1[i]), Math.abs(data2[i])));
  }

  return fromArray(result, a1.shape, dtype);
}

/**
 * Compute LCM using GCD.
 */
function _lcmScalar(a: number, b: number): number {
  if (a === 0 || b === 0) {
    return 0;
  }

  return Math.abs(a * b) / _gcdScalar(a, b);
}

/**
 * Check if dtype is an integer type.
 */
function _isIntegerDtype(dtype: DType): boolean {
  return dtype === DType.Int8 || dtype === DType.Int16 ||
         dtype === DType.Int32 || dtype === DType.Int64 ||
         dtype === DType.UInt8 || dtype === DType.UInt16 ||
         dtype === DType.UInt32 || dtype === DType.UInt64 ||
         dtype === DType.Bool;
}

/**
 * Promote integer types for output.
 */
function _promoteIntegerType(dtype1: DType, dtype2: DType): DType {
  // Return the larger type
  const order = [
    DType.Bool, DType.Int8, DType.UInt8,
    DType.Int16, DType.UInt16,
    DType.Int32, DType.UInt32,
    DType.Int64, DType.UInt64
  ];

  const idx1 = order.indexOf(dtype1);
  const idx2 = order.indexOf(dtype2);

  return order[Math.max(idx1, idx2)];
}
```

---

### 26.3 Special Functions

#### 26.3.1 sinc

**File:** `src/ts/ufunc.ts` (additions)

```typescript
/**
 * Return the normalized sinc function.
 *
 * The sinc function is sin(pi*x) / (pi*x).
 *
 * @param x - Input array
 * @returns sinc(x)
 *
 * @example
 * sinc(0)  // 1.0 (limit as x -> 0)
 * sinc(1)  // 0.0
 * sinc(0.5)  // ~0.637
 *
 * // Useful for signal processing (ideal lowpass filter)
 * sinc(linspace(-3, 3, 7))
 */
export function sinc(x: NDArray | ArrayLike<number>): NDArray {
  const arr = asarray(x);
  const dtype = _floatType(arr.dtype);

  const data = arr.toArray();
  const result: number[] = [];

  for (let i = 0; i < data.length; i++) {
    result.push(_sincScalar(data[i]));
  }

  return fromArray(result, arr.shape, dtype);
}

/**
 * Compute sinc for scalar value.
 */
function _sincScalar(x: number): number {
  if (x === 0) {
    return 1.0;
  }

  const pix = Math.PI * x;
  return Math.sin(pix) / pix;
}
```

#### 26.3.2 heaviside

**File:** `src/ts/ufunc.ts` (additions)

```typescript
/**
 * Compute the Heaviside step function.
 *
 * H(x) = 0 if x < 0
 * H(x) = h0 if x == 0
 * H(x) = 1 if x > 0
 *
 * @param x1 - Input values
 * @param x2 - Value of the function at x1 == 0
 * @returns Heaviside step function values
 *
 * @example
 * heaviside([-1, 0, 1], 0.5)  // [0, 0.5, 1]
 * heaviside([-1, 0, 1], 1.0)  // [0, 1.0, 1]
 */
export function heaviside(
  x1: NDArray | ArrayLike<number>,
  x2: NDArray | ArrayLike<number>
): NDArray {
  const arr1 = asarray(x1);
  const arr2 = asarray(x2);

  const [a1, a2] = broadcastArrays(arr1, arr2);
  const dtype = _floatType(promoteTypes(a1.dtype, a2.dtype));

  const data1 = a1.toArray();
  const data2 = a2.toArray();
  const result: number[] = [];

  for (let i = 0; i < data1.length; i++) {
    const x = data1[i];
    const h0 = data2[i];

    if (isNaN(x)) {
      result.push(NaN);
    } else if (x < 0) {
      result.push(0);
    } else if (x === 0) {
      result.push(h0);
    } else {
      result.push(1);
    }
  }

  return fromArray(result, a1.shape, dtype);
}
```

#### 26.3.3 divmod

**File:** `src/ts/ufunc.ts` (additions)

```typescript
/**
 * Return element-wise quotient and remainder simultaneously.
 *
 * Equivalent to (x1 // x2, x1 % x2), but computes both simultaneously
 * for efficiency.
 *
 * @param x1 - Dividend array
 * @param x2 - Divisor array
 * @returns Tuple of (quotient, remainder)
 *
 * @example
 * const [q, r] = divmod(10, 3);
 * // q: 3, r: 1
 *
 * const [q, r] = divmod([10, 11, 12], 3);
 * // q: [3, 3, 4], r: [1, 2, 0]
 *
 * // For floats, uses floor division
 * const [q, r] = divmod(7.5, 2.5);
 * // q: 3.0, r: 0.0
 */
export function divmod(
  x1: NDArray | ArrayLike<number>,
  x2: NDArray | ArrayLike<number>
): [NDArray, NDArray] {
  const arr1 = asarray(x1);
  const arr2 = asarray(x2);

  const [a1, a2] = broadcastArrays(arr1, arr2);
  const dtype = promoteTypes(a1.dtype, a2.dtype);

  const quotient = floor_divide(a1, a2);
  const remainder = mod(a1, a2);

  return [quotient, remainder];
}
```

---

### 26.4 Bitwise Functions

#### 26.4.1 bitwise_count

**File:** `src/ts/ufunc.ts` (additions)

```typescript
/**
 * Computes the number of 1-bits in the absolute value of x.
 *
 * Also known as popcount or population count.
 *
 * @param x - Input array of integers
 * @returns Number of 1-bits in each element
 *
 * @example
 * bitwise_count(7)  // 3 (binary: 111)
 * bitwise_count([0, 1, 2, 3, 4, 5, 6, 7])
 * // [0, 1, 1, 2, 1, 2, 2, 3]
 *
 * bitwise_count(-1)  // All bits set for the type
 */
export function bitwise_count(x: NDArray | ArrayLike<number>): NDArray {
  const arr = asarray(x);

  if (!_isIntegerDtype(arr.dtype)) {
    throw new TypeError('bitwise_count requires integer inputs');
  }

  const data = arr.toArray();
  const result: number[] = [];

  for (let i = 0; i < data.length; i++) {
    result.push(_popcountScalar(data[i], arr.dtype));
  }

  return fromArray(result, arr.shape, DType.UInt8);
}

/**
 * Count number of set bits in an integer.
 */
function _popcountScalar(x: number, dtype: DType): number {
  // Handle negative numbers based on dtype
  let value: number;

  if (x < 0) {
    // Get bit width for this dtype
    const bitWidth = _getBitWidth(dtype);
    // Convert to unsigned representation
    value = (1 << bitWidth) + x;
  } else {
    value = Math.floor(Math.abs(x));
  }

  // Brian Kernighan's algorithm
  let count = 0;
  while (value > 0) {
    value &= value - 1;
    count++;
  }

  return count;
}

/**
 * Get bit width for integer dtype.
 */
function _getBitWidth(dtype: DType): number {
  switch (dtype) {
    case DType.Int8:
    case DType.UInt8:
    case DType.Bool:
      return 8;
    case DType.Int16:
    case DType.UInt16:
      return 16;
    case DType.Int32:
    case DType.UInt32:
      return 32;
    case DType.Int64:
    case DType.UInt64:
      return 64;
    default:
      return 64;
  }
}
```

---

### 26.5 WASM Acceleration

**File:** `src/wasm/misc_ufunc.c` (new file)

```c
#ifndef NUMJS_MISC_UFUNC_H
#define NUMJS_MISC_UFUNC_H

#include "ndarray.h"
#include <math.h>

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#define EXPORT EMSCRIPTEN_KEEPALIVE
#else
#define EXPORT
#endif

/* ============ Floating Point Utilities ============ */

/**
 * Decompose float into mantissa and exponent.
 */
EXPORT void frexp_f64(int32_t n, const double* x,
                       double* mantissa, int32_t* exponent);

EXPORT void frexp_f32(int32_t n, const float* x,
                       float* mantissa, int32_t* exponent);

/**
 * Compute x1 * 2^x2.
 */
EXPORT void ldexp_f64(int32_t n, const double* x1, const int32_t* x2,
                       double* result);

/* ============ Rational Functions ============ */

/**
 * Greatest common divisor.
 */
EXPORT void gcd_i32(int32_t n, const int32_t* x1, const int32_t* x2,
                     int32_t* result);

EXPORT void gcd_i64(int32_t n, const int64_t* x1, const int64_t* x2,
                     int64_t* result);

/**
 * Lowest common multiple.
 */
EXPORT void lcm_i32(int32_t n, const int32_t* x1, const int32_t* x2,
                     int32_t* result);

EXPORT void lcm_i64(int32_t n, const int64_t* x1, const int64_t* x2,
                     int64_t* result);

/* ============ Bitwise ============ */

/**
 * Population count (number of 1-bits).
 */
EXPORT void bitwise_count_i32(int32_t n, const int32_t* x, uint8_t* result);
EXPORT void bitwise_count_i64(int32_t n, const int64_t* x, uint8_t* result);

#endif /* NUMJS_MISC_UFUNC_H */
```

**File:** `src/wasm/misc_ufunc.c` (implementation)

```c
#include "misc_ufunc.h"

/* ============ Floating Point Utilities ============ */

EXPORT void frexp_f64(int32_t n, const double* x,
                       double* mantissa, int32_t* exponent) {
    for (int32_t i = 0; i < n; i++) {
        int exp;
        mantissa[i] = frexp(x[i], &exp);
        exponent[i] = exp;
    }
}

EXPORT void frexp_f32(int32_t n, const float* x,
                       float* mantissa, int32_t* exponent) {
    for (int32_t i = 0; i < n; i++) {
        int exp;
        mantissa[i] = frexpf(x[i], &exp);
        exponent[i] = exp;
    }
}

EXPORT void ldexp_f64(int32_t n, const double* x1, const int32_t* x2,
                       double* result) {
    for (int32_t i = 0; i < n; i++) {
        result[i] = ldexp(x1[i], x2[i]);
    }
}

/* ============ Rational Functions ============ */

static int32_t gcd_scalar_i32(int32_t a, int32_t b) {
    a = a < 0 ? -a : a;
    b = b < 0 ? -b : b;

    while (b != 0) {
        int32_t t = b;
        b = a % b;
        a = t;
    }
    return a;
}

static int64_t gcd_scalar_i64(int64_t a, int64_t b) {
    a = a < 0 ? -a : a;
    b = b < 0 ? -b : b;

    while (b != 0) {
        int64_t t = b;
        b = a % b;
        a = t;
    }
    return a;
}

EXPORT void gcd_i32(int32_t n, const int32_t* x1, const int32_t* x2,
                     int32_t* result) {
    for (int32_t i = 0; i < n; i++) {
        result[i] = gcd_scalar_i32(x1[i], x2[i]);
    }
}

EXPORT void gcd_i64(int32_t n, const int64_t* x1, const int64_t* x2,
                     int64_t* result) {
    for (int32_t i = 0; i < n; i++) {
        result[i] = gcd_scalar_i64(x1[i], x2[i]);
    }
}

EXPORT void lcm_i32(int32_t n, const int32_t* x1, const int32_t* x2,
                     int32_t* result) {
    for (int32_t i = 0; i < n; i++) {
        int32_t a = x1[i] < 0 ? -x1[i] : x1[i];
        int32_t b = x2[i] < 0 ? -x2[i] : x2[i];

        if (a == 0 || b == 0) {
            result[i] = 0;
        } else {
            result[i] = (a / gcd_scalar_i32(a, b)) * b;
        }
    }
}

EXPORT void lcm_i64(int32_t n, const int64_t* x1, const int64_t* x2,
                     int64_t* result) {
    for (int32_t i = 0; i < n; i++) {
        int64_t a = x1[i] < 0 ? -x1[i] : x1[i];
        int64_t b = x2[i] < 0 ? -x2[i] : x2[i];

        if (a == 0 || b == 0) {
            result[i] = 0;
        } else {
            result[i] = (a / gcd_scalar_i64(a, b)) * b;
        }
    }
}

/* ============ Bitwise ============ */

/* Population count using Brian Kernighan's algorithm */
static uint8_t popcount_i32(int32_t x) {
    uint32_t v = (uint32_t)x;
    uint8_t count = 0;
    while (v) {
        v &= v - 1;
        count++;
    }
    return count;
}

static uint8_t popcount_i64(int64_t x) {
    uint64_t v = (uint64_t)x;
    uint8_t count = 0;
    while (v) {
        v &= v - 1;
        count++;
    }
    return count;
}

EXPORT void bitwise_count_i32(int32_t n, const int32_t* x, uint8_t* result) {
    for (int32_t i = 0; i < n; i++) {
        result[i] = popcount_i32(x[i]);
    }
}

EXPORT void bitwise_count_i64(int32_t n, const int64_t* x, uint8_t* result) {
    for (int32_t i = 0; i < n; i++) {
        result[i] = popcount_i64(x[i]);
    }
}
```

---

## File Changes Summary

### New Files to Create

```
src/wasm/
└── misc_ufunc.c         # WASM implementations

tests/ts/
└── misc-ufunc.test.ts   # Test suite
```

### Files to Modify

```
src/ts/ufunc.ts
├── Add frexp, ldexp
├── Add nextafter, spacing
├── Add modf
├── Add gcd, lcm
├── Add sinc, heaviside
├── Add divmod
└── Add bitwise_count

src/ts/index.ts
├── Export frexp, ldexp
├── Export nextafter, spacing
├── Export modf
├── Export gcd, lcm
├── Export sinc, heaviside
├── Export divmod
└── Export bitwise_count

src/ts/types.ts
└── Add WASM function declarations

scripts/build-wasm.sh
├── Add misc_ufunc.c to compilation
└── Add EXPORTED_FUNCTIONS
```

---

## Build Script Updates

Add to `scripts/build-wasm.sh` EXPORTED_FUNCTIONS:

```bash
"_frexp_f64",
"_frexp_f32",
"_ldexp_f64",
"_gcd_i32",
"_gcd_i64",
"_lcm_i32",
"_lcm_i64",
"_bitwise_count_i32",
"_bitwise_count_i64"
```

---

## Implementation Order

```
Phase 26.1: Floating Point (Day 1-2)
├── Day 1: frexp, ldexp
│   ├── TypeScript implementation
│   ├── WASM acceleration
│   └── Tests
│
└── Day 2: nextafter, spacing, modf
    ├── Bit manipulation for nextafter
    ├── Tests with edge cases
    └── Tests

Phase 26.2: Rational Functions (Day 3)
├── gcd implementation
├── lcm implementation
├── WASM acceleration
└── Tests

Phase 26.3: Special Functions (Day 4)
├── sinc
├── heaviside
├── divmod
└── Tests

Phase 26.4: Bitwise Functions (Day 5)
├── bitwise_count
├── WASM acceleration
└── Tests

Phase 26.5: Polish (Day 6)
├── Edge cases (NaN, Inf, zeros)
├── NumPy comparison tests
└── Documentation
```

---

## Verification Plan

After Phase 26 completion, verify:

```bash
# Build
npm run build

# Run tests
npm test

# Phase 26 specific tests:

# frexp/ldexp
✓ ldexp(frexp(x)[0], frexp(x)[1]) ≈ x
✓ frexp(0) = (0, 0)
✓ frexp(Inf) = (Inf, 0)

# nextafter/spacing
✓ nextafter(1, 2) > 1
✓ nextafter(1, 0) < 1
✓ spacing(1) ≈ 2.22e-16

# modf
✓ modf(3.5) = (0.5, 3)
✓ modf(-2.7) = (-0.7, -2)

# gcd/lcm
✓ gcd(12, 8) = 4
✓ lcm(12, 8) = 24
✓ gcd(0, 5) = 5
✓ lcm(0, 5) = 0

# sinc
✓ sinc(0) = 1
✓ sinc(1) ≈ 0
✓ sinc(0.5) ≈ 0.637

# heaviside
✓ heaviside(-1, 0.5) = 0
✓ heaviside(0, 0.5) = 0.5
✓ heaviside(1, 0.5) = 1

# divmod
✓ divmod(10, 3) = (3, 1)

# bitwise_count
✓ bitwise_count(7) = 3
✓ bitwise_count(0) = 0
```

Generate NumPy comparison vectors:

```python
import numpy as np
import json

tests = {
    "frexp": {
        "input": [1.0, 2.0, 4.0, 8.0, 0.5],
        "mantissa": [float(m) for m, _ in [np.frexp(x) for x in [1.0, 2.0, 4.0, 8.0, 0.5]]],
        "exponent": [int(e) for _, e in [np.frexp(x) for x in [1.0, 2.0, 4.0, 8.0, 0.5]]]
    },
    "ldexp": {
        "mantissa": [0.5, 0.5, 0.5],
        "exponent": [1, 2, 3],
        "expected": np.ldexp([0.5, 0.5, 0.5], [1, 2, 3]).tolist()
    },
    "gcd": {
        "x1": [12, 15, 20, 0],
        "x2": [8, 10, 15, 5],
        "expected": np.gcd([12, 15, 20, 0], [8, 10, 15, 5]).tolist()
    },
    "lcm": {
        "x1": [12, 4, 6],
        "x2": [8, 8, 9],
        "expected": np.lcm([12, 4, 6], [8, 8, 9]).tolist()
    },
    "sinc": {
        "input": [0, 0.5, 1, -1],
        "expected": np.sinc([0, 0.5, 1, -1]).tolist()
    },
    "heaviside": {
        "x": [-1, 0, 1],
        "h0": 0.5,
        "expected": np.heaviside([-1, 0, 1], 0.5).tolist()
    }
}

with open("tests/fixtures/misc_ufunc_vectors.json", "w") as f:
    json.dump(tests, f, indent=2)
```

---

## API Compatibility Notes

### NumPy Signature Match

```typescript
frexp(x) → [mantissa, exponent]
ldexp(x1, x2)
nextafter(x1, x2)
spacing(x)
modf(x) → [fractional, integral]
gcd(x1, x2)
lcm(x1, x2)
sinc(x)
heaviside(x1, x2)
divmod(x1, x2) → [quotient, remainder]
bitwise_count(x)
```

### Differences from NumPy

1. **Return types**: Functions returning tuples (frexp, modf, divmod) return arrays in NumJS, matching NumPy behavior.

2. **Integer precision**: GCD/LCM may have different overflow behavior for very large numbers due to JavaScript's number representation.

3. **nextafter precision**: Implementation uses Float64Array for bit manipulation, matching IEEE 754 double precision.
