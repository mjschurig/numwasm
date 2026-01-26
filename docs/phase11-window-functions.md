# Phase 11: Window Functions Implementation Plan

Phase 11 implements NumPy's window functions for signal processing applications. These functions generate window arrays used for smoothing, tapering, and spectral analysis. This is a TypeScript-only phase with optional WASM optimization for the Bessel function used by `kaiser()`.

---

## Overview

Window functions are essential for digital signal processing. They are used to:
- Reduce spectral leakage in FFT analysis
- Smooth transitions at signal boundaries (apodization/tapering)
- Filter design and convolution operations

All window functions follow the same pattern:
1. Accept window length `M` as input
2. Return empty array for `M < 1`
3. Return `[1.0]` for `M == 1`
4. Generate symmetric window centered at `(M-1)/2`

---

## NumPy Reference

**Source:** `numpy/lib/_function_base_impl.py`

| Function | Lines | Formula |
|----------|-------|---------|
| `blackman(M)` | 3037-3133 | `0.42 + 0.5*cos(π*n/(M-1)) + 0.08*cos(2π*n/(M-1))` |
| `bartlett(M)` | 3137-3240 | Triangular: `1 - |n|/(M-1)` |
| `hanning(M)` | 3244-3342 | `0.5 + 0.5*cos(π*n/(M-1))` |
| `hamming(M)` | 3346-3441 | `0.54 + 0.46*cos(π*n/(M-1))` |
| `kaiser(M, beta)` | 3596-3723 | `I₀(β*sqrt(1-(2n/(M-1))²)) / I₀(β)` |

**Supporting Functions:**
| Function | Lines | Description |
|----------|-------|-------------|
| `i0(x)` | 3533-3590 | Modified Bessel function of first kind, order 0 |
| `_chbevl(x, vals)` | 3508-3517 | Chebyshev polynomial evaluation |
| `_i0_1(x)` | 3520-3521 | Bessel I₀ for |x| ≤ 8 |
| `_i0_2(x)` | 3524-3525 | Bessel I₀ for |x| > 8 |

---

## Current State

```
src/ts/
├── NDArray.ts         # Core array class with arange, cos operations needed
├── indexing.ts        # where() function needed for bartlett
├── types.ts           # DType definitions
└── index.ts           # Public exports

Dependencies needed (must exist before Phase 11):
├── ✅ arange(start, stop, step)     # Generates n values
├── ❌ cos(arr) → unary ufunc        # Phase 4 - NOT STARTED
├── ❌ sqrt(arr) → unary ufunc       # Phase 4 - NOT STARTED
├── ❌ abs(arr) → unary ufunc        # Phase 4 - NOT STARTED
├── ❌ exp(arr) → unary ufunc        # Phase 4 - NOT STARTED
├── ✅ where(cond, x, y)             # Conditional selection
├── ❌ less_equal(x1, x2)            # Phase 4 - comparison ufunc
├── ❌ piecewise(x, condlist, funclist) # Phase 10 - NOT STARTED
└── ❌ multiply, add, divide         # Phase 4 - binary ufuncs
```

**Note:** Phase 11 has dependencies on Phase 4 (Ufuncs). Window functions need:
- Element-wise math operations: `cos`, `sqrt`, `abs`, `exp`
- Element-wise arithmetic: `multiply`, `add`, `divide`, `subtract`
- Comparison operations: `less_equal`

---

## Implementation Tree

```
PHASE 11: WINDOW FUNCTIONS
│
├── 11.1 Dependencies Check
│   ├── 11.1.1 Verify ufunc infrastructure exists (Phase 4)
│   ├── 11.1.2 Verify cos, sqrt, abs, exp ufuncs exist
│   ├── 11.1.3 Verify arithmetic ufuncs exist
│   └── 11.1.4 Verify comparison ufuncs exist
│
├── 11.2 Bessel Function I₀ (TypeScript)
│   ├── 11.2.1 _chbevl() - Chebyshev polynomial evaluation
│   ├── 11.2.2 _i0A, _i0B coefficient arrays
│   ├── 11.2.3 _i0_1(x) - Bessel I₀ for |x| ≤ 8
│   ├── 11.2.4 _i0_2(x) - Bessel I₀ for |x| > 8
│   └── 11.2.5 i0(x) - Main Bessel function
│
├── 11.3 Cosine-Based Windows (TypeScript)
│   ├── 11.3.1 blackman(M) - Near-optimal window
│   ├── 11.3.2 hanning(M) - Cosine bell / Hann window
│   └── 11.3.3 hamming(M) - Similar to Hanning with raised endpoints
│
├── 11.4 Other Windows (TypeScript)
│   ├── 11.4.1 bartlett(M) - Triangular window
│   └── 11.4.2 kaiser(M, beta) - Bessel-based window
│
└── 11.5 Optional WASM Optimization
    └── 11.5.1 i0_wasm() - WASM-accelerated Bessel function (for large arrays)
```

---

## Detailed Implementation Specifications

### 11.2 Bessel Function I₀ (TypeScript)

**File:** `src/ts/window.ts` (new file)

```typescript
/**
 * Window functions for signal processing.
 *
 * Implements NumPy-compatible window functions:
 * - blackman(M)
 * - bartlett(M)
 * - hanning(M)
 * - hamming(M)
 * - kaiser(M, beta)
 *
 * Reference: numpy/lib/_function_base_impl.py
 */

import { NDArray } from './NDArray.js';
import { DType } from './types.js';

// Chebyshev coefficients for I₀ computation (from Cephes library)
// Used for |x| <= 8
const _i0A: number[] = [
  -4.41534164647933937950e-18,
   3.33079451882223809783e-17,
  -2.43127984654795469359e-16,
   1.71539128555513303061e-15,
  -1.16853328779934516808e-14,
   7.67618549860493561688e-14,
  -4.85644678311192946090e-13,
   2.95505266312963983461e-12,
  -1.72682629144155570723e-11,
   9.67580903537323691224e-11,
  -5.18979560163526290666e-10,
   2.65982372468238665035e-9,
  -1.30002500998624804212e-8,
   6.04699502254191894932e-8,
  -2.67079385394061173391e-7,
   1.11738753912010371815e-6,
  -4.41673835845875056359e-6,
   1.64484480707288970893e-5,
  -5.75419501008210370398e-5,
   1.88502885095841655729e-4,
  -5.76375574538582365885e-4,
   1.63947561694133579842e-3,
  -4.32430999505057594430e-3,
   1.05464603945949983183e-2,
  -2.37374148058994688156e-2,
   4.93052842396707084878e-2,
  -9.49010970480476444210e-2,
   1.71620901522208775349e-1,
  -3.04682672343198398683e-1,
   6.76795274409476084995e-1,
];

// Chebyshev coefficients for I₀ computation (from Cephes library)
// Used for |x| > 8
const _i0B: number[] = [
  -7.23318048787475395456e-18,
  -4.83050448594418207126e-18,
   4.46562142029675999901e-17,
   3.46122286769746109310e-17,
  -2.82762398051658348494e-16,
  -3.42548561967721913462e-16,
   1.77256013305652638360e-15,
   3.81168066935262242075e-15,
  -9.55484669882830764870e-15,
  -4.15056934728722208663e-14,
   1.54008621752140982691e-14,
   3.85277838274214270114e-13,
   7.18012445138366623367e-13,
  -1.79417853150680611778e-12,
  -1.32158118404477131188e-11,
  -3.14991652796324136454e-11,
   1.18891471078464383424e-11,
   4.94060238822496958910e-10,
   3.39623202570838634515e-9,
   2.26666899049817806459e-8,
   2.04891858946906374183e-7,
   2.89137052083475648297e-6,
   6.88975834691682398426e-5,
   3.36911647825569408990e-3,
   8.04490411014108831608e-1,
];

/**
 * Evaluate Chebyshev polynomial using Clenshaw recurrence.
 * @internal
 */
function _chbevl(x: number, vals: number[]): number {
  let b0 = vals[0];
  let b1 = 0.0;
  let b2: number;

  for (let i = 1; i < vals.length; i++) {
    b2 = b1;
    b1 = b0;
    b0 = x * b1 - b2 + vals[i];
  }

  return 0.5 * (b0 - b2);
}

/**
 * Compute I₀(x) for |x| <= 8 using Chebyshev approximation.
 * @internal
 */
function _i0_1(x: number): number {
  return Math.exp(x) * _chbevl(x / 2.0 - 2, _i0A);
}

/**
 * Compute I₀(x) for |x| > 8 using asymptotic Chebyshev approximation.
 * @internal
 */
function _i0_2(x: number): number {
  return Math.exp(x) * _chbevl(32.0 / x - 2.0, _i0B) / Math.sqrt(x);
}

/**
 * Compute I₀(x) for a single scalar value.
 * @internal
 */
function _i0_scalar(x: number): number {
  const absX = Math.abs(x);
  return absX <= 8.0 ? _i0_1(absX) : _i0_2(absX);
}

/**
 * Modified Bessel function of the first kind, order 0.
 *
 * Usually denoted I₀. This implementation uses Chebyshev polynomial
 * approximations from the Cephes library, partitioning the domain
 * into [0, 8] and (8, ∞).
 *
 * @param x - Input array or scalar
 * @returns Modified Bessel function I₀ evaluated at x
 *
 * @example
 * ```typescript
 * i0(0)           // 1.0
 * i0([0, 1, 2])   // [1.0, 1.266..., 2.279...]
 * ```
 *
 * @remarks
 * Complex values are not supported and will throw an error.
 * Relative error on [0, 30] is documented as having a peak of 5.8e-16.
 *
 * Reference: numpy/lib/_function_base_impl.py lines 3533-3590
 */
export function i0(x: NDArray | number | number[]): NDArray {
  // Handle scalar input
  if (typeof x === 'number') {
    const result = NDArray.zeros([1], DType.Float64);
    result.set(_i0_scalar(x), 0);
    return result.reshape([]);  // Return 0-d array
  }

  // Handle array input
  if (Array.isArray(x)) {
    x = NDArray.fromArray(x, undefined, DType.Float64);
  }

  // Validate dtype
  if (x.dtype === DType.Complex64 || x.dtype === DType.Complex128) {
    throw new TypeError('i0 not supported for complex values');
  }

  // Create output array
  const result = NDArray.empty(x.shape, DType.Float64);

  // Compute I₀ for each element
  const size = x.size;
  for (let i = 0; i < size; i++) {
    const val = x.getFlat(i);
    result.setFlat(i, _i0_scalar(val));
  }

  return result;
}
```

---

### 11.3 Cosine-Based Windows

**File:** `src/ts/window.ts` (continued)

```typescript
/**
 * Return the Blackman window.
 *
 * The Blackman window is a taper formed by using the first three
 * terms of a summation of cosines. It was designed to have close to
 * the minimal leakage possible.
 *
 * Formula: w(n) = 0.42 + 0.5*cos(π*n/(M-1)) + 0.08*cos(2π*n/(M-1))
 *
 * @param M - Number of points in the output window. If zero or less,
 *            an empty array is returned.
 * @returns The window, with the maximum value normalized to one
 *          (appears only if M is odd)
 *
 * @example
 * ```typescript
 * blackman(12)
 * // array([-1.39e-17, 0.0326, 0.1599, 0.4144, 0.7360, 0.9670,
 * //         0.9670, 0.7360, 0.4144, 0.1599, 0.0326, -1.39e-17])
 * ```
 *
 * Reference: numpy/lib/_function_base_impl.py lines 3037-3133
 */
export function blackman(M: number): NDArray {
  M = Math.floor(M);

  if (M < 1) {
    return NDArray.empty([0], DType.Float64);
  }
  if (M === 1) {
    return NDArray.ones([1], DType.Float64);
  }

  // n = arange(1 - M, M, 2) gives symmetric indices centered at 0
  // This creates values: [1-M, 3-M, 5-M, ..., M-3, M-1]
  const n = NDArray.arange(1 - M, M, 2, DType.Float64);
  const denom = M - 1;

  // w(n) = 0.42 + 0.5*cos(π*n/(M-1)) + 0.08*cos(2π*n/(M-1))
  // Note: NumPy uses + 0.5*cos and + 0.08*cos because n is centered at 0
  const result = NDArray.empty([M], DType.Float64);
  for (let i = 0; i < M; i++) {
    const ni = n.getFlat(i);
    const term1 = 0.42;
    const term2 = 0.5 * Math.cos(Math.PI * ni / denom);
    const term3 = 0.08 * Math.cos(2.0 * Math.PI * ni / denom);
    result.setFlat(i, term1 + term2 + term3);
  }

  return result;
}

/**
 * Return the Hanning window (also known as Hann window).
 *
 * The Hanning window is a taper formed by using a weighted cosine.
 * Named for Julius von Hann.
 *
 * Formula: w(n) = 0.5 + 0.5*cos(π*n/(M-1))
 *
 * @param M - Number of points in the output window. If zero or less,
 *            an empty array is returned.
 * @returns The window, with the maximum value normalized to one
 *          (appears only if M is odd)
 *
 * @example
 * ```typescript
 * hanning(12)
 * // array([0, 0.0794, 0.2923, 0.5712, 0.8274, 0.9797,
 * //        0.9797, 0.8274, 0.5712, 0.2923, 0.0794, 0])
 * ```
 *
 * Reference: numpy/lib/_function_base_impl.py lines 3244-3342
 */
export function hanning(M: number): NDArray {
  M = Math.floor(M);

  if (M < 1) {
    return NDArray.empty([0], DType.Float64);
  }
  if (M === 1) {
    return NDArray.ones([1], DType.Float64);
  }

  const n = NDArray.arange(1 - M, M, 2, DType.Float64);
  const denom = M - 1;

  // w(n) = 0.5 + 0.5*cos(π*n/(M-1))
  const result = NDArray.empty([M], DType.Float64);
  for (let i = 0; i < M; i++) {
    const ni = n.getFlat(i);
    result.setFlat(i, 0.5 + 0.5 * Math.cos(Math.PI * ni / denom));
  }

  return result;
}

/**
 * Return the Hamming window.
 *
 * The Hamming window is a taper formed by using a weighted cosine,
 * similar to Hanning but with raised endpoints to reduce the first
 * side lobe.
 *
 * Formula: w(n) = 0.54 + 0.46*cos(π*n/(M-1))
 *
 * @param M - Number of points in the output window. If zero or less,
 *            an empty array is returned.
 * @returns The window, with the maximum value normalized to one
 *          (appears only if M is odd)
 *
 * @example
 * ```typescript
 * hamming(12)
 * // array([0.08, 0.153, 0.349, 0.605, 0.841, 0.981,
 * //        0.981, 0.841, 0.605, 0.349, 0.153, 0.08])
 * ```
 *
 * Reference: numpy/lib/_function_base_impl.py lines 3346-3441
 */
export function hamming(M: number): NDArray {
  M = Math.floor(M);

  if (M < 1) {
    return NDArray.empty([0], DType.Float64);
  }
  if (M === 1) {
    return NDArray.ones([1], DType.Float64);
  }

  const n = NDArray.arange(1 - M, M, 2, DType.Float64);
  const denom = M - 1;

  // w(n) = 0.54 + 0.46*cos(π*n/(M-1))
  const result = NDArray.empty([M], DType.Float64);
  for (let i = 0; i < M; i++) {
    const ni = n.getFlat(i);
    result.setFlat(i, 0.54 + 0.46 * Math.cos(Math.PI * ni / denom));
  }

  return result;
}
```

---

### 11.4 Other Windows

**File:** `src/ts/window.ts` (continued)

```typescript
/**
 * Return the Bartlett window.
 *
 * The Bartlett window is very similar to a triangular window, except
 * that the end points are at zero. Also known as the triangular or
 * Fejér window.
 *
 * Formula: w(n) = 1 - |n - (M-1)/2| / ((M-1)/2)
 *               = (2/(M-1)) * ((M-1)/2 - |n - (M-1)/2|)
 *
 * For the symmetric form with n centered at 0:
 *   w(n) = 1 + n/(M-1) for n <= 0
 *   w(n) = 1 - n/(M-1) for n > 0
 *
 * @param M - Number of points in the output window. If zero or less,
 *            an empty array is returned.
 * @returns The triangular window, with the maximum value normalized to one
 *          (appears only if M is odd), with first and last samples equal
 *          to zero.
 *
 * @example
 * ```typescript
 * bartlett(12)
 * // array([0, 0.182, 0.364, 0.545, 0.727, 0.909,
 * //        0.909, 0.727, 0.545, 0.364, 0.182, 0])
 * ```
 *
 * Reference: numpy/lib/_function_base_impl.py lines 3137-3240
 */
export function bartlett(M: number): NDArray {
  M = Math.floor(M);

  if (M < 1) {
    return NDArray.empty([0], DType.Float64);
  }
  if (M === 1) {
    return NDArray.ones([1], DType.Float64);
  }

  const n = NDArray.arange(1 - M, M, 2, DType.Float64);
  const denom = M - 1;

  // w(n) = 1 + n/(M-1) for n <= 0
  // w(n) = 1 - n/(M-1) for n > 0
  const result = NDArray.empty([M], DType.Float64);
  for (let i = 0; i < M; i++) {
    const ni = n.getFlat(i);
    if (ni <= 0) {
      result.setFlat(i, 1 + ni / denom);
    } else {
      result.setFlat(i, 1 - ni / denom);
    }
  }

  return result;
}

/**
 * Return the Kaiser window.
 *
 * The Kaiser window is a taper formed by using a Bessel function.
 * It provides a good approximation to the Digital Prolate Spheroidal
 * Sequence (DPSS/Slepian window).
 *
 * Formula: w(n) = I₀(β * sqrt(1 - (2n/(M-1))²)) / I₀(β)
 *
 * The beta parameter controls the window shape:
 * - β = 0: Rectangular window
 * - β ≈ 5: Similar to Hamming
 * - β ≈ 6: Similar to Hanning
 * - β ≈ 8.6: Similar to Blackman
 * - β = 14: Good starting point for most applications
 *
 * @param M - Number of points in the output window. If zero or less,
 *            an empty array is returned.
 * @param beta - Shape parameter for the window. Larger values produce
 *               narrower main lobe with lower side lobes.
 * @returns The window, with the maximum value normalized to one
 *          (appears only if M is odd)
 *
 * @example
 * ```typescript
 * kaiser(12, 14)
 * // array([7.73e-06, 0.00346, 0.0465, 0.230, 0.600, 0.946,
 * //        0.946, 0.600, 0.230, 0.0465, 0.00346, 7.73e-06])
 * ```
 *
 * @remarks
 * As beta increases, the window narrows. If beta is too large relative
 * to M, NaN values may be returned.
 *
 * Reference: numpy/lib/_function_base_impl.py lines 3596-3723
 */
export function kaiser(M: number, beta: number): NDArray {
  M = Math.floor(M);

  if (M < 1) {
    return NDArray.empty([0], DType.Float64);
  }
  if (M === 1) {
    return NDArray.ones([1], DType.Float64);
  }

  const alpha = (M - 1) / 2.0;
  const i0Beta = _i0_scalar(beta);

  // w(n) = I₀(β * sqrt(1 - ((n - α) / α)²)) / I₀(β)
  const result = NDArray.empty([M], DType.Float64);
  for (let n = 0; n < M; n++) {
    const x = (n - alpha) / alpha;
    const arg = 1 - x * x;
    // Handle numerical issues where arg might be slightly negative
    const sqrtArg = arg >= 0 ? Math.sqrt(arg) : 0;
    result.setFlat(n, _i0_scalar(beta * sqrtArg) / i0Beta);
  }

  return result;
}
```

---

## File Changes Summary

### New Files to Create

```
src/ts/
└── window.ts          # All window functions and i0 Bessel function
```

### Files to Modify

```
src/ts/index.ts
├── Import window functions
└── Export: blackman, bartlett, hanning, hamming, kaiser, i0
```

---

## Export Updates

**File:** `src/ts/index.ts` (additions)

```typescript
// Window functions (Phase 11)
export {
  blackman,
  bartlett,
  hanning,
  hamming,
  kaiser,
  i0,
} from './window.js';
```

---

## Implementation Order

```
Phase 11 Implementation:
│
├── Step 1: Create window.ts with i0 Bessel function
│   ├── Implement _chbevl() Chebyshev evaluation
│   ├── Add _i0A and _i0B coefficient arrays
│   ├── Implement _i0_1() and _i0_2() helpers
│   ├── Implement _i0_scalar() for single values
│   └── Implement i0() main function
│
├── Step 2: Implement cosine-based windows
│   ├── Implement blackman(M)
│   ├── Implement hanning(M)
│   └── Implement hamming(M)
│
├── Step 3: Implement other windows
│   ├── Implement bartlett(M)
│   └── Implement kaiser(M, beta)
│
├── Step 4: Update exports
│   └── Add exports to index.ts
│
├── Step 5: Write tests
│   ├── Test i0() against NumPy values
│   ├── Test each window function for M=0, M=1, M=12, M=51
│   ├── Test symmetry properties
│   └── Test normalization (max value = 1 for odd M)
│
└── Step 6: Documentation
    └── Add JSDoc comments with mathematical formulas
```

---

## Test Vectors

Generate NumPy comparison vectors:

```python
# tests/python/generate_phase11_tests.py
import numpy as np
import json

def to_list(arr):
    """Convert numpy array to JSON-serializable list."""
    return [float(x) for x in arr]

tests = {
    "i0": [
        {"input": 0.0, "expected": 1.0},
        {"input": 1.0, "expected": 1.2660658777520082},
        {"input": 2.0, "expected": 2.2795853023360673},
        {"input": 3.0, "expected": 4.880792585865024},
        {"input": 8.0, "expected": 427.56411572180454},  # Boundary
        {"input": 10.0, "expected": 2815.716628466254},  # > 8 branch
        {"input": -5.0, "expected": 27.239871823604442},  # Negative (uses abs)
    ],
    "blackman": [
        {"M": 0, "expected": []},
        {"M": 1, "expected": [1.0]},
        {"M": 12, "expected": to_list(np.blackman(12))},
        {"M": 51, "expected": to_list(np.blackman(51))},
    ],
    "bartlett": [
        {"M": 0, "expected": []},
        {"M": 1, "expected": [1.0]},
        {"M": 12, "expected": to_list(np.bartlett(12))},
        {"M": 51, "expected": to_list(np.bartlett(51))},
    ],
    "hanning": [
        {"M": 0, "expected": []},
        {"M": 1, "expected": [1.0]},
        {"M": 12, "expected": to_list(np.hanning(12))},
        {"M": 51, "expected": to_list(np.hanning(51))},
    ],
    "hamming": [
        {"M": 0, "expected": []},
        {"M": 1, "expected": [1.0]},
        {"M": 12, "expected": to_list(np.hamming(12))},
        {"M": 51, "expected": to_list(np.hamming(51))},
    ],
    "kaiser": [
        {"M": 0, "beta": 14, "expected": []},
        {"M": 1, "beta": 14, "expected": [1.0]},
        {"M": 12, "beta": 0, "expected": to_list(np.kaiser(12, 0))},  # Rectangular
        {"M": 12, "beta": 5, "expected": to_list(np.kaiser(12, 5))},  # ~Hamming
        {"M": 12, "beta": 6, "expected": to_list(np.kaiser(12, 6))},  # ~Hanning
        {"M": 12, "beta": 8.6, "expected": to_list(np.kaiser(12, 8.6))},  # ~Blackman
        {"M": 12, "beta": 14, "expected": to_list(np.kaiser(12, 14))},
        {"M": 51, "beta": 14, "expected": to_list(np.kaiser(51, 14))},
    ],
}

# Save test vectors
with open("tests/fixtures/phase11_vectors.json", "w") as f:
    json.dump(tests, f, indent=2)

print("Test vectors generated successfully!")
```

---

## Verification Plan

After Phase 11 completion, verify:

```bash
# Build
npm run build

# Run all tests
npm test

# Specific Phase 11 tests should pass:

# i0 Bessel function
✓ i0(0) returns 1.0
✓ i0([0, 1, 2, 3]) matches NumPy values
✓ i0() uses abs(x) for negative inputs
✓ i0() switches algorithm at |x| = 8 boundary
✓ i0() throws for complex input

# Cosine-based windows
✓ blackman(0) returns empty array
✓ blackman(1) returns [1.0]
✓ blackman(12) matches NumPy within 1e-14
✓ blackman(M) is symmetric
✓ hanning(12) matches NumPy within 1e-14
✓ hamming(12) matches NumPy within 1e-14

# Other windows
✓ bartlett(12) matches NumPy within 1e-14
✓ bartlett endpoints are 0
✓ kaiser(12, 14) matches NumPy within 1e-14
✓ kaiser(12, 0) approximates rectangular window
✓ kaiser(12, 5) approximates hamming window
```

---

## Performance Considerations

1. **Pure TypeScript Implementation**: All window functions use pure TypeScript with Math.* functions, which is efficient for typical window sizes (< 10000 points).

2. **No Ufunc Dependency**: Unlike NumPy's implementation which uses array operations, this implementation uses scalar loops. This avoids the Phase 4 ufunc dependency while maintaining numerical accuracy.

3. **Bessel Function**: The I₀ implementation uses the same Chebyshev polynomial approximation as NumPy/Cephes, providing excellent accuracy (relative error < 5.8e-16).

4. **Memory**: Each function allocates a single output array of size M, with no intermediate allocations.

---

## Dependencies for Later Phases

Phase 11 completion enables:

- **numpy.fft (Phase 14)**: Window functions are commonly used with FFT for spectral analysis
- **scipy.signal equivalent**: Windowed filtering and spectral estimation
- **Audio processing**: Window functions are essential for audio analysis

---

## Notes on NumPy Compatibility

1. **Floating Point Precision**: NumPy coerces M to float64 for consistent behavior. This implementation uses `Math.floor(M)` for integer handling.

2. **Empty vs Zero**: NumPy returns empty array for M < 1, not an array of zeros.

3. **Symmetry**: All windows are symmetric around the center point `(M-1)/2`.

4. **Endpoint Values**:
   - Blackman, Hanning: Endpoints approach 0 (exact 0 for hanning)
   - Hamming: Endpoints are 0.08 (raised)
   - Bartlett: Endpoints are exactly 0
   - Kaiser: Endpoints depend on beta

5. **Maximum Value**: The maximum value is 1.0, occurring at the center only when M is odd.
