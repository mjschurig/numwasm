# Phase 22: Cumulative Operations Implementation Plan

Complete implementation roadmap for cumulative sum, product, and their NaN-handling variants.

---

## ⚠️ Implementation Guideline

**The original NumPy source code is available in `/numpy`.**

Key reference files:
- `numpy/_core/fromnumeric.py` - `cumsum`, `cumprod` functions (~4,233 lines)
- `numpy/lib/_nanfunctions_impl.py` - `nancumsum`, `nancumprod` functions
- `numpy/_core/src/multiarray/calculation.c` - C implementation

Implementation should follow NumPy's algorithms and error handling for consistency.

---

## Current State (Pre-Phase 22)

```
Already Implemented:
├── sum(arr, axis, dtype, keepdims)
├── prod(arr, axis, dtype, keepdims)
├── diff(arr, n, axis, prepend, append)
└── ediff1d(arr, to_end, to_begin)

Missing:
├── cumsum(arr, axis, dtype, out)
├── cumprod(arr, axis, dtype, out)
├── nancumsum(arr, axis, dtype)
└── nancumprod(arr, axis, dtype)
```

---

## Phase 22 Dependency Tree

```
PHASE 22: CUMULATIVE OPERATIONS
│
├── 22.1 Core Cumulative Functions (TypeScript + WASM)
│   ├── 22.1.1 cumsum(a, axis, dtype, out)
│   │   ├── Axis handling (None = flatten)
│   │   ├── dtype promotion
│   │   └── Output array allocation
│   │
│   └── 22.1.2 cumprod(a, axis, dtype, out)
│       ├── Axis handling (None = flatten)
│       ├── dtype promotion
│       └── Output array allocation
│
│   Dependencies: NDArray core, dtype system, axis iteration
│
├── 22.2 NaN-handling Variants (TypeScript)
│   ├── 22.2.1 nancumsum(a, axis, dtype)
│   │   └── Replace NaN with 0 before cumsum
│   │
│   └── 22.2.2 nancumprod(a, axis, dtype)
│       └── Replace NaN with 1 before cumprod
│
│   Dependencies: 22.1.* (core cumulative functions)
│
└── 22.3 WASM Acceleration (C)
    ├── 22.3.1 cumsum_contiguous - fast path
    ├── 22.3.2 cumprod_contiguous - fast path
    └── 22.3.3 cumulative_along_axis - strided access

    Dependencies: NDArray WASM core
```

---

## Detailed Implementation Specifications

### 22.1 Core Cumulative Functions

#### 22.1.1 cumsum

**File:** `src/ts/statistics.ts` (additions)

```typescript
/**
 * Return the cumulative sum of the elements along a given axis.
 *
 * @param a - Input array
 * @param axis - Axis along which to compute. None = flatten first
 * @param dtype - Type of output array. If not specified, uses input dtype
 *                (or float64 for integer inputs to avoid overflow)
 * @param out - Output array (must have correct shape)
 * @returns Cumulative sum array with same shape as input (or flattened if axis=null)
 *
 * @example
 * cumsum([1, 2, 3, 4])  // [1, 3, 6, 10]
 *
 * cumsum([[1, 2], [3, 4]], axis=0)  // [[1, 2], [4, 6]]
 * cumsum([[1, 2], [3, 4]], axis=1)  // [[1, 3], [3, 7]]
 */
export function cumsum(
  a: NDArray | ArrayLike<number>,
  axis: number | null = null,
  dtype: DType | null = null,
  out: NDArray | null = null
): NDArray {
  const arr = asarray(a);

  // If axis is null, flatten the array
  if (axis === null) {
    const flat = arr.ravel();
    return _cumsumFlat(flat, dtype, out);
  }

  // Normalize negative axis
  const normalizedAxis = normalizeAxis(axis, arr.ndim);

  // Determine output dtype
  const outDtype = dtype ?? _cumsumDtype(arr.dtype);

  // Allocate or validate output
  const result = out ?? empty(arr.shape, outDtype);

  if (out !== null) {
    if (!arraysEqual(out.shape, arr.shape)) {
      throw new ValueError(
        `output array has wrong shape: expected ${arr.shape}, got ${out.shape}`
      );
    }
  }

  // Call WASM implementation
  _cumsumAxis(arr, result, normalizedAxis);

  return result;
}

/**
 * Determine appropriate dtype for cumsum output.
 * Integer types are promoted to avoid overflow.
 */
function _cumsumDtype(inputDtype: DType): DType {
  // Promote small integer types to larger ones
  switch (inputDtype) {
    case DType.Bool:
    case DType.Int8:
    case DType.Int16:
    case DType.Int32:
      return DType.Int64;
    case DType.UInt8:
    case DType.UInt16:
    case DType.UInt32:
      return DType.UInt64;
    default:
      return inputDtype;
  }
}
```

#### 22.1.2 cumprod

**File:** `src/ts/statistics.ts` (additions)

```typescript
/**
 * Return the cumulative product of elements along a given axis.
 *
 * @param a - Input array
 * @param axis - Axis along which to compute. None = flatten first
 * @param dtype - Type of output array
 * @param out - Output array
 * @returns Cumulative product array
 *
 * @example
 * cumprod([1, 2, 3, 4])  // [1, 2, 6, 24]
 *
 * cumprod([[1, 2], [3, 4]], axis=0)  // [[1, 2], [3, 8]]
 * cumprod([[1, 2], [3, 4]], axis=1)  // [[1, 2], [3, 12]]
 */
export function cumprod(
  a: NDArray | ArrayLike<number>,
  axis: number | null = null,
  dtype: DType | null = null,
  out: NDArray | null = null
): NDArray {
  const arr = asarray(a);

  if (axis === null) {
    const flat = arr.ravel();
    return _cumprodFlat(flat, dtype, out);
  }

  const normalizedAxis = normalizeAxis(axis, arr.ndim);
  const outDtype = dtype ?? _cumprodDtype(arr.dtype);
  const result = out ?? empty(arr.shape, outDtype);

  if (out !== null) {
    if (!arraysEqual(out.shape, arr.shape)) {
      throw new ValueError(
        `output array has wrong shape: expected ${arr.shape}, got ${out.shape}`
      );
    }
  }

  _cumprodAxis(arr, result, normalizedAxis);

  return result;
}
```

---

### 22.2 NaN-handling Variants

**File:** `src/ts/statistics.ts` (additions)

```typescript
/**
 * Return the cumulative sum of array elements over a given axis
 * treating Not a Numbers (NaNs) as zero.
 *
 * @param a - Input array
 * @param axis - Axis along which to compute
 * @param dtype - Type of output array
 * @returns Cumulative sum with NaN replaced by zero
 *
 * @example
 * nancumsum([1, NaN, 3, 4])  // [1, 1, 4, 8]
 */
export function nancumsum(
  a: NDArray | ArrayLike<number>,
  axis: number | null = null,
  dtype: DType | null = null
): NDArray {
  const arr = asarray(a);

  // Create a copy with NaN replaced by 0
  const cleaned = where(isnan(arr), 0, arr);

  return cumsum(cleaned, axis, dtype);
}

/**
 * Return the cumulative product of array elements over a given axis
 * treating Not a Numbers (NaNs) as one.
 *
 * @param a - Input array
 * @param axis - Axis along which to compute
 * @param dtype - Type of output array
 * @returns Cumulative product with NaN replaced by one
 *
 * @example
 * nancumprod([1, NaN, 3, 4])  // [1, 1, 3, 12]
 */
export function nancumprod(
  a: NDArray | ArrayLike<number>,
  axis: number | null = null,
  dtype: DType | null = null
): NDArray {
  const arr = asarray(a);

  // Create a copy with NaN replaced by 1
  const cleaned = where(isnan(arr), 1, arr);

  return cumprod(cleaned, axis, dtype);
}
```

---

### 22.3 WASM Acceleration

**File:** `src/wasm/cumulative.c` (new file)

```c
#ifndef NUMJS_CUMULATIVE_H
#define NUMJS_CUMULATIVE_H

#include "ndarray.h"

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#define EXPORT EMSCRIPTEN_KEEPALIVE
#else
#define EXPORT
#endif

/**
 * Cumulative sum for contiguous 1D array (fast path).
 *
 * @param n      Number of elements
 * @param input  Input array
 * @param output Output array (must be pre-allocated)
 * @param dtype  Data type
 */
EXPORT void cumsum_contiguous_f64(int32_t n, const double* input, double* output);
EXPORT void cumsum_contiguous_f32(int32_t n, const float* input, float* output);
EXPORT void cumsum_contiguous_i32(int32_t n, const int32_t* input, int64_t* output);
EXPORT void cumsum_contiguous_i64(int32_t n, const int64_t* input, int64_t* output);

/**
 * Cumulative product for contiguous 1D array (fast path).
 */
EXPORT void cumprod_contiguous_f64(int32_t n, const double* input, double* output);
EXPORT void cumprod_contiguous_f32(int32_t n, const float* input, float* output);
EXPORT void cumprod_contiguous_i32(int32_t n, const int32_t* input, int64_t* output);
EXPORT void cumprod_contiguous_i64(int32_t n, const int64_t* input, int64_t* output);

/**
 * Cumulative sum along axis with strided access.
 *
 * @param ndim        Number of dimensions
 * @param shape       Array shape
 * @param input       Input data pointer
 * @param in_strides  Input strides (in bytes)
 * @param output      Output data pointer
 * @param out_strides Output strides (in bytes)
 * @param axis        Axis along which to compute
 * @param dtype       Data type
 */
EXPORT void cumsum_axis(
    int32_t ndim,
    const int32_t* shape,
    const void* input,
    const int32_t* in_strides,
    void* output,
    const int32_t* out_strides,
    int32_t axis,
    int32_t dtype
);

/**
 * Cumulative product along axis with strided access.
 */
EXPORT void cumprod_axis(
    int32_t ndim,
    const int32_t* shape,
    const void* input,
    const int32_t* in_strides,
    void* output,
    const int32_t* out_strides,
    int32_t axis,
    int32_t dtype
);

#endif /* NUMJS_CUMULATIVE_H */
```

**File:** `src/wasm/cumulative.c` (implementation)

```c
#include "cumulative.h"
#include <string.h>

/* ============ Contiguous Fast Paths ============ */

EXPORT void cumsum_contiguous_f64(int32_t n, const double* input, double* output) {
    if (n <= 0) return;

    double sum = 0.0;
    for (int32_t i = 0; i < n; i++) {
        sum += input[i];
        output[i] = sum;
    }
}

EXPORT void cumsum_contiguous_f32(int32_t n, const float* input, float* output) {
    if (n <= 0) return;

    float sum = 0.0f;
    for (int32_t i = 0; i < n; i++) {
        sum += input[i];
        output[i] = sum;
    }
}

EXPORT void cumsum_contiguous_i64(int32_t n, const int64_t* input, int64_t* output) {
    if (n <= 0) return;

    int64_t sum = 0;
    for (int32_t i = 0; i < n; i++) {
        sum += input[i];
        output[i] = sum;
    }
}

EXPORT void cumprod_contiguous_f64(int32_t n, const double* input, double* output) {
    if (n <= 0) return;

    double prod = 1.0;
    for (int32_t i = 0; i < n; i++) {
        prod *= input[i];
        output[i] = prod;
    }
}

EXPORT void cumprod_contiguous_f32(int32_t n, const float* input, float* output) {
    if (n <= 0) return;

    float prod = 1.0f;
    for (int32_t i = 0; i < n; i++) {
        prod *= input[i];
        output[i] = prod;
    }
}

/* ============ Strided Axis Implementation ============ */

EXPORT void cumsum_axis(
    int32_t ndim,
    const int32_t* shape,
    const void* input,
    const int32_t* in_strides,
    void* output,
    const int32_t* out_strides,
    int32_t axis,
    int32_t dtype
) {
    /* Calculate the number of "lanes" to process
     * (all dimensions except the cumsum axis) */
    int64_t n_lanes = 1;
    for (int32_t d = 0; d < ndim; d++) {
        if (d != axis) {
            n_lanes *= shape[d];
        }
    }

    int32_t axis_len = shape[axis];
    int32_t axis_in_stride = in_strides[axis];
    int32_t axis_out_stride = out_strides[axis];

    /* Process each lane */
    int32_t* indices = (int32_t*)calloc(ndim, sizeof(int32_t));

    for (int64_t lane = 0; lane < n_lanes; lane++) {
        /* Calculate base offset for this lane */
        int64_t in_offset = 0;
        int64_t out_offset = 0;
        for (int32_t d = 0; d < ndim; d++) {
            if (d != axis) {
                in_offset += indices[d] * in_strides[d];
                out_offset += indices[d] * out_strides[d];
            }
        }

        /* Perform cumsum along axis */
        if (dtype == 10) { /* Float64 */
            double sum = 0.0;
            const char* in_ptr = (const char*)input + in_offset;
            char* out_ptr = (char*)output + out_offset;

            for (int32_t i = 0; i < axis_len; i++) {
                sum += *(const double*)(in_ptr + i * axis_in_stride);
                *(double*)(out_ptr + i * axis_out_stride) = sum;
            }
        } else if (dtype == 9) { /* Float32 */
            float sum = 0.0f;
            const char* in_ptr = (const char*)input + in_offset;
            char* out_ptr = (char*)output + out_offset;

            for (int32_t i = 0; i < axis_len; i++) {
                sum += *(const float*)(in_ptr + i * axis_in_stride);
                *(float*)(out_ptr + i * axis_out_stride) = sum;
            }
        }
        /* Add more dtype cases as needed */

        /* Increment indices (skip axis dimension) */
        for (int32_t d = ndim - 1; d >= 0; d--) {
            if (d == axis) continue;
            indices[d]++;
            if (indices[d] < shape[d]) break;
            indices[d] = 0;
        }
    }

    free(indices);
}

/* Similar implementation for cumprod_axis */
EXPORT void cumprod_axis(
    int32_t ndim,
    const int32_t* shape,
    const void* input,
    const int32_t* in_strides,
    void* output,
    const int32_t* out_strides,
    int32_t axis,
    int32_t dtype
) {
    int64_t n_lanes = 1;
    for (int32_t d = 0; d < ndim; d++) {
        if (d != axis) {
            n_lanes *= shape[d];
        }
    }

    int32_t axis_len = shape[axis];
    int32_t axis_in_stride = in_strides[axis];
    int32_t axis_out_stride = out_strides[axis];

    int32_t* indices = (int32_t*)calloc(ndim, sizeof(int32_t));

    for (int64_t lane = 0; lane < n_lanes; lane++) {
        int64_t in_offset = 0;
        int64_t out_offset = 0;
        for (int32_t d = 0; d < ndim; d++) {
            if (d != axis) {
                in_offset += indices[d] * in_strides[d];
                out_offset += indices[d] * out_strides[d];
            }
        }

        if (dtype == 10) { /* Float64 */
            double prod = 1.0;
            const char* in_ptr = (const char*)input + in_offset;
            char* out_ptr = (char*)output + out_offset;

            for (int32_t i = 0; i < axis_len; i++) {
                prod *= *(const double*)(in_ptr + i * axis_in_stride);
                *(double*)(out_ptr + i * axis_out_stride) = prod;
            }
        } else if (dtype == 9) { /* Float32 */
            float prod = 1.0f;
            const char* in_ptr = (const char*)input + in_offset;
            char* out_ptr = (char*)output + out_offset;

            for (int32_t i = 0; i < axis_len; i++) {
                prod *= *(const float*)(in_ptr + i * axis_in_stride);
                *(float*)(out_ptr + i * axis_out_stride) = prod;
            }
        }

        for (int32_t d = ndim - 1; d >= 0; d--) {
            if (d == axis) continue;
            indices[d]++;
            if (indices[d] < shape[d]) break;
            indices[d] = 0;
        }
    }

    free(indices);
}
```

---

## File Changes Summary

### New Files to Create

```
src/wasm/
└── cumulative.c         # WASM cumulative operations

tests/ts/
└── cumulative.test.ts   # Test suite
```

### Files to Modify

```
src/ts/statistics.ts
├── Add cumsum() function
├── Add cumprod() function
├── Add nancumsum() function
└── Add nancumprod() function

src/ts/index.ts
├── Export cumsum
├── Export cumprod
├── Export nancumsum
└── Export nancumprod

src/ts/types.ts
└── Add WASM function declarations

scripts/build-wasm.sh
├── Add cumulative.c to compilation
└── Add EXPORTED_FUNCTIONS
```

---

## Build Script Updates

Add to `scripts/build-wasm.sh` EXPORTED_FUNCTIONS:

```bash
"_cumsum_contiguous_f64",
"_cumsum_contiguous_f32",
"_cumsum_contiguous_i64",
"_cumprod_contiguous_f64",
"_cumprod_contiguous_f32",
"_cumprod_contiguous_i64",
"_cumsum_axis",
"_cumprod_axis"
```

---

## Implementation Order

```
Phase 22.1: Core Functions (Day 1-2)
├── Day 1: cumsum implementation
│   ├── TypeScript function with axis handling
│   ├── WASM contiguous fast path
│   └── Basic tests
│
└── Day 2: cumprod implementation
    ├── TypeScript function with axis handling
    ├── WASM contiguous fast path
    └── Basic tests

Phase 22.2: NaN Variants (Day 3)
├── nancumsum implementation
├── nancumprod implementation
└── Tests with NaN values

Phase 22.3: WASM Optimization (Day 4)
├── Strided axis implementation
├── All dtype support
└── Performance benchmarks

Phase 22.4: Testing & Polish (Day 5)
├── Edge cases (empty arrays, single element)
├── NumPy comparison tests
└── Documentation
```

---

## Verification Plan

After Phase 22 completion, verify:

```bash
# Build
npm run build

# Run tests
npm test

# Phase 22 specific tests:

# Basic functionality
✓ cumsum([1, 2, 3, 4]) === [1, 3, 6, 10]
✓ cumprod([1, 2, 3, 4]) === [1, 2, 6, 24]

# Axis support
✓ cumsum([[1, 2], [3, 4]], axis=0) === [[1, 2], [4, 6]]
✓ cumsum([[1, 2], [3, 4]], axis=1) === [[1, 3], [3, 7]]

# NaN handling
✓ nancumsum([1, NaN, 3]) === [1, 1, 4]
✓ nancumprod([1, NaN, 3]) === [1, 1, 3]

# Edge cases
✓ cumsum([]) === []
✓ cumsum([5]) === [5]
✓ cumsum on non-contiguous view works correctly

# dtype handling
✓ cumsum of int32 returns int64 (to prevent overflow)
✓ cumsum with explicit dtype=float32 works
```

Generate NumPy comparison vectors:

```python
import numpy as np
import json

tests = {
    "cumsum_1d": {
        "input": [1, 2, 3, 4, 5],
        "expected": np.cumsum([1, 2, 3, 4, 5]).tolist()
    },
    "cumsum_2d_axis0": {
        "input": [[1, 2, 3], [4, 5, 6]],
        "axis": 0,
        "expected": np.cumsum([[1, 2, 3], [4, 5, 6]], axis=0).tolist()
    },
    "cumsum_2d_axis1": {
        "input": [[1, 2, 3], [4, 5, 6]],
        "axis": 1,
        "expected": np.cumsum([[1, 2, 3], [4, 5, 6]], axis=1).tolist()
    },
    "cumprod_1d": {
        "input": [1, 2, 3, 4],
        "expected": np.cumprod([1, 2, 3, 4]).tolist()
    },
    "nancumsum": {
        "input": [1.0, float('nan'), 3.0, 4.0],
        "expected": np.nancumsum([1.0, np.nan, 3.0, 4.0]).tolist()
    },
    "nancumprod": {
        "input": [1.0, float('nan'), 3.0, 4.0],
        "expected": np.nancumprod([1.0, np.nan, 3.0, 4.0]).tolist()
    }
}

with open("tests/fixtures/cumulative_vectors.json", "w") as f:
    json.dump(tests, f, indent=2)
```

---

## API Compatibility Notes

### NumPy Signature Match

```typescript
// NumPy: numpy.cumsum(a, axis=None, dtype=None, out=None)
// NumJS: cumsum(a, axis=null, dtype=null, out=null)

// NumPy: numpy.cumprod(a, axis=None, dtype=None, out=None)
// NumJS: cumprod(a, axis=null, dtype=null, out=null)

// NumPy: numpy.nancumsum(a, axis=None, dtype=None, out=None)
// NumJS: nancumsum(a, axis=null, dtype=null)  // out not supported initially

// NumPy: numpy.nancumprod(a, axis=None, dtype=None, out=None)
// NumJS: nancumprod(a, axis=null, dtype=null)  // out not supported initially
```

### Differences from NumPy

1. **Initial `out` parameter**: The `out` parameter for `nancumsum`/`nancumprod` may be added in a future version.

2. **Integer overflow**: NumPy allows integer overflow; NumJS promotes to int64 by default for safety.
