# Level 3: Universal Functions (Ufuncs) Implementation Plan

Level 3 implements NumPy's Universal Functions (ufuncs) infrastructure, enabling element-wise operations, broadcasting for binary operations, and axis-based reductions. This is the critical foundation for all mathematical operations.

---

## Current State (Level 2 Complete)

```
src/wasm/
├── ndarray.h/c        # Core array with views, slicing, element access
├── dtype.h/c          # DType utilities, type promotion
├── broadcast.h/c      # Broadcasting shape/stride computation
├── indexing.h/c       # take, put, nonzero, where
└── pairwise_sum.h/c   # Accurate pairwise summation algorithm

src/ts/
├── types.ts           # DType enum, WasmModule interface, flags
├── NDArray.ts         # Factory methods, slicing, shape manipulation
├── dtype.ts           # Type promotion (promoteTypes, canCast, etc.)
├── broadcast.ts       # broadcastTo, broadcastArrays, broadcastShapes
├── slice.ts           # Slice class, ellipsis, newaxis
├── indexing.ts        # take, put, nonzero, where, compress, extract
├── iterators.ts       # FlatIterator, nditer, ndenumerate, ndindex
└── index.ts           # Public exports
```

**Existing Infrastructure:**
- Broadcasting: `broadcast_shapes()`, `broadcast_strides()`, `ndarray_broadcast_to()`
- Type promotion: `dtype_promote()`, `promoteTypes()`, `canCast()`
- Pairwise summation: `pairwise_sum_*()` for accurate reduction
- Strided iteration: Element access respects shape/strides
- Flags: OWNDATA, WRITEABLE, C_CONTIGUOUS, F_CONTIGUOUS

---

## Level 3 Implementation Tree

```
LEVEL 3: UNIVERSAL FUNCTIONS (UFUNCS)
│
├── 3.1 Ufunc Infrastructure
│   ├── 3.1.1 C: UfuncDef struct and registry
│   ├── 3.1.2 C: Type resolution system
│   ├── 3.1.3 C: Strided iteration helpers
│   ├── 3.1.4 C: Output array allocation
│   ├── 3.1.5 TS: Ufunc class
│   ├── 3.1.6 TS: ufunc registry object
│   └── 3.1.7 TS: WasmModule interface extensions
│
├── 3.2 Unary Ufuncs (C + TypeScript)
│   ├── 3.2.1 Arithmetic: negative, positive, absolute, sign
│   ├── 3.2.2 Powers: sqrt, square, cbrt, reciprocal
│   ├── 3.2.3 Exponential: exp, exp2, expm1
│   ├── 3.2.4 Logarithmic: log, log2, log10, log1p
│   ├── 3.2.5 Trigonometric: sin, cos, tan, arcsin, arccos, arctan
│   ├── 3.2.6 Hyperbolic: sinh, cosh, tanh, arcsinh, arccosh, arctanh
│   ├── 3.2.7 Rounding: floor, ceil, trunc, rint, round
│   ├── 3.2.8 Angle conversion: degrees, radians, deg2rad, rad2deg
│   └── 3.2.9 Predicates: isnan, isinf, isfinite, signbit
│
│   Dependencies: 3.1.*
│
├── 3.3 Binary Ufuncs (C + TypeScript)
│   ├── 3.3.1 Arithmetic: add, subtract, multiply, divide
│   ├── 3.3.2 Division variants: true_divide, floor_divide, remainder, mod, fmod
│   ├── 3.3.3 Powers: power, float_power
│   ├── 3.3.4 Comparison: greater, greater_equal, less, less_equal, equal, not_equal
│   ├── 3.3.5 Extrema: maximum, minimum, fmax, fmin
│   ├── 3.3.6 Logical: logical_and, logical_or, logical_xor, logical_not
│   ├── 3.3.7 Bitwise: bitwise_and, bitwise_or, bitwise_xor, invert
│   ├── 3.3.8 Shifts: left_shift, right_shift
│   ├── 3.3.9 Special: arctan2, hypot, copysign, nextafter
│   └── 3.3.10 Rational: gcd, lcm
│
│   Dependencies: 3.1.*, broadcasting from Level 2
│
├── 3.4 Reduction Operations (C + TypeScript)
│   ├── 3.4.1 C: ReduceContext struct
│   ├── 3.4.2 C: Axis iteration infrastructure
│   ├── 3.4.3 C: ndarray_reduce() with axis support
│   ├── 3.4.4 sum(arr, axis, dtype, keepdims, initial)
│   ├── 3.4.5 prod(arr, axis, dtype, keepdims, initial)
│   ├── 3.4.6 min/max/amax/amin(arr, axis, keepdims)
│   ├── 3.4.7 mean(arr, axis, dtype, keepdims)
│   ├── 3.4.8 argmin/argmax(arr, axis, keepdims)
│   ├── 3.4.9 all/any(arr, axis, keepdims)
│   └── 3.4.10 nansum, nanprod, nanmin, nanmax, nanmean
│
│   Dependencies: 3.1.*, 3.2.*, 3.3.*
│
├── 3.5 Accumulation Operations (C + TypeScript)
│   ├── 3.5.1 C: ndarray_accumulate() infrastructure
│   ├── 3.5.2 cumsum(arr, axis, dtype)
│   ├── 3.5.3 cumprod(arr, axis, dtype)
│   ├── 3.5.4 nancumsum, nancumprod
│   └── 3.5.5 diff(arr, n, axis, prepend, append)
│
│   Dependencies: 3.4.*
│
└── 3.6 Ufunc Methods (TypeScript)
    ├── 3.6.1 ufunc.reduce(arr, axis, dtype, keepdims, initial)
    ├── 3.6.2 ufunc.accumulate(arr, axis, dtype)
    ├── 3.6.3 ufunc.outer(arr1, arr2)
    ├── 3.6.4 ufunc.at(arr, indices, b) [in-place]
    └── 3.6.5 ufunc.reduceat(arr, indices, axis)

    Dependencies: 3.4.*, 3.5.*
```

---

## NumPy Reference Files

| Component | NumPy Reference |
|-----------|----------------|
| Ufunc object | `numpy/_core/src/umath/ufunc_object.c` |
| Reduction wrapper | `numpy/_core/src/umath/reduction.c` |
| Inner loops | `numpy/_core/src/umath/loops.c.src` |
| Type resolution | `numpy/_core/src/umath/ufunc_type_resolution.c` |
| Math functions | `numpy/_core/src/npymath/npy_math.c` |
| Python interface | `numpy/_core/umath.py` |

---

## Detailed Implementation Specifications

### 3.1 Ufunc Infrastructure

#### 3.1.1 C: UfuncDef Struct and Registry

**File:** `src/wasm/ufunc.h` (new file)

```c
#ifndef NUMJS_UFUNC_H
#define NUMJS_UFUNC_H

#include "ndarray.h"
#include "dtype.h"

/* ============ Ufunc Type Constants ============ */

#define UFUNC_TYPE_UNARY   1
#define UFUNC_TYPE_BINARY  2

/* Identity values for reductions */
#define UFUNC_IDENTITY_ZERO    0
#define UFUNC_IDENTITY_ONE     1
#define UFUNC_IDENTITY_NONE   -1  /* No identity (e.g., max on empty) */
#define UFUNC_IDENTITY_REORDERABLE 2  /* Like min/max - reorderable but no identity */

/* ============ Function Pointer Types ============ */

/**
 * Unary inner loop: applies operation element-wise.
 * Handles strided input/output.
 */
typedef void (*UnaryLoopFunc)(
    const char* input,      /* Input data pointer */
    char* output,           /* Output data pointer */
    size_t count,           /* Number of elements */
    int64_t in_stride,      /* Input stride in bytes */
    int64_t out_stride      /* Output stride in bytes */
);

/**
 * Binary inner loop: applies operation element-wise with broadcasting.
 * Both inputs may have different strides (0 for broadcast dims).
 */
typedef void (*BinaryLoopFunc)(
    const char* input1,     /* First input data */
    const char* input2,     /* Second input data */
    char* output,           /* Output data */
    size_t count,           /* Number of elements */
    int64_t in1_stride,     /* First input stride */
    int64_t in2_stride,     /* Second input stride */
    int64_t out_stride      /* Output stride */
);

/**
 * Reduction inner loop: reduces array along contiguous dimension.
 */
typedef void (*ReduceLoopFunc)(
    const char* input,      /* Input data */
    char* output,           /* Output (single value or accumulator) */
    size_t count,           /* Number of elements to reduce */
    int64_t stride,         /* Input stride in bytes */
    const char* identity    /* Identity value for empty reductions */
);

/* ============ Ufunc Definition ============ */

/**
 * Type-specific implementation for a ufunc.
 */
typedef struct {
    DType in_dtype1;        /* First input dtype */
    DType in_dtype2;        /* Second input dtype (DTYPE_NONE for unary) */
    DType out_dtype;        /* Output dtype */
    void* loop_func;        /* UnaryLoopFunc or BinaryLoopFunc */
} UfuncTypeLoop;

/**
 * Complete ufunc definition.
 */
typedef struct {
    const char* name;       /* Ufunc name: "add", "sin", etc. */
    int type;               /* UFUNC_TYPE_UNARY or UFUNC_TYPE_BINARY */
    int nin;                /* Number of inputs (1 or 2) */
    int nout;               /* Number of outputs (always 1 for now) */
    int identity;           /* Identity for reduction */

    /* Type-specific loops */
    UfuncTypeLoop* loops;   /* Array of type-specific implementations */
    int nloops;             /* Number of loops */

    /* Reduction loop (if applicable) */
    ReduceLoopFunc reduce_loop;
} UfuncDef;

/* ============ Ufunc Registry ============ */

/**
 * Get ufunc definition by name.
 * @return Pointer to static UfuncDef or NULL if not found.
 */
const UfuncDef* ufunc_get(const char* name);

/**
 * Get ufunc definition by ID (faster lookup).
 */
const UfuncDef* ufunc_get_by_id(int id);

/* Ufunc IDs for fast lookup */
#define UFUNC_ID_ADD        0
#define UFUNC_ID_SUBTRACT   1
#define UFUNC_ID_MULTIPLY   2
#define UFUNC_ID_DIVIDE     3
#define UFUNC_ID_NEGATIVE   4
#define UFUNC_ID_ABSOLUTE   5
#define UFUNC_ID_SQRT       6
#define UFUNC_ID_EXP        7
#define UFUNC_ID_LOG        8
#define UFUNC_ID_SIN        9
#define UFUNC_ID_COS        10
#define UFUNC_ID_MAXIMUM    11
#define UFUNC_ID_MINIMUM    12
#define UFUNC_ID_EQUAL      13
#define UFUNC_ID_LESS       14
#define UFUNC_ID_GREATER    15
/* ... more IDs ... */
#define UFUNC_COUNT         64

/* ============ Ufunc Application ============ */

/**
 * Apply unary ufunc to array.
 *
 * @param ufunc_id  Ufunc identifier
 * @param input     Input array
 * @param output    Output array (NULL to allocate)
 * @param out_dtype Desired output dtype (or input dtype if DTYPE_NONE)
 * @return          Result array (may be output if provided)
 */
NDArray* ufunc_apply_unary(int ufunc_id, NDArray* input,
                           NDArray* output, DType out_dtype);

/**
 * Apply binary ufunc to two arrays with broadcasting.
 *
 * @param ufunc_id  Ufunc identifier
 * @param input1    First input array
 * @param input2    Second input array
 * @param output    Output array (NULL to allocate)
 * @param out_dtype Desired output dtype
 * @return          Result array
 */
NDArray* ufunc_apply_binary(int ufunc_id, NDArray* input1, NDArray* input2,
                            NDArray* output, DType out_dtype);

/* ============ Type Resolution ============ */

/**
 * Resolve output dtype for ufunc given input dtypes.
 * Follows NumPy type promotion rules.
 */
DType ufunc_resolve_dtype(int ufunc_id, DType in_dtype1, DType in_dtype2);

/**
 * Find matching loop for given dtypes.
 * @return Loop index or -1 if no match.
 */
int ufunc_find_loop(const UfuncDef* ufunc, DType in1, DType in2, DType out);

#endif /* NUMJS_UFUNC_H */
```

#### 3.1.2 C: Strided Iteration Helpers

**File:** `src/wasm/ufunc.c` (new file, partial)

```c
#include "ufunc.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#define EXPORT EMSCRIPTEN_KEEPALIVE
#else
#define EXPORT
#endif

/* ============ Strided Iteration ============ */

/**
 * Multi-dimensional iterator for ufunc operations.
 * Handles arbitrary strides and broadcasting.
 */
typedef struct {
    int ndim;
    int64_t shape[32];
    int64_t strides[4][32];  /* Up to 4 arrays (2 in + 1 out + 1 spare) */
    char* data[4];
    int64_t coords[32];
    size_t index;
    size_t size;
} UfuncIterator;

/**
 * Initialize iterator for unary operation.
 */
static void ufunc_iter_init_unary(UfuncIterator* iter,
                                   NDArray* input, NDArray* output)
{
    iter->ndim = output->ndim;
    iter->size = 1;

    for (int i = 0; i < iter->ndim; i++) {
        iter->shape[i] = output->shape[i];
        iter->size *= output->shape[i];
        iter->coords[i] = 0;

        /* Input strides (may be 0 if broadcasting, but unary shouldn't broadcast) */
        iter->strides[0][i] = input->strides[i];
        /* Output strides */
        iter->strides[1][i] = output->strides[i];
    }

    iter->data[0] = (char*)input->data;
    iter->data[1] = (char*)output->data;
    iter->index = 0;
}

/**
 * Initialize iterator for binary operation with broadcasting.
 */
static void ufunc_iter_init_binary(UfuncIterator* iter,
                                    NDArray* in1, NDArray* in2, NDArray* out,
                                    const int32_t* broadcast_shape, int broadcast_ndim)
{
    iter->ndim = broadcast_ndim;
    iter->size = 1;

    int in1_offset = broadcast_ndim - in1->ndim;
    int in2_offset = broadcast_ndim - in2->ndim;

    for (int i = 0; i < broadcast_ndim; i++) {
        iter->shape[i] = broadcast_shape[i];
        iter->size *= broadcast_shape[i];
        iter->coords[i] = 0;

        /* Input 1 strides: 0 for prepended or size-1 dims */
        int in1_idx = i - in1_offset;
        if (in1_idx < 0 || in1->shape[in1_idx] == 1) {
            iter->strides[0][i] = 0;
        } else {
            iter->strides[0][i] = in1->strides[in1_idx];
        }

        /* Input 2 strides */
        int in2_idx = i - in2_offset;
        if (in2_idx < 0 || in2->shape[in2_idx] == 1) {
            iter->strides[1][i] = 0;
        } else {
            iter->strides[1][i] = in2->strides[in2_idx];
        }

        /* Output strides (no broadcasting) */
        iter->strides[2][i] = out->strides[i];
    }

    iter->data[0] = (char*)in1->data;
    iter->data[1] = (char*)in2->data;
    iter->data[2] = (char*)out->data;
    iter->index = 0;
}

/**
 * Advance iterator to next position.
 * Updates data pointers based on strides.
 */
static int ufunc_iter_next(UfuncIterator* iter, int num_arrays)
{
    iter->index++;
    if (iter->index >= iter->size) return 0;

    /* Increment coordinates (like an odometer) */
    for (int i = iter->ndim - 1; i >= 0; i--) {
        iter->coords[i]++;
        if (iter->coords[i] < iter->shape[i]) {
            /* Update data pointers */
            for (int j = 0; j < num_arrays; j++) {
                iter->data[j] += iter->strides[j][i];
            }
            return 1;
        }

        /* Carry over: reset this dim, continue to next */
        iter->coords[i] = 0;
        for (int j = 0; j < num_arrays; j++) {
            iter->data[j] -= iter->strides[j][i] * (iter->shape[i] - 1);
        }
    }

    return 0;
}
```

#### 3.1.3 TS: Ufunc Class

**File:** `src/ts/ufunc.ts` (new file)

```typescript
import { NDArray, zeros, empty } from './NDArray.js';
import { DType } from './types.js';
import { promoteTypes } from './dtype.js';
import { broadcastShapes, broadcastArrays } from './broadcast.js';
import { getWasmModule } from './wasm-loader.js';

/**
 * Ufunc identity values for reduction.
 */
export enum UfuncIdentity {
  Zero = 0,
  One = 1,
  None = -1,
  Reorderable = 2,
}

/**
 * Universal Function class.
 * Provides element-wise operations with broadcasting support.
 */
export class Ufunc {
  /** Ufunc ID for WASM dispatch */
  readonly id: number;

  /** Human-readable name */
  readonly name: string;

  /** Number of inputs (1 for unary, 2 for binary) */
  readonly nin: number;

  /** Number of outputs (always 1) */
  readonly nout: number;

  /** Identity value for reduction (null if none) */
  readonly identity: number | null;

  constructor(
    id: number,
    name: string,
    nin: number,
    nout: number,
    identity: UfuncIdentity
  ) {
    this.id = id;
    this.name = name;
    this.nin = nin;
    this.nout = nout;
    this.identity = identity === UfuncIdentity.None ? null : identity;
  }

  /**
   * Apply ufunc to input array(s).
   *
   * @param args Input arrays
   * @param out Optional output array
   * @param dtype Optional output dtype
   * @returns Result array
   */
  call(args: NDArray[], out?: NDArray, dtype?: DType): NDArray {
    if (args.length !== this.nin) {
      throw new Error(
        `${this.name} requires ${this.nin} input(s), got ${args.length}`
      );
    }

    const module = getWasmModule();

    if (this.nin === 1) {
      return this._applyUnary(args[0], out, dtype);
    } else {
      return this._applyBinary(args[0], args[1], out, dtype);
    }
  }

  /**
   * Reduce array along axis using this ufunc.
   */
  reduce(
    arr: NDArray,
    axis?: number | number[] | null,
    dtype?: DType,
    out?: NDArray,
    keepdims: boolean = false,
    initial?: number
  ): NDArray {
    if (this.nin !== 2) {
      throw new Error(`reduce requires binary ufunc, ${this.name} is unary`);
    }

    // Implementation calls WASM reduction
    return this._reduce(arr, axis, dtype, out, keepdims, initial);
  }

  /**
   * Accumulate operation along axis.
   */
  accumulate(
    arr: NDArray,
    axis: number = 0,
    dtype?: DType,
    out?: NDArray
  ): NDArray {
    if (this.nin !== 2) {
      throw new Error(`accumulate requires binary ufunc`);
    }

    return this._accumulate(arr, axis, dtype, out);
  }

  /**
   * Outer operation: apply to all pairs.
   */
  outer(arr1: NDArray, arr2: NDArray, dtype?: DType): NDArray {
    if (this.nin !== 2) {
      throw new Error(`outer requires binary ufunc`);
    }

    return this._outer(arr1, arr2, dtype);
  }

  /* ============ Private Implementation ============ */

  private _applyUnary(input: NDArray, out?: NDArray, dtype?: DType): NDArray {
    const module = getWasmModule();
    const outDtype = dtype ?? input.dtype;

    // Allocate output if not provided
    const output = out ?? empty(input.shape, outDtype);

    // Validate output shape
    if (!arraysEqual(output.shape, input.shape)) {
      throw new Error('Output shape does not match input shape');
    }

    // Call WASM
    const resultPtr = module._ufunc_apply_unary(
      this.id,
      input._ptr,
      output._ptr,
      outDtype
    );

    if (resultPtr === 0) {
      throw new Error(`${this.name} failed`);
    }

    return output;
  }

  private _applyBinary(
    in1: NDArray,
    in2: NDArray,
    out?: NDArray,
    dtype?: DType
  ): NDArray {
    const module = getWasmModule();

    // Compute broadcast shape
    const broadcastShape = broadcastShapes(in1.shape, in2.shape);

    // Determine output dtype
    const outDtype = dtype ?? promoteTypes(in1.dtype, in2.dtype);

    // Allocate output if not provided
    const output = out ?? empty(broadcastShape, outDtype);

    // Validate output shape
    if (!arraysEqual(output.shape, broadcastShape)) {
      throw new Error(
        `Output shape ${output.shape} does not match broadcast shape ${broadcastShape}`
      );
    }

    // Call WASM
    const resultPtr = module._ufunc_apply_binary(
      this.id,
      in1._ptr,
      in2._ptr,
      output._ptr,
      outDtype
    );

    if (resultPtr === 0) {
      throw new Error(`${this.name} failed`);
    }

    return output;
  }

  private _reduce(
    arr: NDArray,
    axis: number | number[] | null | undefined,
    dtype: DType | undefined,
    out: NDArray | undefined,
    keepdims: boolean,
    initial: number | undefined
  ): NDArray {
    const module = getWasmModule();

    // Normalize axis
    const axes = normalizeAxis(axis, arr.ndim);

    // Compute result shape
    const resultShape = computeReduceShape(arr.shape, axes, keepdims);

    // Determine output dtype
    const outDtype = dtype ?? arr.dtype;

    // Allocate output if not provided
    const output = out ?? empty(resultShape, outDtype);

    // Prepare axis flags array
    const axisFlags = new Int32Array(arr.ndim);
    for (const ax of axes) {
      axisFlags[ax] = 1;
    }

    // Call WASM reduction
    // ... implementation details ...

    return output;
  }

  private _accumulate(
    arr: NDArray,
    axis: number,
    dtype: DType | undefined,
    out: NDArray | undefined
  ): NDArray {
    // Implementation
    throw new Error('Not implemented');
  }

  private _outer(arr1: NDArray, arr2: NDArray, dtype?: DType): NDArray {
    // Implementation: reshape arr1 to [...shape1, 1, 1, ...], arr2 to [1, 1, ..., ...shape2]
    // Then apply binary ufunc
    throw new Error('Not implemented');
  }
}

/* ============ Helper Functions ============ */

function arraysEqual(a: number[], b: number[]): boolean {
  if (a.length !== b.length) return false;
  for (let i = 0; i < a.length; i++) {
    if (a[i] !== b[i]) return false;
  }
  return true;
}

function normalizeAxis(
  axis: number | number[] | null | undefined,
  ndim: number
): number[] {
  if (axis === null || axis === undefined) {
    // Reduce over all axes
    return Array.from({ length: ndim }, (_, i) => i);
  }

  if (typeof axis === 'number') {
    const ax = axis < 0 ? ndim + axis : axis;
    if (ax < 0 || ax >= ndim) {
      throw new Error(`axis ${axis} is out of bounds for array of dimension ${ndim}`);
    }
    return [ax];
  }

  // Array of axes
  return axis.map(ax => {
    const normalized = ax < 0 ? ndim + ax : ax;
    if (normalized < 0 || normalized >= ndim) {
      throw new Error(`axis ${ax} is out of bounds for array of dimension ${ndim}`);
    }
    return normalized;
  });
}

function computeReduceShape(
  shape: number[],
  axes: number[],
  keepdims: boolean
): number[] {
  if (keepdims) {
    return shape.map((dim, i) => axes.includes(i) ? 1 : dim);
  } else {
    return shape.filter((_, i) => !axes.includes(i));
  }
}
```

---

### 3.2 Unary Ufuncs

#### 3.2.1 C: Unary Loop Implementations

**File:** `src/wasm/ufunc_unary.h` (new file)

```c
#ifndef NUMJS_UFUNC_UNARY_H
#define NUMJS_UFUNC_UNARY_H

#include "ufunc.h"

/* ============ Arithmetic ============ */

void ufunc_negative_f64(const char* in, char* out, size_t n,
                        int64_t in_stride, int64_t out_stride);
void ufunc_negative_f32(const char* in, char* out, size_t n,
                        int64_t in_stride, int64_t out_stride);
void ufunc_negative_i32(const char* in, char* out, size_t n,
                        int64_t in_stride, int64_t out_stride);

void ufunc_absolute_f64(const char* in, char* out, size_t n,
                        int64_t in_stride, int64_t out_stride);
void ufunc_absolute_f32(const char* in, char* out, size_t n,
                        int64_t in_stride, int64_t out_stride);
void ufunc_absolute_i32(const char* in, char* out, size_t n,
                        int64_t in_stride, int64_t out_stride);

void ufunc_sign_f64(const char* in, char* out, size_t n,
                    int64_t in_stride, int64_t out_stride);

/* ============ Powers ============ */

void ufunc_sqrt_f64(const char* in, char* out, size_t n,
                    int64_t in_stride, int64_t out_stride);
void ufunc_sqrt_f32(const char* in, char* out, size_t n,
                    int64_t in_stride, int64_t out_stride);

void ufunc_square_f64(const char* in, char* out, size_t n,
                      int64_t in_stride, int64_t out_stride);

void ufunc_cbrt_f64(const char* in, char* out, size_t n,
                    int64_t in_stride, int64_t out_stride);

void ufunc_reciprocal_f64(const char* in, char* out, size_t n,
                          int64_t in_stride, int64_t out_stride);

/* ============ Exponential ============ */

void ufunc_exp_f64(const char* in, char* out, size_t n,
                   int64_t in_stride, int64_t out_stride);
void ufunc_exp2_f64(const char* in, char* out, size_t n,
                    int64_t in_stride, int64_t out_stride);
void ufunc_expm1_f64(const char* in, char* out, size_t n,
                     int64_t in_stride, int64_t out_stride);

/* ============ Logarithmic ============ */

void ufunc_log_f64(const char* in, char* out, size_t n,
                   int64_t in_stride, int64_t out_stride);
void ufunc_log2_f64(const char* in, char* out, size_t n,
                    int64_t in_stride, int64_t out_stride);
void ufunc_log10_f64(const char* in, char* out, size_t n,
                     int64_t in_stride, int64_t out_stride);
void ufunc_log1p_f64(const char* in, char* out, size_t n,
                     int64_t in_stride, int64_t out_stride);

/* ============ Trigonometric ============ */

void ufunc_sin_f64(const char* in, char* out, size_t n,
                   int64_t in_stride, int64_t out_stride);
void ufunc_cos_f64(const char* in, char* out, size_t n,
                   int64_t in_stride, int64_t out_stride);
void ufunc_tan_f64(const char* in, char* out, size_t n,
                   int64_t in_stride, int64_t out_stride);
void ufunc_arcsin_f64(const char* in, char* out, size_t n,
                      int64_t in_stride, int64_t out_stride);
void ufunc_arccos_f64(const char* in, char* out, size_t n,
                      int64_t in_stride, int64_t out_stride);
void ufunc_arctan_f64(const char* in, char* out, size_t n,
                      int64_t in_stride, int64_t out_stride);

/* ============ Hyperbolic ============ */

void ufunc_sinh_f64(const char* in, char* out, size_t n,
                    int64_t in_stride, int64_t out_stride);
void ufunc_cosh_f64(const char* in, char* out, size_t n,
                    int64_t in_stride, int64_t out_stride);
void ufunc_tanh_f64(const char* in, char* out, size_t n,
                    int64_t in_stride, int64_t out_stride);
void ufunc_arcsinh_f64(const char* in, char* out, size_t n,
                       int64_t in_stride, int64_t out_stride);
void ufunc_arccosh_f64(const char* in, char* out, size_t n,
                       int64_t in_stride, int64_t out_stride);
void ufunc_arctanh_f64(const char* in, char* out, size_t n,
                       int64_t in_stride, int64_t out_stride);

/* ============ Rounding ============ */

void ufunc_floor_f64(const char* in, char* out, size_t n,
                     int64_t in_stride, int64_t out_stride);
void ufunc_ceil_f64(const char* in, char* out, size_t n,
                    int64_t in_stride, int64_t out_stride);
void ufunc_trunc_f64(const char* in, char* out, size_t n,
                     int64_t in_stride, int64_t out_stride);
void ufunc_rint_f64(const char* in, char* out, size_t n,
                    int64_t in_stride, int64_t out_stride);

/* ============ Angle Conversion ============ */

void ufunc_degrees_f64(const char* in, char* out, size_t n,
                       int64_t in_stride, int64_t out_stride);
void ufunc_radians_f64(const char* in, char* out, size_t n,
                       int64_t in_stride, int64_t out_stride);

/* ============ Predicates (return bool) ============ */

void ufunc_isnan_f64(const char* in, char* out, size_t n,
                     int64_t in_stride, int64_t out_stride);
void ufunc_isinf_f64(const char* in, char* out, size_t n,
                     int64_t in_stride, int64_t out_stride);
void ufunc_isfinite_f64(const char* in, char* out, size_t n,
                        int64_t in_stride, int64_t out_stride);
void ufunc_signbit_f64(const char* in, char* out, size_t n,
                       int64_t in_stride, int64_t out_stride);

#endif /* NUMJS_UFUNC_UNARY_H */
```

**File:** `src/wasm/ufunc_unary.c` (new file)

```c
#include "ufunc_unary.h"
#include <math.h>
#include <stdint.h>

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#define EXPORT EMSCRIPTEN_KEEPALIVE
#else
#define EXPORT
#endif

/* ============ Macro for generating strided loops ============ */

#define UNARY_LOOP_F64(name, op) \
EXPORT void ufunc_##name##_f64(const char* in, char* out, size_t n, \
                                int64_t in_stride, int64_t out_stride) \
{ \
    for (size_t i = 0; i < n; i++) { \
        double x = *(const double*)(in + i * in_stride); \
        *(double*)(out + i * out_stride) = op; \
    } \
}

#define UNARY_LOOP_F32(name, op) \
EXPORT void ufunc_##name##_f32(const char* in, char* out, size_t n, \
                                int64_t in_stride, int64_t out_stride) \
{ \
    for (size_t i = 0; i < n; i++) { \
        float x = *(const float*)(in + i * in_stride); \
        *(float*)(out + i * out_stride) = op; \
    } \
}

#define UNARY_LOOP_I32(name, op) \
EXPORT void ufunc_##name##_i32(const char* in, char* out, size_t n, \
                                int64_t in_stride, int64_t out_stride) \
{ \
    for (size_t i = 0; i < n; i++) { \
        int32_t x = *(const int32_t*)(in + i * in_stride); \
        *(int32_t*)(out + i * out_stride) = op; \
    } \
}

/* ============ Arithmetic ============ */

UNARY_LOOP_F64(negative, -x)
UNARY_LOOP_F32(negative, -x)
UNARY_LOOP_I32(negative, -x)

UNARY_LOOP_F64(absolute, fabs(x))
UNARY_LOOP_F32(absolute, fabsf(x))
UNARY_LOOP_I32(absolute, (x < 0) ? -x : x)

UNARY_LOOP_F64(sign, (x > 0) ? 1.0 : ((x < 0) ? -1.0 : 0.0))

/* ============ Powers ============ */

UNARY_LOOP_F64(sqrt, sqrt(x))
UNARY_LOOP_F32(sqrt, sqrtf(x))

UNARY_LOOP_F64(square, x * x)
UNARY_LOOP_F64(cbrt, cbrt(x))
UNARY_LOOP_F64(reciprocal, 1.0 / x)

/* ============ Exponential ============ */

UNARY_LOOP_F64(exp, exp(x))
UNARY_LOOP_F64(exp2, exp2(x))
UNARY_LOOP_F64(expm1, expm1(x))

/* ============ Logarithmic ============ */

UNARY_LOOP_F64(log, log(x))
UNARY_LOOP_F64(log2, log2(x))
UNARY_LOOP_F64(log10, log10(x))
UNARY_LOOP_F64(log1p, log1p(x))

/* ============ Trigonometric ============ */

UNARY_LOOP_F64(sin, sin(x))
UNARY_LOOP_F64(cos, cos(x))
UNARY_LOOP_F64(tan, tan(x))
UNARY_LOOP_F64(arcsin, asin(x))
UNARY_LOOP_F64(arccos, acos(x))
UNARY_LOOP_F64(arctan, atan(x))

/* ============ Hyperbolic ============ */

UNARY_LOOP_F64(sinh, sinh(x))
UNARY_LOOP_F64(cosh, cosh(x))
UNARY_LOOP_F64(tanh, tanh(x))
UNARY_LOOP_F64(arcsinh, asinh(x))
UNARY_LOOP_F64(arccosh, acosh(x))
UNARY_LOOP_F64(arctanh, atanh(x))

/* ============ Rounding ============ */

UNARY_LOOP_F64(floor, floor(x))
UNARY_LOOP_F64(ceil, ceil(x))
UNARY_LOOP_F64(trunc, trunc(x))
UNARY_LOOP_F64(rint, rint(x))

/* ============ Angle Conversion ============ */

#define RAD_TO_DEG (180.0 / 3.14159265358979323846)
#define DEG_TO_RAD (3.14159265358979323846 / 180.0)

UNARY_LOOP_F64(degrees, x * RAD_TO_DEG)
UNARY_LOOP_F64(radians, x * DEG_TO_RAD)

/* ============ Predicates ============ */

EXPORT void ufunc_isnan_f64(const char* in, char* out, size_t n,
                             int64_t in_stride, int64_t out_stride)
{
    for (size_t i = 0; i < n; i++) {
        double x = *(const double*)(in + i * in_stride);
        *(uint8_t*)(out + i * out_stride) = isnan(x) ? 1 : 0;
    }
}

EXPORT void ufunc_isinf_f64(const char* in, char* out, size_t n,
                             int64_t in_stride, int64_t out_stride)
{
    for (size_t i = 0; i < n; i++) {
        double x = *(const double*)(in + i * in_stride);
        *(uint8_t*)(out + i * out_stride) = isinf(x) ? 1 : 0;
    }
}

EXPORT void ufunc_isfinite_f64(const char* in, char* out, size_t n,
                                int64_t in_stride, int64_t out_stride)
{
    for (size_t i = 0; i < n; i++) {
        double x = *(const double*)(in + i * in_stride);
        *(uint8_t*)(out + i * out_stride) = isfinite(x) ? 1 : 0;
    }
}

EXPORT void ufunc_signbit_f64(const char* in, char* out, size_t n,
                               int64_t in_stride, int64_t out_stride)
{
    for (size_t i = 0; i < n; i++) {
        double x = *(const double*)(in + i * in_stride);
        *(uint8_t*)(out + i * out_stride) = signbit(x) ? 1 : 0;
    }
}
```

---

### 3.3 Binary Ufuncs

#### 3.3.1 C: Binary Loop Implementations

**File:** `src/wasm/ufunc_binary.h` (new file)

```c
#ifndef NUMJS_UFUNC_BINARY_H
#define NUMJS_UFUNC_BINARY_H

#include "ufunc.h"

/* ============ Arithmetic ============ */

void ufunc_add_f64(const char* in1, const char* in2, char* out, size_t n,
                   int64_t s1, int64_t s2, int64_t so);
void ufunc_add_f32(const char* in1, const char* in2, char* out, size_t n,
                   int64_t s1, int64_t s2, int64_t so);
void ufunc_add_i32(const char* in1, const char* in2, char* out, size_t n,
                   int64_t s1, int64_t s2, int64_t so);

void ufunc_subtract_f64(const char* in1, const char* in2, char* out, size_t n,
                        int64_t s1, int64_t s2, int64_t so);
void ufunc_multiply_f64(const char* in1, const char* in2, char* out, size_t n,
                        int64_t s1, int64_t s2, int64_t so);
void ufunc_divide_f64(const char* in1, const char* in2, char* out, size_t n,
                      int64_t s1, int64_t s2, int64_t so);
void ufunc_floor_divide_f64(const char* in1, const char* in2, char* out, size_t n,
                            int64_t s1, int64_t s2, int64_t so);
void ufunc_remainder_f64(const char* in1, const char* in2, char* out, size_t n,
                         int64_t s1, int64_t s2, int64_t so);
void ufunc_power_f64(const char* in1, const char* in2, char* out, size_t n,
                     int64_t s1, int64_t s2, int64_t so);

/* ============ Comparison ============ */

void ufunc_equal_f64(const char* in1, const char* in2, char* out, size_t n,
                     int64_t s1, int64_t s2, int64_t so);
void ufunc_not_equal_f64(const char* in1, const char* in2, char* out, size_t n,
                         int64_t s1, int64_t s2, int64_t so);
void ufunc_less_f64(const char* in1, const char* in2, char* out, size_t n,
                    int64_t s1, int64_t s2, int64_t so);
void ufunc_less_equal_f64(const char* in1, const char* in2, char* out, size_t n,
                          int64_t s1, int64_t s2, int64_t so);
void ufunc_greater_f64(const char* in1, const char* in2, char* out, size_t n,
                       int64_t s1, int64_t s2, int64_t so);
void ufunc_greater_equal_f64(const char* in1, const char* in2, char* out, size_t n,
                             int64_t s1, int64_t s2, int64_t so);

/* ============ Extrema ============ */

void ufunc_maximum_f64(const char* in1, const char* in2, char* out, size_t n,
                       int64_t s1, int64_t s2, int64_t so);
void ufunc_minimum_f64(const char* in1, const char* in2, char* out, size_t n,
                       int64_t s1, int64_t s2, int64_t so);
void ufunc_fmax_f64(const char* in1, const char* in2, char* out, size_t n,
                    int64_t s1, int64_t s2, int64_t so);
void ufunc_fmin_f64(const char* in1, const char* in2, char* out, size_t n,
                    int64_t s1, int64_t s2, int64_t so);

/* ============ Logical ============ */

void ufunc_logical_and_f64(const char* in1, const char* in2, char* out, size_t n,
                           int64_t s1, int64_t s2, int64_t so);
void ufunc_logical_or_f64(const char* in1, const char* in2, char* out, size_t n,
                          int64_t s1, int64_t s2, int64_t so);
void ufunc_logical_xor_f64(const char* in1, const char* in2, char* out, size_t n,
                           int64_t s1, int64_t s2, int64_t so);

/* ============ Bitwise (integer types) ============ */

void ufunc_bitwise_and_i32(const char* in1, const char* in2, char* out, size_t n,
                           int64_t s1, int64_t s2, int64_t so);
void ufunc_bitwise_or_i32(const char* in1, const char* in2, char* out, size_t n,
                          int64_t s1, int64_t s2, int64_t so);
void ufunc_bitwise_xor_i32(const char* in1, const char* in2, char* out, size_t n,
                           int64_t s1, int64_t s2, int64_t so);
void ufunc_left_shift_i32(const char* in1, const char* in2, char* out, size_t n,
                          int64_t s1, int64_t s2, int64_t so);
void ufunc_right_shift_i32(const char* in1, const char* in2, char* out, size_t n,
                           int64_t s1, int64_t s2, int64_t so);

/* ============ Special ============ */

void ufunc_arctan2_f64(const char* in1, const char* in2, char* out, size_t n,
                       int64_t s1, int64_t s2, int64_t so);
void ufunc_hypot_f64(const char* in1, const char* in2, char* out, size_t n,
                     int64_t s1, int64_t s2, int64_t so);
void ufunc_copysign_f64(const char* in1, const char* in2, char* out, size_t n,
                        int64_t s1, int64_t s2, int64_t so);

#endif /* NUMJS_UFUNC_BINARY_H */
```

**File:** `src/wasm/ufunc_binary.c` (new file)

```c
#include "ufunc_binary.h"
#include <math.h>
#include <stdint.h>

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#define EXPORT EMSCRIPTEN_KEEPALIVE
#else
#define EXPORT
#endif

/* ============ Macro for generating binary strided loops ============ */

#define BINARY_LOOP_F64(name, op) \
EXPORT void ufunc_##name##_f64(const char* in1, const char* in2, char* out, \
                                size_t n, int64_t s1, int64_t s2, int64_t so) \
{ \
    for (size_t i = 0; i < n; i++) { \
        double x = *(const double*)(in1 + i * s1); \
        double y = *(const double*)(in2 + i * s2); \
        *(double*)(out + i * so) = op; \
    } \
}

#define BINARY_LOOP_F32(name, op) \
EXPORT void ufunc_##name##_f32(const char* in1, const char* in2, char* out, \
                                size_t n, int64_t s1, int64_t s2, int64_t so) \
{ \
    for (size_t i = 0; i < n; i++) { \
        float x = *(const float*)(in1 + i * s1); \
        float y = *(const float*)(in2 + i * s2); \
        *(float*)(out + i * so) = op; \
    } \
}

#define BINARY_LOOP_I32(name, op) \
EXPORT void ufunc_##name##_i32(const char* in1, const char* in2, char* out, \
                                size_t n, int64_t s1, int64_t s2, int64_t so) \
{ \
    for (size_t i = 0; i < n; i++) { \
        int32_t x = *(const int32_t*)(in1 + i * s1); \
        int32_t y = *(const int32_t*)(in2 + i * s2); \
        *(int32_t*)(out + i * so) = op; \
    } \
}

/* Comparison returns bool (uint8) */
#define BINARY_CMP_F64(name, op) \
EXPORT void ufunc_##name##_f64(const char* in1, const char* in2, char* out, \
                                size_t n, int64_t s1, int64_t s2, int64_t so) \
{ \
    for (size_t i = 0; i < n; i++) { \
        double x = *(const double*)(in1 + i * s1); \
        double y = *(const double*)(in2 + i * s2); \
        *(uint8_t*)(out + i * so) = (op) ? 1 : 0; \
    } \
}

/* ============ Arithmetic ============ */

BINARY_LOOP_F64(add, x + y)
BINARY_LOOP_F32(add, x + y)
BINARY_LOOP_I32(add, x + y)

BINARY_LOOP_F64(subtract, x - y)
BINARY_LOOP_F32(subtract, x - y)
BINARY_LOOP_I32(subtract, x - y)

BINARY_LOOP_F64(multiply, x * y)
BINARY_LOOP_F32(multiply, x * y)
BINARY_LOOP_I32(multiply, x * y)

BINARY_LOOP_F64(divide, x / y)
BINARY_LOOP_F32(divide, x / y)

BINARY_LOOP_F64(floor_divide, floor(x / y))
BINARY_LOOP_F64(remainder, fmod(x, y))
BINARY_LOOP_F64(power, pow(x, y))

/* ============ Comparison ============ */

BINARY_CMP_F64(equal, x == y)
BINARY_CMP_F64(not_equal, x != y)
BINARY_CMP_F64(less, x < y)
BINARY_CMP_F64(less_equal, x <= y)
BINARY_CMP_F64(greater, x > y)
BINARY_CMP_F64(greater_equal, x >= y)

/* ============ Extrema ============ */

BINARY_LOOP_F64(maximum, (x > y) ? x : y)
BINARY_LOOP_F64(minimum, (x < y) ? x : y)

/* fmax/fmin ignore NaN (return other value if one is NaN) */
BINARY_LOOP_F64(fmax, fmax(x, y))
BINARY_LOOP_F64(fmin, fmin(x, y))

/* ============ Logical ============ */

BINARY_CMP_F64(logical_and, (x != 0) && (y != 0))
BINARY_CMP_F64(logical_or, (x != 0) || (y != 0))
BINARY_CMP_F64(logical_xor, ((x != 0) && (y == 0)) || ((x == 0) && (y != 0)))

/* ============ Bitwise ============ */

BINARY_LOOP_I32(bitwise_and, x & y)
BINARY_LOOP_I32(bitwise_or, x | y)
BINARY_LOOP_I32(bitwise_xor, x ^ y)
BINARY_LOOP_I32(left_shift, x << y)
BINARY_LOOP_I32(right_shift, x >> y)

/* ============ Special ============ */

BINARY_LOOP_F64(arctan2, atan2(x, y))
BINARY_LOOP_F64(hypot, hypot(x, y))
BINARY_LOOP_F64(copysign, copysign(x, y))
```

---

### 3.4 Reduction Operations

#### 3.4.1 C: Reduction Infrastructure

**File:** `src/wasm/reduction.h` (new file)

```c
#ifndef NUMJS_REDUCTION_H
#define NUMJS_REDUCTION_H

#include "ndarray.h"
#include "ufunc.h"

/**
 * Reduction context for axis-based reductions.
 */
typedef struct {
    NDArray* operand;           /* Input array */
    NDArray* result;            /* Output array */

    int32_t* axis_flags;        /* Boolean flags: 1 = reduce this axis */
    int keepdims;               /* Keep reduced dimensions as size 1 */

    DType out_dtype;            /* Output dtype */
    double initial;             /* Initial value (or NaN for none) */
    int has_initial;            /* Whether initial was provided */

    /* Ufunc for reduction */
    int ufunc_id;               /* add, multiply, maximum, minimum, etc. */
} ReduceContext;

/**
 * Perform reduction over specified axes.
 *
 * @param arr         Input array
 * @param ufunc_id    Ufunc to use (add for sum, multiply for prod, etc.)
 * @param axis_flags  Boolean array: 1 for axes to reduce (size = arr->ndim)
 * @param keepdims    Keep reduced dims as size 1
 * @param out_dtype   Output dtype
 * @param initial     Initial value (ignored if has_initial=0)
 * @param has_initial Whether to use initial value
 * @return            Result array or NULL on error
 */
NDArray* ndarray_reduce(NDArray* arr, int ufunc_id,
                        const int32_t* axis_flags, int keepdims,
                        DType out_dtype, double initial, int has_initial);

/**
 * Compute sum over specified axes.
 * Uses pairwise summation for numerical accuracy.
 */
NDArray* ndarray_sum(NDArray* arr, const int32_t* axis_flags, int keepdims,
                     DType out_dtype, double initial, int has_initial);

/**
 * Compute product over specified axes.
 */
NDArray* ndarray_prod(NDArray* arr, const int32_t* axis_flags, int keepdims,
                      DType out_dtype, double initial, int has_initial);

/**
 * Compute maximum over specified axes.
 */
NDArray* ndarray_max(NDArray* arr, const int32_t* axis_flags, int keepdims);

/**
 * Compute minimum over specified axes.
 */
NDArray* ndarray_min(NDArray* arr, const int32_t* axis_flags, int keepdims);

/**
 * Compute mean over specified axes.
 */
NDArray* ndarray_mean(NDArray* arr, const int32_t* axis_flags, int keepdims,
                      DType out_dtype);

/**
 * Find index of maximum value along axis.
 */
NDArray* ndarray_argmax(NDArray* arr, int axis, int keepdims);

/**
 * Find index of minimum value along axis.
 */
NDArray* ndarray_argmin(NDArray* arr, int axis, int keepdims);

/**
 * Test if all elements are true along axis.
 */
NDArray* ndarray_all(NDArray* arr, const int32_t* axis_flags, int keepdims);

/**
 * Test if any element is true along axis.
 */
NDArray* ndarray_any(NDArray* arr, const int32_t* axis_flags, int keepdims);

#endif /* NUMJS_REDUCTION_H */
```

#### 3.4.2 C: Reduction Implementation

**File:** `src/wasm/reduction.c` (new file)

```c
#include "reduction.h"
#include "pairwise_sum.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#define EXPORT EMSCRIPTEN_KEEPALIVE
#else
#define EXPORT
#endif

/* ============ Helper: Compute result shape ============ */

static void compute_result_shape(const int32_t* in_shape, int in_ndim,
                                  const int32_t* axis_flags, int keepdims,
                                  int32_t* out_shape, int* out_ndim)
{
    int j = 0;
    for (int i = 0; i < in_ndim; i++) {
        if (axis_flags[i]) {
            if (keepdims) {
                out_shape[j++] = 1;
            }
            /* else: skip this dimension */
        } else {
            out_shape[j++] = in_shape[i];
        }
    }
    *out_ndim = j;

    /* Handle scalar result */
    if (*out_ndim == 0) {
        *out_ndim = 0;  /* 0-dimensional array */
    }
}

/* ============ Sum Reduction ============ */

/**
 * Internal: sum along innermost reduced axis using pairwise summation.
 */
static void sum_along_axis_f64(const double* data, double* result,
                                const int32_t* shape, const int32_t* strides,
                                int ndim, int axis)
{
    /* Compute iteration dimensions */
    size_t outer_size = 1;
    for (int i = 0; i < axis; i++) {
        outer_size *= shape[i];
    }

    size_t inner_size = 1;
    for (int i = axis + 1; i < ndim; i++) {
        inner_size *= shape[i];
    }

    size_t reduce_size = shape[axis];
    int64_t reduce_stride = strides[axis];

    /* For each position in outer × inner dimensions */
    for (size_t outer = 0; outer < outer_size; outer++) {
        for (size_t inner = 0; inner < inner_size; inner++) {
            /* Compute starting offset */
            size_t offset = 0;
            size_t remaining = outer;
            for (int i = axis - 1; i >= 0; i--) {
                size_t idx = remaining % shape[i];
                remaining /= shape[i];
                offset += idx * strides[i];
            }
            remaining = inner;
            for (int i = ndim - 1; i > axis; i--) {
                size_t idx = remaining % shape[i];
                remaining /= shape[i];
                offset += idx * strides[i];
            }

            /* Use pairwise sum along reduce axis */
            const char* start = (const char*)data + offset;
            double sum = pairwise_sum_strided_f64(start, reduce_size, reduce_stride);

            /* Store result */
            result[outer * inner_size + inner] = sum;
        }
    }
}

EXPORT NDArray* ndarray_sum(NDArray* arr, const int32_t* axis_flags, int keepdims,
                             DType out_dtype, double initial, int has_initial)
{
    if (!arr || !axis_flags) return NULL;

    /* Check if reducing all axes */
    int reduce_all = 1;
    for (int i = 0; i < arr->ndim; i++) {
        if (!axis_flags[i]) {
            reduce_all = 0;
            break;
        }
    }

    if (reduce_all) {
        /* Special case: sum all elements */
        double total;
        if (ndarray_is_c_contiguous(arr)) {
            total = pairwise_sum_f64((const double*)arr->data, arr->size);
        } else {
            /* Non-contiguous: iterate with strides */
            total = 0.0;
            /* ... strided iteration ... */
        }

        if (has_initial) {
            total += initial;
        }

        /* Create scalar result */
        int32_t scalar_shape[] = {};
        NDArray* result = ndarray_empty(0, scalar_shape, out_dtype);
        if (result) {
            *(double*)result->data = total;
        }
        return result;
    }

    /* Compute result shape */
    int32_t result_shape[32];
    int result_ndim;
    compute_result_shape(arr->shape, arr->ndim, axis_flags, keepdims,
                         result_shape, &result_ndim);

    /* Allocate result */
    NDArray* result = ndarray_empty(result_ndim, result_shape, out_dtype);
    if (!result) return NULL;

    /* Initialize with initial value if provided */
    if (has_initial) {
        double* out_data = (double*)result->data;
        for (size_t i = 0; i < result->size; i++) {
            out_data[i] = initial;
        }
    }

    /* Reduce along each flagged axis */
    /* ... implementation ... */

    return result;
}

/* ============ Max Reduction ============ */

EXPORT NDArray* ndarray_max(NDArray* arr, const int32_t* axis_flags, int keepdims)
{
    if (!arr || !axis_flags) return NULL;

    /* Compute result shape */
    int32_t result_shape[32];
    int result_ndim;
    compute_result_shape(arr->shape, arr->ndim, axis_flags, keepdims,
                         result_shape, &result_ndim);

    /* Allocate result */
    NDArray* result = ndarray_empty(result_ndim, result_shape, arr->dtype);
    if (!result) return NULL;

    /* Initialize to -infinity */
    double* out_data = (double*)result->data;
    for (size_t i = 0; i < result->size; i++) {
        out_data[i] = -DBL_MAX;
    }

    /* Iterate and find maximum */
    /* ... implementation ... */

    return result;
}

/* ============ Min Reduction ============ */

EXPORT NDArray* ndarray_min(NDArray* arr, const int32_t* axis_flags, int keepdims)
{
    if (!arr || !axis_flags) return NULL;

    int32_t result_shape[32];
    int result_ndim;
    compute_result_shape(arr->shape, arr->ndim, axis_flags, keepdims,
                         result_shape, &result_ndim);

    NDArray* result = ndarray_empty(result_ndim, result_shape, arr->dtype);
    if (!result) return NULL;

    /* Initialize to +infinity */
    double* out_data = (double*)result->data;
    for (size_t i = 0; i < result->size; i++) {
        out_data[i] = DBL_MAX;
    }

    /* ... implementation ... */

    return result;
}

/* ============ Mean Reduction ============ */

EXPORT NDArray* ndarray_mean(NDArray* arr, const int32_t* axis_flags, int keepdims,
                              DType out_dtype)
{
    /* Mean = sum / count */
    NDArray* sum_result = ndarray_sum(arr, axis_flags, keepdims,
                                       DTYPE_FLOAT64, 0.0, 0);
    if (!sum_result) return NULL;

    /* Compute count of reduced elements */
    size_t count = 1;
    for (int i = 0; i < arr->ndim; i++) {
        if (axis_flags[i]) {
            count *= arr->shape[i];
        }
    }

    /* Divide sum by count */
    double* data = (double*)sum_result->data;
    for (size_t i = 0; i < sum_result->size; i++) {
        data[i] /= (double)count;
    }

    return sum_result;
}

/* ============ Argmax / Argmin ============ */

EXPORT NDArray* ndarray_argmax(NDArray* arr, int axis, int keepdims)
{
    if (!arr) return NULL;

    /* Normalize axis */
    if (axis < 0) axis += arr->ndim;
    if (axis < 0 || axis >= arr->ndim) return NULL;

    /* Compute result shape */
    int32_t result_shape[32];
    int result_ndim = 0;
    for (int i = 0; i < arr->ndim; i++) {
        if (i == axis) {
            if (keepdims) {
                result_shape[result_ndim++] = 1;
            }
        } else {
            result_shape[result_ndim++] = arr->shape[i];
        }
    }

    /* Result is int64 indices */
    NDArray* result = ndarray_empty(result_ndim, result_shape, DTYPE_INT64);
    if (!result) return NULL;

    /* ... implementation ... */

    return result;
}

/* ============ All / Any ============ */

EXPORT NDArray* ndarray_all(NDArray* arr, const int32_t* axis_flags, int keepdims)
{
    if (!arr || !axis_flags) return NULL;

    int32_t result_shape[32];
    int result_ndim;
    compute_result_shape(arr->shape, arr->ndim, axis_flags, keepdims,
                         result_shape, &result_ndim);

    /* Result is boolean */
    NDArray* result = ndarray_empty(result_ndim, result_shape, DTYPE_BOOL);
    if (!result) return NULL;

    /* Initialize to true (1) */
    uint8_t* out_data = (uint8_t*)result->data;
    for (size_t i = 0; i < result->size; i++) {
        out_data[i] = 1;
    }

    /* ... implementation ... */

    return result;
}

EXPORT NDArray* ndarray_any(NDArray* arr, const int32_t* axis_flags, int keepdims)
{
    if (!arr || !axis_flags) return NULL;

    int32_t result_shape[32];
    int result_ndim;
    compute_result_shape(arr->shape, arr->ndim, axis_flags, keepdims,
                         result_shape, &result_ndim);

    /* Result is boolean */
    NDArray* result = ndarray_empty(result_ndim, result_shape, DTYPE_BOOL);
    if (!result) return NULL;

    /* Initialize to false (0) */
    uint8_t* out_data = (uint8_t*)result->data;
    for (size_t i = 0; i < result->size; i++) {
        out_data[i] = 0;
    }

    /* ... implementation ... */

    return result;
}
```

---

### 3.5 Accumulation Operations

#### 3.5.1 C: Accumulation Infrastructure

**File:** `src/wasm/accumulate.h` (new file)

```c
#ifndef NUMJS_ACCUMULATE_H
#define NUMJS_ACCUMULATE_H

#include "ndarray.h"

/**
 * Cumulative sum along axis.
 */
NDArray* ndarray_cumsum(NDArray* arr, int axis, DType out_dtype);

/**
 * Cumulative product along axis.
 */
NDArray* ndarray_cumprod(NDArray* arr, int axis, DType out_dtype);

/**
 * Compute n-th discrete difference along axis.
 */
NDArray* ndarray_diff(NDArray* arr, int n, int axis,
                      NDArray* prepend, NDArray* append);

#endif /* NUMJS_ACCUMULATE_H */
```

---

## File Changes Summary

### New Files to Create

```
src/wasm/
├── ufunc.h            # Ufunc infrastructure declarations
├── ufunc.c            # Ufunc registry and application
├── ufunc_unary.h      # Unary ufunc declarations
├── ufunc_unary.c      # Unary ufunc implementations
├── ufunc_binary.h     # Binary ufunc declarations
├── ufunc_binary.c     # Binary ufunc implementations
├── reduction.h        # Reduction declarations
├── reduction.c        # Reduction implementations
├── accumulate.h       # Accumulation declarations
└── accumulate.c       # Accumulation implementations

src/ts/
├── ufunc.ts           # Ufunc class and helpers
├── ufuncs.ts          # Ufunc registry (all exported ufuncs)
└── math.ts            # Convenience math functions (sin, cos, exp, etc.)
```

### Files to Modify

```
src/wasm/ndarray.h
└── Add forward declarations for ufunc functions

src/ts/types.ts
├── Add UfuncIdentity enum
├── Add ufunc-related WASM function declarations
└── Add reduction function declarations

src/ts/NDArray.ts
├── Import ufuncs
├── Add operator methods: add(), sub(), mul(), div()
└── Add math methods: sum(), mean(), max(), min()

src/ts/index.ts
├── Export Ufunc class
├── Export all ufuncs (add, subtract, multiply, etc.)
├── Export math functions (sin, cos, exp, log, etc.)
└── Export reduction functions (sum, mean, max, min)

scripts/build-wasm.sh
├── Add ufunc.c to compilation
├── Add ufunc_unary.c to compilation
├── Add ufunc_binary.c to compilation
├── Add reduction.c to compilation
├── Add accumulate.c to compilation
└── Add new EXPORTED_FUNCTIONS
```

---

## Build Script Updates

Add to `scripts/build-wasm.sh` EXPORTED_FUNCTIONS:

```bash
# Ufunc infrastructure
"_ufunc_apply_unary",
"_ufunc_apply_binary",
"_ufunc_resolve_dtype",

# Unary ufuncs (direct calls for performance)
"_ufunc_negative_f64",
"_ufunc_absolute_f64",
"_ufunc_sqrt_f64",
"_ufunc_exp_f64",
"_ufunc_log_f64",
"_ufunc_sin_f64",
"_ufunc_cos_f64",
# ... more unary functions ...

# Binary ufuncs
"_ufunc_add_f64",
"_ufunc_subtract_f64",
"_ufunc_multiply_f64",
"_ufunc_divide_f64",
"_ufunc_power_f64",
"_ufunc_maximum_f64",
"_ufunc_minimum_f64",
# ... more binary functions ...

# Reductions
"_ndarray_sum",
"_ndarray_prod",
"_ndarray_max",
"_ndarray_min",
"_ndarray_mean",
"_ndarray_argmax",
"_ndarray_argmin",
"_ndarray_all",
"_ndarray_any",

# Accumulations
"_ndarray_cumsum",
"_ndarray_cumprod",
"_ndarray_diff"
```

Add new source files:

```bash
emcc \
    "$SRC_DIR/ndarray.c" \
    "$SRC_DIR/dtype.c" \
    "$SRC_DIR/pairwise_sum.c" \
    "$SRC_DIR/broadcast.c" \
    "$SRC_DIR/indexing.c" \
    "$SRC_DIR/ufunc.c" \
    "$SRC_DIR/ufunc_unary.c" \
    "$SRC_DIR/ufunc_binary.c" \
    "$SRC_DIR/reduction.c" \
    "$SRC_DIR/accumulate.c" \
    -o "$OUT_DIR/numjs.cjs" \
    ...
```

---

## WasmModule Interface Updates

Add to `src/ts/types.ts`:

```typescript
// Ufunc application
_ufunc_apply_unary(ufuncId: number, inputPtr: number,
                   outputPtr: number, outDtype: number): number;
_ufunc_apply_binary(ufuncId: number, in1Ptr: number, in2Ptr: number,
                    outputPtr: number, outDtype: number): number;

// Type resolution
_ufunc_resolve_dtype(ufuncId: number, inDtype1: number, inDtype2: number): number;

// Reductions
_ndarray_sum(arrPtr: number, axisFlagsPtr: number, keepdims: number,
             outDtype: number, initial: number, hasInitial: number): number;
_ndarray_prod(arrPtr: number, axisFlagsPtr: number, keepdims: number,
              outDtype: number, initial: number, hasInitial: number): number;
_ndarray_max(arrPtr: number, axisFlagsPtr: number, keepdims: number): number;
_ndarray_min(arrPtr: number, axisFlagsPtr: number, keepdims: number): number;
_ndarray_mean(arrPtr: number, axisFlagsPtr: number, keepdims: number,
              outDtype: number): number;
_ndarray_argmax(arrPtr: number, axis: number, keepdims: number): number;
_ndarray_argmin(arrPtr: number, axis: number, keepdims: number): number;
_ndarray_all(arrPtr: number, axisFlagsPtr: number, keepdims: number): number;
_ndarray_any(arrPtr: number, axisFlagsPtr: number, keepdims: number): number;

// Accumulations
_ndarray_cumsum(arrPtr: number, axis: number, outDtype: number): number;
_ndarray_cumprod(arrPtr: number, axis: number, outDtype: number): number;
_ndarray_diff(arrPtr: number, n: number, axis: number,
              prependPtr: number, appendPtr: number): number;
```

---

## TypeScript API Surface

### Ufunc Registry (`src/ts/ufuncs.ts`)

```typescript
import { Ufunc, UfuncIdentity } from './ufunc.js';

/* ============ Unary Ufuncs ============ */

// Arithmetic
export const negative = new Ufunc(4, 'negative', 1, 1, UfuncIdentity.None);
export const positive = new Ufunc(5, 'positive', 1, 1, UfuncIdentity.None);
export const absolute = new Ufunc(6, 'absolute', 1, 1, UfuncIdentity.None);
export const abs = absolute;  // Alias
export const fabs = new Ufunc(7, 'fabs', 1, 1, UfuncIdentity.None);
export const sign = new Ufunc(8, 'sign', 1, 1, UfuncIdentity.None);

// Powers
export const sqrt = new Ufunc(9, 'sqrt', 1, 1, UfuncIdentity.None);
export const square = new Ufunc(10, 'square', 1, 1, UfuncIdentity.None);
export const cbrt = new Ufunc(11, 'cbrt', 1, 1, UfuncIdentity.None);
export const reciprocal = new Ufunc(12, 'reciprocal', 1, 1, UfuncIdentity.None);

// Exponential
export const exp = new Ufunc(13, 'exp', 1, 1, UfuncIdentity.None);
export const exp2 = new Ufunc(14, 'exp2', 1, 1, UfuncIdentity.None);
export const expm1 = new Ufunc(15, 'expm1', 1, 1, UfuncIdentity.None);

// Logarithmic
export const log = new Ufunc(16, 'log', 1, 1, UfuncIdentity.None);
export const log2 = new Ufunc(17, 'log2', 1, 1, UfuncIdentity.None);
export const log10 = new Ufunc(18, 'log10', 1, 1, UfuncIdentity.None);
export const log1p = new Ufunc(19, 'log1p', 1, 1, UfuncIdentity.None);

// Trigonometric
export const sin = new Ufunc(20, 'sin', 1, 1, UfuncIdentity.None);
export const cos = new Ufunc(21, 'cos', 1, 1, UfuncIdentity.None);
export const tan = new Ufunc(22, 'tan', 1, 1, UfuncIdentity.None);
export const arcsin = new Ufunc(23, 'arcsin', 1, 1, UfuncIdentity.None);
export const arccos = new Ufunc(24, 'arccos', 1, 1, UfuncIdentity.None);
export const arctan = new Ufunc(25, 'arctan', 1, 1, UfuncIdentity.None);

// Hyperbolic
export const sinh = new Ufunc(26, 'sinh', 1, 1, UfuncIdentity.None);
export const cosh = new Ufunc(27, 'cosh', 1, 1, UfuncIdentity.None);
export const tanh = new Ufunc(28, 'tanh', 1, 1, UfuncIdentity.None);
export const arcsinh = new Ufunc(29, 'arcsinh', 1, 1, UfuncIdentity.None);
export const arccosh = new Ufunc(30, 'arccosh', 1, 1, UfuncIdentity.None);
export const arctanh = new Ufunc(31, 'arctanh', 1, 1, UfuncIdentity.None);

// Rounding
export const floor = new Ufunc(32, 'floor', 1, 1, UfuncIdentity.None);
export const ceil = new Ufunc(33, 'ceil', 1, 1, UfuncIdentity.None);
export const trunc = new Ufunc(34, 'trunc', 1, 1, UfuncIdentity.None);
export const rint = new Ufunc(35, 'rint', 1, 1, UfuncIdentity.None);

// Angle conversion
export const degrees = new Ufunc(36, 'degrees', 1, 1, UfuncIdentity.None);
export const radians = new Ufunc(37, 'radians', 1, 1, UfuncIdentity.None);
export const deg2rad = radians;  // Alias
export const rad2deg = degrees;  // Alias

// Predicates
export const isnan = new Ufunc(38, 'isnan', 1, 1, UfuncIdentity.None);
export const isinf = new Ufunc(39, 'isinf', 1, 1, UfuncIdentity.None);
export const isfinite = new Ufunc(40, 'isfinite', 1, 1, UfuncIdentity.None);
export const signbit = new Ufunc(41, 'signbit', 1, 1, UfuncIdentity.None);

/* ============ Binary Ufuncs ============ */

// Arithmetic
export const add = new Ufunc(0, 'add', 2, 1, UfuncIdentity.Zero);
export const subtract = new Ufunc(1, 'subtract', 2, 1, UfuncIdentity.None);
export const multiply = new Ufunc(2, 'multiply', 2, 1, UfuncIdentity.One);
export const divide = new Ufunc(3, 'divide', 2, 1, UfuncIdentity.None);
export const true_divide = divide;  // Alias
export const floor_divide = new Ufunc(42, 'floor_divide', 2, 1, UfuncIdentity.None);
export const remainder = new Ufunc(43, 'remainder', 2, 1, UfuncIdentity.None);
export const mod = remainder;  // Alias
export const fmod = new Ufunc(44, 'fmod', 2, 1, UfuncIdentity.None);
export const power = new Ufunc(45, 'power', 2, 1, UfuncIdentity.None);
export const float_power = power;  // Alias (different behavior in NumPy, same here)

// Comparison
export const equal = new Ufunc(46, 'equal', 2, 1, UfuncIdentity.None);
export const not_equal = new Ufunc(47, 'not_equal', 2, 1, UfuncIdentity.None);
export const less = new Ufunc(48, 'less', 2, 1, UfuncIdentity.None);
export const less_equal = new Ufunc(49, 'less_equal', 2, 1, UfuncIdentity.None);
export const greater = new Ufunc(50, 'greater', 2, 1, UfuncIdentity.None);
export const greater_equal = new Ufunc(51, 'greater_equal', 2, 1, UfuncIdentity.None);

// Extrema
export const maximum = new Ufunc(52, 'maximum', 2, 1, UfuncIdentity.Reorderable);
export const minimum = new Ufunc(53, 'minimum', 2, 1, UfuncIdentity.Reorderable);
export const fmax = new Ufunc(54, 'fmax', 2, 1, UfuncIdentity.Reorderable);
export const fmin = new Ufunc(55, 'fmin', 2, 1, UfuncIdentity.Reorderable);

// Logical
export const logical_and = new Ufunc(56, 'logical_and', 2, 1, UfuncIdentity.One);
export const logical_or = new Ufunc(57, 'logical_or', 2, 1, UfuncIdentity.Zero);
export const logical_xor = new Ufunc(58, 'logical_xor', 2, 1, UfuncIdentity.Zero);
export const logical_not = new Ufunc(59, 'logical_not', 1, 1, UfuncIdentity.None);

// Bitwise
export const bitwise_and = new Ufunc(60, 'bitwise_and', 2, 1, UfuncIdentity.None);
export const bitwise_or = new Ufunc(61, 'bitwise_or', 2, 1, UfuncIdentity.Zero);
export const bitwise_xor = new Ufunc(62, 'bitwise_xor', 2, 1, UfuncIdentity.Zero);
export const invert = new Ufunc(63, 'invert', 1, 1, UfuncIdentity.None);
export const left_shift = new Ufunc(64, 'left_shift', 2, 1, UfuncIdentity.None);
export const right_shift = new Ufunc(65, 'right_shift', 2, 1, UfuncIdentity.None);

// Special
export const arctan2 = new Ufunc(66, 'arctan2', 2, 1, UfuncIdentity.None);
export const hypot = new Ufunc(67, 'hypot', 2, 1, UfuncIdentity.Zero);
export const copysign = new Ufunc(68, 'copysign', 2, 1, UfuncIdentity.None);
```

### Convenience Functions (`src/ts/math.ts`)

```typescript
import { NDArray } from './NDArray.js';
import * as ufuncs from './ufuncs.js';

/**
 * Convenience wrapper: apply ufunc to single array.
 */
function applyUnary(ufunc: Ufunc, arr: NDArray): NDArray {
  return ufunc.call([arr]);
}

/**
 * Convenience wrapper: apply ufunc to two arrays.
 */
function applyBinary(ufunc: Ufunc, arr1: NDArray, arr2: NDArray): NDArray {
  return ufunc.call([arr1, arr2]);
}

/* ============ Exported Functions ============ */

// Unary
export const negative = (arr: NDArray) => applyUnary(ufuncs.negative, arr);
export const abs = (arr: NDArray) => applyUnary(ufuncs.absolute, arr);
export const sqrt = (arr: NDArray) => applyUnary(ufuncs.sqrt, arr);
export const exp = (arr: NDArray) => applyUnary(ufuncs.exp, arr);
export const log = (arr: NDArray) => applyUnary(ufuncs.log, arr);
export const sin = (arr: NDArray) => applyUnary(ufuncs.sin, arr);
export const cos = (arr: NDArray) => applyUnary(ufuncs.cos, arr);
export const tan = (arr: NDArray) => applyUnary(ufuncs.tan, arr);
// ... more unary functions ...

// Binary
export const add = (a: NDArray, b: NDArray) => applyBinary(ufuncs.add, a, b);
export const subtract = (a: NDArray, b: NDArray) => applyBinary(ufuncs.subtract, a, b);
export const multiply = (a: NDArray, b: NDArray) => applyBinary(ufuncs.multiply, a, b);
export const divide = (a: NDArray, b: NDArray) => applyBinary(ufuncs.divide, a, b);
export const power = (a: NDArray, b: NDArray) => applyBinary(ufuncs.power, a, b);
export const maximum = (a: NDArray, b: NDArray) => applyBinary(ufuncs.maximum, a, b);
export const minimum = (a: NDArray, b: NDArray) => applyBinary(ufuncs.minimum, a, b);
// ... more binary functions ...
```

---

## Implementation Order

```
Week 1: Ufunc Infrastructure
├── Day 1: C: UfuncDef struct, registry skeleton
├── Day 2: C: Strided iteration helpers
├── Day 3: C: Type resolution system
├── Day 4: TS: Ufunc class with call() method
└── Day 5: TS: WasmModule interface, basic integration tests

Week 2: Unary Ufuncs
├── Day 1: C: Arithmetic (negative, absolute, sign)
├── Day 2: C: Powers and exponentials (sqrt, exp, log)
├── Day 3: C: Trigonometric (sin, cos, tan, arc*)
├── Day 4: C: Rounding and angle conversion
├── Day 5: TS: Unary wrappers + tests

Week 3: Binary Ufuncs
├── Day 1: C: Arithmetic (add, subtract, multiply, divide)
├── Day 2: C: Comparison (equal, less, greater, etc.)
├── Day 3: C: Binary with broadcasting integration
├── Day 4: C: Logical and bitwise operations
├── Day 5: TS: Binary wrappers + comprehensive tests

Week 4: Reductions
├── Day 1: C: Reduction infrastructure (axis iteration)
├── Day 2: C: sum, prod with pairwise algorithm
├── Day 3: C: max, min, argmax, argmin
├── Day 4: C: mean, all, any
├── Day 5: TS: Reduction wrappers + tests

Week 5: Accumulations & Polish
├── Day 1: C: cumsum, cumprod
├── Day 2: C: diff
├── Day 3: TS: Ufunc methods (reduce, accumulate, outer)
├── Day 4: Edge case handling and error messages
└── Day 5: Documentation, benchmarks, final cleanup
```

---

## Verification Plan

After Level 3 completion, verify:

```bash
# Build
npm run build

# Run all tests
npm test

# Specific Level 3 tests should pass:

# Unary Ufuncs
✓ negative(arr) returns element-wise negation
✓ abs(arr) returns absolute values
✓ sqrt(arr) returns square roots
✓ exp(arr) returns exponentials
✓ log(arr) returns natural logarithms
✓ sin/cos/tan work correctly
✓ floor/ceil/trunc/rint round correctly
✓ isnan/isinf/isfinite return boolean arrays

# Binary Ufuncs
✓ add(a, b) broadcasts correctly
✓ subtract(a, b) broadcasts correctly
✓ multiply(a, b) broadcasts correctly
✓ divide(a, b) broadcasts correctly
✓ power(a, b) computes element-wise power
✓ comparison ufuncs return boolean arrays
✓ maximum/minimum respect NaN handling
✓ logical operations work on boolean arrays
✓ bitwise operations work on integer arrays

# Broadcasting
✓ (5,) + (1,) → (5,)
✓ (3,1) + (1,4) → (3,4)
✓ (2,3,4) + (4,) → (2,3,4)
✓ incompatible shapes raise error

# Reductions
✓ sum(arr) returns total
✓ sum(arr, axis=0) reduces along axis 0
✓ sum(arr, axis=[0,2]) reduces multiple axes
✓ sum(arr, keepdims=true) preserves dimensions
✓ sum with initial value works
✓ prod(arr) returns product
✓ max(arr, axis) finds maximum along axis
✓ min(arr, axis) finds minimum along axis
✓ mean(arr, axis) computes mean
✓ argmax/argmin return indices
✓ all/any work with boolean logic

# Accumulations
✓ cumsum(arr) computes running sum
✓ cumprod(arr) computes running product
✓ diff(arr, n) computes n-th difference

# Type Promotion
✓ int32 + float64 → float64
✓ float32 + float64 → float64
✓ bool + int32 → int32
✓ divide always returns float
```

Generate NumPy comparison vectors:

```python
# tests/python/generate_level3_tests.py
import numpy as np
import json

tests = {
    "unary": [
        {"name": "negative", "input": [1.0, -2.0, 3.0], "expected": [-1.0, 2.0, -3.0]},
        {"name": "sqrt", "input": [1.0, 4.0, 9.0], "expected": [1.0, 2.0, 3.0]},
        {"name": "exp", "input": [0.0, 1.0], "expected": [1.0, 2.718281828]},
        {"name": "sin", "input": [0.0, 1.5707963], "expected": [0.0, 1.0]},
    ],
    "binary_broadcast": [
        {"op": "add", "a_shape": [3, 1], "b_shape": [1, 4], "expected_shape": [3, 4]},
        {"op": "multiply", "a_shape": [5], "b_shape": [1], "expected_shape": [5]},
    ],
    "reductions": [
        {"op": "sum", "shape": [2, 3], "axis": None, "expected_shape": []},
        {"op": "sum", "shape": [2, 3], "axis": 0, "expected_shape": [3]},
        {"op": "sum", "shape": [2, 3], "axis": 1, "expected_shape": [2]},
        {"op": "sum", "shape": [2, 3], "axis": 0, "keepdims": True, "expected_shape": [1, 3]},
        {"op": "max", "shape": [3, 4], "axis": 1, "expected_shape": [3]},
    ],
}

with open("tests/fixtures/level3_vectors.json", "w") as f:
    json.dump(tests, f, indent=2)
```

---

## Critical Dependencies for Level 4

Level 3 completion enables:

- **Level 4.1 Array Manipulation**: Uses ufuncs for concatenate, stack operations
- **Level 4.2 Sorting**: Comparison ufuncs for sort algorithms
- **Level 4.3 Statistics**: Built on reductions (mean, std, var require sum)
- **Level 5 numpy.linalg**: Matrix operations use binary ufuncs

Level 3 MUST be complete before proceeding to Level 4.

---

## Performance Considerations

### 1. Contiguous Fast Path

For C-contiguous arrays, use specialized loops without stride computation:

```c
if (ndarray_is_c_contiguous(arr)) {
    /* Direct memory access, no stride calculation */
    for (size_t i = 0; i < arr->size; i++) {
        out[i] = op(in[i]);
    }
} else {
    /* Strided access */
    for (size_t i = 0; i < arr->size; i++) {
        /* ... compute offset with strides ... */
    }
}
```

### 2. Loop Unrolling

For inner loops, consider 4-way or 8-way unrolling:

```c
/* 4-way unrolled */
size_t i = 0;
for (; i + 4 <= n; i += 4) {
    out[i] = op(in[i]);
    out[i+1] = op(in[i+1]);
    out[i+2] = op(in[i+2]);
    out[i+3] = op(in[i+3]);
}
for (; i < n; i++) {
    out[i] = op(in[i]);
}
```

### 3. SIMD Hints

While WASM SIMD is available, it requires explicit vectorization. For initial implementation, rely on the compiler's auto-vectorization and consider SIMD optimization as a future enhancement.

### 4. Memory Allocation

Minimize allocations by:
- Reusing output arrays when possible (`out` parameter)
- Pre-allocating iteration buffers
- Using stack allocation for small temporary arrays

---

## Estimated Code Volume

| Component | Lines of Code | Complexity |
|-----------|--------------|-----------|
| Ufunc infrastructure (C) | 600-800 | High |
| Unary ufuncs (C) | 400-500 | Low |
| Binary ufuncs (C) | 500-600 | Medium |
| Reductions (C) | 800-1000 | High |
| Accumulations (C) | 300-400 | Medium |
| TypeScript wrapper | 1000-1200 | Medium |
| Ufunc registry | 300-400 | Low |
| Tests | 2000+ | Low |
| **Total** | **5900-6900** | - |
