# Level 0: Foundation Layer Implementation Plan

The foundation layer must be complete before any higher-level features can be implemented. This layer establishes the core data structures, memory management, and type system.

---

## Current State (What Exists)

```
src/wasm/
├── ndarray.h          # NDArray struct, DType enum (4 types)
├── ndarray.c          # create, from_data, free, sum, fill, accessors
├── pairwise_sum.h     # Pairwise summation declarations
└── pairwise_sum.c     # Pairwise sum for float32/float64

src/ts/
├── types.ts           # DType enum, WasmModule interface
├── NDArray.ts         # zeros, ones, fromArray, arange, sum, fill, toArray
├── wasm-loader.ts     # WASM module loading
└── index.ts           # Public exports
```

**Existing DTypes:** Float32, Float64, Int32, Int64
**Existing Operations:** zeros, ones, fromArray, arange, sum, fill, toArray, dispose

---

## Level 0 Implementation Tree

```
LEVEL 0: FOUNDATION
│
├── 0.1 Extended DType System
│   │
│   ├── 0.1.1 C: New DType Enum Values
│   │   ├── DTYPE_BOOL = 4
│   │   ├── DTYPE_INT8 = 5
│   │   ├── DTYPE_INT16 = 6
│   │   ├── DTYPE_UINT8 = 7
│   │   ├── DTYPE_UINT16 = 8
│   │   ├── DTYPE_UINT32 = 9
│   │   ├── DTYPE_UINT64 = 10
│   │   ├── DTYPE_FLOAT16 = 11 (optional, complex to implement)
│   │   ├── DTYPE_COMPLEX64 = 12
│   │   └── DTYPE_COMPLEX128 = 13
│   │
│   ├── 0.1.2 C: dtype_size() Extension
│   │   └── Update switch statement for all new types
│   │
│   ├── 0.1.3 C: DType Info Structure
│   │   ├── dtype_info struct { size, alignment, name, kind }
│   │   ├── dtype_get_info(dtype) → DTypeInfo*
│   │   ├── dtype_is_integer(dtype) → bool
│   │   ├── dtype_is_floating(dtype) → bool
│   │   ├── dtype_is_complex(dtype) → bool
│   │   ├── dtype_is_signed(dtype) → bool
│   │   └── dtype_is_bool(dtype) → bool
│   │
│   ├── 0.1.4 C: Type Promotion Rules
│   │   ├── dtype_promote(dtype1, dtype2) → DType
│   │   ├── dtype_common_type(dtypes[], count) → DType
│   │   ├── dtype_can_cast(from, to, casting) → bool
│   │   └── CastingKind enum { no, equiv, safe, same_kind, unsafe }
│   │
│   ├── 0.1.5 TypeScript: Extended DType Enum
│   │   └── types.ts: Add all new DType values
│   │
│   ├── 0.1.6 TypeScript: DType Utilities
│   │   ├── DTYPE_SIZES: Record<DType, number>
│   │   ├── DTYPE_NAMES: Record<DType, string>
│   │   ├── dtypeFromString(name: string) → DType
│   │   ├── isIntegerDType(dtype: DType) → boolean
│   │   ├── isFloatDType(dtype: DType) → boolean
│   │   └── isComplexDType(dtype: DType) → boolean
│   │
│   └── 0.1.7 TypeScript: TypedArray Mapping
│       ├── TypedArrayType union (extend with Uint8Array, Int8Array, etc.)
│       ├── dtypeToTypedArray(dtype: DType) → TypedArrayConstructor
│       └── typedArrayToDType(arr: TypedArray) → DType
│
├── 0.2 Element Access
│   │
│   ├── 0.2.1 C: Flat Index Calculation
│   │   ├── ndarray_flat_index(arr, indices[], ndim) → size_t
│   │   │   └── Computes: sum(indices[i] * strides[i]) / element_size
│   │   └── ndarray_check_bounds(arr, indices[], ndim) → bool
│   │
│   ├── 0.2.2 C: Element Getters (by flat index)
│   │   ├── ndarray_get_float64(arr, flat_idx) → double
│   │   ├── ndarray_get_float32(arr, flat_idx) → float
│   │   ├── ndarray_get_int64(arr, flat_idx) → int64_t
│   │   ├── ndarray_get_int32(arr, flat_idx) → int32_t
│   │   ├── ndarray_get_int16(arr, flat_idx) → int16_t
│   │   ├── ndarray_get_int8(arr, flat_idx) → int8_t
│   │   ├── ndarray_get_uint64(arr, flat_idx) → uint64_t
│   │   ├── ndarray_get_uint32(arr, flat_idx) → uint32_t
│   │   ├── ndarray_get_uint16(arr, flat_idx) → uint16_t
│   │   ├── ndarray_get_uint8(arr, flat_idx) → uint8_t
│   │   ├── ndarray_get_bool(arr, flat_idx) → bool
│   │   ├── ndarray_get_complex64_real(arr, flat_idx) → float
│   │   ├── ndarray_get_complex64_imag(arr, flat_idx) → float
│   │   ├── ndarray_get_complex128_real(arr, flat_idx) → double
│   │   └── ndarray_get_complex128_imag(arr, flat_idx) → double
│   │
│   ├── 0.2.3 C: Element Setters (by flat index)
│   │   ├── ndarray_set_float64(arr, flat_idx, value)
│   │   ├── ndarray_set_float32(arr, flat_idx, value)
│   │   ├── ndarray_set_int64(arr, flat_idx, value)
│   │   ├── ndarray_set_int32(arr, flat_idx, value)
│   │   ├── ndarray_set_int16(arr, flat_idx, value)
│   │   ├── ndarray_set_int8(arr, flat_idx, value)
│   │   ├── ndarray_set_uint64(arr, flat_idx, value)
│   │   ├── ndarray_set_uint32(arr, flat_idx, value)
│   │   ├── ndarray_set_uint16(arr, flat_idx, value)
│   │   ├── ndarray_set_uint8(arr, flat_idx, value)
│   │   ├── ndarray_set_bool(arr, flat_idx, value)
│   │   ├── ndarray_set_complex64(arr, flat_idx, real, imag)
│   │   └── ndarray_set_complex128(arr, flat_idx, real, imag)
│   │
│   ├── 0.2.4 C: Generic Element Access (multi-index)
│   │   ├── ndarray_get_item(arr, indices[], ndim) → double
│   │   │   └── Returns value converted to double
│   │   ├── ndarray_set_item(arr, indices[], ndim, value)
│   │   │   └── Converts double to array's dtype
│   │   ├── ndarray_get_item_complex(arr, indices[], ndim, *real, *imag)
│   │   └── ndarray_set_item_complex(arr, indices[], ndim, real, imag)
│   │
│   ├── 0.2.5 TypeScript: Element Access Methods
│   │   ├── get(...indices: number[]) → number
│   │   ├── set(value: number, ...indices: number[])
│   │   ├── getComplex(...indices: number[]) → { real: number, imag: number }
│   │   ├── setComplex(real: number, imag: number, ...indices: number[])
│   │   ├── item() → number  (for 0-d or single-element arrays)
│   │   └── itemset(value: number)
│   │
│   └── 0.2.6 TypeScript: Flat Access
│       ├── flat property → FlatIterator (deferred to Level 1)
│       ├── getFlatIndex(index: number) → number
│       └── setFlatIndex(index: number, value: number)
│
├── 0.3 Array Creation (extend existing)
│   │
│   ├── 0.3.1 C: empty() - Uninitialized Array
│   │   └── ndarray_empty(ndim, shape, dtype) → NDArray*
│   │       └── Like create() but skips memset (faster)
│   │
│   ├── 0.3.2 C: full() - Fill with Value
│   │   └── ndarray_full(ndim, shape, dtype, value) → NDArray*
│   │
│   ├── 0.3.3 C: Scalar Array
│   │   └── ndarray_scalar(value, dtype) → NDArray*
│   │       └── Creates 0-dimensional array
│   │
│   ├── 0.3.4 TypeScript: empty()
│   │   └── static async empty(shape, options?) → NDArray
│   │
│   ├── 0.3.5 TypeScript: full()
│   │   └── static async full(shape, fill_value, options?) → NDArray
│   │
│   ├── 0.3.6 TypeScript: zeros_like, ones_like, empty_like, full_like
│   │   ├── static async zerosLike(arr, options?) → NDArray
│   │   ├── static async onesLike(arr, options?) → NDArray
│   │   ├── static async emptyLike(arr, options?) → NDArray
│   │   └── static async fullLike(arr, fill_value, options?) → NDArray
│   │
│   ├── 0.3.7 TypeScript: linspace
│   │   └── static async linspace(start, stop, num?, endpoint?, options?) → NDArray
│   │
│   ├── 0.3.8 TypeScript: logspace
│   │   └── static async logspace(start, stop, num?, endpoint?, base?, options?) → NDArray
│   │
│   ├── 0.3.9 TypeScript: geomspace
│   │   └── static async geomspace(start, stop, num?, endpoint?, options?) → NDArray
│   │
│   ├── 0.3.10 TypeScript: eye, identity
│   │   ├── static async eye(N, M?, k?, options?) → NDArray
│   │   └── static async identity(n, options?) → NDArray
│   │
│   ├── 0.3.11 TypeScript: diag
│   │   └── static async diag(v, k?) → NDArray
│   │       └── Extract diagonal or construct diagonal array
│   │
│   └── 0.3.12 TypeScript: tri, tril, triu
│       ├── static async tri(N, M?, k?, options?) → NDArray
│       ├── static async tril(arr, k?) → NDArray
│       └── static async triu(arr, k?) → NDArray
│
├── 0.4 Memory & Flags
│   │
│   ├── 0.4.1 C: Extended Flags
│   │   ├── NDARRAY_C_CONTIGUOUS = 0x0004
│   │   ├── NDARRAY_F_CONTIGUOUS = 0x0008
│   │   ├── NDARRAY_ALIGNED = 0x0010
│   │   └── NDARRAY_UPDATEIFCOPY = 0x0020 (for future use)
│   │
│   ├── 0.4.2 C: Contiguity Checks
│   │   ├── ndarray_is_c_contiguous(arr) → bool
│   │   ├── ndarray_is_f_contiguous(arr) → bool
│   │   └── ndarray_update_contiguity_flags(arr)
│   │
│   ├── 0.4.3 C: Base Pointer (for views)
│   │   ├── Add: NDArray* base to struct
│   │   ├── ndarray_get_base(arr) → NDArray*
│   │   └── ndarray_set_base(arr, base)
│   │
│   ├── 0.4.4 TypeScript: Flags Property
│   │   └── get flags() → { c_contiguous, f_contiguous, writeable, owndata, aligned }
│   │
│   └── 0.4.5 TypeScript: Memory Info
│       ├── get nbytes() → number  (total bytes)
│       ├── get itemsize() → number  (bytes per element)
│       └── get strides() → number[]
│
└── 0.5 Copy & Conversion
    │
    ├── 0.5.1 C: Array Copy
    │   ├── ndarray_copy(arr) → NDArray*
    │   │   └── Deep copy, always owns data
    │   └── ndarray_copy_to(src, dst)
    │       └── Copy data from src to dst (must be same size)
    │
    ├── 0.5.2 C: Type Casting
    │   └── ndarray_astype(arr, dtype) → NDArray*
    │       └── Returns new array with converted dtype
    │
    ├── 0.5.3 TypeScript: copy()
    │   └── copy() → Promise<NDArray>
    │
    ├── 0.5.4 TypeScript: astype()
    │   └── astype(dtype: DType) → Promise<NDArray>
    │
    └── 0.5.5 TypeScript: toTypedArray()
        └── toTypedArray() → TypedArray
            └── Returns appropriate TypedArray view/copy
```

---

## File Changes Summary

### New Files to Create

```
src/wasm/
├── dtype.h            # DType info, promotion rules, casting
└── dtype.c            # DType utility implementations

src/ts/
└── dtype.ts           # DType utilities, type guards, conversions
```

### Files to Modify

```
src/wasm/ndarray.h
├── Add new DType enum values (DTYPE_BOOL through DTYPE_COMPLEX128)
├── Add base pointer to NDArray struct
├── Add new flag constants
├── Declare element access functions
├── Declare empty, full, scalar creation functions
└── Declare copy, astype functions

src/wasm/ndarray.c
├── Extend dtype_size() for new types
├── Implement element access functions
├── Implement empty, full, scalar creation
├── Implement copy, astype functions
└── Update contiguity flag handling

src/ts/types.ts
├── Extend DType enum with new values
├── Add DTYPE_NAMES mapping
├── Extend TypedArrayType union
└── Add Complex type { real: number, imag: number }

src/ts/NDArray.ts
├── Add get/set methods for element access
├── Add empty, full, *Like factory methods
├── Add linspace, logspace, geomspace
├── Add eye, identity, diag, tri, tril, triu
├── Add copy, astype methods
├── Add flags, nbytes, itemsize, strides properties
└── Add toTypedArray method
```

---

## Detailed Implementation Specifications

### 0.1.1 Extended DType Enum (C)

```c
// In ndarray.h
typedef enum {
    // Existing
    DTYPE_FLOAT32 = 0,
    DTYPE_FLOAT64 = 1,
    DTYPE_INT32 = 2,
    DTYPE_INT64 = 3,
    // New
    DTYPE_BOOL = 4,
    DTYPE_INT8 = 5,
    DTYPE_INT16 = 6,
    DTYPE_UINT8 = 7,
    DTYPE_UINT16 = 8,
    DTYPE_UINT32 = 9,
    DTYPE_UINT64 = 10,
    DTYPE_FLOAT16 = 11,     // Optional: requires half-precision support
    DTYPE_COMPLEX64 = 12,   // float real + float imag
    DTYPE_COMPLEX128 = 13,  // double real + double imag
    DTYPE_COUNT = 14        // Total count for iteration
} DType;
```

### 0.1.2 Extended dtype_size() (C)

```c
// In ndarray.c or dtype.c
EXPORT size_t dtype_size(DType dtype) {
    switch (dtype) {
        case DTYPE_BOOL:       return sizeof(uint8_t);  // 1 byte
        case DTYPE_INT8:       return sizeof(int8_t);   // 1 byte
        case DTYPE_INT16:      return sizeof(int16_t);  // 2 bytes
        case DTYPE_INT32:      return sizeof(int32_t);  // 4 bytes
        case DTYPE_INT64:      return sizeof(int64_t);  // 8 bytes
        case DTYPE_UINT8:      return sizeof(uint8_t);  // 1 byte
        case DTYPE_UINT16:     return sizeof(uint16_t); // 2 bytes
        case DTYPE_UINT32:     return sizeof(uint32_t); // 4 bytes
        case DTYPE_UINT64:     return sizeof(uint64_t); // 8 bytes
        case DTYPE_FLOAT16:    return 2;                // 2 bytes (half)
        case DTYPE_FLOAT32:    return sizeof(float);    // 4 bytes
        case DTYPE_FLOAT64:    return sizeof(double);   // 8 bytes
        case DTYPE_COMPLEX64:  return 2 * sizeof(float);  // 8 bytes
        case DTYPE_COMPLEX128: return 2 * sizeof(double); // 16 bytes
        default: return 0;
    }
}
```

### 0.1.3 DType Info Structure (C)

```c
// In dtype.h
typedef enum {
    DTYPE_KIND_BOOL = 'b',
    DTYPE_KIND_INT = 'i',
    DTYPE_KIND_UINT = 'u',
    DTYPE_KIND_FLOAT = 'f',
    DTYPE_KIND_COMPLEX = 'c'
} DTypeKind;

typedef struct {
    size_t size;        // Size in bytes
    size_t alignment;   // Alignment requirement
    const char* name;   // String name (e.g., "float64")
    DTypeKind kind;     // Category
    bool is_signed;     // For integers
} DTypeInfo;

// Get info for a dtype
const DTypeInfo* dtype_get_info(DType dtype);

// Type predicates
bool dtype_is_integer(DType dtype);
bool dtype_is_floating(DType dtype);
bool dtype_is_complex(DType dtype);
bool dtype_is_signed(DType dtype);
bool dtype_is_bool(DType dtype);
bool dtype_is_numeric(DType dtype);
```

```c
// In dtype.c
static const DTypeInfo DTYPE_INFO[] = {
    [DTYPE_FLOAT32]    = { 4,  4, "float32",    DTYPE_KIND_FLOAT,   true  },
    [DTYPE_FLOAT64]    = { 8,  8, "float64",    DTYPE_KIND_FLOAT,   true  },
    [DTYPE_INT32]      = { 4,  4, "int32",      DTYPE_KIND_INT,     true  },
    [DTYPE_INT64]      = { 8,  8, "int64",      DTYPE_KIND_INT,     true  },
    [DTYPE_BOOL]       = { 1,  1, "bool",       DTYPE_KIND_BOOL,    false },
    [DTYPE_INT8]       = { 1,  1, "int8",       DTYPE_KIND_INT,     true  },
    [DTYPE_INT16]      = { 2,  2, "int16",      DTYPE_KIND_INT,     true  },
    [DTYPE_UINT8]      = { 1,  1, "uint8",      DTYPE_KIND_UINT,    false },
    [DTYPE_UINT16]     = { 2,  2, "uint16",     DTYPE_KIND_UINT,    false },
    [DTYPE_UINT32]     = { 4,  4, "uint32",     DTYPE_KIND_UINT,    false },
    [DTYPE_UINT64]     = { 8,  8, "uint64",     DTYPE_KIND_UINT,    false },
    [DTYPE_FLOAT16]    = { 2,  2, "float16",    DTYPE_KIND_FLOAT,   true  },
    [DTYPE_COMPLEX64]  = { 8,  4, "complex64",  DTYPE_KIND_COMPLEX, true  },
    [DTYPE_COMPLEX128] = { 16, 8, "complex128", DTYPE_KIND_COMPLEX, true  },
};

EXPORT const DTypeInfo* dtype_get_info(DType dtype) {
    if (dtype < 0 || dtype >= DTYPE_COUNT) return NULL;
    return &DTYPE_INFO[dtype];
}

EXPORT bool dtype_is_integer(DType dtype) {
    const DTypeInfo* info = dtype_get_info(dtype);
    return info && (info->kind == DTYPE_KIND_INT || info->kind == DTYPE_KIND_UINT);
}

EXPORT bool dtype_is_floating(DType dtype) {
    const DTypeInfo* info = dtype_get_info(dtype);
    return info && info->kind == DTYPE_KIND_FLOAT;
}

EXPORT bool dtype_is_complex(DType dtype) {
    const DTypeInfo* info = dtype_get_info(dtype);
    return info && info->kind == DTYPE_KIND_COMPLEX;
}
```

### 0.1.4 Type Promotion Rules (C)

```c
// In dtype.h
typedef enum {
    CASTING_NO = 0,        // No casting allowed
    CASTING_EQUIV = 1,     // Only byte-order changes
    CASTING_SAFE = 2,      // Only casts preserving values
    CASTING_SAME_KIND = 3, // Only safe casts or within a kind
    CASTING_UNSAFE = 4     // Any data conversion
} CastingKind;

// Get result type when combining two dtypes
DType dtype_promote(DType dtype1, DType dtype2);

// Get common type for array of dtypes
DType dtype_common_type(DType* dtypes, size_t count);

// Check if cast is allowed
bool dtype_can_cast(DType from, DType to, CastingKind casting);
```

```c
// In dtype.c
// Simplified promotion table (NumPy compatible)
// Higher priority wins, complex > float > int > uint > bool
static const int DTYPE_PRIORITY[] = {
    [DTYPE_BOOL]       = 0,
    [DTYPE_UINT8]      = 1,
    [DTYPE_UINT16]     = 2,
    [DTYPE_UINT32]     = 3,
    [DTYPE_UINT64]     = 4,
    [DTYPE_INT8]       = 5,
    [DTYPE_INT16]      = 6,
    [DTYPE_INT32]      = 7,
    [DTYPE_INT64]      = 8,
    [DTYPE_FLOAT16]    = 9,
    [DTYPE_FLOAT32]    = 10,
    [DTYPE_FLOAT64]    = 11,
    [DTYPE_COMPLEX64]  = 12,
    [DTYPE_COMPLEX128] = 13,
};

EXPORT DType dtype_promote(DType dtype1, DType dtype2) {
    // Handle same type
    if (dtype1 == dtype2) return dtype1;

    // Complex always promotes to complex
    if (dtype_is_complex(dtype1) || dtype_is_complex(dtype2)) {
        if (dtype1 == DTYPE_COMPLEX128 || dtype2 == DTYPE_COMPLEX128 ||
            dtype1 == DTYPE_FLOAT64 || dtype2 == DTYPE_FLOAT64) {
            return DTYPE_COMPLEX128;
        }
        return DTYPE_COMPLEX64;
    }

    // Float promotion
    if (dtype_is_floating(dtype1) || dtype_is_floating(dtype2)) {
        // int64/uint64 + float32 -> float64 (precision)
        if ((dtype1 == DTYPE_INT64 || dtype1 == DTYPE_UINT64 ||
             dtype2 == DTYPE_INT64 || dtype2 == DTYPE_UINT64) &&
            (dtype1 == DTYPE_FLOAT32 || dtype2 == DTYPE_FLOAT32)) {
            return DTYPE_FLOAT64;
        }
        // Otherwise higher float wins
        if (dtype1 == DTYPE_FLOAT64 || dtype2 == DTYPE_FLOAT64) return DTYPE_FLOAT64;
        if (dtype1 == DTYPE_FLOAT32 || dtype2 == DTYPE_FLOAT32) return DTYPE_FLOAT32;
        return DTYPE_FLOAT16;
    }

    // Integer promotion
    // Signed + unsigned of same size -> next larger signed
    const DTypeInfo* info1 = dtype_get_info(dtype1);
    const DTypeInfo* info2 = dtype_get_info(dtype2);

    if (info1->is_signed != info2->is_signed) {
        // Mixed signed/unsigned
        size_t max_size = info1->size > info2->size ? info1->size : info2->size;
        // Need signed type that can hold both
        if (max_size >= 8) return DTYPE_FLOAT64; // Can't fit in int64
        if (max_size >= 4) return DTYPE_INT64;
        if (max_size >= 2) return DTYPE_INT32;
        return DTYPE_INT16;
    }

    // Same signedness: larger wins
    return DTYPE_PRIORITY[dtype1] > DTYPE_PRIORITY[dtype2] ? dtype1 : dtype2;
}
```

### 0.2.1 Flat Index Calculation (C)

```c
// In ndarray.c
EXPORT size_t ndarray_flat_index(NDArray* arr, int32_t* indices, int32_t ndim) {
    if (!arr || ndim != arr->ndim) return SIZE_MAX; // Error

    size_t flat = 0;
    size_t elem_size = dtype_size(arr->dtype);

    for (int32_t i = 0; i < ndim; i++) {
        // Bounds check
        if (indices[i] < 0 || indices[i] >= arr->shape[i]) {
            return SIZE_MAX; // Out of bounds
        }
        flat += (size_t)indices[i] * ((size_t)arr->strides[i] / elem_size);
    }
    return flat;
}

EXPORT bool ndarray_check_bounds(NDArray* arr, int32_t* indices, int32_t ndim) {
    if (!arr || ndim != arr->ndim) return false;

    for (int32_t i = 0; i < ndim; i++) {
        if (indices[i] < 0 || indices[i] >= arr->shape[i]) {
            return false;
        }
    }
    return true;
}
```

### 0.2.2-0.2.3 Element Getters/Setters (C)

```c
// In ndarray.c - Example implementations
EXPORT double ndarray_get_float64(NDArray* arr, size_t flat_idx) {
    if (!arr || flat_idx >= arr->size) return 0.0;
    return ((double*)arr->data)[flat_idx];
}

EXPORT void ndarray_set_float64(NDArray* arr, size_t flat_idx, double value) {
    if (!arr || flat_idx >= arr->size) return;
    if (!(arr->flags & NDARRAY_WRITEABLE)) return;
    ((double*)arr->data)[flat_idx] = value;
}

EXPORT float ndarray_get_complex64_real(NDArray* arr, size_t flat_idx) {
    if (!arr || flat_idx >= arr->size) return 0.0f;
    // Complex64 stored as [real, imag, real, imag, ...]
    return ((float*)arr->data)[flat_idx * 2];
}

EXPORT float ndarray_get_complex64_imag(NDArray* arr, size_t flat_idx) {
    if (!arr || flat_idx >= arr->size) return 0.0f;
    return ((float*)arr->data)[flat_idx * 2 + 1];
}

EXPORT void ndarray_set_complex64(NDArray* arr, size_t flat_idx, float real, float imag) {
    if (!arr || flat_idx >= arr->size) return;
    if (!(arr->flags & NDARRAY_WRITEABLE)) return;
    ((float*)arr->data)[flat_idx * 2] = real;
    ((float*)arr->data)[flat_idx * 2 + 1] = imag;
}
```

### 0.2.4 Generic Element Access (C)

```c
// In ndarray.c
EXPORT double ndarray_get_item(NDArray* arr, int32_t* indices, int32_t ndim) {
    size_t flat_idx = ndarray_flat_index(arr, indices, ndim);
    if (flat_idx == SIZE_MAX) return 0.0;

    switch (arr->dtype) {
        case DTYPE_FLOAT64:    return ndarray_get_float64(arr, flat_idx);
        case DTYPE_FLOAT32:    return (double)ndarray_get_float32(arr, flat_idx);
        case DTYPE_INT64:      return (double)ndarray_get_int64(arr, flat_idx);
        case DTYPE_INT32:      return (double)ndarray_get_int32(arr, flat_idx);
        case DTYPE_INT16:      return (double)ndarray_get_int16(arr, flat_idx);
        case DTYPE_INT8:       return (double)ndarray_get_int8(arr, flat_idx);
        case DTYPE_UINT64:     return (double)ndarray_get_uint64(arr, flat_idx);
        case DTYPE_UINT32:     return (double)ndarray_get_uint32(arr, flat_idx);
        case DTYPE_UINT16:     return (double)ndarray_get_uint16(arr, flat_idx);
        case DTYPE_UINT8:      return (double)ndarray_get_uint8(arr, flat_idx);
        case DTYPE_BOOL:       return (double)ndarray_get_bool(arr, flat_idx);
        case DTYPE_COMPLEX64:  return (double)ndarray_get_complex64_real(arr, flat_idx);
        case DTYPE_COMPLEX128: return ndarray_get_complex128_real(arr, flat_idx);
        default: return 0.0;
    }
}

EXPORT void ndarray_set_item(NDArray* arr, int32_t* indices, int32_t ndim, double value) {
    size_t flat_idx = ndarray_flat_index(arr, indices, ndim);
    if (flat_idx == SIZE_MAX) return;

    switch (arr->dtype) {
        case DTYPE_FLOAT64:    ndarray_set_float64(arr, flat_idx, value); break;
        case DTYPE_FLOAT32:    ndarray_set_float32(arr, flat_idx, (float)value); break;
        case DTYPE_INT64:      ndarray_set_int64(arr, flat_idx, (int64_t)value); break;
        case DTYPE_INT32:      ndarray_set_int32(arr, flat_idx, (int32_t)value); break;
        case DTYPE_INT16:      ndarray_set_int16(arr, flat_idx, (int16_t)value); break;
        case DTYPE_INT8:       ndarray_set_int8(arr, flat_idx, (int8_t)value); break;
        case DTYPE_UINT64:     ndarray_set_uint64(arr, flat_idx, (uint64_t)value); break;
        case DTYPE_UINT32:     ndarray_set_uint32(arr, flat_idx, (uint32_t)value); break;
        case DTYPE_UINT16:     ndarray_set_uint16(arr, flat_idx, (uint16_t)value); break;
        case DTYPE_UINT8:      ndarray_set_uint8(arr, flat_idx, (uint8_t)value); break;
        case DTYPE_BOOL:       ndarray_set_bool(arr, flat_idx, value != 0.0); break;
        case DTYPE_COMPLEX64:  ndarray_set_complex64(arr, flat_idx, (float)value, 0.0f); break;
        case DTYPE_COMPLEX128: ndarray_set_complex128(arr, flat_idx, value, 0.0); break;
    }
}
```

### 0.2.5 TypeScript Element Access

```typescript
// In NDArray.ts
/**
 * Get element at specified indices.
 * @param indices - Multi-dimensional indices
 * @returns Element value
 */
get(...indices: number[]): number {
    this.ensureNotDisposed();

    if (indices.length !== this.ndim) {
        throw new Error(`Expected ${this.ndim} indices, got ${indices.length}`);
    }

    // Allocate indices array in WASM memory
    const indicesPtr = this._module._malloc(indices.length * 4);
    for (let i = 0; i < indices.length; i++) {
        this._module.setValue(indicesPtr + i * 4, indices[i], 'i32');
    }

    const value = this._module._ndarray_get_item(this._ptr, indicesPtr, indices.length);
    this._module._free(indicesPtr);

    return value;
}

/**
 * Set element at specified indices.
 * @param value - Value to set
 * @param indices - Multi-dimensional indices
 */
set(value: number, ...indices: number[]): void {
    this.ensureNotDisposed();

    if (indices.length !== this.ndim) {
        throw new Error(`Expected ${this.ndim} indices, got ${indices.length}`);
    }

    const indicesPtr = this._module._malloc(indices.length * 4);
    for (let i = 0; i < indices.length; i++) {
        this._module.setValue(indicesPtr + i * 4, indices[i], 'i32');
    }

    this._module._ndarray_set_item(this._ptr, indicesPtr, indices.length, value);
    this._module._free(indicesPtr);
}
```

### 0.3.1 empty() Implementation (C)

```c
// In ndarray.c
EXPORT NDArray* ndarray_empty(int32_t ndim, int32_t* shape, DType dtype) {
    NDArray* arr = (NDArray*)malloc(sizeof(NDArray));
    if (!arr) return NULL;

    arr->ndim = ndim;
    arr->dtype = dtype;
    arr->size = compute_size(ndim, shape);
    arr->flags = NDARRAY_OWNDATA | NDARRAY_WRITEABLE;
    arr->base = NULL;

    if (ndim > 0) {
        arr->shape = (int32_t*)malloc((size_t)ndim * sizeof(int32_t));
        arr->strides = (int32_t*)malloc((size_t)ndim * sizeof(int32_t));
        if (!arr->shape || !arr->strides) {
            free(arr->shape);
            free(arr->strides);
            free(arr);
            return NULL;
        }
        memcpy(arr->shape, shape, (size_t)ndim * sizeof(int32_t));
        compute_strides(ndim, shape, dtype, arr->strides);
    } else {
        arr->shape = NULL;
        arr->strides = NULL;
    }

    // Allocate data WITHOUT initialization (key difference from create)
    size_t data_size = arr->size * dtype_size(dtype);
    if (data_size > 0) {
        arr->data = malloc(data_size);
        if (!arr->data) {
            free(arr->shape);
            free(arr->strides);
            free(arr);
            return NULL;
        }
        // NO memset - leave uninitialized for speed
    } else {
        arr->data = NULL;
    }

    // Update contiguity flags
    ndarray_update_contiguity_flags(arr);

    return arr;
}
```

### 0.4.2 Contiguity Checks (C)

```c
// In ndarray.c
EXPORT bool ndarray_is_c_contiguous(NDArray* arr) {
    if (!arr || arr->ndim == 0) return true;

    size_t elem_size = dtype_size(arr->dtype);
    size_t expected = elem_size;

    for (int32_t i = arr->ndim - 1; i >= 0; i--) {
        if (arr->shape[i] == 0) return true; // Empty array
        if (arr->shape[i] != 1 && (size_t)arr->strides[i] != expected) {
            return false;
        }
        expected *= (size_t)arr->shape[i];
    }
    return true;
}

EXPORT bool ndarray_is_f_contiguous(NDArray* arr) {
    if (!arr || arr->ndim == 0) return true;

    size_t elem_size = dtype_size(arr->dtype);
    size_t expected = elem_size;

    for (int32_t i = 0; i < arr->ndim; i++) {
        if (arr->shape[i] == 0) return true;
        if (arr->shape[i] != 1 && (size_t)arr->strides[i] != expected) {
            return false;
        }
        expected *= (size_t)arr->shape[i];
    }
    return true;
}

EXPORT void ndarray_update_contiguity_flags(NDArray* arr) {
    if (!arr) return;

    arr->flags &= ~(NDARRAY_C_CONTIGUOUS | NDARRAY_F_CONTIGUOUS);

    if (ndarray_is_c_contiguous(arr)) {
        arr->flags |= NDARRAY_C_CONTIGUOUS;
    }
    if (ndarray_is_f_contiguous(arr)) {
        arr->flags |= NDARRAY_F_CONTIGUOUS;
    }
}
```

### 0.4.3 Base Pointer for Views (C)

```c
// In ndarray.h - Update struct
typedef struct NDArray_s {
    void* data;
    int32_t ndim;
    int32_t* shape;
    int32_t* strides;
    DType dtype;
    int32_t flags;
    size_t size;
    struct NDArray_s* base;  // NEW: Points to original array if this is a view
} NDArray;

// In ndarray.c
EXPORT NDArray* ndarray_get_base(NDArray* arr) {
    return arr ? arr->base : NULL;
}

EXPORT void ndarray_set_base(NDArray* arr, NDArray* base) {
    if (arr) arr->base = base;
}
```

### 0.5.1 Array Copy (C)

```c
// In ndarray.c
EXPORT NDArray* ndarray_copy(NDArray* arr) {
    if (!arr) return NULL;

    NDArray* copy = ndarray_empty(arr->ndim, arr->shape, arr->dtype);
    if (!copy) return NULL;

    // Copy data
    if (ndarray_is_c_contiguous(arr)) {
        // Fast path: memcpy for contiguous arrays
        memcpy(copy->data, arr->data, arr->size * dtype_size(arr->dtype));
    } else {
        // Slow path: element-by-element for non-contiguous
        ndarray_copy_to(arr, copy);
    }

    return copy;
}

EXPORT void ndarray_copy_to(NDArray* src, NDArray* dst) {
    if (!src || !dst || src->size != dst->size) return;

    // Element-by-element copy (handles different strides)
    for (size_t i = 0; i < src->size; i++) {
        double val;
        switch (src->dtype) {
            case DTYPE_FLOAT64: val = ((double*)src->data)[i]; break;
            case DTYPE_FLOAT32: val = ((float*)src->data)[i]; break;
            case DTYPE_INT32:   val = ((int32_t*)src->data)[i]; break;
            // ... other types
            default: val = 0.0;
        }

        switch (dst->dtype) {
            case DTYPE_FLOAT64: ((double*)dst->data)[i] = val; break;
            case DTYPE_FLOAT32: ((float*)dst->data)[i] = (float)val; break;
            case DTYPE_INT32:   ((int32_t*)dst->data)[i] = (int32_t)val; break;
            // ... other types
        }
    }
}
```

### 0.5.2 Type Casting (C)

```c
// In ndarray.c
EXPORT NDArray* ndarray_astype(NDArray* arr, DType dtype) {
    if (!arr) return NULL;

    // Same dtype: just copy
    if (arr->dtype == dtype) {
        return ndarray_copy(arr);
    }

    NDArray* result = ndarray_empty(arr->ndim, arr->shape, dtype);
    if (!result) return NULL;

    // Convert each element
    for (size_t i = 0; i < arr->size; i++) {
        double val = ndarray_get_item_flat(arr, i);
        ndarray_set_item_flat(result, i, val);
    }

    return result;
}
```

### TypeScript: linspace, eye, etc.

```typescript
// In NDArray.ts

/**
 * Create array of evenly spaced numbers over interval.
 */
static async linspace(
    start: number,
    stop: number,
    num: number = 50,
    endpoint: boolean = true,
    options: NDArrayOptions = {}
): Promise<NDArray> {
    const n = endpoint ? num - 1 : num;
    const step = n > 0 ? (stop - start) / n : 0;

    const data: number[] = [];
    for (let i = 0; i < num; i++) {
        data.push(start + i * step);
    }

    return NDArray.fromArray(data, [num], options);
}

/**
 * Create 2D identity matrix.
 */
static async eye(
    N: number,
    M?: number,
    k: number = 0,
    options: NDArrayOptions = {}
): Promise<NDArray> {
    M = M ?? N;
    const arr = await NDArray.zeros([N, M], options);

    // Set diagonal
    const start = k >= 0 ? k : -k * M;
    const diagLen = Math.min(N, M - Math.max(0, k), N + Math.min(0, k));

    for (let i = 0; i < diagLen; i++) {
        const row = k >= 0 ? i : i - k;
        const col = k >= 0 ? i + k : i;
        if (row >= 0 && row < N && col >= 0 && col < M) {
            arr.set(1, row, col);
        }
    }

    return arr;
}

/**
 * Create array like another array.
 */
static async zerosLike(arr: NDArray, options: NDArrayOptions = {}): Promise<NDArray> {
    return NDArray.zeros(arr.shape, { dtype: options.dtype ?? arr.dtype });
}

static async onesLike(arr: NDArray, options: NDArrayOptions = {}): Promise<NDArray> {
    return NDArray.ones(arr.shape, { dtype: options.dtype ?? arr.dtype });
}

static async emptyLike(arr: NDArray, options: NDArrayOptions = {}): Promise<NDArray> {
    return NDArray.empty(arr.shape, { dtype: options.dtype ?? arr.dtype });
}
```

---

## Implementation Order (Dependency-Driven)

```
Week 1: DType System
├── Day 1-2: C dtype enum extension + dtype_size update
├── Day 3-4: C dtype info structure + utility functions
├── Day 5: TypeScript DType enum + utilities
└── Tests: dtype_size, dtype_is_*, dtype_promote

Week 2: Element Access
├── Day 1-2: C flat_index calculation + bounds checking
├── Day 3-4: C element getters/setters (all types)
├── Day 5: C generic get_item/set_item
└── Tests: get/set for all dtypes, bounds errors

Week 3: TypeScript Element Access + Creation
├── Day 1-2: TS get/set methods
├── Day 3-4: TS empty, full, *Like methods
├── Day 5: TS linspace, logspace, geomspace
└── Tests: element access, creation methods

Week 4: Memory, Flags & Copy
├── Day 1-2: C flags extension + contiguity checks
├── Day 3: C base pointer for views
├── Day 4: C/TS copy and astype
├── Day 5: TS eye, identity, diag, tri
└── Tests: flags, copy, astype, matrix creation

Week 5: Integration & Polish
├── Day 1-2: Complex number support (c64, c128)
├── Day 3-4: Full test coverage
├── Day 5: Documentation + benchmark updates
└── Final: Ensure all Level 0 functions work together
```

---

## Verification Plan

After Level 0 completion, verify:

```bash
# Run all tests
npm test

# Specific Level 0 tests should pass:
✓ DType enum has all 14 types
✓ dtype_size returns correct sizes for all types
✓ Type promotion follows NumPy rules
✓ Element get/set works for all dtypes
✓ Element access with multi-dimensional indices
✓ Bounds checking throws appropriate errors
✓ empty() creates uninitialized array
✓ full() creates array with fill value
✓ *Like() methods copy shape/dtype correctly
✓ linspace/logspace/geomspace produce correct sequences
✓ eye/identity/diag/tri produce correct matrices
✓ copy() creates independent deep copy
✓ astype() converts between all compatible dtypes
✓ flags property reports correct contiguity
✓ Complex number get/set works correctly
```

Generate NumPy comparison vectors:

```python
# tests/python/generate_level0_tests.py
import numpy as np
import json

tests = {
    "dtype_sizes": {
        "bool": 1, "int8": 1, "int16": 2, "int32": 4, "int64": 8,
        "uint8": 1, "uint16": 2, "uint32": 4, "uint64": 8,
        "float32": 4, "float64": 8, "complex64": 8, "complex128": 16
    },
    "promotion": [
        {"a": "int32", "b": "float32", "result": "float64"},
        {"a": "int64", "b": "float32", "result": "float64"},
        {"a": "uint8", "b": "int8", "result": "int16"},
        {"a": "float32", "b": "complex64", "result": "complex64"},
    ],
    "linspace": [
        {"start": 0, "stop": 10, "num": 5, "result": [0, 2.5, 5, 7.5, 10]},
        {"start": 0, "stop": 1, "num": 5, "endpoint": False, "result": [0, 0.2, 0.4, 0.6, 0.8]},
    ],
    "eye": [
        {"N": 3, "result": [[1,0,0],[0,1,0],[0,0,1]]},
        {"N": 3, "M": 4, "result": [[1,0,0,0],[0,1,0,0],[0,0,1,0]]},
        {"N": 3, "k": 1, "result": [[0,1,0],[0,0,1],[0,0,0]]},
    ],
}

with open("tests/fixtures/level0_vectors.json", "w") as f:
    json.dump(tests, f, indent=2)
```

---

## Critical Dependencies for Level 1

Level 0 completion enables:

- **Level 1.1 Iterators**: Requires element access (0.2)
- **Level 1.2 Shape Manipulation**: Requires flags, base pointer (0.4)
- **Level 1.3 Views**: Requires base pointer, contiguity checks (0.4)

Level 0 MUST be complete before proceeding to Level 1.
