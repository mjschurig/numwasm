#ifndef NUMJS_DTYPE_H
#define NUMJS_DTYPE_H

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

/*
 * NumJS DType System
 *
 * Provides data type information, promotion rules, and casting utilities.
 * Inspired by NumPy's dtype system.
 */

/* ============ DType Enum ============ */

/*
 * Data type enumeration - matches TypeScript DType enum.
 * Values must stay synchronized with src/ts/types.ts
 */
typedef enum {
    /* Original types (maintain backward compatibility) */
    DTYPE_FLOAT32 = 0,
    DTYPE_FLOAT64 = 1,
    DTYPE_INT32 = 2,
    DTYPE_INT64 = 3,

    /* Boolean */
    DTYPE_BOOL = 4,

    /* Additional integers */
    DTYPE_INT8 = 5,
    DTYPE_INT16 = 6,
    DTYPE_UINT8 = 7,
    DTYPE_UINT16 = 8,
    DTYPE_UINT32 = 9,
    DTYPE_UINT64 = 10,

    /* Half precision (optional - complex to implement fully) */
    DTYPE_FLOAT16 = 11,

    /* Complex types */
    DTYPE_COMPLEX64 = 12,   /* float real + float imag */
    DTYPE_COMPLEX128 = 13,  /* double real + double imag */

    /* Count for iteration */
    DTYPE_COUNT = 14
} DType;

/* ============ DType Kind ============ */

/*
 * Category of data type.
 * Uses NumPy's character codes for compatibility.
 */
typedef enum {
    DTYPE_KIND_BOOL = 'b',
    DTYPE_KIND_INT = 'i',
    DTYPE_KIND_UINT = 'u',
    DTYPE_KIND_FLOAT = 'f',
    DTYPE_KIND_COMPLEX = 'c'
} DTypeKind;

/* ============ Casting Kind ============ */

/*
 * Casting safety levels (matches NumPy).
 */
typedef enum {
    CASTING_NO = 0,        /* No casting allowed */
    CASTING_EQUIV = 1,     /* Only byte-order changes */
    CASTING_SAFE = 2,      /* Only casts preserving values */
    CASTING_SAME_KIND = 3, /* Safe casts or within same kind */
    CASTING_UNSAFE = 4     /* Any data conversion */
} CastingKind;

/* ============ DType Info ============ */

/*
 * Information about a data type.
 */
typedef struct {
    size_t size;           /* Size in bytes */
    size_t alignment;      /* Alignment requirement */
    const char* name;      /* String name (e.g., "float64") */
    DTypeKind kind;        /* Category */
    bool is_signed;        /* True for signed types */
} DTypeInfo;

/* ============ DType Functions ============ */

/*
 * Get size of a dtype in bytes.
 */
size_t dtype_size(DType dtype);

/*
 * Get information about a dtype.
 * Returns NULL for invalid dtypes.
 */
const DTypeInfo* dtype_get_info(DType dtype);

/*
 * Get the name of a dtype as a string.
 * Returns "unknown" for invalid dtypes.
 */
const char* dtype_name(DType dtype);

/* ============ Type Predicates ============ */

/*
 * Check if dtype is an integer type (signed or unsigned).
 */
bool dtype_is_integer(DType dtype);

/*
 * Check if dtype is a floating point type.
 */
bool dtype_is_floating(DType dtype);

/*
 * Check if dtype is a complex type.
 */
bool dtype_is_complex(DType dtype);

/*
 * Check if dtype is signed (integers and floats).
 */
bool dtype_is_signed(DType dtype);

/*
 * Check if dtype is boolean.
 */
bool dtype_is_bool(DType dtype);

/*
 * Check if dtype is numeric (not bool).
 */
bool dtype_is_numeric(DType dtype);

/* ============ Type Promotion ============ */

/*
 * Get result type when combining two dtypes.
 * Follows NumPy promotion rules.
 */
DType dtype_promote(DType dtype1, DType dtype2);

/*
 * Get common type for an array of dtypes.
 */
DType dtype_common_type(DType* dtypes, size_t count);

/*
 * Check if casting from one dtype to another is allowed.
 */
bool dtype_can_cast(DType from, DType to, CastingKind casting);

/*
 * Get the minimum dtype that can hold the given integer value.
 */
DType dtype_min_int(int64_t value);

/*
 * Get the minimum unsigned dtype that can hold the given value.
 */
DType dtype_min_uint(uint64_t value);

#endif /* NUMJS_DTYPE_H */
