/*
 * NumJS DType Implementation
 *
 * Data type utilities, promotion rules, and casting logic.
 */

#include "dtype.h"
#include <stdint.h>
#include <limits.h>

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#define EXPORT EMSCRIPTEN_KEEPALIVE
#else
#define EXPORT
#endif

/* ============ DType Info Table ============ */

/*
 * Static table of dtype information.
 * Order must match DType enum values.
 */
static const DTypeInfo DTYPE_INFO[DTYPE_COUNT] = {
    /* DTYPE_FLOAT32 = 0 */
    { .size = 4,  .alignment = 4, .name = "float32",    .kind = DTYPE_KIND_FLOAT,   .is_signed = true  },
    /* DTYPE_FLOAT64 = 1 */
    { .size = 8,  .alignment = 8, .name = "float64",    .kind = DTYPE_KIND_FLOAT,   .is_signed = true  },
    /* DTYPE_INT32 = 2 */
    { .size = 4,  .alignment = 4, .name = "int32",      .kind = DTYPE_KIND_INT,     .is_signed = true  },
    /* DTYPE_INT64 = 3 */
    { .size = 8,  .alignment = 8, .name = "int64",      .kind = DTYPE_KIND_INT,     .is_signed = true  },
    /* DTYPE_BOOL = 4 */
    { .size = 1,  .alignment = 1, .name = "bool",       .kind = DTYPE_KIND_BOOL,    .is_signed = false },
    /* DTYPE_INT8 = 5 */
    { .size = 1,  .alignment = 1, .name = "int8",       .kind = DTYPE_KIND_INT,     .is_signed = true  },
    /* DTYPE_INT16 = 6 */
    { .size = 2,  .alignment = 2, .name = "int16",      .kind = DTYPE_KIND_INT,     .is_signed = true  },
    /* DTYPE_UINT8 = 7 */
    { .size = 1,  .alignment = 1, .name = "uint8",      .kind = DTYPE_KIND_UINT,    .is_signed = false },
    /* DTYPE_UINT16 = 8 */
    { .size = 2,  .alignment = 2, .name = "uint16",     .kind = DTYPE_KIND_UINT,    .is_signed = false },
    /* DTYPE_UINT32 = 9 */
    { .size = 4,  .alignment = 4, .name = "uint32",     .kind = DTYPE_KIND_UINT,    .is_signed = false },
    /* DTYPE_UINT64 = 10 */
    { .size = 8,  .alignment = 8, .name = "uint64",     .kind = DTYPE_KIND_UINT,    .is_signed = false },
    /* DTYPE_FLOAT16 = 11 */
    { .size = 2,  .alignment = 2, .name = "float16",    .kind = DTYPE_KIND_FLOAT,   .is_signed = true  },
    /* DTYPE_COMPLEX64 = 12 */
    { .size = 8,  .alignment = 4, .name = "complex64",  .kind = DTYPE_KIND_COMPLEX, .is_signed = true  },
    /* DTYPE_COMPLEX128 = 13 */
    { .size = 16, .alignment = 8, .name = "complex128", .kind = DTYPE_KIND_COMPLEX, .is_signed = true  },
};

/* ============ Basic Functions ============ */

EXPORT size_t dtype_size(DType dtype)
{
    if (dtype < 0 || dtype >= DTYPE_COUNT) {
        return 0;
    }
    return DTYPE_INFO[dtype].size;
}

EXPORT const DTypeInfo* dtype_get_info(DType dtype)
{
    if (dtype < 0 || dtype >= DTYPE_COUNT) {
        return NULL;
    }
    return &DTYPE_INFO[dtype];
}

EXPORT const char* dtype_name(DType dtype)
{
    const DTypeInfo* info = dtype_get_info(dtype);
    return info ? info->name : "unknown";
}

/* ============ Type Predicates ============ */

EXPORT bool dtype_is_integer(DType dtype)
{
    const DTypeInfo* info = dtype_get_info(dtype);
    if (!info) return false;
    return info->kind == DTYPE_KIND_INT || info->kind == DTYPE_KIND_UINT;
}

EXPORT bool dtype_is_floating(DType dtype)
{
    const DTypeInfo* info = dtype_get_info(dtype);
    if (!info) return false;
    return info->kind == DTYPE_KIND_FLOAT;
}

EXPORT bool dtype_is_complex(DType dtype)
{
    const DTypeInfo* info = dtype_get_info(dtype);
    if (!info) return false;
    return info->kind == DTYPE_KIND_COMPLEX;
}

EXPORT bool dtype_is_signed(DType dtype)
{
    const DTypeInfo* info = dtype_get_info(dtype);
    if (!info) return false;
    return info->is_signed;
}

EXPORT bool dtype_is_bool(DType dtype)
{
    const DTypeInfo* info = dtype_get_info(dtype);
    if (!info) return false;
    return info->kind == DTYPE_KIND_BOOL;
}

EXPORT bool dtype_is_numeric(DType dtype)
{
    return dtype_is_integer(dtype) || dtype_is_floating(dtype) || dtype_is_complex(dtype);
}

/* ============ Type Promotion ============ */

/*
 * Priority values for type promotion.
 * Higher priority wins when combining types.
 */
static const int DTYPE_PRIORITY[DTYPE_COUNT] = {
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

EXPORT DType dtype_promote(DType dtype1, DType dtype2)
{
    /* Invalid types */
    if (dtype1 < 0 || dtype1 >= DTYPE_COUNT ||
        dtype2 < 0 || dtype2 >= DTYPE_COUNT) {
        return DTYPE_FLOAT64; /* Default fallback */
    }

    /* Same type: no promotion needed */
    if (dtype1 == dtype2) {
        return dtype1;
    }

    const DTypeInfo* info1 = &DTYPE_INFO[dtype1];
    const DTypeInfo* info2 = &DTYPE_INFO[dtype2];

    /* Complex always promotes to complex */
    if (info1->kind == DTYPE_KIND_COMPLEX || info2->kind == DTYPE_KIND_COMPLEX) {
        /* If either is complex128 or float64, result is complex128 */
        if (dtype1 == DTYPE_COMPLEX128 || dtype2 == DTYPE_COMPLEX128 ||
            dtype1 == DTYPE_FLOAT64 || dtype2 == DTYPE_FLOAT64 ||
            dtype1 == DTYPE_INT64 || dtype2 == DTYPE_INT64 ||
            dtype1 == DTYPE_UINT64 || dtype2 == DTYPE_UINT64) {
            return DTYPE_COMPLEX128;
        }
        return DTYPE_COMPLEX64;
    }

    /* Float promotion */
    if (info1->kind == DTYPE_KIND_FLOAT || info2->kind == DTYPE_KIND_FLOAT) {
        /* int64/uint64 + float32 -> float64 (precision) */
        bool has_large_int = (dtype1 == DTYPE_INT64 || dtype1 == DTYPE_UINT64 ||
                              dtype2 == DTYPE_INT64 || dtype2 == DTYPE_UINT64 ||
                              dtype1 == DTYPE_INT32 || dtype1 == DTYPE_UINT32 ||
                              dtype2 == DTYPE_INT32 || dtype2 == DTYPE_UINT32);
        bool has_float32 = (dtype1 == DTYPE_FLOAT32 || dtype2 == DTYPE_FLOAT32);

        if (has_large_int && has_float32) {
            return DTYPE_FLOAT64;
        }

        /* Otherwise higher float wins */
        if (dtype1 == DTYPE_FLOAT64 || dtype2 == DTYPE_FLOAT64) {
            return DTYPE_FLOAT64;
        }
        if (dtype1 == DTYPE_FLOAT32 || dtype2 == DTYPE_FLOAT32) {
            return DTYPE_FLOAT32;
        }
        return DTYPE_FLOAT16;
    }

    /* Bool promotion: bool + any numeric -> that numeric */
    if (info1->kind == DTYPE_KIND_BOOL) {
        return dtype2;
    }
    if (info2->kind == DTYPE_KIND_BOOL) {
        return dtype1;
    }

    /* Integer promotion */
    /* Mixed signed/unsigned: need type that can hold both ranges */
    if (info1->is_signed != info2->is_signed) {
        size_t max_size = info1->size > info2->size ? info1->size : info2->size;

        /* Unsigned larger or equal to signed: need next larger signed */
        DType unsigned_dt = info1->is_signed ? dtype2 : dtype1;
        DType signed_dt = info1->is_signed ? dtype1 : dtype2;
        size_t unsigned_size = DTYPE_INFO[unsigned_dt].size;
        size_t signed_size = DTYPE_INFO[signed_dt].size;

        if (unsigned_size >= signed_size) {
            /* Need larger signed type to hold unsigned range */
            if (unsigned_size >= 8) {
                /* uint64 can't fit in int64, use float64 */
                return DTYPE_FLOAT64;
            }
            if (unsigned_size >= 4) {
                return DTYPE_INT64;
            }
            if (unsigned_size >= 2) {
                return DTYPE_INT32;
            }
            return DTYPE_INT16;
        } else {
            /* Signed is larger, can hold unsigned */
            return signed_dt;
        }
    }

    /* Same signedness: larger size wins */
    return DTYPE_PRIORITY[dtype1] > DTYPE_PRIORITY[dtype2] ? dtype1 : dtype2;
}

EXPORT DType dtype_common_type(DType* dtypes, size_t count)
{
    if (count == 0) {
        return DTYPE_FLOAT64; /* Default */
    }

    DType result = dtypes[0];
    for (size_t i = 1; i < count; i++) {
        result = dtype_promote(result, dtypes[i]);
    }
    return result;
}

/* ============ Casting ============ */

EXPORT bool dtype_can_cast(DType from, DType to, CastingKind casting)
{
    if (from < 0 || from >= DTYPE_COUNT || to < 0 || to >= DTYPE_COUNT) {
        return false;
    }

    /* Same type always allowed */
    if (from == to) {
        return true;
    }

    switch (casting) {
        case CASTING_NO:
            return false;

        case CASTING_EQUIV:
            /* Only allow equivalent types (same size, same kind) */
            return DTYPE_INFO[from].size == DTYPE_INFO[to].size &&
                   DTYPE_INFO[from].kind == DTYPE_INFO[to].kind;

        case CASTING_SAFE:
            /* Allow if promotion would yield 'to' type */
            return dtype_promote(from, to) == to;

        case CASTING_SAME_KIND:
            /* Allow within same kind, or safe casts */
            if (DTYPE_INFO[from].kind == DTYPE_INFO[to].kind) {
                return true;
            }
            return dtype_promote(from, to) == to;

        case CASTING_UNSAFE:
            return true;
    }

    return false;
}

/* ============ Minimum Type for Value ============ */

EXPORT DType dtype_min_int(int64_t value)
{
    if (value >= INT8_MIN && value <= INT8_MAX) {
        return DTYPE_INT8;
    }
    if (value >= INT16_MIN && value <= INT16_MAX) {
        return DTYPE_INT16;
    }
    if (value >= INT32_MIN && value <= INT32_MAX) {
        return DTYPE_INT32;
    }
    return DTYPE_INT64;
}

EXPORT DType dtype_min_uint(uint64_t value)
{
    if (value <= UINT8_MAX) {
        return DTYPE_UINT8;
    }
    if (value <= UINT16_MAX) {
        return DTYPE_UINT16;
    }
    if (value <= UINT32_MAX) {
        return DTYPE_UINT32;
    }
    return DTYPE_UINT64;
}
