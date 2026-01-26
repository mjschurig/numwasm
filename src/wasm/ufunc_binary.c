/**
 * NumJS-WASM Binary Universal Functions Implementation
 *
 * Implements ~25 binary element-wise operations with broadcasting support.
 */

#include "ufunc_binary.h"
#include "ufunc.h"
#include "dtype.h"
#include <math.h>
#include <stdint.h>

/* ============ Arithmetic Loops ============ */

/* add: x + y */
BINARY_LOOP_F64(add, x + y)
BINARY_LOOP_F32(add, x + y)
BINARY_LOOP_I32(add, x + y)

/* subtract: x - y */
BINARY_LOOP_F64(subtract, x - y)
BINARY_LOOP_F32(subtract, x - y)
BINARY_LOOP_I32(subtract, x - y)

/* multiply: x * y */
BINARY_LOOP_F64(multiply, x * y)
BINARY_LOOP_F32(multiply, x * y)
BINARY_LOOP_I32(multiply, x * y)

/* divide: x / y (true division) */
BINARY_LOOP_F64(divide, x / y)
BINARY_LOOP_F32(divide, x / y)
/* Integer division returns float for true divide */
static void ufunc_loop_divide_i32(const char* in1, const char* in2, char* out,
                                   size_t n, int64_t s1, int64_t s2, int64_t so) {
    for (size_t i = 0; i < n; i++) {
        int32_t x = *(const int32_t*)(in1 + i * s1);
        int32_t y = *(const int32_t*)(in2 + i * s2);
        *(double*)(out + i * so) = (y != 0) ? ((double)x / (double)y) : (x >= 0 ? INFINITY : -INFINITY);
    }
}

/* floor_divide: floor(x / y) */
BINARY_LOOP_F64(floor_divide, floor(x / y))
BINARY_LOOP_F32(floor_divide, floorf(x / y))
static void ufunc_loop_floor_divide_i32(const char* in1, const char* in2, char* out,
                                         size_t n, int64_t s1, int64_t s2, int64_t so) {
    for (size_t i = 0; i < n; i++) {
        int32_t x = *(const int32_t*)(in1 + i * s1);
        int32_t y = *(const int32_t*)(in2 + i * s2);
        if (y == 0) {
            *(int32_t*)(out + i * so) = 0;  /* Undefined behavior, return 0 */
        } else {
            /* Floor division for integers (Python-style) */
            int32_t q = x / y;
            int32_t r = x % y;
            /* Adjust for negative remainders */
            if ((r != 0) && ((r ^ y) < 0)) {
                q -= 1;
            }
            *(int32_t*)(out + i * so) = q;
        }
    }
}

/* remainder (modulo): x % y (Python-style) */
BINARY_LOOP_F64(remainder, fmod(x, y) + ((fmod(x, y) != 0 && ((x < 0) != (y < 0))) ? y : 0))
BINARY_LOOP_F32(remainder, fmodf(x, y) + ((fmodf(x, y) != 0 && ((x < 0) != (y < 0))) ? y : 0))
static void ufunc_loop_remainder_i32(const char* in1, const char* in2, char* out,
                                      size_t n, int64_t s1, int64_t s2, int64_t so) {
    for (size_t i = 0; i < n; i++) {
        int32_t x = *(const int32_t*)(in1 + i * s1);
        int32_t y = *(const int32_t*)(in2 + i * s2);
        if (y == 0) {
            *(int32_t*)(out + i * so) = 0;
        } else {
            int32_t r = x % y;
            /* Python-style modulo: result has same sign as divisor */
            if ((r != 0) && ((r ^ y) < 0)) {
                r += y;
            }
            *(int32_t*)(out + i * so) = r;
        }
    }
}

/* fmod: C-style remainder (same sign as dividend) */
BINARY_LOOP_F64(fmod, fmod(x, y))
BINARY_LOOP_F32(fmod, fmodf(x, y))
static void ufunc_loop_fmod_i32(const char* in1, const char* in2, char* out,
                                 size_t n, int64_t s1, int64_t s2, int64_t so) {
    for (size_t i = 0; i < n; i++) {
        int32_t x = *(const int32_t*)(in1 + i * s1);
        int32_t y = *(const int32_t*)(in2 + i * s2);
        *(int32_t*)(out + i * so) = (y != 0) ? (x % y) : 0;
    }
}

/* power: x^y */
BINARY_LOOP_F64(power, pow(x, y))
BINARY_LOOP_F32(power, powf(x, y))
static void ufunc_loop_power_i32(const char* in1, const char* in2, char* out,
                                  size_t n, int64_t s1, int64_t s2, int64_t so) {
    for (size_t i = 0; i < n; i++) {
        int32_t x = *(const int32_t*)(in1 + i * s1);
        int32_t y = *(const int32_t*)(in2 + i * s2);
        *(int32_t*)(out + i * so) = (int32_t)pow((double)x, (double)y);
    }
}

/* ============ Comparison Loops ============ */

/* equal: x == y */
BINARY_CMP_LOOP_F64(equal, x == y)
BINARY_CMP_LOOP_F32(equal, x == y)
BINARY_CMP_LOOP_I32(equal, x == y)

/* not_equal: x != y */
BINARY_CMP_LOOP_F64(not_equal, x != y)
BINARY_CMP_LOOP_F32(not_equal, x != y)
BINARY_CMP_LOOP_I32(not_equal, x != y)

/* less: x < y */
BINARY_CMP_LOOP_F64(less, x < y)
BINARY_CMP_LOOP_F32(less, x < y)
BINARY_CMP_LOOP_I32(less, x < y)

/* less_equal: x <= y */
BINARY_CMP_LOOP_F64(less_equal, x <= y)
BINARY_CMP_LOOP_F32(less_equal, x <= y)
BINARY_CMP_LOOP_I32(less_equal, x <= y)

/* greater: x > y */
BINARY_CMP_LOOP_F64(greater, x > y)
BINARY_CMP_LOOP_F32(greater, x > y)
BINARY_CMP_LOOP_I32(greater, x > y)

/* greater_equal: x >= y */
BINARY_CMP_LOOP_F64(greater_equal, x >= y)
BINARY_CMP_LOOP_F32(greater_equal, x >= y)
BINARY_CMP_LOOP_I32(greater_equal, x >= y)

/* ============ Extrema Loops ============ */

/* maximum: max(x, y), propagates NaN */
BINARY_LOOP_F64(maximum, (isnan(x) || isnan(y)) ? NAN : ((x > y) ? x : y))
BINARY_LOOP_F32(maximum, (isnan(x) || isnan(y)) ? NAN : ((x > y) ? x : y))
BINARY_LOOP_I32(maximum, (x > y) ? x : y)

/* minimum: min(x, y), propagates NaN */
BINARY_LOOP_F64(minimum, (isnan(x) || isnan(y)) ? NAN : ((x < y) ? x : y))
BINARY_LOOP_F32(minimum, (isnan(x) || isnan(y)) ? NAN : ((x < y) ? x : y))
BINARY_LOOP_I32(minimum, (x < y) ? x : y)

/* fmax: max ignoring NaN (use fmax which handles NaN) */
BINARY_LOOP_F64(fmax, fmax(x, y))
BINARY_LOOP_F32(fmax, fmaxf(x, y))
BINARY_LOOP_I32(fmax, (x > y) ? x : y)

/* fmin: min ignoring NaN */
BINARY_LOOP_F64(fmin, fmin(x, y))
BINARY_LOOP_F32(fmin, fminf(x, y))
BINARY_LOOP_I32(fmin, (x < y) ? x : y)

/* ============ Logical Loops ============ */

/* logical_and: x && y (output is bool) */
BINARY_CMP_LOOP_F64(logical_and, (x != 0) && (y != 0))
BINARY_CMP_LOOP_F32(logical_and, (x != 0) && (y != 0))
BINARY_CMP_LOOP_I32(logical_and, (x != 0) && (y != 0))

/* logical_or: x || y */
BINARY_CMP_LOOP_F64(logical_or, (x != 0) || (y != 0))
BINARY_CMP_LOOP_F32(logical_or, (x != 0) || (y != 0))
BINARY_CMP_LOOP_I32(logical_or, (x != 0) || (y != 0))

/* logical_xor: x ^ y (logical) */
BINARY_CMP_LOOP_F64(logical_xor, ((x != 0) != (y != 0)))
BINARY_CMP_LOOP_F32(logical_xor, ((x != 0) != (y != 0)))
BINARY_CMP_LOOP_I32(logical_xor, ((x != 0) != (y != 0)))

/* ============ Bitwise Loops (Integer only) ============ */

static void ufunc_loop_bitwise_and_i32(const char* in1, const char* in2, char* out,
                                        size_t n, int64_t s1, int64_t s2, int64_t so) {
    for (size_t i = 0; i < n; i++) {
        int32_t x = *(const int32_t*)(in1 + i * s1);
        int32_t y = *(const int32_t*)(in2 + i * s2);
        *(int32_t*)(out + i * so) = x & y;
    }
}

static void ufunc_loop_bitwise_or_i32(const char* in1, const char* in2, char* out,
                                       size_t n, int64_t s1, int64_t s2, int64_t so) {
    for (size_t i = 0; i < n; i++) {
        int32_t x = *(const int32_t*)(in1 + i * s1);
        int32_t y = *(const int32_t*)(in2 + i * s2);
        *(int32_t*)(out + i * so) = x | y;
    }
}

static void ufunc_loop_bitwise_xor_i32(const char* in1, const char* in2, char* out,
                                        size_t n, int64_t s1, int64_t s2, int64_t so) {
    for (size_t i = 0; i < n; i++) {
        int32_t x = *(const int32_t*)(in1 + i * s1);
        int32_t y = *(const int32_t*)(in2 + i * s2);
        *(int32_t*)(out + i * so) = x ^ y;
    }
}

static void ufunc_loop_left_shift_i32(const char* in1, const char* in2, char* out,
                                       size_t n, int64_t s1, int64_t s2, int64_t so) {
    for (size_t i = 0; i < n; i++) {
        int32_t x = *(const int32_t*)(in1 + i * s1);
        int32_t y = *(const int32_t*)(in2 + i * s2);
        *(int32_t*)(out + i * so) = (y >= 0 && y < 32) ? (x << y) : 0;
    }
}

static void ufunc_loop_right_shift_i32(const char* in1, const char* in2, char* out,
                                        size_t n, int64_t s1, int64_t s2, int64_t so) {
    for (size_t i = 0; i < n; i++) {
        int32_t x = *(const int32_t*)(in1 + i * s1);
        int32_t y = *(const int32_t*)(in2 + i * s2);
        *(int32_t*)(out + i * so) = (y >= 0 && y < 32) ? (x >> y) : 0;
    }
}

/* ============ Special Math Loops ============ */

/* arctan2: atan2(x, y) */
BINARY_LOOP_F64(arctan2, atan2(x, y))
BINARY_LOOP_F32(arctan2, atan2f(x, y))
static void ufunc_loop_arctan2_i32(const char* in1, const char* in2, char* out,
                                    size_t n, int64_t s1, int64_t s2, int64_t so) {
    for (size_t i = 0; i < n; i++) {
        int32_t x = *(const int32_t*)(in1 + i * s1);
        int32_t y = *(const int32_t*)(in2 + i * s2);
        *(double*)(out + i * so) = atan2((double)x, (double)y);
    }
}

/* hypot: sqrt(x^2 + y^2) */
BINARY_LOOP_F64(hypot, hypot(x, y))
BINARY_LOOP_F32(hypot, hypotf(x, y))
static void ufunc_loop_hypot_i32(const char* in1, const char* in2, char* out,
                                  size_t n, int64_t s1, int64_t s2, int64_t so) {
    for (size_t i = 0; i < n; i++) {
        int32_t x = *(const int32_t*)(in1 + i * s1);
        int32_t y = *(const int32_t*)(in2 + i * s2);
        *(double*)(out + i * so) = hypot((double)x, (double)y);
    }
}

/* copysign: sign of y applied to x */
BINARY_LOOP_F64(copysign, copysign(x, y))
BINARY_LOOP_F32(copysign, copysignf(x, y))
static void ufunc_loop_copysign_i32(const char* in1, const char* in2, char* out,
                                     size_t n, int64_t s1, int64_t s2, int64_t so) {
    for (size_t i = 0; i < n; i++) {
        int32_t x = *(const int32_t*)(in1 + i * s1);
        int32_t y = *(const int32_t*)(in2 + i * s2);
        int32_t ax = (x < 0) ? -x : x;
        *(int32_t*)(out + i * so) = (y < 0) ? -ax : ax;
    }
}

/* logaddexp: log(exp(x) + exp(y)) */
BINARY_LOOP_F64(logaddexp, log(exp(x) + exp(y)))
BINARY_LOOP_F32(logaddexp, logf(expf(x) + expf(y)))
static void ufunc_loop_logaddexp_i32(const char* in1, const char* in2, char* out,
                                      size_t n, int64_t s1, int64_t s2, int64_t so) {
    for (size_t i = 0; i < n; i++) {
        int32_t x = *(const int32_t*)(in1 + i * s1);
        int32_t y = *(const int32_t*)(in2 + i * s2);
        *(double*)(out + i * so) = log(exp((double)x) + exp((double)y));
    }
}

/* logaddexp2: log2(2^x + 2^y) */
BINARY_LOOP_F64(logaddexp2, log2(exp2(x) + exp2(y)))
BINARY_LOOP_F32(logaddexp2, log2f(exp2f(x) + exp2f(y)))
static void ufunc_loop_logaddexp2_i32(const char* in1, const char* in2, char* out,
                                       size_t n, int64_t s1, int64_t s2, int64_t so) {
    for (size_t i = 0; i < n; i++) {
        int32_t x = *(const int32_t*)(in1 + i * s1);
        int32_t y = *(const int32_t*)(in2 + i * s2);
        *(double*)(out + i * so) = log2(exp2((double)x) + exp2((double)y));
    }
}

/* ============ Exported Functions ============ */

/* Helper macro for binary ufuncs that preserve dtype */
#define IMPL_BINARY_UFUNC(name) \
NDArray* ufunc_##name(NDArray* x1, NDArray* x2) { \
    if (x1 == NULL || x2 == NULL) return NULL; \
    DType out_dtype = dtype_promote(x1->dtype, x2->dtype); \
    BinaryLoopFunc loop = ufunc_select_binary_loop( \
        out_dtype, \
        ufunc_loop_##name##_f64, \
        ufunc_loop_##name##_f32, \
        ufunc_loop_##name##_i32); \
    return ufunc_apply_binary(x1, x2, out_dtype, loop); \
}

/* Helper macro for comparison ufuncs (always return bool) */
#define IMPL_BINARY_CMP(name) \
NDArray* ufunc_##name(NDArray* x1, NDArray* x2) { \
    if (x1 == NULL || x2 == NULL) return NULL; \
    DType common_dtype = dtype_promote(x1->dtype, x2->dtype); \
    BinaryLoopFunc loop = ufunc_select_binary_loop( \
        common_dtype, \
        ufunc_loop_##name##_f64, \
        ufunc_loop_##name##_f32, \
        ufunc_loop_##name##_i32); \
    return ufunc_apply_binary_cmp(x1, x2, loop); \
}

/* Arithmetic */
IMPL_BINARY_UFUNC(add)
IMPL_BINARY_UFUNC(subtract)
IMPL_BINARY_UFUNC(multiply)

/* divide: for integers, output is float64 */
NDArray* ufunc_divide(NDArray* x1, NDArray* x2) {
    if (x1 == NULL || x2 == NULL) return NULL;
    DType common_dtype = dtype_promote(x1->dtype, x2->dtype);
    DType out_dtype;
    BinaryLoopFunc loop;

    if (dtype_is_integer(common_dtype)) {
        /* True division of integers produces float64 */
        out_dtype = DTYPE_FLOAT64;
        loop = ufunc_loop_divide_i32;
    } else {
        out_dtype = common_dtype;
        loop = ufunc_select_binary_loop(common_dtype,
            ufunc_loop_divide_f64, ufunc_loop_divide_f32, ufunc_loop_divide_i32);
    }
    return ufunc_apply_binary(x1, x2, out_dtype, loop);
}

NDArray* ufunc_true_divide(NDArray* x1, NDArray* x2) { return ufunc_divide(x1, x2); }

IMPL_BINARY_UFUNC(floor_divide)
IMPL_BINARY_UFUNC(remainder)
IMPL_BINARY_UFUNC(fmod)
IMPL_BINARY_UFUNC(power)

NDArray* ufunc_mod(NDArray* x1, NDArray* x2) { return ufunc_remainder(x1, x2); }

/* Comparisons */
IMPL_BINARY_CMP(equal)
IMPL_BINARY_CMP(not_equal)
IMPL_BINARY_CMP(less)
IMPL_BINARY_CMP(less_equal)
IMPL_BINARY_CMP(greater)
IMPL_BINARY_CMP(greater_equal)

/* Extrema */
IMPL_BINARY_UFUNC(maximum)
IMPL_BINARY_UFUNC(minimum)
IMPL_BINARY_UFUNC(fmax)
IMPL_BINARY_UFUNC(fmin)

/* Logical (return bool) */
IMPL_BINARY_CMP(logical_and)
IMPL_BINARY_CMP(logical_or)
IMPL_BINARY_CMP(logical_xor)

/* Bitwise (integer only) */
NDArray* ufunc_bitwise_and(NDArray* x1, NDArray* x2) {
    if (x1 == NULL || x2 == NULL) return NULL;
    return ufunc_apply_binary(x1, x2, DTYPE_INT32, ufunc_loop_bitwise_and_i32);
}

NDArray* ufunc_bitwise_or(NDArray* x1, NDArray* x2) {
    if (x1 == NULL || x2 == NULL) return NULL;
    return ufunc_apply_binary(x1, x2, DTYPE_INT32, ufunc_loop_bitwise_or_i32);
}

NDArray* ufunc_bitwise_xor(NDArray* x1, NDArray* x2) {
    if (x1 == NULL || x2 == NULL) return NULL;
    return ufunc_apply_binary(x1, x2, DTYPE_INT32, ufunc_loop_bitwise_xor_i32);
}

NDArray* ufunc_left_shift(NDArray* x1, NDArray* x2) {
    if (x1 == NULL || x2 == NULL) return NULL;
    return ufunc_apply_binary(x1, x2, DTYPE_INT32, ufunc_loop_left_shift_i32);
}

NDArray* ufunc_right_shift(NDArray* x1, NDArray* x2) {
    if (x1 == NULL || x2 == NULL) return NULL;
    return ufunc_apply_binary(x1, x2, DTYPE_INT32, ufunc_loop_right_shift_i32);
}

/* Special math functions */
/* arctan2: output is always float */
NDArray* ufunc_arctan2(NDArray* x1, NDArray* x2) {
    if (x1 == NULL || x2 == NULL) return NULL;
    DType common_dtype = dtype_promote(x1->dtype, x2->dtype);
    DType out_dtype = dtype_is_integer(common_dtype) ? DTYPE_FLOAT64 : common_dtype;
    BinaryLoopFunc loop = ufunc_select_binary_loop(common_dtype,
        ufunc_loop_arctan2_f64, ufunc_loop_arctan2_f32, ufunc_loop_arctan2_i32);
    return ufunc_apply_binary(x1, x2, out_dtype, loop);
}

/* hypot: output is always float */
NDArray* ufunc_hypot(NDArray* x1, NDArray* x2) {
    if (x1 == NULL || x2 == NULL) return NULL;
    DType common_dtype = dtype_promote(x1->dtype, x2->dtype);
    DType out_dtype = dtype_is_integer(common_dtype) ? DTYPE_FLOAT64 : common_dtype;
    BinaryLoopFunc loop = ufunc_select_binary_loop(common_dtype,
        ufunc_loop_hypot_f64, ufunc_loop_hypot_f32, ufunc_loop_hypot_i32);
    return ufunc_apply_binary(x1, x2, out_dtype, loop);
}

IMPL_BINARY_UFUNC(copysign)

/* logaddexp: output is always float */
NDArray* ufunc_logaddexp(NDArray* x1, NDArray* x2) {
    if (x1 == NULL || x2 == NULL) return NULL;
    DType common_dtype = dtype_promote(x1->dtype, x2->dtype);
    DType out_dtype = dtype_is_integer(common_dtype) ? DTYPE_FLOAT64 : common_dtype;
    BinaryLoopFunc loop = ufunc_select_binary_loop(common_dtype,
        ufunc_loop_logaddexp_f64, ufunc_loop_logaddexp_f32, ufunc_loop_logaddexp_i32);
    return ufunc_apply_binary(x1, x2, out_dtype, loop);
}

NDArray* ufunc_logaddexp2(NDArray* x1, NDArray* x2) {
    if (x1 == NULL || x2 == NULL) return NULL;
    DType common_dtype = dtype_promote(x1->dtype, x2->dtype);
    DType out_dtype = dtype_is_integer(common_dtype) ? DTYPE_FLOAT64 : common_dtype;
    BinaryLoopFunc loop = ufunc_select_binary_loop(common_dtype,
        ufunc_loop_logaddexp2_f64, ufunc_loop_logaddexp2_f32, ufunc_loop_logaddexp2_i32);
    return ufunc_apply_binary(x1, x2, out_dtype, loop);
}

/* signbit: unary function returning bool */
NDArray* ufunc_signbit(NDArray* x) {
    if (x == NULL) return NULL;

    NDArray* output = ndarray_empty(x->ndim, x->shape, DTYPE_BOOL);
    if (output == NULL) return NULL;

    size_t n = x->size;
    size_t elem_size = dtype_size(x->dtype);

    for (size_t i = 0; i < n; i++) {
        size_t offset = 0;
        size_t out_offset = 0;

        if (x->ndim > 0) {
            /* Compute offset for strided access */
            size_t idx = i;
            for (int32_t d = x->ndim - 1; d >= 0; d--) {
                size_t coord = idx % x->shape[d];
                offset += coord * x->strides[d];
                out_offset += coord * output->strides[d];
                idx /= x->shape[d];
            }
        }

        uint8_t* out_ptr = (uint8_t*)output->data + out_offset;

        switch (x->dtype) {
            case DTYPE_FLOAT64: {
                double val = *(const double*)((const char*)x->data + offset);
                *out_ptr = signbit(val) ? 1 : 0;
                break;
            }
            case DTYPE_FLOAT32: {
                float val = *(const float*)((const char*)x->data + offset);
                *out_ptr = signbit(val) ? 1 : 0;
                break;
            }
            case DTYPE_INT32: {
                int32_t val = *(const int32_t*)((const char*)x->data + offset);
                *out_ptr = (val < 0) ? 1 : 0;
                break;
            }
            default: {
                double val = ndarray_get_flat(x, i);
                *out_ptr = (val < 0) ? 1 : 0;
                break;
            }
        }
    }

    return output;
}
