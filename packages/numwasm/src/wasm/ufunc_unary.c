/**
 * NumJS-WASM Unary Universal Functions Implementation
 *
 * Implements ~30 unary element-wise operations using type-specific loops.
 */

#include "ufunc_unary.h"
#include "ufunc.h"
#include "dtype.h"
#include <math.h>
#include <stdint.h>

/* ============ Constants ============ */

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define DEG_TO_RAD (M_PI / 180.0)
#define RAD_TO_DEG (180.0 / M_PI)

/* ============ Arithmetic Loops ============ */

/* negative: -x */
UNARY_LOOP_F64(negative, -x)
UNARY_LOOP_F32(negative, -x)
UNARY_LOOP_I32(negative, -x)

/* positive: +x (identity/copy) */
UNARY_LOOP_F64(positive, x)
UNARY_LOOP_F32(positive, x)
UNARY_LOOP_I32(positive, x)

/* absolute: |x| */
UNARY_LOOP_F64(absolute, fabs(x))
UNARY_LOOP_F32(absolute, fabsf(x))
UNARY_LOOP_I32(absolute, (x < 0 ? -x : x))

/* sign: -1, 0, or +1 */
UNARY_LOOP_F64(sign, (x > 0) ? 1.0 : ((x < 0) ? -1.0 : (x == 0 ? 0.0 : x)))  /* NaN -> NaN */
UNARY_LOOP_F32(sign, (x > 0) ? 1.0f : ((x < 0) ? -1.0f : (x == 0 ? 0.0f : x)))
UNARY_LOOP_I32(sign, (x > 0) ? 1 : ((x < 0) ? -1 : 0))

/* ============ Powers and Roots Loops ============ */

/* sqrt */
UNARY_LOOP_F64(sqrt, sqrt(x))
UNARY_LOOP_F32(sqrt, sqrtf(x))
UNARY_LOOP_I32(sqrt, (int32_t)sqrt((double)x))

/* square: x^2 */
UNARY_LOOP_F64(square, x * x)
UNARY_LOOP_F32(square, x * x)
UNARY_LOOP_I32(square, x * x)

/* cbrt: cube root */
UNARY_LOOP_F64(cbrt, cbrt(x))
UNARY_LOOP_F32(cbrt, cbrtf(x))
UNARY_LOOP_I32(cbrt, (int32_t)cbrt((double)x))

/* reciprocal: 1/x */
UNARY_LOOP_F64(reciprocal, 1.0 / x)
UNARY_LOOP_F32(reciprocal, 1.0f / x)
UNARY_LOOP_I32(reciprocal, (x != 0) ? (1 / x) : 0)

/* ============ Exponential Loops ============ */

/* exp: e^x */
UNARY_LOOP_F64(exp, exp(x))
UNARY_LOOP_F32(exp, expf(x))
UNARY_LOOP_I32(exp, (int32_t)exp((double)x))

/* exp2: 2^x */
UNARY_LOOP_F64(exp2, exp2(x))
UNARY_LOOP_F32(exp2, exp2f(x))
UNARY_LOOP_I32(exp2, (int32_t)exp2((double)x))

/* expm1: e^x - 1 */
UNARY_LOOP_F64(expm1, expm1(x))
UNARY_LOOP_F32(expm1, expm1f(x))
UNARY_LOOP_I32(expm1, (int32_t)expm1((double)x))

/* ============ Logarithmic Loops ============ */

/* log: natural log */
UNARY_LOOP_F64(log, log(x))
UNARY_LOOP_F32(log, logf(x))
UNARY_LOOP_I32(log, (int32_t)log((double)x))

/* log2 */
UNARY_LOOP_F64(log2, log2(x))
UNARY_LOOP_F32(log2, log2f(x))
UNARY_LOOP_I32(log2, (int32_t)log2((double)x))

/* log10 */
UNARY_LOOP_F64(log10, log10(x))
UNARY_LOOP_F32(log10, log10f(x))
UNARY_LOOP_I32(log10, (int32_t)log10((double)x))

/* log1p: log(1 + x) */
UNARY_LOOP_F64(log1p, log1p(x))
UNARY_LOOP_F32(log1p, log1pf(x))
UNARY_LOOP_I32(log1p, (int32_t)log1p((double)x))

/* ============ Trigonometric Loops ============ */

/* sin */
UNARY_LOOP_F64(sin, sin(x))
UNARY_LOOP_F32(sin, sinf(x))
UNARY_LOOP_I32(sin, (int32_t)sin((double)x))

/* cos */
UNARY_LOOP_F64(cos, cos(x))
UNARY_LOOP_F32(cos, cosf(x))
UNARY_LOOP_I32(cos, (int32_t)cos((double)x))

/* tan */
UNARY_LOOP_F64(tan, tan(x))
UNARY_LOOP_F32(tan, tanf(x))
UNARY_LOOP_I32(tan, (int32_t)tan((double)x))

/* arcsin (asin) */
UNARY_LOOP_F64(arcsin, asin(x))
UNARY_LOOP_F32(arcsin, asinf(x))
UNARY_LOOP_I32(arcsin, (int32_t)asin((double)x))

/* arccos (acos) */
UNARY_LOOP_F64(arccos, acos(x))
UNARY_LOOP_F32(arccos, acosf(x))
UNARY_LOOP_I32(arccos, (int32_t)acos((double)x))

/* arctan (atan) */
UNARY_LOOP_F64(arctan, atan(x))
UNARY_LOOP_F32(arctan, atanf(x))
UNARY_LOOP_I32(arctan, (int32_t)atan((double)x))

/* ============ Hyperbolic Loops ============ */

/* sinh */
UNARY_LOOP_F64(sinh, sinh(x))
UNARY_LOOP_F32(sinh, sinhf(x))
UNARY_LOOP_I32(sinh, (int32_t)sinh((double)x))

/* cosh */
UNARY_LOOP_F64(cosh, cosh(x))
UNARY_LOOP_F32(cosh, coshf(x))
UNARY_LOOP_I32(cosh, (int32_t)cosh((double)x))

/* tanh */
UNARY_LOOP_F64(tanh, tanh(x))
UNARY_LOOP_F32(tanh, tanhf(x))
UNARY_LOOP_I32(tanh, (int32_t)tanh((double)x))

/* arcsinh (asinh) */
UNARY_LOOP_F64(arcsinh, asinh(x))
UNARY_LOOP_F32(arcsinh, asinhf(x))
UNARY_LOOP_I32(arcsinh, (int32_t)asinh((double)x))

/* arccosh (acosh) */
UNARY_LOOP_F64(arccosh, acosh(x))
UNARY_LOOP_F32(arccosh, acoshf(x))
UNARY_LOOP_I32(arccosh, (int32_t)acosh((double)x))

/* arctanh (atanh) */
UNARY_LOOP_F64(arctanh, atanh(x))
UNARY_LOOP_F32(arctanh, atanhf(x))
UNARY_LOOP_I32(arctanh, (int32_t)atanh((double)x))

/* ============ Rounding Loops ============ */

/* floor */
UNARY_LOOP_F64(floor, floor(x))
UNARY_LOOP_F32(floor, floorf(x))
UNARY_LOOP_I32(floor, x)  /* no-op for integers */

/* ceil */
UNARY_LOOP_F64(ceil, ceil(x))
UNARY_LOOP_F32(ceil, ceilf(x))
UNARY_LOOP_I32(ceil, x)

/* trunc */
UNARY_LOOP_F64(trunc, trunc(x))
UNARY_LOOP_F32(trunc, truncf(x))
UNARY_LOOP_I32(trunc, x)

/* rint: round to nearest integer */
UNARY_LOOP_F64(rint, rint(x))
UNARY_LOOP_F32(rint, rintf(x))
UNARY_LOOP_I32(rint, x)

/* round: round to nearest, with ties going to nearest even */
UNARY_LOOP_F64(round, round(x))
UNARY_LOOP_F32(round, roundf(x))
UNARY_LOOP_I32(round, x)

/* ============ Angle Conversion Loops ============ */

/* degrees: radians to degrees */
UNARY_LOOP_F64(degrees, x * RAD_TO_DEG)
UNARY_LOOP_F32(degrees, x * (float)RAD_TO_DEG)
UNARY_LOOP_I32(degrees, (int32_t)(x * RAD_TO_DEG))

/* radians: degrees to radians */
UNARY_LOOP_F64(radians, x * DEG_TO_RAD)
UNARY_LOOP_F32(radians, x * (float)DEG_TO_RAD)
UNARY_LOOP_I32(radians, (int32_t)(x * DEG_TO_RAD))

/* ============ Logical Loops ============ */

/* logical_not: !x (output is uint8) */
static void ufunc_loop_logical_not_f64(const char* in, char* out, size_t n,
                                        int64_t in_stride, int64_t out_stride) {
    for (size_t i = 0; i < n; i++) {
        double x = *(const double*)(in + i * in_stride);
        *(uint8_t*)(out + i * out_stride) = (x == 0.0) ? 1 : 0;
    }
}

static void ufunc_loop_logical_not_f32(const char* in, char* out, size_t n,
                                        int64_t in_stride, int64_t out_stride) {
    for (size_t i = 0; i < n; i++) {
        float x = *(const float*)(in + i * in_stride);
        *(uint8_t*)(out + i * out_stride) = (x == 0.0f) ? 1 : 0;
    }
}

static void ufunc_loop_logical_not_i32(const char* in, char* out, size_t n,
                                        int64_t in_stride, int64_t out_stride) {
    for (size_t i = 0; i < n; i++) {
        int32_t x = *(const int32_t*)(in + i * in_stride);
        *(uint8_t*)(out + i * out_stride) = (x == 0) ? 1 : 0;
    }
}

/* ============ Bitwise Loops (Integer only) ============ */

static void ufunc_loop_invert_i32(const char* in, char* out, size_t n,
                                   int64_t in_stride, int64_t out_stride) {
    for (size_t i = 0; i < n; i++) {
        int32_t x = *(const int32_t*)(in + i * in_stride);
        *(int32_t*)(out + i * out_stride) = ~x;
    }
}

/* ============ Exported Functions ============ */

/* Helper macro for implementing unary functions */
#define IMPL_UNARY_UFUNC(name) \
NDArray* ufunc_##name(NDArray* input) { \
    if (input == NULL) return NULL; \
    UnaryLoopFunc loop = ufunc_select_unary_loop( \
        input->dtype, \
        ufunc_loop_##name##_f64, \
        ufunc_loop_##name##_f32, \
        ufunc_loop_##name##_i32); \
    return ufunc_apply_unary(input, input->dtype, loop); \
}

/* Arithmetic */
IMPL_UNARY_UFUNC(negative)
IMPL_UNARY_UFUNC(positive)
IMPL_UNARY_UFUNC(absolute)
IMPL_UNARY_UFUNC(sign)

NDArray* ufunc_abs(NDArray* input) { return ufunc_absolute(input); }

/* Powers and Roots */
IMPL_UNARY_UFUNC(sqrt)
IMPL_UNARY_UFUNC(square)
IMPL_UNARY_UFUNC(cbrt)
IMPL_UNARY_UFUNC(reciprocal)

/* Exponential */
IMPL_UNARY_UFUNC(exp)
IMPL_UNARY_UFUNC(exp2)
IMPL_UNARY_UFUNC(expm1)

/* Logarithmic */
IMPL_UNARY_UFUNC(log)
IMPL_UNARY_UFUNC(log2)
IMPL_UNARY_UFUNC(log10)
IMPL_UNARY_UFUNC(log1p)

/* Trigonometric */
IMPL_UNARY_UFUNC(sin)
IMPL_UNARY_UFUNC(cos)
IMPL_UNARY_UFUNC(tan)
IMPL_UNARY_UFUNC(arcsin)
IMPL_UNARY_UFUNC(arccos)
IMPL_UNARY_UFUNC(arctan)

/* Hyperbolic */
IMPL_UNARY_UFUNC(sinh)
IMPL_UNARY_UFUNC(cosh)
IMPL_UNARY_UFUNC(tanh)
IMPL_UNARY_UFUNC(arcsinh)
IMPL_UNARY_UFUNC(arccosh)
IMPL_UNARY_UFUNC(arctanh)

/* Rounding */
IMPL_UNARY_UFUNC(floor)
IMPL_UNARY_UFUNC(ceil)
IMPL_UNARY_UFUNC(trunc)
IMPL_UNARY_UFUNC(rint)
IMPL_UNARY_UFUNC(round)

/* Angle conversion */
IMPL_UNARY_UFUNC(degrees)
IMPL_UNARY_UFUNC(radians)

NDArray* ufunc_rad2deg(NDArray* input) { return ufunc_degrees(input); }
NDArray* ufunc_deg2rad(NDArray* input) { return ufunc_radians(input); }

/* Logical not (output is bool) */
NDArray* ufunc_logical_not(NDArray* input) {
    if (input == NULL) return NULL;
    UnaryLoopFunc loop = ufunc_select_unary_loop(
        input->dtype,
        ufunc_loop_logical_not_f64,
        ufunc_loop_logical_not_f32,
        ufunc_loop_logical_not_i32);
    return ufunc_apply_unary(input, DTYPE_BOOL, loop);
}

/* Bitwise invert (integer only) */
NDArray* ufunc_invert(NDArray* input) {
    if (input == NULL) return NULL;
    if (!dtype_is_integer(input->dtype) && input->dtype != DTYPE_BOOL) {
        /* For non-integer types, convert to int32, apply, stay int32 */
        NDArray* as_int = ndarray_astype(input, DTYPE_INT32);
        if (as_int == NULL) return NULL;
        NDArray* result = ufunc_apply_unary(as_int, DTYPE_INT32, ufunc_loop_invert_i32);
        ndarray_free(as_int);
        return result;
    }
    return ufunc_apply_unary(input, input->dtype, ufunc_loop_invert_i32);
}

NDArray* ufunc_bitwise_not(NDArray* input) { return ufunc_invert(input); }
