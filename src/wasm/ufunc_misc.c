/**
 * NumJS-WASM Miscellaneous Universal Functions Implementation
 * Phase 26: frexp, ldexp, nextafter, modf, gcd, lcm, heaviside, divmod, bitwise_count
 */

#include "ufunc.h"
#include "ndarray.h"
#include "dtype.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

/* ============ Tuple Result Structure ============ */

/**
 * Structure for functions returning two arrays (frexp, modf, divmod).
 */
typedef struct {
    NDArray* out1;
    NDArray* out2;
} TupleResult;

/* ============ ldexp Loops: x1 * 2^x2 ============ */
/* NOTE: ldexp is special - x1 is float, x2 is always integer exponent */

static void ufunc_loop_ldexp_f64_i32(const char* in1, const char* in2, char* out,
                                      size_t n, int64_t s1, int64_t s2, int64_t so) {
    for (size_t i = 0; i < n; i++) {
        double x = *(const double*)(in1 + i * s1);
        int32_t exp = *(const int32_t*)(in2 + i * s2);
        *(double*)(out + i * so) = ldexp(x, exp);
    }
}

static void ufunc_loop_ldexp_f32_i32(const char* in1, const char* in2, char* out,
                                      size_t n, int64_t s1, int64_t s2, int64_t so) {
    for (size_t i = 0; i < n; i++) {
        float x = *(const float*)(in1 + i * s1);
        int32_t exp = *(const int32_t*)(in2 + i * s2);
        *(float*)(out + i * so) = ldexpf(x, exp);
    }
}

static void ufunc_loop_ldexp_f64_f64(const char* in1, const char* in2, char* out,
                                      size_t n, int64_t s1, int64_t s2, int64_t so) {
    for (size_t i = 0; i < n; i++) {
        double x = *(const double*)(in1 + i * s1);
        double exp_f = *(const double*)(in2 + i * s2);
        *(double*)(out + i * so) = ldexp(x, (int)exp_f);
    }
}

static void ufunc_loop_ldexp_i32_i32(const char* in1, const char* in2, char* out,
                                      size_t n, int64_t s1, int64_t s2, int64_t so) {
    for (size_t i = 0; i < n; i++) {
        int32_t x = *(const int32_t*)(in1 + i * s1);
        int32_t exp = *(const int32_t*)(in2 + i * s2);
        *(double*)(out + i * so) = ldexp((double)x, exp);
    }
}

/* ============ nextafter Loops ============ */

static void ufunc_loop_nextafter_f64(const char* in1, const char* in2, char* out,
                                      size_t n, int64_t s1, int64_t s2, int64_t so) {
    for (size_t i = 0; i < n; i++) {
        double x = *(const double*)(in1 + i * s1);
        double y = *(const double*)(in2 + i * s2);
        *(double*)(out + i * so) = nextafter(x, y);
    }
}

static void ufunc_loop_nextafter_f32(const char* in1, const char* in2, char* out,
                                      size_t n, int64_t s1, int64_t s2, int64_t so) {
    for (size_t i = 0; i < n; i++) {
        float x = *(const float*)(in1 + i * s1);
        float y = *(const float*)(in2 + i * s2);
        *(float*)(out + i * so) = nextafterf(x, y);
    }
}

static void ufunc_loop_nextafter_i32(const char* in1, const char* in2, char* out,
                                      size_t n, int64_t s1, int64_t s2, int64_t so) {
    for (size_t i = 0; i < n; i++) {
        int32_t x = *(const int32_t*)(in1 + i * s1);
        int32_t y = *(const int32_t*)(in2 + i * s2);
        *(double*)(out + i * so) = nextafter((double)x, (double)y);
    }
}

/* ============ heaviside Loops ============ */

static void ufunc_loop_heaviside_f64(const char* in1, const char* in2, char* out,
                                      size_t n, int64_t s1, int64_t s2, int64_t so) {
    for (size_t i = 0; i < n; i++) {
        double x = *(const double*)(in1 + i * s1);
        double h0 = *(const double*)(in2 + i * s2);
        double result;
        if (isnan(x)) {
            result = NAN;
        } else if (x < 0) {
            result = 0.0;
        } else if (x == 0) {
            result = h0;
        } else {
            result = 1.0;
        }
        *(double*)(out + i * so) = result;
    }
}

static void ufunc_loop_heaviside_f32(const char* in1, const char* in2, char* out,
                                      size_t n, int64_t s1, int64_t s2, int64_t so) {
    for (size_t i = 0; i < n; i++) {
        float x = *(const float*)(in1 + i * s1);
        float h0 = *(const float*)(in2 + i * s2);
        float result;
        if (isnan(x)) {
            result = NAN;
        } else if (x < 0) {
            result = 0.0f;
        } else if (x == 0) {
            result = h0;
        } else {
            result = 1.0f;
        }
        *(float*)(out + i * so) = result;
    }
}

static void ufunc_loop_heaviside_i32(const char* in1, const char* in2, char* out,
                                      size_t n, int64_t s1, int64_t s2, int64_t so) {
    for (size_t i = 0; i < n; i++) {
        int32_t x = *(const int32_t*)(in1 + i * s1);
        int32_t h0 = *(const int32_t*)(in2 + i * s2);
        double result;
        if (x < 0) {
            result = 0.0;
        } else if (x == 0) {
            result = (double)h0;
        } else {
            result = 1.0;
        }
        *(double*)(out + i * so) = result;
    }
}

/* ============ GCD Scalar Implementation ============ */

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

/* ============ GCD Loops ============ */

static void ufunc_loop_gcd_i32(const char* in1, const char* in2, char* out,
                                size_t n, int64_t s1, int64_t s2, int64_t so) {
    for (size_t i = 0; i < n; i++) {
        int32_t x = *(const int32_t*)(in1 + i * s1);
        int32_t y = *(const int32_t*)(in2 + i * s2);
        *(int32_t*)(out + i * so) = gcd_scalar_i32(x, y);
    }
}

static void ufunc_loop_gcd_i64(const char* in1, const char* in2, char* out,
                                size_t n, int64_t s1, int64_t s2, int64_t so) {
    for (size_t i = 0; i < n; i++) {
        int64_t x = *(const int64_t*)(in1 + i * s1);
        int64_t y = *(const int64_t*)(in2 + i * s2);
        *(int64_t*)(out + i * so) = gcd_scalar_i64(x, y);
    }
}

/* ============ LCM Loops ============ */

static void ufunc_loop_lcm_i32(const char* in1, const char* in2, char* out,
                                size_t n, int64_t s1, int64_t s2, int64_t so) {
    for (size_t i = 0; i < n; i++) {
        int32_t a = *(const int32_t*)(in1 + i * s1);
        int32_t b = *(const int32_t*)(in2 + i * s2);
        a = a < 0 ? -a : a;
        b = b < 0 ? -b : b;
        if (a == 0 || b == 0) {
            *(int32_t*)(out + i * so) = 0;
        } else {
            *(int32_t*)(out + i * so) = (a / gcd_scalar_i32(a, b)) * b;
        }
    }
}

static void ufunc_loop_lcm_i64(const char* in1, const char* in2, char* out,
                                size_t n, int64_t s1, int64_t s2, int64_t so) {
    for (size_t i = 0; i < n; i++) {
        int64_t a = *(const int64_t*)(in1 + i * s1);
        int64_t b = *(const int64_t*)(in2 + i * s2);
        a = a < 0 ? -a : a;
        b = b < 0 ? -b : b;
        if (a == 0 || b == 0) {
            *(int64_t*)(out + i * so) = 0;
        } else {
            *(int64_t*)(out + i * so) = (a / gcd_scalar_i64(a, b)) * b;
        }
    }
}

/* ============ Bitwise Count (Popcount) ============ */

static uint8_t popcount32(uint32_t v) {
    uint8_t count = 0;
    while (v) {
        v &= v - 1;  /* Brian Kernighan's algorithm */
        count++;
    }
    return count;
}

static void ufunc_loop_bitwise_count_i32(const char* in, char* out, size_t n,
                                          int64_t in_stride, int64_t out_stride) {
    for (size_t i = 0; i < n; i++) {
        int32_t x = *(const int32_t*)(in + i * in_stride);
        *(uint8_t*)(out + i * out_stride) = popcount32((uint32_t)x);
    }
}

static void ufunc_loop_bitwise_count_i8(const char* in, char* out, size_t n,
                                         int64_t in_stride, int64_t out_stride) {
    for (size_t i = 0; i < n; i++) {
        int8_t x = *(const int8_t*)(in + i * in_stride);
        *(uint8_t*)(out + i * out_stride) = popcount32((uint8_t)x);
    }
}

static void ufunc_loop_bitwise_count_u8(const char* in, char* out, size_t n,
                                         int64_t in_stride, int64_t out_stride) {
    for (size_t i = 0; i < n; i++) {
        uint8_t x = *(const uint8_t*)(in + i * in_stride);
        *(uint8_t*)(out + i * out_stride) = popcount32(x);
    }
}

/* ============ spacing Loop: distance to next representable number ============ */

static void ufunc_loop_spacing_f64(const char* in, char* out, size_t n,
                                    int64_t in_stride, int64_t out_stride) {
    for (size_t i = 0; i < n; i++) {
        double x = *(const double*)(in + i * in_stride);
        double result;
        if (!isfinite(x)) {
            result = NAN;
        } else if (x == 0.0) {
            /* Smallest positive subnormal */
            result = 5e-324;  /* ~= DBL_TRUE_MIN */
        } else {
            double ax = fabs(x);
            double next = nextafter(ax, INFINITY);
            result = next - ax;
        }
        *(double*)(out + i * out_stride) = result;
    }
}

static void ufunc_loop_spacing_f32(const char* in, char* out, size_t n,
                                    int64_t in_stride, int64_t out_stride) {
    for (size_t i = 0; i < n; i++) {
        float x = *(const float*)(in + i * in_stride);
        float result;
        if (!isfinite(x)) {
            result = NAN;
        } else if (x == 0.0f) {
            result = 1e-45f;  /* ~= FLT_TRUE_MIN */
        } else {
            float ax = fabsf(x);
            float next = nextafterf(ax, INFINITY);
            result = next - ax;
        }
        *(float*)(out + i * out_stride) = result;
    }
}

static void ufunc_loop_spacing_i32(const char* in, char* out, size_t n,
                                    int64_t in_stride, int64_t out_stride) {
    for (size_t i = 0; i < n; i++) {
        int32_t x = *(const int32_t*)(in + i * in_stride);
        double dx = (double)x;
        double result;
        if (x == 0) {
            result = 5e-324;
        } else {
            double ax = fabs(dx);
            double next = nextafter(ax, INFINITY);
            result = next - ax;
        }
        *(double*)(out + i * out_stride) = result;
    }
}

/* ============ sinc Loop: sin(pi*x) / (pi*x) ============ */

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static void ufunc_loop_sinc_f64(const char* in, char* out, size_t n,
                                 int64_t in_stride, int64_t out_stride) {
    for (size_t i = 0; i < n; i++) {
        double x = *(const double*)(in + i * in_stride);
        double result;
        if (x == 0.0) {
            result = 1.0;
        } else {
            double pix = M_PI * x;
            result = sin(pix) / pix;
        }
        *(double*)(out + i * out_stride) = result;
    }
}

static void ufunc_loop_sinc_f32(const char* in, char* out, size_t n,
                                 int64_t in_stride, int64_t out_stride) {
    for (size_t i = 0; i < n; i++) {
        float x = *(const float*)(in + i * in_stride);
        float result;
        if (x == 0.0f) {
            result = 1.0f;
        } else {
            float pix = (float)M_PI * x;
            result = sinf(pix) / pix;
        }
        *(float*)(out + i * out_stride) = result;
    }
}

static void ufunc_loop_sinc_i32(const char* in, char* out, size_t n,
                                 int64_t in_stride, int64_t out_stride) {
    for (size_t i = 0; i < n; i++) {
        int32_t x = *(const int32_t*)(in + i * in_stride);
        double result;
        if (x == 0) {
            result = 1.0;
        } else {
            double pix = M_PI * (double)x;
            result = sin(pix) / pix;
        }
        *(double*)(out + i * out_stride) = result;
    }
}

/* ============ Exported Functions ============ */

/* ldexp: x1 * 2^x2 - special case where x2 should be integer */
NDArray* ufunc_ldexp(NDArray* x1, NDArray* x2) {
    if (x1 == NULL || x2 == NULL) return NULL;

    DType out_dtype;
    BinaryLoopFunc loop;

    /* Select loop based on actual input dtypes, not promoted dtype */
    bool x1_int = dtype_is_integer(x1->dtype);
    bool x2_int = dtype_is_integer(x2->dtype);

    if (x1_int && x2_int) {
        /* Both integers: output float64 */
        out_dtype = DTYPE_FLOAT64;
        loop = ufunc_loop_ldexp_i32_i32;
    } else if (x1->dtype == DTYPE_FLOAT32 && x2_int) {
        /* Float32 mantissa, int exponent */
        out_dtype = DTYPE_FLOAT32;
        loop = ufunc_loop_ldexp_f32_i32;
    } else if ((x1->dtype == DTYPE_FLOAT64 || x1_int) && x2_int) {
        /* Float64 mantissa (or promoted int), int exponent */
        out_dtype = DTYPE_FLOAT64;
        loop = ufunc_loop_ldexp_f64_i32;
    } else {
        /* Fallback: both treated as float64 (e.g., when exponent is float) */
        out_dtype = DTYPE_FLOAT64;
        loop = ufunc_loop_ldexp_f64_f64;
    }

    return ufunc_apply_binary(x1, x2, out_dtype, loop);
}

/* nextafter: next float toward direction */
NDArray* ufunc_nextafter(NDArray* x1, NDArray* x2) {
    if (x1 == NULL || x2 == NULL) return NULL;

    DType common_dtype = dtype_promote(x1->dtype, x2->dtype);
    DType out_dtype;
    BinaryLoopFunc loop;

    if (dtype_is_integer(common_dtype)) {
        out_dtype = DTYPE_FLOAT64;
        loop = ufunc_loop_nextafter_i32;
    } else if (common_dtype == DTYPE_FLOAT32) {
        out_dtype = DTYPE_FLOAT32;
        loop = ufunc_loop_nextafter_f32;
    } else {
        out_dtype = DTYPE_FLOAT64;
        loop = ufunc_loop_nextafter_f64;
    }

    return ufunc_apply_binary(x1, x2, out_dtype, loop);
}

/* heaviside: step function */
NDArray* ufunc_heaviside(NDArray* x1, NDArray* x2) {
    if (x1 == NULL || x2 == NULL) return NULL;

    DType common_dtype = dtype_promote(x1->dtype, x2->dtype);
    DType out_dtype;
    BinaryLoopFunc loop;

    if (dtype_is_integer(common_dtype)) {
        out_dtype = DTYPE_FLOAT64;
        loop = ufunc_loop_heaviside_i32;
    } else if (common_dtype == DTYPE_FLOAT32) {
        out_dtype = DTYPE_FLOAT32;
        loop = ufunc_loop_heaviside_f32;
    } else {
        out_dtype = DTYPE_FLOAT64;
        loop = ufunc_loop_heaviside_f64;
    }

    return ufunc_apply_binary(x1, x2, out_dtype, loop);
}

/* gcd: greatest common divisor (integer only) */
NDArray* ufunc_gcd(NDArray* x1, NDArray* x2) {
    if (x1 == NULL || x2 == NULL) return NULL;

    /* Convert to int32 if needed, use i64 loop for int64 */
    DType common_dtype = dtype_promote(x1->dtype, x2->dtype);

    if (common_dtype == DTYPE_INT64 || common_dtype == DTYPE_UINT64) {
        return ufunc_apply_binary(x1, x2, DTYPE_INT64, ufunc_loop_gcd_i64);
    }

    return ufunc_apply_binary(x1, x2, DTYPE_INT32, ufunc_loop_gcd_i32);
}

/* lcm: lowest common multiple (integer only) */
NDArray* ufunc_lcm(NDArray* x1, NDArray* x2) {
    if (x1 == NULL || x2 == NULL) return NULL;

    DType common_dtype = dtype_promote(x1->dtype, x2->dtype);

    if (common_dtype == DTYPE_INT64 || common_dtype == DTYPE_UINT64) {
        return ufunc_apply_binary(x1, x2, DTYPE_INT64, ufunc_loop_lcm_i64);
    }

    return ufunc_apply_binary(x1, x2, DTYPE_INT32, ufunc_loop_lcm_i32);
}

/* spacing: distance to next representable number */
NDArray* ufunc_spacing(NDArray* input) {
    if (input == NULL) return NULL;

    DType out_dtype;
    UnaryLoopFunc loop;

    if (dtype_is_integer(input->dtype)) {
        out_dtype = DTYPE_FLOAT64;
        loop = ufunc_loop_spacing_i32;
    } else if (input->dtype == DTYPE_FLOAT32) {
        out_dtype = DTYPE_FLOAT32;
        loop = ufunc_loop_spacing_f32;
    } else {
        out_dtype = DTYPE_FLOAT64;
        loop = ufunc_loop_spacing_f64;
    }

    return ufunc_apply_unary(input, out_dtype, loop);
}

/* sinc: sin(pi*x) / (pi*x) */
NDArray* ufunc_sinc(NDArray* input) {
    if (input == NULL) return NULL;

    DType out_dtype;
    UnaryLoopFunc loop;

    if (dtype_is_integer(input->dtype)) {
        out_dtype = DTYPE_FLOAT64;
        loop = ufunc_loop_sinc_i32;
    } else if (input->dtype == DTYPE_FLOAT32) {
        out_dtype = DTYPE_FLOAT32;
        loop = ufunc_loop_sinc_f32;
    } else {
        out_dtype = DTYPE_FLOAT64;
        loop = ufunc_loop_sinc_f64;
    }

    return ufunc_apply_unary(input, out_dtype, loop);
}

/* bitwise_count: population count (integer only) */
NDArray* ufunc_bitwise_count(NDArray* input) {
    if (input == NULL) return NULL;

    UnaryLoopFunc loop;

    switch (input->dtype) {
        case DTYPE_INT8:
            loop = ufunc_loop_bitwise_count_i8;
            break;
        case DTYPE_UINT8:
        case DTYPE_BOOL:
            loop = ufunc_loop_bitwise_count_u8;
            break;
        default:
            /* All other integer types use i32 loop */
            loop = ufunc_loop_bitwise_count_i32;
            break;
    }

    return ufunc_apply_unary(input, DTYPE_UINT8, loop);
}

/* ============ Tuple-Returning Functions ============ */

/* frexp: decompose into mantissa and exponent */
TupleResult* ufunc_frexp(NDArray* input) {
    if (input == NULL) return NULL;

    TupleResult* result = (TupleResult*)malloc(sizeof(TupleResult));
    if (result == NULL) return NULL;

    /* Mantissa is same float type, exponent is always int32 */
    DType mantissa_dtype = dtype_is_floating(input->dtype) ? input->dtype : DTYPE_FLOAT64;

    result->out1 = ndarray_empty(input->ndim, input->shape, mantissa_dtype);
    result->out2 = ndarray_empty(input->ndim, input->shape, DTYPE_INT32);

    if (result->out1 == NULL || result->out2 == NULL) {
        if (result->out1) ndarray_free(result->out1);
        if (result->out2) ndarray_free(result->out2);
        free(result);
        return NULL;
    }

    size_t n = input->size;

    if (input->dtype == DTYPE_FLOAT32) {
        for (size_t i = 0; i < n; i++) {
            float x = ((float*)input->data)[i];
            int exp;
            float mantissa = frexpf(x, &exp);
            ((float*)result->out1->data)[i] = mantissa;
            ((int32_t*)result->out2->data)[i] = exp;
        }
    } else {
        /* Float64 or convert from integer */
        for (size_t i = 0; i < n; i++) {
            double x;
            if (dtype_is_integer(input->dtype)) {
                x = (double)ndarray_get_flat(input, i);
            } else {
                x = ((double*)input->data)[i];
            }
            int exp;
            double mantissa = frexp(x, &exp);
            ((double*)result->out1->data)[i] = mantissa;
            ((int32_t*)result->out2->data)[i] = exp;
        }
    }

    return result;
}

/* modf: separate into fractional and integral parts */
TupleResult* ufunc_modf(NDArray* input) {
    if (input == NULL) return NULL;

    TupleResult* result = (TupleResult*)malloc(sizeof(TupleResult));
    if (result == NULL) return NULL;

    DType out_dtype = dtype_is_floating(input->dtype) ? input->dtype : DTYPE_FLOAT64;

    result->out1 = ndarray_empty(input->ndim, input->shape, out_dtype);  /* fractional */
    result->out2 = ndarray_empty(input->ndim, input->shape, out_dtype);  /* integral */

    if (result->out1 == NULL || result->out2 == NULL) {
        if (result->out1) ndarray_free(result->out1);
        if (result->out2) ndarray_free(result->out2);
        free(result);
        return NULL;
    }

    size_t n = input->size;

    if (input->dtype == DTYPE_FLOAT32) {
        for (size_t i = 0; i < n; i++) {
            float x = ((float*)input->data)[i];
            float intpart;
            float fracpart = modff(x, &intpart);
            ((float*)result->out1->data)[i] = fracpart;
            ((float*)result->out2->data)[i] = intpart;
        }
    } else {
        for (size_t i = 0; i < n; i++) {
            double x;
            if (dtype_is_integer(input->dtype)) {
                x = (double)ndarray_get_flat(input, i);
            } else {
                x = ((double*)input->data)[i];
            }
            double intpart;
            double fracpart = modf(x, &intpart);
            ((double*)result->out1->data)[i] = fracpart;
            ((double*)result->out2->data)[i] = intpart;
        }
    }

    return result;
}

/* Forward declarations for floor_divide and remainder from ufunc_binary.c */
extern NDArray* ufunc_floor_divide(NDArray* x1, NDArray* x2);
extern NDArray* ufunc_remainder(NDArray* x1, NDArray* x2);

/* divmod: quotient and remainder */
TupleResult* ufunc_divmod(NDArray* x1, NDArray* x2) {
    if (x1 == NULL || x2 == NULL) return NULL;

    TupleResult* result = (TupleResult*)malloc(sizeof(TupleResult));
    if (result == NULL) return NULL;

    result->out1 = ufunc_floor_divide(x1, x2);
    result->out2 = ufunc_remainder(x1, x2);

    if (result->out1 == NULL || result->out2 == NULL) {
        if (result->out1) ndarray_free(result->out1);
        if (result->out2) ndarray_free(result->out2);
        free(result);
        return NULL;
    }

    return result;
}

/* ============ Tuple Result Helpers ============ */

void ufunc_tuple_result_free(TupleResult* result) {
    if (result != NULL) {
        /* Don't free the arrays - TypeScript manages their lifecycle */
        free(result);
    }
}

NDArray* ufunc_tuple_get_first(TupleResult* result) {
    return result ? result->out1 : NULL;
}

NDArray* ufunc_tuple_get_second(TupleResult* result) {
    return result ? result->out2 : NULL;
}
