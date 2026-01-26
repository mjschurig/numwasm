/**
 * NumJS-WASM Universal Functions (Ufuncs) Infrastructure
 *
 * Provides element-wise operations on NDArrays with broadcasting support.
 * Inspired by NumPy's ufunc implementation in numpy/_core/src/umath/ufunc_object.c
 *
 * Design decisions:
 * - Direct function exports (not registry pattern) for simplicity and tree-shaking
 * - Type-specific loops via macros for f64, f32, i32
 * - Optimized strided iteration with pre-computed byte strides
 * - C-side output allocation (TypeScript manages lifecycle)
 */

#ifndef NUMJS_UFUNC_H
#define NUMJS_UFUNC_H

#include "ndarray.h"
#include "broadcast.h"
#include <stdint.h>

/* ============ Identity Values for Reductions ============ */

#define UFUNC_IDENTITY_ZERO   0
#define UFUNC_IDENTITY_ONE    1
#define UFUNC_IDENTITY_NONE  -1  /* No identity (e.g., for max/min) */

/* ============ Loop Function Types ============ */

/**
 * Unary loop function signature.
 * Processes n elements from input to output with given strides.
 *
 * @param in         Input data pointer
 * @param out        Output data pointer
 * @param n          Number of elements to process
 * @param in_stride  Stride for input (in bytes)
 * @param out_stride Stride for output (in bytes)
 */
typedef void (*UnaryLoopFunc)(const char* in, char* out, size_t n,
                              int64_t in_stride, int64_t out_stride);

/**
 * Binary loop function signature.
 * Processes n elements from two inputs to output with given strides.
 *
 * @param in1        First input data pointer
 * @param in2        Second input data pointer
 * @param out        Output data pointer
 * @param n          Number of elements to process
 * @param stride1    Stride for first input (in bytes)
 * @param stride2    Stride for second input (in bytes)
 * @param out_stride Stride for output (in bytes)
 */
typedef void (*BinaryLoopFunc)(const char* in1, const char* in2, char* out,
                               size_t n, int64_t stride1, int64_t stride2,
                               int64_t out_stride);

/* ============ Unary Loop Macros ============ */

/**
 * Generate a unary loop for float64.
 * Usage: UNARY_LOOP_F64(sqrt, sqrt(x))
 */
#define UNARY_LOOP_F64(name, op) \
static void ufunc_loop_##name##_f64(const char* in, char* out, size_t n, \
                                    int64_t in_stride, int64_t out_stride) { \
    for (size_t i = 0; i < n; i++) { \
        double x = *(const double*)(in + i * in_stride); \
        *(double*)(out + i * out_stride) = op; \
    } \
}

/**
 * Generate a unary loop for float32.
 */
#define UNARY_LOOP_F32(name, op) \
static void ufunc_loop_##name##_f32(const char* in, char* out, size_t n, \
                                    int64_t in_stride, int64_t out_stride) { \
    for (size_t i = 0; i < n; i++) { \
        float x = *(const float*)(in + i * in_stride); \
        *(float*)(out + i * out_stride) = (float)(op); \
    } \
}

/**
 * Generate a unary loop for int32.
 */
#define UNARY_LOOP_I32(name, op) \
static void ufunc_loop_##name##_i32(const char* in, char* out, size_t n, \
                                    int64_t in_stride, int64_t out_stride) { \
    for (size_t i = 0; i < n; i++) { \
        int32_t x = *(const int32_t*)(in + i * in_stride); \
        *(int32_t*)(out + i * out_stride) = (int32_t)(op); \
    } \
}

/**
 * Generate all numeric type loops for a unary operation.
 * Creates f64, f32, and i32 variants.
 */
#define UNARY_LOOPS_NUMERIC(name, op) \
    UNARY_LOOP_F64(name, op) \
    UNARY_LOOP_F32(name, op) \
    UNARY_LOOP_I32(name, op)

/* ============ Binary Loop Macros ============ */

/**
 * Generate a binary loop for float64.
 * Usage: BINARY_LOOP_F64(add, x + y)
 */
#define BINARY_LOOP_F64(name, op) \
static void ufunc_loop_##name##_f64(const char* in1, const char* in2, char* out, \
                                    size_t n, int64_t s1, int64_t s2, int64_t so) { \
    for (size_t i = 0; i < n; i++) { \
        double x = *(const double*)(in1 + i * s1); \
        double y = *(const double*)(in2 + i * s2); \
        *(double*)(out + i * so) = op; \
    } \
}

/**
 * Generate a binary loop for float32.
 */
#define BINARY_LOOP_F32(name, op) \
static void ufunc_loop_##name##_f32(const char* in1, const char* in2, char* out, \
                                    size_t n, int64_t s1, int64_t s2, int64_t so) { \
    for (size_t i = 0; i < n; i++) { \
        float x = *(const float*)(in1 + i * s1); \
        float y = *(const float*)(in2 + i * s2); \
        *(float*)(out + i * so) = (float)(op); \
    } \
}

/**
 * Generate a binary loop for int32.
 */
#define BINARY_LOOP_I32(name, op) \
static void ufunc_loop_##name##_i32(const char* in1, const char* in2, char* out, \
                                    size_t n, int64_t s1, int64_t s2, int64_t so) { \
    for (size_t i = 0; i < n; i++) { \
        int32_t x = *(const int32_t*)(in1 + i * s1); \
        int32_t y = *(const int32_t*)(in2 + i * s2); \
        *(int32_t*)(out + i * so) = (int32_t)(op); \
    } \
}

/**
 * Generate all numeric type loops for a binary operation.
 */
#define BINARY_LOOPS_NUMERIC(name, op) \
    BINARY_LOOP_F64(name, op) \
    BINARY_LOOP_F32(name, op) \
    BINARY_LOOP_I32(name, op)

/**
 * Generate a binary comparison loop (output is always uint8/bool).
 */
#define BINARY_CMP_LOOP_F64(name, op) \
static void ufunc_loop_##name##_f64(const char* in1, const char* in2, char* out, \
                                    size_t n, int64_t s1, int64_t s2, int64_t so) { \
    for (size_t i = 0; i < n; i++) { \
        double x = *(const double*)(in1 + i * s1); \
        double y = *(const double*)(in2 + i * s2); \
        *(uint8_t*)(out + i * so) = (op) ? 1 : 0; \
    } \
}

#define BINARY_CMP_LOOP_F32(name, op) \
static void ufunc_loop_##name##_f32(const char* in1, const char* in2, char* out, \
                                    size_t n, int64_t s1, int64_t s2, int64_t so) { \
    for (size_t i = 0; i < n; i++) { \
        float x = *(const float*)(in1 + i * s1); \
        float y = *(const float*)(in2 + i * s2); \
        *(uint8_t*)(out + i * so) = (op) ? 1 : 0; \
    } \
}

#define BINARY_CMP_LOOP_I32(name, op) \
static void ufunc_loop_##name##_i32(const char* in1, const char* in2, char* out, \
                                    size_t n, int64_t s1, int64_t s2, int64_t so) { \
    for (size_t i = 0; i < n; i++) { \
        int32_t x = *(const int32_t*)(in1 + i * s1); \
        int32_t y = *(const int32_t*)(in2 + i * s2); \
        *(uint8_t*)(out + i * so) = (op) ? 1 : 0; \
    } \
}

#define BINARY_CMP_LOOPS(name, op) \
    BINARY_CMP_LOOP_F64(name, op) \
    BINARY_CMP_LOOP_F32(name, op) \
    BINARY_CMP_LOOP_I32(name, op)

/* ============ High-Level Apply Functions ============ */

/**
 * Apply a unary loop function to an array.
 *
 * Handles:
 * - Contiguous fast path (direct memory scan)
 * - Non-contiguous iteration via strides
 * - Output allocation with same shape
 *
 * @param input    Input array
 * @param out_dtype Output dtype (use input->dtype for same type)
 * @param loop     Loop function to apply
 * @return         New array with result, or NULL on error
 */
NDArray* ufunc_apply_unary(NDArray* input, DType out_dtype, UnaryLoopFunc loop);

/**
 * Apply a binary loop function to two arrays.
 *
 * Handles:
 * - Broadcasting to common shape
 * - Stride calculation for broadcast dimensions (stride=0)
 * - Contiguous fast path when possible
 * - Output allocation with broadcast shape
 *
 * @param in1       First input array
 * @param in2       Second input array
 * @param out_dtype Output dtype
 * @param loop      Loop function to apply
 * @return          New array with result, or NULL on error
 */
NDArray* ufunc_apply_binary(NDArray* in1, NDArray* in2, DType out_dtype, BinaryLoopFunc loop);

/**
 * Apply a binary comparison loop (output dtype is always BOOL).
 *
 * @param in1   First input array
 * @param in2   Second input array
 * @param loop  Comparison loop function
 * @return      Boolean array with result, or NULL on error
 */
NDArray* ufunc_apply_binary_cmp(NDArray* in1, NDArray* in2, BinaryLoopFunc loop);

/* ============ Strided Iterator Helpers ============ */

/**
 * Get contiguous inner loop size.
 * Returns the number of elements that can be processed contiguously
 * from the innermost dimension.
 *
 * @param arr   Array to check
 * @return      Number of contiguous elements (at least 1)
 */
size_t ufunc_get_inner_loop_size(NDArray* arr);

/**
 * Check if array can use fast contiguous path.
 *
 * @param arr   Array to check
 * @return      true if C-contiguous
 */
bool ufunc_is_contiguous(NDArray* arr);

/**
 * Check if two arrays are both contiguous and have same shape.
 * Used to enable fast binary loop path.
 */
bool ufunc_binary_contiguous(NDArray* arr1, NDArray* arr2);

/* ============ Type Selection ============ */

/**
 * Select the appropriate loop function based on dtype.
 * Returns NULL if dtype is not supported.
 */
UnaryLoopFunc ufunc_select_unary_loop(DType dtype,
                                       UnaryLoopFunc f64_loop,
                                       UnaryLoopFunc f32_loop,
                                       UnaryLoopFunc i32_loop);

BinaryLoopFunc ufunc_select_binary_loop(DType dtype,
                                         BinaryLoopFunc f64_loop,
                                         BinaryLoopFunc f32_loop,
                                         BinaryLoopFunc i32_loop);

/**
 * Determine output dtype for type promotion.
 * For binary operations, promotes to common type.
 */
DType ufunc_result_dtype(DType dtype1, DType dtype2);

#endif /* NUMJS_UFUNC_H */
