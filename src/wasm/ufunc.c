/**
 * NumJS-WASM Universal Functions Implementation
 *
 * Provides high-level apply functions for ufuncs that handle:
 * - Broadcasting for binary operations
 * - Strided iteration for non-contiguous arrays
 * - Fast path for contiguous arrays
 * - Output array allocation
 */

#include "ufunc.h"
#include "ndarray.h"
#include "broadcast.h"
#include "dtype.h"
#include <stdlib.h>
#include <string.h>

/* ============ Strided Iterator Helpers ============ */

size_t ufunc_get_inner_loop_size(NDArray* arr) {
    if (arr == NULL || arr->ndim == 0) {
        return 1;
    }

    /* Check if innermost dimension is contiguous */
    size_t elem_size = dtype_size(arr->dtype);
    if (arr->strides[arr->ndim - 1] == (int32_t)elem_size) {
        return arr->shape[arr->ndim - 1];
    }

    return 1;
}

bool ufunc_is_contiguous(NDArray* arr) {
    if (arr == NULL) return false;
    return ndarray_is_c_contiguous(arr);
}

bool ufunc_binary_contiguous(NDArray* arr1, NDArray* arr2) {
    if (!ufunc_is_contiguous(arr1) || !ufunc_is_contiguous(arr2)) {
        return false;
    }

    /* Must have same shape for fast path */
    if (arr1->ndim != arr2->ndim) {
        return false;
    }

    for (int32_t i = 0; i < arr1->ndim; i++) {
        if (arr1->shape[i] != arr2->shape[i]) {
            return false;
        }
    }

    return true;
}

/* ============ Type Selection ============ */

UnaryLoopFunc ufunc_select_unary_loop(DType dtype,
                                       UnaryLoopFunc f64_loop,
                                       UnaryLoopFunc f32_loop,
                                       UnaryLoopFunc i32_loop) {
    switch (dtype) {
        case DTYPE_FLOAT64:
            return f64_loop;
        case DTYPE_FLOAT32:
            return f32_loop;
        case DTYPE_INT32:
            return i32_loop;
        default:
            /* Other types fall back to f64 */
            return f64_loop;
    }
}

BinaryLoopFunc ufunc_select_binary_loop(DType dtype,
                                         BinaryLoopFunc f64_loop,
                                         BinaryLoopFunc f32_loop,
                                         BinaryLoopFunc i32_loop) {
    switch (dtype) {
        case DTYPE_FLOAT64:
            return f64_loop;
        case DTYPE_FLOAT32:
            return f32_loop;
        case DTYPE_INT32:
            return i32_loop;
        default:
            return f64_loop;
    }
}

DType ufunc_result_dtype(DType dtype1, DType dtype2) {
    return dtype_promote(dtype1, dtype2);
}

/* ============ Multi-dimensional Iteration ============ */

/**
 * Internal: iterate over all elements using multi-dimensional indices.
 * Used for non-contiguous arrays.
 */
static void unary_strided_loop(NDArray* input, NDArray* output, UnaryLoopFunc loop) {
    size_t n = input->size;
    size_t in_elem_size = dtype_size(input->dtype);
    size_t out_elem_size = dtype_size(output->dtype);

    if (n == 0) return;

    /* For 0-d arrays */
    if (input->ndim == 0) {
        loop((const char*)input->data, (char*)output->data, 1, 0, 0);
        return;
    }

    /* Check for contiguous fast path */
    if (ndarray_is_c_contiguous(input) && ndarray_is_c_contiguous(output)) {
        loop((const char*)input->data, (char*)output->data,
             n, in_elem_size, out_elem_size);
        return;
    }

    /* Strided iteration using multi-index */
    int32_t* indices = (int32_t*)calloc(input->ndim, sizeof(int32_t));
    if (indices == NULL) return;

    for (size_t i = 0; i < n; i++) {
        /* Compute byte offsets */
        size_t in_offset = 0;
        size_t out_offset = 0;
        for (int32_t d = 0; d < input->ndim; d++) {
            in_offset += indices[d] * input->strides[d];
            out_offset += indices[d] * output->strides[d];
        }

        /* Process single element */
        loop((const char*)input->data + in_offset,
             (char*)output->data + out_offset,
             1, 0, 0);

        /* Increment multi-index (row-major order) */
        for (int32_t d = input->ndim - 1; d >= 0; d--) {
            indices[d]++;
            if (indices[d] < input->shape[d]) {
                break;
            }
            indices[d] = 0;
        }
    }

    free(indices);
}

/**
 * Internal: iterate over two arrays with broadcasting.
 */
static void binary_strided_loop(NDArray* in1, NDArray* in2, NDArray* output,
                                 int32_t* strides1, int32_t* strides2,
                                 BinaryLoopFunc loop) {
    size_t n = output->size;
    int32_t ndim = output->ndim;
    size_t out_elem_size = dtype_size(output->dtype);

    if (n == 0) return;

    /* For 0-d arrays */
    if (ndim == 0) {
        loop((const char*)in1->data, (const char*)in2->data,
             (char*)output->data, 1, 0, 0, 0);
        return;
    }

    /* Check for contiguous fast path (same shape, no broadcasting) */
    if (ufunc_binary_contiguous(in1, in2) && ndarray_is_c_contiguous(output)) {
        size_t elem_size = dtype_size(in1->dtype);
        loop((const char*)in1->data, (const char*)in2->data,
             (char*)output->data, n, elem_size, elem_size, out_elem_size);
        return;
    }

    /* Strided iteration with broadcasting */
    int32_t* indices = (int32_t*)calloc(ndim, sizeof(int32_t));
    if (indices == NULL) return;

    for (size_t i = 0; i < n; i++) {
        /* Compute byte offsets using broadcast strides */
        size_t off1 = 0, off2 = 0, out_off = 0;
        for (int32_t d = 0; d < ndim; d++) {
            off1 += indices[d] * strides1[d];
            off2 += indices[d] * strides2[d];
            out_off += indices[d] * output->strides[d];
        }

        /* Process single element */
        loop((const char*)in1->data + off1,
             (const char*)in2->data + off2,
             (char*)output->data + out_off,
             1, 0, 0, 0);

        /* Increment multi-index */
        for (int32_t d = ndim - 1; d >= 0; d--) {
            indices[d]++;
            if (indices[d] < output->shape[d]) {
                break;
            }
            indices[d] = 0;
        }
    }

    free(indices);
}

/* ============ High-Level Apply Functions ============ */

NDArray* ufunc_apply_unary(NDArray* input, DType out_dtype, UnaryLoopFunc loop) {
    if (input == NULL || loop == NULL) {
        return NULL;
    }

    /* Create output array with same shape */
    NDArray* output = ndarray_empty(input->ndim, input->shape, out_dtype);
    if (output == NULL) {
        return NULL;
    }

    /* Apply the loop */
    unary_strided_loop(input, output, loop);

    return output;
}

NDArray* ufunc_apply_binary(NDArray* in1, NDArray* in2, DType out_dtype, BinaryLoopFunc loop) {
    if (in1 == NULL || in2 == NULL || loop == NULL) {
        return NULL;
    }

    /* Compute broadcast shape */
    int32_t out_shape[NPY_MAXDIMS];
    int32_t out_ndim;

    if (broadcast_shapes(in1->shape, in1->ndim, in2->shape, in2->ndim,
                         out_shape, &out_ndim) != 0) {
        /* Shapes not broadcastable */
        return NULL;
    }

    /* Compute broadcast strides for each input */
    int32_t strides1[NPY_MAXDIMS];
    int32_t strides2[NPY_MAXDIMS];

    if (broadcast_strides(in1, out_shape, out_ndim, strides1) != 0 ||
        broadcast_strides(in2, out_shape, out_ndim, strides2) != 0) {
        return NULL;
    }

    /* Create output array */
    NDArray* output = ndarray_empty(out_ndim, out_shape, out_dtype);
    if (output == NULL) {
        return NULL;
    }

    /* Apply the loop with broadcasting */
    binary_strided_loop(in1, in2, output, strides1, strides2, loop);

    return output;
}

NDArray* ufunc_apply_binary_cmp(NDArray* in1, NDArray* in2, BinaryLoopFunc loop) {
    /* Comparison always outputs bool (stored as uint8) */
    return ufunc_apply_binary(in1, in2, DTYPE_BOOL, loop);
}
