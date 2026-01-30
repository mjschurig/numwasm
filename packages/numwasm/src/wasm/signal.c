/**
 * NumWasm Signal Processing Functions
 *
 * Implements 1D convolution and correlation operations.
 */

#include "signal.h"
#include <stdlib.h>
#include <math.h>

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#define EXPORT EMSCRIPTEN_KEEPALIVE
#else
#define EXPORT
#endif

/**
 * 1D discrete convolution: (a * v)[n] = sum_m( a[m] * v[n - m] )
 *
 * Convolution is equivalent to correlation with the second array reversed.
 *
 * Output sizes by mode:
 *   - FULL:  N + M - 1
 *   - SAME:  max(N, M)
 *   - VALID: |N - M| + 1 (only where both signals overlap completely)
 */
EXPORT NDArray* ndarray_convolve(NDArray* a, NDArray* v, int32_t mode) {
    if (!a || !a->data || !v || !v->data) return NULL;
    if (a->ndim != 1 || v->ndim != 1) return NULL;  /* 1D only */

    int32_t na = a->shape[0];
    int32_t nv = v->shape[0];

    if (na == 0 || nv == 0) {
        int32_t out_shape[1] = { 0 };
        return ndarray_create(1, out_shape, DTYPE_FLOAT64);
    }

    /* Swap if v is longer (optimization and consistent behavior) */
    NDArray* arr_long = a;
    NDArray* arr_short = v;
    int32_t n_long = na;
    int32_t n_short = nv;

    if (nv > na) {
        arr_long = v;
        arr_short = a;
        n_long = nv;
        n_short = na;
    }

    int32_t full_len = n_long + n_short - 1;
    int32_t out_len;
    int32_t out_start = 0;

    switch (mode) {
        case CONVOLVE_MODE_FULL:
            out_len = full_len;
            out_start = 0;
            break;
        case CONVOLVE_MODE_SAME:
            out_len = n_long;  /* max(na, nv) but we swapped so n_long >= n_short */
            out_start = (n_short - 1) / 2;
            break;
        case CONVOLVE_MODE_VALID:
            out_len = n_long - n_short + 1;
            out_start = n_short - 1;
            break;
        default:
            return NULL;
    }

    if (out_len <= 0) {
        int32_t out_shape[1] = { 0 };
        return ndarray_create(1, out_shape, DTYPE_FLOAT64);
    }

    int32_t out_shape[1] = { out_len };
    NDArray* result = ndarray_create(1, out_shape, DTYPE_FLOAT64);
    if (!result) return NULL;

    /* Direct convolution: O(N*M) - efficient for small kernels
     * For each output position k (in full output), compute:
     *   result[k] = sum_j( arr_long[k - n_short + 1 + j] * arr_short[n_short - 1 - j] )
     * Note: arr_short is effectively reversed for convolution
     */
    for (int32_t k = 0; k < out_len; k++) {
        int32_t full_k = k + out_start;
        double sum = 0.0;

        for (int32_t j = 0; j < n_short; j++) {
            /* Convolution reverses the kernel */
            int32_t i = full_k - (n_short - 1 - j);
            if (i >= 0 && i < n_long) {
                double a_val = ndarray_get_flat(arr_long, i);
                double v_val = ndarray_get_flat(arr_short, j);
                sum += a_val * v_val;
            }
        }
        ndarray_set_flat(result, k, sum);
    }

    return result;
}

/**
 * 1D cross-correlation: (a <*> v)[n] = sum_m( a[n + m] * conj(v[m]) )
 *
 * For real arrays, this is like convolution but without reversing v.
 * Default mode is VALID (unlike convolution's FULL).
 *
 * Output sizes by mode:
 *   - FULL:  N + M - 1
 *   - SAME:  max(N, M)
 *   - VALID: |N - M| + 1
 */
EXPORT NDArray* ndarray_correlate(NDArray* a, NDArray* v, int32_t mode) {
    if (!a || !a->data || !v || !v->data) return NULL;
    if (a->ndim != 1 || v->ndim != 1) return NULL;  /* 1D only */

    int32_t na = a->shape[0];
    int32_t nv = v->shape[0];

    if (na == 0 || nv == 0) {
        int32_t out_shape[1] = { 0 };
        return ndarray_create(1, out_shape, DTYPE_FLOAT64);
    }

    int32_t full_len = na + nv - 1;
    int32_t out_len;
    int32_t out_start = 0;

    switch (mode) {
        case CONVOLVE_MODE_FULL:
            out_len = full_len;
            out_start = 0;
            break;
        case CONVOLVE_MODE_SAME:
            out_len = (na > nv) ? na : nv;
            out_start = (nv - 1) / 2;
            break;
        case CONVOLVE_MODE_VALID:
        default:
            /* VALID is the default for correlate */
            if (na >= nv) {
                out_len = na - nv + 1;
                out_start = nv - 1;
            } else {
                out_len = nv - na + 1;
                out_start = na - 1;
            }
            break;
    }

    if (out_len <= 0) {
        int32_t out_shape[1] = { 0 };
        return ndarray_create(1, out_shape, DTYPE_FLOAT64);
    }

    int32_t out_shape[1] = { out_len };
    NDArray* result = ndarray_create(1, out_shape, DTYPE_FLOAT64);
    if (!result) return NULL;

    /* Direct correlation: O(N*M)
     * For each output position k (in full output), compute:
     *   result[k] = sum_j( a[k - nv + 1 + j] * v[j] )
     * Note: v is NOT reversed (unlike convolution)
     */
    for (int32_t k = 0; k < out_len; k++) {
        int32_t full_k = k + out_start;
        double sum = 0.0;

        for (int32_t j = 0; j < nv; j++) {
            int32_t i = full_k - nv + 1 + j;
            if (i >= 0 && i < na) {
                double a_val = ndarray_get_flat(a, i);
                double v_val = ndarray_get_flat(v, j);
                sum += a_val * v_val;
            }
        }
        ndarray_set_flat(result, k, sum);
    }

    return result;
}
