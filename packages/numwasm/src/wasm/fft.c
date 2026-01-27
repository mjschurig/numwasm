/**
 * NumJS FFT Implementation
 *
 * Implements:
 * - Cooley-Tukey radix-2 FFT for power-of-2 sizes
 * - Bluestein's algorithm (chirp-z transform) for arbitrary sizes
 * - Real FFT with Hermitian symmetry optimization
 */

#include "fft.h"
#include "ndarray.h"
#include "dtype.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#define EXPORT EMSCRIPTEN_KEEPALIVE
#else
#define EXPORT
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ============ Utility Functions ============ */

EXPORT int32_t fft_is_power_of_2(int32_t n)
{
    return (n > 0) && ((n & (n - 1)) == 0);
}

EXPORT int32_t fft_next_power_of_2(int32_t n)
{
    if (n <= 0) return 1;
    int32_t p = 1;
    while (p < n) p <<= 1;
    return p;
}

EXPORT int32_t fft_bit_reverse(int32_t i, int32_t log2n)
{
    int32_t rev = 0;
    for (int32_t j = 0; j < log2n; j++) {
        rev = (rev << 1) | (i & 1);
        i >>= 1;
    }
    return rev;
}

EXPORT void fft_bit_reverse_permute(double* data, int32_t n)
{
    /* Compute log2(n) */
    int32_t log2n = 0;
    int32_t temp = n;
    while (temp > 1) {
        temp >>= 1;
        log2n++;
    }

    /* Swap elements with their bit-reversed indices */
    for (int32_t i = 0; i < n; i++) {
        int32_t j = fft_bit_reverse(i, log2n);
        if (j > i) {
            /* Swap complex elements (2 doubles each) */
            double temp_re = data[2 * i];
            double temp_im = data[2 * i + 1];
            data[2 * i]     = data[2 * j];
            data[2 * i + 1] = data[2 * j + 1];
            data[2 * j]     = temp_re;
            data[2 * j + 1] = temp_im;
        }
    }
}

/* ============ Radix-2 Cooley-Tukey FFT ============ */

EXPORT int32_t fft_radix2(double* data, int32_t n, int32_t inverse)
{
    if (!fft_is_power_of_2(n)) {
        return -1;  /* Error: n must be power of 2 */
    }

    if (n <= 1) return 0;

    /* Bit-reversal permutation */
    fft_bit_reverse_permute(data, n);

    /* Direction: -1 for forward, +1 for inverse */
    double sign = inverse ? 1.0 : -1.0;

    /* Cooley-Tukey iterative algorithm (decimation-in-time) */
    for (int32_t len = 2; len <= n; len <<= 1) {
        double theta = sign * 2.0 * M_PI / len;
        double wpr = cos(theta);
        double wpi = sin(theta);

        for (int32_t i = 0; i < n; i += len) {
            double wr = 1.0;
            double wi = 0.0;

            for (int32_t j = 0; j < len / 2; j++) {
                int32_t idx1 = i + j;
                int32_t idx2 = i + j + len / 2;

                /* Butterfly operation:
                 * data[idx1] = data[idx1] + W * data[idx2]
                 * data[idx2] = data[idx1] - W * data[idx2]
                 */
                double t_re = wr * data[2 * idx2]     - wi * data[2 * idx2 + 1];
                double t_im = wr * data[2 * idx2 + 1] + wi * data[2 * idx2];

                data[2 * idx2]     = data[2 * idx1]     - t_re;
                data[2 * idx2 + 1] = data[2 * idx1 + 1] - t_im;
                data[2 * idx1]     = data[2 * idx1]     + t_re;
                data[2 * idx1 + 1] = data[2 * idx1 + 1] + t_im;

                /* Update twiddle factor: W = W * W_len */
                double temp_w = wr;
                wr = wr * wpr - wi * wpi;
                wi = temp_w * wpi + wi * wpr;
            }
        }
    }

    return 0;
}

/* ============ Bluestein's Algorithm ============ */

EXPORT int32_t fft_bluestein(double* data, int32_t n, int32_t inverse, double* work)
{
    if (n <= 0 || work == NULL) return -1;
    if (n == 1) return 0;

    /* Convolution size: smallest power of 2 >= 2*n - 1 */
    int32_t m = fft_next_power_of_2(2 * n - 1);

    double sign = inverse ? 1.0 : -1.0;
    double theta = sign * M_PI / n;

    /* Workspace layout:
     * work[0..2m-1]       = chirp signal
     * work[2m..4m-1]      = conv_a (input convolved)
     * work[4m..6m-1]      = conv_b (chirp for convolution)
     */
    double* chirp  = work;
    double* conv_a = work + 2 * m;
    double* conv_b = work + 4 * m;

    /* Initialize chirp signal: w[k] = exp(sign * i * pi * k^2 / n) */
    for (int32_t k = 0; k < n; k++) {
        double angle = theta * (double)k * (double)k;
        chirp[2 * k]     = cos(angle);
        chirp[2 * k + 1] = sin(angle);
    }

    /* Build conv_a[k] = x[k] * conj(chirp[k]) for k < n, zero for k >= n */
    for (int32_t k = 0; k < n; k++) {
        double x_re = data[2 * k];
        double x_im = data[2 * k + 1];
        double c_re = chirp[2 * k];
        double c_im = -chirp[2 * k + 1];  /* Conjugate */
        conv_a[2 * k]     = x_re * c_re - x_im * c_im;
        conv_a[2 * k + 1] = x_re * c_im + x_im * c_re;
    }
    /* Zero-pad */
    for (int32_t k = n; k < m; k++) {
        conv_a[2 * k]     = 0.0;
        conv_a[2 * k + 1] = 0.0;
    }

    /* Build conv_b: chirp[k] for k < n, chirp[m-k] for k > m-n, zero elsewhere */
    for (int32_t k = 0; k < m; k++) {
        conv_b[2 * k]     = 0.0;
        conv_b[2 * k + 1] = 0.0;
    }
    for (int32_t k = 0; k < n; k++) {
        conv_b[2 * k]     = chirp[2 * k];
        conv_b[2 * k + 1] = chirp[2 * k + 1];
    }
    for (int32_t k = 1; k < n; k++) {
        conv_b[2 * (m - k)]     = chirp[2 * k];
        conv_b[2 * (m - k) + 1] = chirp[2 * k + 1];
    }

    /* Circular convolution via FFT:
     * result = IFFT(FFT(conv_a) * FFT(conv_b))
     */
    fft_radix2(conv_a, m, 0);  /* Forward FFT */
    fft_radix2(conv_b, m, 0);  /* Forward FFT */

    /* Point-wise complex multiplication */
    for (int32_t k = 0; k < m; k++) {
        double a_re = conv_a[2 * k];
        double a_im = conv_a[2 * k + 1];
        double b_re = conv_b[2 * k];
        double b_im = conv_b[2 * k + 1];
        conv_a[2 * k]     = a_re * b_re - a_im * b_im;
        conv_a[2 * k + 1] = a_re * b_im + a_im * b_re;
    }

    /* Inverse FFT */
    fft_radix2(conv_a, m, 1);

    /* Scale by 1/m and multiply by chirp */
    double scale = 1.0 / m;
    for (int32_t k = 0; k < n; k++) {
        double c_re = conv_a[2 * k] * scale;
        double c_im = conv_a[2 * k + 1] * scale;
        double w_re = chirp[2 * k];
        double w_im = -chirp[2 * k + 1];  /* Conjugate */
        data[2 * k]     = c_re * w_re - c_im * w_im;
        data[2 * k + 1] = c_re * w_im + c_im * w_re;
    }

    return 0;
}

/* ============ General FFT Dispatcher ============ */

EXPORT int32_t fft_complex(double* data, int32_t n, int32_t inverse, double* work)
{
    if (n <= 0) return -1;
    if (n == 1) return 0;

    if (fft_is_power_of_2(n)) {
        return fft_radix2(data, n, inverse);
    } else {
        if (work == NULL) return -1;
        return fft_bluestein(data, n, inverse, work);
    }
}

/* ============ Real FFT Functions ============ */

EXPORT int32_t fft_rfft(const double* real_in, double* out, int32_t n, double* work)
{
    if (n <= 0 || real_in == NULL || out == NULL || work == NULL) return -1;

    /* Pack real input as complex with zero imaginary */
    double* complex_data = work;
    for (int32_t i = 0; i < n; i++) {
        complex_data[2 * i]     = real_in[i];
        complex_data[2 * i + 1] = 0.0;
    }

    /* Perform complex FFT */
    double* fft_work = work + 2 * n;
    int32_t result = fft_complex(complex_data, n, 0, fft_work);
    if (result != 0) return result;

    /* Copy first n/2+1 values to output (Hermitian symmetry) */
    int32_t out_size = n / 2 + 1;
    for (int32_t i = 0; i < out_size; i++) {
        out[2 * i]     = complex_data[2 * i];
        out[2 * i + 1] = complex_data[2 * i + 1];
    }

    return 0;
}

EXPORT int32_t fft_irfft(const double* complex_in, double* out, int32_t n, double* work)
{
    if (n <= 0 || complex_in == NULL || out == NULL || work == NULL) return -1;

    int32_t in_size = n / 2 + 1;

    /* Reconstruct full complex spectrum using Hermitian symmetry:
     * X[n-k] = conj(X[k]) for k > 0
     */
    double* complex_data = work;
    for (int32_t i = 0; i < in_size; i++) {
        complex_data[2 * i]     = complex_in[2 * i];
        complex_data[2 * i + 1] = complex_in[2 * i + 1];
    }

    /* Fill negative frequencies using Hermitian symmetry */
    for (int32_t i = in_size; i < n; i++) {
        int32_t j = n - i;
        complex_data[2 * i]     =  complex_in[2 * j];      /* Real part */
        complex_data[2 * i + 1] = -complex_in[2 * j + 1];  /* Conjugate: negate imag */
    }

    /* Perform inverse FFT */
    double* fft_work = work + 2 * n;
    int32_t result = fft_complex(complex_data, n, 1, fft_work);
    if (result != 0) return result;

    /* Extract real parts and scale by 1/n */
    double scale = 1.0 / n;
    for (int32_t i = 0; i < n; i++) {
        out[i] = complex_data[2 * i] * scale;
    }

    return 0;
}

/* ============ NDArray-Level Operations ============ */

/**
 * Get workspace size needed for FFT of size n.
 */
static size_t fft_workspace_size(int32_t n)
{
    if (fft_is_power_of_2(n)) {
        return 0;  /* Radix-2 is in-place */
    }
    /* Bluestein needs 6 * m complex doubles where m = next_power_of_2(2n-1) */
    int32_t m = fft_next_power_of_2(2 * n - 1);
    return 6 * m * sizeof(double) * 2;
}

/**
 * Get workspace size for rfft/irfft of size n.
 */
static size_t rfft_workspace_size(int32_t n)
{
    /* Need 2*n for complex data + fft_workspace_size(n) */
    size_t fft_ws = fft_workspace_size(n);
    if (fft_ws == 0) {
        /* For power-of-2, we still need space for complex data */
        return 2 * n * sizeof(double);
    }
    return 2 * n * sizeof(double) + fft_ws;
}

EXPORT NDArray* ndarray_fft(NDArray* arr, int32_t n, int32_t axis, int32_t inverse)
{
    if (arr == NULL || arr->data == NULL) return NULL;

    /* Normalize axis */
    if (axis < 0) axis += arr->ndim;
    if (axis < 0 || axis >= arr->ndim) return NULL;

    /* Determine FFT size */
    int32_t axis_size = arr->shape[axis];
    if (n < 0) n = axis_size;

    /* Output shape is same as input (or with n along axis if different) */
    int32_t* out_shape = (int32_t*)malloc(arr->ndim * sizeof(int32_t));
    if (out_shape == NULL) return NULL;

    for (int32_t i = 0; i < arr->ndim; i++) {
        out_shape[i] = (i == axis) ? n : arr->shape[i];
    }

    /* Create output array (always complex128) */
    NDArray* result = ndarray_create(arr->ndim, out_shape, DTYPE_COMPLEX128);
    free(out_shape);
    if (result == NULL) return NULL;

    /* Allocate workspace */
    size_t ws_size = fft_workspace_size(n);
    double* workspace = NULL;
    if (ws_size > 0) {
        workspace = (double*)malloc(ws_size);
        if (workspace == NULL) {
            ndarray_free(result);
            return NULL;
        }
    }

    /* Compute iteration parameters:
     * outer = product of dimensions before axis
     * inner = product of dimensions after axis
     */
    size_t outer = 1, inner = 1;
    for (int32_t i = 0; i < axis; i++) outer *= arr->shape[i];
    for (int32_t i = axis + 1; i < arr->ndim; i++) inner *= arr->shape[i];

    /* Allocate buffer for 1D slice */
    double* slice_buf = (double*)malloc(n * 2 * sizeof(double));
    if (slice_buf == NULL) {
        if (workspace) free(workspace);
        ndarray_free(result);
        return NULL;
    }

    int is_complex = dtype_is_complex(arr->dtype);

    /* Process each 1D slice along axis */
    for (size_t o = 0; o < outer; o++) {
        for (size_t in = 0; in < inner; in++) {
            /* Extract 1D slice into buffer */
            for (int32_t k = 0; k < axis_size && k < n; k++) {
                size_t flat_idx = o * axis_size * inner + k * inner + in;
                if (is_complex) {
                    slice_buf[2 * k]     = ndarray_get_complex_real(arr, flat_idx);
                    slice_buf[2 * k + 1] = ndarray_get_complex_imag(arr, flat_idx);
                } else {
                    slice_buf[2 * k]     = ndarray_get_flat(arr, flat_idx);
                    slice_buf[2 * k + 1] = 0.0;
                }
            }
            /* Zero-pad if n > axis_size */
            for (int32_t k = axis_size; k < n; k++) {
                slice_buf[2 * k]     = 0.0;
                slice_buf[2 * k + 1] = 0.0;
            }

            /* Apply FFT to slice */
            fft_complex(slice_buf, n, inverse, workspace);

            /* Apply scaling for inverse FFT */
            if (inverse) {
                double scale = 1.0 / n;
                for (int32_t k = 0; k < n; k++) {
                    slice_buf[2 * k]     *= scale;
                    slice_buf[2 * k + 1] *= scale;
                }
            }

            /* Write results back to output */
            for (int32_t k = 0; k < n; k++) {
                size_t out_idx = o * n * inner + k * inner + in;
                ndarray_set_complex(result, out_idx, slice_buf[2 * k], slice_buf[2 * k + 1]);
            }
        }
    }

    free(slice_buf);
    if (workspace) free(workspace);

    return result;
}

EXPORT NDArray* ndarray_rfft(NDArray* arr, int32_t n, int32_t axis)
{
    if (arr == NULL || arr->data == NULL) return NULL;

    /* Real FFT should only be applied to real arrays */
    if (dtype_is_complex(arr->dtype)) {
        /* For complex input, just use regular fft */
        return ndarray_fft(arr, n, axis, 0);
    }

    /* Normalize axis */
    if (axis < 0) axis += arr->ndim;
    if (axis < 0 || axis >= arr->ndim) return NULL;

    /* Determine FFT size */
    int32_t axis_size = arr->shape[axis];
    if (n < 0) n = axis_size;

    /* Output size along axis is n/2 + 1 */
    int32_t out_axis_size = n / 2 + 1;

    /* Output shape */
    int32_t* out_shape = (int32_t*)malloc(arr->ndim * sizeof(int32_t));
    if (out_shape == NULL) return NULL;

    for (int32_t i = 0; i < arr->ndim; i++) {
        out_shape[i] = (i == axis) ? out_axis_size : arr->shape[i];
    }

    /* Create output array (complex128) */
    NDArray* result = ndarray_create(arr->ndim, out_shape, DTYPE_COMPLEX128);
    free(out_shape);
    if (result == NULL) return NULL;

    /* Allocate workspace */
    size_t ws_size = rfft_workspace_size(n);
    double* workspace = (double*)malloc(ws_size);
    if (workspace == NULL) {
        ndarray_free(result);
        return NULL;
    }

    /* Compute iteration parameters */
    size_t outer = 1, inner = 1;
    for (int32_t i = 0; i < axis; i++) outer *= arr->shape[i];
    for (int32_t i = axis + 1; i < arr->ndim; i++) inner *= arr->shape[i];

    /* Allocate buffers */
    double* real_buf = (double*)malloc(n * sizeof(double));
    double* out_buf = (double*)malloc(out_axis_size * 2 * sizeof(double));
    if (real_buf == NULL || out_buf == NULL) {
        if (real_buf) free(real_buf);
        if (out_buf) free(out_buf);
        free(workspace);
        ndarray_free(result);
        return NULL;
    }

    /* Process each 1D slice */
    for (size_t o = 0; o < outer; o++) {
        for (size_t in = 0; in < inner; in++) {
            /* Extract real 1D slice */
            for (int32_t k = 0; k < axis_size && k < n; k++) {
                size_t flat_idx = o * axis_size * inner + k * inner + in;
                real_buf[k] = ndarray_get_flat(arr, flat_idx);
            }
            /* Zero-pad if n > axis_size */
            for (int32_t k = axis_size; k < n; k++) {
                real_buf[k] = 0.0;
            }

            /* Apply rfft */
            fft_rfft(real_buf, out_buf, n, workspace);

            /* Write results */
            for (int32_t k = 0; k < out_axis_size; k++) {
                size_t out_idx = o * out_axis_size * inner + k * inner + in;
                ndarray_set_complex(result, out_idx, out_buf[2 * k], out_buf[2 * k + 1]);
            }
        }
    }

    free(real_buf);
    free(out_buf);
    free(workspace);

    return result;
}

EXPORT NDArray* ndarray_irfft(NDArray* arr, int32_t n, int32_t axis)
{
    if (arr == NULL || arr->data == NULL) return NULL;

    /* Normalize axis */
    if (axis < 0) axis += arr->ndim;
    if (axis < 0 || axis >= arr->ndim) return NULL;

    /* Input has n/2+1 complex values along axis */
    int32_t in_axis_size = arr->shape[axis];

    /* Determine output size */
    if (n < 0) {
        n = 2 * (in_axis_size - 1);  /* Default: reconstruct original size */
    }

    /* Output shape (real output) */
    int32_t* out_shape = (int32_t*)malloc(arr->ndim * sizeof(int32_t));
    if (out_shape == NULL) return NULL;

    for (int32_t i = 0; i < arr->ndim; i++) {
        out_shape[i] = (i == axis) ? n : arr->shape[i];
    }

    /* Create output array (float64 - real) */
    NDArray* result = ndarray_create(arr->ndim, out_shape, DTYPE_FLOAT64);
    free(out_shape);
    if (result == NULL) return NULL;

    /* Allocate workspace */
    size_t ws_size = rfft_workspace_size(n);
    double* workspace = (double*)malloc(ws_size);
    if (workspace == NULL) {
        ndarray_free(result);
        return NULL;
    }

    /* Compute iteration parameters */
    size_t outer = 1, inner = 1;
    for (int32_t i = 0; i < axis; i++) outer *= arr->shape[i];
    for (int32_t i = axis + 1; i < arr->ndim; i++) inner *= arr->shape[i];

    /* Allocate buffers */
    int32_t expected_in_size = n / 2 + 1;
    double* complex_buf = (double*)malloc(expected_in_size * 2 * sizeof(double));
    double* real_out = (double*)malloc(n * sizeof(double));
    if (complex_buf == NULL || real_out == NULL) {
        if (complex_buf) free(complex_buf);
        if (real_out) free(real_out);
        free(workspace);
        ndarray_free(result);
        return NULL;
    }

    int is_complex = dtype_is_complex(arr->dtype);

    /* Process each 1D slice */
    for (size_t o = 0; o < outer; o++) {
        for (size_t in = 0; in < inner; in++) {
            /* Extract complex 1D slice */
            for (int32_t k = 0; k < in_axis_size && k < expected_in_size; k++) {
                size_t flat_idx = o * in_axis_size * inner + k * inner + in;
                if (is_complex) {
                    complex_buf[2 * k]     = ndarray_get_complex_real(arr, flat_idx);
                    complex_buf[2 * k + 1] = ndarray_get_complex_imag(arr, flat_idx);
                } else {
                    complex_buf[2 * k]     = ndarray_get_flat(arr, flat_idx);
                    complex_buf[2 * k + 1] = 0.0;
                }
            }
            /* Zero-pad if needed */
            for (int32_t k = in_axis_size; k < expected_in_size; k++) {
                complex_buf[2 * k]     = 0.0;
                complex_buf[2 * k + 1] = 0.0;
            }

            /* Apply irfft */
            fft_irfft(complex_buf, real_out, n, workspace);

            /* Write results */
            for (int32_t k = 0; k < n; k++) {
                size_t out_idx = o * n * inner + k * inner + in;
                ndarray_set_flat(result, out_idx, real_out[k]);
            }
        }
    }

    free(complex_buf);
    free(real_out);
    free(workspace);

    return result;
}
