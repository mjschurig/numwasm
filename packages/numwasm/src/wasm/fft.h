/**
 * NumJS FFT (Fast Fourier Transform) Module
 *
 * Provides FFT algorithms for complex and real arrays.
 * Implements Cooley-Tukey radix-2 for power-of-2 sizes and
 * Bluestein's algorithm for arbitrary sizes.
 */

#ifndef NUMJS_FFT_H
#define NUMJS_FFT_H

#include <stdint.h>
#include <stddef.h>
#include "ndarray.h"

/* ============ Utility Functions ============ */

/**
 * Check if n is a power of 2.
 *
 * @param n  Number to check
 * @return   1 if power of 2, 0 otherwise
 */
int32_t fft_is_power_of_2(int32_t n);

/**
 * Find the smallest power of 2 >= n.
 *
 * @param n  Input number
 * @return   Smallest power of 2 >= n
 */
int32_t fft_next_power_of_2(int32_t n);

/**
 * Compute bit-reversed index.
 *
 * @param i      Input index
 * @param log2n  log2(n) where n is FFT size
 * @return       Bit-reversed index
 */
int32_t fft_bit_reverse(int32_t i, int32_t log2n);

/**
 * Perform in-place bit-reversal permutation on complex array.
 *
 * @param data  Complex array (interleaved real/imag, size 2*n doubles)
 * @param n     Number of complex elements
 */
void fft_bit_reverse_permute(double* data, int32_t n);

/* ============ Core FFT Algorithms ============ */

/**
 * In-place radix-2 Cooley-Tukey FFT.
 * Input must be power-of-2 length.
 *
 * @param data     Complex array (interleaved real/imag, size 2*n doubles)
 * @param n        Number of complex elements (must be power of 2)
 * @param inverse  0 for forward FFT, 1 for inverse FFT
 * @return         0 on success, -1 if n is not power of 2
 */
int32_t fft_radix2(double* data, int32_t n, int32_t inverse);

/**
 * FFT for arbitrary size using Bluestein's algorithm (chirp-z transform).
 * Converts arbitrary-size FFT to power-of-2 convolution.
 *
 * @param data     Complex input/output array (interleaved, size 2*n doubles)
 * @param n        Number of complex elements (any positive integer)
 * @param inverse  0 for forward, 1 for inverse
 * @param work     Workspace (size >= 6 * next_power_of_2(2*n-1) doubles)
 * @return         0 on success, -1 on error
 */
int32_t fft_bluestein(double* data, int32_t n, int32_t inverse, double* work);

/**
 * General FFT dispatcher - calls radix-2 or Bluestein based on size.
 *
 * @param data     Complex input/output array
 * @param n        Number of complex elements
 * @param inverse  0 for forward, 1 for inverse
 * @param work     Workspace (can be NULL for power-of-2 sizes)
 * @return         0 on success, -1 on error
 */
int32_t fft_complex(double* data, int32_t n, int32_t inverse, double* work);

/* ============ Real FFT Functions ============ */

/**
 * Real-to-complex FFT exploiting Hermitian symmetry.
 * Output is n/2 + 1 complex values.
 *
 * @param real_in  Real input array (size n doubles)
 * @param out      Complex output (size 2*(n/2+1) doubles, interleaved)
 * @param n        Number of real input samples
 * @param work     Workspace (size >= 2*n + 6*next_power_of_2(2*n-1) doubles)
 * @return         0 on success, -1 on error
 */
int32_t fft_rfft(const double* real_in, double* out, int32_t n, double* work);

/**
 * Complex-to-real inverse FFT.
 * Input is n/2 + 1 complex values with Hermitian symmetry.
 *
 * @param complex_in  Complex input (size 2*(n/2+1) doubles)
 * @param out         Real output (size n doubles)
 * @param n           Number of real output samples
 * @param work        Workspace (size >= 2*n + 6*next_power_of_2(2*n-1) doubles)
 * @return            0 on success, -1 on error
 */
int32_t fft_irfft(const double* complex_in, double* out, int32_t n, double* work);

/* ============ NDArray-Level Operations ============ */

/**
 * Compute FFT of NDArray along specified axis.
 *
 * @param arr      Input array (real or complex)
 * @param n        FFT size (-1 to use axis size)
 * @param axis     Axis along which to compute FFT
 * @param inverse  0 for forward FFT, 1 for inverse FFT
 * @return         New complex array with FFT result, or NULL on error
 */
NDArray* ndarray_fft(NDArray* arr, int32_t n, int32_t axis, int32_t inverse);

/**
 * Compute real-to-complex FFT of NDArray along specified axis.
 * Output has size n/2+1 along the transformed axis.
 *
 * @param arr   Input array (real-valued)
 * @param n     FFT size (-1 to use axis size)
 * @param axis  Axis along which to compute FFT
 * @return      New complex array with FFT result, or NULL on error
 */
NDArray* ndarray_rfft(NDArray* arr, int32_t n, int32_t axis);

/**
 * Compute complex-to-real inverse FFT of NDArray along specified axis.
 * Input should have Hermitian symmetry (from rfft).
 *
 * @param arr   Input array (complex, from rfft)
 * @param n     Output size (-1 to infer from input: 2*(axis_size-1))
 * @param axis  Axis along which to compute inverse FFT
 * @return      New real array with inverse FFT result, or NULL on error
 */
NDArray* ndarray_irfft(NDArray* arr, int32_t n, int32_t axis);

#endif /* NUMJS_FFT_H */
