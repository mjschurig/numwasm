#ifndef NUMJS_PAIRWISE_SUM_H
#define NUMJS_PAIRWISE_SUM_H

#include <stdint.h>
#include <stddef.h>

/*
 * Cutoff blocksize for pairwise summation.
 * Ported from NumPy: numpy/_core/src/umath/loops_utils.h.src line 63
 *
 * Decreasing it decreases errors slightly as more pairs are summed but
 * also lowers performance. Since the inner loop is unrolled 8 times,
 * the effective blocksize is 16.
 */
#define PW_BLOCKSIZE 128

/*
 * Pairwise summation for float64, rounding error O(lg n) instead of O(n).
 * The recursion depth is O(lg n) as well.
 *
 * Ported from NumPy: numpy/_core/src/umath/loops_utils.h.src lines 80-145
 * (DOUBLE_pairwise_sum variant)
 *
 * @param data   Pointer to contiguous array of doubles
 * @param n      Number of elements
 * @param stride Stride between elements (in number of doubles, not bytes)
 * @return       Sum of all elements
 */
double pairwise_sum_f64(const double* data, size_t n, size_t stride);

/*
 * Pairwise summation for float32
 * Same algorithm as f64 version but for single precision
 */
float pairwise_sum_f32(const float* data, size_t n, size_t stride);

#endif /* NUMJS_PAIRWISE_SUM_H */
