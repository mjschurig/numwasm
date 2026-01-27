/*
 * Pairwise Summation Implementation
 *
 * Direct port of NumPy's pairwise summation algorithm.
 * Original source: numpy/_core/src/umath/loops_utils.h.src lines 80-145
 *
 * Key features:
 * - O(lg n) rounding error instead of O(n) for naive summation
 * - 8-way unrolling for better CPU pipelining and vectorization
 * - Recursive divide-and-conquer for large arrays
 */

#include "pairwise_sum.h"

/*
 * Pairwise summation for float64
 *
 * The algorithm works in three cases:
 * 1. n < 8: Simple loop (base case)
 * 2. 8 <= n <= 128: 8-accumulator unrolled loop with tree reduction
 * 3. n > 128: Recursive divide and conquer
 */
double pairwise_sum_f64(const double* data, size_t n, size_t stride)
{
    if (n < 8) {
        /*
         * Base case: simple summation for small arrays.
         * Start with -0.0 to preserve -0 values.
         * Reason: summing only -0 should return -0, but 0 + -0 == 0
         * while -0 + -0 == -0.
         * (NumPy comment from line 86-88)
         */
        double res = -0.0;
        for (size_t i = 0; i < n; i++) {
            res += data[i * stride];
        }
        return res;
    }
    else if (n <= PW_BLOCKSIZE) {
        /*
         * Block case: 8-accumulator unrolled loop.
         * 8 times unroll reduces blocksize to 16 and allows vectorization
         * with AVX without changing summation ordering.
         * (NumPy comment from lines 100-104)
         */
        double r[8];
        double res;
        size_t i;

        /* Initialize 8 accumulators (lines 105-112) */
        r[0] = data[0 * stride];
        r[1] = data[1 * stride];
        r[2] = data[2 * stride];
        r[3] = data[3 * stride];
        r[4] = data[4 * stride];
        r[5] = data[5 * stride];
        r[6] = data[6 * stride];
        r[7] = data[7 * stride];

        /* 8-way unrolled accumulation (lines 114-125) */
        for (i = 8; i < n - (n % 8); i += 8) {
            r[0] += data[(i + 0) * stride];
            r[1] += data[(i + 1) * stride];
            r[2] += data[(i + 2) * stride];
            r[3] += data[(i + 3) * stride];
            r[4] += data[(i + 4) * stride];
            r[5] += data[(i + 5) * stride];
            r[6] += data[(i + 6) * stride];
            r[7] += data[(i + 7) * stride];
        }

        /*
         * Accumulate now to avoid stack spills for single peel loop.
         * Tree-structured accumulation minimizes rounding error.
         * (NumPy lines 127-129)
         */
        res = ((r[0] + r[1]) + (r[2] + r[3])) +
              ((r[4] + r[5]) + (r[6] + r[7]));

        /* Handle remainder elements (lines 131-134) */
        for (; i < n; i++) {
            res += data[i * stride];
        }
        return res;
    }
    else {
        /*
         * Recursive case: divide and conquer.
         * Divide by two but avoid non-multiples of unroll factor.
         * (NumPy lines 137-144)
         */
        size_t n2 = n / 2;
        n2 -= n2 % 8;
        return pairwise_sum_f64(data, n2, stride) +
               pairwise_sum_f64(data + n2 * stride, n - n2, stride);
    }
}

/*
 * Pairwise summation for float32
 * Same algorithm as f64 but for single precision floats.
 */
float pairwise_sum_f32(const float* data, size_t n, size_t stride)
{
    if (n < 8) {
        float res = -0.0f;
        for (size_t i = 0; i < n; i++) {
            res += data[i * stride];
        }
        return res;
    }
    else if (n <= PW_BLOCKSIZE) {
        float r[8];
        float res;
        size_t i;

        r[0] = data[0 * stride];
        r[1] = data[1 * stride];
        r[2] = data[2 * stride];
        r[3] = data[3 * stride];
        r[4] = data[4 * stride];
        r[5] = data[5 * stride];
        r[6] = data[6 * stride];
        r[7] = data[7 * stride];

        for (i = 8; i < n - (n % 8); i += 8) {
            r[0] += data[(i + 0) * stride];
            r[1] += data[(i + 1) * stride];
            r[2] += data[(i + 2) * stride];
            r[3] += data[(i + 3) * stride];
            r[4] += data[(i + 4) * stride];
            r[5] += data[(i + 5) * stride];
            r[6] += data[(i + 6) * stride];
            r[7] += data[(i + 7) * stride];
        }

        res = ((r[0] + r[1]) + (r[2] + r[3])) +
              ((r[4] + r[5]) + (r[6] + r[7]));

        for (; i < n; i++) {
            res += data[i * stride];
        }
        return res;
    }
    else {
        size_t n2 = n / 2;
        n2 -= n2 % 8;
        return pairwise_sum_f32(data, n2, stride) +
               pairwise_sum_f32(data + n2 * stride, n - n2, stride);
    }
}
