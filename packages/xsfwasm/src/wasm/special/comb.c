/**
 * Combinatorial functions for sciwasm: binomial coefficients and permutations
 * Based on scipy.special implementations
 */

#include "comb.h"
#include <math.h>
#include <stdint.h>
#include <limits.h>

// lgamma is available from standard C math library

/**
 * Binomial coefficient (non-exact mode) using lgamma formula:
 * C(n, k) = exp(lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1))
 */
double wasm_binom(double n, double k) {
    // Handle edge cases
    if (k > n || k < 0 || n < 0) {
        return 0.0;
    }

    if (k == 0 || k == n) {
        return 1.0;
    }

    // Use symmetry for optimization: C(n, k) = C(n, n-k)
    if (k > n - k) {
        k = n - k;
    }

    // Compute using lgamma formula
    // lgamma(x) = log(|Gamma(x)|)
    double result = lgamma(n + 1.0) - lgamma(k + 1.0) - lgamma(n - k + 1.0);
    return exp(result);
}

/**
 * Binomial coefficient (exact mode) using integer arithmetic
 * Based on scipy's _comb_int_long implementation
 * Returns -1.0 to signal overflow or invalid input
 */
double wasm_binom_exact(double n, double k) {
    // Convert to unsigned long long
    unsigned long long N = (unsigned long long)n;
    unsigned long long K = (unsigned long long)k;

    // Validate that inputs are actually integers
    if (n != (double)N || k != (double)K) {
        return -1.0;  // Non-integer input
    }

    // Handle negative inputs
    if (n < 0 || k < 0) {
        return 0.0;
    }

    // Handle edge cases
    if (K > N) {
        return 0.0;
    }

    if (K == 0 || K == N) {
        return 1.0;
    }

    // Optimize: use min(k, n-k) to reduce iterations
    unsigned long long nterms = K;
    if (K > N - K) {
        nterms = N - K;
    }

    unsigned long long M = N + 1;
    unsigned long long val = 1;

    // Interleave multiplication and division to reduce overflow
    // Compute: (N * (N-1) * ... * (N-nterms+1)) / (1 * 2 * ... * nterms)
    for (unsigned long long j = 1; j <= nterms; j++) {
        // Check for overflow before multiplication
        if (val > ULLONG_MAX / (M - j)) {
            return -1.0;  // Signal overflow
        }

        val *= (M - j);
        val /= j;
    }

    return (double)val;
}

/**
 * Pochhammer symbol (rising factorial): (x)_m = x * (x+1) * ... * (x+m-1)
 * Computed using gamma ratio: poch(x, m) = Gamma(x+m) / Gamma(x)
 * With lgamma: poch(x, m) = exp(lgamma(x+m) - lgamma(x))
 */
double wasm_poch(double x, double m) {
    // Handle edge cases
    if (m == 0) {
        return 1.0;
    }

    if (m < 0) {
        return 0.0;  // Undefined for negative m
    }

    // Use gamma ratio formula with lgamma
    double result = lgamma(x + m) - lgamma(x);
    return exp(result);
}

/**
 * Permutations (exact mode): P(n, k) = n! / (n-k)! = n * (n-1) * ... * (n-k+1)
 * Based on scipy's perm exact implementation
 * Returns -1.0 to signal overflow or invalid input
 */
double wasm_perm_exact(double n, double k) {
    // Convert to long long for validation
    long long N = (long long)n;
    long long K = (long long)k;

    // Validate that inputs are actually integers
    if (n != (double)N || k != (double)K) {
        return -1.0;  // Non-integer input
    }

    // Handle edge cases
    if (K > N || N < 0 || K < 0) {
        return 0.0;
    }

    if (K == 0) {
        return 1.0;
    }

    unsigned long long val = 1;

    // Direct multiplication: N * (N-1) * ... * (N-K+1)
    for (long long i = N - K + 1; i <= N; i++) {
        // Check for overflow before multiplication
        if (val > ULLONG_MAX / (unsigned long long)i) {
            return -1.0;  // Signal overflow
        }
        val *= (unsigned long long)i;
    }

    return (double)val;
}
