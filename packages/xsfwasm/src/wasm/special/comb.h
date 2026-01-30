#ifndef SCIWASM_SPECIAL_COMB_H
#define SCIWASM_SPECIAL_COMB_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Binomial coefficient (non-exact mode): C(n, k) = n! / (k! * (n-k)!)
 * Uses lgamma for floating-point computation.
 *
 * @param n Number of items
 * @param k Number of items to choose
 * @return Binomial coefficient as double
 */
double wasm_binom(double n, double k);

/**
 * Binomial coefficient (exact mode): C(n, k) using integer arithmetic
 * Uses unsigned long long with overflow detection.
 *
 * @param n Number of items (must be integer)
 * @param k Number of items to choose (must be integer)
 * @return Binomial coefficient as double, or -1.0 on overflow/error
 */
double wasm_binom_exact(double n, double k);

/**
 * Pochhammer symbol (rising factorial): (x)_m = x * (x+1) * ... * (x+m-1)
 * Used for computing permutations in non-exact mode.
 *
 * @param x Starting value
 * @param m Number of terms
 * @return Pochhammer symbol as double
 */
double wasm_poch(double x, double m);

/**
 * Permutations (exact mode): P(n, k) = n! / (n-k)!
 * Uses unsigned long long with overflow detection.
 *
 * @param n Number of items (must be integer)
 * @param k Number of items to arrange (must be integer)
 * @return Number of permutations as double, or -1.0 on overflow/error
 */
double wasm_perm_exact(double n, double k);

#ifdef __cplusplus
}
#endif

#endif // SCIWASM_SPECIAL_COMB_H
