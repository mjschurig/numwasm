/**
 * NumJS Random - Distribution Algorithms
 *
 * Implements various probability distributions using the BitGenerator interface.
 * Algorithms include:
 * - Ziggurat for Normal and Exponential distributions
 * - Marsaglia-Tsang for Gamma distribution
 * - Johnk's algorithm for Beta distribution
 * - BTPE for Binomial distribution
 * - PTRS for Poisson distribution
 * - Lemire's algorithm for bounded integers
 */

#ifndef NUMJS_DISTRIBUTIONS_H
#define NUMJS_DISTRIBUTIONS_H

#include "bitgen.h"
#include <stdint.h>

/* ============ Uniform Distributions ============ */

/**
 * Generate uniform random double in [0, 1).
 */
double random_standard_uniform(bitgen_t *state);

/**
 * Generate uniform random float in [0, 1).
 */
float random_standard_uniform_f(bitgen_t *state);

/**
 * Fill array with uniform random doubles in [0, 1).
 */
void random_uniform_fill(bitgen_t *state, int32_t cnt, double *out);

/**
 * Fill array with uniform random floats in [0, 1).
 */
void random_uniform_fill_f(bitgen_t *state, int32_t cnt, float *out);

/* ============ Normal Distribution (Ziggurat) ============ */

/**
 * Generate standard normal variate N(0, 1) using Ziggurat algorithm.
 * ~98.9% of samples return on the first iteration.
 */
double random_standard_normal(bitgen_t *state);

/**
 * Generate standard normal variate (single precision).
 */
float random_standard_normal_f(bitgen_t *state);

/**
 * Fill array with standard normal variates.
 */
void random_standard_normal_fill(bitgen_t *state, int32_t cnt, double *out);

/**
 * Fill array with standard normal variates (single precision).
 */
void random_standard_normal_fill_f(bitgen_t *state, int32_t cnt, float *out);

/**
 * Generate normal variate with given location and scale.
 * Returns: loc + scale * N(0, 1)
 */
double random_normal(bitgen_t *state, double loc, double scale);

/* ============ Exponential Distribution (Ziggurat) ============ */

/**
 * Generate standard exponential variate Exp(1) using Ziggurat algorithm.
 */
double random_standard_exponential(bitgen_t *state);

/**
 * Generate standard exponential variate (single precision).
 */
float random_standard_exponential_f(bitgen_t *state);

/**
 * Fill array with standard exponential variates.
 */
void random_standard_exponential_fill(bitgen_t *state, int32_t cnt, double *out);

/**
 * Fill array with standard exponential variates (single precision).
 */
void random_standard_exponential_fill_f(bitgen_t *state, int32_t cnt, float *out);

/**
 * Fill array using inverse CDF method (alternative to Ziggurat).
 */
void random_standard_exponential_inv_fill(bitgen_t *state, int32_t cnt, double *out);

/**
 * Generate exponential variate with given scale.
 * Returns: scale * Exp(1)
 */
double random_exponential(bitgen_t *state, double scale);

/* ============ Gamma Distribution ============ */

/**
 * Generate standard gamma variate using Marsaglia-Tsang method.
 * For shape < 1, uses Ahrens-Dieter transformation.
 * For shape >= 1, uses Marsaglia-Tsang squeeze method.
 */
double random_standard_gamma(bitgen_t *state, double shape);

/**
 * Generate standard gamma variate (single precision).
 */
float random_standard_gamma_f(bitgen_t *state, float shape);

/**
 * Fill array with standard gamma variates.
 */
void random_standard_gamma_fill(bitgen_t *state, int32_t cnt, double shape, double *out);

/**
 * Generate gamma variate with given shape and scale.
 * Returns: scale * Gamma(shape)
 */
double random_gamma(bitgen_t *state, double shape, double scale);

/* ============ Beta Distribution ============ */

/**
 * Generate beta variate using appropriate algorithm:
 * - Johnk's algorithm for a <= 1 and b <= 1
 * - Gamma ratio method otherwise
 */
double random_beta(bitgen_t *state, double a, double b);

/**
 * Fill array with beta variates.
 */
void random_beta_fill(bitgen_t *state, int32_t cnt, double a, double b, double *out);

/* ============ Chi-Square Distribution ============ */

/**
 * Generate chi-square variate with given degrees of freedom.
 * Equivalent to 2 * Gamma(df/2).
 */
double random_chisquare(bitgen_t *state, double df);

/**
 * Generate noncentral chi-square variate.
 */
double random_noncentral_chisquare(bitgen_t *state, double df, double nonc);

/* ============ F Distribution ============ */

/**
 * Generate F-distributed variate.
 */
double random_f(bitgen_t *state, double dfnum, double dfden);

/**
 * Generate noncentral F-distributed variate.
 */
double random_noncentral_f(bitgen_t *state, double dfnum, double dfden, double nonc);

/* ============ Student's t Distribution ============ */

/**
 * Generate Student's t-distributed variate.
 */
double random_standard_t(bitgen_t *state, double df);

/* ============ Other Continuous Distributions ============ */

/**
 * Generate standard Cauchy variate (Lorentz distribution).
 */
double random_standard_cauchy(bitgen_t *state);

/**
 * Generate Pareto II (Lomax) variate with shape parameter a.
 */
double random_pareto(bitgen_t *state, double a);

/**
 * Generate Weibull variate with shape parameter a.
 */
double random_weibull(bitgen_t *state, double a);

/**
 * Generate power-distributed variate with exponent a.
 */
double random_power(bitgen_t *state, double a);

/**
 * Generate Laplace (double exponential) variate.
 */
double random_laplace(bitgen_t *state, double loc, double scale);

/**
 * Generate Gumbel (extreme value type I) variate.
 */
double random_gumbel(bitgen_t *state, double loc, double scale);

/**
 * Generate logistic-distributed variate.
 */
double random_logistic(bitgen_t *state, double loc, double scale);

/**
 * Generate log-normal variate.
 * Returns exp(N(mean, sigma)).
 */
double random_lognormal(bitgen_t *state, double mean, double sigma);

/**
 * Generate Rayleigh-distributed variate.
 */
double random_rayleigh(bitgen_t *state, double scale);

/**
 * Generate Wald (inverse Gaussian) variate.
 */
double random_wald(bitgen_t *state, double mean, double scale);

/**
 * Generate triangular-distributed variate.
 */
double random_triangular(bitgen_t *state, double left, double mode, double right);

/**
 * Generate von Mises distributed variate (circular).
 */
double random_vonmises(bitgen_t *state, double mu, double kappa);

/* ============ Discrete Distributions ============ */

/**
 * State structure for binomial distribution.
 * Caches values for repeated calls with same parameters.
 */
typedef struct {
    int initialized;
    int64_t nsave;
    double psave;
    double r;
    double q;
    double fm;
    int64_t m;
    double p1;
    double xm;
    double xl;
    double xr;
    double c;
    double laml;
    double lamr;
    double p2;
    double p3;
    double p4;
} binomial_t;

/**
 * Generate binomial variate using BTPE algorithm for large n*p,
 * or inversion method for small n*p.
 */
int64_t random_binomial(bitgen_t *state, double p, int64_t n);

/**
 * Generate binomial variate using BTPE (BinomialTriangleParetoExponential).
 * More efficient for large n*min(p, 1-p).
 */
int64_t random_binomial_btpe(bitgen_t *state, binomial_t *binomial, int64_t n, double p);

/**
 * Generate binomial variate using inversion method.
 * Efficient for small n*min(p, 1-p).
 */
int64_t random_binomial_inversion(bitgen_t *state, int64_t n, double p);

/**
 * Generate negative binomial variate.
 */
int64_t random_negative_binomial(bitgen_t *state, double n, double p);

/**
 * Generate Poisson variate using PTRS algorithm for large lambda,
 * or direct method for small lambda.
 */
int64_t random_poisson(bitgen_t *state, double lam);

/**
 * Generate geometric variate.
 */
int64_t random_geometric(bitgen_t *state, double p);

/**
 * Generate hypergeometric variate.
 */
int64_t random_hypergeometric(bitgen_t *state, int64_t ngood, int64_t nbad, int64_t nsample);

/**
 * Generate logarithmic series variate.
 */
int64_t random_logseries(bitgen_t *state, double p);

/**
 * Generate Zipf (zeta) variate.
 */
int64_t random_zipf(bitgen_t *state, double a);

/* ============ Bounded Integers ============ */

/**
 * Generate unbiased bounded uint64 using Lemire's algorithm.
 * Returns a value in [off, off + rng].
 * mask is a bitmask >= rng for optimization.
 */
uint64_t random_bounded_uint64(bitgen_t *state, uint64_t off, uint64_t rng, uint64_t mask);

/**
 * Generate unbiased bounded uint32 using Lemire's algorithm.
 */
uint32_t random_bounded_uint32(bitgen_t *state, uint32_t off, uint32_t rng, uint32_t mask);

/**
 * Generate random integer in [low, high].
 */
int64_t random_integers(bitgen_t *state, int64_t low, int64_t high);

/**
 * Fill array with random integers in [low, high].
 */
void random_integers_fill(bitgen_t *state, int32_t cnt, int64_t low, int64_t high, int64_t *out);

#endif /* NUMJS_DISTRIBUTIONS_H */
