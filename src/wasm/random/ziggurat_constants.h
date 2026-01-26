/**
 * NumJS Random - Ziggurat Algorithm Constants
 *
 * Pre-computed tables for the Ziggurat algorithm used in
 * normal and exponential distribution sampling.
 *
 * The Ziggurat algorithm uses 256 layers (rectangles) to cover the
 * probability density function. ~98.9% of samples return on the
 * first iteration, making it very efficient.
 *
 * Tables are provided for both double (64-bit) and float (32-bit) precision.
 */

#ifndef NUMJS_ZIGGURAT_CONSTANTS_H
#define NUMJS_ZIGGURAT_CONSTANTS_H

#include <stdint.h>

/* ============ Normal Distribution Constants ============ */

/* Tail cutoff for normal distribution */
#define ZIGGURAT_NOR_R       3.6541528853610088
#define ZIGGURAT_NOR_INV_R   0.27366123732975828

/* Float versions */
#define ZIGGURAT_NOR_R_F     3.6541528f
#define ZIGGURAT_NOR_INV_R_F 0.27366123f

/**
 * Normal distribution Ziggurat tables (double precision).
 *
 * ki_double[i]: Threshold - if random < ki_double[i], accept x*wi_double[i]
 * wi_double[i]: Width scale factor
 * fi_double[i]: f(x_i) = exp(-0.5 * x_i^2) at layer boundary
 */
extern const uint64_t ki_double[256];
extern const double wi_double[256];
extern const double fi_double[256];

/**
 * Normal distribution Ziggurat tables (single precision).
 */
extern const uint32_t ki_float[256];
extern const float wi_float[256];
extern const float fi_float[256];

/* ============ Exponential Distribution Constants ============ */

/* Tail cutoff for exponential distribution */
#define ZIGGURAT_EXP_R       7.6971174701310497
#define ZIGGURAT_EXP_R_F     7.6971175f

/**
 * Exponential distribution Ziggurat tables (double precision).
 *
 * ke_double[i]: Threshold for fast acceptance
 * we_double[i]: Width scale factor
 * fe_double[i]: f(x_i) = exp(-x_i) at layer boundary
 */
extern const uint64_t ke_double[256];
extern const double we_double[256];
extern const double fe_double[256];

/**
 * Exponential distribution Ziggurat tables (single precision).
 */
extern const uint32_t ke_float[256];
extern const float we_float[256];
extern const float fe_float[256];

#endif /* NUMJS_ZIGGURAT_CONSTANTS_H */
