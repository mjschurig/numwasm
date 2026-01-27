/**
 * NumJS Random - BitGenerator Interface
 *
 * Defines the interface for all bit generators (PRNGs).
 * Based on NumPy's bitgen_t structure.
 */

#ifndef NUMJS_BITGEN_H
#define NUMJS_BITGEN_H

#include <stdint.h>

/**
 * BitGenerator interface structure.
 * All bit generators must implement these function pointers.
 */
typedef struct {
    void *state;                                    /* Generator-specific state */
    uint64_t (*next_uint64)(void *state);          /* Generate 64-bit unsigned */
    uint32_t (*next_uint32)(void *state);          /* Generate 32-bit unsigned */
    double (*next_double)(void *state);            /* Generate [0, 1) double */
    uint64_t (*next_raw)(void *state);             /* Raw output (same as uint64) */
} bitgen_t;

/**
 * Convert uint64 to double in [0, 1).
 * Uses upper 53 bits for IEEE 754 double precision.
 *
 * The conversion divides by 2^53 to get the range [0, 1).
 * This gives uniform spacing of ~1.1e-16 between possible values.
 */
static inline double uint64_to_double(uint64_t x) {
    return (x >> 11) * (1.0 / 9007199254740992.0);  /* 2^53 = 9007199254740992 */
}

/**
 * Convert uint32 to float in [0, 1).
 * Uses upper 24 bits for IEEE 754 single precision.
 *
 * Float has 23 mantissa bits + 1 implicit, so we use 24 bits.
 */
static inline float uint32_to_float(uint32_t x) {
    return (x >> 8) * (1.0f / 16777216.0f);  /* 2^24 = 16777216 */
}

/**
 * Convert uint64 to double in (0, 1].
 * Useful for algorithms that need to avoid 0 (e.g., log).
 */
static inline double uint64_to_double_open(uint64_t x) {
    return ((x >> 11) + 1) * (1.0 / 9007199254740992.0);
}

/**
 * Convert uint64 to double in (0, 1).
 * Useful for algorithms that need to avoid both 0 and 1.
 */
static inline double uint64_to_double_open_open(uint64_t x) {
    return ((x >> 12) + 0.5) * (1.0 / 4503599627370496.0);  /* 2^52 */
}

#endif /* NUMJS_BITGEN_H */
