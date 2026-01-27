/**
 * NumJS Random - PCG64 BitGenerator
 *
 * PCG64 is the default random number generator for NumJS.
 * It uses a 128-bit Linear Congruential Generator with the XSL-RR output function.
 *
 * Based on Melissa O'Neill's PCG family of generators.
 * Reference: https://www.pcg-random.org
 *
 * Properties:
 * - Period: 2^128
 * - State: 128 bits (state) + 128 bits (increment/stream)
 * - Output: 64 bits
 * - Statistically excellent (passes BigCrush, PractRand)
 * - Supports jump-ahead in O(log n) time
 */

#ifndef NUMJS_PCG64_H
#define NUMJS_PCG64_H

#include "bitgen.h"
#include <stdint.h>

/**
 * 128-bit unsigned integer for PCG state.
 * Uses struct for platforms without native __uint128_t (like WASM).
 */
typedef struct {
    uint64_t high;
    uint64_t low;
} pcg128_t;

/**
 * PCG64 internal state structure.
 * Uses Linear Congruential Generator formula:
 *   state = state * multiplier + increment
 */
typedef struct {
    pcg128_t state;    /* Current LCG state */
    pcg128_t inc;      /* Increment (stream selector, must be odd) */
} pcg_state_t;

/**
 * PCG64 wrapper with buffered 32-bit output.
 * When generating 32-bit values, we save the upper half
 * of a 64-bit generation for the next call.
 */
typedef struct {
    pcg_state_t *pcg_state;  /* Pointer to internal state */
    int has_uint32;          /* Flag: buffered 32-bit value available */
    uint32_t uinteger;       /* Buffered upper 32 bits */
} pcg64_state;

/* PCG64 multiplier constant (from the original PCG paper) */
#define PCG_DEFAULT_MULTIPLIER_HIGH 2549297995355413924ULL
#define PCG_DEFAULT_MULTIPLIER_LOW  4865540595714422341ULL

/* ============ 128-bit Arithmetic ============ */

/**
 * Multiply two 128-bit integers.
 */
pcg128_t pcg128_mult(pcg128_t a, pcg128_t b);

/**
 * Add two 128-bit integers.
 */
pcg128_t pcg128_add(pcg128_t a, pcg128_t b);

/* ============ PCG64 Functions ============ */

/**
 * Create a new PCG64 state structure.
 * Returns pointer to allocated state, must be freed with pcg64_free().
 */
pcg64_state* pcg64_create(void);

/**
 * Free a PCG64 state structure.
 */
void pcg64_free(pcg64_state *state);

/**
 * Seed the PCG64 generator.
 *
 * @param state     PCG64 state structure
 * @param seed_high Upper 64 bits of seed
 * @param seed_low  Lower 64 bits of seed
 * @param inc_high  Upper 64 bits of increment (stream selector)
 * @param inc_low   Lower 64 bits of increment (will be made odd)
 */
void pcg64_seed(pcg64_state *state,
                uint64_t seed_high, uint64_t seed_low,
                uint64_t inc_high, uint64_t inc_low);

/**
 * Generate a 64-bit unsigned random integer.
 */
uint64_t pcg64_next64(pcg64_state *state);

/**
 * Generate a 32-bit unsigned random integer.
 * Uses buffering to avoid wasting bits.
 */
uint32_t pcg64_next32(pcg64_state *state);

/**
 * Generate a random double in [0, 1).
 */
double pcg64_next_double(pcg64_state *state);

/**
 * Advance the generator by delta steps.
 * Uses O(log n) algorithm.
 *
 * @param state      PCG64 state
 * @param delta_high Upper 64 bits of step count
 * @param delta_low  Lower 64 bits of step count
 */
void pcg64_advance(pcg64_state *state, uint64_t delta_high, uint64_t delta_low);

/**
 * Get the current state as an array of uint64.
 * Output array must have at least 6 elements:
 *   [0]: state.high, [1]: state.low
 *   [2]: inc.high,   [3]: inc.low
 *   [4]: has_uint32, [5]: uinteger
 */
void pcg64_get_state(pcg64_state *state, uint64_t *out);

/**
 * Set the state from an array of uint64.
 * Input array format same as get_state output.
 */
void pcg64_set_state(pcg64_state *state, uint64_t *in);

/**
 * Initialize a bitgen_t interface for PCG64.
 */
void pcg64_init_bitgen(bitgen_t *bitgen, pcg64_state *state);

#endif /* NUMJS_PCG64_H */
