/**
 * NumJS Random - SFC64 BitGenerator
 *
 * SFC64 (Small Fast Chaotic 64) pseudo-random number generator.
 * A fast generator with a 256-bit state (4 × 64-bit values).
 *
 * Properties:
 * - Period: approximately 2^255
 * - State: 256 bits (a, b, c, counter)
 * - Output: 64 bits
 * - Very fast generation with good statistical properties
 *
 * Reference: PractRand by Chris Doty-Humphrey
 */

#ifndef NUMJS_SFC64_H
#define NUMJS_SFC64_H

#include "bitgen.h"
#include <stdint.h>

/**
 * SFC64 internal state structure.
 * Uses 4 × 64-bit values: a, b, c, and a Weyl counter.
 */
typedef struct {
    uint64_t a;
    uint64_t b;
    uint64_t c;
    uint64_t w;  /* Weyl sequence counter */
} sfc64_internal_state;

/**
 * SFC64 wrapper with buffered 32-bit output.
 */
typedef struct {
    sfc64_internal_state *state;
    int has_uint32;      /* Flag: buffered 32-bit value available */
    uint32_t uinteger;   /* Buffered upper 32 bits */
} sfc64_state;

/* ============ SFC64 Functions ============ */

/**
 * Create a new SFC64 state structure.
 * Returns pointer to allocated state, must be freed with sfc64_free().
 */
sfc64_state* sfc64_create(void);

/**
 * Free a SFC64 state structure.
 */
void sfc64_free(sfc64_state *state);

/**
 * Seed the SFC64 generator.
 *
 * @param state  SFC64 state structure
 * @param seed   Array of 4 uint64 values (a, b, c, w)
 */
void sfc64_seed(sfc64_state *state, uint64_t *seed);

/**
 * Seed using 32-bit parts (for JavaScript compatibility).
 *
 * @param state  SFC64 state structure
 * @param parts  Array of 8 uint32 values representing 4 uint64 values
 */
void sfc64_seed_parts(sfc64_state *state, uint32_t *parts);

/**
 * Generate a 64-bit unsigned random integer.
 */
uint64_t sfc64_next64(sfc64_state *state);

/**
 * Generate a 32-bit unsigned random integer.
 * Uses buffering to avoid wasting bits.
 */
uint32_t sfc64_next32(sfc64_state *state);

/**
 * Generate a random double in [0, 1).
 */
double sfc64_next_double(sfc64_state *state);

/**
 * Get next 64-bit value and return parts (for JavaScript).
 *
 * @param state    SFC64 state
 * @param high_out Output for upper 32 bits
 * @return Lower 32 bits
 */
uint32_t sfc64_next64_parts(sfc64_state *state, uint32_t *high_out);

/**
 * Get the current state as an array of uint64.
 * Output array must have at least 6 elements:
 *   [0]: a, [1]: b, [2]: c, [3]: w, [4]: has_uint32, [5]: uinteger
 */
void sfc64_get_state(sfc64_state *state, uint64_t *out);

/**
 * Set the state from an array of uint64.
 * Input array format same as get_state output.
 */
void sfc64_set_state(sfc64_state *state, uint64_t *in);

/**
 * Initialize a bitgen_t interface for SFC64.
 */
void sfc64_init_bitgen(bitgen_t *bitgen, sfc64_state *state);

/**
 * Fill an array with random uint64 values.
 */
void sfc64_fill_uint64(sfc64_state *state, uint64_t *out, int64_t count);

/**
 * Fill an array with random doubles in [0, 1).
 */
void sfc64_fill_double(sfc64_state *state, double *out, int64_t count);

#endif /* NUMJS_SFC64_H */
