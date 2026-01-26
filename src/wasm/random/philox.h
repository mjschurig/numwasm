/**
 * NumJS Random - Philox BitGenerator
 *
 * Philox (Philox4x64-10) counter-based pseudo-random number generator.
 * A fast, parallelizable PRNG well-suited for GPU and distributed computing.
 *
 * Properties:
 * - Period: 2^256
 * - State: 4 × 64-bit counter + 2 × 64-bit key
 * - Output: 4 × 64-bit values per counter increment
 * - Supports efficient jump and advance operations
 *
 * Reference: John K. Salmon et al., "Parallel Random Numbers: As Easy as 1, 2, 3",
 * SC'11 Proceedings, 2011
 */

#ifndef NUMJS_PHILOX_H
#define NUMJS_PHILOX_H

#include "bitgen.h"
#include <stdint.h>

/* Philox constants */
#define PHILOX_ROUNDS 10
#define PHILOX_BUFFER_SIZE 4

/* Multiplication constants for Philox4x64 */
#define PHILOX_M0 0xD2E7470EE14C6C93ULL
#define PHILOX_M1 0xCA5A826395121157ULL

/* Weyl sequence constants (key bump) */
#define PHILOX_W0 0x9E3779B97F4A7C15ULL
#define PHILOX_W1 0xBB67AE8584CAA73BULL

/**
 * Philox 4x64 counter structure.
 */
typedef struct {
    uint64_t v[4];
} philox_ctr_t;

/**
 * Philox 2x64 key structure.
 */
typedef struct {
    uint64_t v[2];
} philox_key_t;

/**
 * Philox state structure with output buffer.
 */
typedef struct {
    philox_ctr_t ctr;                    /* 4x64-bit counter */
    philox_key_t key;                    /* 2x64-bit key */
    int buffer_pos;                      /* Position in output buffer (0-4) */
    uint64_t buffer[PHILOX_BUFFER_SIZE]; /* Output buffer */
    int has_uint32;                      /* Buffered 32-bit value flag */
    uint32_t uinteger;                   /* Buffered 32-bit value */
} philox_state;

/* ============ Philox Functions ============ */

/**
 * Create a new Philox state structure.
 * Returns pointer to allocated state, must be freed with philox_free().
 */
philox_state* philox_create(void);

/**
 * Free a Philox state structure.
 */
void philox_free(philox_state *state);

/**
 * Seed the Philox generator with a 2x64-bit key.
 * Counter is reset to zero.
 *
 * @param state  Philox state structure
 * @param key0   First 64 bits of key
 * @param key1   Second 64 bits of key
 */
void philox_seed(philox_state *state, uint64_t key0, uint64_t key1);

/**
 * Seed using 32-bit parts (for JavaScript compatibility).
 *
 * @param state  Philox state structure
 * @param parts  Array of 4 uint32 values representing key (key0_hi, key0_lo, key1_hi, key1_lo)
 */
void philox_seed_parts(philox_state *state, uint32_t *parts);

/**
 * Set both counter and key directly.
 *
 * @param state   Philox state structure
 * @param counter Array of 4 uint64 values for counter
 * @param key     Array of 2 uint64 values for key
 */
void philox_set_counter(philox_state *state, uint64_t *counter, uint64_t *key);

/**
 * Generate a 64-bit unsigned random integer.
 */
uint64_t philox_next64(philox_state *state);

/**
 * Generate a 32-bit unsigned random integer.
 * Uses buffering to avoid wasting bits.
 */
uint32_t philox_next32(philox_state *state);

/**
 * Generate a random double in [0, 1).
 */
double philox_next_double(philox_state *state);

/**
 * Get next 64-bit value and return parts (for JavaScript).
 *
 * @param state    Philox state
 * @param high_out Output for upper 32 bits
 * @return Lower 32 bits
 */
uint32_t philox_next64_parts(philox_state *state, uint32_t *high_out);

/**
 * Jump ahead by 2^128 draws.
 * This is equivalent to incrementing counter[2] by 1.
 */
void philox_jump(philox_state *state);

/**
 * Advance the counter by an arbitrary number of steps.
 * Each step generates 4 values, so advancing by n increments counter by n.
 *
 * @param state Philox state
 * @param step  Array of 4 uint64 values representing 256-bit step
 */
void philox_advance(philox_state *state, uint64_t *step);

/**
 * Advance by a single 64-bit step count (for convenience).
 */
void philox_advance64(philox_state *state, uint64_t step);

/**
 * Get the current state.
 *
 * @param state      Philox state
 * @param ctr_out    Array of 4 uint64 for counter
 * @param key_out    Array of 2 uint64 for key
 * @param buffer_pos_out  Output for buffer position
 * @param buffer_out Array of 4 uint64 for buffer
 * @param has_uint32_out  Output for has_uint32 flag
 * @param uinteger_out    Output for buffered value
 */
void philox_get_state(philox_state *state,
                       uint64_t *ctr_out,
                       uint64_t *key_out,
                       int *buffer_pos_out,
                       uint64_t *buffer_out,
                       int *has_uint32_out,
                       uint32_t *uinteger_out);

/**
 * Set the state.
 */
void philox_set_state(philox_state *state,
                       uint64_t *ctr_in,
                       uint64_t *key_in,
                       int buffer_pos,
                       uint64_t *buffer_in,
                       int has_uint32,
                       uint32_t uinteger);

/**
 * Initialize a bitgen_t interface for Philox.
 */
void philox_init_bitgen(bitgen_t *bitgen, philox_state *state);

/**
 * Fill an array with random uint64 values.
 */
void philox_fill_uint64(philox_state *state, uint64_t *out, int64_t count);

/**
 * Fill an array with random doubles in [0, 1).
 */
void philox_fill_double(philox_state *state, double *out, int64_t count);

#endif /* NUMJS_PHILOX_H */
