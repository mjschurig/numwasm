/**
 * NumJS Random - MT19937 BitGenerator
 *
 * MT19937 (Mersenne Twister) pseudo-random number generator.
 * The classic PRNG with period 2^19937 - 1.
 *
 * Properties:
 * - Period: 2^19937 - 1
 * - State: 624 × 32-bit words + position index
 * - Output: 32 bits (native), 64 bits via combining
 * - Widely used for compatibility with legacy systems
 *
 * Reference: Makoto Matsumoto and Takuji Nishimura,
 * "Mersenne Twister: A 623-Dimensionally Equidistributed Uniform
 * Pseudo-Random Number Generator", ACM TOMACS 1998
 */

#ifndef NUMJS_MT19937_H
#define NUMJS_MT19937_H

#include "bitgen.h"
#include <stdint.h>

/* MT19937 constants */
#define MT19937_N 624         /* State array size */
#define MT19937_M 397         /* Middle point for twist */
#define MT19937_MATRIX_A 0x9908b0dfUL    /* Constant vector a */
#define MT19937_UPPER_MASK 0x80000000UL  /* Most significant w-r bits */
#define MT19937_LOWER_MASK 0x7fffffffUL  /* Least significant r bits */

/**
 * MT19937 internal state structure.
 */
typedef struct {
    uint32_t key[MT19937_N];  /* State array */
    int pos;                  /* Current position (0 to N) */
} mt19937_internal_state;

/**
 * MT19937 wrapper with buffered 32-bit output.
 * Since MT19937 is natively 32-bit, buffering works opposite to PCG64:
 * we buffer 32 bits when generating 64-bit values.
 */
typedef struct {
    mt19937_internal_state *state;
    int has_uint32;      /* Flag: buffered 32-bit value available */
    uint32_t uinteger;   /* Buffered 32-bit value */
} mt19937_state;

/* ============ MT19937 Functions ============ */

/**
 * Create a new MT19937 state structure.
 * Returns pointer to allocated state, must be freed with mt19937_free().
 */
mt19937_state* mt19937_create(void);

/**
 * Free a MT19937 state structure.
 */
void mt19937_free(mt19937_state *state);

/**
 * Seed the MT19937 generator with a single 32-bit value.
 *
 * @param state  MT19937 state structure
 * @param seed   32-bit seed value
 */
void mt19937_seed(mt19937_state *state, uint32_t seed);

/**
 * Seed the MT19937 generator with an array of 32-bit values.
 * This is the init_by_array initialization from the reference implementation.
 *
 * @param state      MT19937 state structure
 * @param init_key   Array of seed values
 * @param key_length Length of init_key array
 */
void mt19937_seed_array(mt19937_state *state, uint32_t *init_key, int key_length);

/**
 * Generate a 32-bit unsigned random integer.
 * This is the native output size for MT19937.
 */
uint32_t mt19937_next32(mt19937_state *state);

/**
 * Generate a 64-bit unsigned random integer.
 * Combines two 32-bit values.
 */
uint64_t mt19937_next64(mt19937_state *state);

/**
 * Generate a random double in [0, 1).
 * Uses the MT19937 convention: (a * 67108864.0 + b) / 9007199254740992.0
 * where a and b are 27 and 26 bits respectively.
 */
double mt19937_next_double(mt19937_state *state);

/**
 * Get next 64-bit value and return parts (for JavaScript).
 *
 * @param state    MT19937 state
 * @param high_out Output for upper 32 bits
 * @return Lower 32 bits
 */
uint32_t mt19937_next64_parts(mt19937_state *state, uint32_t *high_out);

/**
 * Get the current state.
 * Output array must have space for 624 uint32 values (key) + 1 int (pos)
 * + 1 int (has_uint32) + 1 uint32 (uinteger).
 *
 * @param state   MT19937 state
 * @param key_out Array of at least 624 uint32 for the key
 * @param pos_out Pointer to store position
 * @param has_uint32_out Pointer to store has_uint32 flag
 * @param uinteger_out Pointer to store buffered value
 */
void mt19937_get_state(mt19937_state *state, uint32_t *key_out,
                        int *pos_out, int *has_uint32_out, uint32_t *uinteger_out);

/**
 * Set the state.
 *
 * @param state   MT19937 state
 * @param key_in  Array of 624 uint32 values
 * @param pos     Position value
 * @param has_uint32 Buffered flag
 * @param uinteger Buffered value
 */
void mt19937_set_state(mt19937_state *state, uint32_t *key_in,
                        int pos, int has_uint32, uint32_t uinteger);

/**
 * Initialize a bitgen_t interface for MT19937.
 */
void mt19937_init_bitgen(bitgen_t *bitgen, mt19937_state *state);

/**
 * Fill an array with random uint32 values.
 */
void mt19937_fill_uint32(mt19937_state *state, uint32_t *out, int64_t count);

/**
 * Fill an array with random uint64 values.
 */
void mt19937_fill_uint64(mt19937_state *state, uint64_t *out, int64_t count);

/**
 * Fill an array with random doubles in [0, 1).
 */
void mt19937_fill_double(mt19937_state *state, double *out, int64_t count);

#endif /* NUMJS_MT19937_H */
