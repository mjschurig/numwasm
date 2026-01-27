/**
 * NumJS Random - SeedSequence
 *
 * Provides reproducible entropy management for seeding BitGenerators.
 * Based on Melissa E. O'Neill's C++11 seed_seq design.
 *
 * Key features:
 * - Mixes entropy from multiple sources
 * - Supports spawn() for creating independent child sequences
 * - Generates high-quality seed state for any BitGenerator
 */

#ifndef NUMJS_SEEDSEQ_H
#define NUMJS_SEEDSEQ_H

#include <stdint.h>

/**
 * SeedSequence state structure.
 */
typedef struct {
    uint32_t *pool;       /* Entropy pool */
    int32_t pool_size;    /* Pool size in uint32 words (default 4 = 128 bits) */
    int32_t n_children;   /* Number of children spawned */
} seed_sequence_t;

/* Mixing constants (from NumPy) */
#define SEED_SEQ_INIT_A  0x43b0d7e5u
#define SEED_SEQ_MULT_A  0x931e8875u
#define SEED_SEQ_INIT_B  0x8b51f9ddu
#define SEED_SEQ_MULT_B  0x58f38dadu

/**
 * Initialize a SeedSequence with entropy and optional spawn key.
 *
 * @param seq          SeedSequence structure to initialize
 * @param entropy      Array of entropy values
 * @param entropy_len  Number of entropy values
 * @param spawn_key    Array of spawn key values (for child sequences)
 * @param spawn_key_len Number of spawn key values
 * @param pool_size    Size of entropy pool in uint32 words (default 4)
 */
void seed_seq_init(seed_sequence_t *seq,
                   const uint32_t *entropy, int32_t entropy_len,
                   const uint32_t *spawn_key, int32_t spawn_key_len,
                   int32_t pool_size);

/**
 * Generate seed state words.
 *
 * @param seq      Initialized SeedSequence
 * @param out      Output array for generated words
 * @param n_words  Number of words to generate
 */
void seed_seq_generate(seed_sequence_t *seq, uint32_t *out, int32_t n_words);

/**
 * Free resources allocated by seed_seq_init.
 */
void seed_seq_free(seed_sequence_t *seq);

/**
 * Convenience function to generate words from entropy directly.
 * Combines init, generate, and free in one call.
 */
void seed_seq_generate_words(const uint32_t *entropy, int32_t entropy_len,
                             const uint32_t *spawn_key, int32_t spawn_key_len,
                             int32_t pool_size,
                             uint32_t *out, int32_t n_words);

#endif /* NUMJS_SEEDSEQ_H */
