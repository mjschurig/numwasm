/**
 * Utility functions for benchmarks.
 * Includes seeded PRNG for reproducible random data generation.
 */

/**
 * Seeded PRNG (Mulberry32) for reproducibility.
 * This ensures consistent random data across runs.
 */
export function createSeededRandom(seed: number): () => number {
  let s = seed;
  return function () {
    let t = (s += 0x6d2b79f5);
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

/**
 * Box-Muller transform for normal distribution.
 * Generates numbers similar to numpy.random.randn().
 */
export function randn(random: () => number): number {
  const u1 = random();
  const u2 = random();
  return Math.sqrt(-2 * Math.log(u1)) * Math.cos(2 * Math.PI * u2);
}

/**
 * Generate an array of random normal values with a specific seed.
 */
export function generateRandomData(size: number, seed: number = 42): number[] {
  const random = createSeededRandom(seed);
  return Array.from({ length: size }, () => randn(random));
}

/**
 * Generate an array of uniform random values [0, 1) with a specific seed.
 */
export function generateUniformData(size: number, seed: number = 42): number[] {
  const random = createSeededRandom(seed);
  return Array.from({ length: size }, () => random());
}

/**
 * Format a size number for display (e.g., 1000000 -> "1M").
 */
export function formatSize(size: number): string {
  if (size >= 1_000_000) {
    return `${size / 1_000_000}M`;
  } else if (size >= 1_000) {
    return `${size / 1_000}K`;
  }
  return size.toString();
}

/**
 * Format milliseconds for display.
 */
export function formatMs(ms: number): string {
  if (ms < 0.001) {
    return `${(ms * 1000).toFixed(3)} μs`;
  } else if (ms < 1) {
    return `${ms.toFixed(4)} ms`;
  } else {
    return `${ms.toFixed(2)} ms`;
  }
}
