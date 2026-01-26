/**
 * Random Module Tests (Phase 15)
 *
 * Tests for numpy.random compatible random number generation.
 */

import { describe, it, expect, beforeAll } from 'vitest';
import {
  NDArray,
  loadWasmModule,
  Generator,
  PCG64,
  SeedSequence,
  default_rng,
  random,
  seed,
  randn,
  randint,
  DType,
} from '../../dist/numjs.mjs';

// Helper to check if value is in range [min, max)
function inRange(val: number, min: number, max: number): boolean {
  return val >= min && val < max;
}

// Helper to compute mean of array
function computeMean(arr: NDArray | number[]): number {
  if (Array.isArray(arr)) {
    return arr.reduce((a, b) => a + b, 0) / arr.length;
  }
  let sum = 0;
  for (let i = 0; i < arr.size; i++) {
    sum += arr.get(i) as number;
  }
  return sum / arr.size;
}

// Helper to compute variance of array
function computeVariance(arr: NDArray | number[], mean?: number): number {
  const m = mean ?? computeMean(arr);
  let sumSq = 0;
  const size = Array.isArray(arr) ? arr.length : arr.size;
  for (let i = 0; i < size; i++) {
    const val = Array.isArray(arr) ? arr[i] : arr.get(i) as number;
    sumSq += (val - m) ** 2;
  }
  return sumSq / size;
}

// Helper for statistical tests - check if value is within expected range
function withinTolerance(actual: number, expected: number, tolerance: number): boolean {
  return Math.abs(actual - expected) <= tolerance;
}

describe('Random Module', () => {
  beforeAll(async () => {
    await loadWasmModule();
  });

  describe('SeedSequence', () => {
    it('creates from integer seed', () => {
      const ss = new SeedSequence(12345);
      expect(ss).toBeDefined();
      expect(ss.entropy.length).toBeGreaterThan(0);
    });

    it('creates from array seed', () => {
      const ss = new SeedSequence([1, 2, 3, 4]);
      expect(ss).toBeDefined();
    });

    it('creates with OS entropy when no seed', () => {
      const ss = new SeedSequence();
      expect(ss).toBeDefined();
      expect(ss.entropy.length).toBeGreaterThan(0);
    });

    it('generates state words', () => {
      const ss = new SeedSequence(12345);
      const state = ss.generateState(4);
      expect(state).toBeInstanceOf(Uint32Array);
      expect(state.length).toBe(4);
    });

    it('spawns independent sequences', () => {
      const ss = new SeedSequence(12345);
      const [child1, child2] = ss.spawn(2);

      expect(child1).toBeInstanceOf(SeedSequence);
      expect(child2).toBeInstanceOf(SeedSequence);

      // Children should have different spawn keys
      expect(child1.spawnKey).not.toEqual(child2.spawnKey);
    });

    it('produces reproducible state', () => {
      const ss1 = new SeedSequence(12345);
      const ss2 = new SeedSequence(12345);

      const state1 = ss1.generateState(4);
      const state2 = ss2.generateState(4);

      expect(Array.from(state1)).toEqual(Array.from(state2));
    });
  });

  describe('PCG64 BitGenerator', () => {
    it('creates with integer seed', () => {
      const pcg = new PCG64(12345);
      expect(pcg).toBeDefined();
      pcg.dispose();
    });

    it('creates with SeedSequence', () => {
      const ss = new SeedSequence(12345);
      const pcg = new PCG64(ss);
      expect(pcg).toBeDefined();
      pcg.dispose();
    });

    it('generates uint64 values', () => {
      const pcg = new PCG64(12345);
      const val = pcg.next_uint64();
      expect(typeof val).toBe('bigint');
      pcg.dispose();
    });

    it('generates uint32 values', () => {
      const pcg = new PCG64(12345);
      const val = pcg.next_uint32();
      expect(typeof val).toBe('number');
      expect(val).toBeGreaterThanOrEqual(0);
      expect(val).toBeLessThanOrEqual(0xFFFFFFFF);
      pcg.dispose();
    });

    it('generates double values in [0, 1)', () => {
      const pcg = new PCG64(12345);
      for (let i = 0; i < 100; i++) {
        const val = pcg.next_double();
        expect(val).toBeGreaterThanOrEqual(0);
        expect(val).toBeLessThan(1);
      }
      pcg.dispose();
    });

    it('produces reproducible sequences', () => {
      const pcg1 = new PCG64(12345);
      const pcg2 = new PCG64(12345);

      for (let i = 0; i < 10; i++) {
        expect(pcg1.next_uint64()).toBe(pcg2.next_uint64());
      }

      pcg1.dispose();
      pcg2.dispose();
    });

    it('can get and set state', () => {
      const pcg1 = new PCG64(12345);

      // Advance some steps
      for (let i = 0; i < 100; i++) {
        pcg1.next_uint64();
      }

      // Save state
      const state = pcg1.getState();
      expect(state.bit_generator).toBe('PCG64');

      // Continue and record values
      const values1: bigint[] = [];
      for (let i = 0; i < 10; i++) {
        values1.push(pcg1.next_uint64());
      }

      // Create new generator and restore state
      const pcg2 = new PCG64(0);
      pcg2.setState(state);

      // Should produce same values
      for (let i = 0; i < 10; i++) {
        expect(pcg2.next_uint64()).toBe(values1[i]);
      }

      pcg1.dispose();
      pcg2.dispose();
    });

    it('spawns independent generators', () => {
      const pcg = new PCG64(12345);
      const [child1, child2] = pcg.spawn(2);

      // Children should produce different sequences
      const val1 = child1.next_uint64();
      const val2 = child2.next_uint64();
      expect(val1).not.toBe(val2);

      pcg.dispose();
      child1.dispose();
      child2.dispose();
    });
  });

  describe('Generator class', () => {
    describe('random()', () => {
      it('returns scalar without size argument', () => {
        const rng = default_rng(12345);
        const val = rng.random();
        expect(typeof val).toBe('number');
        expect(val).toBeGreaterThanOrEqual(0);
        expect(val).toBeLessThan(1);
        rng.dispose();
      });

      it('returns array with size argument', async () => {
        const rng = default_rng(12345);
        const arr = rng.random(10) as NDArray;
        expect(arr.shape).toEqual([10]);
        for (let i = 0; i < 10; i++) {
          const val = arr.get(i) as number;
          expect(val).toBeGreaterThanOrEqual(0);
          expect(val).toBeLessThan(1);
        }
        arr.dispose();
        rng.dispose();
      });

      it('respects shape argument', async () => {
        const rng = default_rng(12345);
        const arr = rng.random([3, 4]) as NDArray;
        expect(arr.shape).toEqual([3, 4]);
        arr.dispose();
        rng.dispose();
      });
    });

    describe('integers()', () => {
      it('generates integers in range', () => {
        const rng = default_rng(12345);
        for (let i = 0; i < 100; i++) {
          const val = rng.integers(0, 10) as number;
          expect(Number.isInteger(val)).toBe(true);
          expect(val).toBeGreaterThanOrEqual(0);
          expect(val).toBeLessThan(10);
        }
        rng.dispose();
      });

      it('handles single argument as upper bound', () => {
        const rng = default_rng(12345);
        for (let i = 0; i < 100; i++) {
          const val = rng.integers(5) as number;
          expect(val).toBeGreaterThanOrEqual(0);
          expect(val).toBeLessThan(5);
        }
        rng.dispose();
      });

      it('returns array with size argument', async () => {
        const rng = default_rng(12345);
        const arr = rng.integers(0, 100, 20) as NDArray;
        expect(arr.shape).toEqual([20]);
        for (let i = 0; i < 20; i++) {
          const val = arr.get(i) as number;
          expect(val).toBeGreaterThanOrEqual(0);
          expect(val).toBeLessThan(100);
        }
        arr.dispose();
        rng.dispose();
      });

      it('supports endpoint option', () => {
        const rng = default_rng(12345);
        let foundMax = false;
        for (let i = 0; i < 1000; i++) {
          const val = rng.integers(0, 5, null, DType.Int64, true) as number;
          if (val === 5) foundMax = true;
        }
        expect(foundMax).toBe(true);
        rng.dispose();
      });
    });

    describe('uniform()', () => {
      it('generates values in specified range', () => {
        const rng = default_rng(12345);
        for (let i = 0; i < 100; i++) {
          const val = rng.uniform(10, 20) as number;
          expect(val).toBeGreaterThanOrEqual(10);
          expect(val).toBeLessThan(20);
        }
        rng.dispose();
      });

      it('returns array with size argument', async () => {
        const rng = default_rng(12345);
        const arr = rng.uniform(-1, 1, 50) as NDArray;
        expect(arr.shape).toEqual([50]);
        for (let i = 0; i < 50; i++) {
          const val = arr.get(i) as number;
          expect(val).toBeGreaterThanOrEqual(-1);
          expect(val).toBeLessThan(1);
        }
        arr.dispose();
        rng.dispose();
      });
    });

    describe('standard_normal()', () => {
      it('generates values from standard normal distribution', async () => {
        const rng = default_rng(12345);
        const arr = rng.standard_normal(10000) as NDArray;

        // Check statistical properties (mean ~0, variance ~1)
        const mean = computeMean(arr);
        const variance = computeVariance(arr, mean);

        expect(withinTolerance(mean, 0, 0.1)).toBe(true);
        expect(withinTolerance(variance, 1, 0.1)).toBe(true);

        arr.dispose();
        rng.dispose();
      });
    });

    describe('normal()', () => {
      it('generates values with specified mean and scale', async () => {
        const rng = default_rng(12345);
        const arr = rng.normal(5, 2, 10000) as NDArray;

        const mean = computeMean(arr);
        const variance = computeVariance(arr, mean);

        expect(withinTolerance(mean, 5, 0.2)).toBe(true);
        expect(withinTolerance(variance, 4, 0.4)).toBe(true);

        arr.dispose();
        rng.dispose();
      });

      it('throws on negative scale', () => {
        const rng = default_rng(12345);
        expect(() => rng.normal(0, -1)).toThrow('scale must be non-negative');
        rng.dispose();
      });
    });

    describe('exponential()', () => {
      it('generates non-negative values', async () => {
        const rng = default_rng(12345);
        const arr = rng.exponential(1, 1000) as NDArray;

        for (let i = 0; i < arr.size; i++) {
          expect(arr.get(i)).toBeGreaterThanOrEqual(0);
        }

        // Mean should be ~scale
        const mean = computeMean(arr);
        expect(withinTolerance(mean, 1, 0.2)).toBe(true);

        arr.dispose();
        rng.dispose();
      });

      it('respects scale parameter', async () => {
        const rng = default_rng(12345);
        const arr = rng.exponential(5, 5000) as NDArray;

        const mean = computeMean(arr);
        expect(withinTolerance(mean, 5, 0.5)).toBe(true);

        arr.dispose();
        rng.dispose();
      });
    });

    describe('gamma()', () => {
      it('generates non-negative values', async () => {
        const rng = default_rng(12345);
        const arr = rng.gamma(2, 1, 1000) as NDArray;

        for (let i = 0; i < arr.size; i++) {
          expect(arr.get(i)).toBeGreaterThanOrEqual(0);
        }

        // Mean should be shape * scale
        const mean = computeMean(arr);
        expect(withinTolerance(mean, 2, 0.3)).toBe(true);

        arr.dispose();
        rng.dispose();
      });
    });

    describe('beta()', () => {
      it('generates values in [0, 1]', async () => {
        const rng = default_rng(12345);
        const arr = rng.beta(2, 5, 1000) as NDArray;

        for (let i = 0; i < arr.size; i++) {
          const val = arr.get(i) as number;
          expect(val).toBeGreaterThanOrEqual(0);
          expect(val).toBeLessThanOrEqual(1);
        }

        // Mean should be a / (a + b) = 2/7 ~ 0.286
        const mean = computeMean(arr);
        expect(withinTolerance(mean, 2 / 7, 0.05)).toBe(true);

        arr.dispose();
        rng.dispose();
      });

      it('throws on non-positive parameters', () => {
        const rng = default_rng(12345);
        expect(() => rng.beta(0, 1)).toThrow('a must be positive');
        expect(() => rng.beta(1, 0)).toThrow('b must be positive');
        rng.dispose();
      });
    });

    describe('binomial()', () => {
      it('generates integers in [0, n]', async () => {
        const rng = default_rng(12345);
        const arr = rng.binomial(10, 0.5, 1000) as NDArray;

        for (let i = 0; i < arr.size; i++) {
          const val = arr.get(i) as number;
          expect(Number.isInteger(val)).toBe(true);
          expect(val).toBeGreaterThanOrEqual(0);
          expect(val).toBeLessThanOrEqual(10);
        }

        // Mean should be n * p = 5
        const mean = computeMean(arr);
        expect(withinTolerance(mean, 5, 0.3)).toBe(true);

        arr.dispose();
        rng.dispose();
      });
    });

    describe('poisson()', () => {
      it('generates non-negative integers', async () => {
        const rng = default_rng(12345);
        const arr = rng.poisson(5, 1000) as NDArray;

        for (let i = 0; i < arr.size; i++) {
          const val = arr.get(i) as number;
          expect(Number.isInteger(val)).toBe(true);
          expect(val).toBeGreaterThanOrEqual(0);
        }

        // Mean should be ~lambda
        const mean = computeMean(arr);
        expect(withinTolerance(mean, 5, 0.3)).toBe(true);

        arr.dispose();
        rng.dispose();
      });
    });

    describe('geometric()', () => {
      it('generates positive integers', async () => {
        const rng = default_rng(12345);
        const arr = rng.geometric(0.5, 1000) as NDArray;

        for (let i = 0; i < arr.size; i++) {
          const val = arr.get(i) as number;
          expect(Number.isInteger(val)).toBe(true);
          expect(val).toBeGreaterThanOrEqual(1);
        }

        // Mean should be 1/p = 2
        const mean = computeMean(arr);
        expect(withinTolerance(mean, 2, 0.3)).toBe(true);

        arr.dispose();
        rng.dispose();
      });
    });

    describe('bytes()', () => {
      it('returns Uint8Array of specified length', () => {
        const rng = default_rng(12345);
        const bytes = rng.bytes(16);
        expect(bytes).toBeInstanceOf(Uint8Array);
        expect(bytes.length).toBe(16);
        rng.dispose();
      });
    });

    describe('spawn()', () => {
      it('creates independent generators', () => {
        const rng = default_rng(12345);
        const [child1, child2] = rng.spawn(2);

        // Children should produce different sequences
        const val1 = child1.random();
        const val2 = child2.random();
        expect(val1).not.toBe(val2);

        rng.dispose();
        child1.dispose();
        child2.dispose();
      });
    });
  });

  describe('default_rng() factory', () => {
    it('creates generator with integer seed', () => {
      const rng = default_rng(12345);
      expect(rng).toBeInstanceOf(Generator);
      rng.dispose();
    });

    it('creates generator with SeedSequence', () => {
      const ss = new SeedSequence(12345);
      const rng = default_rng(ss);
      expect(rng).toBeInstanceOf(Generator);
      rng.dispose();
    });

    it('creates generator with PCG64', () => {
      const pcg = new PCG64(12345);
      const rng = default_rng(pcg);
      expect(rng).toBeInstanceOf(Generator);
      rng.dispose();
    });

    it('creates generator with OS entropy when no seed', () => {
      const rng = default_rng();
      expect(rng).toBeInstanceOf(Generator);
      rng.dispose();
    });
  });

  describe('Module-level functions', () => {
    beforeAll(() => {
      seed(12345);
    });

    describe('seed()', () => {
      it('makes random() reproducible', () => {
        seed(12345);
        const val1 = random();
        seed(12345);
        const val2 = random();
        expect(val1).toBe(val2);
      });
    });

    describe('random()', () => {
      it('returns scalar without size', () => {
        const val = random();
        expect(typeof val).toBe('number');
        expect(val).toBeGreaterThanOrEqual(0);
        expect(val).toBeLessThan(1);
      });

      it('returns array with size', async () => {
        const arr = random(10) as NDArray;
        expect(arr.shape).toEqual([10]);
        arr.dispose();
      });
    });

    describe('randn()', () => {
      it('returns scalar without arguments', () => {
        const val = randn();
        expect(typeof val).toBe('number');
      });

      it('returns array with shape arguments', async () => {
        const arr = randn(3, 4) as NDArray;
        expect(arr.shape).toEqual([3, 4]);
        arr.dispose();
      });
    });

    describe('randint()', () => {
      it('generates integers in range', () => {
        for (let i = 0; i < 100; i++) {
          const val = randint(0, 10) as number;
          expect(Number.isInteger(val)).toBe(true);
          expect(val).toBeGreaterThanOrEqual(0);
          expect(val).toBeLessThan(10);
        }
      });
    });
  });

  describe('Reproducibility', () => {
    it('same seed produces same sequence', async () => {
      const rng1 = default_rng(42);
      const rng2 = default_rng(42);

      // Uniform
      expect(rng1.random()).toBe(rng2.random());

      // Normal
      expect(rng1.standard_normal()).toBe(rng2.standard_normal());

      // Integers
      expect(rng1.integers(0, 1000)).toBe(rng2.integers(0, 1000));

      // Arrays
      const arr1 = rng1.random(10) as NDArray;
      const arr2 = rng2.random(10) as NDArray;
      for (let i = 0; i < 10; i++) {
        expect(arr1.get(i)).toBe(arr2.get(i));
      }

      arr1.dispose();
      arr2.dispose();
      rng1.dispose();
      rng2.dispose();
    });

    it('different seeds produce different sequences', () => {
      const rng1 = default_rng(12345);
      const rng2 = default_rng(54321);

      expect(rng1.random()).not.toBe(rng2.random());

      rng1.dispose();
      rng2.dispose();
    });
  });
});
