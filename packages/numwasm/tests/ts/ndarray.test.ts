/**
 * NDArray Unit Tests
 *
 * Tests for the core NDArray class functionality.
 */

import { describe, it, expect, beforeAll, afterEach } from 'vitest';
import { NDArray, DType, loadWasmModule } from 'numwasm';

describe('NDArray', () => {
  // Load WASM module before all tests
  beforeAll(async () => {
    await loadWasmModule();
  });

  // Track arrays for cleanup
  const arrays: NDArray[] = [];
  afterEach(() => {
    for (const arr of arrays) {
      if (!arr.isDisposed) {
        arr.dispose();
      }
    }
    arrays.length = 0;
  });

  const track = (arr: NDArray) => {
    arrays.push(arr);
    return arr;
  };

  describe('zeros()', () => {
    it('should create a 1D array of zeros', async () => {
      const arr = track(await NDArray.zeros([5]));
      expect(arr.shape).toEqual([5]);
      expect(arr.ndim).toBe(1);
      expect(arr.size).toBe(5);
      expect(arr.toArray()).toEqual([0, 0, 0, 0, 0]);
    });

    it('should create a 2D array of zeros', async () => {
      const arr = track(await NDArray.zeros([3, 4]));
      expect(arr.shape).toEqual([3, 4]);
      expect(arr.ndim).toBe(2);
      expect(arr.size).toBe(12);
      expect(arr.toArray().every((v) => v === 0)).toBe(true);
    });

    it('should create a 3D array of zeros', async () => {
      const arr = track(await NDArray.zeros([2, 3, 4]));
      expect(arr.shape).toEqual([2, 3, 4]);
      expect(arr.ndim).toBe(3);
      expect(arr.size).toBe(24);
    });

    it('should support Float32 dtype', async () => {
      const arr = track(await NDArray.zeros([5], { dtype: DType.Float32 }));
      expect(arr.dtype).toBe(DType.Float32);
      expect(arr.toArray()).toEqual([0, 0, 0, 0, 0]);
    });

    it('should support Int32 dtype', async () => {
      const arr = track(await NDArray.zeros([5], { dtype: DType.Int32 }));
      expect(arr.dtype).toBe(DType.Int32);
      expect(arr.toArray()).toEqual([0, 0, 0, 0, 0]);
    });

    it('should throw on empty shape', async () => {
      await expect(NDArray.zeros([])).rejects.toThrow();
    });

    it('should throw on negative dimension', async () => {
      await expect(NDArray.zeros([-1])).rejects.toThrow();
    });
  });

  describe('ones()', () => {
    it('should create an array of ones', async () => {
      const arr = track(await NDArray.ones([5]));
      expect(arr.toArray()).toEqual([1, 1, 1, 1, 1]);
    });

    it('should create a 2D array of ones', async () => {
      const arr = track(await NDArray.ones([2, 3]));
      expect(arr.size).toBe(6);
      expect(arr.toArray().every((v) => v === 1)).toBe(true);
    });
  });

  describe('fromArray()', () => {
    it('should create from number array', async () => {
      const data = [1, 2, 3, 4, 5];
      const arr = track(await NDArray.fromArray(data));
      expect(arr.toArray()).toEqual(data);
      expect(arr.shape).toEqual([5]);
    });

    it('should create from Float64Array', async () => {
      const data = new Float64Array([1.5, 2.5, 3.5]);
      const arr = track(await NDArray.fromArray(data));
      expect(arr.toArray()).toEqual([1.5, 2.5, 3.5]);
    });

    it('should create with explicit shape', async () => {
      const arr = track(await NDArray.fromArray([1, 2, 3, 4, 5, 6], [2, 3]));
      expect(arr.shape).toEqual([2, 3]);
      expect(arr.size).toBe(6);
    });

    it('should throw on shape mismatch', async () => {
      await expect(NDArray.fromArray([1, 2, 3], [2, 2])).rejects.toThrow();
    });

    it('should preserve floating point values', async () => {
      const arr = track(await NDArray.fromArray([1.1, 2.2, 3.3]));
      const result = arr.toArray();
      expect(result[0]).toBeCloseTo(1.1, 10);
      expect(result[1]).toBeCloseTo(2.2, 10);
      expect(result[2]).toBeCloseTo(3.3, 10);
    });
  });

  describe('arange()', () => {
    it('should create range with single argument', async () => {
      const arr = track(await NDArray.arange(5));
      expect(arr.toArray()).toEqual([0, 1, 2, 3, 4]);
    });

    it('should create range with start and end', async () => {
      const arr = track(await NDArray.arange(2, 7));
      expect(arr.toArray()).toEqual([2, 3, 4, 5, 6]);
    });

    it('should create range with step', async () => {
      const arr = track(await NDArray.arange(0, 10, 2));
      expect(arr.toArray()).toEqual([0, 2, 4, 6, 8]);
    });

    it('should handle negative step', async () => {
      const arr = track(await NDArray.arange(5, 0, -1));
      expect(arr.toArray()).toEqual([5, 4, 3, 2, 1]);
    });

    it('should throw on zero step', async () => {
      await expect(NDArray.arange(0, 5, 0)).rejects.toThrow();
    });
  });

  describe('sum()', () => {
    it('should sum simple array', async () => {
      const arr = track(await NDArray.fromArray([1, 2, 3, 4, 5]));
      expect(arr.sum()).toBe(15);
    });

    it('should sum floating point array', async () => {
      const arr = track(await NDArray.fromArray([1.5, 2.5, 3.0]));
      expect(arr.sum()).toBeCloseTo(7.0, 10);
    });

    it('should sum zeros', async () => {
      const arr = track(await NDArray.zeros([100]));
      expect(arr.sum()).toBe(0);
    });

    it('should sum ones', async () => {
      const arr = track(await NDArray.ones([100]));
      expect(arr.sum()).toBe(100);
    });

    it('should handle negative values', async () => {
      const arr = track(await NDArray.fromArray([-1, -2, 3, 4]));
      expect(arr.sum()).toBe(4);
    });

    it('should use pairwise summation for large arrays', async () => {
      // Sum of 0..999 = 999 * 1000 / 2 = 499500
      const data = Array.from({ length: 1000 }, (_, i) => i);
      const arr = track(await NDArray.fromArray(data));
      expect(arr.sum()).toBe(499500);
    });

    it('should handle array larger than PW_BLOCKSIZE (128)', async () => {
      const data = Array.from({ length: 200 }, (_, i) => i + 1);
      const arr = track(await NDArray.fromArray(data));
      // Sum of 1..200 = 200 * 201 / 2 = 20100
      expect(arr.sum()).toBe(20100);
    });
  });

  describe('fill()', () => {
    it('should fill array with value', async () => {
      const arr = track(await NDArray.zeros([5]));
      arr.fill(42);
      expect(arr.toArray()).toEqual([42, 42, 42, 42, 42]);
    });

    it('should fill with floating point value', async () => {
      const arr = track(await NDArray.zeros([3]));
      arr.fill(3.14);
      const result = arr.toArray();
      expect(result[0]).toBeCloseTo(3.14, 10);
      expect(result[1]).toBeCloseTo(3.14, 10);
      expect(result[2]).toBeCloseTo(3.14, 10);
    });
  });

  describe('properties', () => {
    it('should return correct shape', async () => {
      const arr = track(await NDArray.zeros([3, 4, 5]));
      expect(arr.shape).toEqual([3, 4, 5]);
    });

    it('should return correct ndim', async () => {
      const arr = track(await NDArray.zeros([3, 4, 5]));
      expect(arr.ndim).toBe(3);
    });

    it('should return correct size', async () => {
      const arr = track(await NDArray.zeros([3, 4, 5]));
      expect(arr.size).toBe(60);
    });

    it('should return correct dtype', async () => {
      const arr = track(await NDArray.zeros([5], { dtype: DType.Float32 }));
      expect(arr.dtype).toBe(DType.Float32);
    });
  });

  describe('memory management', () => {
    it('should track disposed state', async () => {
      const arr = await NDArray.zeros([5]);
      expect(arr.isDisposed).toBe(false);
      arr.dispose();
      expect(arr.isDisposed).toBe(true);
    });

    it('should throw when accessing disposed array', async () => {
      const arr = await NDArray.zeros([5]);
      arr.dispose();
      expect(() => arr.sum()).toThrow('disposed');
      expect(() => arr.toArray()).toThrow('disposed');
      expect(() => arr.fill(1)).toThrow('disposed');
    });

    it('should be safe to dispose twice', async () => {
      const arr = await NDArray.zeros([5]);
      arr.dispose();
      expect(() => arr.dispose()).not.toThrow();
    });
  });
});
