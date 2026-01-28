/**
 * Tests for perm() function
 * Based on scipy tests from test_basic.py:1605-1629
 */

import { describe, it, expect, beforeAll } from 'vitest';
import { perm } from '../../../src/ts/special/perm.js';
import { loadWasmModule } from '../../../src/ts/wasm-loader.js';
import { NDArray } from 'numwasm';

describe('perm', () => {
  beforeAll(async () => {
    await loadWasmModule();
  });

  describe('non-exact mode', () => {
    it('should handle scalar inputs', async () => {
      const result = await perm(10, 3);
      expect(result).toBeCloseTo(720, 1);
    });

    it('should handle array inputs', async () => {
      const N = await NDArray.array([10, 10]);
      const k = await NDArray.array([3, 4]);
      const result = await perm(N, k);

      expect(result).toBeInstanceOf(NDArray);
      const data = await (result as NDArray).getData();
      expect(data[0]).toBeCloseTo(720, 1);
      expect(data[1]).toBeCloseTo(5040, 1);
    });

    it('should return 0 for k > N', async () => {
      const result = await perm(2, 3);
      expect(result).toBeCloseTo(0, 1);
    });

    it('should return 0 for negative inputs', async () => {
      const result1 = await perm(-1, 3);
      const result2 = await perm(2, -1);
      expect(result1).toBeCloseTo(0, 1);
      expect(result2).toBeCloseTo(0, 1);
    });

    it('should handle mixed valid/invalid array inputs', async () => {
      const N = await NDArray.array([2, -1, 2, 10]);
      const k = await NDArray.array([3, 3, -1, 3]);
      const result = await perm(N, k);

      const data = await (result as NDArray).getData();
      expect(data[0]).toBeCloseTo(0, 1);  // k > N
      expect(data[1]).toBeCloseTo(0, 1);  // N < 0
      expect(data[2]).toBeCloseTo(0, 1);  // k < 0
      expect(data[3]).toBeCloseTo(720, 1);
    });
  });

  describe('exact mode', () => {
    it('should compute exact permutations', async () => {
      const result = await perm(10, 3, { exact: true });
      expect(result).toBe(720);
    });

    it('should handle small exact values', async () => {
      expect(await perm(10, 3, { exact: true })).toBe(720);
      expect(await perm(10, 4, { exact: true })).toBe(5040);
    });

    it('should return 0 for k > N in exact mode', async () => {
      expect(await perm(2, 3, { exact: true })).toBe(0);
    });

    it('should return 0 for negative inputs in exact mode', async () => {
      expect(await perm(-1, 3, { exact: true })).toBe(0);
      expect(await perm(2, -1, { exact: true })).toBe(0);
    });

    it('should throw error for non-integer inputs', async () => {
      await expect(perm(4.6, 6, { exact: true })).rejects.toThrow(/Non-integer/);
      await expect(perm(-4.6, 3, { exact: true })).rejects.toThrow(/Non-integer/);
      await expect(perm(4, -3.9, { exact: true })).rejects.toThrow(/Non-integer/);
      await expect(perm(6.0, 4.6, { exact: true })).rejects.toThrow(/Non-integer/);
    });

    it('should throw error for array inputs in exact mode', async () => {
      const N = await NDArray.array([1, 2]);
      const k = await NDArray.array([4, 5]);
      await expect(perm(N, k, { exact: true })).rejects.toThrow(/scalar integers/);
    });
  });
});
