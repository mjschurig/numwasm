/**
 * Tests for comb() function
 * Based on scipy tests from test_basic.py:1550-1604
 */

import { describe, it, expect, beforeAll } from 'vitest';
import { comb } from '../../../src/ts/special/comb.js';
import { loadWasmModule } from '../../../src/ts/wasm-loader.js';
import { NDArray } from 'numwasm';

describe('comb', () => {
  beforeAll(async () => {
    await loadWasmModule();
  });

  describe('non-exact mode', () => {
    it('should handle scalar inputs', async () => {
      const result = await comb(10, 3);
      expect(result).toBeCloseTo(120, 1);
    });

    it('should handle array inputs', async () => {
      const N = await NDArray.array([10, 10]);
      const k = await NDArray.array([3, 4]);
      const result = await comb(N, k);

      expect(result).toBeInstanceOf(NDArray);
      const data = await (result as NDArray).getData();
      expect(data[0]).toBeCloseTo(120, 1);
      expect(data[1]).toBeCloseTo(210, 1);
    });

    it('should return 0 for k > N', async () => {
      const result = await comb(2, 3);
      expect(result).toBe(0);
    });

    it('should return 0 for negative inputs', async () => {
      const result1 = await comb(-1, 3);
      const result2 = await comb(2, -1);
      expect(result1).toBe(0);
      expect(result2).toBe(0);
    });

    it('should handle mixed valid/invalid array inputs', async () => {
      const N = await NDArray.array([2, -1, 2, 10]);
      const k = await NDArray.array([3, 3, -1, 3]);
      const result = await comb(N, k);

      const data = await (result as NDArray).getData();
      expect(data[0]).toBeCloseTo(0, 1);  // k > N
      expect(data[1]).toBeCloseTo(0, 1);  // N < 0
      expect(data[2]).toBeCloseTo(0, 1);  // k < 0
      expect(data[3]).toBeCloseTo(120, 1);
    });
  });

  describe('exact mode', () => {
    it('should compute exact binomial coefficients', async () => {
      const result = await comb(10, 3, { exact: true });
      expect(result).toBe(120);
    });

    it('should handle small exact values', async () => {
      expect(await comb(10, 3, { exact: true })).toBe(120);
      expect(await comb(20, 10, { exact: true })).toBe(184756);
    });

    it('should return 0 for k > N in exact mode', async () => {
      expect(await comb(2, 3, { exact: true })).toBe(0);
    });

    it('should return 0 for negative inputs in exact mode', async () => {
      expect(await comb(-1, 3, { exact: true })).toBe(0);
      expect(await comb(2, -1, { exact: true })).toBe(0);
    });

    it('should throw error for non-integer inputs', async () => {
      await expect(comb(3.4, 4, { exact: true })).rejects.toThrow(/Non-integer/);
      await expect(comb(3, 4.4, { exact: true })).rejects.toThrow(/Non-integer/);
    });

    it('should throw error for array inputs in exact mode', async () => {
      const N = await NDArray.array([10, 10]);
      const k = await NDArray.array([3, 4]);
      await expect(comb(N, k, { exact: true })).rejects.toThrow(/scalar integers/);
    });
  });

  describe('repetition mode', () => {
    it('should compute combinations with repetition', async () => {
      const result = await comb(10, 3, { exact: true, repetition: true });
      expect(result).toBe(220);  // C(10+3-1, 3) = C(12, 3) = 220
    });

    it('should handle k=0 with repetition', async () => {
      expect(await comb(0, 0, { exact: true, repetition: true })).toBe(1);
      expect(await comb(5, 0, { exact: true, repetition: true })).toBe(1);
      expect(await comb(10, 0, { exact: true, repetition: true })).toBe(1);
    });

    it('should handle k=0 with repetition in non-exact mode', async () => {
      const N = await NDArray.array([0, 5, 10]);
      const k = await NDArray.array([0, 0, 0]);
      const result = await comb(N, k, { repetition: true });

      const data = await (result as NDArray).getData();
      expect(data[0]).toBeCloseTo(1.0, 1);
      expect(data[1]).toBeCloseTo(1.0, 1);
      expect(data[2]).toBeCloseTo(1.0, 1);
    });
  });
});
