/**
 * Tests for factorial functions
 * Based on scipy.special.tests.test_basic.py
 */
import { describe, it, expect } from 'vitest';
import { special } from '../../../src/ts/index.js';

describe('special.factorial', () => {
  it('should compute factorials accurately for n=0 to 20', async () => {
    const expected = [
      1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880,
      3628800, 39916800, 479001600, 6227020800, 87178291200,
      1307674368000, 20922789888000, 355687428096000,
      6402373705728000, 121645100408832000, 2432902008176640000
    ];

    for (let n = 0; n <= 20; n++) {
      const result = await special.factorial(n, { exact: false });
      expect(result).toBeCloseTo(expected[n], 10);
    }
  });

  it('should compute exact integer factorials', async () => {
    const result5 = await special.factorial(5, { exact: true });
    expect(result5).toBe(120);

    const result10 = await special.factorial(10, { exact: true });
    expect(result10).toBe(3628800);
  });

  it('should handle large factorials with BigInt', async () => {
    const result25 = await special.factorial(25, { exact: true });
    expect(result25).toBe(15511210043330985984000000n);
  });

  it('should return 0 for negative inputs with extend="zero"', async () => {
    const result = await special.factorial(-1, { exact: false });
    expect(result).toBe(0);

    const resultExact = await special.factorial(-5, { exact: true });
    expect(resultExact).toBe(0);
  });

  it('should return scalar for scalar input', async () => {
    const result = await special.factorial(5, { exact: false });
    expect(typeof result).toBe('number');
    expect(Array.isArray(result)).toBe(false);
  });

  it('should handle array inputs', async () => {
    const result = await special.factorial([3, 4, 5], { exact: false });
    expect(Array.isArray(result)).toBe(true);
    expect(result).toHaveLength(3);
    expect(result[0]).toBeCloseTo(6, 10);
    expect(result[1]).toBeCloseTo(24, 10);
    expect(result[2]).toBeCloseTo(120, 10);
  });

  it('should handle array inputs with exact mode', async () => {
    const result = await special.factorial([3, 4, 5], { exact: true });
    expect(Array.isArray(result)).toBe(true);
    expect(result[0]).toBe(6);
    expect(result[1]).toBe(24);
    expect(result[2]).toBe(120);
  });
});

describe('special.factorial2', () => {
  it('should compute double factorial correctly', async () => {
    // 7!! = 7 * 5 * 3 * 1 = 105
    const result7 = await special.factorial2(7, { exact: true });
    expect(result7).toBe(105);

    // 8!! = 8 * 6 * 4 * 2 = 384
    const result8 = await special.factorial2(8, { exact: true });
    expect(result8).toBe(384);
  });

  it('should compute double factorial approximately', async () => {
    const result = await special.factorial2(7, { exact: false });
    expect(result).toBeCloseTo(105, 8);
  });

  it('should handle array inputs', async () => {
    const result = await special.factorial2([3, 5, 7], { exact: true });
    expect(result[0]).toBe(3);
    expect(result[1]).toBe(15);
    expect(result[2]).toBe(105);
  });
});

describe('special.factorialk', () => {
  it('should compute multifactorial correctly', async () => {
    // 17!!!! = 17 * 13 * 9 * 5 * 1 = 8721
    const result = await special.factorialk(17, 4, { exact: true });
    expect(result).toBe(8721);
  });

  it('should match factorial when k=1', async () => {
    const fact5 = await special.factorial(5, { exact: true });
    const factk5 = await special.factorialk(5, 1, { exact: true });
    expect(factk5).toBe(fact5);
  });

  it('should match factorial2 when k=2', async () => {
    const fact27 = await special.factorial2(7, { exact: true });
    const factk7 = await special.factorialk(7, 2, { exact: true });
    expect(factk7).toBe(fact27);
  });

  it('should handle array inputs', async () => {
    const result = await special.factorialk([5, 7, 9], 3, { exact: true });
    expect(result[0]).toBe(10);
    expect(result[1]).toBe(28);
    expect(result[2]).toBe(162);
  });

  it('should throw error for invalid k', async () => {
    await expect(special.factorialk(5, 0)).rejects.toThrow('k must be a positive integer');
    await expect(special.factorialk(5, -1)).rejects.toThrow('k must be a positive integer');
  });
});
