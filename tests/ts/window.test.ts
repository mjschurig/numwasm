/**
 * Tests for Window Functions (Phase 11)
 *
 * Reference: NumPy window functions in numpy/lib/_function_base_impl.py
 */

import { describe, it, expect, beforeAll, afterEach } from 'vitest';
import {
  loadWasmModule,
  NDArray,
  blackman,
  bartlett,
  hanning,
  hamming,
  kaiser,
  i0,
} from '../../dist/numjs.mjs';

// Track resources for cleanup
const resources: NDArray[] = [];
function track<T extends NDArray>(arr: T): T {
  resources.push(arr);
  return arr;
}

// Helper to compare arrays with tolerance
function expectClose(
  actual: number[],
  expected: number[],
  rtol = 1e-7,
  atol = 1e-10
) {
  expect(actual.length).toBe(expected.length);
  for (let i = 0; i < actual.length; i++) {
    const diff = Math.abs(actual[i] - expected[i]);
    const tol = atol + rtol * Math.abs(expected[i]);
    expect(diff).toBeLessThanOrEqual(tol);
  }
}

describe('Window Function Tests', () => {
  beforeAll(async () => {
    await loadWasmModule();
  });

  afterEach(() => {
    for (const arr of resources) {
      if (!arr.isDisposed) {
        arr.dispose();
      }
    }
    resources.length = 0;
  });

  // =========================================================================
  // i0 (Modified Bessel function) Tests
  // =========================================================================
  describe('i0 - Modified Bessel function', () => {
    it('returns 1.0 for i0(0)', async () => {
      const result = track(await i0(0));
      expect(result.item()).toBeCloseTo(1.0, 14);
    });

    it('computes correct values for small inputs (|x| <= 8)', async () => {
      // NumPy reference values
      const inputs = [0, 1, 2, 3, 4, 5, 6, 7, 8];
      const expected = [
        1.0,
        1.2660658777520082,
        2.2795853023360673,
        4.880792585865024,
        11.301921952136328,
        27.239871823604442,
        67.23440697647798,
        168.59390851028544,
        427.56411572180454,
      ];

      const result = track(await i0(inputs));
      expectClose(result.toArray(), expected);
    });

    it('computes correct values for large inputs (|x| > 8)', async () => {
      // NumPy reference values
      const inputs = [9, 10, 15, 20];
      const expected = [
        1093.5883545113747,
        2815.7166284662544,
        339649.37329596426,
        4.355828227818201e7,
      ];

      const result = track(await i0(inputs));
      expectClose(result.toArray(), expected, 1e-8);
    });

    it('handles negative inputs using absolute value', async () => {
      const positive = track(await i0(5));
      const negative = track(await i0(-5));
      expect(positive.item()).toBeCloseTo(negative.item(), 14);
    });

    it('handles array input', async () => {
      const arr = track(await NDArray.fromArray([0, 1, 2]));
      const result = track(await i0(arr));
      expect(result.shape).toEqual([3]);
      expect(result.getFlat(0)).toBeCloseTo(1.0, 14);
      expect(result.getFlat(1)).toBeCloseTo(1.2660658777520082, 10);
      expect(result.getFlat(2)).toBeCloseTo(2.2795853023360673, 10);
    });
  });

  // =========================================================================
  // Blackman Window Tests
  // =========================================================================
  describe('blackman', () => {
    it('returns empty array for M < 1', async () => {
      const result0 = track(await blackman(0));
      const resultNeg = track(await blackman(-5));
      expect(result0.size).toBe(0);
      expect(resultNeg.size).toBe(0);
    });

    it('returns [1.0] for M = 1', async () => {
      const result = track(await blackman(1));
      expect(result.toArray()).toEqual([1.0]);
    });

    it('computes correct values for M = 12', async () => {
      // NumPy: np.blackman(12)
      const expected = [
        -1.3877787807814457e-17,
        0.032606434624560324,
        0.159903634783434,
        0.4143979812474828,
        0.7360451799107798,
        0.9670467694337431,
        0.9670467694337431,
        0.7360451799107798,
        0.4143979812474828,
        0.159903634783434,
        0.032606434624560324,
        -1.3877787807814457e-17,
      ];

      const result = track(await blackman(12));
      expectClose(result.toArray(), expected);
    });

    it('is symmetric', async () => {
      const result = track(await blackman(11));
      const arr = result.toArray();
      for (let i = 0; i < 5; i++) {
        expect(arr[i]).toBeCloseTo(arr[10 - i], 14);
      }
    });

    it('has maximum value of 1 for odd M', async () => {
      const result = track(await blackman(11));
      const arr = result.toArray();
      const maxVal = Math.max(...arr);
      expect(maxVal).toBeCloseTo(1.0, 10);
    });
  });

  // =========================================================================
  // Hanning Window Tests
  // =========================================================================
  describe('hanning', () => {
    it('returns empty array for M < 1', async () => {
      const result = track(await hanning(0));
      expect(result.size).toBe(0);
    });

    it('returns [1.0] for M = 1', async () => {
      const result = track(await hanning(1));
      expect(result.toArray()).toEqual([1.0]);
    });

    it('computes correct values for M = 12', async () => {
      // NumPy: np.hanning(12)
      const expected = [
        0.0,
        0.07937323293697949,
        0.29229249349905033,
        0.5711574182308893,
        0.8274303757572243,
        0.9797464868072289,
        0.9797464868072289,
        0.8274303757572243,
        0.5711574182308893,
        0.29229249349905033,
        0.07937323293697949,
        0.0,
      ];

      const result = track(await hanning(12));
      expectClose(result.toArray(), expected);
    });

    it('has zero endpoints', async () => {
      const result = track(await hanning(12));
      const arr = result.toArray();
      expect(arr[0]).toBeCloseTo(0, 14);
      expect(arr[11]).toBeCloseTo(0, 14);
    });

    it('is symmetric', async () => {
      const result = track(await hanning(10));
      const arr = result.toArray();
      for (let i = 0; i < 5; i++) {
        expect(arr[i]).toBeCloseTo(arr[9 - i], 14);
      }
    });
  });

  // =========================================================================
  // Hamming Window Tests
  // =========================================================================
  describe('hamming', () => {
    it('returns empty array for M < 1', async () => {
      const result = track(await hamming(0));
      expect(result.size).toBe(0);
    });

    it('returns [1.0] for M = 1', async () => {
      const result = track(await hamming(1));
      expect(result.toArray()).toEqual([1.0]);
    });

    it('computes correct values for M = 12', async () => {
      // NumPy: np.hamming(12)
      const expected = [
        0.08,
        0.15302337489435634,
        0.3489090909090909,
        0.6054648259178303,
        0.8412359409678934,
        0.9813667678238129,
        0.9813667678238129,
        0.8412359409678934,
        0.6054648259178303,
        0.3489090909090909,
        0.15302337489435634,
        0.08,
      ];

      const result = track(await hamming(12));
      expectClose(result.toArray(), expected);
    });

    it('has raised endpoints (0.08)', async () => {
      const result = track(await hamming(12));
      const arr = result.toArray();
      expect(arr[0]).toBeCloseTo(0.08, 14);
      expect(arr[11]).toBeCloseTo(0.08, 14);
    });

    it('is symmetric', async () => {
      const result = track(await hamming(10));
      const arr = result.toArray();
      for (let i = 0; i < 5; i++) {
        expect(arr[i]).toBeCloseTo(arr[9 - i], 14);
      }
    });
  });

  // =========================================================================
  // Bartlett Window Tests
  // =========================================================================
  describe('bartlett', () => {
    it('returns empty array for M < 1', async () => {
      const result = track(await bartlett(0));
      expect(result.size).toBe(0);
    });

    it('returns [1.0] for M = 1', async () => {
      const result = track(await bartlett(1));
      expect(result.toArray()).toEqual([1.0]);
    });

    it('computes correct values for M = 12', async () => {
      // NumPy: np.bartlett(12)
      const expected = [
        0.0,
        0.18181818181818182,
        0.36363636363636365,
        0.5454545454545454,
        0.7272727272727273,
        0.9090909090909091,
        0.9090909090909091,
        0.7272727272727273,
        0.5454545454545454,
        0.36363636363636365,
        0.18181818181818182,
        0.0,
      ];

      const result = track(await bartlett(12));
      expectClose(result.toArray(), expected);
    });

    it('has zero endpoints', async () => {
      const result = track(await bartlett(12));
      const arr = result.toArray();
      expect(arr[0]).toBeCloseTo(0, 14);
      expect(arr[11]).toBeCloseTo(0, 14);
    });

    it('forms a triangular shape', async () => {
      const result = track(await bartlett(11));
      const arr = result.toArray();
      // Should increase to middle, then decrease
      for (let i = 1; i <= 5; i++) {
        expect(arr[i]).toBeGreaterThan(arr[i - 1]);
      }
      for (let i = 6; i <= 10; i++) {
        expect(arr[i]).toBeLessThan(arr[i - 1]);
      }
    });

    it('is symmetric', async () => {
      const result = track(await bartlett(10));
      const arr = result.toArray();
      for (let i = 0; i < 5; i++) {
        expect(arr[i]).toBeCloseTo(arr[9 - i], 14);
      }
    });
  });

  // =========================================================================
  // Kaiser Window Tests
  // =========================================================================
  describe('kaiser', () => {
    it('returns empty array for M < 1', async () => {
      const result = track(await kaiser(0, 14));
      expect(result.size).toBe(0);
    });

    it('returns [1.0] for M = 1', async () => {
      const result = track(await kaiser(1, 14));
      expect(result.toArray()).toEqual([1.0]);
    });

    it('approximates rectangular window when beta = 0', async () => {
      // NumPy: np.kaiser(12, 0)
      const expected = [
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
      ];

      const result = track(await kaiser(12, 0));
      expectClose(result.toArray(), expected);
    });

    it('computes correct values for M = 12, beta = 14', async () => {
      // NumPy: np.kaiser(12, 14)
      const expected = [
        7.726866835270368e-6,
        0.003460091937889429,
        0.04652001885143685,
        0.22973712022098938,
        0.5998853159552917,
        0.9456748983666177,
        0.9456748983666177,
        0.5998853159552917,
        0.22973712022098938,
        0.04652001885143685,
        0.003460091937889429,
        7.726866835270368e-6,
      ];

      const result = track(await kaiser(12, 14));
      expectClose(result.toArray(), expected, 1e-8);
    });

    it('is symmetric', async () => {
      const result = track(await kaiser(11, 14));
      const arr = result.toArray();
      for (let i = 0; i < 5; i++) {
        expect(arr[i]).toBeCloseTo(arr[10 - i], 10);
      }
    });

    it('has maximum value of 1 at center for odd M', async () => {
      const result = track(await kaiser(11, 14));
      const arr = result.toArray();
      const maxVal = Math.max(...arr);
      expect(maxVal).toBeCloseTo(1.0, 10);
      expect(arr[5]).toBeCloseTo(1.0, 10); // Center index
    });

    it('produces narrower window as beta increases', async () => {
      const beta5 = track(await kaiser(12, 5));
      const beta14 = track(await kaiser(12, 14));

      // With higher beta, endpoints should be smaller (more concentrated)
      expect(beta14.getFlat(0)).toBeLessThan(beta5.getFlat(0));
      expect(beta14.getFlat(1)).toBeLessThan(beta5.getFlat(1));
    });
  });

  // =========================================================================
  // Edge Cases
  // =========================================================================
  describe('Edge cases', () => {
    it('handles fractional M by truncating', async () => {
      const result = track(await blackman(5.9));
      expect(result.size).toBe(5);
    });

    it('handles very small windows', async () => {
      const b2 = track(await blackman(2));
      const h2 = track(await hanning(2));
      const m2 = track(await hamming(2));
      const t2 = track(await bartlett(2));
      const k2 = track(await kaiser(2, 14));

      expect(b2.size).toBe(2);
      expect(h2.size).toBe(2);
      expect(m2.size).toBe(2);
      expect(t2.size).toBe(2);
      expect(k2.size).toBe(2);
    });

    it('handles large window sizes', async () => {
      const result = track(await blackman(1000));
      expect(result.size).toBe(1000);
      // Check symmetry
      expect(result.getFlat(0)).toBeCloseTo(result.getFlat(999), 10);
      expect(result.getFlat(100)).toBeCloseTo(result.getFlat(899), 10);
    });
  });
});
