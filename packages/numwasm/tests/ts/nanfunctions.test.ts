/**
 * Tests for NaN-Handling Functions (Phase 23)
 *
 * Tests for nansum, nanprod, nanmin, nanmax, nanargmin, nanargmax,
 * nanmean, nanvar, nanstd, nanmedian, nanquantile, nanpercentile, nan_to_num.
 *
 * Reference: NumPy lib/_nanfunctions_impl.py
 */

import { describe, it, expect, beforeAll, afterEach } from 'vitest';
import {
  loadWasmModule,
  NDArray,
  DType,
  nansum,
  nanprod,
  nanmin,
  nanmax,
  nanargmin,
  nanargmax,
  nanmean,
  nanvar,
  nanstd,
  nanmedian,
  nanquantile,
  nanpercentile,
  nan_to_num,
} from '../../dist/numjs.mjs';

// Track resources for cleanup
const resources: NDArray[] = [];
function track<T extends NDArray>(arr: T): T {
  resources.push(arr);
  return arr;
}

// Helper to compare arrays with tolerance
function expectClose(actual: number, expected: number, rtol = 1e-7, atol = 1e-10) {
  const diff = Math.abs(actual - expected);
  const tol = atol + rtol * Math.abs(expected);
  expect(diff).toBeLessThanOrEqual(tol);
}

// Helper to compare array contents
function expectArrayClose(
  actual: number[],
  expected: number[],
  rtol = 1e-7,
  atol = 1e-10
) {
  expect(actual.length).toBe(expected.length);
  for (let i = 0; i < actual.length; i++) {
    if (Number.isNaN(expected[i])) {
      expect(Number.isNaN(actual[i])).toBe(true);
    } else {
      const diff = Math.abs(actual[i] - expected[i]);
      const tol = atol + rtol * Math.abs(expected[i]);
      expect(diff).toBeLessThanOrEqual(tol);
    }
  }
}

describe('NaN-Handling Functions (Phase 23)', () => {
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
  // 23.1 Basic NaN Aggregations
  // =========================================================================
  describe('23.1 Basic NaN Aggregations', () => {
    describe('nansum', () => {
      it('should sum treating NaN as zero', async () => {
        const arr = track(await NDArray.fromArray([1, NaN, 3]));
        const result = await nansum(arr);
        expect(result).toBe(4);
      });

      it('should return 0 for all-NaN array', async () => {
        const arr = track(await NDArray.fromArray([NaN, NaN, NaN]));
        const result = await nansum(arr);
        expect(result).toBe(0);
      });

      it('should handle 2D array with axis=0', async () => {
        const arr = track(await NDArray.fromArray([[1, NaN], [3, 4]]));
        const result = track(await nansum(arr, 0) as NDArray);
        expectArrayClose(result.toArray() as number[], [4, 4]);
      });

      it('should handle 2D array with axis=1', async () => {
        const arr = track(await NDArray.fromArray([[1, NaN], [3, 4]]));
        const result = track(await nansum(arr, 1) as NDArray);
        expectArrayClose(result.toArray() as number[], [1, 7]);
      });

      it('should respect keepdims parameter', async () => {
        const arr = track(await NDArray.fromArray([[1, NaN], [3, 4]]));
        const result = track(await nansum(arr, 1, true) as NDArray);
        expect(result.shape).toEqual([2, 1]);
      });

      it('should handle array with no NaN values', async () => {
        const arr = track(await NDArray.fromArray([1, 2, 3, 4, 5]));
        const result = await nansum(arr);
        expect(result).toBe(15);
      });
    });

    describe('nanprod', () => {
      it('should compute product treating NaN as one', async () => {
        const arr = track(await NDArray.fromArray([1, NaN, 3]));
        const result = await nanprod(arr);
        expect(result).toBe(3);
      });

      it('should return 1 for all-NaN array', async () => {
        const arr = track(await NDArray.fromArray([NaN, NaN, NaN]));
        const result = await nanprod(arr);
        expect(result).toBe(1);
      });

      it('should handle 2D array with axis=1', async () => {
        const arr = track(await NDArray.fromArray([[1, NaN], [3, 4]]));
        const result = track(await nanprod(arr, 1) as NDArray);
        expectArrayClose(result.toArray() as number[], [1, 12]);
      });

      it('should handle array with no NaN values', async () => {
        const arr = track(await NDArray.fromArray([1, 2, 3, 4]));
        const result = await nanprod(arr);
        expect(result).toBe(24);
      });
    });

    describe('nanmin', () => {
      it('should find minimum ignoring NaN', async () => {
        const arr = track(await NDArray.fromArray([NaN, 2, 1, NaN]));
        const result = await nanmin(arr);
        expect(result).toBe(1);
      });

      it('should return NaN for all-NaN array', async () => {
        const arr = track(await NDArray.fromArray([NaN, NaN, NaN]));
        const result = await nanmin(arr);
        expect(Number.isNaN(result as number)).toBe(true);
      });

      it('should handle 2D array with axis=1', async () => {
        const arr = track(await NDArray.fromArray([[1, NaN], [NaN, 4]]));
        const result = track(await nanmin(arr, 1) as NDArray);
        expectArrayClose(result.toArray() as number[], [1, 4]);
      });

      it('should handle negative values', async () => {
        const arr = track(await NDArray.fromArray([-5, NaN, -2, NaN, 3]));
        const result = await nanmin(arr);
        expect(result).toBe(-5);
      });
    });

    describe('nanmax', () => {
      it('should find maximum ignoring NaN', async () => {
        const arr = track(await NDArray.fromArray([NaN, 2, 3, NaN]));
        const result = await nanmax(arr);
        expect(result).toBe(3);
      });

      it('should return NaN for all-NaN array', async () => {
        const arr = track(await NDArray.fromArray([NaN, NaN, NaN]));
        const result = await nanmax(arr);
        expect(Number.isNaN(result as number)).toBe(true);
      });

      it('should handle 2D array with axis=1', async () => {
        const arr = track(await NDArray.fromArray([[1, NaN], [NaN, 4]]));
        const result = track(await nanmax(arr, 1) as NDArray);
        expectArrayClose(result.toArray() as number[], [1, 4]);
      });

      it('should handle negative values', async () => {
        const arr = track(await NDArray.fromArray([-5, NaN, -2, NaN, 3]));
        const result = await nanmax(arr);
        expect(result).toBe(3);
      });
    });

    describe('nanargmin', () => {
      it('should find index of minimum ignoring NaN', async () => {
        const arr = track(await NDArray.fromArray([NaN, 2, 1, NaN]));
        const result = await nanargmin(arr);
        expect(result).toBe(2);
      });

      it('should handle 2D array with axis=1', async () => {
        const arr = track(await NDArray.fromArray([[3, NaN, 1], [NaN, 4, 2]]));
        const result = track(await nanargmin(arr, 1) as NDArray);
        expectArrayClose(result.toArray() as number[], [2, 2]);
      });
    });

    describe('nanargmax', () => {
      it('should find index of maximum ignoring NaN', async () => {
        const arr = track(await NDArray.fromArray([NaN, 2, 3, NaN]));
        const result = await nanargmax(arr);
        expect(result).toBe(2);
      });

      it('should handle 2D array with axis=1', async () => {
        const arr = track(await NDArray.fromArray([[3, NaN, 1], [NaN, 4, 2]]));
        const result = track(await nanargmax(arr, 1) as NDArray);
        expectArrayClose(result.toArray() as number[], [0, 1]);
      });
    });
  });

  // =========================================================================
  // 23.2 NaN Statistics
  // =========================================================================
  describe('23.2 NaN Statistics', () => {
    describe('nanmean', () => {
      it('should compute mean ignoring NaN', async () => {
        const arr = track(await NDArray.fromArray([1, NaN, 3]));
        const result = await nanmean(arr);
        expect(result).toBe(2);
      });

      it('should return NaN for all-NaN array', async () => {
        const arr = track(await NDArray.fromArray([NaN, NaN, NaN]));
        const result = await nanmean(arr);
        expect(Number.isNaN(result as number)).toBe(true);
      });

      it('should handle 2D array with axis=1', async () => {
        const arr = track(await NDArray.fromArray([[1, NaN], [3, 4]]));
        const result = track(await nanmean(arr, 1) as NDArray);
        expectArrayClose(result.toArray() as number[], [1, 3.5]);
      });

      it('should handle array with no NaN values', async () => {
        const arr = track(await NDArray.fromArray([1, 2, 3, 4, 5]));
        const result = await nanmean(arr);
        expect(result).toBe(3);
      });

      it('should respect keepdims parameter', async () => {
        const arr = track(await NDArray.fromArray([[1, NaN], [3, 4]]));
        const result = track(await nanmean(arr, 1, true) as NDArray);
        expect(result.shape).toEqual([2, 1]);
      });
    });

    describe('nanvar', () => {
      it('should compute variance ignoring NaN', async () => {
        const arr = track(await NDArray.fromArray([1, NaN, 3]));
        const result = await nanvar(arr);
        // Variance of [1, 3] with mean 2: ((1-2)^2 + (3-2)^2) / 2 = 1
        expectClose(result as number, 1.0);
      });

      it('should return NaN for all-NaN array', async () => {
        const arr = track(await NDArray.fromArray([NaN, NaN, NaN]));
        const result = await nanvar(arr);
        expect(Number.isNaN(result as number)).toBe(true);
      });

      it('should handle ddof parameter', async () => {
        const arr = track(await NDArray.fromArray([1, NaN, 3]));
        const result = await nanvar(arr, null, 1);
        // Sample variance with ddof=1: ((1-2)^2 + (3-2)^2) / 1 = 2
        expectClose(result as number, 2.0);
      });

      it('should handle 2D array with axis', async () => {
        const arr = track(await NDArray.fromArray([[1, NaN, 3], [4, 5, NaN]]));
        const result = track(await nanvar(arr, 1) as NDArray);
        // Row 0: variance of [1, 3] = 1
        // Row 1: variance of [4, 5] = 0.25
        const expected = [1.0, 0.25];
        expectArrayClose(result.toArray() as number[], expected);
      });
    });

    describe('nanstd', () => {
      it('should compute standard deviation ignoring NaN', async () => {
        const arr = track(await NDArray.fromArray([1, NaN, 3]));
        const result = await nanstd(arr);
        // std of [1, 3] = sqrt(1) = 1
        expectClose(result as number, 1.0);
      });

      it('should return NaN for all-NaN array', async () => {
        const arr = track(await NDArray.fromArray([NaN, NaN, NaN]));
        const result = await nanstd(arr);
        expect(Number.isNaN(result as number)).toBe(true);
      });

      it('should handle ddof parameter', async () => {
        const arr = track(await NDArray.fromArray([1, NaN, 3]));
        const result = await nanstd(arr, null, 1);
        // Sample std with ddof=1: sqrt(2) ≈ 1.414
        expectClose(result as number, Math.sqrt(2));
      });
    });

    describe('nanmedian', () => {
      it('should compute median ignoring NaN', async () => {
        const arr = track(await NDArray.fromArray([1, NaN, 3, 4]));
        const result = await nanmedian(arr);
        // Median of [1, 3, 4] = 3
        expect(result).toBe(3);
      });

      it('should return NaN for all-NaN array', async () => {
        const arr = track(await NDArray.fromArray([NaN, NaN, NaN]));
        const result = await nanmedian(arr);
        expect(Number.isNaN(result as number)).toBe(true);
      });

      it('should handle even number of non-NaN values', async () => {
        const arr = track(await NDArray.fromArray([1, NaN, 3, 5, NaN]));
        const result = await nanmedian(arr);
        // Median of [1, 3, 5] = 3
        expect(result).toBe(3);
      });

      it('should handle 2D array with axis', async () => {
        const arr = track(await NDArray.fromArray([[1, NaN, 3], [4, 5, NaN]]));
        const result = track(await nanmedian(arr, 1) as NDArray);
        // Row 0: median of [1, 3] = 2
        // Row 1: median of [4, 5] = 4.5
        const expected = [2, 4.5];
        expectArrayClose(result.toArray() as number[], expected);
      });
    });
  });

  // =========================================================================
  // 23.3 NaN Quantiles
  // =========================================================================
  describe('23.3 NaN Quantiles', () => {
    describe('nanquantile', () => {
      it('should compute quantile ignoring NaN', async () => {
        const arr = track(await NDArray.fromArray([1, 2, NaN, 4]));
        const result = await nanquantile(arr, 0.5);
        // Median of [1, 2, 4] = 2
        expect(result).toBe(2);
      });

      it('should return NaN for all-NaN array', async () => {
        const arr = track(await NDArray.fromArray([NaN, NaN, NaN]));
        const result = await nanquantile(arr, 0.5);
        expect(Number.isNaN(result as number)).toBe(true);
      });

      it('should throw for invalid quantile values', async () => {
        const arr = track(await NDArray.fromArray([1, 2, 3]));
        await expect(nanquantile(arr, 1.5)).rejects.toThrow();
        await expect(nanquantile(arr, -0.5)).rejects.toThrow();
      });

      it('should handle multiple quantiles', async () => {
        const arr = track(await NDArray.fromArray([1, 2, NaN, 4, 5]));
        const result = track(await nanquantile(arr, [0, 0.5, 1]) as NDArray);
        // Values: [1, 2, 4, 5]
        // q=0: 1, q=0.5: (2+4)/2=3, q=1: 5
        expectArrayClose(result.toArray() as number[], [1, 3, 5]);
      });

      it('should support different interpolation methods', async () => {
        const arr = track(await NDArray.fromArray([1, NaN, 3]));

        const lower = await nanquantile(arr, 0.5, null, 'lower');
        expect(lower).toBe(1);

        const higher = await nanquantile(arr, 0.5, null, 'higher');
        expect(higher).toBe(3);

        const midpoint = await nanquantile(arr, 0.5, null, 'midpoint');
        expect(midpoint).toBe(2);
      });
    });

    describe('nanpercentile', () => {
      it('should compute percentile ignoring NaN', async () => {
        const arr = track(await NDArray.fromArray([1, 2, NaN, 4]));
        const result = await nanpercentile(arr, 50);
        // Same as nanquantile(arr, 0.5)
        expect(result).toBe(2);
      });

      it('should throw for invalid percentile values', async () => {
        const arr = track(await NDArray.fromArray([1, 2, 3]));
        await expect(nanpercentile(arr, 150)).rejects.toThrow();
        await expect(nanpercentile(arr, -50)).rejects.toThrow();
      });

      it('should handle 2D array with axis', async () => {
        const arr = track(await NDArray.fromArray([[1, NaN], [3, 4]]));
        const result = track(await nanpercentile(arr, 50, 1) as NDArray);
        // Row 0: percentile of [1] = 1
        // Row 1: percentile of [3, 4] = 3.5
        expectArrayClose(result.toArray() as number[], [1, 3.5]);
      });
    });
  });

  // =========================================================================
  // 23.4 NaN Utilities
  // =========================================================================
  describe('23.4 NaN Utilities', () => {
    describe('nan_to_num', () => {
      it('should replace NaN with 0 by default', async () => {
        const arr = track(await NDArray.fromArray([1, NaN, 3]));
        const result = track(await nan_to_num(arr));
        expectArrayClose(result.toArray() as number[], [1, 0, 3]);
      });

      it('should replace Infinity with large finite numbers', async () => {
        const arr = track(await NDArray.fromArray([1, Infinity, -Infinity]));
        const result = track(await nan_to_num(arr));
        const data = result.toArray() as number[];
        expect(data[0]).toBe(1);
        expect(Number.isFinite(data[1])).toBe(true);
        expect(data[1]).toBeGreaterThan(1e300);
        expect(Number.isFinite(data[2])).toBe(true);
        expect(data[2]).toBeLessThan(-1e300);
      });

      it('should use custom replacement values', async () => {
        const arr = track(await NDArray.fromArray([NaN, Infinity, -Infinity]));
        const result = track(await nan_to_num(arr, -1, 999, -999));
        expectArrayClose(result.toArray() as number[], [-1, 999, -999]);
      });

      it('should handle array with no special values', async () => {
        const arr = track(await NDArray.fromArray([1, 2, 3, 4, 5]));
        const result = track(await nan_to_num(arr));
        expectArrayClose(result.toArray() as number[], [1, 2, 3, 4, 5]);
      });

      it('should handle mixed special values', async () => {
        const arr = track(await NDArray.fromArray([NaN, 1, Infinity, 2, -Infinity, NaN]));
        const result = track(await nan_to_num(arr, 0, 100, -100));
        expectArrayClose(result.toArray() as number[], [0, 1, 100, 2, -100, 0]);
      });
    });
  });

  // =========================================================================
  // Edge Cases
  // =========================================================================
  describe('Edge Cases', () => {
    it('should handle single element array', async () => {
      const arr = track(await NDArray.fromArray([5]));
      expect(await nansum(arr)).toBe(5);
      expect(await nanprod(arr)).toBe(5);
      expect(await nanmin(arr)).toBe(5);
      expect(await nanmax(arr)).toBe(5);
      expect(await nanmean(arr)).toBe(5);
    });

    it('should handle single NaN element', async () => {
      const arr = track(await NDArray.fromArray([NaN]));
      expect(await nansum(arr)).toBe(0);
      expect(await nanprod(arr)).toBe(1);
      expect(Number.isNaN(await nanmin(arr) as number)).toBe(true);
      expect(Number.isNaN(await nanmax(arr) as number)).toBe(true);
      expect(Number.isNaN(await nanmean(arr) as number)).toBe(true);
    });

    it('should handle mixed NaN and valid values', async () => {
      const arr = track(await NDArray.fromArray([NaN, 1, NaN, 2, NaN, 3, NaN]));
      expect(await nansum(arr)).toBe(6);
      expect(await nanprod(arr)).toBe(6);
      expect(await nanmin(arr)).toBe(1);
      expect(await nanmax(arr)).toBe(3);
      expect(await nanmean(arr)).toBe(2);
    });

    it('should handle negative axis indexing', async () => {
      const arr = track(await NDArray.fromArray([[1, NaN], [3, 4]]));
      const result = track(await nansum(arr, -1) as NDArray);
      expectArrayClose(result.toArray() as number[], [1, 7]);
    });
  });
});
