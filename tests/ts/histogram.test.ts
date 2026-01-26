/**
 * Histogram Functions Tests (Phase 24)
 *
 * Tests for histogram computation, binning, and counting functions.
 */

import { describe, it, expect, beforeAll } from 'vitest';
import {
  NDArray,
  loadWasmModule,
  bincount,
  digitize,
  histogram_bin_edges,
  histogram,
  histogram2d,
  histogramdd,
  HistogramError,
  DType,
} from '../../dist/numjs.mjs';

// Helper to check if arrays are close
function allClose(
  actual: number[],
  expected: number[],
  rtol = 1e-5,
  atol = 1e-8
): boolean {
  if (actual.length !== expected.length) return false;
  for (let i = 0; i < actual.length; i++) {
    const diff = Math.abs(actual[i] - expected[i]);
    const tolerance = atol + rtol * Math.abs(expected[i]);
    if (diff > tolerance && !Number.isNaN(expected[i])) return false;
  }
  return true;
}

describe('Histogram Functions (Phase 24)', () => {
  beforeAll(async () => {
    await loadWasmModule();
  });

  // ========================================================================
  // bincount tests
  // ========================================================================
  describe('bincount', () => {
    it('counts occurrences of integers', async () => {
      const result = await bincount([0, 1, 1, 3, 2, 1, 7]);
      const data = result.toArray() as number[];
      expect(data).toEqual([1, 3, 1, 1, 0, 0, 0, 1]);
      result.dispose();
    });

    it('handles empty array', async () => {
      const result = await bincount([]);
      expect(result.size).toBe(0);
      result.dispose();
    });

    it('handles empty array with minlength', async () => {
      const result = await bincount([], null, 5);
      const data = result.toArray() as number[];
      expect(data).toEqual([0, 0, 0, 0, 0]);
      result.dispose();
    });

    it('respects minlength parameter', async () => {
      const result = await bincount([0, 1, 2], null, 10);
      const data = result.toArray() as number[];
      expect(data.length).toBe(10);
      expect(data[0]).toBe(1);
      expect(data[1]).toBe(1);
      expect(data[2]).toBe(1);
      expect(data.slice(3)).toEqual([0, 0, 0, 0, 0, 0, 0]);
      result.dispose();
    });

    it('handles weighted counting', async () => {
      const result = await bincount([0, 1, 1, 2], [0.5, 1.0, 0.5, 2.0]);
      const data = result.toArray() as number[];
      expect(allClose(data, [0.5, 1.5, 2.0])).toBe(true);
      result.dispose();
    });

    it('handles single element', async () => {
      const result = await bincount([5]);
      const data = result.toArray() as number[];
      expect(data).toEqual([0, 0, 0, 0, 0, 1]);
      result.dispose();
    });

    it('handles all zeros', async () => {
      const result = await bincount([0, 0, 0, 0]);
      const data = result.toArray() as number[];
      expect(data).toEqual([4]);
      result.dispose();
    });

    it('throws on negative values', async () => {
      await expect(bincount([-1, 0, 1])).rejects.toThrow(HistogramError);
    });

    it('throws on non-integers', async () => {
      await expect(bincount([0, 1.5, 2])).rejects.toThrow(HistogramError);
    });

    it('throws on weights shape mismatch', async () => {
      await expect(bincount([0, 1, 2], [1.0, 2.0])).rejects.toThrow(
        HistogramError
      );
    });

    it('accepts NDArray input', async () => {
      const arr = await NDArray.fromArray([0, 1, 1, 2]);
      const result = await bincount(arr);
      const data = result.toArray() as number[];
      expect(data).toEqual([1, 2, 1]);
      arr.dispose();
      result.dispose();
    });
  });

  // ========================================================================
  // digitize tests
  // ========================================================================
  describe('digitize', () => {
    it('assigns bins correctly for increasing bins', async () => {
      const x = [0.2, 6.4, 3.0, 1.6];
      const bins = [0.0, 1.0, 2.5, 4.0, 10.0];
      const result = await digitize(x, bins);
      const data = result.toArray() as number[];
      // np.digitize([0.2, 6.4, 3.0, 1.6], [0, 1, 2.5, 4, 10]) = [1, 4, 3, 2]
      expect(data).toEqual([1, 4, 3, 2]);
      result.dispose();
    });

    it('handles values at bin edges', async () => {
      const x = [0.0, 1.0, 2.5, 4.0, 10.0];
      const bins = [0.0, 1.0, 2.5, 4.0, 10.0];
      const result = await digitize(x, bins);
      const data = result.toArray() as number[];
      // Values at edges go to the right bin (since right=false by default)
      expect(data).toEqual([1, 2, 3, 4, 5]);
      result.dispose();
    });

    it('respects right parameter', async () => {
      const x = [0.5, 1.0, 1.5];
      const bins = [0.0, 1.0, 2.0];

      const resultLeft = await digitize(x, bins, false);
      const dataLeft = resultLeft.toArray() as number[];
      // right=false: bins[i-1] <= x < bins[i]
      expect(dataLeft).toEqual([1, 2, 2]);
      resultLeft.dispose();

      const resultRight = await digitize(x, bins, true);
      const dataRight = resultRight.toArray() as number[];
      // right=true: bins[i-1] < x <= bins[i]
      expect(dataRight).toEqual([1, 1, 2]);
      resultRight.dispose();
    });

    it('handles decreasing bins', async () => {
      const x = [0.2, 6.4, 3.0, 1.6];
      const bins = [10.0, 4.0, 2.5, 1.0, 0.0];
      const result = await digitize(x, bins);
      const data = result.toArray() as number[];
      // For decreasing bins, indices are adjusted
      expect(data).toEqual([4, 1, 2, 3]);
      result.dispose();
    });

    it('handles single bin edge', async () => {
      const x = [0.5, 1.5, 2.5];
      const bins = [1.0];
      const result = await digitize(x, bins);
      const data = result.toArray() as number[];
      expect(data).toEqual([0, 1, 1]);
      result.dispose();
    });

    it('throws on non-monotonic bins', async () => {
      await expect(digitize([1, 2, 3], [0, 2, 1, 3])).rejects.toThrow(
        HistogramError
      );
    });

    it('throws on empty bins', async () => {
      await expect(digitize([1, 2, 3], [])).rejects.toThrow(HistogramError);
    });

    it('accepts NDArray inputs', async () => {
      const x = await NDArray.fromArray([0.5, 1.5, 2.5]);
      const bins = await NDArray.fromArray([0.0, 1.0, 2.0, 3.0]);
      const result = await digitize(x, bins);
      const data = result.toArray() as number[];
      expect(data).toEqual([1, 2, 3]);
      x.dispose();
      bins.dispose();
      result.dispose();
    });
  });

  // ========================================================================
  // histogram_bin_edges tests
  // ========================================================================
  describe('histogram_bin_edges', () => {
    it('computes edges for integer bin count', async () => {
      const edges = await histogram_bin_edges([1, 2, 3, 4, 5], 4);
      const data = edges.toArray() as number[];
      expect(data.length).toBe(5); // 4 bins = 5 edges
      expect(data[0]).toBe(1);
      expect(data[4]).toBe(5);
      edges.dispose();
    });

    it('passes through explicit edges', async () => {
      const edges = await histogram_bin_edges([1, 2, 3], [0, 1, 2, 3, 4]);
      const data = edges.toArray() as number[];
      expect(data).toEqual([0, 1, 2, 3, 4]);
      edges.dispose();
    });

    it('respects range parameter', async () => {
      const edges = await histogram_bin_edges([1, 2, 3, 4, 5], 2, [0, 10]);
      const data = edges.toArray() as number[];
      expect(data.length).toBe(3);
      expect(data[0]).toBe(0);
      expect(data[2]).toBe(10);
      edges.dispose();
    });

    it('handles equal min/max by expanding range', async () => {
      const edges = await histogram_bin_edges([5, 5, 5], 2);
      const data = edges.toArray() as number[];
      expect(data.length).toBe(3);
      expect(data[0]).toBe(4.5);
      expect(data[2]).toBe(5.5);
      edges.dispose();
    });

    it('handles empty array', async () => {
      const edges = await histogram_bin_edges([], 5);
      const data = edges.toArray() as number[];
      expect(data.length).toBe(6); // 5 bins + 1
      expect(data[0]).toBe(0);
      expect(data[5]).toBe(1);
      edges.dispose();
    });

    it('handles auto bin method', async () => {
      const edges = await histogram_bin_edges([1, 2, 3, 4, 5, 6, 7, 8, 9, 10], 'auto');
      const data = edges.toArray() as number[];
      expect(data.length).toBeGreaterThan(1);
      expect(data[0]).toBe(1);
      expect(data[data.length - 1]).toBe(10);
      edges.dispose();
    });

    it('handles sqrt bin method', async () => {
      const arr = Array.from({ length: 100 }, (_, i) => i);
      const edges = await histogram_bin_edges(arr, 'sqrt');
      const data = edges.toArray() as number[];
      // sqrt(100) = 10 bins = 11 edges
      expect(data.length).toBe(11);
      edges.dispose();
    });

    it('handles sturges bin method', async () => {
      const arr = Array.from({ length: 100 }, (_, i) => i);
      const edges = await histogram_bin_edges(arr, 'sturges');
      const data = edges.toArray() as number[];
      // sturges: ceil(log2(100) + 1) = ceil(7.64) = 8 bins = 9 edges
      expect(data.length).toBe(9);
      edges.dispose();
    });

    it('throws on invalid range', async () => {
      await expect(
        histogram_bin_edges([1, 2, 3], 5, [5, 0])
      ).rejects.toThrow(HistogramError);
    });

    it('throws on non-monotonic explicit edges', async () => {
      await expect(
        histogram_bin_edges([1, 2, 3], [0, 2, 1, 3])
      ).rejects.toThrow(HistogramError);
    });
  });

  // ========================================================================
  // histogram tests
  // ========================================================================
  describe('histogram', () => {
    it('computes basic histogram', async () => {
      const { hist, bin_edges } = await histogram([1, 2, 1, 3, 2, 2, 3], 3);
      const histData = hist.toArray() as number[];
      const edgesData = bin_edges.toArray() as number[];

      expect(histData.length).toBe(3);
      // Values: 1,1 in first bin, 2,2,2 in second, 3,3 in third
      // Note: bins are [1, 1.67), [1.67, 2.33), [2.33, 3]
      expect(histData[0]).toBe(2); // 1, 1
      expect(histData[1]).toBe(3); // 2, 2, 2
      expect(histData[2]).toBe(2); // 3, 3

      expect(edgesData.length).toBe(4);
      expect(edgesData[0]).toBe(1);
      expect(edgesData[3]).toBe(3);

      hist.dispose();
      bin_edges.dispose();
    });

    it('handles explicit bin edges', async () => {
      const { hist, bin_edges } = await histogram([1, 2, 3, 4], [0, 2, 4]);
      const histData = hist.toArray() as number[];

      expect(histData.length).toBe(2);
      expect(histData[0]).toBe(1); // value 1 in [0, 2)
      expect(histData[1]).toBe(3); // values 2, 3, 4 in [2, 4]

      hist.dispose();
      bin_edges.dispose();
    });

    it('includes values exactly on last edge', async () => {
      const { hist } = await histogram([1, 2, 3, 4, 5], [1, 3, 5]);
      const histData = hist.toArray() as number[];

      // 5 should be included in last bin [3, 5]
      expect(histData[1]).toBe(3); // 3, 4, 5

      hist.dispose();
    });

    it('excludes values outside range', async () => {
      const { hist } = await histogram([0, 1, 2, 3, 10], [1, 2, 3]);
      const histData = hist.toArray() as number[];

      // 0 and 10 should be excluded
      expect(histData[0]).toBe(1); // 1
      expect(histData[1]).toBe(2); // 2, 3 (3 on last edge)

      hist.dispose();
    });

    it('handles empty array', async () => {
      const { hist, bin_edges } = await histogram([], 3);
      const histData = hist.toArray() as number[];

      expect(histData).toEqual([0, 0, 0]);

      hist.dispose();
      bin_edges.dispose();
    });

    it('handles weighted histogram', async () => {
      const { hist } = await histogram(
        [0, 1, 2, 3],
        [0, 2, 4],
        null,
        false,
        [1.0, 2.0, 3.0, 4.0]
      );
      const histData = hist.toArray() as number[];

      expect(histData[0]).toBeCloseTo(3.0); // weights 1.0 + 2.0
      expect(histData[1]).toBeCloseTo(7.0); // weights 3.0 + 4.0

      hist.dispose();
    });

    it('handles density normalization', async () => {
      const { hist, bin_edges } = await histogram([1, 2, 3, 4], 2, null, true);
      const histData = hist.toArray() as number[];
      const edgesData = bin_edges.toArray() as number[];

      // Density should integrate to 1
      let integral = 0;
      for (let i = 0; i < histData.length; i++) {
        const binWidth = edgesData[i + 1] - edgesData[i];
        integral += histData[i] * binWidth;
      }
      expect(integral).toBeCloseTo(1.0);

      hist.dispose();
      bin_edges.dispose();
    });

    it('handles range parameter', async () => {
      const { hist, bin_edges } = await histogram([1, 2, 3, 4, 5], 5, [0, 10]);
      const edgesData = bin_edges.toArray() as number[];

      expect(edgesData[0]).toBe(0);
      expect(edgesData[5]).toBe(10);

      hist.dispose();
      bin_edges.dispose();
    });

    it('accepts NDArray input', async () => {
      const arr = await NDArray.fromArray([1, 2, 2, 3]);
      const { hist } = await histogram(arr, 3);
      const histData = hist.toArray() as number[];

      expect(histData[0]).toBe(1); // 1
      expect(histData[1]).toBe(2); // 2, 2
      expect(histData[2]).toBe(1); // 3

      arr.dispose();
      hist.dispose();
    });
  });

  // ========================================================================
  // histogram2d tests
  // ========================================================================
  describe('histogram2d', () => {
    it('computes 2D histogram', async () => {
      const { H, xedges, yedges } = await histogram2d(
        [0, 1, 1, 2],
        [0, 0, 1, 1],
        2
      );

      expect(H.shape).toEqual([2, 2]);
      const hData = H.toArray() as number[];

      // With 2 bins for x in [0,2] and y in [0,1]:
      // x edges: [0, 1, 2], y edges: [0, 0.5, 1]
      // Point (0,0): bin [0,0]
      // Point (1,0): bin [1,0]
      // Point (1,1): bin [1,1]
      // Point (2,1): bin [1,1] (2 is on right edge)

      H.dispose();
      xedges.dispose();
      yedges.dispose();
    });

    it('handles per-axis bin counts', async () => {
      const { H, xedges, yedges } = await histogram2d(
        [0, 1, 2, 3],
        [0, 1, 2, 3],
        [2, 4]
      );

      expect(H.shape).toEqual([2, 4]);
      expect(xedges.size).toBe(3); // 2 bins + 1
      expect(yedges.size).toBe(5); // 4 bins + 1

      H.dispose();
      xedges.dispose();
      yedges.dispose();
    });

    it('handles explicit bin edges', async () => {
      const { H, xedges, yedges } = await histogram2d(
        [0.5, 1.5, 2.5],
        [0.5, 1.5, 2.5],
        [[0, 1, 2, 3], [0, 1, 2, 3]]
      );

      expect(H.shape).toEqual([3, 3]);

      // Each point should be in diagonal bins
      const hData = H.toArray() as number[];

      H.dispose();
      xedges.dispose();
      yedges.dispose();
    });

    it('handles range parameter', async () => {
      const { H, xedges, yedges } = await histogram2d(
        [0, 1, 2],
        [0, 1, 2],
        2,
        [[0, 4], [0, 4]]
      );

      const xData = xedges.toArray() as number[];
      const yData = yedges.toArray() as number[];

      expect(xData[0]).toBe(0);
      expect(xData[2]).toBe(4);
      expect(yData[0]).toBe(0);
      expect(yData[2]).toBe(4);

      H.dispose();
      xedges.dispose();
      yedges.dispose();
    });

    it('handles weighted 2D histogram', async () => {
      const { H } = await histogram2d(
        [0, 1],
        [0, 0],
        2,
        null,
        false,
        [1.0, 2.0]
      );

      const hData = H.toArray() as number[];
      const total = hData.reduce((a, b) => a + b, 0);
      expect(total).toBeCloseTo(3.0);

      H.dispose();
    });

    it('throws on mismatched x and y lengths', async () => {
      await expect(histogram2d([0, 1, 2], [0, 1])).rejects.toThrow(
        HistogramError
      );
    });
  });

  // ========================================================================
  // histogramdd tests
  // ========================================================================
  describe('histogramdd', () => {
    it('computes N-dimensional histogram', async () => {
      const sample = [
        [0, 0],
        [1, 1],
        [1, 1],
        [2, 1],
      ];
      const { H, edges } = await histogramdd(sample, 2);

      expect(H.ndim).toBe(2);
      expect(edges.length).toBe(2);

      H.dispose();
      edges.forEach((e) => e.dispose());
    });

    it('handles 1D input', async () => {
      const { H, edges } = await histogramdd([1, 2, 3, 4, 5], 3);

      expect(H.ndim).toBe(1);
      expect(H.size).toBe(3);
      expect(edges.length).toBe(1);

      H.dispose();
      edges.forEach((e) => e.dispose());
    });

    it('handles per-dimension bin counts', async () => {
      const sample = [
        [0, 0, 0],
        [1, 1, 1],
        [2, 2, 2],
      ];
      const { H, edges } = await histogramdd(sample, [2, 3, 4]);

      expect(H.shape).toEqual([2, 3, 4]);
      expect(edges[0].size).toBe(3); // 2 bins
      expect(edges[1].size).toBe(4); // 3 bins
      expect(edges[2].size).toBe(5); // 4 bins

      H.dispose();
      edges.forEach((e) => e.dispose());
    });

    it('handles range parameter', async () => {
      const sample = [
        [0, 0],
        [1, 1],
      ];
      const { H, edges } = await histogramdd(sample, 2, [
        [0, 4],
        [0, 4],
      ]);

      const e0 = edges[0].toArray() as number[];
      const e1 = edges[1].toArray() as number[];

      expect(e0[0]).toBe(0);
      expect(e0[2]).toBe(4);
      expect(e1[0]).toBe(0);
      expect(e1[2]).toBe(4);

      H.dispose();
      edges.forEach((e) => e.dispose());
    });

    it('handles density normalization', async () => {
      const sample = [
        [0, 0],
        [1, 1],
        [2, 2],
        [3, 3],
      ];
      const { H, edges } = await histogramdd(sample, 2, null, true);

      // Density should sum to 1 when multiplied by bin volumes
      const hData = H.toArray() as number[];
      const e0 = edges[0].toArray() as number[];
      const e1 = edges[1].toArray() as number[];

      let integral = 0;
      for (let i = 0; i < 2; i++) {
        for (let j = 0; j < 2; j++) {
          const volume = (e0[i + 1] - e0[i]) * (e1[j + 1] - e1[j]);
          integral += hData[i * 2 + j] * volume;
        }
      }
      expect(integral).toBeCloseTo(1.0);

      H.dispose();
      edges.forEach((e) => e.dispose());
    });

    it('handles weighted histogram', async () => {
      const sample = [
        [0, 0],
        [1, 1],
      ];
      const { H } = await histogramdd(sample, 2, null, false, [1.0, 3.0]);

      const hData = H.toArray() as number[];
      const total = hData.reduce((a, b) => a + b, 0);
      expect(total).toBeCloseTo(4.0);

      H.dispose();
    });

    it('throws on invalid sample dimensions', async () => {
      const sample = await NDArray.fromTypedArray(
        new Float64Array(24),
        [2, 3, 4],
        DType.Float64
      );
      await expect(histogramdd(sample)).rejects.toThrow(HistogramError);
      sample.dispose();
    });

    it('throws on bins length mismatch', async () => {
      const sample = [
        [0, 0],
        [1, 1],
      ];
      await expect(histogramdd(sample, [2, 3, 4])).rejects.toThrow(
        HistogramError
      );
    });
  });

  // ========================================================================
  // HistogramError tests
  // ========================================================================
  describe('HistogramError', () => {
    it('is a proper Error subclass', () => {
      const err = new HistogramError('test error');
      expect(err).toBeInstanceOf(Error);
      expect(err).toBeInstanceOf(HistogramError);
      expect(err.name).toBe('HistogramError');
      expect(err.message).toBe('test error');
    });
  });
});
