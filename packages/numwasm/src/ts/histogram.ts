/**
 * Histogram Functions (Phase 24)
 *
 * Provides histogram computation, binning, and counting functions
 * compatible with NumPy's histogram API.
 */

import { NDArray } from "./_core/NDArray.js";
import { DType } from "./types.js";
import { searchsorted, min, max, mean, std } from "./statistics.js";
import { sort } from "./sorting.js";

/**
 * Error class for histogram-related errors.
 */
export class HistogramError extends Error {
  constructor(message: string) {
    super(message);
    this.name = "HistogramError";
  }
}

/**
 * Supported bin estimation methods.
 */
export type BinMethod =
  | "auto"
  | "fd"
  | "doane"
  | "scott"
  | "stone"
  | "rice"
  | "sturges"
  | "sqrt";

/**
 * Result of histogram computation.
 */
export interface HistogramResult {
  /** Bin counts or density values */
  hist: NDArray;
  /** Bin edge values (length = n_bins + 1) */
  bin_edges: NDArray;
}

/**
 * Result of 2D histogram computation.
 */
export interface Histogram2DResult {
  /** 2D histogram (shape: [nx, ny]) */
  H: NDArray;
  /** Bin edges along x */
  xedges: NDArray;
  /** Bin edges along y */
  yedges: NDArray;
}

/**
 * Result of N-dimensional histogram computation.
 */
export interface HistogramDDResult {
  /** N-dimensional histogram */
  H: NDArray;
  /** List of bin edges for each dimension */
  edges: NDArray[];
}

// ============================================================================
// Helper Functions
// ============================================================================

/**
 * Check if array is monotonically increasing.
 */
function _isMonotonicIncreasing(arr: number[]): boolean {
  for (let i = 1; i < arr.length; i++) {
    if (arr[i] < arr[i - 1]) {
      return false;
    }
  }
  return true;
}

/**
 * Check if array is monotonically decreasing.
 */
function _isMonotonicDecreasing(arr: number[]): boolean {
  for (let i = 1; i < arr.length; i++) {
    if (arr[i] > arr[i - 1]) {
      return false;
    }
  }
  return true;
}

/**
 * Compute interquartile range.
 */
async function _computeIQR(arr: NDArray): Promise<number> {
  const sorted = await sort(arr.flatten());
  const n = sorted.size;

  if (n === 0) {
    sorted.dispose();
    return 0;
  }

  const q1Idx = Math.floor(n * 0.25);
  const q3Idx = Math.floor(n * 0.75);

  const iqr = sorted.getFlat(q3Idx) - sorted.getFlat(q1Idx);
  sorted.dispose();
  return iqr;
}

/**
 * Compute skewness.
 */
function _computeSkewness(
  data: number[],
  meanVal: number,
  stdVal: number,
): number {
  const n = data.length;
  if (n === 0 || stdVal === 0) return 0;

  let sumCubed = 0;
  for (let i = 0; i < n; i++) {
    sumCubed += Math.pow((data[i] - meanVal) / stdVal, 3);
  }

  return sumCubed / n;
}

/**
 * Find bin index for a value using binary search.
 * Returns index i where edges[i] <= val < edges[i+1]
 * Returns -1 if val < edges[0]
 * Returns nBins if val > edges[nBins]
 */
function _findBin(val: number, edges: number[]): number {
  const nBins = edges.length - 1;

  // Handle out of range
  if (val < edges[0]) return -1;
  if (val > edges[nBins]) return nBins;

  // Binary search
  let lo = 0;
  let hi = nBins;

  while (lo < hi) {
    const mid = Math.floor((lo + hi) / 2);
    if (edges[mid] <= val) {
      lo = mid + 1;
    } else {
      hi = mid;
    }
  }

  return lo - 1;
}

/**
 * Generate evenly spaced numbers over a specified interval.
 * Simple implementation for histogram bin edges.
 */
async function _linspace(
  start: number,
  stop: number,
  num: number,
): Promise<NDArray> {
  const result = new Float64Array(num);
  if (num === 1) {
    result[0] = start;
  } else {
    const step = (stop - start) / (num - 1);
    for (let i = 0; i < num; i++) {
      result[i] = start + i * step;
    }
    // Ensure last value is exactly stop (avoid floating point errors)
    result[num - 1] = stop;
  }
  return NDArray.fromTypedArray(result, [num], DType.Float64);
}

// ============================================================================
// Core Binning Functions
// ============================================================================

/**
 * Count number of occurrences of each value in array of non-negative ints.
 *
 * @param x - Input array of non-negative integers
 * @param weights - Weights, array of the same shape as x (optional)
 * @param minlength - Minimum number of bins for the output array (default: 0)
 * @returns Array of counts (or weighted counts)
 *
 * @example
 * ```typescript
 * const counts = await bincount([0, 1, 1, 3, 2, 1, 7]);
 * // [1, 3, 1, 1, 0, 0, 0, 1]
 *
 * const weighted = await bincount([0, 1, 1, 2], [0.5, 1.0, 0.5, 2.0]);
 * // [0.5, 1.5, 2.0]
 * ```
 */
export async function bincount(
  x: NDArray | number[],
  weights: NDArray | number[] | null = null,
  minlength: number = 0,
): Promise<NDArray> {
  // Convert to NDArray if needed
  const arr = Array.isArray(x) ? await NDArray.fromArray(x) : x;
  const disposeArr = Array.isArray(x);

  try {
    // Validate input is 1D
    if (arr.ndim !== 1) {
      throw new HistogramError("bincount only works on 1-D arrays");
    }

    // Get data as array
    const xData = arr.toArray() as number[];

    // Check for non-negative integers
    let maxVal = -1;
    for (let i = 0; i < xData.length; i++) {
      const val = xData[i];
      if (val < 0) {
        throw new HistogramError(
          `'x' must contain non-negative integers, got ${val} at index ${i}`,
        );
      }
      if (!Number.isInteger(val)) {
        throw new HistogramError(
          `'x' must contain non-negative integers, got ${val} at index ${i}`,
        );
      }
      if (val > maxVal) maxVal = val;
    }

    // Determine output size
    const size = Math.max(minlength, maxVal + 1);

    // Handle weights
    if (weights !== null) {
      const weightsArr = Array.isArray(weights)
        ? await NDArray.fromArray(weights)
        : weights;
      const disposeWeights = Array.isArray(weights);

      try {
        if (weightsArr.size !== arr.size) {
          throw new HistogramError(
            `weights and x must have the same size, got ${weightsArr.size} and ${arr.size}`,
          );
        }

        const wData = weightsArr.toArray() as number[];
        const result = new Float64Array(size);

        // Weighted counting
        for (let i = 0; i < xData.length; i++) {
          const idx = xData[i];
          result[idx] += wData[i];
        }

        return NDArray.fromTypedArray(result, [size], DType.Float64);
      } finally {
        if (disposeWeights) weightsArr.dispose();
      }
    }

    // Unweighted counting
    const result = new Float64Array(size);
    for (let i = 0; i < xData.length; i++) {
      const idx = xData[i];
      result[idx] += 1;
    }

    return NDArray.fromTypedArray(result, [size], DType.Float64);
  } finally {
    if (disposeArr) arr.dispose();
  }
}

/**
 * Return the indices of the bins to which each value in input array belongs.
 *
 * @param x - Input array to be binned
 * @param bins - Array of bins. Must be 1-dimensional and monotonic.
 * @param right - If true, bins[i-1] < x <= bins[i]. If false (default), bins[i-1] <= x < bins[i].
 * @returns Array of bin indices
 *
 * @example
 * ```typescript
 * const x = [0.2, 6.4, 3.0, 1.6];
 * const bins = [0, 1, 2.5, 4, 10];
 * const indices = await digitize(x, bins);
 * // [1, 4, 3, 2]
 *
 * const indicesRight = await digitize(x, bins, true);
 * // [1, 4, 3, 2]
 * ```
 */
export async function digitize(
  x: NDArray | number[],
  bins: NDArray | number[],
  right: boolean = false,
): Promise<NDArray> {
  // Convert to NDArray if needed
  const xArr = Array.isArray(x) ? await NDArray.fromArray(x) : x;
  const binsArr = Array.isArray(bins) ? await NDArray.fromArray(bins) : bins;
  const disposeX = Array.isArray(x);
  const disposeBins = Array.isArray(bins);

  try {
    // Validate bins
    if (binsArr.ndim !== 1) {
      throw new HistogramError("bins must be 1-dimensional");
    }

    if (binsArr.size < 1) {
      throw new HistogramError("bins must have at least one element");
    }

    // Check monotonicity
    const binsData = binsArr.toArray() as number[];
    const isIncreasing = _isMonotonicIncreasing(binsData);
    const isDecreasing = _isMonotonicDecreasing(binsData);

    if (!isIncreasing && !isDecreasing) {
      throw new HistogramError(
        "bins must be monotonically increasing or decreasing",
      );
    }

    // For decreasing bins (but not increasing - single element is both), reverse and adjust
    if (isDecreasing && !isIncreasing) {
      // Reverse the bins
      const reversedBins = binsData.slice().reverse();
      const reversedBinsArr = await NDArray.fromArray(reversedBins);

      try {
        // Use searchsorted on reversed bins
        const side = right ? "left" : "right";
        const indices = (await searchsorted(
          reversedBinsArr,
          xArr,
          side,
        )) as NDArray;

        // Adjust indices for reversed bins
        const indicesData = indices.toArray() as number[];
        const result = new Float64Array(indicesData.length);
        const nBins = binsArr.size;

        for (let i = 0; i < indicesData.length; i++) {
          result[i] = nBins - indicesData[i];
        }

        indices.dispose();
        return NDArray.fromTypedArray(result, xArr.shape, DType.Float64);
      } finally {
        reversedBinsArr.dispose();
      }
    }

    // For increasing bins, use searchsorted directly
    const side = right ? "left" : "right";
    const result = await searchsorted(binsArr, xArr, side);

    if (typeof result === "number") {
      return NDArray.fromTypedArray(
        new Float64Array([result]),
        [1],
        DType.Float64,
      );
    }
    return result;
  } finally {
    if (disposeX) xArr.dispose();
    if (disposeBins) binsArr.dispose();
  }
}

// ============================================================================
// Bin Edge Computation
// ============================================================================

/**
 * Estimate optimal number of bins using various methods.
 */
async function _estimateBins(
  arr: NDArray,
  method: BinMethod,
  dataRange: number,
): Promise<number> {
  const n = arr.size;

  if (n === 0) {
    return 1;
  }

  if (dataRange === 0) {
    return 1;
  }

  switch (method) {
    case "sqrt":
      // Square root rule
      return Math.ceil(Math.sqrt(n));

    case "sturges":
      // Sturges' formula
      return Math.ceil(Math.log2(n) + 1);

    case "rice":
      // Rice rule
      return Math.ceil(2 * Math.cbrt(n));

    case "scott": {
      // Scott's rule
      const stdVal = (await std(arr)) as number;
      if (stdVal === 0) return 1;
      const scottBinWidth = (3.5 * stdVal) / Math.cbrt(n);
      return Math.max(1, Math.ceil(dataRange / scottBinWidth));
    }

    case "fd": {
      // Freedman-Diaconis rule
      const iqr = await _computeIQR(arr);
      if (iqr === 0) return Math.ceil(Math.sqrt(n)); // fallback
      const fdBinWidth = (2 * iqr) / Math.cbrt(n);
      return Math.max(1, Math.ceil(dataRange / fdBinWidth));
    }

    case "doane": {
      // Doane's formula
      const stdDev = (await std(arr)) as number;
      const meanVal = (await mean(arr)) as number;
      if (stdDev === 0) return 1;

      // Compute skewness
      const data = arr.toArray() as number[];
      const skewness = _computeSkewness(data, meanVal, stdDev);
      const sigmaSg = Math.sqrt((6 * (n - 2)) / ((n + 1) * (n + 3)));
      return Math.ceil(
        1 + Math.log2(n) + Math.log2(1 + Math.abs(skewness) / sigmaSg),
      );
    }

    case "stone": {
      // Stone's rule (simplified)
      return Math.max(1, Math.ceil(Math.pow(n, 1 / 3) * 2));
    }

    case "auto": {
      // Auto: choose between Sturges and FD based on data
      const fdBins = await _estimateBins(arr, "fd", dataRange);
      const sturgesBins = await _estimateBins(arr, "sturges", dataRange);
      return Math.max(fdBins, sturgesBins);
    }

    default:
      throw new HistogramError(`Unknown bin method: ${method}`);
  }
}

/**
 * Compute histogram bin edges.
 *
 * @param a - Input data
 * @param bins - Number of bins, bin method string, or explicit bin edges (default: 'auto')
 * @param range - Lower and upper range of bins
 * @param weights - Weights (used for some bin estimation methods)
 * @returns Array of bin edges (length = n_bins + 1)
 *
 * @example
 * ```typescript
 * const edges = await histogram_bin_edges([1, 2, 3, 4, 5], 3);
 * // [1, 2.333, 3.667, 5]
 *
 * const autoEdges = await histogram_bin_edges([1, 2, 3, 4, 5], 'auto');
 * // Automatically determined edges
 * ```
 */
export async function histogram_bin_edges(
  a: NDArray | number[],
  bins: number | BinMethod | NDArray | number[] = "auto",
  range: [number, number] | null = null,
  weights: NDArray | number[] | null = null,
): Promise<NDArray> {
  // Convert to NDArray and flatten
  const arr = Array.isArray(a) ? await NDArray.fromArray(a) : a;
  const disposeArr = Array.isArray(a);

  try {
    const flat = arr.flatten();
    const disposeFlat = true;

    try {
      // Handle explicit bin edges
      if (Array.isArray(bins) || bins instanceof NDArray) {
        const binsArr = Array.isArray(bins)
          ? await NDArray.fromArray(bins)
          : bins;
        const disposeBins = Array.isArray(bins);

        try {
          const binsData = binsArr.toArray() as number[];
          if (!_isMonotonicIncreasing(binsData)) {
            throw new HistogramError("bins must be monotonically increasing");
          }
          // Return a copy
          return NDArray.fromTypedArray(
            new Float64Array(binsData),
            [binsData.length],
            DType.Float64,
          );
        } finally {
          if (disposeBins) binsArr.dispose();
        }
      }

      // Determine range
      let minVal: number;
      let maxVal: number;

      if (range !== null) {
        [minVal, maxVal] = range;
        if (!Number.isFinite(minVal) || !Number.isFinite(maxVal)) {
          throw new HistogramError("range must contain finite values");
        }
        if (minVal > maxVal) {
          throw new HistogramError(
            `range must have first_edge <= last_edge, got ${minVal} > ${maxVal}`,
          );
        }
      } else {
        if (flat.size === 0) {
          minVal = 0;
          maxVal = 1;
        } else {
          minVal = (await min(flat)) as number;
          maxVal = (await max(flat)) as number;
        }
      }

      // Handle edge case where min == max
      if (minVal === maxVal) {
        minVal = minVal - 0.5;
        maxVal = maxVal + 0.5;
      }

      const dataRange = maxVal - minVal;

      // Determine number of bins
      let nBins: number;

      if (typeof bins === "number") {
        if (bins < 1) {
          throw new HistogramError("bins must be a positive integer");
        }
        nBins = Math.round(bins);
      } else if (typeof bins === "string") {
        // String method - weights not supported
        if (weights !== null) {
          throw new HistogramError(
            `Automated bin estimation is not supported for weighted data. Use explicit bin edges or integer bin count.`,
          );
        }
        nBins = await _estimateBins(flat, bins as BinMethod, dataRange);
      } else {
        throw new HistogramError(`Invalid bins specification: ${bins}`);
      }

      // Ensure at least 1 bin
      nBins = Math.max(1, nBins);

      // Generate bin edges
      return _linspace(minVal, maxVal, nBins + 1);
    } finally {
      if (disposeFlat) flat.dispose();
    }
  } finally {
    if (disposeArr) arr.dispose();
  }
}

// ============================================================================
// 1D Histogram
// ============================================================================

/**
 * Compute the histogram of a dataset.
 *
 * @param a - Input data
 * @param bins - Number of bins, method string, or explicit bin edges (default: 10)
 * @param range - Lower and upper range of bins
 * @param density - If true, return probability density instead of counts (default: false)
 * @param weights - Weights for each data point
 * @returns Object with hist (bin counts/density) and bin_edges
 *
 * @example
 * ```typescript
 * const { hist, bin_edges } = await histogram([1, 2, 1, 3, 2, 2, 3], 3);
 * // hist: [2, 3, 2]
 * // bin_edges: [1, 1.667, 2.333, 3]
 *
 * const { hist: h2 } = await histogram([1, 2, 3, 4], [0, 2, 4]);
 * // h2: [1, 3]
 * ```
 */
export async function histogram(
  a: NDArray | number[],
  bins: number | BinMethod | NDArray | number[] = 10,
  range: [number, number] | null = null,
  density: boolean = false,
  weights: NDArray | number[] | null = null,
): Promise<HistogramResult> {
  // Convert to NDArray and flatten
  const arr = Array.isArray(a) ? await NDArray.fromArray(a) : a;
  const disposeArr = Array.isArray(a);

  try {
    const flat = arr.flatten();

    try {
      // Compute bin edges
      const bin_edges = await histogram_bin_edges(flat, bins, range, weights);
      const nBins = bin_edges.size - 1;

      // Initialize histogram
      const histData = new Float64Array(nBins);

      if (flat.size === 0) {
        const hist = await NDArray.fromTypedArray(
          histData,
          [nBins],
          DType.Float64,
        );
        return { hist, bin_edges };
      }

      // Get bin edges as array
      const edges = bin_edges.toArray() as number[];
      const data = flat.toArray() as number[];

      // Handle weights
      let weightsData: number[] | null = null;
      let weightsArr: NDArray | null = null;
      let disposeWeights = false;

      if (weights !== null) {
        weightsArr = Array.isArray(weights)
          ? await NDArray.fromArray(weights)
          : weights;
        disposeWeights = Array.isArray(weights);

        const weightsFlat = weightsArr.flatten();
        if (weightsFlat.size !== flat.size) {
          weightsFlat.dispose();
          if (disposeWeights && weightsArr) weightsArr.dispose();
          throw new HistogramError(
            `weights must have the same shape as a, got ${weightsArr.size} vs ${arr.size}`,
          );
        }
        weightsData = weightsFlat.toArray() as number[];
        weightsFlat.dispose();
      }

      try {
        // Bin the data
        for (let i = 0; i < data.length; i++) {
          const val = data[i];

          // Find bin index
          let binIdx = _findBin(val, edges);

          // Handle edge cases
          if (binIdx < 0) continue; // below range
          if (binIdx >= nBins) {
            // Check if exactly on upper edge (include in last bin)
            if (val === edges[nBins]) {
              binIdx = nBins - 1;
            } else {
              continue; // above range
            }
          }

          if (weightsData !== null) {
            histData[binIdx] += weightsData[i];
          } else {
            histData[binIdx] += 1;
          }
        }

        // Create histogram array
        let hist = await NDArray.fromTypedArray(
          histData,
          [nBins],
          DType.Float64,
        );

        // Convert to density if requested
        if (density) {
          // Compute bin widths
          const binWidths = new Float64Array(nBins);
          for (let i = 0; i < nBins; i++) {
            binWidths[i] = edges[i + 1] - edges[i];
          }

          // Compute total count
          let totalCount = 0;
          for (let i = 0; i < nBins; i++) {
            totalCount += histData[i];
          }

          if (totalCount > 0) {
            // Normalize by bin width and total count
            const densityData = new Float64Array(nBins);
            for (let i = 0; i < nBins; i++) {
              densityData[i] = histData[i] / (binWidths[i] * totalCount);
            }
            hist.dispose();
            hist = await NDArray.fromTypedArray(
              densityData,
              [nBins],
              DType.Float64,
            );
          }
        }

        return { hist, bin_edges };
      } finally {
        if (disposeWeights && weightsArr) weightsArr.dispose();
      }
    } finally {
      flat.dispose();
    }
  } finally {
    if (disposeArr) arr.dispose();
  }
}

// ============================================================================
// Multi-dimensional Histograms
// ============================================================================

/**
 * Compute the multidimensional histogram of some data.
 *
 * @param sample - Data to histogram (shape: [N, D] for N samples in D dimensions,
 *                 or 1D array for single dimension)
 * @param bins - Number of bins per dimension or list of bin edges (default: 10)
 * @param range - Ranges for each dimension [[min, max], ...]
 * @param density - If true, return probability density (default: false)
 * @param weights - Weights for each sample
 * @returns Object with H (N-dimensional histogram) and edges (list of bin edges)
 *
 * @example
 * ```typescript
 * const sample = [[0, 0], [1, 1], [1, 1], [2, 1]];
 * const { H, edges } = await histogramdd(sample, 2);
 * // H shape: [2, 2]
 * ```
 */
export async function histogramdd(
  sample: NDArray | number[][] | number[],
  bins: number | number[] | NDArray[] = 10,
  range: Array<[number, number]> | null = null,
  density: boolean = false,
  weights: NDArray | number[] | null = null,
): Promise<HistogramDDResult> {
  // Convert sample to NDArray
  const sampleArr = Array.isArray(sample)
    ? await NDArray.fromArray(sample)
    : sample;
  const disposeSample = Array.isArray(sample);

  try {
    // Determine data shape and dimensions
    let data: NDArray;
    let ndim: number;
    let nSamples: number;

    if (sampleArr.ndim === 1) {
      // 1D input - treat as single dimension
      data = sampleArr.reshape([sampleArr.size, 1]);
      ndim = 1;
      nSamples = sampleArr.size;
    } else if (sampleArr.ndim === 2) {
      data = sampleArr;
      ndim = sampleArr.shape[1];
      nSamples = sampleArr.shape[0];
    } else {
      throw new HistogramError("sample must be 1D or 2D array");
    }

    // Process bins specification for each dimension
    const binsPerDim: (number | NDArray)[] = [];
    if (typeof bins === "number") {
      for (let d = 0; d < ndim; d++) {
        binsPerDim.push(bins);
      }
    } else if (Array.isArray(bins)) {
      if (bins.length !== ndim) {
        throw new HistogramError(
          `bins must have ${ndim} elements, got ${bins.length}`,
        );
      }
      for (let d = 0; d < ndim; d++) {
        binsPerDim.push(bins[d]);
      }
    } else {
      throw new HistogramError("Invalid bins specification");
    }

    // Compute bin edges for each dimension
    const edges: NDArray[] = [];
    const nBinsPerDim: number[] = [];
    const edgesArrays: number[][] = [];

    for (let d = 0; d < ndim; d++) {
      // Extract data for this dimension
      const dimDataArr = new Float64Array(nSamples);
      for (let i = 0; i < nSamples; i++) {
        dimDataArr[i] = data.get(i, d);
      }
      const dimData = await NDArray.fromTypedArray(
        dimDataArr,
        [nSamples],
        DType.Float64,
      );

      const dimRange = range ? range[d] : null;
      const dimBins = binsPerDim[d];

      const dimEdges = await histogram_bin_edges(dimData, dimBins, dimRange);
      edges.push(dimEdges);
      nBinsPerDim.push(dimEdges.size - 1);
      edgesArrays.push(dimEdges.toArray() as number[]);

      dimData.dispose();
    }

    // Compute total number of bins
    let totalBins = 1;
    for (let d = 0; d < ndim; d++) {
      totalBins *= nBinsPerDim[d];
    }

    // Initialize histogram (flat array, will reshape later)
    const histData = new Float64Array(totalBins);

    // Handle weights
    let weightsData: number[] | null = null;
    let weightsArr: NDArray | null = null;
    let disposeWeights = false;

    if (weights !== null) {
      weightsArr = Array.isArray(weights)
        ? await NDArray.fromArray(weights)
        : weights;
      disposeWeights = Array.isArray(weights);

      if (weightsArr.size !== nSamples) {
        if (disposeWeights) weightsArr.dispose();
        throw new HistogramError(
          `weights must have the same length as samples, got ${weightsArr.size} vs ${nSamples}`,
        );
      }
      weightsData = weightsArr.toArray() as number[];
    }

    try {
      // Bin each sample
      for (let i = 0; i < nSamples; i++) {
        const indices: number[] = [];
        let valid = true;

        for (let d = 0; d < ndim; d++) {
          const val = data.get(i, d);
          const dimEdges = edgesArrays[d];
          let binIdx = _findBin(val, dimEdges);

          if (binIdx < 0 || binIdx >= nBinsPerDim[d]) {
            // Handle edge case: exactly on upper bound
            if (
              binIdx === nBinsPerDim[d] &&
              val === dimEdges[dimEdges.length - 1]
            ) {
              binIdx = nBinsPerDim[d] - 1;
            } else {
              valid = false;
              break;
            }
          }
          indices.push(binIdx);
        }

        if (valid) {
          // Convert multi-dimensional index to flat index
          let flatIdx = 0;
          let multiplier = 1;
          for (let d = ndim - 1; d >= 0; d--) {
            flatIdx += indices[d] * multiplier;
            multiplier *= nBinsPerDim[d];
          }

          const increment = weightsData ? weightsData[i] : 1;
          histData[flatIdx] += increment;
        }
      }

      // Create histogram array with correct shape
      let H = await NDArray.fromTypedArray(
        histData,
        nBinsPerDim,
        DType.Float64,
      );

      // Apply density normalization
      if (density) {
        // Compute bin volumes
        let totalCount = 0;
        for (let i = 0; i < histData.length; i++) {
          totalCount += histData[i];
        }

        if (totalCount > 0) {
          // Compute volume for each bin
          const densityData = new Float64Array(totalBins);

          for (let flatIdx = 0; flatIdx < totalBins; flatIdx++) {
            // Convert flat index to multi-dimensional indices
            const indices: number[] = [];
            let remainder = flatIdx;
            for (let d = ndim - 1; d >= 0; d--) {
              indices.unshift(remainder % nBinsPerDim[d]);
              remainder = Math.floor(remainder / nBinsPerDim[d]);
            }

            // Compute bin volume
            let volume = 1;
            for (let d = 0; d < ndim; d++) {
              const dimEdges = edgesArrays[d];
              const binIdx = indices[d];
              volume *= dimEdges[binIdx + 1] - dimEdges[binIdx];
            }

            densityData[flatIdx] = histData[flatIdx] / (volume * totalCount);
          }

          H.dispose();
          H = await NDArray.fromTypedArray(
            densityData,
            nBinsPerDim,
            DType.Float64,
          );
        }
      }

      return { H, edges };
    } finally {
      if (disposeWeights && weightsArr) weightsArr.dispose();
    }
  } finally {
    if (disposeSample) sampleArr.dispose();
  }
}

/**
 * Compute the bi-dimensional histogram of two data samples.
 *
 * @param x - Array of x coordinates
 * @param y - Array of y coordinates
 * @param bins - Number of bins or [nx, ny] or [xedges, yedges] (default: 10)
 * @param range - [[xmin, xmax], [ymin, ymax]]
 * @param density - If true, return probability density (default: false)
 * @param weights - Weights for each sample
 * @returns Object with H (2D histogram), xedges, and yedges
 *
 * @example
 * ```typescript
 * const { H, xedges, yedges } = await histogram2d(
 *   [0, 1, 1, 2],
 *   [0, 0, 1, 1],
 *   2
 * );
 * // H: [[1, 0], [1, 2]]
 * ```
 */
export async function histogram2d(
  x: NDArray | number[],
  y: NDArray | number[],
  bins:
    | number
    | [number, number]
    | [NDArray | number[], NDArray | number[]] = 10,
  range: [[number, number], [number, number]] | null = null,
  density: boolean = false,
  weights: NDArray | number[] | null = null,
): Promise<Histogram2DResult> {
  // Convert to NDArrays if needed
  const xArr = Array.isArray(x) ? await NDArray.fromArray(x) : x;
  const yArr = Array.isArray(y) ? await NDArray.fromArray(y) : y;
  const disposeX = Array.isArray(x);
  const disposeY = Array.isArray(y);

  try {
    // Flatten inputs
    const xFlat = xArr.flatten();
    const yFlat = yArr.flatten();

    try {
      if (xFlat.size !== yFlat.size) {
        throw new HistogramError("x and y must have the same length");
      }

      const nSamples = xFlat.size;

      // Stack into 2D sample array [N, 2]
      const sampleData = new Float64Array(nSamples * 2);
      const xData = xFlat.toArray() as number[];
      const yData = yFlat.toArray() as number[];

      for (let i = 0; i < nSamples; i++) {
        sampleData[i * 2] = xData[i];
        sampleData[i * 2 + 1] = yData[i];
      }

      const sample = await NDArray.fromTypedArray(
        sampleData,
        [nSamples, 2],
        DType.Float64,
      );

      try {
        // Handle bins specification
        let binsSpec: number | number[] | NDArray[];

        if (typeof bins === "number") {
          binsSpec = [bins, bins];
        } else if (Array.isArray(bins)) {
          if (bins.length !== 2) {
            throw new HistogramError(
              "bins must have 2 elements for 2D histogram",
            );
          }
          // Check if bins are numbers or arrays
          const [b0, b1] = bins;
          if (typeof b0 === "number" && typeof b1 === "number") {
            binsSpec = [b0, b1];
          } else {
            // Convert any number[] to NDArray
            const processedBins: NDArray[] = [];
            for (const b of bins) {
              if (Array.isArray(b)) {
                processedBins.push(await NDArray.fromArray(b));
              } else if (b instanceof NDArray) {
                processedBins.push(b);
              } else {
                throw new HistogramError("Invalid bins specification");
              }
            }
            binsSpec = processedBins;
          }
        } else {
          throw new HistogramError("Invalid bins specification");
        }

        // Call histogramdd
        const result = await histogramdd(
          sample,
          binsSpec,
          range,
          density,
          weights,
        );

        return {
          H: result.H,
          xedges: result.edges[0],
          yedges: result.edges[1],
        };
      } finally {
        sample.dispose();
      }
    } finally {
      xFlat.dispose();
      yFlat.dispose();
    }
  } finally {
    if (disposeX) xArr.dispose();
    if (disposeY) yArr.dispose();
  }
}
