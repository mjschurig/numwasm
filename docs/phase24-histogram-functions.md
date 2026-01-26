# Phase 24: Histogram Functions Implementation Plan

Complete implementation roadmap for histogram computation, binning, and counting functions.

---

## ⚠️ Implementation Guideline

**The original NumPy source code is available in `/numpy`.**

Key reference files:
- `numpy/lib/_histograms_impl.py` - Histogram functions (~1,100 lines)
- `numpy/_core/fromnumeric.py` - bincount implementation
- `numpy/_core/src/multiarray/compiled_base.c` - C implementation

Implementation should follow NumPy's algorithms and error handling for consistency.

---

## Current State (Pre-Phase 24)

```
Already Implemented:
├── searchsorted(a, v, side, sorter)
├── sort, argsort
├── unique with return_counts
└── Basic statistics (min, max, mean)

Missing:
├── histogram(a, bins, range, density, weights)
├── histogram2d(x, y, bins, range, density, weights)
├── histogramdd(sample, bins, range, density, weights)
├── histogram_bin_edges(a, bins, range, weights)
├── bincount(x, weights, minlength)
└── digitize(x, bins, right)
```

---

## Phase 24 Dependency Tree

```
PHASE 24: HISTOGRAM FUNCTIONS
│
├── 24.1 Core Binning (TypeScript + WASM)
│   ├── 24.1.1 bincount(x, weights, minlength)
│   │   ├── Fast integer binning
│   │   ├── Weighted counting
│   │   └── WASM acceleration
│   │
│   └── 24.1.2 digitize(x, bins, right)
│       ├── Bin assignment via searchsorted
│       └── Monotonicity handling
│
│   Dependencies: searchsorted, asarray
│
├── 24.2 Histogram Bin Edge Computation (TypeScript)
│   └── 24.2.1 histogram_bin_edges(a, bins, range, weights)
│       ├── String bin strategies ('auto', 'fd', 'scott', etc.)
│       ├── Integer bin count
│       └── Explicit bin edges
│
│   Dependencies: min, max, std, unique
│
├── 24.3 1D Histogram (TypeScript + WASM)
│   └── 24.3.1 histogram(a, bins, range, density, weights)
│       ├── Bin edge computation
│       ├── Sample counting
│       └── Density normalization
│
│   Dependencies: 24.1.*, 24.2.*
│
├── 24.4 Multi-dimensional Histograms (TypeScript)
│   ├── 24.4.1 histogram2d(x, y, bins, range, density, weights)
│   │   └── 2D binning as special case of histogramdd
│   │
│   └── 24.4.2 histogramdd(sample, bins, range, density, weights)
│       ├── N-dimensional binning
│       └── Per-dimension bin specification
│
│   Dependencies: 24.3.* (1D histogram)
│
└── 24.5 WASM Acceleration (C)
    ├── 24.5.1 bincount_i32, bincount_weighted
    └── 24.5.2 histogram_count (fast path)

    Dependencies: NDArray WASM core
```

---

## Detailed Implementation Specifications

### 24.1 Core Binning

#### 24.1.1 bincount

**File:** `src/ts/histogram.ts` (new file)

```typescript
import { NDArray } from './NDArray.js';
import { DType } from './types.js';
import { asarray, zeros, max as arrayMax } from './index.js';

/**
 * Count number of occurrences of each value in array of non-negative ints.
 *
 * @param x - Input array of non-negative integers
 * @param weights - Weights, array of the same shape as x (optional)
 * @param minlength - Minimum number of bins for the output array
 * @returns Array of counts (or weighted counts)
 *
 * @example
 * bincount([0, 1, 1, 3, 2, 1, 7])
 * // [1, 3, 1, 1, 0, 0, 0, 1]
 *
 * bincount([0, 1, 1, 2], weights=[0.5, 1.0, 0.5, 2.0])
 * // [0.5, 1.5, 2.0]
 */
export function bincount(
  x: NDArray | ArrayLike<number>,
  weights: NDArray | ArrayLike<number> | null = null,
  minlength: number = 0
): NDArray {
  const arr = asarray(x);

  // Validate input
  if (arr.ndim !== 1) {
    throw new ValueError('bincount only works on 1-D arrays');
  }

  // Check for non-negative integers
  const xData = arr.toArray();
  for (let i = 0; i < xData.length; i++) {
    if (xData[i] < 0 || !Number.isInteger(xData[i])) {
      throw new ValueError(
        `'x' must contain non-negative integers, got ${xData[i]} at index ${i}`
      );
    }
  }

  // Determine output size
  const maxVal = xData.length > 0 ? Math.max(...xData) : -1;
  const size = Math.max(minlength, maxVal + 1);

  // Handle weights
  if (weights !== null) {
    const weightsArr = asarray(weights);

    if (weightsArr.size !== arr.size) {
      throw new ValueError(
        `weights and x must have the same size, got ${weightsArr.size} and ${arr.size}`
      );
    }

    const wData = weightsArr.toArray();
    const result = zeros([size], DType.Float64);

    // Weighted counting (use WASM for large arrays)
    if (arr.size > 1000) {
      return _bincountWeightedWASM(xData, wData, size);
    }

    // TypeScript fallback
    for (let i = 0; i < xData.length; i++) {
      const idx = xData[i];
      result.setFlat(idx, result.getFlat(idx) + wData[i]);
    }

    return result;
  }

  // Unweighted counting
  const result = zeros([size], DType.Int64);

  // Use WASM for large arrays
  if (arr.size > 1000) {
    return _bincountWASM(xData, size);
  }

  // TypeScript fallback
  for (let i = 0; i < xData.length; i++) {
    const idx = xData[i];
    result.setFlat(idx, result.getFlat(idx) + 1);
  }

  return result;
}

/**
 * WASM-accelerated bincount for large arrays.
 */
function _bincountWASM(x: number[], size: number): NDArray {
  const result = zeros([size], DType.Int64);

  // Call WASM function
  const module = getWasmModule();
  const xPtr = module._malloc(x.length * 4);
  const resultPtr = result.dataPtr;

  // Copy input to WASM memory
  const xHeap = new Int32Array(module.HEAP32.buffer, xPtr, x.length);
  xHeap.set(x);

  module._bincount_i32(x.length, xPtr, size, resultPtr);

  module._free(xPtr);

  return result;
}

/**
 * WASM-accelerated weighted bincount.
 */
function _bincountWeightedWASM(
  x: number[],
  weights: number[],
  size: number
): NDArray {
  const result = zeros([size], DType.Float64);

  const module = getWasmModule();
  const xPtr = module._malloc(x.length * 4);
  const wPtr = module._malloc(weights.length * 8);
  const resultPtr = result.dataPtr;

  const xHeap = new Int32Array(module.HEAP32.buffer, xPtr, x.length);
  xHeap.set(x);

  const wHeap = new Float64Array(module.HEAPF64.buffer, wPtr / 8, weights.length);
  wHeap.set(weights);

  module._bincount_weighted_f64(x.length, xPtr, wPtr, size, resultPtr);

  module._free(xPtr);
  module._free(wPtr);

  return result;
}
```

#### 24.1.2 digitize

**File:** `src/ts/histogram.ts` (additions)

```typescript
/**
 * Return the indices of the bins to which each value in input array belongs.
 *
 * @param x - Input array to be binned
 * @param bins - Array of bins. Must be 1-dimensional and monotonic.
 * @param right - If True, bins[i-1] < x <= bins[i]. If False, bins[i-1] <= x < bins[i].
 * @returns Array of bin indices
 *
 * @example
 * x = [0.2, 6.4, 3.0, 1.6]
 * bins = [0, 1, 2.5, 4, 10]
 * digitize(x, bins)  // [1, 4, 3, 2]
 *
 * digitize(x, bins, right=true)  // [1, 4, 3, 2]
 */
export function digitize(
  x: NDArray | ArrayLike<number>,
  bins: NDArray | ArrayLike<number>,
  right: boolean = false
): NDArray {
  const xArr = asarray(x);
  const binsArr = asarray(bins);

  // Validate bins
  if (binsArr.ndim !== 1) {
    throw new ValueError('bins must be 1-dimensional');
  }

  if (binsArr.size < 1) {
    throw new ValueError('bins must have at least one element');
  }

  // Check monotonicity
  const binsData = binsArr.toArray();
  const isIncreasing = _isMonotonicIncreasing(binsData);
  const isDecreasing = _isMonotonicDecreasing(binsData);

  if (!isIncreasing && !isDecreasing) {
    throw new ValueError('bins must be monotonically increasing or decreasing');
  }

  // For decreasing bins, reverse and adjust
  if (isDecreasing) {
    const reversedBins = binsArr.slice([null, null, -1]); // reverse
    const indices = searchsorted(reversedBins, xArr, right ? 'left' : 'right');

    // Adjust indices for reversed bins
    return subtract(binsArr.size, indices);
  }

  // For increasing bins, use searchsorted directly
  return searchsorted(binsArr, xArr, right ? 'left' : 'right');
}

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
```

---

### 24.2 Histogram Bin Edge Computation

**File:** `src/ts/histogram.ts` (additions)

```typescript
/**
 * Supported bin estimation methods.
 */
type BinMethod = 'auto' | 'fd' | 'doane' | 'scott' | 'stone' | 'rice' | 'sturges' | 'sqrt';

/**
 * Compute histogram bin edges.
 *
 * @param a - Input data
 * @param bins - Number of bins, bin method string, or explicit bin edges
 * @param range - Lower and upper range of bins
 * @param weights - Weights (used for weighted bin estimation)
 * @returns Array of bin edges (length = n_bins + 1)
 *
 * @example
 * histogram_bin_edges([1, 2, 3, 4, 5], bins=3)
 * // [1, 2.333, 3.667, 5]
 *
 * histogram_bin_edges([1, 2, 3, 4, 5], bins='auto')
 * // Automatically determined edges
 */
export function histogram_bin_edges(
  a: NDArray | ArrayLike<number>,
  bins: number | BinMethod | NDArray | ArrayLike<number> = 'auto',
  range: [number, number] | null = null,
  weights: NDArray | ArrayLike<number> | null = null
): NDArray {
  const arr = asarray(a).ravel();

  // Handle explicit bin edges
  if (Array.isArray(bins) || (bins instanceof NDArray)) {
    const binsArr = asarray(bins);
    if (!_isMonotonicIncreasing(binsArr.toArray())) {
      throw new ValueError('bins must be monotonically increasing');
    }
    return binsArr;
  }

  // Determine range
  let [minVal, maxVal] = range ?? [NaN, NaN];

  if (isNaN(minVal)) {
    minVal = arr.size > 0 ? Number(min(arr)) : 0;
  }
  if (isNaN(maxVal)) {
    maxVal = arr.size > 0 ? Number(max(arr)) : 1;
  }

  // Handle edge case where min == max
  if (minVal === maxVal) {
    minVal = minVal - 0.5;
    maxVal = maxVal + 0.5;
  }

  // Determine number of bins
  let nBins: number;

  if (typeof bins === 'number') {
    nBins = bins;
  } else if (typeof bins === 'string') {
    nBins = _estimateBins(arr, bins as BinMethod, range);
  } else {
    throw new ValueError(`Invalid bins specification: ${bins}`);
  }

  // Ensure at least 1 bin
  nBins = Math.max(1, Math.round(nBins));

  // Generate bin edges
  return linspace(minVal, maxVal, nBins + 1);
}

/**
 * Estimate optimal number of bins using various methods.
 */
function _estimateBins(
  arr: NDArray,
  method: BinMethod,
  range: [number, number] | null
): number {
  const n = arr.size;

  if (n === 0) {
    return 1;
  }

  const [minVal, maxVal] = range ?? [Number(min(arr)), Number(max(arr))];
  const dataRange = maxVal - minVal;

  if (dataRange === 0) {
    return 1;
  }

  switch (method) {
    case 'sqrt':
      // Square root rule
      return Math.ceil(Math.sqrt(n));

    case 'sturges':
      // Sturges' formula
      return Math.ceil(Math.log2(n) + 1);

    case 'rice':
      // Rice rule
      return Math.ceil(2 * Math.cbrt(n));

    case 'scott':
      // Scott's rule
      const stdVal = Number(std(arr));
      if (stdVal === 0) return 1;
      const scottBinWidth = 3.5 * stdVal / Math.cbrt(n);
      return Math.ceil(dataRange / scottBinWidth);

    case 'fd':
      // Freedman-Diaconis rule
      const iqr = _computeIQR(arr);
      if (iqr === 0) return Math.ceil(Math.sqrt(n)); // fallback
      const fdBinWidth = 2 * iqr / Math.cbrt(n);
      return Math.ceil(dataRange / fdBinWidth);

    case 'doane':
      // Doane's formula
      const stdDev = Number(std(arr));
      const meanVal = Number(mean(arr));
      if (stdDev === 0) return 1;

      // Compute skewness
      const skewness = _computeSkewness(arr, meanVal, stdDev);
      const sigmaSg = Math.sqrt(6 * (n - 2) / ((n + 1) * (n + 3)));
      return Math.ceil(1 + Math.log2(n) + Math.log2(1 + Math.abs(skewness) / sigmaSg));

    case 'stone':
      // Stone's rule (optimal for normal distributions)
      // Uses cross-validation approach - simplified here
      const stoneN = Math.max(1, Math.ceil(n ** (1/3) * 2));
      return stoneN;

    case 'auto':
      // Auto: choose between Sturges and FD based on data
      const fdBins = _estimateBins(arr, 'fd', range);
      const sturgesBins = _estimateBins(arr, 'sturges', range);
      return Math.max(fdBins, sturgesBins);

    default:
      throw new ValueError(`Unknown bin method: ${method}`);
  }
}

/**
 * Compute interquartile range.
 */
function _computeIQR(arr: NDArray): number {
  const sorted = sort(arr.ravel());
  const n = sorted.size;

  const q1Idx = Math.floor(n * 0.25);
  const q3Idx = Math.floor(n * 0.75);

  return sorted.getFlat(q3Idx) - sorted.getFlat(q1Idx);
}

/**
 * Compute skewness.
 */
function _computeSkewness(arr: NDArray, meanVal: number, stdVal: number): number {
  const n = arr.size;
  const data = arr.toArray();

  let sumCubed = 0;
  for (let i = 0; i < n; i++) {
    sumCubed += Math.pow((data[i] - meanVal) / stdVal, 3);
  }

  return sumCubed / n;
}
```

---

### 24.3 1D Histogram

**File:** `src/ts/histogram.ts` (additions)

```typescript
/**
 * Result of histogram computation.
 */
export interface HistogramResult {
  hist: NDArray;      // Bin counts or density
  bin_edges: NDArray; // Bin edge values
}

/**
 * Compute the histogram of a dataset.
 *
 * @param a - Input data
 * @param bins - Number of bins, method string, or explicit bin edges
 * @param range - Lower and upper range of bins
 * @param density - If True, return probability density instead of counts
 * @param weights - Weights for each data point
 * @returns Tuple of (histogram values, bin edges)
 *
 * @example
 * const { hist, bin_edges } = histogram([1, 2, 1, 3, 2, 2, 3], bins=3);
 * // hist: [2, 3, 2]
 * // bin_edges: [1, 1.667, 2.333, 3]
 *
 * histogram([1, 2, 3, 4], bins=[0, 2, 4])
 * // hist: [1, 3], bin_edges: [0, 2, 4]
 */
export function histogram(
  a: NDArray | ArrayLike<number>,
  bins: number | BinMethod | NDArray | ArrayLike<number> = 10,
  range: [number, number] | null = null,
  density: boolean = false,
  weights: NDArray | ArrayLike<number> | null = null
): HistogramResult {
  const arr = asarray(a).ravel();

  // Compute bin edges
  const bin_edges = histogram_bin_edges(arr, bins, range, weights);
  const nBins = bin_edges.size - 1;

  // Initialize histogram
  const weighted = weights !== null;
  const hist = weighted
    ? zeros([nBins], DType.Float64)
    : zeros([nBins], DType.Int64);

  if (arr.size === 0) {
    return { hist, bin_edges };
  }

  // Get bin edges as array
  const edges = bin_edges.toArray();
  const data = arr.toArray();
  const weightsData = weighted ? asarray(weights!).toArray() : null;

  // Use WASM for large arrays
  if (arr.size > 1000) {
    return _histogramWASM(data, edges, weightsData, density);
  }

  // TypeScript implementation
  for (let i = 0; i < data.length; i++) {
    const val = data[i];

    // Find bin index using binary search
    let binIdx = _findBin(val, edges);

    // Handle edge cases
    if (binIdx < 0) continue; // below range
    if (binIdx >= nBins) {
      // Check if exactly on upper edge
      if (val === edges[nBins]) {
        binIdx = nBins - 1; // Include in last bin
      } else {
        continue; // above range
      }
    }

    if (weighted) {
      hist.setFlat(binIdx, hist.getFlat(binIdx) + weightsData![i]);
    } else {
      hist.setFlat(binIdx, hist.getFlat(binIdx) + 1);
    }
  }

  // Convert to density if requested
  if (density) {
    const binWidths = diff(bin_edges);
    const totalArea = Number(sum(multiply(hist, binWidths)));

    if (totalArea > 0) {
      return {
        hist: divide(hist, multiply(binWidths, totalArea)) as NDArray,
        bin_edges
      };
    }
  }

  return { hist, bin_edges };
}

/**
 * Find bin index for a value using binary search.
 * Returns index i where edges[i] <= val < edges[i+1]
 */
function _findBin(val: number, edges: number[]): number {
  // Handle out of range
  if (val < edges[0]) return -1;
  if (val > edges[edges.length - 1]) return edges.length;

  // Binary search
  let lo = 0;
  let hi = edges.length - 1;

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
 * WASM-accelerated histogram computation.
 */
function _histogramWASM(
  data: number[],
  edges: number[],
  weights: number[] | null,
  density: boolean
): HistogramResult {
  const module = getWasmModule();
  const nBins = edges.length - 1;

  // Allocate memory
  const dataPtr = module._malloc(data.length * 8);
  const edgesPtr = module._malloc(edges.length * 8);
  const histPtr = module._malloc(nBins * 8);

  // Copy data
  const dataHeap = new Float64Array(module.HEAPF64.buffer, dataPtr / 8, data.length);
  dataHeap.set(data);

  const edgesHeap = new Float64Array(module.HEAPF64.buffer, edgesPtr / 8, edges.length);
  edgesHeap.set(edges);

  if (weights !== null) {
    const weightsPtr = module._malloc(weights.length * 8);
    const weightsHeap = new Float64Array(module.HEAPF64.buffer, weightsPtr / 8, weights.length);
    weightsHeap.set(weights);

    module._histogram_weighted_f64(data.length, dataPtr, edgesPtr, nBins, weightsPtr, histPtr);
    module._free(weightsPtr);
  } else {
    module._histogram_f64(data.length, dataPtr, edgesPtr, nBins, histPtr);
  }

  // Read result
  const histHeap = new Float64Array(module.HEAPF64.buffer, histPtr / 8, nBins);
  const hist = fromArray(Array.from(histHeap), [nBins], DType.Float64);
  const bin_edges = fromArray(edges);

  // Cleanup
  module._free(dataPtr);
  module._free(edgesPtr);
  module._free(histPtr);

  // Apply density normalization if needed
  if (density) {
    const binWidths = diff(bin_edges);
    const totalArea = Number(sum(multiply(hist, binWidths)));
    if (totalArea > 0) {
      return {
        hist: divide(hist, multiply(binWidths, totalArea)) as NDArray,
        bin_edges
      };
    }
  }

  return { hist, bin_edges };
}
```

---

### 24.4 Multi-dimensional Histograms

**File:** `src/ts/histogram.ts` (additions)

```typescript
/**
 * Result of 2D histogram computation.
 */
export interface Histogram2DResult {
  H: NDArray;        // 2D histogram (shape: [nx, ny])
  xedges: NDArray;   // Bin edges along x
  yedges: NDArray;   // Bin edges along y
}

/**
 * Compute the bi-dimensional histogram of two data samples.
 *
 * @param x - Array of x coordinates
 * @param y - Array of y coordinates
 * @param bins - Number of bins or [nx, ny] or [xedges, yedges]
 * @param range - [[xmin, xmax], [ymin, ymax]]
 * @param density - If True, return probability density
 * @param weights - Weights for each sample
 * @returns 2D histogram and bin edges
 *
 * @example
 * const { H, xedges, yedges } = histogram2d(
 *   [0, 1, 1, 2],
 *   [0, 0, 1, 1],
 *   bins=2
 * );
 * // H: [[1, 0], [1, 2]]
 */
export function histogram2d(
  x: NDArray | ArrayLike<number>,
  y: NDArray | ArrayLike<number>,
  bins: number | [number, number] | [NDArray, NDArray] = 10,
  range: [[number, number], [number, number]] | null = null,
  density: boolean = false,
  weights: NDArray | ArrayLike<number> | null = null
): Histogram2DResult {
  const xArr = asarray(x).ravel();
  const yArr = asarray(y).ravel();

  if (xArr.size !== yArr.size) {
    throw new ValueError('x and y must have the same length');
  }

  // Stack into 2D sample array
  const sample = stack([xArr, yArr], 1); // shape: [N, 2]

  // Handle bins specification
  let binsSpec: [number | NDArray, number | NDArray];
  if (typeof bins === 'number') {
    binsSpec = [bins, bins];
  } else {
    binsSpec = bins as [number | NDArray, number | NDArray];
  }

  // Call histogramdd
  const result = histogramdd(sample, binsSpec, range, density, weights);

  return {
    H: result.H,
    xedges: result.edges[0],
    yedges: result.edges[1]
  };
}

/**
 * Result of N-dimensional histogram computation.
 */
export interface HistogramDDResult {
  H: NDArray;         // N-dimensional histogram
  edges: NDArray[];   // List of bin edges for each dimension
}

/**
 * Compute the multidimensional histogram of some data.
 *
 * @param sample - Data to histogram (shape: [N, D] for N samples in D dimensions)
 * @param bins - Number of bins per dimension or list of bin edges
 * @param range - Ranges for each dimension [[min, max], ...]
 * @param density - If True, return probability density
 * @param weights - Weights for each sample
 * @returns N-dimensional histogram and bin edges
 *
 * @example
 * const sample = [[0, 0], [1, 1], [1, 1], [2, 1]];
 * const { H, edges } = histogramdd(sample, bins=2);
 */
export function histogramdd(
  sample: NDArray | ArrayLike<number[]>,
  bins: number | number[] | NDArray[] = 10,
  range: Array<[number, number]> | null = null,
  density: boolean = false,
  weights: NDArray | ArrayLike<number> | null = null
): HistogramDDResult {
  const sampleArr = asarray(sample);

  // Handle 1D input
  let data: NDArray;
  let ndim: number;

  if (sampleArr.ndim === 1) {
    data = sampleArr.reshape([sampleArr.size, 1]);
    ndim = 1;
  } else if (sampleArr.ndim === 2) {
    data = sampleArr;
    ndim = sampleArr.shape[1];
  } else {
    throw new ValueError('sample must be 1D or 2D array');
  }

  const nSamples = data.shape[0];

  // Process bins for each dimension
  const binsPerDim: (number | NDArray)[] = [];
  if (typeof bins === 'number') {
    for (let d = 0; d < ndim; d++) {
      binsPerDim.push(bins);
    }
  } else if (Array.isArray(bins)) {
    if (bins.length !== ndim) {
      throw new ValueError(`bins must have ${ndim} elements, got ${bins.length}`);
    }
    for (let d = 0; d < ndim; d++) {
      binsPerDim.push(bins[d]);
    }
  }

  // Compute bin edges for each dimension
  const edges: NDArray[] = [];
  const nBinsPerDim: number[] = [];

  for (let d = 0; d < ndim; d++) {
    const dimData = data.slice([null, d]).ravel();
    const dimRange = range ? range[d] : null;

    const dimEdges = histogram_bin_edges(dimData, binsPerDim[d], dimRange);
    edges.push(dimEdges);
    nBinsPerDim.push(dimEdges.size - 1);
  }

  // Initialize histogram
  const H = zeros(nBinsPerDim, weights ? DType.Float64 : DType.Int64);

  // Bin each sample
  const weightsData = weights ? asarray(weights).toArray() : null;

  for (let i = 0; i < nSamples; i++) {
    const indices: number[] = [];
    let valid = true;

    for (let d = 0; d < ndim; d++) {
      const val = data.get(i, d);
      const dimEdges = edges[d].toArray();
      const binIdx = _findBin(val, dimEdges);

      if (binIdx < 0 || binIdx >= nBinsPerDim[d]) {
        // Handle edge case: exactly on upper bound
        if (binIdx === nBinsPerDim[d] && val === dimEdges[dimEdges.length - 1]) {
          indices.push(nBinsPerDim[d] - 1);
        } else {
          valid = false;
          break;
        }
      } else {
        indices.push(binIdx);
      }
    }

    if (valid) {
      const currentVal = H.get(...indices);
      const increment = weightsData ? weightsData[i] : 1;
      H.set(...indices, currentVal + increment);
    }
  }

  // Apply density normalization
  if (density) {
    // Compute bin volumes
    let totalVolume = 1;
    const binVolumes = ones(nBinsPerDim);

    for (let d = 0; d < ndim; d++) {
      const dimWidths = diff(edges[d]);

      // Broadcast widths to full histogram shape
      const shape = nBinsPerDim.map((_, i) => i === d ? nBinsPerDim[d] : 1);
      const broadcastWidths = dimWidths.reshape(shape);

      binVolumes.multiply(broadcastTo(broadcastWidths, nBinsPerDim));
    }

    const totalCount = sum(H);
    if (Number(totalCount) > 0) {
      return {
        H: divide(H, multiply(binVolumes, totalCount)) as NDArray,
        edges
      };
    }
  }

  return { H, edges };
}
```

---

### 24.5 WASM Acceleration

**File:** `src/wasm/histogram.c` (new file)

```c
#ifndef NUMJS_HISTOGRAM_H
#define NUMJS_HISTOGRAM_H

#include "ndarray.h"

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#define EXPORT EMSCRIPTEN_KEEPALIVE
#else
#define EXPORT
#endif

/**
 * Count occurrences of non-negative integers.
 *
 * @param n      Number of elements in x
 * @param x      Input array of non-negative integers
 * @param size   Output array size
 * @param result Output array (must be pre-allocated and zeroed)
 */
EXPORT void bincount_i32(int32_t n, const int32_t* x, int32_t size, int64_t* result);

/**
 * Weighted count of non-negative integers.
 */
EXPORT void bincount_weighted_f64(
    int32_t n,
    const int32_t* x,
    const double* weights,
    int32_t size,
    double* result
);

/**
 * Compute 1D histogram.
 *
 * @param n       Number of data points
 * @param data    Input data array
 * @param edges   Bin edge array (length = n_bins + 1)
 * @param n_bins  Number of bins
 * @param result  Output histogram (must be pre-allocated and zeroed)
 */
EXPORT void histogram_f64(
    int32_t n,
    const double* data,
    const double* edges,
    int32_t n_bins,
    double* result
);

/**
 * Compute weighted 1D histogram.
 */
EXPORT void histogram_weighted_f64(
    int32_t n,
    const double* data,
    const double* edges,
    int32_t n_bins,
    const double* weights,
    double* result
);

#endif /* NUMJS_HISTOGRAM_H */
```

**File:** `src/wasm/histogram.c` (implementation)

```c
#include "histogram.h"

/* Binary search to find bin index */
static int32_t find_bin(double val, const double* edges, int32_t n_bins) {
    if (val < edges[0]) return -1;
    if (val > edges[n_bins]) return n_bins;

    int32_t lo = 0;
    int32_t hi = n_bins;

    while (lo < hi) {
        int32_t mid = (lo + hi) / 2;
        if (edges[mid] <= val) {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }

    return lo - 1;
}

EXPORT void bincount_i32(int32_t n, const int32_t* x, int32_t size, int64_t* result) {
    for (int32_t i = 0; i < n; i++) {
        int32_t idx = x[i];
        if (idx >= 0 && idx < size) {
            result[idx]++;
        }
    }
}

EXPORT void bincount_weighted_f64(
    int32_t n,
    const int32_t* x,
    const double* weights,
    int32_t size,
    double* result
) {
    for (int32_t i = 0; i < n; i++) {
        int32_t idx = x[i];
        if (idx >= 0 && idx < size) {
            result[idx] += weights[i];
        }
    }
}

EXPORT void histogram_f64(
    int32_t n,
    const double* data,
    const double* edges,
    int32_t n_bins,
    double* result
) {
    for (int32_t i = 0; i < n; i++) {
        int32_t bin_idx = find_bin(data[i], edges, n_bins);

        /* Handle right edge inclusion */
        if (bin_idx == n_bins && data[i] == edges[n_bins]) {
            bin_idx = n_bins - 1;
        }

        if (bin_idx >= 0 && bin_idx < n_bins) {
            result[bin_idx] += 1.0;
        }
    }
}

EXPORT void histogram_weighted_f64(
    int32_t n,
    const double* data,
    const double* edges,
    int32_t n_bins,
    const double* weights,
    double* result
) {
    for (int32_t i = 0; i < n; i++) {
        int32_t bin_idx = find_bin(data[i], edges, n_bins);

        if (bin_idx == n_bins && data[i] == edges[n_bins]) {
            bin_idx = n_bins - 1;
        }

        if (bin_idx >= 0 && bin_idx < n_bins) {
            result[bin_idx] += weights[i];
        }
    }
}
```

---

## File Changes Summary

### New Files to Create

```
src/ts/
└── histogram.ts         # All histogram functions

src/wasm/
└── histogram.c          # WASM histogram operations

tests/ts/
└── histogram.test.ts    # Test suite
```

### Files to Modify

```
src/ts/index.ts
├── Export bincount
├── Export digitize
├── Export histogram_bin_edges
├── Export histogram
├── Export histogram2d
└── Export histogramdd

src/ts/types.ts
└── Add WASM function declarations

scripts/build-wasm.sh
├── Add histogram.c to compilation
└── Add EXPORTED_FUNCTIONS
```

---

## Build Script Updates

Add to `scripts/build-wasm.sh` EXPORTED_FUNCTIONS:

```bash
"_bincount_i32",
"_bincount_weighted_f64",
"_histogram_f64",
"_histogram_weighted_f64"
```

---

## Implementation Order

```
Phase 24.1: Core Binning (Day 1-2)
├── Day 1: bincount
│   ├── TypeScript implementation
│   ├── WASM fast path
│   └── Tests
│
└── Day 2: digitize
    ├── Monotonicity handling
    ├── searchsorted integration
    └── Tests

Phase 24.2: Bin Edge Computation (Day 3)
├── histogram_bin_edges
├── All bin estimation methods
└── Tests

Phase 24.3: 1D Histogram (Day 4)
├── histogram function
├── Density normalization
├── WASM acceleration
└── Tests

Phase 24.4: Multi-dimensional (Day 5-6)
├── Day 5: histogram2d
│   └── Tests
│
└── Day 6: histogramdd
    ├── N-dimensional support
    └── Tests

Phase 24.5: Polish (Day 7)
├── Edge case handling
├── NumPy comparison tests
└── Documentation
```

---

## Verification Plan

After Phase 24 completion, verify:

```bash
# Build
npm run build

# Run tests
npm test

# Phase 24 specific tests:

# bincount
✓ bincount([0, 1, 1, 3, 2, 1, 7]) === [1, 3, 1, 1, 0, 0, 0, 1]
✓ bincount([0, 1, 1, 2], weights=[0.5, 1.0, 0.5, 2.0]) === [0.5, 1.5, 2.0]
✓ bincount([], minlength=5) === [0, 0, 0, 0, 0]

# digitize
✓ digitize([0.2, 6.4, 3.0], [0, 1, 2.5, 4, 10]) === [1, 4, 3]
✓ digitize handles decreasing bins correctly

# histogram_bin_edges
✓ histogram_bin_edges(data, 'auto') produces reasonable bins
✓ histogram_bin_edges(data, 10) produces 11 edges

# histogram
✓ histogram([1, 2, 1], bins=3).hist sums to 3
✓ histogram with density=true integrates to 1
✓ histogram with weights works correctly

# histogram2d
✓ histogram2d produces correct 2D shape
✓ Sum of histogram equals number of samples (or sum of weights)

# histogramdd
✓ histogramdd produces correct N-dimensional histogram
✓ Edge cases handled properly
```

Generate NumPy comparison vectors:

```python
import numpy as np
import json

tests = {
    "bincount_basic": {
        "input": [0, 1, 1, 3, 2, 1, 7],
        "expected": np.bincount([0, 1, 1, 3, 2, 1, 7]).tolist()
    },
    "bincount_weighted": {
        "input": [0, 1, 1, 2],
        "weights": [0.5, 1.0, 0.5, 2.0],
        "expected": np.bincount([0, 1, 1, 2], weights=[0.5, 1.0, 0.5, 2.0]).tolist()
    },
    "digitize": {
        "x": [0.2, 6.4, 3.0, 1.6],
        "bins": [0.0, 1.0, 2.5, 4.0, 10.0],
        "expected": np.digitize([0.2, 6.4, 3.0, 1.6], [0, 1, 2.5, 4, 10]).tolist()
    },
    "histogram": {
        "data": [1, 2, 1, 3, 2, 2, 3],
        "bins": 3,
        "hist": np.histogram([1, 2, 1, 3, 2, 2, 3], bins=3)[0].tolist(),
        "edges": np.histogram([1, 2, 1, 3, 2, 2, 3], bins=3)[1].tolist()
    },
    "histogram2d": {
        "x": [0, 1, 1, 2],
        "y": [0, 0, 1, 1],
        "bins": 2,
        "H": np.histogram2d([0, 1, 1, 2], [0, 0, 1, 1], bins=2)[0].tolist()
    }
}

with open("tests/fixtures/histogram_vectors.json", "w") as f:
    json.dump(tests, f, indent=2)
```

---

## API Compatibility Notes

### NumPy Signature Match

```typescript
bincount(x, weights=None, minlength=0)
digitize(x, bins, right=False)
histogram_bin_edges(a, bins='auto', range=None, weights=None)
histogram(a, bins=10, range=None, density=False, weights=None)
histogram2d(x, y, bins=10, range=None, density=False, weights=None)
histogramdd(sample, bins=10, range=None, density=False, weights=None)
```

### Return Value Differences

```typescript
// NumPy: returns tuple (hist, bin_edges)
// np.histogram([1, 2, 3], bins=3)
// Returns: (array([1, 1, 1]), array([1., 1.67, 2.33, 3.]))

// NumJS: returns object { hist, bin_edges }
// histogram([1, 2, 3], 3)
// Returns: { hist: NDArray, bin_edges: NDArray }
```
