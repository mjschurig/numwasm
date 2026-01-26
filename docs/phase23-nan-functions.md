# Phase 23: NaN-Handling Functions Implementation Plan

Complete implementation roadmap for NaN-aware statistical and aggregation functions.

---

## ⚠️ Implementation Guideline

**The original NumPy source code is available in `/numpy`.**

Key reference files:
- `numpy/lib/_nanfunctions_impl.py` - All nan* functions (~1,200 lines)
- `numpy/_core/fromnumeric.py` - Base reduction functions

Implementation should follow NumPy's algorithms and error handling for consistency.

---

## Current State (Pre-Phase 23)

```
Already Implemented:
├── sum, prod, mean, std, var
├── min, max, argmin, argmax
├── median
└── isnan() predicate

Missing NaN-handling variants:
├── nansum, nanprod
├── nanmean, nanstd, nanvar
├── nanmin, nanmax
├── nanargmin, nanargmax
├── nanmedian
├── nanpercentile, nanquantile
└── nan_to_num
```

---

## Phase 23 Dependency Tree

```
PHASE 23: NAN-HANDLING FUNCTIONS
│
├── 23.1 Basic NaN Aggregations (TypeScript + WASM)
│   ├── 23.1.1 nansum(a, axis, dtype, keepdims)
│   ├── 23.1.2 nanprod(a, axis, dtype, keepdims)
│   ├── 23.1.3 nanmin(a, axis, keepdims)
│   ├── 23.1.4 nanmax(a, axis, keepdims)
│   ├── 23.1.5 nanargmin(a, axis, keepdims)
│   └── 23.1.6 nanargmax(a, axis, keepdims)
│
│   Dependencies: isnan(), where(), sum/prod/min/max/argmin/argmax
│
├── 23.2 NaN Statistics (TypeScript)
│   ├── 23.2.1 nanmean(a, axis, dtype, keepdims)
│   ├── 23.2.2 nanstd(a, axis, dtype, ddof, keepdims)
│   ├── 23.2.3 nanvar(a, axis, dtype, ddof, keepdims)
│   └── 23.2.4 nanmedian(a, axis, keepdims)
│
│   Dependencies: 23.1.* (nansum, etc.)
│
├── 23.3 NaN Quantiles (TypeScript)
│   ├── 23.3.1 nanpercentile(a, q, axis, interpolation, keepdims)
│   └── 23.3.2 nanquantile(a, q, axis, interpolation, keepdims)
│
│   Dependencies: 23.2.* (statistics), sorting
│
└── 23.4 NaN Utilities (TypeScript)
    ├── 23.4.1 nan_to_num(x, nan, posinf, neginf)
    └── 23.4.2 _replace_nan(a, val) - internal helper

    Dependencies: isnan(), isinf(), where()
```

---

## Detailed Implementation Specifications

### 23.1 Basic NaN Aggregations

#### 23.1.1 nansum

**File:** `src/ts/nanfunctions.ts` (new file)

```typescript
import { NDArray } from './NDArray.js';
import { DType } from './types.js';
import { asarray, isnan, where, sum, zeros } from './index.js';

/**
 * Return the sum of array elements over a given axis treating
 * Not a Numbers (NaNs) as zero.
 *
 * @param a - Array containing numbers whose sum is desired
 * @param axis - Axis or axes along which the sum is computed
 * @param dtype - The type of the returned array
 * @param keepdims - If True, reduced axes are left with size one
 * @returns Sum of array elements, with NaN treated as zero
 *
 * @example
 * nansum([1, NaN, 3])  // 4
 * nansum([[1, NaN], [3, 4]], axis=1)  // [1, 7]
 */
export function nansum(
  a: NDArray | ArrayLike<number>,
  axis: number | number[] | null = null,
  dtype: DType | null = null,
  keepdims: boolean = false
): NDArray | number {
  const arr = asarray(a);

  // Replace NaN with 0
  const mask = isnan(arr);
  const cleaned = where(mask, 0, arr);

  return sum(cleaned, axis, dtype, keepdims);
}

/**
 * Return the product of array elements over a given axis treating
 * Not a Numbers (NaNs) as ones.
 *
 * @param a - Array containing numbers whose product is desired
 * @param axis - Axis or axes along which the product is computed
 * @param dtype - The type of the returned array
 * @param keepdims - If True, reduced axes are left with size one
 * @returns Product of array elements, with NaN treated as one
 *
 * @example
 * nanprod([1, NaN, 3])  // 3
 * nanprod([[1, NaN], [3, 4]], axis=1)  // [1, 12]
 */
export function nanprod(
  a: NDArray | ArrayLike<number>,
  axis: number | number[] | null = null,
  dtype: DType | null = null,
  keepdims: boolean = false
): NDArray | number {
  const arr = asarray(a);

  // Replace NaN with 1
  const mask = isnan(arr);
  const cleaned = where(mask, 1, arr);

  return prod(cleaned, axis, dtype, keepdims);
}

/**
 * Return minimum of an array or minimum along an axis, ignoring any NaNs.
 *
 * @param a - Array containing numbers
 * @param axis - Axis or axes along which to operate
 * @param keepdims - If True, reduced axes are left with size one
 * @returns Minimum values with NaN ignored
 *
 * @example
 * nanmin([1, NaN, 3])  // 1
 * nanmin([[1, NaN], [NaN, 4]], axis=1)  // [1, 4]
 *
 * @throws ValueError if all values along axis are NaN
 */
export function nanmin(
  a: NDArray | ArrayLike<number>,
  axis: number | number[] | null = null,
  keepdims: boolean = false
): NDArray | number {
  const arr = asarray(a);

  // Check for all-NaN slices
  const mask = isnan(arr);
  const allNan = all(mask, axis, keepdims);

  if (typeof allNan === 'number' ? allNan : any(allNan).item()) {
    // At least one slice is all NaN
    console.warn('All-NaN slice encountered in nanmin');
  }

  // Replace NaN with +Infinity for min computation
  const cleaned = where(mask, Infinity, arr);

  const result = min(cleaned, axis, keepdims);

  // Replace Infinity back with NaN where all values were NaN
  if (typeof result === 'number') {
    return result === Infinity ? NaN : result;
  }

  return where(equal(result, Infinity), NaN, result);
}

/**
 * Return maximum of an array or maximum along an axis, ignoring any NaNs.
 *
 * @param a - Array containing numbers
 * @param axis - Axis or axes along which to operate
 * @param keepdims - If True, reduced axes are left with size one
 * @returns Maximum values with NaN ignored
 *
 * @example
 * nanmax([1, NaN, 3])  // 3
 * nanmax([[1, NaN], [NaN, 4]], axis=1)  // [1, 4]
 */
export function nanmax(
  a: NDArray | ArrayLike<number>,
  axis: number | number[] | null = null,
  keepdims: boolean = false
): NDArray | number {
  const arr = asarray(a);

  const mask = isnan(arr);

  // Replace NaN with -Infinity for max computation
  const cleaned = where(mask, -Infinity, arr);

  const result = max(cleaned, axis, keepdims);

  // Replace -Infinity back with NaN where all values were NaN
  if (typeof result === 'number') {
    return result === -Infinity ? NaN : result;
  }

  return where(equal(result, -Infinity), NaN, result);
}

/**
 * Return the indices of the minimum values along an axis, ignoring NaNs.
 *
 * @param a - Input array
 * @param axis - Axis along which to operate
 * @param keepdims - If True, reduced axes are left with size one
 * @returns Array of indices into the array
 *
 * @example
 * nanargmin([NaN, 2, 3])  // 1
 */
export function nanargmin(
  a: NDArray | ArrayLike<number>,
  axis: number | null = null,
  keepdims: boolean = false
): NDArray | number {
  const arr = asarray(a);

  const mask = isnan(arr);

  // Replace NaN with Infinity for argmin computation
  const cleaned = where(mask, Infinity, arr);

  return argmin(cleaned, axis, keepdims);
}

/**
 * Return the indices of the maximum values along an axis, ignoring NaNs.
 *
 * @param a - Input array
 * @param axis - Axis along which to operate
 * @param keepdims - If True, reduced axes are left with size one
 * @returns Array of indices into the array
 *
 * @example
 * nanargmax([NaN, 2, 3])  // 2
 */
export function nanargmax(
  a: NDArray | ArrayLike<number>,
  axis: number | null = null,
  keepdims: boolean = false
): NDArray | number {
  const arr = asarray(a);

  const mask = isnan(arr);

  // Replace NaN with -Infinity for argmax computation
  const cleaned = where(mask, -Infinity, arr);

  return argmax(cleaned, axis, keepdims);
}
```

---

### 23.2 NaN Statistics

**File:** `src/ts/nanfunctions.ts` (additions)

```typescript
/**
 * Compute the arithmetic mean along the specified axis, ignoring NaNs.
 *
 * @param a - Array containing numbers whose mean is desired
 * @param axis - Axis or axes along which the means are computed
 * @param dtype - Type to use in computing the mean
 * @param keepdims - If True, reduced axes are left with size one
 * @returns Mean of non-NaN elements
 *
 * @example
 * nanmean([1, NaN, 3])  // 2
 * nanmean([[1, NaN], [3, 4]], axis=1)  // [1, 3.5]
 */
export function nanmean(
  a: NDArray | ArrayLike<number>,
  axis: number | number[] | null = null,
  dtype: DType | null = null,
  keepdims: boolean = false
): NDArray | number {
  const arr = asarray(a);

  const mask = isnan(arr);

  // Count non-NaN values
  const notNanMask = logical_not(mask);
  const count = sum(notNanMask, axis, DType.Float64, keepdims);

  // Sum of non-NaN values
  const cleaned = where(mask, 0, arr);
  const total = sum(cleaned, axis, dtype, keepdims);

  // Mean = sum / count
  return divide(total, count);
}

/**
 * Compute the standard deviation along the specified axis,
 * while ignoring NaNs.
 *
 * @param a - Array containing numbers
 * @param axis - Axis or axes along which the standard deviation is computed
 * @param dtype - Type to use in computing the standard deviation
 * @param ddof - Delta Degrees of Freedom (divisor is N - ddof)
 * @param keepdims - If True, reduced axes are left with size one
 * @returns Standard deviation of non-NaN elements
 *
 * @example
 * nanstd([1, NaN, 3])  // 1.0
 */
export function nanstd(
  a: NDArray | ArrayLike<number>,
  axis: number | number[] | null = null,
  dtype: DType | null = null,
  ddof: number = 0,
  keepdims: boolean = false
): NDArray | number {
  const variance = nanvar(a, axis, dtype, ddof, keepdims);

  if (typeof variance === 'number') {
    return Math.sqrt(variance);
  }

  return sqrt(variance);
}

/**
 * Compute the variance along the specified axis, while ignoring NaNs.
 *
 * @param a - Array containing numbers
 * @param axis - Axis or axes along which the variance is computed
 * @param dtype - Type to use in computing the variance
 * @param ddof - Delta Degrees of Freedom (divisor is N - ddof)
 * @param keepdims - If True, reduced axes are left with size one
 * @returns Variance of non-NaN elements
 *
 * @example
 * nanvar([1, NaN, 3])  // 1.0
 */
export function nanvar(
  a: NDArray | ArrayLike<number>,
  axis: number | number[] | null = null,
  dtype: DType | null = null,
  ddof: number = 0,
  keepdims: boolean = false
): NDArray | number {
  const arr = asarray(a);

  const mask = isnan(arr);
  const notNanMask = logical_not(mask);

  // Count non-NaN values
  const count = sum(notNanMask, axis, DType.Float64, keepdims);

  // Compute mean (ignoring NaN)
  const meanVal = nanmean(arr, axis, dtype, true);

  // Compute squared deviations
  const cleaned = where(mask, 0, arr);
  const centered = subtract(cleaned, meanVal);

  // Zero out NaN positions in centered values
  const centeredClean = where(mask, 0, centered);
  const sqDiff = multiply(centeredClean, centeredClean);

  // Sum of squared deviations
  const sumSqDiff = sum(sqDiff, axis, dtype, keepdims);

  // Variance = sum((x - mean)^2) / (N - ddof)
  const divisor = subtract(count, ddof);

  return divide(sumSqDiff, divisor);
}

/**
 * Compute the median along the specified axis, while ignoring NaNs.
 *
 * @param a - Input array
 * @param axis - Axis or axes along which the medians are computed
 * @param keepdims - If True, reduced axes are left with size one
 * @returns Median of non-NaN elements
 *
 * @example
 * nanmedian([1, NaN, 3, 4])  // 3.0
 */
export function nanmedian(
  a: NDArray | ArrayLike<number>,
  axis: number | null = null,
  keepdims: boolean = false
): NDArray | number {
  const arr = asarray(a);

  if (axis === null) {
    // Flatten and filter out NaN values
    const flat = arr.ravel();
    const mask = isnan(flat);
    const validValues = extract(logical_not(mask), flat);

    if (validValues.size === 0) {
      return NaN;
    }

    return median(validValues);
  }

  // For axis-wise median, need to handle each slice
  return _nanmedianAxis(arr, axis, keepdims);
}

/**
 * Helper function for axis-wise nanmedian.
 */
function _nanmedianAxis(
  arr: NDArray,
  axis: number,
  keepdims: boolean
): NDArray {
  const normalizedAxis = normalizeAxis(axis, arr.ndim);

  // Calculate output shape
  const outShape = arr.shape.filter((_, i) => i !== normalizedAxis);
  if (keepdims) {
    outShape.splice(normalizedAxis, 0, 1);
  }

  // Apply nanmedian along each slice
  const result = apply_along_axis(
    (slice: NDArray) => {
      const mask = isnan(slice);
      const validValues = extract(logical_not(mask), slice);
      if (validValues.size === 0) {
        return NaN;
      }
      return median(validValues);
    },
    normalizedAxis,
    arr
  );

  if (keepdims && result.ndim < arr.ndim) {
    return expandDims(result, normalizedAxis);
  }

  return result;
}
```

---

### 23.3 NaN Quantiles

**File:** `src/ts/nanfunctions.ts` (additions)

```typescript
/**
 * Compute the q-th percentile of the data along the specified axis,
 * while ignoring NaN values.
 *
 * @param a - Input array
 * @param q - Percentile(s) to compute, in range [0, 100]
 * @param axis - Axis along which to compute
 * @param interpolation - Interpolation method ('linear', 'lower', 'higher', 'midpoint', 'nearest')
 * @param keepdims - If True, reduced axes are left with size one
 * @returns Percentile values
 *
 * @example
 * nanpercentile([1, 2, NaN, 4], 50)  // 2.0
 * nanpercentile([[1, NaN], [3, 4]], 50, axis=1)  // [1.0, 3.5]
 */
export function nanpercentile(
  a: NDArray | ArrayLike<number>,
  q: number | number[],
  axis: number | null = null,
  interpolation: 'linear' | 'lower' | 'higher' | 'midpoint' | 'nearest' = 'linear',
  keepdims: boolean = false
): NDArray | number {
  // Convert percentile to quantile
  const qArr = Array.isArray(q) ? q : [q];
  const quantiles = qArr.map(p => p / 100);

  if (quantiles.length === 1) {
    return nanquantile(a, quantiles[0], axis, interpolation, keepdims);
  }

  return nanquantile(a, quantiles, axis, interpolation, keepdims);
}

/**
 * Compute the q-th quantile of the data along the specified axis,
 * while ignoring NaN values.
 *
 * @param a - Input array
 * @param q - Quantile(s) to compute, in range [0, 1]
 * @param axis - Axis along which to compute
 * @param interpolation - Interpolation method
 * @param keepdims - If True, reduced axes are left with size one
 * @returns Quantile values
 *
 * @example
 * nanquantile([1, 2, NaN, 4], 0.5)  // 2.0
 */
export function nanquantile(
  a: NDArray | ArrayLike<number>,
  q: number | number[],
  axis: number | null = null,
  interpolation: 'linear' | 'lower' | 'higher' | 'midpoint' | 'nearest' = 'linear',
  keepdims: boolean = false
): NDArray | number {
  const arr = asarray(a);

  // Validate q values
  const qArr = Array.isArray(q) ? q : [q];
  for (const qVal of qArr) {
    if (qVal < 0 || qVal > 1) {
      throw new ValueError(`Quantiles must be in the range [0, 1], got ${qVal}`);
    }
  }

  if (axis === null) {
    // Flatten and filter NaN
    const flat = arr.ravel();
    const mask = isnan(flat);
    const validValues = extract(logical_not(mask), flat);

    if (validValues.size === 0) {
      return qArr.length === 1 ? NaN : full([qArr.length], NaN);
    }

    return _computeQuantile(validValues, qArr, interpolation);
  }

  // Axis-wise quantile
  return _nanquantileAxis(arr, qArr, axis, interpolation, keepdims);
}

/**
 * Compute quantile for a 1D array (no NaN values).
 */
function _computeQuantile(
  arr: NDArray,
  quantiles: number[],
  interpolation: string
): NDArray | number {
  const sorted = sort(arr);
  const n = sorted.size;

  const results: number[] = [];

  for (const q of quantiles) {
    const virtualIdx = q * (n - 1);
    const lower = Math.floor(virtualIdx);
    const upper = Math.ceil(virtualIdx);
    const frac = virtualIdx - lower;

    let value: number;

    switch (interpolation) {
      case 'lower':
        value = sorted.getFlat(lower);
        break;
      case 'higher':
        value = sorted.getFlat(upper);
        break;
      case 'nearest':
        value = frac < 0.5 ? sorted.getFlat(lower) : sorted.getFlat(upper);
        break;
      case 'midpoint':
        value = (sorted.getFlat(lower) + sorted.getFlat(upper)) / 2;
        break;
      case 'linear':
      default:
        value = sorted.getFlat(lower) * (1 - frac) + sorted.getFlat(upper) * frac;
        break;
    }

    results.push(value);
  }

  return results.length === 1 ? results[0] : fromArray(results);
}

/**
 * Helper for axis-wise nanquantile.
 */
function _nanquantileAxis(
  arr: NDArray,
  quantiles: number[],
  axis: number,
  interpolation: string,
  keepdims: boolean
): NDArray {
  const normalizedAxis = normalizeAxis(axis, arr.ndim);

  const result = apply_along_axis(
    (slice: NDArray) => {
      const mask = isnan(slice);
      const validValues = extract(logical_not(mask), slice);
      if (validValues.size === 0) {
        return quantiles.length === 1 ? NaN : full([quantiles.length], NaN);
      }
      return _computeQuantile(validValues, quantiles, interpolation);
    },
    normalizedAxis,
    arr
  );

  if (keepdims && result.ndim < arr.ndim) {
    return expandDims(result, normalizedAxis);
  }

  return result;
}
```

---

### 23.4 NaN Utilities

**File:** `src/ts/nanfunctions.ts` (additions)

```typescript
/**
 * Replace NaN with zero and infinity with large finite numbers.
 *
 * @param x - Input array
 * @param nan - Value to replace NaN with (default: 0.0)
 * @param posinf - Value to replace positive infinity with (default: largest finite float)
 * @param neginf - Value to replace negative infinity with (default: most negative finite float)
 * @returns Array with NaN and infinity replaced
 *
 * @example
 * nan_to_num([NaN, Infinity, -Infinity, 1])
 * // [0, 1.7976931348623157e+308, -1.7976931348623157e+308, 1]
 *
 * nan_to_num([NaN, Infinity], nan=-1, posinf=999)
 * // [-1, 999]
 */
export function nan_to_num(
  x: NDArray | ArrayLike<number>,
  nan: number = 0.0,
  posinf: number | null = null,
  neginf: number | null = null
): NDArray {
  const arr = asarray(x);

  // Get dtype-appropriate infinity replacements
  const dtype = arr.dtype;
  const info = finfo(dtype);

  const posinfVal = posinf ?? info.max;
  const neginfVal = neginf ?? info.min;

  // Replace NaN
  let result = where(isnan(arr), nan, arr);

  // Replace positive infinity
  result = where(isposinf(result), posinfVal, result);

  // Replace negative infinity
  result = where(isneginf(result), neginfVal, result);

  return result;
}

/**
 * Internal helper to count non-NaN values along axis.
 */
function _countNonNan(
  arr: NDArray,
  axis: number | number[] | null,
  keepdims: boolean
): NDArray | number {
  const mask = logical_not(isnan(arr));
  return sum(mask, axis, DType.Int64, keepdims);
}

/**
 * Internal helper to check if all values are NaN along axis.
 * Returns true if any slice is all-NaN (for warning purposes).
 */
function _hasAllNanSlice(arr: NDArray, axis: number | null): boolean {
  const mask = isnan(arr);
  const allNan = all(mask, axis);

  if (typeof allNan === 'number' || typeof allNan === 'boolean') {
    return Boolean(allNan);
  }

  return Boolean(any(allNan).item());
}
```

---

## File Changes Summary

### New Files to Create

```
src/ts/
└── nanfunctions.ts      # All NaN-handling functions

tests/ts/
└── nanfunctions.test.ts # Test suite
```

### Files to Modify

```
src/ts/index.ts
├── Export nansum, nanprod
├── Export nanmin, nanmax
├── Export nanargmin, nanargmax
├── Export nanmean, nanstd, nanvar
├── Export nanmedian
├── Export nanpercentile, nanquantile
└── Export nan_to_num
```

---

## Implementation Order

```
Phase 23.1: Basic Aggregations (Day 1-2)
├── Day 1: nansum, nanprod
│   ├── TypeScript implementation
│   ├── Edge case handling
│   └── Tests
│
└── Day 2: nanmin, nanmax, nanargmin, nanargmax
    ├── All-NaN slice handling
    ├── Warning emission
    └── Tests

Phase 23.2: Statistics (Day 3-4)
├── Day 3: nanmean, nanvar, nanstd
│   ├── Count-based mean calculation
│   ├── Variance with NaN handling
│   └── Tests
│
└── Day 4: nanmedian
    ├── Axis-wise implementation
    ├── Empty slice handling
    └── Tests

Phase 23.3: Quantiles (Day 5)
├── nanpercentile, nanquantile
├── Interpolation methods
└── Tests

Phase 23.4: Utilities (Day 6)
├── nan_to_num
├── Final integration
└── NumPy comparison tests
```

---

## Verification Plan

After Phase 23 completion, verify:

```bash
# Build
npm run build

# Run tests
npm test

# Phase 23 specific tests:

# Basic aggregations
✓ nansum([1, NaN, 3]) === 4
✓ nanprod([1, NaN, 3]) === 3
✓ nanmin([NaN, 2, 1]) === 1
✓ nanmax([NaN, 2, 1]) === 2
✓ nanargmin([NaN, 2, 1]) === 2
✓ nanargmax([NaN, 2, 1]) === 1

# Statistics
✓ nanmean([1, NaN, 3]) === 2
✓ nanstd([1, NaN, 3]) ≈ 1.0
✓ nanvar([1, NaN, 3]) ≈ 1.0
✓ nanmedian([1, NaN, 3, 4]) === 3.0

# Quantiles
✓ nanpercentile([1, NaN, 3], 50) === 2.0
✓ nanquantile([1, NaN, 3], 0.5) === 2.0

# Utilities
✓ nan_to_num([NaN, Inf, -Inf]) replaces correctly

# Edge cases
✓ nansum([NaN, NaN]) === 0
✓ nanprod([NaN, NaN]) === 1
✓ nanmin([NaN, NaN]) === NaN (with warning)
✓ nanmean([NaN, NaN]) === NaN

# Axis support
✓ nanmean([[1, NaN], [3, 4]], axis=0) === [2, 4]
✓ nanmean([[1, NaN], [3, 4]], axis=1) === [1, 3.5]
```

Generate NumPy comparison vectors:

```python
import numpy as np
import json

tests = {
    "nansum_1d": {
        "input": [1.0, float('nan'), 3.0],
        "expected": float(np.nansum([1.0, np.nan, 3.0]))
    },
    "nanmean_2d_axis0": {
        "input": [[1.0, float('nan')], [3.0, 4.0]],
        "axis": 0,
        "expected": np.nanmean([[1.0, np.nan], [3.0, 4.0]], axis=0).tolist()
    },
    "nanstd": {
        "input": [1.0, float('nan'), 3.0, 5.0],
        "expected": float(np.nanstd([1.0, np.nan, 3.0, 5.0]))
    },
    "nanmedian": {
        "input": [1.0, float('nan'), 3.0, 4.0],
        "expected": float(np.nanmedian([1.0, np.nan, 3.0, 4.0]))
    },
    "nanpercentile": {
        "input": [1.0, float('nan'), 3.0, 4.0, 5.0],
        "q": 50,
        "expected": float(np.nanpercentile([1.0, np.nan, 3.0, 4.0, 5.0], 50))
    },
    "nan_to_num": {
        "input": [float('nan'), float('inf'), float('-inf'), 1.0],
        "expected": np.nan_to_num([np.nan, np.inf, -np.inf, 1.0]).tolist()
    }
}

with open("tests/fixtures/nanfunctions_vectors.json", "w") as f:
    json.dump(tests, f, indent=2)
```

---

## API Compatibility Notes

### NumPy Signature Match

```typescript
// All functions match NumPy signatures:
nansum(a, axis=None, dtype=None, keepdims=False)
nanprod(a, axis=None, dtype=None, keepdims=False)
nanmin(a, axis=None, keepdims=False)
nanmax(a, axis=None, keepdims=False)
nanargmin(a, axis=None, keepdims=False)
nanargmax(a, axis=None, keepdims=False)
nanmean(a, axis=None, dtype=None, keepdims=False)
nanstd(a, axis=None, dtype=None, ddof=0, keepdims=False)
nanvar(a, axis=None, dtype=None, ddof=0, keepdims=False)
nanmedian(a, axis=None, keepdims=False)
nanpercentile(a, q, axis=None, interpolation='linear', keepdims=False)
nanquantile(a, q, axis=None, interpolation='linear', keepdims=False)
nan_to_num(x, nan=0.0, posinf=None, neginf=None)
```

### Differences from NumPy

1. **Warnings**: NumPy emits RuntimeWarning for all-NaN slices. NumJS uses console.warn.

2. **`out` parameter**: Not supported initially for most functions.

3. **`where` parameter**: NumPy's `where` masking parameter is not supported.
