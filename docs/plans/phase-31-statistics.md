# Phase 31: Statistics Functions

## Overview

Implement remaining statistics functions that are commonly used in NumPy but not yet available. This includes weighted averages, percentiles/quantiles, numerical calculus functions, and signal correlation.

## Functions to Implement

### 1. `average(a, axis?, weights?, returned?)` - Weighted average
```typescript
async function average(
  a: NDArray,
  axis?: number | null,
  weights?: NDArray | number[] | null,
  returned?: boolean
): Promise<NDArray | [NDArray, NDArray]>
```
- Computes weighted average along axis
- If weights=null, equivalent to mean()
- If returned=true, also returns sum of weights
- Different from `ma.average()` which handles masked arrays

**NumPy Reference:** `numpy/lib/_function_base_impl.py` line ~550

### 2. `ptp(a, axis?, keepdims?)` - Peak to peak (range)
```typescript
async function ptp(
  a: NDArray,
  axis?: number | null,
  keepdims?: boolean
): Promise<NDArray | number>
```
- Returns max - min along axis
- "Peak to peak" value
- Simple but commonly used

**NumPy Reference:** `numpy/_core/fromnumeric.py` line ~2700

### 3. `percentile(a, q, axis?, ...)` - Percentile
```typescript
async function percentile(
  a: NDArray,
  q: number | number[],
  axis?: number | null,
  interpolation?: 'linear' | 'lower' | 'higher' | 'midpoint' | 'nearest',
  keepdims?: boolean
): Promise<NDArray | number>
```
- Computes q-th percentile(s) of data
- q should be in range [0, 100]
- `nanpercentile` already exists, this is the non-nan version

**NumPy Reference:** `numpy/lib/_function_base_impl.py` line ~4100

### 4. `quantile(a, q, axis?, ...)` - Quantile
```typescript
async function quantile(
  a: NDArray,
  q: number | number[],
  axis?: number | null,
  interpolation?: 'linear' | 'lower' | 'higher' | 'midpoint' | 'nearest',
  keepdims?: boolean
): Promise<NDArray | number>
```
- Computes q-th quantile(s) of data
- q should be in range [0, 1]
- Similar to percentile but different scale
- `nanquantile` already exists

**NumPy Reference:** `numpy/lib/_function_base_impl.py` line ~4300

### 5. `gradient(f, *varargs, axis?, edge_order?)` - Numerical gradient
```typescript
async function gradient(
  f: NDArray,
  ...varargs: (number | NDArray)[]
): Promise<NDArray | NDArray[]>
// With options object:
async function gradient(
  f: NDArray,
  options?: { spacing?: number | number[], axis?: number | number[], edge_order?: 1 | 2 }
): Promise<NDArray | NDArray[]>
```
- Computes N-dimensional gradient using finite differences
- Returns gradient for each axis (or single array if axis specified)
- Supports non-uniform spacing
- edge_order: 1 for first-order, 2 for second-order accurate at boundaries

**NumPy Reference:** `numpy/lib/_function_base_impl.py` line ~1000

### 6. `trapezoid(y, x?, dx?, axis?)` - Trapezoidal integration
```typescript
async function trapezoid(
  y: NDArray,
  x?: NDArray | null,
  dx?: number,
  axis?: number
): Promise<NDArray | number>
```
- Integrates along axis using trapezoidal rule
- Named `trapezoid` in NumPy 2.0 (was `trapz`)
- Essential for numerical integration

**NumPy Reference:** `numpy/lib/_function_base_impl.py` line ~5000

### 7. `correlate(a, v, mode?)` - Cross-correlation
```typescript
async function correlate(
  a: NDArray,
  v: NDArray,
  mode?: 'valid' | 'same' | 'full'
): Promise<NDArray>
```
- Computes cross-correlation of two 1D arrays
- mode='valid': Output only where arrays fully overlap
- mode='same': Output same size as first input
- mode='full': Full cross-correlation
- Foundation for signal processing

**NumPy Reference:** `numpy/_core/numeric.py` line ~700

## Implementation Strategy

### File Location
- Add to `src/ts/statistics.ts`
- Export from `src/ts/index.ts`

### Dependencies
- Existing reduction functions (sum, mean, min, max)
- Existing nanpercentile, nanquantile for reference
- No WASM changes required initially (pure TypeScript)

## Implementation Details

### average() Implementation

```typescript
/**
 * Compute weighted average along specified axis.
 *
 * @param a - Input array
 * @param axis - Axis along which to average (null for flattened)
 * @param weights - Weights for each element (same shape as a, or 1D for axis)
 * @param returned - If true, return tuple of (average, sum_of_weights)
 */
export async function average(
  a: NDArray,
  axis: number | null = null,
  weights: NDArray | number[] | null = null,
  returned: boolean = false
): Promise<NDArray | number | [NDArray | number, NDArray | number]> {
  if (weights === null) {
    // No weights: equivalent to mean
    const avg = await mean(a, axis);
    if (returned) {
      const wsum = axis === null ? a.size : a.shape[axis];
      return [avg, wsum];
    }
    return avg;
  }

  const w = Array.isArray(weights) ? fromArray(weights) : weights;

  // Validate weights shape
  if (axis === null) {
    // Weights must match total size
    if (w.size !== a.size) {
      throw new Error('weights must have same size as a when axis=null');
    }
    const flat = await ravel(a);
    const wflat = await ravel(w);
    const weighted = await multiply(flat, wflat);
    const wsum = await sum(wflat);
    const avg = await divide(await sum(weighted), wsum);
    if (returned) {
      return [avg.item(), wsum.item()];
    }
    return avg.item();
  }

  // Weights along axis
  if (w.ndim === 1 && w.size === a.shape[axis]) {
    // Broadcast weights to match array shape along axis
    const broadcastShape = a.shape.map((s, i) => i === axis ? s : 1);
    const wBroad = reshape(w, broadcastShape);
    const weighted = await multiply(a, wBroad);
    const wsum = await sum(wBroad, axis);
    const avg = await divide(await sum(weighted, axis), wsum);
    if (returned) {
      return [avg, wsum];
    }
    return avg;
  }

  throw new Error('weights must be 1D array matching axis dimension');
}
```

### ptp() Implementation

```typescript
/**
 * Peak to peak (maximum - minimum) along axis.
 */
export async function ptp(
  a: NDArray,
  axis: number | null = null,
  keepdims: boolean = false
): Promise<NDArray | number> {
  const maxVal = await amax(a, axis, keepdims);
  const minVal = await amin(a, axis, keepdims);
  return subtract(maxVal, minVal);
}
```

### percentile() / quantile() Implementation

```typescript
/**
 * Compute q-th percentile of data along axis.
 */
export async function percentile(
  a: NDArray,
  q: number | number[],
  axis: number | null = null,
  interpolation: 'linear' | 'lower' | 'higher' | 'midpoint' | 'nearest' = 'linear',
  keepdims: boolean = false
): Promise<NDArray | number> {
  // Convert percentile (0-100) to quantile (0-1)
  const qNorm = Array.isArray(q) ? q.map(v => v / 100) : q / 100;
  return quantile(a, qNorm, axis, interpolation, keepdims);
}

/**
 * Compute q-th quantile of data along axis.
 */
export async function quantile(
  a: NDArray,
  q: number | number[],
  axis: number | null = null,
  interpolation: 'linear' | 'lower' | 'higher' | 'midpoint' | 'nearest' = 'linear',
  keepdims: boolean = false
): Promise<NDArray | number> {
  // Similar to nanquantile but without NaN filtering
  // Can reuse internal logic from nanquantile

  const qArr = Array.isArray(q) ? q : [q];

  // Validate q values
  for (const qi of qArr) {
    if (qi < 0 || qi > 1) {
      throw new Error('quantile values must be between 0 and 1');
    }
  }

  if (axis === null) {
    const sorted = await sort(ravel(a));
    return computeQuantiles(sorted, qArr, interpolation);
  }

  // Apply along axis
  return applyAlongAxis(
    async (arr1d) => computeQuantiles(await sort(arr1d), qArr, interpolation),
    axis,
    a,
    keepdims
  );
}

function computeQuantiles(
  sorted: NDArray,
  q: number[],
  interpolation: string
): NDArray | number {
  const n = sorted.size;
  const results: number[] = [];

  for (const qi of q) {
    const idx = qi * (n - 1);
    const lo = Math.floor(idx);
    const hi = Math.ceil(idx);
    const frac = idx - lo;

    let val: number;
    switch (interpolation) {
      case 'lower':
        val = sorted.getFlat(lo);
        break;
      case 'higher':
        val = sorted.getFlat(hi);
        break;
      case 'nearest':
        val = sorted.getFlat(frac < 0.5 ? lo : hi);
        break;
      case 'midpoint':
        val = (sorted.getFlat(lo) + sorted.getFlat(hi)) / 2;
        break;
      case 'linear':
      default:
        val = sorted.getFlat(lo) * (1 - frac) + sorted.getFlat(hi) * frac;
    }
    results.push(val);
  }

  return results.length === 1 ? results[0] : fromArray(results);
}
```

### gradient() Implementation

```typescript
/**
 * Compute N-dimensional gradient using finite differences.
 */
export async function gradient(
  f: NDArray,
  options: {
    spacing?: number | number[] | NDArray[],
    axis?: number | number[] | null,
    edge_order?: 1 | 2
  } = {}
): Promise<NDArray | NDArray[]> {
  const { spacing = 1, axis = null, edge_order = 1 } = options;

  const axes = axis === null
    ? Array.from({ length: f.ndim }, (_, i) => i)
    : Array.isArray(axis) ? axis : [axis];

  const results: NDArray[] = [];

  for (const ax of axes) {
    const n = f.shape[ax];
    const h = typeof spacing === 'number' ? spacing :
              Array.isArray(spacing) ? spacing[axes.indexOf(ax)] : 1;

    const grad = zeros(f.shape, f.dtype);

    // Interior points: central difference
    // grad[i] = (f[i+1] - f[i-1]) / (2h)

    // Edge points: forward/backward difference
    // edge_order=1: first-order accurate
    // edge_order=2: second-order accurate

    // Implementation involves slicing along axis...
    // This is complex - see NumPy source for full implementation
  }

  return results.length === 1 ? results[0] : results;
}
```

### trapezoid() Implementation

```typescript
/**
 * Integrate along axis using trapezoidal rule.
 */
export async function trapezoid(
  y: NDArray,
  x: NDArray | null = null,
  dx: number = 1.0,
  axis: number = -1
): Promise<NDArray | number> {
  // Normalize axis
  if (axis < 0) axis += y.ndim;

  if (x !== null) {
    // Non-uniform spacing
    // dx[i] = x[i+1] - x[i]
    // integral ≈ sum((y[i] + y[i+1]) / 2 * dx[i])
    const d = await diff(x, 1, axis);
    const ySlice1 = sliceArray(y, axis, 0, -1);  // y[:-1]
    const ySlice2 = sliceArray(y, axis, 1, null); // y[1:]
    const yAvg = await divide(await add(ySlice1, ySlice2), 2);
    return sum(await multiply(yAvg, d), axis);
  }

  // Uniform spacing
  // integral ≈ dx * (y[0]/2 + y[1] + y[2] + ... + y[n-1]/2)
  // = dx * (sum(y) - (y[0] + y[n-1])/2)
  const total = await sum(y, axis);
  const first = sliceArray(y, axis, 0, 1);
  const last = sliceArray(y, axis, -1, null);
  const correction = await divide(await add(first, last), 2);
  return multiply(await subtract(total, await sum(correction, axis)), dx);
}
```

### correlate() Implementation

```typescript
/**
 * Cross-correlation of two 1D arrays.
 */
export async function correlate(
  a: NDArray,
  v: NDArray,
  mode: 'valid' | 'same' | 'full' = 'valid'
): Promise<NDArray> {
  if (a.ndim !== 1 || v.ndim !== 1) {
    throw new Error('correlate only supports 1D arrays');
  }

  const n = a.size;
  const m = v.size;

  // Reverse v for correlation (vs convolution)
  const vRev = await flip(v);

  let outSize: number;
  let offset: number;

  switch (mode) {
    case 'full':
      outSize = n + m - 1;
      offset = m - 1;
      break;
    case 'same':
      outSize = Math.max(n, m);
      offset = Math.floor(m / 2);
      break;
    case 'valid':
      outSize = Math.max(n, m) - Math.min(n, m) + 1;
      offset = Math.min(n, m) - 1;
      break;
  }

  const result = zeros([outSize], promoteTypes(a.dtype, v.dtype));

  // Direct computation (can be optimized with FFT for large arrays)
  for (let i = 0; i < outSize; i++) {
    let sum = 0;
    for (let j = 0; j < m; j++) {
      const aIdx = i - offset + j;
      if (aIdx >= 0 && aIdx < n) {
        sum += a.getFlat(aIdx) * vRev.getFlat(j);
      }
    }
    result.setFlat(i, sum);
  }

  return result;
}
```

## Testing

Create or extend `tests/ts/statistics.test.ts`:

```typescript
describe('Statistics functions', () => {
  describe('average', () => {
    it('computes unweighted average', async () => {
      const a = fromArray([1, 2, 3, 4]);
      expect(await average(a)).toBe(2.5);
    });

    it('computes weighted average', async () => {
      const a = fromArray([1, 2, 3, 4]);
      const w = fromArray([1, 1, 1, 1]);
      expect(await average(a, null, w)).toBe(2.5);

      const w2 = fromArray([4, 3, 2, 1]);
      expect(await average(a, null, w2)).toBe(2.0);
    });

    it('returns weights sum when returned=true', async () => {
      const a = fromArray([1, 2, 3, 4]);
      const w = fromArray([1, 2, 3, 4]);
      const [avg, wsum] = await average(a, null, w, true);
      expect(wsum).toBe(10);
    });
  });

  describe('ptp', () => {
    it('computes peak to peak', async () => {
      const a = fromArray([1, 5, 2, 8, 3]);
      expect(await ptp(a)).toBe(7);  // 8 - 1
    });

    it('works along axis', async () => {
      const a = fromArray([[1, 5], [2, 8]], [2, 2]);
      const result = await ptp(a, 0);
      expect(await result.toArray()).toEqual([1, 3]);
    });
  });

  describe('percentile', () => {
    it('computes median as 50th percentile', async () => {
      const a = fromArray([1, 2, 3, 4, 5]);
      expect(await percentile(a, 50)).toBe(3);
    });

    it('computes multiple percentiles', async () => {
      const a = fromArray([1, 2, 3, 4, 5]);
      const result = await percentile(a, [25, 50, 75]);
      // Linear interpolation
    });
  });

  describe('quantile', () => {
    it('computes 0.5 quantile as median', async () => {
      const a = fromArray([1, 2, 3, 4, 5]);
      expect(await quantile(a, 0.5)).toBe(3);
    });
  });

  describe('gradient', () => {
    it('computes 1D gradient', async () => {
      const f = fromArray([1, 2, 4, 7, 11]);
      const g = await gradient(f);
      // Central differences: [1, 1.5, 2.5, 3.5, 4]
    });

    it('handles non-uniform spacing', async () => {
      const f = fromArray([1, 4, 9]);
      const x = fromArray([0, 1, 3]);
      const g = await gradient(f, { spacing: x });
    });
  });

  describe('trapezoid', () => {
    it('integrates constant function', async () => {
      const y = fromArray([1, 1, 1, 1, 1]);
      expect(await trapezoid(y)).toBe(4);  // 5 points, dx=1
    });

    it('integrates linear function', async () => {
      const y = fromArray([0, 1, 2, 3, 4]);
      expect(await trapezoid(y)).toBe(8);  // Area of triangle
    });

    it('handles non-uniform x', async () => {
      const y = fromArray([1, 2, 3]);
      const x = fromArray([0, 1, 3]);
      const result = await trapezoid(y, x);
    });
  });

  describe('correlate', () => {
    it('computes valid correlation', async () => {
      const a = fromArray([1, 2, 3]);
      const v = fromArray([0, 1, 0.5]);
      const result = await correlate(a, v, 'valid');
    });

    it('computes full correlation', async () => {
      const a = fromArray([1, 2, 3]);
      const v = fromArray([1, 1]);
      const result = await correlate(a, v, 'full');
      expect(result.size).toBe(4);  // n + m - 1
    });
  });
});
```

## Exports to Add

In `src/ts/index.ts`:
```typescript
export { average, ptp, percentile, quantile, gradient, trapezoid, correlate } from './statistics';
```

## Priority

Medium:
- `average` - Commonly needed for weighted means
- `percentile/quantile` - Essential statistics (nan versions exist)
- `ptp` - Simple, commonly used
- `gradient` - Important for numerical computing
- `trapezoid` - Essential for integration
- `correlate` - Foundation for signal processing

## Estimated Scope

- ~400 lines of TypeScript
- ~200 lines of tests
- No WASM changes required initially
- `gradient` is the most complex (finite difference formulas)

## Notes

- `percentile` and `quantile` can share implementation with nanpercentile/nanquantile
- `correlate` can be optimized with FFT for large arrays (future enhancement)
- Consider adding `convolve` alongside `correlate` (they're related)
- `trapezoid` replaces deprecated `trapz` in NumPy 2.0
