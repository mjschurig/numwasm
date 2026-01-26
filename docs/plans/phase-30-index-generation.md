# Phase 30: Index Generation Functions

## Overview

Implement advanced index generation functions that provide slice-based grid generation and mask-based indexing.

## Functions to Implement

### 1. `mgrid` - Dense meshgrid from slices
```typescript
// NumPy uses slice notation: np.mgrid[0:5, 0:3]
// TypeScript equivalent using function call
function mgrid(...slices: SliceSpec[]): NDArray
```
- Creates dense multi-dimensional coordinate arrays
- Equivalent to `meshgrid(*xi, indexing='ij')` with `sparse=False`
- Returns stacked arrays where each has shape `(n0, n1, ..., nk)`

**NumPy Reference:** `numpy/lib/_index_tricks_impl.py` class `MGridClass`

### 2. `ogrid` - Open (sparse) meshgrid from slices
```typescript
function ogrid(...slices: SliceSpec[]): NDArray[]
```
- Creates open (sparse) multi-dimensional coordinate arrays
- Equivalent to `meshgrid(*xi, indexing='ij')` with `sparse=True`
- Returns arrays where each has shape with 1s except on its own axis

**NumPy Reference:** `numpy/lib/_index_tricks_impl.py` class `OGridClass`

### 3. `mask_indices(n, mask_func, k?)` - Indices from mask function
```typescript
function mask_indices(
  n: number,
  mask_func: (m: NDArray, k: number) => NDArray,
  k: number = 0
): [NDArray, NDArray]
```
- Returns indices to access elements where mask_func returns True
- `mask_func` should be a function like `triu` or `tril`
- Returns tuple of row and column indices

**NumPy Reference:** `numpy/lib/_index_tricks_impl.py` line ~800

## Implementation Strategy

### TypeScript API Design

Since TypeScript doesn't support Python's slice subscript notation (`mgrid[0:5:1]`), we need an alternative API:

**Option A: Function with slice specs**
```typescript
// Use objects for slice specs
mgrid({ start: 0, stop: 5 }, { start: 0, stop: 3 })

// Or helper function
mgrid(slice(0, 5), slice(0, 3))
```

**Option B: Variadic arguments with arrays**
```typescript
// [start, stop] or [start, stop, step]
mgrid([0, 5], [0, 3])
mgrid([0, 5, 0.5], [0, 3, 1])  // with step
```

**Recommended: Option B** - More concise and natural for TypeScript

### File Location
- Add to `src/ts/indexing.ts`
- Export from `src/ts/index.ts`

### Dependencies
- Existing `meshgrid()`, `indices()`, `arange()`, `linspace()`
- Existing `triu()`, `tril()` for mask_indices examples

## Implementation Details

### mgrid Implementation

```typescript
type GridSlice = [number, number] | [number, number, number];

/**
 * Create a dense multi-dimensional meshgrid.
 *
 * Similar to meshgrid with indexing='ij', but accepts slice-like specifications.
 *
 * @param slices - Array of [start, stop] or [start, stop, step] tuples
 * @returns Stacked array of shape (ndim, n0, n1, ..., nk-1)
 *
 * @example
 * // Create 2D grid from 0-2 and 0-3
 * const grid = mgrid([0, 3], [0, 4]);
 * // grid[0] contains row indices: [[0,0,0,0],[1,1,1,1],[2,2,2,2]]
 * // grid[1] contains col indices: [[0,1,2,3],[0,1,2,3],[0,1,2,3]]
 */
export function mgrid(...slices: GridSlice[]): NDArray {
  if (slices.length === 0) {
    throw new Error('mgrid requires at least one slice specification');
  }

  // Generate 1D arrays for each dimension
  const arrays: NDArray[] = slices.map(spec => {
    const [start, stop, step = 1] = spec;
    if (Number.isInteger(step) && Number.isInteger(start) && Number.isInteger(stop)) {
      return arange(start, stop, step);
    } else {
      // For float steps, use number of points
      const num = Math.ceil((stop - start) / step);
      return linspace(start, stop - step, num, false);
    }
  });

  // Create dense meshgrid with 'ij' indexing
  const grids = meshgrid(...arrays, { indexing: 'ij' });

  // Stack into single array
  return stack(grids, 0);
}
```

### ogrid Implementation

```typescript
/**
 * Create an open (sparse) multi-dimensional meshgrid.
 *
 * Each output array has shape with 1s except along its own axis.
 *
 * @param slices - Array of [start, stop] or [start, stop, step] tuples
 * @returns Array of sparse coordinate arrays
 *
 * @example
 * // Create sparse 2D grid
 * const [x, y] = ogrid([0, 3], [0, 4]);
 * // x.shape = [3, 1], y.shape = [1, 4]
 */
export function ogrid(...slices: GridSlice[]): NDArray[] {
  if (slices.length === 0) {
    throw new Error('ogrid requires at least one slice specification');
  }

  // Generate 1D arrays for each dimension
  const arrays: NDArray[] = slices.map(spec => {
    const [start, stop, step = 1] = spec;
    if (Number.isInteger(step) && Number.isInteger(start) && Number.isInteger(stop)) {
      return arange(start, stop, step);
    } else {
      const num = Math.ceil((stop - start) / step);
      return linspace(start, stop - step, num, false);
    }
  });

  // Create sparse meshgrid with 'ij' indexing
  return meshgrid(...arrays, { indexing: 'ij', sparse: true });
}
```

### mask_indices Implementation

```typescript
/**
 * Return indices to access (n, n) arrays based on a mask function.
 *
 * @param n - Size of the square array
 * @param mask_func - Function that takes (arr, k) and returns boolean mask
 * @param k - Diagonal offset passed to mask_func (default: 0)
 * @returns Tuple of [row_indices, col_indices]
 *
 * @example
 * // Get upper triangle indices
 * const [row, col] = mask_indices(3, triu_mask);
 * // row = [0, 0, 0, 1, 1, 2]
 * // col = [0, 1, 2, 1, 2, 2]
 */
export async function mask_indices(
  n: number,
  mask_func: (arr: NDArray, k: number) => NDArray | Promise<NDArray>,
  k: number = 0
): Promise<[NDArray, NDArray]> {
  // Create an n x n array of ones
  const m = ones([n, n]);

  // Apply mask function
  const mask = await mask_func(m, k);

  // Get indices where mask is non-zero
  return nonzero(mask) as Promise<[NDArray, NDArray]>;
}
```

### Helper: triu_mask and tril_mask

```typescript
// These can be used with mask_indices
function triu_mask(arr: NDArray, k: number = 0): NDArray {
  return triu(ones(arr.shape), k);
}

function tril_mask(arr: NDArray, k: number = 0): NDArray {
  return tril(ones(arr.shape), k);
}
```

## Complex Number Support in Slices

For compatibility with NumPy's complex step notation (`1j` means use number of points):

```typescript
type GridSlice =
  | [number, number]           // [start, stop], step=1
  | [number, number, number]   // [start, stop, step]
  | [number, number, number, 'j'];  // [start, stop, num, 'j'] for linspace-style

// Example: mgrid([0, 5, 10, 'j']) → linspace(0, 5, 10)
```

## Testing

Create or extend `tests/ts/indexing.test.ts`:

```typescript
describe('Grid generation', () => {
  describe('mgrid', () => {
    it('creates 1D dense grid', () => {
      const grid = mgrid([0, 5]);
      expect(grid.shape).toEqual([1, 5]);
      expect(await grid.toArray()).toEqual([[0, 1, 2, 3, 4]]);
    });

    it('creates 2D dense grid', () => {
      const grid = mgrid([0, 3], [0, 4]);
      expect(grid.shape).toEqual([2, 3, 4]);
      // grid[0] = row indices broadcasted
      // grid[1] = col indices broadcasted
    });

    it('supports custom step', () => {
      const grid = mgrid([0, 1, 0.5]);  // 0, 0.5
      expect(grid.shape).toEqual([1, 2]);
    });

    it('matches NumPy behavior', async () => {
      // Compare with known NumPy output
      const grid = mgrid([0, 2], [0, 3]);
      const expected = [
        [[0, 0, 0], [1, 1, 1]],  // row indices
        [[0, 1, 2], [0, 1, 2]]   // col indices
      ];
    });
  });

  describe('ogrid', () => {
    it('creates sparse grid', () => {
      const [x, y] = ogrid([0, 3], [0, 4]);
      expect(x.shape).toEqual([3, 1]);
      expect(y.shape).toEqual([1, 4]);
    });

    it('can broadcast to dense grid', async () => {
      const [x, y] = ogrid([0, 2], [0, 3]);
      const dense = add(x, multiply(y, 0));  // broadcast x
      expect(dense.shape).toEqual([2, 3]);
    });
  });

  describe('mask_indices', () => {
    it('returns upper triangle indices', async () => {
      const [row, col] = await mask_indices(3, (m, k) => triu(m, k));
      // Upper triangle of 3x3: (0,0), (0,1), (0,2), (1,1), (1,2), (2,2)
      expect(await row.toArray()).toEqual([0, 0, 0, 1, 1, 2]);
      expect(await col.toArray()).toEqual([0, 1, 2, 1, 2, 2]);
    });

    it('respects k parameter', async () => {
      const [row, col] = await mask_indices(3, (m, k) => triu(m, k), 1);
      // Strict upper triangle (k=1)
      expect(await row.toArray()).toEqual([0, 0, 1]);
      expect(await col.toArray()).toEqual([1, 2, 2]);
    });

    it('works with lower triangle', async () => {
      const [row, col] = await mask_indices(3, (m, k) => tril(m, k));
      // Lower triangle of 3x3
    });
  });
});
```

## Exports to Add

In `src/ts/index.ts`:
```typescript
export { mgrid, ogrid, mask_indices } from './indexing';
```

## API Compatibility Note

NumPy's `mgrid` and `ogrid` use special indexer objects with `__getitem__`:
```python
np.mgrid[0:5, 0:3]
np.ogrid[0:5:1j*10]  # complex step for linspace
```

Our TypeScript API uses function calls instead:
```typescript
mgrid([0, 5], [0, 3])
mgrid([0, 5, 10, 'j'])  // linspace-style
```

Document this difference clearly in API docs.

## Priority

Low - These are convenience functions. `meshgrid()` and `indices()` already cover most use cases.

## Estimated Scope

- ~100 lines of TypeScript
- ~80 lines of tests
- No WASM changes required
