# Phase 29: Advanced Sorting Functions

## Overview

Implement specialized sorting functions that extend the basic `sort()` and `argsort()` already available.

## Functions to Implement

### 1. `lexsort(keys, axis?)` - Indirect lexicographic sort
```typescript
function lexsort(keys: NDArray[], axis: number = -1): NDArray
```
- Performs indirect stable sort using a sequence of keys
- Last key is the primary sort key, first key is the last resort
- Returns indices that would sort the arrays
- Essential for sorting structured data (like database records)

**NumPy Reference:** `numpy/_core/fromnumeric.py` line ~1100

### 2. `sort_complex(a)` - Sort complex numbers
```typescript
function sort_complex(a: NDArray): NDArray
```
- Sorts complex numbers first by real part, then by imaginary part
- Returns a sorted copy (not in-place)
- Uses lexsort internally: `lexsort([imag(a), real(a)])`

**NumPy Reference:** `numpy/lib/_function_base_impl.py` line ~700

### 3. `msort(a)` - Sort along first axis
```typescript
function msort(a: NDArray): NDArray
```
- Equivalent to `sort(a, axis=0)`
- Returns a sorted copy
- Deprecated in NumPy 2.0 but still used

**NumPy Reference:** `numpy/lib/_function_base_impl.py` line ~730

## Implementation Strategy

### File Location
- Add to `src/ts/sorting.ts`
- Export from `src/ts/index.ts`

### Dependencies
- Existing `sort()`, `argsort()` in `src/ts/sorting.ts`
- `real()`, `imag()` from Phase 28 (for sort_complex)
- No WASM changes required

## Implementation Details

### lexsort() Implementation

```typescript
/**
 * Perform indirect stable sort using sequence of keys.
 *
 * The last key in the sequence is used for the primary sort order,
 * ties are broken by the second-last key, etc.
 *
 * @param keys - Sequence of arrays to sort by
 * @param axis - Axis along which to sort (default: -1, last axis)
 * @returns Indices that would sort the arrays
 */
export async function lexsort(
  keys: NDArray[],
  axis: number = -1
): Promise<NDArray> {
  if (keys.length === 0) {
    throw new Error('lexsort requires at least one key');
  }

  // All keys must have same shape
  const shape = keys[0].shape;
  for (const key of keys) {
    if (!arraysEqual(key.shape, shape)) {
      throw new Error('all keys must have the same shape');
    }
  }

  // Normalize axis
  const ndim = keys[0].ndim;
  if (axis < 0) axis += ndim;

  const n = shape[axis];

  // For 1D arrays, simple approach
  if (ndim === 1) {
    // Create array of indices
    const indices = Array.from({ length: n }, (_, i) => i);

    // Sort indices using comparison function
    indices.sort((i, j) => {
      // Compare using keys from last to first
      for (let k = keys.length - 1; k >= 0; k--) {
        const a = keys[k].getFlat(i);
        const b = keys[k].getFlat(j);
        if (a < b) return -1;
        if (a > b) return 1;
      }
      return 0;  // All keys equal, maintain original order (stable)
    });

    return fromArray(indices, [n], 'int64');
  }

  // For N-D arrays, apply along axis
  // ... (more complex implementation)
}
```

### sort_complex() Implementation

```typescript
/**
 * Sort a complex array by real part, then imaginary part.
 *
 * @param a - Complex array to sort
 * @returns Sorted copy of the array
 */
export async function sort_complex(a: NDArray): Promise<NDArray> {
  const flat = await ravel(a);

  // Get real and imaginary parts
  const re = real(flat);
  const im = imag(flat);

  // Use lexsort: sort by imaginary (secondary), then real (primary)
  const indices = await lexsort([im, re]);

  // Apply indices to get sorted array
  const sorted = await take(flat, indices);

  // Reshape to original shape if needed
  if (a.ndim > 1) {
    return reshape(sorted, a.shape);
  }

  return sorted;
}
```

### msort() Implementation

```typescript
/**
 * Sort array along first axis.
 *
 * Equivalent to np.sort(a, axis=0).
 *
 * @param a - Array to sort
 * @returns Sorted copy of array
 */
export async function msort(a: NDArray): Promise<NDArray> {
  return sort(a, 0);
}
```

## Edge Cases

### lexsort
- Empty keys array → Error
- Keys with different shapes → Error
- Single key → Equivalent to argsort
- NaN handling: NaN values sort to the end

### sort_complex
- Real-valued input: Works (imaginary part is 0)
- NaN in real or imaginary part: Follows NumPy NaN handling

### msort
- 0-D array → Return copy
- 1-D array → Equivalent to sort(a)

## Testing

Create or extend `tests/ts/sorting.test.ts`:

```typescript
describe('Advanced sorting', () => {
  describe('lexsort', () => {
    it('sorts by multiple keys', async () => {
      // Sort by last name (primary), first name (secondary)
      const firstNames = fromArray(['John', 'Jane', 'John', 'Jane']);
      const lastNames = fromArray(['Doe', 'Doe', 'Smith', 'Smith']);

      const indices = await lexsort([firstNames, lastNames]);
      // Expected order: Jane Doe, John Doe, Jane Smith, John Smith
      // Indices: [1, 0, 3, 2]
    });

    it('uses last key as primary', async () => {
      const a = fromArray([1, 2, 1, 2]);
      const b = fromArray([4, 3, 2, 1]);

      const indices = await lexsort([a, b]);
      // Sort by b (primary), then a
      // b: [4, 3, 2, 1] → sorted: [1, 2, 3, 4] → indices: [3, 2, 1, 0]
    });

    it('maintains stability', async () => {
      const a = fromArray([1, 1, 1]);
      const b = fromArray([2, 2, 2]);

      const indices = await lexsort([a, b]);
      // All equal, should maintain original order: [0, 1, 2]
    });
  });

  describe('sort_complex', () => {
    it('sorts by real then imaginary', async () => {
      // [3+1i, 1+2i, 1+1i, 2+0i]
      // Sorted: [1+1i, 1+2i, 2+0i, 3+1i]
    });

    it('handles real arrays', async () => {
      const a = fromArray([3, 1, 2]);
      const sorted = await sort_complex(a);
      expect(await sorted.toArray()).toEqual([1, 2, 3]);
    });
  });

  describe('msort', () => {
    it('sorts along first axis', async () => {
      const a = fromArray([[3, 1], [2, 4]], [2, 2]);
      const sorted = await msort(a);
      // Each column sorted: [[2, 1], [3, 4]]
    });

    it('is equivalent to sort(a, 0)', async () => {
      const a = fromArray([[3, 1, 2], [1, 2, 3]], [2, 3]);
      const m = await msort(a);
      const s = await sort(a, 0);
      expect(await m.toArray()).toEqual(await s.toArray());
    });
  });
});
```

## Exports to Add

In `src/ts/index.ts`:
```typescript
export { lexsort, sort_complex, msort } from './sorting';
```

## Priority

Medium-Low:
- `lexsort` is useful for multi-key sorting
- `sort_complex` needed once complex accessors are added
- `msort` is deprecated in NumPy but included for compatibility

## Estimated Scope

- ~150 lines of TypeScript
- ~100 lines of tests
- Depends on Phase 28 for sort_complex (can stub or skip initially)

## Notes

- `lexsort` should use a stable sort algorithm (mergesort)
- `sort_complex` depends on `real()` and `imag()` from Phase 28
- Consider deprecation warning for `msort` to match NumPy 2.0
