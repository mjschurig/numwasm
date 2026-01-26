/**
 * NumJS Set Operations Tests
 *
 * Tests for unique, union1d, intersect1d, setdiff1d, setxor1d, isin, in1d, ediff1d
 */

import { describe, it, expect, beforeAll, afterEach } from 'vitest';
import {
  NDArray,
  loadWasmModule,
  unique,
  uniqueValues,
  uniqueIndex,
  uniqueInverse,
  uniqueCounts,
  uniqueAll,
  union1d,
  intersect1d,
  setdiff1d,
  setxor1d,
  isin,
  in1d,
  ediff1d,
} from 'numjs-wasm';

// Initialize WASM module
beforeAll(async () => {
  await loadWasmModule();
});

describe('unique', () => {
  it('should find unique values', async () => {
    const arr = await NDArray.fromArray([1, 2, 2, 3, 1, 3, 4]);
    const result = await uniqueValues(arr);

    expect(result.shape).toEqual([4]);
    expect(await result.toArray()).toEqual([1, 2, 3, 4]);

    arr.dispose();
    result.dispose();
  });

  it('should return unique values from 2D array (flattened)', async () => {
    const arr = await NDArray.fromArray([[1, 2], [2, 3]]);
    const result = await uniqueValues(arr);

    expect(result.shape).toEqual([3]);
    expect(await result.toArray()).toEqual([1, 2, 3]);

    arr.dispose();
    result.dispose();
  });

  it('should handle empty array', async () => {
    const arr = await NDArray.fromArray([]);
    const result = await uniqueValues(arr);

    expect(result.shape).toEqual([0]);
    expect(result.size).toBe(0);

    arr.dispose();
    result.dispose();
  });

  it('should return indices with returnIndex', async () => {
    const arr = await NDArray.fromArray([3, 1, 2, 1, 3]);
    const { values, indices } = await uniqueIndex(arr);

    expect(await values.toArray()).toEqual([1, 2, 3]);
    // First occurrence indices
    expect(await indices.toArray()).toEqual([1, 2, 0]);

    arr.dispose();
    values.dispose();
    indices.dispose();
  });

  it('should return inverse indices', async () => {
    const arr = await NDArray.fromArray([3, 1, 2, 1, 3]);
    const { values, inverse } = await uniqueInverse(arr);

    expect(await values.toArray()).toEqual([1, 2, 3]);
    // inverse[i] is the index in values that reconstructs arr[i]
    expect(await inverse.toArray()).toEqual([2, 0, 1, 0, 2]);

    arr.dispose();
    values.dispose();
    inverse.dispose();
  });

  it('should return counts', async () => {
    const arr = await NDArray.fromArray([3, 1, 2, 1, 3, 1]);
    const { values, counts } = await uniqueCounts(arr);

    expect(await values.toArray()).toEqual([1, 2, 3]);
    expect(await counts.toArray()).toEqual([3, 1, 2]); // 1 appears 3x, 2 appears 1x, 3 appears 2x

    arr.dispose();
    values.dispose();
    counts.dispose();
  });

  it('should return all arrays with uniqueAll', async () => {
    const arr = await NDArray.fromArray([3, 1, 2, 1, 3]);
    const result = await uniqueAll(arr);

    expect(await result.values.toArray()).toEqual([1, 2, 3]);
    expect(result.indices).toBeDefined();
    expect(result.inverse).toBeDefined();
    expect(result.counts).toBeDefined();

    arr.dispose();
    result.values.dispose();
    result.indices!.dispose();
    result.inverse!.dispose();
    result.counts!.dispose();
  });
});

describe('union1d', () => {
  it('should find union of two arrays', async () => {
    const a = await NDArray.fromArray([1, 2, 3]);
    const b = await NDArray.fromArray([2, 3, 4]);
    const result = await union1d(a, b);

    expect(result.shape).toEqual([4]);
    expect(await result.toArray()).toEqual([1, 2, 3, 4]);

    a.dispose();
    b.dispose();
    result.dispose();
  });

  it('should handle non-overlapping arrays', async () => {
    const a = await NDArray.fromArray([1, 2]);
    const b = await NDArray.fromArray([3, 4]);
    const result = await union1d(a, b);

    expect(await result.toArray()).toEqual([1, 2, 3, 4]);

    a.dispose();
    b.dispose();
    result.dispose();
  });

  it('should handle empty arrays', async () => {
    const a = await NDArray.fromArray([1, 2, 3]);
    const b = await NDArray.fromArray([]);
    const result = await union1d(a, b);

    expect(await result.toArray()).toEqual([1, 2, 3]);

    a.dispose();
    b.dispose();
    result.dispose();
  });
});

describe('intersect1d', () => {
  it('should find intersection of two arrays', async () => {
    const a = await NDArray.fromArray([1, 3, 4, 3]);
    const b = await NDArray.fromArray([3, 1, 2, 1]);
    const result = (await intersect1d(a, b)) as NDArray;

    expect(result.shape).toEqual([2]);
    expect(await result.toArray()).toEqual([1, 3]);

    a.dispose();
    b.dispose();
    result.dispose();
  });

  it('should handle arrays with no common elements', async () => {
    const a = await NDArray.fromArray([1, 2]);
    const b = await NDArray.fromArray([3, 4]);
    const result = (await intersect1d(a, b)) as NDArray;

    expect(result.shape).toEqual([0]);

    a.dispose();
    b.dispose();
    result.dispose();
  });

  it('should return indices when requested', async () => {
    const a = await NDArray.fromArray([1, 3, 4, 3]);
    const b = await NDArray.fromArray([3, 1, 2, 1]);
    const result = (await intersect1d(a, b, { returnIndices: true })) as {
      values: NDArray;
      indices1: NDArray;
      indices2: NDArray;
    };

    expect(await result.values.toArray()).toEqual([1, 3]);
    // indices1: where 1 and 3 first appear in a
    // indices2: where 1 and 3 first appear in b

    result.values.dispose();
    result.indices1.dispose();
    result.indices2.dispose();
    a.dispose();
    b.dispose();
  });
});

describe('setdiff1d', () => {
  it('should find set difference', async () => {
    const a = await NDArray.fromArray([1, 2, 3, 4]);
    const b = await NDArray.fromArray([2, 4]);
    const result = await setdiff1d(a, b);

    expect(await result.toArray()).toEqual([1, 3]);

    a.dispose();
    b.dispose();
    result.dispose();
  });

  it('should handle empty result', async () => {
    const a = await NDArray.fromArray([1, 2]);
    const b = await NDArray.fromArray([1, 2, 3, 4]);
    const result = await setdiff1d(a, b);

    expect(result.shape).toEqual([0]);

    a.dispose();
    b.dispose();
    result.dispose();
  });

  it('should handle empty second array', async () => {
    const a = await NDArray.fromArray([1, 2, 3]);
    const b = await NDArray.fromArray([]);
    const result = await setdiff1d(a, b);

    expect(await result.toArray()).toEqual([1, 2, 3]);

    a.dispose();
    b.dispose();
    result.dispose();
  });
});

describe('setxor1d', () => {
  it('should find symmetric difference', async () => {
    const a = await NDArray.fromArray([1, 2, 3]);
    const b = await NDArray.fromArray([2, 3, 4]);
    const result = await setxor1d(a, b);

    expect(await result.toArray()).toEqual([1, 4]);

    a.dispose();
    b.dispose();
    result.dispose();
  });

  it('should handle identical arrays', async () => {
    const a = await NDArray.fromArray([1, 2, 3]);
    const b = await NDArray.fromArray([1, 2, 3]);
    const result = await setxor1d(a, b);

    expect(result.shape).toEqual([0]);

    a.dispose();
    b.dispose();
    result.dispose();
  });

  it('should handle non-overlapping arrays', async () => {
    const a = await NDArray.fromArray([1, 2]);
    const b = await NDArray.fromArray([3, 4]);
    const result = await setxor1d(a, b);

    expect(await result.toArray()).toEqual([1, 2, 3, 4]);

    a.dispose();
    b.dispose();
    result.dispose();
  });
});

describe('isin', () => {
  it('should test membership', async () => {
    const arr = await NDArray.fromArray([2, 0, 1, 4]);
    const test = await NDArray.fromArray([0, 2]);
    const result = await isin(arr, test);

    expect(result.shape).toEqual([4]);
    // [2 in [0,2], 0 in [0,2], 1 in [0,2], 4 in [0,2]]
    // [true, true, false, false]
    const data = await result.toArray();
    expect(data).toEqual([1, 1, 0, 0]); // Boolean as 1/0

    arr.dispose();
    test.dispose();
    result.dispose();
  });

  it('should preserve shape of input', async () => {
    const arr = await NDArray.fromArray([[0, 2], [4, 6]]);
    const test = await NDArray.fromArray([1, 2, 4, 8]);
    const result = await isin(arr, test);

    expect(result.shape).toEqual([2, 2]);
    // [[0 in test, 2 in test], [4 in test, 6 in test]]
    // [[false, true], [true, false]]

    arr.dispose();
    test.dispose();
    result.dispose();
  });

  it('should invert result when requested', async () => {
    const arr = await NDArray.fromArray([2, 0, 1, 4]);
    const test = await NDArray.fromArray([0, 2]);
    const result = await isin(arr, test, { invert: true });

    // Inverted: [false, false, true, true]
    const data = await result.toArray();
    expect(data).toEqual([0, 0, 1, 1]);

    arr.dispose();
    test.dispose();
    result.dispose();
  });

  it('should handle empty test array', async () => {
    const arr = await NDArray.fromArray([1, 2, 3]);
    const test = await NDArray.fromArray([]);
    const result = await isin(arr, test);

    // Nothing is in empty array
    const data = await result.toArray();
    expect(data).toEqual([0, 0, 0]);

    arr.dispose();
    test.dispose();
    result.dispose();
  });
});

describe('in1d', () => {
  it('should test membership with flattened output', async () => {
    const arr = await NDArray.fromArray([[0, 2], [4, 6]]);
    const test = await NDArray.fromArray([1, 2, 4, 8]);
    const result = await in1d(arr, test);

    // Flattened: [0, 2, 4, 6] in [1, 2, 4, 8]
    // [false, true, true, false]
    expect(result.shape).toEqual([4]);

    arr.dispose();
    test.dispose();
    result.dispose();
  });
});

describe('ediff1d', () => {
  it('should compute differences between consecutive elements', async () => {
    const arr = await NDArray.fromArray([1, 2, 4, 7]);
    const result = await ediff1d(arr);

    expect(await result.toArray()).toEqual([1, 2, 3]);

    arr.dispose();
    result.dispose();
  });

  it('should handle prepend and append', async () => {
    const arr = await NDArray.fromArray([1, 2, 4, 7]);
    const result = await ediff1d(arr, { toPrepend: [0], toAppend: [10] });

    expect(await result.toArray()).toEqual([0, 1, 2, 3, 10]);

    arr.dispose();
    result.dispose();
  });

  it('should handle empty array', async () => {
    const arr = await NDArray.fromArray([]);
    const result = await ediff1d(arr);

    expect(result.size).toBe(0);

    arr.dispose();
    result.dispose();
  });

  it('should handle single element array', async () => {
    const arr = await NDArray.fromArray([5]);
    const result = await ediff1d(arr);

    expect(result.size).toBe(0); // n-1 differences for n elements

    arr.dispose();
    result.dispose();
  });
});
