/**
 * Level 1 Tests: Element Access, Shape Manipulation, Iterators
 */

import { describe, it, expect, beforeEach, afterEach } from 'vitest';
import {
  NDArray,
  DType,
  FlatIterator,
  nditer,
  ndenumerate,
  ndindex,
} from 'numjs-wasm';

// Track arrays for cleanup
const tracked: NDArray[] = [];

function track<T extends NDArray>(arr: T): T {
  tracked.push(arr);
  return arr;
}

beforeEach(() => {
  tracked.length = 0;
});

afterEach(() => {
  for (const arr of tracked) {
    if (!arr.isDisposed) {
      arr.dispose();
    }
  }
});

/* ============ Element Access Tests ============ */

describe('Element Access', () => {
  describe('get()', () => {
    it('should get element from 1D array', async () => {
      const arr = track(await NDArray.fromArray([1, 2, 3, 4, 5]));
      expect(arr.get(0)).toBe(1);
      expect(arr.get(2)).toBe(3);
      expect(arr.get(4)).toBe(5);
    });

    it('should get element from 2D array', async () => {
      const arr = track(await NDArray.fromArray([1, 2, 3, 4, 5, 6], [2, 3]));
      expect(arr.get(0, 0)).toBe(1);
      expect(arr.get(0, 1)).toBe(2);
      expect(arr.get(0, 2)).toBe(3);
      expect(arr.get(1, 0)).toBe(4);
      expect(arr.get(1, 1)).toBe(5);
      expect(arr.get(1, 2)).toBe(6);
    });

    it('should get element from 3D array', async () => {
      // 2x2x3 array
      const data = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12];
      const arr = track(await NDArray.fromArray(data, [2, 2, 3]));
      expect(arr.get(0, 0, 0)).toBe(1);
      expect(arr.get(0, 0, 2)).toBe(3);
      expect(arr.get(0, 1, 0)).toBe(4);
      expect(arr.get(1, 0, 0)).toBe(7);
      expect(arr.get(1, 1, 2)).toBe(12);
    });

    it('should throw on out-of-bounds index', async () => {
      const arr = track(await NDArray.fromArray([1, 2, 3]));
      expect(() => arr.get(5)).toThrow();
      expect(() => arr.get(-1)).toThrow();
    });

    it('should throw on wrong number of indices', async () => {
      const arr = track(await NDArray.fromArray([1, 2, 3, 4], [2, 2]));
      expect(() => arr.get(0)).toThrow();
      expect(() => arr.get(0, 0, 0)).toThrow();
    });
  });

  describe('set()', () => {
    it('should set element in 1D array', async () => {
      const arr = track(await NDArray.zeros([5]));
      arr.set(42, 2);
      expect(arr.get(2)).toBe(42);
      expect(arr.get(0)).toBe(0);
    });

    it('should set element in 2D array', async () => {
      const arr = track(await NDArray.zeros([2, 3]));
      arr.set(99, 1, 2);
      expect(arr.get(1, 2)).toBe(99);
      expect(arr.get(0, 0)).toBe(0);
    });

    it('should throw on out-of-bounds index', async () => {
      const arr = track(await NDArray.zeros([3]));
      expect(() => arr.set(42, 5)).toThrow();
    });
  });

  describe('item()', () => {
    it('should return value for single-element array', async () => {
      const arr = track(await NDArray.fromArray([42]));
      expect(arr.item()).toBe(42);
    });

    it('should work for 2D single-element array', async () => {
      const arr = track(await NDArray.fromArray([42], [1, 1]));
      expect(arr.item()).toBe(42);
    });

    it('should throw for multi-element array', async () => {
      const arr = track(await NDArray.fromArray([1, 2, 3]));
      expect(() => arr.item()).toThrow();
    });
  });

  describe('getFlat() / setFlat()', () => {
    it('should get element by flat index', async () => {
      const arr = track(await NDArray.fromArray([1, 2, 3, 4, 5, 6], [2, 3]));
      expect(arr.getFlat(0)).toBe(1);
      expect(arr.getFlat(3)).toBe(4);
      expect(arr.getFlat(5)).toBe(6);
    });

    it('should set element by flat index', async () => {
      const arr = track(await NDArray.zeros([2, 3]));
      arr.setFlat(3, 99);
      expect(arr.get(1, 0)).toBe(99);
    });

    it('should throw on out-of-bounds flat index', async () => {
      const arr = track(await NDArray.zeros([6]));
      expect(() => arr.getFlat(10)).toThrow();
      expect(() => arr.setFlat(-1, 0)).toThrow();
    });
  });

  describe('strides property', () => {
    it('should return correct strides for 1D array', async () => {
      const arr = track(await NDArray.zeros([5]));
      const strides = arr.strides;
      expect(strides).toHaveLength(1);
      expect(strides[0]).toBe(8); // Float64 = 8 bytes
    });

    it('should return correct strides for 2D array', async () => {
      const arr = track(await NDArray.zeros([2, 3]));
      const strides = arr.strides;
      expect(strides).toHaveLength(2);
      expect(strides[0]).toBe(24); // 3 * 8 bytes
      expect(strides[1]).toBe(8); // 8 bytes
    });

    it('should return correct strides for Float32', async () => {
      const arr = track(await NDArray.zeros([2, 3], { dtype: DType.Float32 }));
      const strides = arr.strides;
      expect(strides[0]).toBe(12); // 3 * 4 bytes
      expect(strides[1]).toBe(4); // 4 bytes
    });
  });
});

/* ============ Shape Manipulation Tests ============ */

describe('Shape Manipulation', () => {
  describe('reshape()', () => {
    it('should reshape 1D to 2D', async () => {
      const arr = track(await NDArray.fromArray([1, 2, 3, 4, 5, 6]));
      const reshaped = track(arr.reshape([2, 3]));

      expect(reshaped.shape).toEqual([2, 3]);
      expect(reshaped.ndim).toBe(2);
      expect(reshaped.size).toBe(6);

      // Check values
      expect(reshaped.get(0, 0)).toBe(1);
      expect(reshaped.get(0, 2)).toBe(3);
      expect(reshaped.get(1, 0)).toBe(4);
      expect(reshaped.get(1, 2)).toBe(6);
    });

    it('should reshape 2D to 1D', async () => {
      const arr = track(await NDArray.fromArray([1, 2, 3, 4, 5, 6], [2, 3]));
      const reshaped = track(arr.reshape([6]));

      expect(reshaped.shape).toEqual([6]);
      expect(reshaped.toArray()).toEqual([1, 2, 3, 4, 5, 6]);
    });

    it('should reshape with -1 for auto dimension', async () => {
      const arr = track(await NDArray.fromArray([1, 2, 3, 4, 5, 6]));

      const r1 = track(arr.reshape([2, -1]));
      expect(r1.shape).toEqual([2, 3]);

      const r2 = track(arr.reshape([-1, 2]));
      expect(r2.shape).toEqual([3, 2]);

      const r3 = track(arr.reshape([-1]));
      expect(r3.shape).toEqual([6]);
    });

    it('should be a view (share data)', async () => {
      const arr = track(await NDArray.fromArray([1, 2, 3, 4]));
      const reshaped = track(arr.reshape([2, 2]));

      expect(reshaped.isView).toBe(true);

      // Modify original
      arr.set(99, 0);
      expect(reshaped.get(0, 0)).toBe(99);
    });

    it('should throw on incompatible shape', async () => {
      const arr = track(await NDArray.fromArray([1, 2, 3, 4, 5]));
      expect(() => arr.reshape([2, 3])).toThrow();
    });

    it('should throw on multiple -1', async () => {
      const arr = track(await NDArray.fromArray([1, 2, 3, 4, 5, 6]));
      expect(() => arr.reshape([-1, -1])).toThrow();
    });
  });

  describe('T property (transpose)', () => {
    it('should transpose 2D array', async () => {
      const arr = track(await NDArray.fromArray([1, 2, 3, 4, 5, 6], [2, 3]));
      const t = track(arr.T);

      expect(t.shape).toEqual([3, 2]);

      // Original: [[1, 2, 3], [4, 5, 6]]
      // Transposed: [[1, 4], [2, 5], [3, 6]]
      expect(t.get(0, 0)).toBe(1);
      expect(t.get(0, 1)).toBe(4);
      expect(t.get(1, 0)).toBe(2);
      expect(t.get(1, 1)).toBe(5);
      expect(t.get(2, 0)).toBe(3);
      expect(t.get(2, 1)).toBe(6);
    });

    it('should be a view (share data)', async () => {
      const arr = track(await NDArray.fromArray([1, 2, 3, 4], [2, 2]));
      const t = track(arr.T);

      expect(t.isView).toBe(true);

      // Modify original
      arr.set(99, 0, 1);
      expect(t.get(1, 0)).toBe(99);
    });

    it('should transpose 3D array (reverse axes)', async () => {
      const arr = track(await NDArray.zeros([2, 3, 4]));
      const t = track(arr.T);

      expect(t.shape).toEqual([4, 3, 2]);
    });
  });

  describe('transpose() with axes', () => {
    it('should transpose with custom axes', async () => {
      const arr = track(await NDArray.zeros([2, 3, 4]));
      const t = track(arr.transpose([2, 0, 1]));

      expect(t.shape).toEqual([4, 2, 3]);
    });

    it('should throw on invalid axes', async () => {
      const arr = track(await NDArray.zeros([2, 3]));
      expect(() => arr.transpose([0, 0])).toThrow(); // Duplicate
      expect(() => arr.transpose([0, 5])).toThrow(); // Out of range
    });

    it('should throw on wrong number of axes', async () => {
      const arr = track(await NDArray.zeros([2, 3]));
      expect(() => arr.transpose([0])).toThrow();
    });
  });

  describe('ravel()', () => {
    it('should flatten to 1D', async () => {
      const arr = track(await NDArray.fromArray([1, 2, 3, 4, 5, 6], [2, 3]));
      const flat = track(arr.ravel());

      expect(flat.shape).toEqual([6]);
      expect(flat.toArray()).toEqual([1, 2, 3, 4, 5, 6]);
    });

    it('should be a view for contiguous arrays', async () => {
      const arr = track(await NDArray.fromArray([1, 2, 3, 4], [2, 2]));
      const flat = track(arr.ravel());

      expect(flat.isView).toBe(true);

      arr.set(99, 0, 0);
      expect(flat.get(0)).toBe(99);
    });
  });

  describe('flatten()', () => {
    it('should flatten to 1D', async () => {
      const arr = track(await NDArray.fromArray([1, 2, 3, 4, 5, 6], [2, 3]));
      const flat = track(arr.flatten());

      expect(flat.shape).toEqual([6]);
      expect(flat.toArray()).toEqual([1, 2, 3, 4, 5, 6]);
    });

    it('should create a copy (not view)', async () => {
      const arr = track(await NDArray.fromArray([1, 2, 3, 4], [2, 2]));
      const flat = track(arr.flatten());

      expect(flat.isView).toBe(false);

      arr.set(99, 0, 0);
      expect(flat.get(0)).toBe(1); // Unchanged
    });
  });

  describe('squeeze()', () => {
    it('should remove all size-1 dimensions', async () => {
      const arr = track(await NDArray.zeros([1, 3, 1, 4, 1]));
      const squeezed = track(arr.squeeze());

      expect(squeezed.shape).toEqual([3, 4]);
    });

    it('should squeeze specific axis', async () => {
      const arr = track(await NDArray.zeros([1, 3, 1, 4]));
      const squeezed = track(arr.squeeze(0));

      expect(squeezed.shape).toEqual([3, 1, 4]);
    });

    it('should handle negative axis', async () => {
      const arr = track(await NDArray.zeros([3, 1]));
      const squeezed = track(arr.squeeze(-1));

      expect(squeezed.shape).toEqual([3]);
    });

    it('should throw when axis is not size 1', async () => {
      const arr = track(await NDArray.zeros([2, 3]));
      expect(() => arr.squeeze(0)).toThrow();
    });

    it('should be a view', async () => {
      const arr = track(await NDArray.fromArray([1, 2, 3], [1, 3]));
      const squeezed = track(arr.squeeze());

      expect(squeezed.isView).toBe(true);
    });
  });

  describe('expandDims()', () => {
    it('should add dimension at start', async () => {
      const arr = track(await NDArray.zeros([3, 4]));
      const expanded = track(arr.expandDims(0));

      expect(expanded.shape).toEqual([1, 3, 4]);
    });

    it('should add dimension in middle', async () => {
      const arr = track(await NDArray.zeros([3, 4]));
      const expanded = track(arr.expandDims(1));

      expect(expanded.shape).toEqual([3, 1, 4]);
    });

    it('should add dimension at end', async () => {
      const arr = track(await NDArray.zeros([3, 4]));
      const expanded = track(arr.expandDims(2));

      expect(expanded.shape).toEqual([3, 4, 1]);
    });

    it('should handle negative axis', async () => {
      const arr = track(await NDArray.zeros([3, 4]));
      const expanded = track(arr.expandDims(-1));

      expect(expanded.shape).toEqual([3, 4, 1]);
    });

    it('should be a view', async () => {
      const arr = track(await NDArray.fromArray([1, 2, 3]));
      const expanded = track(arr.expandDims(0));

      expect(expanded.isView).toBe(true);
    });
  });

  describe('swapaxes()', () => {
    it('should swap two axes', async () => {
      const arr = track(await NDArray.zeros([2, 3, 4]));
      const swapped = track(arr.swapaxes(0, 2));

      expect(swapped.shape).toEqual([4, 3, 2]);
    });

    it('should handle negative axes', async () => {
      const arr = track(await NDArray.zeros([2, 3, 4]));
      const swapped = track(arr.swapaxes(-1, -3));

      expect(swapped.shape).toEqual([4, 3, 2]);
    });

    it('should be a view', async () => {
      const arr = track(await NDArray.fromArray([1, 2, 3, 4, 5, 6], [2, 3]));
      const swapped = track(arr.swapaxes(0, 1));

      expect(swapped.isView).toBe(true);
    });

    it('should throw on invalid axis', async () => {
      const arr = track(await NDArray.zeros([2, 3]));
      expect(() => arr.swapaxes(0, 5)).toThrow();
    });
  });

  describe('copy()', () => {
    it('should create independent copy', async () => {
      const arr = track(await NDArray.fromArray([1, 2, 3, 4], [2, 2]));
      const copy = track(arr.copy());

      expect(copy.shape).toEqual([2, 2]);
      expect(copy.toArray()).toEqual([1, 2, 3, 4]);
      expect(copy.isView).toBe(false);

      // Modify original
      arr.set(99, 0, 0);
      expect(copy.get(0, 0)).toBe(1); // Unchanged
    });
  });
});

/* ============ Iterator Tests ============ */

describe('Iterators', () => {
  describe('flat property (FlatIterator)', () => {
    it('should iterate over elements in row-major order', async () => {
      const arr = track(await NDArray.fromArray([1, 2, 3, 4, 5, 6], [2, 3]));
      const values: number[] = [];

      for (const x of arr.flat) {
        values.push(x);
      }

      expect(values).toEqual([1, 2, 3, 4, 5, 6]);
    });

    it('should work with 3D array', async () => {
      // 2x2x2 = 8 elements
      const data = [1, 2, 3, 4, 5, 6, 7, 8];
      const arr = track(await NDArray.fromArray(data, [2, 2, 2]));
      const values: number[] = [];

      for (const x of arr.flat) {
        values.push(x);
      }

      expect(values).toEqual(data);
    });

    it('should have index and reset methods', async () => {
      const arr = track(await NDArray.fromArray([1, 2, 3]));
      const flat = arr.flat;

      flat.next();
      expect(flat.index).toBe(1);

      flat.next();
      expect(flat.index).toBe(2);

      flat.reset();
      expect(flat.index).toBe(0);
    });
  });

  describe('nditer()', () => {
    it('should iterate over elements', async () => {
      const arr = track(await NDArray.fromArray([1, 2, 3, 4], [2, 2]));
      const values = [...nditer(arr)];

      expect(values).toEqual([1, 2, 3, 4]);
    });
  });

  describe('ndenumerate()', () => {
    it('should yield index-value pairs', async () => {
      const arr = track(await NDArray.fromArray([1, 2, 3, 4], [2, 2]));
      const pairs = [...ndenumerate(arr)];

      expect(pairs).toEqual([
        [[0, 0], 1],
        [[0, 1], 2],
        [[1, 0], 3],
        [[1, 1], 4],
      ]);
    });

    it('should work with 1D array', async () => {
      const arr = track(await NDArray.fromArray([10, 20, 30]));
      const pairs = [...ndenumerate(arr)];

      expect(pairs).toEqual([
        [[0], 10],
        [[1], 20],
        [[2], 30],
      ]);
    });
  });

  describe('ndindex()', () => {
    it('should generate indices for shape', () => {
      const indices = [...ndindex(2, 3)];

      expect(indices).toEqual([
        [0, 0],
        [0, 1],
        [0, 2],
        [1, 0],
        [1, 1],
        [1, 2],
      ]);
    });

    it('should handle 1D shape', () => {
      const indices = [...ndindex(4)];

      expect(indices).toEqual([[0], [1], [2], [3]]);
    });

    it('should handle 3D shape', () => {
      const indices = [...ndindex(2, 2, 2)];

      expect(indices).toHaveLength(8);
      expect(indices[0]).toEqual([0, 0, 0]);
      expect(indices[7]).toEqual([1, 1, 1]);
    });

    it('should handle empty shape', () => {
      const indices = [...ndindex()];
      expect(indices).toEqual([]);
    });
  });
});

/* ============ View Data Sharing Tests ============ */

describe('View Data Sharing', () => {
  it('reshape view should share data with original', async () => {
    const arr = track(await NDArray.fromArray([1, 2, 3, 4, 5, 6]));
    const view = track(arr.reshape([2, 3]));

    arr.set(100, 0);
    expect(view.get(0, 0)).toBe(100);

    view.set(200, 1, 2);
    expect(arr.get(5)).toBe(200);
  });

  it('transpose view should share data with original', async () => {
    const arr = track(await NDArray.fromArray([1, 2, 3, 4], [2, 2]));
    const view = track(arr.T);

    arr.set(100, 0, 0);
    expect(view.get(0, 0)).toBe(100);

    // Original [1,2; 3,4] -> Transposed [1,3; 2,4]
    // arr[0,1] = 2 -> view[1,0] = 2
    arr.set(200, 0, 1);
    expect(view.get(1, 0)).toBe(200);
  });

  it('chained views should all share data', async () => {
    const arr = track(await NDArray.fromArray([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], [2, 2, 3]));
    const view1 = track(arr.reshape([4, 3]));
    const view2 = track(view1.T);

    arr.set(999, 0, 0, 0);
    expect(view1.get(0, 0)).toBe(999);
    expect(view2.get(0, 0)).toBe(999);
  });
});

/* ============ NumPy Compatibility Tests ============ */
/* Tests ported from NumPy's test_numeric.py to ensure behavioral compatibility */

describe('NumPy Compatibility', () => {
  describe('ravel() - from test_numeric.py:test_ravel', () => {
    it('should flatten nested array to 1D', async () => {
      // NumPy: a = [[1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12]]
      // tgt = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
      const arr = track(await NDArray.fromArray([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], [4, 3]));
      const raveled = track(arr.ravel());

      expect(raveled.shape).toEqual([12]);
      expect(raveled.toArray()).toEqual([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]);
    });
  });

  describe('reshape() - from test_numeric.py:test_reshape', () => {
    it('should reshape array', async () => {
      // NumPy: arr = [[1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12]]
      // tgt = [[1, 2, 3, 4, 5, 6], [7, 8, 9, 10, 11, 12]]
      const arr = track(await NDArray.fromArray([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], [4, 3]));
      const reshaped = track(arr.reshape([2, 6]));

      expect(reshaped.shape).toEqual([2, 6]);
      // Row 0: [1, 2, 3, 4, 5, 6]
      expect(reshaped.get(0, 0)).toBe(1);
      expect(reshaped.get(0, 5)).toBe(6);
      // Row 1: [7, 8, 9, 10, 11, 12]
      expect(reshaped.get(1, 0)).toBe(7);
      expect(reshaped.get(1, 5)).toBe(12);
    });

    it('should share memory with original (view)', async () => {
      // From test_reshape_copy_arg: assert np.shares_memory(np.reshape(arr, shape), arr)
      const arr = track(await NDArray.arange(24));
      const reshaped = track(arr.reshape([2, 3, 4]));

      expect(reshaped.isView).toBe(true);

      // Modifying one affects the other
      arr.set(999, 0);
      expect(reshaped.get(0, 0, 0)).toBe(999);
    });
  });

  describe('squeeze() - from test_numeric.py:test_squeeze', () => {
    it('should remove all size-1 dimensions', async () => {
      // NumPy: A = [[[1, 1, 1], [2, 2, 2], [3, 3, 3]]]
      // np.squeeze(A).shape == (3, 3)
      const arr = track(await NDArray.fromArray([1, 1, 1, 2, 2, 2, 3, 3, 3], [1, 3, 3]));
      const squeezed = track(arr.squeeze());

      expect(squeezed.shape).toEqual([3, 3]);
    });

    it('should remove size-1 from [1, 3, 1]', async () => {
      // np.squeeze(np.zeros((1, 3, 1))).shape == (3,)
      const arr = track(await NDArray.zeros([1, 3, 1]));
      const squeezed = track(arr.squeeze());

      expect(squeezed.shape).toEqual([3]);
    });

    it('should squeeze specific axis', async () => {
      // np.squeeze(np.zeros((1, 3, 1)), axis=0).shape == (3, 1)
      const arr = track(await NDArray.zeros([1, 3, 1]));
      const squeezed = track(arr.squeeze(0));

      expect(squeezed.shape).toEqual([3, 1]);
    });

    it('should squeeze last axis with negative index', async () => {
      // np.squeeze(np.zeros((1, 3, 1)), axis=-1).shape == (1, 3)
      const arr = track(await NDArray.zeros([1, 3, 1]));
      const squeezed = track(arr.squeeze(-1));

      expect(squeezed.shape).toEqual([1, 3]);
    });

    it('should squeeze axis 2', async () => {
      // np.squeeze(np.zeros((1, 3, 1)), axis=2).shape == (1, 3)
      const arr = track(await NDArray.zeros([1, 3, 1]));
      const squeezed = track(arr.squeeze(2));

      expect(squeezed.shape).toEqual([1, 3]);
    });
  });

  describe('swapaxes() - from test_numeric.py:test_swapaxes', () => {
    it('should swap axes correctly', async () => {
      // NumPy: tgt = [[[0, 4], [2, 6]], [[1, 5], [3, 7]]]
      // a = [[[0, 1], [2, 3]], [[4, 5], [6, 7]]]
      // out = np.swapaxes(a, 0, 2)
      const arr = track(await NDArray.fromArray([0, 1, 2, 3, 4, 5, 6, 7], [2, 2, 2]));
      const swapped = track(arr.swapaxes(0, 2));

      // Verify shape changed
      expect(swapped.shape).toEqual([2, 2, 2]);

      // Verify values match NumPy's tgt
      // tgt[0][0][0] = 0, tgt[0][0][1] = 4
      expect(swapped.get(0, 0, 0)).toBe(0);
      expect(swapped.get(0, 0, 1)).toBe(4);
      // tgt[0][1][0] = 2, tgt[0][1][1] = 6
      expect(swapped.get(0, 1, 0)).toBe(2);
      expect(swapped.get(0, 1, 1)).toBe(6);
      // tgt[1][0][0] = 1, tgt[1][0][1] = 5
      expect(swapped.get(1, 0, 0)).toBe(1);
      expect(swapped.get(1, 0, 1)).toBe(5);
      // tgt[1][1][0] = 3, tgt[1][1][1] = 7
      expect(swapped.get(1, 1, 0)).toBe(3);
      expect(swapped.get(1, 1, 1)).toBe(7);
    });
  });

  describe('transpose() - from test_numeric.py:test_transpose', () => {
    it('should transpose 2D array', async () => {
      // NumPy: arr = [[1, 2], [3, 4], [5, 6]]
      // tgt = [[1, 3, 5], [2, 4, 6]]
      // np.transpose(arr, (1, 0)) == tgt
      const arr = track(await NDArray.fromArray([1, 2, 3, 4, 5, 6], [3, 2]));
      const transposed = track(arr.transpose([1, 0]));

      expect(transposed.shape).toEqual([2, 3]);
      // Row 0: [1, 3, 5]
      expect(transposed.get(0, 0)).toBe(1);
      expect(transposed.get(0, 1)).toBe(3);
      expect(transposed.get(0, 2)).toBe(5);
      // Row 1: [2, 4, 6]
      expect(transposed.get(1, 0)).toBe(2);
      expect(transposed.get(1, 1)).toBe(4);
      expect(transposed.get(1, 2)).toBe(6);
    });

    it('should handle negative axis indices', async () => {
      // np.transpose(arr, (-1, -2)) should equal np.transpose(arr, (1, 0))
      const arr = track(await NDArray.fromArray([1, 2, 3, 4, 5, 6], [3, 2]));
      const transposed = track(arr.transpose([-1, -2]));

      expect(transposed.shape).toEqual([2, 3]);
      expect(transposed.get(0, 0)).toBe(1);
      expect(transposed.get(0, 1)).toBe(3);
    });
  });

  describe('flatten() vs ravel()', () => {
    it('flatten should always return a copy', async () => {
      const arr = track(await NDArray.fromArray([1, 2, 3, 4], [2, 2]));
      const flattened = track(arr.flatten());

      expect(flattened.isView).toBe(false);

      // Modifying flattened should not affect original
      flattened.set(999, 0);
      expect(arr.get(0, 0)).toBe(1);
    });

    it('ravel should return a view when contiguous', async () => {
      const arr = track(await NDArray.fromArray([1, 2, 3, 4], [2, 2]));
      const raveled = track(arr.ravel());

      expect(raveled.isView).toBe(true);

      // Modifying raveled should affect original
      raveled.set(999, 0);
      expect(arr.get(0, 0)).toBe(999);
    });
  });

  describe('reshape with -1 dimension', () => {
    it('should infer dimension from -1', async () => {
      // arr = np.arange(12)
      // arr.reshape(-1) -> [12]
      // arr.reshape(3, -1) -> [3, 4]
      // arr.reshape(-1, 4) -> [3, 4]
      // arr.reshape(2, -1, 2) -> [2, 3, 2]
      const arr = track(await NDArray.arange(12));

      const r1 = track(arr.reshape([-1]));
      expect(r1.shape).toEqual([12]);

      const r2 = track(arr.reshape([3, -1]));
      expect(r2.shape).toEqual([3, 4]);

      const r3 = track(arr.reshape([-1, 4]));
      expect(r3.shape).toEqual([3, 4]);

      const r4 = track(arr.reshape([2, -1, 2]));
      expect(r4.shape).toEqual([2, 3, 2]);
    });

    it('should throw for invalid -1 reshape', async () => {
      const arr = track(await NDArray.arange(12));

      // Cannot have multiple -1
      expect(() => arr.reshape([-1, -1])).toThrow();

      // Size must be divisible
      expect(() => arr.reshape([5, -1])).toThrow();
    });
  });

  describe('expandDims() - np.expand_dims equivalent', () => {
    it('should add dimension at axis 0', async () => {
      const arr = track(await NDArray.fromArray([1, 2, 3]));
      const expanded = track(arr.expandDims(0));

      expect(expanded.shape).toEqual([1, 3]);
    });

    it('should add dimension at axis 1', async () => {
      const arr = track(await NDArray.fromArray([1, 2, 3]));
      const expanded = track(arr.expandDims(1));

      expect(expanded.shape).toEqual([3, 1]);
    });

    it('should add dimension at negative axis', async () => {
      const arr = track(await NDArray.fromArray([1, 2, 3]));
      const expanded = track(arr.expandDims(-1));

      expect(expanded.shape).toEqual([3, 1]);
    });
  });
});
