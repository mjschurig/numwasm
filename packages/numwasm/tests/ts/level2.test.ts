/**
 * Level 2 Tests: Views, Slicing & Broadcasting
 */

import { describe, it, expect, beforeAll } from 'vitest';
import {
  NDArray,
  DType,
  slice,
  ellipsis,
  newaxis,
  Slice,
  broadcastShapes,
  broadcastTo,
  broadcastArrays,
  take,
  put,
  nonzero,
  where,
  countNonzero,
  extract,
  compress,
  diagonal,
  CLIP_WRAP,
  CLIP_CLIP,
  // New Phase 2 & 3 functions
  atleast_1d,
  atleast_2d,
  atleast_3d,
  argwhere,
  indices,
  ix_,
  diag_indices,
  tril_indices,
  triu_indices,
  take_along_axis,
  put_along_axis,
  putmask,
  place,
  select,
} from '../../dist/numjs.mjs';

describe('Level 2: Slicing', () => {
  describe('Slice class', () => {
    it('should create slice with start, stop, step', () => {
      const s = slice(1, 5, 2);
      expect(s.start).toBe(1);
      expect(s.stop).toBe(5);
      expect(s.step).toBe(2);
    });

    it('should handle null values', () => {
      const s = slice(null, 5);
      expect(s.start).toBe(null);
      expect(s.stop).toBe(5);
      expect(s.step).toBe(null);
    });

    it('should throw on step=0', () => {
      expect(() => slice(0, 5, 0)).toThrow('step cannot be zero');
    });

    it('should compute indices correctly', () => {
      const s = new Slice(1, 5, 2);
      const [start, stop, step, len] = s.indices(10);
      expect(start).toBe(1);
      expect(stop).toBe(5);
      expect(step).toBe(2);
      expect(len).toBe(2); // indices 1, 3
    });

    it('should handle negative indices', () => {
      const s = new Slice(-3, null);
      const [start, stop, step, len] = s.indices(10);
      expect(start).toBe(7);
      expect(stop).toBe(10);
      expect(step).toBe(1);
      expect(len).toBe(3);
    });

    it('should handle reverse iteration', () => {
      const s = new Slice(null, null, -1);
      const [start, stop, step, len] = s.indices(5);
      expect(start).toBe(4);
      expect(step).toBe(-1);
      expect(len).toBe(5);
    });
  });

  describe('NDArray.at()', () => {
    it('should get row from 2D array', async () => {
      const arr = await NDArray.fromArray([
        [1, 2, 3],
        [4, 5, 6],
      ]);
      const row = arr.at(0);
      expect(row.toArray()).toEqual([1, 2, 3]);
      expect(row.shape).toEqual([3]);
      arr.dispose();
      row.dispose();
    });

    it('should handle negative index', async () => {
      const arr = await NDArray.fromArray([
        [1, 2, 3],
        [4, 5, 6],
      ]);
      const lastRow = arr.at(-1);
      expect(lastRow.toArray()).toEqual([4, 5, 6]);
      arr.dispose();
      lastRow.dispose();
    });

    it('should throw on out of bounds', async () => {
      const arr = await NDArray.fromArray([[1, 2], [3, 4]]);
      expect(() => arr.at(5)).toThrow();
      arr.dispose();
    });
  });

  describe('NDArray.slice()', () => {
    it('should slice with integer index', async () => {
      const arr = await NDArray.fromArray([
        [1, 2, 3],
        [4, 5, 6],
      ]);
      const row = arr.slice([0]);
      expect(row.toArray()).toEqual([1, 2, 3]);
      arr.dispose();
      row.dispose();
    });

    it('should slice with Slice object', async () => {
      const arr = await NDArray.fromArray([1, 2, 3, 4, 5]);
      const sliced = arr.slice([slice(1, 4)]);
      expect(sliced.toArray()).toEqual([2, 3, 4]);
      arr.dispose();
      sliced.dispose();
    });

    it('should slice with step', async () => {
      const arr = await NDArray.fromArray([1, 2, 3, 4, 5]);
      const sliced = arr.slice([slice(null, null, 2)]);
      expect(sliced.toArray()).toEqual([1, 3, 5]);
      arr.dispose();
      sliced.dispose();
    });

    it('should handle multi-dimensional slicing', async () => {
      const arr = await NDArray.fromArray([
        [1, 2, 3],
        [4, 5, 6],
        [7, 8, 9],
      ]);
      // Get first column
      const col = arr.slice([slice(null), 0]);
      expect(col.toArray()).toEqual([1, 4, 7]);
      arr.dispose();
      col.dispose();
    });

    it('should handle newaxis', async () => {
      const arr = await NDArray.fromArray([1, 2, 3]);
      const expanded = arr.slice([newaxis]);
      expect(expanded.shape).toEqual([1, 3]);
      arr.dispose();
      expanded.dispose();
    });
  });
});

describe('Level 2: Broadcasting', () => {
  describe('broadcastShapes', () => {
    it('should broadcast compatible shapes', async () => {
      const result = await broadcastShapes([3, 1], [1, 4]);
      expect(result).toEqual([3, 4]);
    });

    it('should broadcast with different ndims', async () => {
      const result = await broadcastShapes([2, 3], [3]);
      expect(result).toEqual([2, 3]);
    });

    it('should return null for incompatible shapes', async () => {
      const result = await broadcastShapes([2, 3], [4]);
      expect(result).toBeNull();
    });

    it('should handle scalar broadcasting', async () => {
      const result = await broadcastShapes([5, 3], []);
      expect(result).toEqual([5, 3]);
    });
  });

  describe('broadcastTo', () => {
    it('should broadcast 1D to 2D', async () => {
      const arr = await NDArray.fromArray([1, 2, 3]);
      const broadcast = await broadcastTo(arr, [2, 3]);
      expect(broadcast.shape).toEqual([2, 3]);
      // Both rows should be the same
      expect(broadcast.get(0, 0)).toBe(1);
      expect(broadcast.get(1, 0)).toBe(1);
      expect(broadcast.get(0, 2)).toBe(3);
      expect(broadcast.get(1, 2)).toBe(3);
      arr.dispose();
      broadcast.dispose();
    });

    it('should throw for incompatible shapes', async () => {
      const arr = await NDArray.fromArray([1, 2, 3]);
      await expect(broadcastTo(arr, [2, 4])).rejects.toThrow();
      arr.dispose();
    });
  });

  describe('broadcastArrays', () => {
    it('should broadcast multiple arrays', async () => {
      const a = await NDArray.fromArray([1, 2, 3]);
      const b = await NDArray.fromArray([[1], [2]]);
      const [aBc, bBc] = await broadcastArrays(a, b);

      expect(aBc.shape).toEqual([2, 3]);
      expect(bBc.shape).toEqual([2, 3]);

      // a is broadcast along rows
      expect(aBc.get(0, 0)).toBe(1);
      expect(aBc.get(1, 0)).toBe(1);

      // b is broadcast along columns
      expect(bBc.get(0, 0)).toBe(1);
      expect(bBc.get(0, 2)).toBe(1);
      expect(bBc.get(1, 0)).toBe(2);

      a.dispose();
      b.dispose();
      aBc.dispose();
      bBc.dispose();
    });
  });
});

describe('Level 2: Index Functions', () => {
  describe('take', () => {
    it('should take elements along axis 0', async () => {
      const arr = await NDArray.fromArray([
        [1, 2],
        [3, 4],
        [5, 6],
      ]);
      const indices = await NDArray.fromArray([0, 2], [2], { dtype: DType.Int32 });
      const taken = await take(arr, indices, 0);
      expect(taken.shape).toEqual([2, 2]);
      expect(taken.toArray()).toEqual([1, 2, 5, 6]);
      arr.dispose();
      indices.dispose();
      taken.dispose();
    });

    it('should handle wrap mode', async () => {
      const arr = await NDArray.fromArray([10, 20, 30]);
      const indices = await NDArray.fromArray([0, 5], [2], { dtype: DType.Int32 }); // 5 is out of bounds
      const taken = await take(arr, indices, 0, CLIP_WRAP);
      expect(taken.get(1)).toBe(30); // 5 % 3 = 2 -> arr[2] = 30
      arr.dispose();
      indices.dispose();
      taken.dispose();
    });
  });

  describe('put', () => {
    it('should put values at indices', async () => {
      const arr = await NDArray.zeros([3, 3]);
      const indices = await NDArray.fromArray([0, 4, 8], [3], { dtype: DType.Int32 });
      const values = await NDArray.fromArray([1, 1, 1]);
      await put(arr, indices, values);

      // Check diagonal is set
      expect(arr.get(0, 0)).toBe(1);
      expect(arr.get(1, 1)).toBe(1);
      expect(arr.get(2, 2)).toBe(1);
      // Off-diagonal should be 0
      expect(arr.get(0, 1)).toBe(0);

      arr.dispose();
      indices.dispose();
      values.dispose();
    });
  });

  describe('nonzero', () => {
    it('should find nonzero indices', async () => {
      const arr = await NDArray.fromArray([
        [1, 0],
        [0, 2],
      ]);
      const indices = await nonzero(arr);
      expect(indices.shape[0]).toBe(2); // 2 nonzero elements
      expect(indices.shape[1]).toBe(2); // 2D array

      // First nonzero at (0, 0)
      expect(indices.get(0, 0)).toBe(0);
      expect(indices.get(0, 1)).toBe(0);

      // Second nonzero at (1, 1)
      expect(indices.get(1, 0)).toBe(1);
      expect(indices.get(1, 1)).toBe(1);

      arr.dispose();
      indices.dispose();
    });
  });

  describe('countNonzero', () => {
    it('should count nonzero elements', async () => {
      const arr = await NDArray.fromArray([1, 0, 2, 0, 3]);
      const count = await countNonzero(arr);
      expect(count).toBe(3);
      arr.dispose();
    });
  });

  describe('where', () => {
    it('should select from x or y based on condition', async () => {
      const cond = await NDArray.fromArray([1, 0, 1], [3], { dtype: DType.Bool });
      const x = await NDArray.fromArray([1, 2, 3]);
      const y = await NDArray.fromArray([10, 20, 30]);
      const result = await where(cond, x, y);

      expect(result.toArray()).toEqual([1, 20, 3]);

      cond.dispose();
      x.dispose();
      y.dispose();
      result.dispose();
    });
  });

  describe('extract', () => {
    it('should extract elements where condition is true', async () => {
      const arr = await NDArray.fromArray([
        [1, 2],
        [3, 4],
      ]);
      const cond = await NDArray.fromArray(
        [
          [1, 0],
          [0, 1],
        ],
        [2, 2],
        { dtype: DType.Bool }
      );
      const result = await extract(cond, arr);

      expect(result.toArray()).toEqual([1, 4]);

      arr.dispose();
      cond.dispose();
      result.dispose();
    });
  });

  describe('compress', () => {
    it('should compress along axis 0', async () => {
      const arr = await NDArray.fromArray([
        [1, 2],
        [3, 4],
        [5, 6],
      ]);
      const cond = await NDArray.fromArray([1, 0, 1], [3], { dtype: DType.Bool });
      const result = await compress(cond, arr, 0);

      expect(result.shape).toEqual([2, 2]);
      expect(result.toArray()).toEqual([1, 2, 5, 6]);

      arr.dispose();
      cond.dispose();
      result.dispose();
    });
  });

  describe('diagonal', () => {
    it('should get main diagonal', async () => {
      const arr = await NDArray.fromArray([
        [1, 2, 3],
        [4, 5, 6],
        [7, 8, 9],
      ]);
      const diag = await diagonal(arr);

      expect(diag.toArray()).toEqual([1, 5, 9]);

      arr.dispose();
      diag.dispose();
    });

    it('should get upper diagonal', async () => {
      const arr = await NDArray.fromArray([
        [1, 2, 3],
        [4, 5, 6],
        [7, 8, 9],
      ]);
      const diag = await diagonal(arr, 1);

      expect(diag.toArray()).toEqual([2, 6]);

      arr.dispose();
      diag.dispose();
    });

    it('should get lower diagonal', async () => {
      const arr = await NDArray.fromArray([
        [1, 2, 3],
        [4, 5, 6],
        [7, 8, 9],
      ]);
      const diag = await diagonal(arr, -1);

      expect(diag.toArray()).toEqual([4, 8]);

      arr.dispose();
      diag.dispose();
    });
  });
});

describe('Level 2: View Extensions', () => {
  describe('ascontiguousarray', () => {
    it('should return view for already contiguous array', async () => {
      const arr = await NDArray.fromArray([1, 2, 3, 4, 5, 6], [2, 3]);
      const contiguous = arr.ascontiguousarray();
      expect(contiguous.flags.c_contiguous).toBe(true);
      expect(contiguous.toArray()).toEqual([1, 2, 3, 4, 5, 6]);
      arr.dispose();
      contiguous.dispose();
    });

    it('should copy transposed array to make contiguous', async () => {
      const arr = await NDArray.fromArray([1, 2, 3, 4, 5, 6], [2, 3]);
      const transposed = arr.T;
      expect(transposed.flags.c_contiguous).toBe(false);

      const contiguous = transposed.ascontiguousarray();
      expect(contiguous.flags.c_contiguous).toBe(true);

      arr.dispose();
      transposed.dispose();
      contiguous.dispose();
    });
  });
});

/* ============ Phase 2 & 3 Missing Features ============ */

describe('Phase 2: Shape Manipulation - moveaxis', () => {
  it('should move single axis to new position', async () => {
    const arr = await NDArray.zeros([3, 4, 5]);
    const moved = arr.moveaxis(0, -1);
    expect(moved.shape).toEqual([4, 5, 3]);
    arr.dispose();
    moved.dispose();
  });

  it('should move multiple axes', async () => {
    const arr = await NDArray.zeros([3, 4, 5]);
    const moved = arr.moveaxis([0, 1], [-1, -2]);
    expect(moved.shape).toEqual([5, 4, 3]);
    arr.dispose();
    moved.dispose();
  });

  it('should be a view (share data)', async () => {
    const arr = await NDArray.fromArray([1, 2, 3, 4, 5, 6], [2, 3]);
    const moved = arr.moveaxis(0, 1);
    expect(moved.isView).toBe(true);
    expect(moved.shape).toEqual([3, 2]);
    arr.dispose();
    moved.dispose();
  });

  it('should throw on invalid axis', async () => {
    const arr = await NDArray.zeros([3, 4]);
    expect(() => arr.moveaxis(5, 0)).toThrow();
    arr.dispose();
  });

  it('should throw on mismatched source/destination lengths', async () => {
    const arr = await NDArray.zeros([3, 4, 5]);
    expect(() => arr.moveaxis([0, 1], [2])).toThrow();
    arr.dispose();
  });
});

describe('Phase 2: atleast_1d, atleast_2d, atleast_3d', () => {
  describe('atleast_1d', () => {
    it('should convert 0-d to 1-d', async () => {
      const arr = await NDArray.fromArray(5, []);
      const result = atleast_1d(arr) as NDArray;
      expect(result.shape).toEqual([1]);
      expect(result.item()).toBe(5);
      arr.dispose();
      result.dispose();
    });

    it('should keep 1-d unchanged', async () => {
      const arr = await NDArray.fromArray([1, 2, 3]);
      const result = atleast_1d(arr) as NDArray;
      expect(result.shape).toEqual([3]);
      arr.dispose();
      result.dispose();
    });

    it('should return tuple for multiple inputs', async () => {
      const a = await NDArray.fromArray([1, 2]);
      const b = await NDArray.fromArray([[1, 2], [3, 4]]);
      const results = atleast_1d(a, b) as NDArray[];
      expect(results.length).toBe(2);
      expect(results[0].shape).toEqual([2]);
      expect(results[1].shape).toEqual([2, 2]);
      a.dispose();
      b.dispose();
      results.forEach(r => r.dispose());
    });
  });

  describe('atleast_2d', () => {
    it('should convert 1-d to 2-d', async () => {
      const arr = await NDArray.fromArray([1, 2, 3]);
      const result = atleast_2d(arr) as NDArray;
      expect(result.shape).toEqual([1, 3]);
      arr.dispose();
      result.dispose();
    });

    it('should keep 2-d unchanged', async () => {
      const arr = await NDArray.fromArray([[1, 2], [3, 4]]);
      const result = atleast_2d(arr) as NDArray;
      expect(result.shape).toEqual([2, 2]);
      arr.dispose();
      result.dispose();
    });
  });

  describe('atleast_3d', () => {
    it('should convert 1-d to 3-d', async () => {
      const arr = await NDArray.fromArray([1, 2, 3]);
      const result = atleast_3d(arr) as NDArray;
      expect(result.shape).toEqual([1, 3, 1]);
      arr.dispose();
      result.dispose();
    });

    it('should convert 2-d to 3-d', async () => {
      const arr = await NDArray.fromArray([[1, 2], [3, 4]]);
      const result = atleast_3d(arr) as NDArray;
      expect(result.shape).toEqual([2, 2, 1]);
      arr.dispose();
      result.dispose();
    });

    it('should keep 3-d unchanged', async () => {
      const arr = await NDArray.zeros([2, 3, 4]);
      const result = atleast_3d(arr) as NDArray;
      expect(result.shape).toEqual([2, 3, 4]);
      arr.dispose();
      result.dispose();
    });
  });
});

describe('Phase 3: Index Generation', () => {
  describe('argwhere', () => {
    it('should find indices of nonzero elements', async () => {
      const arr = await NDArray.fromArray([[1, 0], [0, 2]]);
      const result = await argwhere(arr);
      expect(result.shape[0]).toBe(2);  // 2 nonzero elements
      expect(result.shape[1]).toBe(2);  // 2D array
      // First nonzero at (0, 0)
      expect(result.get(0, 0)).toBe(0);
      expect(result.get(0, 1)).toBe(0);
      // Second nonzero at (1, 1)
      expect(result.get(1, 0)).toBe(1);
      expect(result.get(1, 1)).toBe(1);
      arr.dispose();
      result.dispose();
    });

    it('should return empty for all-zero array', async () => {
      const arr = await NDArray.zeros([3, 3]);
      const result = await argwhere(arr);
      expect(result.shape[0]).toBe(0);
      arr.dispose();
      result.dispose();
    });
  });

  describe('indices', () => {
    it('should generate dense index grid', async () => {
      const result = await indices([2, 3]) as NDArray;
      expect(result.shape).toEqual([2, 2, 3]);
      // First slice contains row indices
      expect(result.get(0, 0, 0)).toBe(0);
      expect(result.get(0, 1, 0)).toBe(1);
      // Second slice contains column indices
      expect(result.get(1, 0, 0)).toBe(0);
      expect(result.get(1, 0, 1)).toBe(1);
      expect(result.get(1, 0, 2)).toBe(2);
      result.dispose();
    });

    it('should generate sparse index grid', async () => {
      const result = await indices([2, 3], DType.Int32, true) as NDArray[];
      expect(result.length).toBe(2);
      expect(result[0].shape).toEqual([2, 1]);
      expect(result[1].shape).toEqual([1, 3]);
      result.forEach(r => r.dispose());
    });
  });

  describe('ix_', () => {
    it('should create open mesh from sequences', async () => {
      const result = await ix_([0, 1], [2, 3, 4]);
      expect(result.length).toBe(2);
      expect(result[0].shape).toEqual([2, 1]);
      expect(result[1].shape).toEqual([1, 3]);
      expect(result[0].get(0, 0)).toBe(0);
      expect(result[0].get(1, 0)).toBe(1);
      expect(result[1].get(0, 0)).toBe(2);
      expect(result[1].get(0, 1)).toBe(3);
      expect(result[1].get(0, 2)).toBe(4);
      result.forEach(r => r.dispose());
    });

    it('should work with NDArray inputs', async () => {
      const a = await NDArray.fromArray([0, 1, 2]);
      const b = await NDArray.fromArray([3, 4]);
      const result = await ix_(a, b);
      expect(result[0].shape).toEqual([3, 1]);
      expect(result[1].shape).toEqual([1, 2]);
      a.dispose();
      b.dispose();
      result.forEach(r => r.dispose());
    });
  });

  describe('diag_indices', () => {
    it('should return diagonal index arrays', async () => {
      const [i, j] = await diag_indices(3);
      expect(i.toArray()).toEqual([0, 1, 2]);
      expect(j.toArray()).toEqual([0, 1, 2]);
      i.dispose();
      j.dispose();
    });

    it('should support higher dimensions', async () => {
      const result = await diag_indices(2, 3);
      expect(result.length).toBe(3);
      result.forEach(r => {
        expect(r.toArray()).toEqual([0, 1]);
        r.dispose();
      });
    });
  });

  describe('tril_indices', () => {
    it('should return lower triangle indices', async () => {
      const [rows, cols] = await tril_indices(3);
      expect(rows.toArray()).toEqual([0, 1, 1, 2, 2, 2]);
      expect(cols.toArray()).toEqual([0, 0, 1, 0, 1, 2]);
      rows.dispose();
      cols.dispose();
    });

    it('should support diagonal offset', async () => {
      const [rows, cols] = await tril_indices(3, 1);
      // k=1 includes first super-diagonal (lower tri + 1 diagonal above main)
      // For 3x3: [0,0], [0,1], [1,0], [1,1], [1,2], [2,0], [2,1], [2,2] = 8 elements
      expect(rows.size).toBe(8);
      rows.dispose();
      cols.dispose();
    });
  });

  describe('triu_indices', () => {
    it('should return upper triangle indices', async () => {
      const [rows, cols] = await triu_indices(3);
      expect(rows.toArray()).toEqual([0, 0, 0, 1, 1, 2]);
      expect(cols.toArray()).toEqual([0, 1, 2, 1, 2, 2]);
      rows.dispose();
      cols.dispose();
    });
  });
});

describe('Phase 3: Advanced Indexing', () => {
  describe('take_along_axis', () => {
    it('should take values along axis', async () => {
      const arr = await NDArray.fromArray([[10, 20, 30], [40, 50, 60]]);
      const idx = await NDArray.fromArray([[0, 2], [1, 0]], [2, 2], { dtype: DType.Int32 });
      const result = await take_along_axis(arr, idx, 1);
      expect(result.shape).toEqual([2, 2]);
      expect(result.toArray()).toEqual([10, 30, 50, 40]);
      arr.dispose();
      idx.dispose();
      result.dispose();
    });

    it('should work with axis 0', async () => {
      const arr = await NDArray.fromArray([[1, 2], [3, 4]]);
      const idx = await NDArray.fromArray([[1, 0], [0, 1]], [2, 2], { dtype: DType.Int32 });
      const result = await take_along_axis(arr, idx, 0);
      expect(result.toArray()).toEqual([3, 2, 1, 4]);
      arr.dispose();
      idx.dispose();
      result.dispose();
    });
  });

  describe('put_along_axis', () => {
    it('should put values along axis', async () => {
      const arr = await NDArray.zeros([2, 3]);
      const idx = await NDArray.fromArray([[0, 2], [1, 0]], [2, 2], { dtype: DType.Int32 });
      const vals = await NDArray.fromArray([[1, 2], [3, 4]]);
      await put_along_axis(arr, idx, vals, 1);
      // Row 0: positions 0 and 2 get 1 and 2
      expect(arr.get(0, 0)).toBe(1);
      expect(arr.get(0, 2)).toBe(2);
      // Row 1: positions 1 and 0 get 3 and 4
      expect(arr.get(1, 1)).toBe(3);
      expect(arr.get(1, 0)).toBe(4);
      arr.dispose();
      idx.dispose();
      vals.dispose();
    });
  });

  describe('putmask', () => {
    it('should put values where mask is true', async () => {
      const arr = await NDArray.fromArray([1, 2, 3, 4]);
      const mask = await NDArray.fromArray([1, 0, 1, 0], [4], { dtype: DType.Bool });
      const vals = await NDArray.fromArray([10, 20]);
      await putmask(arr, mask, vals);
      expect(arr.toArray()).toEqual([10, 2, 20, 4]);
      arr.dispose();
      mask.dispose();
      vals.dispose();
    });

    it('should cycle values if needed', async () => {
      const arr = await NDArray.fromArray([1, 2, 3, 4, 5]);
      const mask = await NDArray.fromArray([1, 1, 1, 1, 1], [5], { dtype: DType.Bool });
      const vals = await NDArray.fromArray([10, 20]);
      await putmask(arr, mask, vals);
      expect(arr.toArray()).toEqual([10, 20, 10, 20, 10]);
      arr.dispose();
      mask.dispose();
      vals.dispose();
    });
  });

  describe('select', () => {
    it('should select based on conditions', async () => {
      const x = await NDArray.fromArray([0, 1, 2, 3, 4]);
      const cond1 = await NDArray.fromArray([0, 0, 1, 1, 0], [5], { dtype: DType.Bool });
      const cond2 = await NDArray.fromArray([1, 0, 0, 0, 1], [5], { dtype: DType.Bool });
      const choice1 = await NDArray.fromArray([10, 10, 10, 10, 10]);
      const choice2 = await NDArray.fromArray([20, 20, 20, 20, 20]);
      const result = await select([cond1, cond2], [choice1, choice2], 0);
      expect(result.toArray()).toEqual([20, 0, 10, 10, 20]);
      x.dispose();
      cond1.dispose();
      cond2.dispose();
      choice1.dispose();
      choice2.dispose();
      result.dispose();
    });

    it('should use default for no match', async () => {
      const cond = await NDArray.fromArray([0, 0, 0], [3], { dtype: DType.Bool });
      const choice = await NDArray.fromArray([10, 20, 30]);
      const result = await select([cond], [choice], -1);
      expect(result.toArray()).toEqual([-1, -1, -1]);
      cond.dispose();
      choice.dispose();
      result.dispose();
    });
  });
});
