/**
 * Tests for Phase 5: Array Manipulation Functions
 */

import { describe, it, expect } from 'vitest';
import {
  NDArray,
  DType,
  // Joining
  concatenate,
  stack,
  vstack,
  hstack,
  dstack,
  column_stack,
  block,
  append,
  // Splitting
  split,
  array_split,
  vsplit,
  hsplit,
  dsplit,
  unstack,
  // Tiling
  tile,
  repeat,
  pad,
  // Rearranging
  flip,
  fliplr,
  flipud,
  roll,
  rot90,
  resize,
  trim_zeros,
  // Insert/Delete
  insert,
  deleteArr,
  // Copying
  copyto,
  asarray,
  slice,
} from '../../dist/numjs.mjs';

describe('Joining Functions', () => {
  describe('concatenate', () => {
    it('should concatenate 1D arrays along axis 0', async () => {
      const a = await NDArray.fromArray([1, 2, 3]);
      const b = await NDArray.fromArray([4, 5, 6]);
      const result = concatenate([a, b], 0);

      expect(result.shape).toEqual([6]);
      expect(result.toArray()).toEqual([1, 2, 3, 4, 5, 6]);

      a.dispose();
      b.dispose();
      result.dispose();
    });

    it('should concatenate 2D arrays along axis 0', async () => {
      const a = await NDArray.fromArray([[1, 2], [3, 4]], [2, 2]);
      const b = await NDArray.fromArray([[5, 6]], [1, 2]);
      const result = concatenate([a, b], 0);

      expect(result.shape).toEqual([3, 2]);
      expect(result.toArray()).toEqual([1, 2, 3, 4, 5, 6]);

      a.dispose();
      b.dispose();
      result.dispose();
    });

    it('should concatenate 2D arrays along axis 1', async () => {
      const a = await NDArray.fromArray([[1, 2], [3, 4]], [2, 2]);
      const b = await NDArray.fromArray([[5], [6]], [2, 1]);
      const result = concatenate([a, b], 1);

      expect(result.shape).toEqual([2, 3]);
      expect(result.toArray()).toEqual([1, 2, 5, 3, 4, 6]);

      a.dispose();
      b.dispose();
      result.dispose();
    });

    it('should support negative axis', async () => {
      const a = await NDArray.fromArray([[1, 2]], [1, 2]);
      const b = await NDArray.fromArray([[3, 4]], [1, 2]);
      const result = concatenate([a, b], -2); // Same as axis 0

      expect(result.shape).toEqual([2, 2]);

      a.dispose();
      b.dispose();
      result.dispose();
    });

    it('should throw on shape mismatch', async () => {
      const a = await NDArray.fromArray([[1, 2], [3, 4]], [2, 2]);
      const b = await NDArray.fromArray([[5, 6, 7]], [1, 3]);

      expect(() => concatenate([a, b], 0)).toThrow();

      a.dispose();
      b.dispose();
    });
  });

  describe('stack', () => {
    it('should stack arrays along a new axis', async () => {
      const a = await NDArray.fromArray([1, 2, 3]);
      const b = await NDArray.fromArray([4, 5, 6]);
      const result = stack([a, b], 0);

      expect(result.shape).toEqual([2, 3]);
      expect(result.toArray()).toEqual([1, 2, 3, 4, 5, 6]);

      a.dispose();
      b.dispose();
      result.dispose();
    });

    it('should stack along axis 1', async () => {
      const a = await NDArray.fromArray([1, 2, 3]);
      const b = await NDArray.fromArray([4, 5, 6]);
      const result = stack([a, b], 1);

      expect(result.shape).toEqual([3, 2]);
      expect(result.toArray()).toEqual([1, 4, 2, 5, 3, 6]);

      a.dispose();
      b.dispose();
      result.dispose();
    });

    it('should throw on shape mismatch', async () => {
      const a = await NDArray.fromArray([1, 2, 3]);
      const b = await NDArray.fromArray([4, 5]);

      expect(() => stack([a, b], 0)).toThrow();

      a.dispose();
      b.dispose();
    });
  });

  describe('vstack', () => {
    it('should stack 1D arrays as rows', async () => {
      const a = await NDArray.fromArray([1, 2, 3]);
      const b = await NDArray.fromArray([4, 5, 6]);
      const result = vstack([a, b]);

      expect(result.shape).toEqual([2, 3]);
      expect(result.toArray()).toEqual([1, 2, 3, 4, 5, 6]);

      a.dispose();
      b.dispose();
      result.dispose();
    });

    it('should stack 2D arrays vertically', async () => {
      const a = await NDArray.fromArray([[1, 2], [3, 4]], [2, 2]);
      const b = await NDArray.fromArray([[5, 6]], [1, 2]);
      const result = vstack([a, b]);

      expect(result.shape).toEqual([3, 2]);

      a.dispose();
      b.dispose();
      result.dispose();
    });
  });

  describe('hstack', () => {
    it('should concatenate 1D arrays', async () => {
      const a = await NDArray.fromArray([1, 2, 3]);
      const b = await NDArray.fromArray([4, 5, 6]);
      const result = hstack([a, b]);

      expect(result.shape).toEqual([6]);
      expect(result.toArray()).toEqual([1, 2, 3, 4, 5, 6]);

      a.dispose();
      b.dispose();
      result.dispose();
    });

    it('should stack 2D arrays horizontally', async () => {
      const a = await NDArray.fromArray([[1, 2], [3, 4]], [2, 2]);
      const b = await NDArray.fromArray([[5], [6]], [2, 1]);
      const result = hstack([a, b]);

      expect(result.shape).toEqual([2, 3]);

      a.dispose();
      b.dispose();
      result.dispose();
    });
  });

  describe('dstack', () => {
    it('should stack along the third axis', async () => {
      const a = await NDArray.fromArray([[1, 2], [3, 4]], [2, 2]);
      const b = await NDArray.fromArray([[5, 6], [7, 8]], [2, 2]);
      const result = dstack([a, b]);

      expect(result.shape).toEqual([2, 2, 2]);

      a.dispose();
      b.dispose();
      result.dispose();
    });
  });

  describe('column_stack', () => {
    it('should stack 1D arrays as columns', async () => {
      const a = await NDArray.fromArray([1, 2, 3]);
      const b = await NDArray.fromArray([4, 5, 6]);
      const result = column_stack([a, b]);

      expect(result.shape).toEqual([3, 2]);
      expect(result.toArray()).toEqual([1, 4, 2, 5, 3, 6]);

      a.dispose();
      b.dispose();
      result.dispose();
    });
  });

  describe('append', () => {
    it('should append values to flattened array', async () => {
      const a = await NDArray.fromArray([1, 2, 3]);
      const b = await NDArray.fromArray([4, 5, 6]);
      const result = append(a, b);

      expect(result.shape).toEqual([6]);
      expect(result.toArray()).toEqual([1, 2, 3, 4, 5, 6]);

      a.dispose();
      b.dispose();
      result.dispose();
    });

    it('should append along specified axis', async () => {
      const a = await NDArray.fromArray([[1, 2]], [1, 2]);
      const b = await NDArray.fromArray([[3, 4]], [1, 2]);
      const result = append(a, b, 0);

      expect(result.shape).toEqual([2, 2]);

      a.dispose();
      b.dispose();
      result.dispose();
    });
  });
});

describe('Splitting Functions', () => {
  describe('split', () => {
    it('should split into equal sections', async () => {
      const arr = await NDArray.fromArray([1, 2, 3, 4, 5, 6]);
      const result = split(arr, 3);

      expect(result.length).toBe(3);
      expect(result[0].shape).toEqual([2]);
      expect(result[0].toArray()).toEqual([1, 2]);
      expect(result[1].toArray()).toEqual([3, 4]);
      expect(result[2].toArray()).toEqual([5, 6]);

      arr.dispose();
      result.forEach(r => r.dispose());
    });

    it('should throw on unequal division', async () => {
      const arr = await NDArray.fromArray([1, 2, 3, 4, 5]);

      expect(() => split(arr, 3)).toThrow();

      arr.dispose();
    });

    it('should split at indices', async () => {
      const arr = await NDArray.fromArray([1, 2, 3, 4, 5, 6]);
      const result = split(arr, [2, 4]);

      expect(result.length).toBe(3);
      expect(result[0].toArray()).toEqual([1, 2]);
      expect(result[1].toArray()).toEqual([3, 4]);
      expect(result[2].toArray()).toEqual([5, 6]);

      arr.dispose();
      result.forEach(r => r.dispose());
    });
  });

  describe('array_split', () => {
    it('should allow unequal division', async () => {
      const arr = await NDArray.fromArray([1, 2, 3, 4, 5]);
      const result = array_split(arr, 3);

      expect(result.length).toBe(3);
      // First section gets extra element
      expect(result[0].toArray()).toEqual([1, 2]);
      expect(result[1].toArray()).toEqual([3, 4]);
      expect(result[2].toArray()).toEqual([5]);

      arr.dispose();
      result.forEach(r => r.dispose());
    });
  });

  describe('vsplit', () => {
    it('should split vertically', async () => {
      const arr = await NDArray.fromArray([[1, 2], [3, 4], [5, 6], [7, 8]], [4, 2]);
      const result = vsplit(arr, 2);

      expect(result.length).toBe(2);
      expect(result[0].shape).toEqual([2, 2]);
      expect(result[1].shape).toEqual([2, 2]);

      arr.dispose();
      result.forEach(r => r.dispose());
    });
  });

  describe('hsplit', () => {
    it('should split horizontally', async () => {
      const arr = await NDArray.fromArray([[1, 2, 3, 4], [5, 6, 7, 8]], [2, 4]);
      const result = hsplit(arr, 2);

      expect(result.length).toBe(2);
      expect(result[0].shape).toEqual([2, 2]);
      expect(result[1].shape).toEqual([2, 2]);

      arr.dispose();
      result.forEach(r => r.dispose());
    });
  });

  describe('unstack', () => {
    it('should unstack along axis 0', async () => {
      const arr = await NDArray.fromArray([[1, 2], [3, 4]], [2, 2]);
      const result = unstack(arr, 0);

      expect(result.length).toBe(2);
      expect(result[0].shape).toEqual([2]);
      expect(result[0].toArray()).toEqual([1, 2]);
      expect(result[1].toArray()).toEqual([3, 4]);

      arr.dispose();
      result.forEach(r => r.dispose());
    });
  });
});

describe('Rearranging Functions', () => {
  describe('flip', () => {
    it('should flip along axis 0', async () => {
      const arr = await NDArray.fromArray([[1, 2], [3, 4]], [2, 2]);
      const result = flip(arr, 0);

      expect(result.shape).toEqual([2, 2]);
      expect(result.toArray()).toEqual([3, 4, 1, 2]);

      arr.dispose();
      result.dispose();
    });

    it('should flip along axis 1', async () => {
      const arr = await NDArray.fromArray([[1, 2], [3, 4]], [2, 2]);
      const result = flip(arr, 1);

      expect(result.shape).toEqual([2, 2]);
      expect(result.toArray()).toEqual([2, 1, 4, 3]);

      arr.dispose();
      result.dispose();
    });

    it('should flip all axes when axis is undefined', async () => {
      const arr = await NDArray.fromArray([[1, 2], [3, 4]], [2, 2]);
      const result = flip(arr);

      expect(result.toArray()).toEqual([4, 3, 2, 1]);

      arr.dispose();
      result.dispose();
    });
  });

  describe('fliplr', () => {
    it('should flip left-right', async () => {
      const arr = await NDArray.fromArray([[1, 2], [3, 4]], [2, 2]);
      const result = fliplr(arr);

      expect(result.toArray()).toEqual([2, 1, 4, 3]);

      arr.dispose();
      result.dispose();
    });
  });

  describe('flipud', () => {
    it('should flip up-down', async () => {
      const arr = await NDArray.fromArray([[1, 2], [3, 4]], [2, 2]);
      const result = flipud(arr);

      expect(result.toArray()).toEqual([3, 4, 1, 2]);

      arr.dispose();
      result.dispose();
    });
  });

  describe('roll', () => {
    it('should roll elements', async () => {
      const arr = await NDArray.fromArray([1, 2, 3, 4, 5]);
      const result = roll(arr, 2);

      expect(result.toArray()).toEqual([4, 5, 1, 2, 3]);

      arr.dispose();
      result.dispose();
    });

    it('should roll with negative shift', async () => {
      const arr = await NDArray.fromArray([1, 2, 3, 4, 5]);
      const result = roll(arr, -2);

      expect(result.toArray()).toEqual([3, 4, 5, 1, 2]);

      arr.dispose();
      result.dispose();
    });

    it('should roll along axis', async () => {
      const arr = await NDArray.fromArray([[1, 2, 3], [4, 5, 6]], [2, 3]);
      const result = roll(arr, 1, 1);

      expect(result.toArray()).toEqual([3, 1, 2, 6, 4, 5]);

      arr.dispose();
      result.dispose();
    });
  });

  describe('rot90', () => {
    it('should rotate 90 degrees', async () => {
      const arr = await NDArray.fromArray([[1, 2], [3, 4]], [2, 2]);
      const result = rot90(arr, 1);

      expect(result.shape).toEqual([2, 2]);
      expect(result.toArray()).toEqual([2, 4, 1, 3]);

      arr.dispose();
      result.dispose();
    });

    it('should rotate 180 degrees', async () => {
      const arr = await NDArray.fromArray([[1, 2], [3, 4]], [2, 2]);
      const result = rot90(arr, 2);

      expect(result.toArray()).toEqual([4, 3, 2, 1]);

      arr.dispose();
      result.dispose();
    });
  });

  describe('resize', () => {
    it('should resize with repetition', async () => {
      const arr = await NDArray.fromArray([1, 2, 3]);
      const result = await resize(arr, [5]);

      expect(result.toArray()).toEqual([1, 2, 3, 1, 2]);

      arr.dispose();
      result.dispose();
    });

    it('should resize with truncation', async () => {
      const arr = await NDArray.fromArray([1, 2, 3, 4, 5]);
      const result = await resize(arr, [3]);

      expect(result.toArray()).toEqual([1, 2, 3]);

      arr.dispose();
      result.dispose();
    });
  });

  describe('trim_zeros', () => {
    it('should trim leading and trailing zeros', async () => {
      const arr = await NDArray.fromArray([0, 0, 1, 2, 3, 0, 0]);
      const result = trim_zeros(arr);

      expect(result.toArray()).toEqual([1, 2, 3]);

      arr.dispose();
      result.dispose();
    });

    it('should trim only front', async () => {
      const arr = await NDArray.fromArray([0, 0, 1, 2, 3, 0, 0]);
      const result = trim_zeros(arr, 'f');

      expect(result.toArray()).toEqual([1, 2, 3, 0, 0]);

      arr.dispose();
      result.dispose();
    });

    it('should trim only back', async () => {
      const arr = await NDArray.fromArray([0, 0, 1, 2, 3, 0, 0]);
      const result = trim_zeros(arr, 'b');

      expect(result.toArray()).toEqual([0, 0, 1, 2, 3]);

      arr.dispose();
      result.dispose();
    });
  });
});

describe('Tiling Functions', () => {
  describe('tile', () => {
    it('should tile 1D array', async () => {
      const arr = await NDArray.fromArray([1, 2]);
      const result = await tile(arr, 3);

      expect(result.toArray()).toEqual([1, 2, 1, 2, 1, 2]);

      arr.dispose();
      result.dispose();
    });

    it('should tile 2D array', async () => {
      const arr = await NDArray.fromArray([[1, 2], [3, 4]], [2, 2]);
      const result = await tile(arr, [2, 1]);

      expect(result.shape).toEqual([4, 2]);
      expect(result.toArray()).toEqual([1, 2, 3, 4, 1, 2, 3, 4]);

      arr.dispose();
      result.dispose();
    });
  });

  describe('repeat', () => {
    it('should repeat elements', async () => {
      const arr = await NDArray.fromArray([1, 2, 3]);
      const result = await repeat(arr, 2);

      expect(result.toArray()).toEqual([1, 1, 2, 2, 3, 3]);

      arr.dispose();
      result.dispose();
    });

    it('should repeat along axis', async () => {
      const arr = await NDArray.fromArray([[1, 2], [3, 4]], [2, 2]);
      const result = await repeat(arr, 2, 0);

      expect(result.shape).toEqual([4, 2]);
      expect(result.toArray()).toEqual([1, 2, 1, 2, 3, 4, 3, 4]);

      arr.dispose();
      result.dispose();
    });
  });

  describe('pad', () => {
    it('should pad with constant', async () => {
      const arr = await NDArray.fromArray([1, 2, 3]);
      const result = await pad(arr, 2, 'constant', 0);

      expect(result.toArray()).toEqual([0, 0, 1, 2, 3, 0, 0]);

      arr.dispose();
      result.dispose();
    });

    it('should pad 2D array', async () => {
      const arr = await NDArray.fromArray([[1, 2], [3, 4]], [2, 2]);
      const result = await pad(arr, [[1, 1], [1, 1]], 'constant', 0);

      expect(result.shape).toEqual([4, 4]);

      arr.dispose();
      result.dispose();
    });
  });
});

describe('Insert/Delete Functions', () => {
  describe('insert', () => {
    it('should insert values', async () => {
      const arr = await NDArray.fromArray([1, 2, 3, 4, 5]);
      const values = await NDArray.fromArray([10]);
      const result = insert(arr, 2, values);

      expect(result.toArray()).toEqual([1, 2, 10, 3, 4, 5]);

      arr.dispose();
      values.dispose();
      result.dispose();
    });
  });

  describe('deleteArr', () => {
    it('should delete elements', async () => {
      const arr = await NDArray.fromArray([1, 2, 3, 4, 5]);
      const result = await deleteArr(arr, [1, 3]);

      expect(result.toArray()).toEqual([1, 3, 5]);

      arr.dispose();
      result.dispose();
    });
  });
});

describe('Copying Functions', () => {
  describe('copyto', () => {
    it('should copy data', async () => {
      const dst = await NDArray.fromArray([0, 0, 0]);
      const src = await NDArray.fromArray([1, 2, 3]);

      copyto(dst, src);

      expect(dst.toArray()).toEqual([1, 2, 3]);

      dst.dispose();
      src.dispose();
    });

    it('should copy with mask', async () => {
      const dst = await NDArray.fromArray([0, 0, 0]);
      const src = await NDArray.fromArray([1, 2, 3]);
      const where = await NDArray.fromArray([1, 0, 1]);

      copyto(dst, src, 'same_kind', where);

      expect(dst.toArray()).toEqual([1, 0, 3]);

      dst.dispose();
      src.dispose();
      where.dispose();
    });
  });

  describe('asarray', () => {
    it('should return NDArray as-is', async () => {
      const arr = await NDArray.fromArray([1, 2, 3]);
      const result = await asarray(arr);

      expect(result).toBe(arr);

      arr.dispose();
    });

    it('should convert JavaScript array', async () => {
      const result = await asarray([1, 2, 3]);

      expect(result.toArray()).toEqual([1, 2, 3]);

      result.dispose();
    });
  });
});
