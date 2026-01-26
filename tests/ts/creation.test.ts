/**
 * Array Creation Function Tests
 *
 * Tests for array creation functions: empty, full, linspace, eye, etc.
 */

import { describe, it, expect } from 'vitest';
import { NDArray, DType } from '../../dist/numjs.mjs';

describe('Array Creation Functions', () => {
  describe('NDArray.empty()', () => {
    it('should create array with specified shape', async () => {
      const arr = await NDArray.empty([3, 4]);
      expect(arr.shape).toEqual([3, 4]);
      expect(arr.size).toBe(12);
      arr.dispose();
    });

    it('should create array with specified dtype', async () => {
      const arr = await NDArray.empty([5], { dtype: DType.Int32 });
      expect(arr.dtype).toBe(DType.Int32);
      arr.dispose();
    });

    it('should create 1D array', async () => {
      const arr = await NDArray.empty([10]);
      expect(arr.shape).toEqual([10]);
      expect(arr.ndim).toBe(1);
      arr.dispose();
    });
  });

  describe('NDArray.full()', () => {
    it('should create array filled with specified value', async () => {
      const arr = await NDArray.full([2, 3], 7);
      expect(arr.shape).toEqual([2, 3]);
      for (let i = 0; i < arr.size; i++) {
        expect(arr.getFlat(i)).toBe(7);
      }
      arr.dispose();
    });

    it('should work with float values', async () => {
      const arr = await NDArray.full([3], 3.14);
      expect(arr.getFlat(0)).toBeCloseTo(3.14, 10);
      expect(arr.getFlat(1)).toBeCloseTo(3.14, 10);
      expect(arr.getFlat(2)).toBeCloseTo(3.14, 10);
      arr.dispose();
    });

    it('should work with negative values', async () => {
      const arr = await NDArray.full([2], -5);
      expect(arr.getFlat(0)).toBe(-5);
      expect(arr.getFlat(1)).toBe(-5);
      arr.dispose();
    });

    it('should respect dtype', async () => {
      const arr = await NDArray.full([3], 3.7, { dtype: DType.Int32 });
      expect(arr.dtype).toBe(DType.Int32);
      expect(arr.getFlat(0)).toBe(3); // Truncated
      arr.dispose();
    });
  });

  describe('NDArray.linspace()', () => {
    it('should create evenly spaced values', async () => {
      const arr = await NDArray.linspace(0, 10, 5);
      expect(arr.size).toBe(5);
      expect(arr.getFlat(0)).toBeCloseTo(0, 10);
      expect(arr.getFlat(1)).toBeCloseTo(2.5, 10);
      expect(arr.getFlat(2)).toBeCloseTo(5, 10);
      expect(arr.getFlat(3)).toBeCloseTo(7.5, 10);
      expect(arr.getFlat(4)).toBeCloseTo(10, 10);
      arr.dispose();
    });

    it('should default to 50 samples', async () => {
      const arr = await NDArray.linspace(0, 1);
      expect(arr.size).toBe(50);
      arr.dispose();
    });

    it('should exclude endpoint when specified', async () => {
      const arr = await NDArray.linspace(0, 10, 5, false);
      expect(arr.size).toBe(5);
      expect(arr.getFlat(0)).toBeCloseTo(0, 10);
      expect(arr.getFlat(4)).toBeCloseTo(8, 10); // Not 10
      arr.dispose();
    });

    it('should handle negative ranges', async () => {
      const arr = await NDArray.linspace(-5, 5, 11);
      expect(arr.getFlat(0)).toBeCloseTo(-5, 10);
      expect(arr.getFlat(5)).toBeCloseTo(0, 10);
      expect(arr.getFlat(10)).toBeCloseTo(5, 10);
      arr.dispose();
    });

    it('should handle single point', async () => {
      const arr = await NDArray.linspace(5, 5, 1);
      expect(arr.size).toBe(1);
      expect(arr.getFlat(0)).toBeCloseTo(5, 10);
      arr.dispose();
    });
  });

  describe('NDArray.logspace()', () => {
    it('should create logarithmically spaced values', async () => {
      const arr = await NDArray.logspace(0, 2, 3);
      expect(arr.size).toBe(3);
      expect(arr.getFlat(0)).toBeCloseTo(1, 10); // 10^0
      expect(arr.getFlat(1)).toBeCloseTo(10, 10); // 10^1
      expect(arr.getFlat(2)).toBeCloseTo(100, 10); // 10^2
      arr.dispose();
    });

    it('should handle custom base', async () => {
      const arr = await NDArray.logspace(0, 3, 4, true, 2);
      expect(arr.getFlat(0)).toBeCloseTo(1, 10); // 2^0
      expect(arr.getFlat(1)).toBeCloseTo(2, 10); // 2^1
      expect(arr.getFlat(2)).toBeCloseTo(4, 10); // 2^2
      expect(arr.getFlat(3)).toBeCloseTo(8, 10); // 2^3
      arr.dispose();
    });

    it('should default to 50 samples', async () => {
      const arr = await NDArray.logspace(0, 1);
      expect(arr.size).toBe(50);
      arr.dispose();
    });
  });

  describe('NDArray.geomspace()', () => {
    it('should create geometrically spaced values', async () => {
      const arr = await NDArray.geomspace(1, 1000, 4);
      expect(arr.size).toBe(4);
      expect(arr.getFlat(0)).toBeCloseTo(1, 10);
      expect(arr.getFlat(1)).toBeCloseTo(10, 10);
      expect(arr.getFlat(2)).toBeCloseTo(100, 10);
      expect(arr.getFlat(3)).toBeCloseTo(1000, 10);
      arr.dispose();
    });

    it('should handle non-power-of-10 endpoints', async () => {
      const arr = await NDArray.geomspace(1, 8, 4);
      expect(arr.getFlat(0)).toBeCloseTo(1, 10);
      expect(arr.getFlat(1)).toBeCloseTo(2, 10);
      expect(arr.getFlat(2)).toBeCloseTo(4, 10);
      expect(arr.getFlat(3)).toBeCloseTo(8, 10);
      arr.dispose();
    });
  });

  describe('NDArray.eye()', () => {
    it('should create identity-like matrix', async () => {
      const arr = await NDArray.eye(3);
      expect(arr.shape).toEqual([3, 3]);
      expect(arr.get(0, 0)).toBe(1);
      expect(arr.get(1, 1)).toBe(1);
      expect(arr.get(2, 2)).toBe(1);
      expect(arr.get(0, 1)).toBe(0);
      expect(arr.get(1, 0)).toBe(0);
      arr.dispose();
    });

    it('should create rectangular matrix', async () => {
      const arr = await NDArray.eye(2, 4);
      expect(arr.shape).toEqual([2, 4]);
      expect(arr.get(0, 0)).toBe(1);
      expect(arr.get(1, 1)).toBe(1);
      expect(arr.get(0, 2)).toBe(0);
      expect(arr.get(0, 3)).toBe(0);
      arr.dispose();
    });

    it('should handle positive diagonal offset', async () => {
      const arr = await NDArray.eye(3, 4, 1);
      expect(arr.get(0, 0)).toBe(0);
      expect(arr.get(0, 1)).toBe(1);
      expect(arr.get(1, 2)).toBe(1);
      expect(arr.get(2, 3)).toBe(1);
      arr.dispose();
    });

    it('should handle negative diagonal offset', async () => {
      const arr = await NDArray.eye(4, 3, -1);
      expect(arr.get(0, 0)).toBe(0);
      expect(arr.get(1, 0)).toBe(1);
      expect(arr.get(2, 1)).toBe(1);
      expect(arr.get(3, 2)).toBe(1);
      arr.dispose();
    });
  });

  describe('NDArray.identity()', () => {
    it('should create square identity matrix', async () => {
      const arr = await NDArray.identity(4);
      expect(arr.shape).toEqual([4, 4]);
      for (let i = 0; i < 4; i++) {
        for (let j = 0; j < 4; j++) {
          expect(arr.get(i, j)).toBe(i === j ? 1 : 0);
        }
      }
      arr.dispose();
    });

    it('should support different dtypes', async () => {
      const arr = await NDArray.identity(3, { dtype: DType.Float32 });
      expect(arr.dtype).toBe(DType.Float32);
      expect(arr.get(0, 0)).toBeCloseTo(1, 5);
      arr.dispose();
    });
  });

  describe('NDArray.diag()', () => {
    it('should create diagonal matrix from 1D array', async () => {
      const v = await NDArray.fromArray([1, 2, 3]);
      const arr = await NDArray.diag(v);
      expect(arr.shape).toEqual([3, 3]);
      expect(arr.get(0, 0)).toBe(1);
      expect(arr.get(1, 1)).toBe(2);
      expect(arr.get(2, 2)).toBe(3);
      expect(arr.get(0, 1)).toBe(0);
      v.dispose();
      arr.dispose();
    });

    it('should handle diagonal offset', async () => {
      const v = await NDArray.fromArray([1, 2]);
      const arr = await NDArray.diag(v, 1);
      expect(arr.shape).toEqual([3, 3]);
      expect(arr.get(0, 1)).toBe(1);
      expect(arr.get(1, 2)).toBe(2);
      expect(arr.get(0, 0)).toBe(0);
      v.dispose();
      arr.dispose();
    });

    it('should extract diagonal from 2D array', async () => {
      const mat = await NDArray.fromArray([[1, 2, 3], [4, 5, 6], [7, 8, 9]]);
      const diag = await NDArray.diag(mat);
      expect(diag.shape).toEqual([3]);
      expect(diag.get(0)).toBe(1);
      expect(diag.get(1)).toBe(5);
      expect(diag.get(2)).toBe(9);
      mat.dispose();
      diag.dispose();
    });
  });

  describe('NDArray.tri()', () => {
    it('should create lower triangular matrix of ones', async () => {
      const arr = await NDArray.tri(3);
      expect(arr.shape).toEqual([3, 3]);
      expect(arr.get(0, 0)).toBe(1);
      expect(arr.get(1, 0)).toBe(1);
      expect(arr.get(1, 1)).toBe(1);
      expect(arr.get(2, 0)).toBe(1);
      expect(arr.get(2, 1)).toBe(1);
      expect(arr.get(2, 2)).toBe(1);
      expect(arr.get(0, 1)).toBe(0);
      expect(arr.get(0, 2)).toBe(0);
      arr.dispose();
    });

    it('should handle rectangular shape', async () => {
      const arr = await NDArray.tri(3, 4);
      expect(arr.shape).toEqual([3, 4]);
      arr.dispose();
    });

    it('should handle diagonal offset', async () => {
      const arr = await NDArray.tri(3, 3, 1);
      expect(arr.get(0, 0)).toBe(1);
      expect(arr.get(0, 1)).toBe(1);
      expect(arr.get(0, 2)).toBe(0);
      arr.dispose();
    });
  });

  describe('NDArray.tril()', () => {
    it('should extract lower triangle', async () => {
      const mat = await NDArray.fromArray([[1, 2, 3], [4, 5, 6], [7, 8, 9]]);
      const lower = await NDArray.tril(mat);
      expect(lower.get(0, 0)).toBe(1);
      expect(lower.get(0, 1)).toBe(0);
      expect(lower.get(0, 2)).toBe(0);
      expect(lower.get(1, 0)).toBe(4);
      expect(lower.get(1, 1)).toBe(5);
      expect(lower.get(1, 2)).toBe(0);
      expect(lower.get(2, 0)).toBe(7);
      expect(lower.get(2, 1)).toBe(8);
      expect(lower.get(2, 2)).toBe(9);
      mat.dispose();
      lower.dispose();
    });

    it('should handle diagonal offset', async () => {
      const mat = await NDArray.fromArray([[1, 2, 3], [4, 5, 6], [7, 8, 9]]);
      const lower = await NDArray.tril(mat, 1);
      expect(lower.get(0, 1)).toBe(2);
      expect(lower.get(0, 2)).toBe(0);
      mat.dispose();
      lower.dispose();
    });
  });

  describe('NDArray.triu()', () => {
    it('should extract upper triangle', async () => {
      const mat = await NDArray.fromArray([[1, 2, 3], [4, 5, 6], [7, 8, 9]]);
      const upper = await NDArray.triu(mat);
      expect(upper.get(0, 0)).toBe(1);
      expect(upper.get(0, 1)).toBe(2);
      expect(upper.get(0, 2)).toBe(3);
      expect(upper.get(1, 0)).toBe(0);
      expect(upper.get(1, 1)).toBe(5);
      expect(upper.get(1, 2)).toBe(6);
      expect(upper.get(2, 0)).toBe(0);
      expect(upper.get(2, 1)).toBe(0);
      expect(upper.get(2, 2)).toBe(9);
      mat.dispose();
      upper.dispose();
    });
  });

  describe('Like functions', () => {
    describe('NDArray.zerosLike()', () => {
      it('should create zeros with same shape and dtype', async () => {
        const template = await NDArray.fromArray([[1, 2], [3, 4]], undefined, { dtype: DType.Int32 });
        const arr = await NDArray.zerosLike(template);
        expect(arr.shape).toEqual([2, 2]);
        expect(arr.dtype).toBe(DType.Int32);
        for (let i = 0; i < arr.size; i++) {
          expect(arr.getFlat(i)).toBe(0);
        }
        template.dispose();
        arr.dispose();
      });
    });

    describe('NDArray.onesLike()', () => {
      it('should create ones with same shape and dtype', async () => {
        const template = await NDArray.fromArray([1, 2, 3], undefined, { dtype: DType.Float32 });
        const arr = await NDArray.onesLike(template);
        expect(arr.shape).toEqual([3]);
        expect(arr.dtype).toBe(DType.Float32);
        for (let i = 0; i < arr.size; i++) {
          expect(arr.getFlat(i)).toBeCloseTo(1, 5);
        }
        template.dispose();
        arr.dispose();
      });
    });

    describe('NDArray.emptyLike()', () => {
      it('should create array with same shape and dtype', async () => {
        const template = await NDArray.fromArray([[1, 2, 3], [4, 5, 6]]);
        const arr = await NDArray.emptyLike(template);
        expect(arr.shape).toEqual([2, 3]);
        expect(arr.dtype).toBe(template.dtype);
        template.dispose();
        arr.dispose();
      });
    });

    describe('NDArray.fullLike()', () => {
      it('should create filled array with same shape and dtype', async () => {
        const template = await NDArray.fromArray([1, 2, 3, 4]);
        const arr = await NDArray.fullLike(template, 42);
        expect(arr.shape).toEqual([4]);
        for (let i = 0; i < arr.size; i++) {
          expect(arr.getFlat(i)).toBe(42);
        }
        template.dispose();
        arr.dispose();
      });
    });
  });

  describe('Copy and Type Conversion', () => {
    describe('NDArray.copy()', () => {
      it('should create independent copy', async () => {
        const original = await NDArray.fromArray([1, 2, 3]);
        const copy = await original.copy();

        expect(copy.shape).toEqual(original.shape);
        expect(copy.dtype).toBe(original.dtype);

        // Modify copy
        copy.set(99, 0);

        // Original should be unchanged
        expect(original.get(0)).toBe(1);
        expect(copy.get(0)).toBe(99);

        original.dispose();
        copy.dispose();
      });
    });

    describe('NDArray.astype()', () => {
      it('should convert dtype', async () => {
        const f64 = await NDArray.fromArray([1.5, 2.7, 3.9]);
        const i32 = await f64.astype(DType.Int32);

        expect(i32.dtype).toBe(DType.Int32);
        expect(i32.get(0)).toBe(1);
        expect(i32.get(1)).toBe(2);
        expect(i32.get(2)).toBe(3);

        f64.dispose();
        i32.dispose();
      });

      it('should preserve values when converting between float types', async () => {
        const f64 = await NDArray.fromArray([1.5, 2.5, 3.5]);
        const f32 = await f64.astype(DType.Float32);

        expect(f32.dtype).toBe(DType.Float32);
        expect(f32.get(0)).toBeCloseTo(1.5, 5);
        expect(f32.get(1)).toBeCloseTo(2.5, 5);
        expect(f32.get(2)).toBeCloseTo(3.5, 5);

        f64.dispose();
        f32.dispose();
      });

      it('should return self if dtype is same', async () => {
        const arr = await NDArray.fromArray([1, 2, 3]);
        const same = await arr.astype(DType.Float64);

        // Values should be same
        expect(same.get(0)).toBe(1);
        expect(same.get(1)).toBe(2);
        expect(same.get(2)).toBe(3);

        arr.dispose();
        same.dispose();
      });
    });
  });

  describe('toTypedArray()', () => {
    it('should export to Float64Array', async () => {
      const arr = await NDArray.fromArray([1, 2, 3, 4, 5]);
      const typed = arr.toTypedArray();

      expect(typed).toBeInstanceOf(Float64Array);
      expect(typed.length).toBe(5);
      expect(typed[0]).toBe(1);
      expect(typed[4]).toBe(5);

      arr.dispose();
    });

    it('should export to appropriate typed array for dtype', async () => {
      const f32 = await NDArray.fromArray([1, 2, 3], undefined, { dtype: DType.Float32 });
      const typed = f32.toTypedArray();

      expect(typed).toBeInstanceOf(Float32Array);

      f32.dispose();
    });
  });
});
