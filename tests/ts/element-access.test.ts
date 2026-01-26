/**
 * Element Access Tests
 *
 * Tests for getting and setting individual array elements.
 */

import { describe, it, expect } from 'vitest';
import { NDArray, DType } from '../../dist/numjs.mjs';

describe('Element Access', () => {
  describe('get() - Multi-dimensional indexing', () => {
    it('should get elements from 1D array', async () => {
      const arr = await NDArray.fromArray([1, 2, 3, 4, 5]);
      expect(arr.get(0)).toBe(1);
      expect(arr.get(2)).toBe(3);
      expect(arr.get(4)).toBe(5);
    });

    it('should get elements from 2D array', async () => {
      const arr = await NDArray.fromArray([[1, 2, 3], [4, 5, 6]]);
      expect(arr.get(0, 0)).toBe(1);
      expect(arr.get(0, 2)).toBe(3);
      expect(arr.get(1, 0)).toBe(4);
      expect(arr.get(1, 2)).toBe(6);
    });

    it('should get elements from 3D array', async () => {
      const arr = await NDArray.fromArray([
        [[1, 2], [3, 4]],
        [[5, 6], [7, 8]],
      ]);
      expect(arr.get(0, 0, 0)).toBe(1);
      expect(arr.get(0, 0, 1)).toBe(2);
      expect(arr.get(0, 1, 0)).toBe(3);
      expect(arr.get(1, 0, 0)).toBe(5);
      expect(arr.get(1, 1, 1)).toBe(8);
    });

    it('should throw for wrong number of indices', async () => {
      const arr = await NDArray.fromArray([[1, 2], [3, 4]]);
      expect(() => arr.get(0)).toThrow('Expected 2 indices');
      expect(() => arr.get(0, 0, 0)).toThrow('Expected 2 indices');
    });

    it('should throw for out-of-bounds indices', async () => {
      const arr = await NDArray.fromArray([[1, 2], [3, 4]]);
      expect(() => arr.get(2, 0)).toThrow('out of bounds');
      expect(() => arr.get(0, 3)).toThrow('out of bounds');
      expect(() => arr.get(-1, 0)).toThrow('out of bounds');
    });
  });

  describe('set() - Multi-dimensional indexing', () => {
    it('should set elements in 1D array', async () => {
      const arr = await NDArray.zeros([5]);
      arr.set(10, 0);
      arr.set(20, 2);
      arr.set(30, 4);
      expect(arr.get(0)).toBe(10);
      expect(arr.get(2)).toBe(20);
      expect(arr.get(4)).toBe(30);
    });

    it('should set elements in 2D array', async () => {
      const arr = await NDArray.zeros([2, 3]);
      arr.set(99, 0, 1);
      arr.set(88, 1, 2);
      expect(arr.get(0, 1)).toBe(99);
      expect(arr.get(1, 2)).toBe(88);
    });

    it('should throw for wrong number of indices', async () => {
      const arr = await NDArray.zeros([2, 3]);
      expect(() => arr.set(1, 0)).toThrow('Expected 2 indices');
    });

    it('should throw for out-of-bounds indices', async () => {
      const arr = await NDArray.zeros([2, 3]);
      expect(() => arr.set(1, 2, 3)).toThrow('out of bounds');
    });
  });

  describe('getFlat() / setFlat() - Flat indexing', () => {
    it('should get elements by flat index', async () => {
      const arr = await NDArray.fromArray([[1, 2, 3], [4, 5, 6]]);
      expect(arr.getFlat(0)).toBe(1);
      expect(arr.getFlat(3)).toBe(4);
      expect(arr.getFlat(5)).toBe(6);
    });

    it('should set elements by flat index', async () => {
      const arr = await NDArray.zeros([2, 3]);
      arr.setFlat(0, 10);
      arr.setFlat(5, 60);
      expect(arr.getFlat(0)).toBe(10);
      expect(arr.getFlat(5)).toBe(60);
    });

    it('should throw for invalid flat index', async () => {
      const arr = await NDArray.zeros([2, 3]);
      expect(() => arr.getFlat(-1)).toThrow('out of bounds');
      expect(() => arr.getFlat(6)).toThrow('out of bounds');
      expect(() => arr.setFlat(-1, 0)).toThrow('out of bounds');
      expect(() => arr.setFlat(6, 0)).toThrow('out of bounds');
    });
  });

  describe('Different dtype element access', () => {
    it('should handle Float32 dtype', async () => {
      const arr = await NDArray.zeros([3], { dtype: DType.Float32 });
      arr.set(1.5, 0);
      arr.set(2.5, 1);
      expect(arr.get(0)).toBeCloseTo(1.5, 5);
      expect(arr.get(1)).toBeCloseTo(2.5, 5);
    });

    it('should handle Int32 dtype', async () => {
      const arr = await NDArray.zeros([3], { dtype: DType.Int32 });
      arr.set(100, 0);
      arr.set(-50, 1);
      expect(arr.get(0)).toBe(100);
      expect(arr.get(1)).toBe(-50);
    });

    it('should truncate floats for integer dtype', async () => {
      const arr = await NDArray.zeros([3], { dtype: DType.Int32 });
      arr.set(3.7, 0);
      arr.set(-2.9, 1);
      expect(arr.get(0)).toBe(3);
      expect(arr.get(1)).toBe(-2);
    });

    it('should handle Int8 dtype with overflow', async () => {
      const arr = await NDArray.zeros([3], { dtype: DType.Int8 });
      arr.set(100, 0);
      arr.set(-100, 1);
      expect(arr.get(0)).toBe(100);
      expect(arr.get(1)).toBe(-100);
    });

    it('should handle Uint8 dtype', async () => {
      const arr = await NDArray.zeros([3], { dtype: DType.Uint8 });
      arr.set(200, 0);
      arr.set(0, 1);
      expect(arr.get(0)).toBe(200);
      expect(arr.get(1)).toBe(0);
    });

    it('should handle Bool dtype', async () => {
      const arr = await NDArray.zeros([4], { dtype: DType.Bool });
      arr.set(1, 0);
      arr.set(0, 1);
      arr.set(5, 2); // Non-zero becomes 1
      arr.set(-1, 3); // Non-zero becomes 1
      expect(arr.get(0)).toBe(1);
      expect(arr.get(1)).toBe(0);
      expect(arr.get(2)).toBe(1);
      expect(arr.get(3)).toBe(1);
    });

    it('should handle Int16 dtype', async () => {
      const arr = await NDArray.zeros([2], { dtype: DType.Int16 });
      arr.set(1000, 0);
      arr.set(-1000, 1);
      expect(arr.get(0)).toBe(1000);
      expect(arr.get(1)).toBe(-1000);
    });

    it('should handle Uint16 dtype', async () => {
      const arr = await NDArray.zeros([2], { dtype: DType.Uint16 });
      arr.set(50000, 0);
      expect(arr.get(0)).toBe(50000);
    });

    it('should handle Uint32 dtype', async () => {
      const arr = await NDArray.zeros([2], { dtype: DType.Uint32 });
      arr.set(3000000000, 0);
      expect(arr.get(0)).toBe(3000000000);
    });
  });

  describe('item() and itemset()', () => {
    it('should get single element with item()', async () => {
      const arr = await NDArray.fromArray([42]);
      expect(arr.item()).toBe(42);
    });

    it('should throw for multi-element array with item()', async () => {
      const arr = await NDArray.fromArray([1, 2, 3]);
      expect(() => arr.item()).toThrow('single-element');
    });

    it('should set single element with itemset()', async () => {
      const arr = await NDArray.zeros([1]);
      arr.itemset(99);
      expect(arr.item()).toBe(99);
    });

    it('should throw for multi-element array with itemset()', async () => {
      const arr = await NDArray.zeros([3]);
      expect(() => arr.itemset(77)).toThrow('single-element');
    });
  });

  describe('Complex number access', () => {
    it('should get complex values with getComplex()', async () => {
      const arr = await NDArray.zeros([2], { dtype: DType.Complex128 });
      // Set using setComplex
      arr.setComplex(3, 4, 0);
      arr.setComplex(1, -2, 1);

      const c0 = arr.getComplex(0);
      expect(c0.real).toBeCloseTo(3, 10);
      expect(c0.imag).toBeCloseTo(4, 10);

      const c1 = arr.getComplex(1);
      expect(c1.real).toBeCloseTo(1, 10);
      expect(c1.imag).toBeCloseTo(-2, 10);
    });

    it('should set complex values with setComplex()', async () => {
      const arr = await NDArray.zeros([2], { dtype: DType.Complex64 });
      arr.setComplex(5.5, 2.5, 0);
      arr.setComplex(-1.0, 0.0, 1);

      const c0 = arr.getComplex(0);
      expect(c0.real).toBeCloseTo(5.5, 5);
      expect(c0.imag).toBeCloseTo(2.5, 5);

      const c1 = arr.getComplex(1);
      expect(c1.real).toBeCloseTo(-1.0, 5);
      expect(c1.imag).toBeCloseTo(0.0, 5);
    });

    it('should throw for non-complex dtype with getComplex()', async () => {
      const arr = await NDArray.zeros([2], { dtype: DType.Float64 });
      expect(() => arr.getComplex(0)).toThrow('complex');
    });

    it('should throw for non-complex dtype with setComplex()', async () => {
      const arr = await NDArray.zeros([2], { dtype: DType.Float64 });
      expect(() => arr.setComplex(1, 2, 0)).toThrow('complex');
    });
  });

  describe('Array properties', () => {
    it('should report correct nbytes', async () => {
      const f64 = await NDArray.zeros([10], { dtype: DType.Float64 });
      expect(f64.nbytes).toBe(80); // 10 * 8 bytes

      const f32 = await NDArray.zeros([10], { dtype: DType.Float32 });
      expect(f32.nbytes).toBe(40); // 10 * 4 bytes

      const i8 = await NDArray.zeros([10], { dtype: DType.Int8 });
      expect(i8.nbytes).toBe(10); // 10 * 1 byte

      const c128 = await NDArray.zeros([5], { dtype: DType.Complex128 });
      expect(c128.nbytes).toBe(80); // 5 * 16 bytes
    });

    it('should report correct itemsize', async () => {
      const f64 = await NDArray.zeros([1], { dtype: DType.Float64 });
      expect(f64.itemsize).toBe(8);

      const f32 = await NDArray.zeros([1], { dtype: DType.Float32 });
      expect(f32.itemsize).toBe(4);

      const i32 = await NDArray.zeros([1], { dtype: DType.Int32 });
      expect(i32.itemsize).toBe(4);

      const i8 = await NDArray.zeros([1], { dtype: DType.Int8 });
      expect(i8.itemsize).toBe(1);
    });

    it('should report correct strides', async () => {
      const arr2d = await NDArray.zeros([3, 4], { dtype: DType.Float64 });
      const strides = arr2d.strides;
      // C-contiguous: strides are [4*8, 8] = [32, 8]
      expect(strides).toEqual([32, 8]);

      const arr3d = await NDArray.zeros([2, 3, 4], { dtype: DType.Float32 });
      const strides3d = arr3d.strides;
      // C-contiguous: [3*4*4, 4*4, 4] = [48, 16, 4]
      expect(strides3d).toEqual([48, 16, 4]);
    });

    it('should report correct flags', async () => {
      const arr = await NDArray.zeros([3, 4]);
      const flags = arr.flags;
      expect(flags.c_contiguous).toBe(true);
      expect(flags.owndata).toBe(true);
      expect(flags.writeable).toBe(true);
    });

    it('should report correct dtypeName', async () => {
      const f64 = await NDArray.zeros([1], { dtype: DType.Float64 });
      expect(f64.dtypeName).toBe('float64');

      const i32 = await NDArray.zeros([1], { dtype: DType.Int32 });
      expect(i32.dtypeName).toBe('int32');

      const c64 = await NDArray.zeros([1], { dtype: DType.Complex64 });
      expect(c64.dtypeName).toBe('complex64');
    });
  });
});
