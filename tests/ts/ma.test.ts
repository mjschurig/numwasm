/**
 * Tests for NumJS Masked Arrays Module (Phase 16d)
 */

import { describe, it, expect, beforeAll } from 'vitest';
import {
  NDArray,
  DType,
  loadWasmModule,
  ma,
  MaskedArray,
} from '../../dist/numjs.mjs';

// Load WASM module before tests
beforeAll(async () => {
  await loadWasmModule();
});

describe('Types and Constants', () => {
  describe('nomask', () => {
    it('is a symbol', () => {
      expect(typeof ma.nomask).toBe('symbol');
    });
  });

  describe('masked', () => {
    it('is a frozen object', () => {
      expect(Object.isFrozen(ma.masked)).toBe(true);
    });

    it('has __masked__ property', () => {
      expect((ma.masked as any).__masked__).toBe(true);
    });

    it('toString returns "--"', () => {
      expect(ma.masked.toString()).toBe('--');
    });

    it('valueOf returns NaN', () => {
      expect(Number.isNaN(ma.masked.valueOf())).toBe(true);
    });
  });

  describe('isMaskedConstant', () => {
    it('returns true for masked constant', () => {
      expect(ma.isMaskedConstant(ma.masked)).toBe(true);
    });

    it('returns false for other values', () => {
      expect(ma.isMaskedConstant(null)).toBe(false);
      expect(ma.isMaskedConstant(undefined)).toBe(false);
      expect(ma.isMaskedConstant(1)).toBe(false);
      expect(ma.isMaskedConstant({})).toBe(false);
    });
  });

  describe('defaultFillValues', () => {
    it('has correct fill values', () => {
      expect(ma.getDefaultFillValue(DType.Float64)).toBe(1e20);
      expect(ma.getDefaultFillValue(DType.Int32)).toBe(999999);
      expect(ma.getDefaultFillValue(DType.Bool)).toBe(true);
    });
  });

  describe('maximumFillValue', () => {
    it('returns appropriate max values', () => {
      expect(ma.maximumFillValue(DType.Int8)).toBe(127);
      expect(ma.maximumFillValue(DType.Float64)).toBe(Number.POSITIVE_INFINITY);
    });
  });

  describe('minimumFillValue', () => {
    it('returns appropriate min values', () => {
      expect(ma.minimumFillValue(DType.Int8)).toBe(-128);
      expect(ma.minimumFillValue(DType.Float64)).toBe(Number.NEGATIVE_INFINITY);
    });
  });
});

describe('MaskedArray Creation', () => {
  describe('constructor with NDArray', () => {
    it('creates from NDArray without mask', async () => {
      const data = await NDArray.fromArray([1, 2, 3, 4, 5]);
      const arr = new MaskedArray(data);
      expect(arr.shape).toEqual([5]);
      expect(arr.size).toBe(5);
      expect(arr._mask).toBe(ma.nomask);
    });

    it('creates from NDArray with NDArray mask', async () => {
      const data = await NDArray.fromArray([1, 2, 3, 4, 5]);
      const mask = await NDArray.fromArray([0, 1, 0, 0, 1], undefined, { dtype: DType.Bool });
      const arr = new MaskedArray(data, mask);
      expect(arr._mask).not.toBe(ma.nomask);
      expect(arr.getFlat(1)).toBe(ma.masked);
      expect(arr.getFlat(4)).toBe(ma.masked);
    });
  });

  describe('MaskedArray.create (async factory)', () => {
    it('creates from plain array', async () => {
      const arr = await MaskedArray.create([1, 2, 3, 4, 5], [false, true, false, false, true]);
      expect(arr.shape).toEqual([5]);
      expect(arr.getFlat(1)).toBe(ma.masked);
      expect(arr.getFlat(0)).toBe(1);
    });

    it('shrinks all-false mask to nomask', async () => {
      const arr = await MaskedArray.create([1, 2, 3], [false, false, false]);
      expect(arr._mask).toBe(ma.nomask);
    });

    it('respects shrink=false', async () => {
      const arr = await MaskedArray.create([1, 2, 3], [false, false, false], DType.Float64, null, false, false);
      expect(arr._mask).not.toBe(ma.nomask);
    });
  });
});

describe('Element Access', () => {
  describe('getFlat', () => {
    it('returns value for non-masked elements', async () => {
      const arr = await MaskedArray.create([1, 2, 3, 4, 5], [false, true, false, false, true]);
      expect(arr.getFlat(0)).toBe(1);
      expect(arr.getFlat(2)).toBe(3);
    });

    it('returns masked constant for masked elements', async () => {
      const arr = await MaskedArray.create([1, 2, 3, 4, 5], [false, true, false, false, true]);
      expect(arr.getFlat(1)).toBe(ma.masked);
      expect(arr.getFlat(4)).toBe(ma.masked);
    });
  });

  describe('setFlat', () => {
    it('sets value for non-masked elements', async () => {
      const arr = await MaskedArray.create([1, 2, 3, 4, 5], [false, true, false, false, true]);
      arr.setFlat(0, 10);
      expect(arr.getFlat(0)).toBe(10);
    });
  });
});

describe('Mask Operations', () => {
  describe('count', () => {
    it('returns total count for nomask', async () => {
      const arr = await MaskedArray.create([1, 2, 3, 4, 5]);
      expect(arr.count()).toBe(5);
    });

    it('excludes masked elements', async () => {
      const arr = await MaskedArray.create([1, 2, 3, 4, 5], [false, true, false, false, true]);
      expect(arr.count()).toBe(3);
    });
  });

  describe('filled', () => {
    it('returns copy for nomask', async () => {
      const arr = await MaskedArray.create([1, 2, 3, 4, 5]);
      const filled = arr.filled(0);
      expect(filled.getFlat(0)).toBe(1);
    });

    it('replaces masked with fill value', async () => {
      const arr = await MaskedArray.create([1, 2, 3, 4, 5], [false, true, false, false, true]);
      const filled = arr.filled(0);
      expect(filled.getFlat(0)).toBe(1);
      expect(filled.getFlat(1)).toBe(0);
      expect(filled.getFlat(2)).toBe(3);
      expect(filled.getFlat(4)).toBe(0);
    });

    it('uses array fill_value by default', async () => {
      const arr = await MaskedArray.create([1, 2, 3], [false, true, false], DType.Float64, -999);
      const filled = arr.filled();
      expect(filled.getFlat(1)).toBe(-999);
    });
  });

  describe('compressed', () => {
    it('returns all values for nomask', async () => {
      const arr = await MaskedArray.create([1, 2, 3, 4, 5]);
      const compressed = await arr.compressed();
      expect(compressed.size).toBe(5);
    });

    it('excludes masked values', async () => {
      const arr = await MaskedArray.create([1, 2, 3, 4, 5], [false, true, false, false, true]);
      const compressed = await arr.compressed();
      expect(compressed.size).toBe(3);
      expect(compressed.getFlat(0)).toBe(1);
      expect(compressed.getFlat(1)).toBe(3);
      expect(compressed.getFlat(2)).toBe(4);
    });
  });

  describe('shrink_mask', () => {
    it('reduces all-false mask to nomask', async () => {
      const data = await NDArray.fromArray([1, 2, 3]);
      const mask = await NDArray.zeros([3], { dtype: DType.Bool });
      const arr = new MaskedArray(data, mask);
      expect(arr._mask).not.toBe(ma.nomask);
      arr.shrink_mask();
      expect(arr._mask).toBe(ma.nomask);
    });
  });

  describe('harden_mask/soften_mask', () => {
    it('toggles hardmask', async () => {
      const arr = await MaskedArray.create([1, 2, 3], [false, true, false]);
      expect(arr.hardmask).toBe(false);
      arr.harden_mask();
      expect(arr.hardmask).toBe(true);
      arr.soften_mask();
      expect(arr.hardmask).toBe(false);
    });
  });
});

describe('Arithmetic Operations', () => {
  describe('add', () => {
    it('adds scalar', async () => {
      const arr = await MaskedArray.create([1, 2, 3]);
      const result = await arr.add(10);
      expect(result.getFlat(0)).toBe(11);
      expect(result.getFlat(1)).toBe(12);
    });

    it('propagates mask', async () => {
      const arr = await MaskedArray.create([1, 2, 3], [false, true, false]);
      const result = await arr.add(10);
      expect(result.getFlat(0)).toBe(11);
      expect(result.getFlat(1)).toBe(ma.masked);
    });

    it('adds two masked arrays', async () => {
      const a = await MaskedArray.create([1, 2, 3], [false, true, false]);
      const b = await MaskedArray.create([10, 20, 30], [false, false, true]);
      const result = await a.add(b);
      expect(result.getFlat(0)).toBe(11);
      expect(result.getFlat(1)).toBe(ma.masked); // masked in a
      expect(result.getFlat(2)).toBe(ma.masked); // masked in b
    });
  });

  describe('divide', () => {
    it('masks division by zero', async () => {
      const arr = await MaskedArray.create([10, 20, 30]);
      const divisor = await MaskedArray.create([2, 0, 5]);
      const result = await arr.divide(divisor);
      expect(result.getFlat(0)).toBe(5);
      expect(result.getFlat(1)).toBe(ma.masked);
      expect(result.getFlat(2)).toBe(6);
    });
  });

  describe('neg', () => {
    it('negates values', async () => {
      const arr = await MaskedArray.create([1, -2, 3], [false, true, false]);
      const result = await arr.neg();
      expect(result.getFlat(0)).toBe(-1);
      expect(result.getFlat(1)).toBe(ma.masked);
      expect(result.getFlat(2)).toBe(-3);
    });
  });
});

describe('Comparison Operations', () => {
  describe('equal', () => {
    it('compares with scalar', async () => {
      const arr = await MaskedArray.create([1, 2, 3, 2, 1]);
      const result = await arr.equal(2);
      expect(result.getFlat(0)).toBe(0);
      expect(result.getFlat(1)).toBe(1);
      expect(result.getFlat(3)).toBe(1);
    });
  });

  describe('greater', () => {
    it('compares with scalar', async () => {
      const arr = await MaskedArray.create([1, 2, 3, 4, 5]);
      const result = await arr.greater(3);
      expect(result.getFlat(2)).toBe(0);
      expect(result.getFlat(3)).toBe(1);
    });
  });
});

describe('Math Functions', () => {
  describe('sqrt', () => {
    it('computes square root', async () => {
      const arr = await MaskedArray.create([1, 4, 9, 16]);
      const result = await arr.sqrt();
      expect(result.getFlat(0)).toBe(1);
      expect(result.getFlat(1)).toBe(2);
      expect(result.getFlat(2)).toBe(3);
    });

    it('masks negative values', async () => {
      const arr = await MaskedArray.create([4, -1, 9]);
      const result = await arr.sqrt();
      expect(result.getFlat(0)).toBe(2);
      expect(result.getFlat(1)).toBe(ma.masked);
      expect(result.getFlat(2)).toBe(3);
    });
  });

  describe('log', () => {
    it('masks non-positive values', async () => {
      const arr = await MaskedArray.create([1, 0, -1, Math.E]);
      const result = await arr.log();
      expect(result.getFlat(0)).toBeCloseTo(0);
      expect(result.getFlat(1)).toBe(ma.masked);
      expect(result.getFlat(2)).toBe(ma.masked);
      expect(result.getFlat(3)).toBeCloseTo(1);
    });
  });
});

describe('Reductions', () => {
  describe('sum', () => {
    it('sums all elements for nomask', async () => {
      const arr = await MaskedArray.create([1, 2, 3, 4, 5]);
      expect(arr.sum()).toBe(15);
    });

    it('excludes masked elements', async () => {
      const arr = await MaskedArray.create([1, 2, 3, 4, 5], [false, true, false, false, true]);
      expect(arr.sum()).toBe(8); // 1 + 3 + 4
    });
  });

  describe('mean', () => {
    it('computes mean for nomask', async () => {
      const arr = await MaskedArray.create([1, 2, 3, 4, 5]);
      expect(arr.mean()).toBe(3);
    });

    it('excludes masked elements', async () => {
      const arr = await MaskedArray.create([1, 2, 3, 4, 5], [false, true, false, false, true]);
      expect(arr.mean()).toBeCloseTo(8 / 3);
    });
  });

  describe('min/max', () => {
    it('finds min of non-masked', async () => {
      const arr = await MaskedArray.create([5, 1, 3, 2, 4], [false, true, false, false, true]);
      expect(arr.min()).toBe(2);
    });

    it('finds max of non-masked', async () => {
      const arr = await MaskedArray.create([5, 1, 3, 2, 4], [true, false, false, false, true]);
      expect(arr.max()).toBe(3);
    });
  });

  describe('var/std', () => {
    it('computes variance', async () => {
      const arr = await MaskedArray.create([2, 4, 4, 4, 5, 5, 7, 9]);
      expect(arr.var()).toBeCloseTo(4);
    });

    it('computes std', async () => {
      const arr = await MaskedArray.create([2, 4, 4, 4, 5, 5, 7, 9]);
      expect(arr.std()).toBeCloseTo(2);
    });
  });

  describe('prod', () => {
    it('computes product', async () => {
      const arr = await MaskedArray.create([1, 2, 3, 4], [false, true, false, false]);
      expect(arr.prod()).toBe(12); // 1 * 3 * 4
    });
  });

  describe('all/any', () => {
    it('all returns true if all non-masked truthy', async () => {
      const arr = await MaskedArray.create([1, 0, 2], [false, true, false]);
      expect(arr.all()).toBe(true);
    });

    it('any returns true if any non-masked truthy', async () => {
      const arr = await MaskedArray.create([0, 1, 0], [true, false, false]);
      expect(arr.any()).toBe(true);
    });
  });

  describe('argmin/argmax', () => {
    it('finds index of min', async () => {
      const arr = await MaskedArray.create([5, 1, 3], [false, true, false]);
      expect(arr.argmin()).toBe(2);
    });

    it('finds index of max', async () => {
      const arr = await MaskedArray.create([5, 9, 3], [false, true, false]);
      expect(arr.argmax()).toBe(0);
    });
  });
});

describe('Shape Manipulation', () => {
  describe('reshape', () => {
    it('reshapes data and mask', async () => {
      const arr = await MaskedArray.create([1, 2, 3, 4, 5, 6], [false, true, false, false, true, false]);
      const reshaped = arr.reshape([2, 3]);
      expect(reshaped.shape).toEqual([2, 3]);
      expect(reshaped.getFlat(1)).toBe(ma.masked);
    });
  });

  describe('transpose', () => {
    it('transposes 2D array', async () => {
      const arr = await MaskedArray.create([1, 2, 3, 4], [false, true, false, false]);
      const reshaped = arr.reshape([2, 2]);
      const transposed = reshaped.T;
      expect(transposed.shape).toEqual([2, 2]);
    });
  });

  describe('flatten/ravel', () => {
    it('flattens to 1D', async () => {
      const arr = await MaskedArray.create([1, 2, 3, 4]);
      const reshaped = arr.reshape([2, 2]);
      const flat = reshaped.flatten();
      expect(flat.shape).toEqual([4]);
    });
  });
});

describe('Creation Functions', () => {
  describe('masked_array', () => {
    it('creates masked array', async () => {
      const arr = await ma.masked_array([1, 2, 3], [false, true, false]);
      expect(arr.count()).toBe(2);
    });
  });

  describe('masked_equal', () => {
    it('masks elements equal to value', async () => {
      const arr = await ma.masked_equal([1, 2, 3, 2, 1], 2);
      expect(arr.getFlat(0)).toBe(1);
      expect(arr.getFlat(1)).toBe(ma.masked);
      expect(arr.getFlat(3)).toBe(ma.masked);
    });
  });

  describe('masked_greater', () => {
    it('masks elements greater than value', async () => {
      const arr = await ma.masked_greater([1, 2, 3, 4, 5], 3);
      expect(arr.getFlat(2)).toBe(3);
      expect(arr.getFlat(3)).toBe(ma.masked);
      expect(arr.getFlat(4)).toBe(ma.masked);
    });
  });

  describe('masked_inside', () => {
    it('masks elements inside interval', async () => {
      const arr = await ma.masked_inside([1, 2, 3, 4, 5], 2, 4);
      expect(arr.getFlat(0)).toBe(1);
      expect(arr.getFlat(1)).toBe(ma.masked);
      expect(arr.getFlat(2)).toBe(ma.masked);
      expect(arr.getFlat(3)).toBe(ma.masked);
      expect(arr.getFlat(4)).toBe(5);
    });
  });

  describe('masked_invalid', () => {
    it('masks NaN and Inf', async () => {
      const arr = await ma.masked_invalid([1, NaN, 3, Infinity, -Infinity]);
      expect(arr.getFlat(0)).toBe(1);
      expect(arr.getFlat(1)).toBe(ma.masked);
      expect(arr.getFlat(2)).toBe(3);
      expect(arr.getFlat(3)).toBe(ma.masked);
      expect(arr.getFlat(4)).toBe(ma.masked);
    });
  });

  describe('masked_where', () => {
    it('masks where condition is true', async () => {
      const arr = await ma.masked_where([true, false, true], [1, 2, 3]);
      expect(arr.getFlat(0)).toBe(ma.masked);
      expect(arr.getFlat(1)).toBe(2);
      expect(arr.getFlat(2)).toBe(ma.masked);
    });
  });

  describe('zeros/ones/empty', () => {
    it('creates zeros', async () => {
      const arr = await ma.zeros([3]);
      expect(arr.size).toBe(3);
      expect(arr.getFlat(0)).toBe(0);
    });

    it('creates ones', async () => {
      const arr = await ma.ones([3]);
      expect(arr.getFlat(0)).toBe(1);
    });
  });

  describe('masked_all', () => {
    it('creates fully masked array', async () => {
      const arr = await ma.masked_all([3]);
      expect(arr.count()).toBe(0);
      expect(arr.getFlat(0)).toBe(ma.masked);
    });
  });
});

describe('Extras', () => {
  describe('average', () => {
    it('computes simple average', async () => {
      const arr = await MaskedArray.create([1, 2, 3, 4], [false, true, false, false]);
      const avg = await ma.average(arr);
      expect(avg).toBeCloseTo(8 / 3);
    });

    it('computes weighted average', async () => {
      const arr = await MaskedArray.create([1, 2, 3]);
      const avg = await ma.average(arr, null, [1, 2, 1]);
      expect(avg).toBeCloseTo((1 + 4 + 3) / 4);
    });
  });

  describe('median', () => {
    it('computes median', async () => {
      const arr = await MaskedArray.create([1, 2, 3, 4, 5], [false, true, false, false, false]);
      const med = ma.median(arr);
      expect(med).toBe(3.5); // median of [1, 3, 4, 5]
    });
  });

  describe('notmasked_edges', () => {
    it('finds first and last non-masked', async () => {
      const arr = await MaskedArray.create([1, 2, 3, 4], [true, false, false, true]);
      const edges = ma.notmasked_edges(arr);
      expect(edges).toEqual([1, 2]);
    });

    it('returns null for all-masked', async () => {
      const arr = await MaskedArray.create([1, 2], [true, true]);
      const edges = ma.notmasked_edges(arr);
      expect(edges).toBeNull();
    });
  });

  describe('flatnotmasked_contiguous', () => {
    it('finds contiguous non-masked regions', async () => {
      const arr = await MaskedArray.create([1, 2, 3, 4, 5], [true, false, false, true, false]);
      const slices = ma.flatnotmasked_contiguous(arr);
      expect(slices).toEqual([
        { start: 1, stop: 3 },
        { start: 4, stop: 5 },
      ]);
    });
  });

  describe('clump_masked', () => {
    it('finds contiguous masked regions', async () => {
      const arr = await MaskedArray.create([1, 2, 3, 4, 5], [true, false, false, true, true]);
      const slices = ma.clump_masked(arr);
      expect(slices).toEqual([
        { start: 0, stop: 1 },
        { start: 3, stop: 5 },
      ]);
    });
  });
});

describe('toArray', () => {
  it('converts to array with nulls for masked', async () => {
    const arr = await MaskedArray.create([1, 2, 3], [false, true, false]);
    expect(arr.toArray()).toEqual([1, null, 3]);
  });
});

describe('toString', () => {
  it('produces string representation', async () => {
    const arr = await MaskedArray.create([1, 2, 3], [false, true, false]);
    const str = arr.toString();
    expect(str).toContain('masked_array');
    expect(str).toContain('--');
  });
});
