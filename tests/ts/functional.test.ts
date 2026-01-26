/**
 * Functional Programming Tests (Level 10)
 *
 * Tests for applyAlongAxis, applyOverAxes, Vectorize/vectorize,
 * frompyfunc, and piecewise.
 */

import { describe, it, expect } from 'vitest';
import {
  NDArray,
  DType,
  applyAlongAxis,
  applyOverAxes,
  Vectorize,
  vectorize,
  frompyfunc,
  piecewise,
  sum,
  less,
  greater_equal,
  negative,
} from '../../dist/numjs.mjs';

describe('Level 10: Functional Programming', () => {
  describe('applyAlongAxis', () => {
    it('should apply function along axis 0 (columns)', async () => {
      const arr = await NDArray.fromArray([
        [1, 2, 3],
        [4, 5, 6],
      ]);

      // Sum each column
      const result = await applyAlongAxis(
        (x) => x.toArray().reduce((a, b) => a + b, 0),
        0,
        arr
      );

      expect(result.shape).toEqual([3]);
      expect(result.toArray()).toEqual([5, 7, 9]);

      arr.dispose();
      result.dispose();
    });

    it('should apply function along axis 1 (rows)', async () => {
      const arr = await NDArray.fromArray([
        [1, 2, 3],
        [4, 5, 6],
      ]);

      // Sum each row
      const result = await applyAlongAxis(
        (x) => x.toArray().reduce((a, b) => a + b, 0),
        1,
        arr
      );

      expect(result.shape).toEqual([2]);
      expect(result.toArray()).toEqual([6, 15]);

      arr.dispose();
      result.dispose();
    });

    it('should apply async function returning array', async () => {
      const arr = await NDArray.fromArray([
        [3, 1, 2],
        [6, 4, 5],
      ]);

      // Sort each row
      const result = await applyAlongAxis(
        async (x) => {
          const data = [...x.toArray()].sort((a, b) => a - b);
          return NDArray.fromArray(data);
        },
        1,
        arr
      );

      expect(result.shape).toEqual([2, 3]);
      expect(result.toArray()).toEqual([1, 2, 3, 4, 5, 6]);

      arr.dispose();
      result.dispose();
    });

    it('should handle negative axis', async () => {
      const arr = await NDArray.fromArray([
        [1, 2],
        [3, 4],
      ]);

      // Sum along last axis (axis=-1 is axis=1)
      const result = await applyAlongAxis(
        (x) => x.toArray().reduce((a, b) => a + b, 0),
        -1,
        arr
      );

      expect(result.shape).toEqual([2]);
      expect(result.toArray()).toEqual([3, 7]);

      arr.dispose();
      result.dispose();
    });

    it('should handle 1-D array', async () => {
      const arr = await NDArray.fromArray([1, 2, 3, 4]);

      const result = await applyAlongAxis(
        (x) => x.toArray().reduce((a, b) => a + b, 0),
        0,
        arr
      );

      expect(result.toArray()).toEqual([10]);

      arr.dispose();
      result.dispose();
    });

    it('should pass additional arguments to function', async () => {
      const arr = await NDArray.fromArray([
        [1, 2],
        [3, 4],
      ]);

      const result = await applyAlongAxis(
        (x, multiplier) =>
          x.toArray().reduce((a, b) => a + b, 0) * (multiplier as number),
        0,
        arr,
        10 // multiplier
      );

      expect(result.shape).toEqual([2]);
      expect(result.toArray()).toEqual([40, 60]); // [4*10, 6*10]

      arr.dispose();
      result.dispose();
    });

    it('should throw on 0-D array', async () => {
      const arr = await NDArray.fromArray(5, []);

      await expect(
        applyAlongAxis((x) => x.item(), 0, arr)
      ).rejects.toThrow('at least 1-dimensional');

      arr.dispose();
    });

    it('should throw on invalid axis', async () => {
      const arr = await NDArray.fromArray([1, 2, 3]);

      await expect(
        applyAlongAxis((x) => x.item(), 5, arr)
      ).rejects.toThrow('out of bounds');

      arr.dispose();
    });
  });

  describe('applyOverAxes', () => {
    it('should apply function over single axis', async () => {
      const arr = await NDArray.fromArray([
        [1, 2, 3],
        [4, 5, 6],
      ]);

      const result = await applyOverAxes(
        async (a, axis) => sum(a, axis),
        arr,
        0
      );

      // After summing axis 0, shape is [3], expanded to [1, 3]
      expect(result.ndim).toBe(2);
      expect(result.shape).toEqual([1, 3]);
      expect(result.toArray()).toEqual([5, 7, 9]);

      arr.dispose();
      result.dispose();
    });

    it('should apply function over multiple axes', async () => {
      const arr = await NDArray.fromArray([
        [[1, 2], [3, 4]],
        [[5, 6], [7, 8]],
      ]);

      const result = await applyOverAxes(
        async (a, axis) => sum(a, axis),
        arr,
        [0, 2]
      );

      // Sum axis 0: shape [2, 2] -> expanded to [1, 2, 2]
      // Sum axis 2: shape [1, 2] -> expanded to [1, 2, 1]
      expect(result.ndim).toBe(3);
      expect(result.shape).toEqual([1, 2, 1]);

      arr.dispose();
      result.dispose();
    });

    it('should return view for empty axes', async () => {
      const arr = await NDArray.fromArray([1, 2, 3]);

      const result = await applyOverAxes(async (a, axis) => sum(a, axis), arr, []);

      expect(result.shape).toEqual([3]);
      expect(result.toArray()).toEqual([1, 2, 3]);

      arr.dispose();
      result.dispose();
    });

    it('should preserve dimensionality', async () => {
      const arr = await NDArray.fromArray([
        [1, 2],
        [3, 4],
      ]);

      const result = await applyOverAxes(
        async (a, axis) => sum(a, axis),
        arr,
        [0, 1]
      );

      // After all reductions, should still be 2D
      expect(result.ndim).toBe(2);
      expect(result.shape).toEqual([1, 1]);
      expect(result.item()).toBe(10);

      arr.dispose();
      result.dispose();
    });
  });

  describe('Vectorize class', () => {
    it('should vectorize simple function', async () => {
      const vmin = new Vectorize((a: unknown, b: unknown) =>
        Math.min(a as number, b as number)
      );

      const arr1 = await NDArray.fromArray([1, 4, 3]);
      const arr2 = await NDArray.fromArray([2, 2, 5]);

      const result = (await vmin.call(arr1, arr2)) as NDArray;

      expect(result.toArray()).toEqual([1, 2, 3]);

      arr1.dispose();
      arr2.dispose();
      result.dispose();
    });

    it('should handle broadcasting', async () => {
      const vadd = new Vectorize((a: unknown, b: unknown) =>
        (a as number) + (b as number)
      );

      const arr1 = await NDArray.fromArray([
        [1, 2],
        [3, 4],
      ]);
      const arr2 = await NDArray.fromArray([10, 20]);

      const result = (await vadd.call(arr1, arr2)) as NDArray;

      expect(result.shape).toEqual([2, 2]);
      expect(result.toArray()).toEqual([11, 22, 13, 24]);

      arr1.dispose();
      arr2.dispose();
      result.dispose();
    });

    it('should handle scalar inputs', async () => {
      const vmul = new Vectorize((a: unknown, b: unknown) =>
        (a as number) * (b as number)
      );

      const arr = await NDArray.fromArray([1, 2, 3]);

      const result = (await vmul.call(arr, 10)) as NDArray;

      expect(result.toArray()).toEqual([10, 20, 30]);

      arr.dispose();
      result.dispose();
    });

    it('should respect otypes option', async () => {
      const visPositive = new Vectorize((x: unknown) => (x as number) > 0, {
        otypes: [DType.Bool],
      });

      const arr = await NDArray.fromArray([-1, 0, 1, 2]);

      const result = (await visPositive.call(arr)) as NDArray;

      expect(result.dtype).toBe(DType.Bool);
      expect(result.toArray()).toEqual([0, 0, 1, 1]);

      arr.dispose();
      result.dispose();
    });

    it('should handle excluded arguments', async () => {
      // Multiply array by a scalar factor, but don't vectorize the factor
      const vscale = new Vectorize(
        (x: unknown, factor: unknown) => (x as number) * (factor as number),
        { excluded: new Set([1]) }
      );

      const arr = await NDArray.fromArray([1, 2, 3]);

      const result = (await vscale.call(arr, 5)) as NDArray;

      expect(result.toArray()).toEqual([5, 10, 15]);

      arr.dispose();
      result.dispose();
    });

    it('should handle multiple outputs', async () => {
      const vdivmod = new Vectorize((x: unknown, y: unknown) => {
        const xn = x as number;
        const yn = y as number;
        return [Math.floor(xn / yn), xn % yn];
      });

      const arr1 = await NDArray.fromArray([10, 11, 12]);
      const arr2 = await NDArray.fromArray([3]);

      const result = (await vdivmod.call(arr1, arr2)) as NDArray[];

      expect(result.length).toBe(2);
      expect(result[0].toArray()).toEqual([3, 3, 4]);
      expect(result[1].toArray()).toEqual([1, 2, 0]);

      arr1.dispose();
      arr2.dispose();
      result[0].dispose();
      result[1].dispose();
    });

    it('should cache first result when cache=true', async () => {
      let callCount = 0;
      const vcached = new Vectorize(
        (x: unknown) => {
          callCount++;
          return (x as number) * 2;
        },
        { cache: true }
      );

      const arr = await NDArray.fromArray([1, 2, 3]);

      // First call - infers types and caches first result
      const result1 = (await vcached.call(arr)) as NDArray;
      const firstCallCount = callCount;

      // Reset counter
      callCount = 0;

      // Second call on same instance should use cached type info
      const result2 = (await vcached.call(arr)) as NDArray;

      expect(result1.toArray()).toEqual([2, 4, 6]);
      expect(result2.toArray()).toEqual([2, 4, 6]);

      arr.dispose();
      result1.dispose();
      result2.dispose();
    });
  });

  describe('vectorize factory', () => {
    it('should create Vectorize instance', async () => {
      const vsquare = vectorize((x: unknown) => (x as number) ** 2);

      const arr = await NDArray.fromArray([1, 2, 3, 4]);

      const result = (await vsquare.call(arr)) as NDArray;

      expect(result.toArray()).toEqual([1, 4, 9, 16]);

      arr.dispose();
      result.dispose();
    });
  });

  describe('frompyfunc', () => {
    it('should create ufunc-like with 2 inputs, 1 output', async () => {
      const addone = frompyfunc((x, y) => x + y + 1, 2, 1);

      expect(addone.nin).toBe(2);
      expect(addone.nout).toBe(1);

      const arr1 = await NDArray.fromArray([1, 2, 3]);
      const arr2 = await NDArray.fromArray([10, 20, 30]);

      const result = (await addone(arr1, arr2)) as NDArray;

      expect(result.toArray()).toEqual([12, 23, 34]);

      arr1.dispose();
      arr2.dispose();
      result.dispose();
    });

    it('should create ufunc-like with multiple outputs', async () => {
      const divmod = frompyfunc((x, y) => [Math.floor(x / y), x % y], 2, 2);

      expect(divmod.nin).toBe(2);
      expect(divmod.nout).toBe(2);

      const arr1 = await NDArray.fromArray([10, 11, 12]);
      const arr2 = await NDArray.fromArray([3]);

      const [quot, rem] = (await divmod(arr1, arr2)) as NDArray[];

      expect(quot.toArray()).toEqual([3, 3, 4]);
      expect(rem.toArray()).toEqual([1, 2, 0]);

      arr1.dispose();
      arr2.dispose();
      quot.dispose();
      rem.dispose();
    });

    it('should handle broadcasting', async () => {
      const add = frompyfunc((x, y) => x + y, 2, 1);

      const arr1 = await NDArray.fromArray([
        [1, 2],
        [3, 4],
      ]);
      const arr2 = await NDArray.fromArray([10]);

      const result = (await add(arr1, arr2)) as NDArray;

      expect(result.shape).toEqual([2, 2]);
      expect(result.toArray()).toEqual([11, 12, 13, 14]);

      arr1.dispose();
      arr2.dispose();
      result.dispose();
    });

    it('should handle scalar inputs', async () => {
      const add = frompyfunc((x, y) => x + y, 2, 1);

      const arr = await NDArray.fromArray([1, 2, 3]);

      const result = (await add(arr, 100)) as NDArray;

      expect(result.toArray()).toEqual([101, 102, 103]);

      arr.dispose();
      result.dispose();
    });

    it('should throw on wrong number of arguments', async () => {
      const add = frompyfunc((x, y) => x + y, 2, 1);

      const arr = await NDArray.fromArray([1, 2, 3]);

      await expect(add(arr)).rejects.toThrow('Expected 2 arguments');

      arr.dispose();
    });

    it('should store identity value', () => {
      const add = frompyfunc((x, y) => x + y, 2, 1, 0);
      expect(add.identity).toBe(0);

      const mul = frompyfunc((x, y) => x * y, 2, 1, 1);
      expect(mul.identity).toBe(1);
    });

    it('should throw on invalid nin/nout', () => {
      expect(() => frompyfunc((x) => x, 0, 1)).toThrow('nin must be at least 1');
      expect(() => frompyfunc((x) => x, 1, 0)).toThrow('nout must be at least 1');
    });
  });

  describe('piecewise', () => {
    it('should evaluate with constant values', async () => {
      const x = await NDArray.fromArray([-2, -1, 0, 1, 2]);

      const cond1 = await less(x, await NDArray.fromArray([0]));
      const cond2 = await greater_equal(x, await NDArray.fromArray([0]));

      const result = await piecewise(x, [cond1, cond2], [-1, 1]);

      expect(result.toArray()).toEqual([-1, -1, 1, 1, 1]);

      x.dispose();
      cond1.dispose();
      cond2.dispose();
      result.dispose();
    });

    it('should evaluate with callable functions', async () => {
      const x = await NDArray.fromArray([-2, -1, 0, 1, 2]);

      const cond1 = await less(x, await NDArray.fromArray([0]));
      const cond2 = await greater_equal(x, await NDArray.fromArray([0]));

      // Absolute value using piecewise
      const result = await piecewise(
        x,
        [cond1, cond2],
        [
          async (vals) => negative(vals), // negate for x < 0
          async (vals) => vals.copy(), // identity for x >= 0
        ]
      );

      expect(result.toArray()).toEqual([2, 1, 0, 1, 2]);

      x.dispose();
      cond1.dispose();
      cond2.dispose();
      result.dispose();
    });

    it('should handle default (otherwise) case', async () => {
      const x = await NDArray.fromArray([1, 2, 3, 4, 5]);

      // Only one condition for x == 3
      const cond = await NDArray.fromArray([0, 0, 1, 0, 0], [5], {
        dtype: DType.Bool,
      });

      // Two functions: one for condition, one for default
      const result = await piecewise(
        x,
        [cond],
        [
          100, // where x == 3
          0, // everywhere else (default)
        ]
      );

      expect(result.toArray()).toEqual([0, 0, 100, 0, 0]);

      x.dispose();
      cond.dispose();
      result.dispose();
    });

    it('should handle mixed constants and functions', async () => {
      const x = await NDArray.fromArray([1, 2, 3, 4, 5]);

      const cond1 = await NDArray.fromArray([1, 1, 0, 0, 0], [5], {
        dtype: DType.Bool,
      });
      const cond2 = await NDArray.fromArray([0, 0, 1, 1, 1], [5], {
        dtype: DType.Bool,
      });

      const result = await piecewise(
        x,
        [cond1, cond2],
        [
          -1, // constant for cond1
          async (vals) => {
            // double for cond2
            const data = vals.toArray().map((v) => v * 2);
            return NDArray.fromArray(data);
          },
        ]
      );

      expect(result.toArray()).toEqual([-1, -1, 6, 8, 10]);

      x.dispose();
      cond1.dispose();
      cond2.dispose();
      result.dispose();
    });

    it('should throw on empty condlist', async () => {
      const x = await NDArray.fromArray([1, 2, 3]);

      await expect(piecewise(x, [], [0])).rejects.toThrow('at least one condition');

      x.dispose();
    });

    it('should throw on mismatched lengths', async () => {
      const x = await NDArray.fromArray([1, 2, 3]);
      const cond = await NDArray.fromArray([1, 0, 1], [3], { dtype: DType.Bool });

      // 3 functions for 1 condition is invalid (should be 1 or 2)
      await expect(piecewise(x, [cond], [1, 2, 3])).rejects.toThrow(
        'same length as condlist'
      );

      x.dispose();
      cond.dispose();
    });

    it('should handle no matching conditions (all zero)', async () => {
      const x = await NDArray.fromArray([1, 2, 3]);
      const cond = await NDArray.fromArray([0, 0, 0], [3], { dtype: DType.Bool });

      const result = await piecewise(x, [cond], [99]);

      // No condition matches, output is zeros (default initialization)
      expect(result.toArray()).toEqual([0, 0, 0]);

      x.dispose();
      cond.dispose();
      result.dispose();
    });

    it('should pass additional arguments to functions', async () => {
      const x = await NDArray.fromArray([1, 2, 3, 4]);
      const cond = await NDArray.fromArray([1, 1, 1, 1], [4], { dtype: DType.Bool });

      const result = await piecewise(
        x,
        [cond],
        [
          async (vals, multiplier) => {
            const data = vals.toArray().map((v) => v * (multiplier as number));
            return NDArray.fromArray(data);
          },
        ],
        10 // extra argument: multiplier
      );

      expect(result.toArray()).toEqual([10, 20, 30, 40]);

      x.dispose();
      cond.dispose();
      result.dispose();
    });
  });
});
