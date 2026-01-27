/**
 * Tests for Universal Functions (Ufuncs)
 */

import { describe, it, expect, beforeAll, afterEach } from 'vitest';
import {
  loadWasmModule,
  NDArray,
  DType,
  // Unary ufuncs
  negative,
  positive,
  absolute,
  abs,
  sign,
  sqrt,
  square,
  cbrt,
  reciprocal,
  exp,
  exp2,
  expm1,
  log,
  log2,
  log10,
  log1p,
  sin,
  cos,
  tan,
  arcsin,
  arccos,
  arctan,
  sinh,
  cosh,
  tanh,
  arcsinh,
  arccosh,
  arctanh,
  floor,
  ceil,
  trunc,
  rint,
  degrees,
  radians,
  logical_not,
  // Binary ufuncs
  add,
  subtract,
  multiply,
  divide,
  floor_divide,
  remainder,
  power,
  equal,
  not_equal,
  less,
  less_equal,
  greater,
  greater_equal,
  maximum,
  minimum,
  logical_and,
  logical_or,
  arctan2,
  hypot,
} from '../../dist/numjs.mjs';

// Track resources for cleanup
const resources: NDArray[] = [];
function track<T extends NDArray>(arr: T): T {
  resources.push(arr);
  return arr;
}

describe('Ufunc Tests', () => {
  beforeAll(async () => {
    await loadWasmModule();
  });

  afterEach(() => {
    for (const arr of resources) {
      if (!arr.isDisposed) {
        arr.dispose();
      }
    }
    resources.length = 0;
  });

  describe('Unary Arithmetic', () => {
    it('negative returns element-wise negation', async () => {
      const arr = track(await NDArray.fromArray([1, -2, 3, -4]));
      const result = track(negative(arr));
      expect(result.toArray()).toEqual([-1, 2, -3, 4]);
    });

    it('positive returns a copy', async () => {
      const arr = track(await NDArray.fromArray([1, -2, 3]));
      const result = track(positive(arr));
      expect(result.toArray()).toEqual([1, -2, 3]);
    });

    it('absolute returns element-wise absolute value', async () => {
      const arr = track(await NDArray.fromArray([-1, 2, -3, 4, 0]));
      const result = track(absolute(arr));
      expect(result.toArray()).toEqual([1, 2, 3, 4, 0]);
    });

    it('abs is an alias for absolute', async () => {
      const arr = track(await NDArray.fromArray([-5, 5]));
      const result = track(abs(arr));
      expect(result.toArray()).toEqual([5, 5]);
    });

    it('sign returns -1, 0, or +1', async () => {
      const arr = track(await NDArray.fromArray([-5, 0, 5]));
      const result = track(sign(arr));
      expect(result.toArray()).toEqual([-1, 0, 1]);
    });
  });

  describe('Powers and Roots', () => {
    it('sqrt computes element-wise square root', async () => {
      const arr = track(await NDArray.fromArray([0, 1, 4, 9, 16]));
      const result = track(sqrt(arr));
      expect(result.toArray()).toEqual([0, 1, 2, 3, 4]);
    });

    it('square computes element-wise x^2', async () => {
      const arr = track(await NDArray.fromArray([0, 1, 2, 3, -2]));
      const result = track(square(arr));
      expect(result.toArray()).toEqual([0, 1, 4, 9, 4]);
    });

    it('cbrt computes element-wise cube root', async () => {
      const arr = track(await NDArray.fromArray([0, 1, 8, 27]));
      const result = track(cbrt(arr));
      const expected = [0, 1, 2, 3];
      result.toArray().forEach((val, i) => {
        expect(val).toBeCloseTo(expected[i], 10);
      });
    });

    it('reciprocal computes element-wise 1/x', async () => {
      const arr = track(await NDArray.fromArray([1, 2, 4, 0.5]));
      const result = track(reciprocal(arr));
      expect(result.toArray()).toEqual([1, 0.5, 0.25, 2]);
    });
  });

  describe('Exponential and Logarithmic', () => {
    it('exp computes e^x', async () => {
      const arr = track(await NDArray.fromArray([0, 1]));
      const result = track(exp(arr));
      expect(result.toArray()[0]).toBeCloseTo(1, 10);
      expect(result.toArray()[1]).toBeCloseTo(Math.E, 10);
    });

    it('log computes natural logarithm', async () => {
      const arr = track(await NDArray.fromArray([1, Math.E, Math.E * Math.E]));
      const result = track(log(arr));
      expect(result.toArray()[0]).toBeCloseTo(0, 10);
      expect(result.toArray()[1]).toBeCloseTo(1, 10);
      expect(result.toArray()[2]).toBeCloseTo(2, 10);
    });

    it('log2 computes base-2 logarithm', async () => {
      const arr = track(await NDArray.fromArray([1, 2, 4, 8]));
      const result = track(log2(arr));
      expect(result.toArray()).toEqual([0, 1, 2, 3]);
    });

    it('log10 computes base-10 logarithm', async () => {
      const arr = track(await NDArray.fromArray([1, 10, 100, 1000]));
      const result = track(log10(arr));
      result.toArray().forEach((val, i) => {
        expect(val).toBeCloseTo(i, 10);
      });
    });
  });

  describe('Trigonometric', () => {
    it('sin computes element-wise sine', async () => {
      const arr = track(await NDArray.fromArray([0, Math.PI / 2, Math.PI]));
      const result = track(sin(arr));
      expect(result.toArray()[0]).toBeCloseTo(0, 10);
      expect(result.toArray()[1]).toBeCloseTo(1, 10);
      expect(result.toArray()[2]).toBeCloseTo(0, 10);
    });

    it('cos computes element-wise cosine', async () => {
      const arr = track(await NDArray.fromArray([0, Math.PI / 2, Math.PI]));
      const result = track(cos(arr));
      expect(result.toArray()[0]).toBeCloseTo(1, 10);
      expect(result.toArray()[1]).toBeCloseTo(0, 10);
      expect(result.toArray()[2]).toBeCloseTo(-1, 10);
    });

    it('arcsin is inverse of sin', async () => {
      const arr = track(await NDArray.fromArray([0, 0.5, 1]));
      const result = track(arcsin(arr));
      expect(result.toArray()[0]).toBeCloseTo(0, 10);
      expect(result.toArray()[1]).toBeCloseTo(Math.PI / 6, 10);
      expect(result.toArray()[2]).toBeCloseTo(Math.PI / 2, 10);
    });
  });

  describe('Rounding', () => {
    it('floor rounds down', async () => {
      const arr = track(await NDArray.fromArray([1.7, -1.7, 2.0, -2.0]));
      const result = track(floor(arr));
      expect(result.toArray()).toEqual([1, -2, 2, -2]);
    });

    it('ceil rounds up', async () => {
      const arr = track(await NDArray.fromArray([1.1, -1.1, 2.0, -2.0]));
      const result = track(ceil(arr));
      expect(result.toArray()).toEqual([2, -1, 2, -2]);
    });

    it('trunc truncates toward zero', async () => {
      const arr = track(await NDArray.fromArray([1.7, -1.7, 2.0, -2.0]));
      const result = track(trunc(arr));
      expect(result.toArray()).toEqual([1, -1, 2, -2]);
    });
  });

  describe('Angle Conversion', () => {
    it('degrees converts radians to degrees', async () => {
      const arr = track(await NDArray.fromArray([0, Math.PI / 2, Math.PI, 2 * Math.PI]));
      const result = track(degrees(arr));
      expect(result.toArray()[0]).toBeCloseTo(0, 10);
      expect(result.toArray()[1]).toBeCloseTo(90, 10);
      expect(result.toArray()[2]).toBeCloseTo(180, 10);
      expect(result.toArray()[3]).toBeCloseTo(360, 10);
    });

    it('radians converts degrees to radians', async () => {
      const arr = track(await NDArray.fromArray([0, 90, 180, 360]));
      const result = track(radians(arr));
      expect(result.toArray()[0]).toBeCloseTo(0, 10);
      expect(result.toArray()[1]).toBeCloseTo(Math.PI / 2, 10);
      expect(result.toArray()[2]).toBeCloseTo(Math.PI, 10);
      expect(result.toArray()[3]).toBeCloseTo(2 * Math.PI, 10);
    });
  });

  describe('Binary Arithmetic', () => {
    it('add performs element-wise addition', async () => {
      const a = track(await NDArray.fromArray([1, 2, 3]));
      const b = track(await NDArray.fromArray([4, 5, 6]));
      const result = track(add(a, b));
      expect(result.toArray()).toEqual([5, 7, 9]);
    });

    it('subtract performs element-wise subtraction', async () => {
      const a = track(await NDArray.fromArray([5, 5, 5]));
      const b = track(await NDArray.fromArray([1, 2, 3]));
      const result = track(subtract(a, b));
      expect(result.toArray()).toEqual([4, 3, 2]);
    });

    it('multiply performs element-wise multiplication', async () => {
      const a = track(await NDArray.fromArray([1, 2, 3]));
      const b = track(await NDArray.fromArray([2, 3, 4]));
      const result = track(multiply(a, b));
      expect(result.toArray()).toEqual([2, 6, 12]);
    });

    it('divide performs element-wise division', async () => {
      const a = track(await NDArray.fromArray([6, 8, 10]));
      const b = track(await NDArray.fromArray([2, 4, 5]));
      const result = track(divide(a, b));
      expect(result.toArray()).toEqual([3, 2, 2]);
    });

    it('power raises to element-wise power', async () => {
      const a = track(await NDArray.fromArray([2, 3, 4]));
      const b = track(await NDArray.fromArray([2, 2, 2]));
      const result = track(power(a, b));
      expect(result.toArray()).toEqual([4, 9, 16]);
    });

    it('floor_divide performs integer division', async () => {
      const a = track(await NDArray.fromArray([7, 8, 9]));
      const b = track(await NDArray.fromArray([2, 3, 4]));
      const result = track(floor_divide(a, b));
      expect(result.toArray()).toEqual([3, 2, 2]);
    });

    it('remainder computes modulo', async () => {
      const a = track(await NDArray.fromArray([7, 8, 9]));
      const b = track(await NDArray.fromArray([2, 3, 4]));
      const result = track(remainder(a, b));
      expect(result.toArray()).toEqual([1, 2, 1]);
    });
  });

  describe('Comparisons', () => {
    it('equal returns element-wise equality', async () => {
      const a = track(await NDArray.fromArray([1, 2, 3]));
      const b = track(await NDArray.fromArray([1, 0, 3]));
      const result = track(equal(a, b));
      expect(result.toArray()).toEqual([1, 0, 1]);
    });

    it('less returns element-wise less than', async () => {
      const a = track(await NDArray.fromArray([1, 2, 3]));
      const b = track(await NDArray.fromArray([2, 2, 2]));
      const result = track(less(a, b));
      expect(result.toArray()).toEqual([1, 0, 0]);
    });

    it('greater returns element-wise greater than', async () => {
      const a = track(await NDArray.fromArray([1, 2, 3]));
      const b = track(await NDArray.fromArray([2, 2, 2]));
      const result = track(greater(a, b));
      expect(result.toArray()).toEqual([0, 0, 1]);
    });
  });

  describe('Extrema', () => {
    it('maximum returns element-wise max', async () => {
      const a = track(await NDArray.fromArray([1, 5, 3]));
      const b = track(await NDArray.fromArray([4, 2, 3]));
      const result = track(maximum(a, b));
      expect(result.toArray()).toEqual([4, 5, 3]);
    });

    it('minimum returns element-wise min', async () => {
      const a = track(await NDArray.fromArray([1, 5, 3]));
      const b = track(await NDArray.fromArray([4, 2, 3]));
      const result = track(minimum(a, b));
      expect(result.toArray()).toEqual([1, 2, 3]);
    });
  });

  describe('Broadcasting', () => {
    it('add broadcasts scalar-like array', async () => {
      const a = track(await NDArray.fromArray([[1, 2], [3, 4]]));
      const b = track(await NDArray.fromArray([10, 20]));
      const result = track(add(a, b));
      expect(result.shape).toEqual([2, 2]);
      expect(result.toArray()).toEqual([11, 22, 13, 24]);
    });

    it('multiply broadcasts different shapes', async () => {
      const a = track(await NDArray.fromArray([[1], [2], [3]]));  // 3x1
      const b = track(await NDArray.fromArray([1, 2, 3, 4]));     // 4
      const result = track(multiply(a, b));
      expect(result.shape).toEqual([3, 4]);
      expect(result.toArray()).toEqual([
        1, 2, 3, 4,
        2, 4, 6, 8,
        3, 6, 9, 12
      ]);
    });
  });

  describe('Special Math Functions', () => {
    it('arctan2 computes atan2(y, x)', async () => {
      const y = track(await NDArray.fromArray([0, 1, 0, -1]));
      const x = track(await NDArray.fromArray([1, 0, -1, 0]));
      const result = track(arctan2(y, x));
      expect(result.toArray()[0]).toBeCloseTo(0, 10);
      expect(result.toArray()[1]).toBeCloseTo(Math.PI / 2, 10);
      expect(result.toArray()[2]).toBeCloseTo(Math.PI, 10);
      expect(result.toArray()[3]).toBeCloseTo(-Math.PI / 2, 10);
    });

    it('hypot computes sqrt(x^2 + y^2)', async () => {
      const a = track(await NDArray.fromArray([3, 5, 8]));
      const b = track(await NDArray.fromArray([4, 12, 15]));
      const result = track(hypot(a, b));
      expect(result.toArray()).toEqual([5, 13, 17]);
    });
  });

  describe('Logical Operations', () => {
    it('logical_not inverts boolean values', async () => {
      const arr = track(await NDArray.fromArray([0, 1, 0, 5, -3]));
      const result = track(logical_not(arr));
      expect(result.toArray()).toEqual([1, 0, 1, 0, 0]);
    });

    it('logical_and performs element-wise AND', async () => {
      const a = track(await NDArray.fromArray([1, 1, 0, 0]));
      const b = track(await NDArray.fromArray([1, 0, 1, 0]));
      const result = track(logical_and(a, b));
      expect(result.toArray()).toEqual([1, 0, 0, 0]);
    });

    it('logical_or performs element-wise OR', async () => {
      const a = track(await NDArray.fromArray([1, 1, 0, 0]));
      const b = track(await NDArray.fromArray([1, 0, 1, 0]));
      const result = track(logical_or(a, b));
      expect(result.toArray()).toEqual([1, 1, 1, 0]);
    });
  });

  describe('Non-contiguous Arrays', () => {
    it('sqrt works on transposed arrays', async () => {
      const arr = track(await NDArray.fromArray([[1, 4], [9, 16]]));
      const transposed = track(arr.T);
      const result = track(sqrt(transposed));
      expect(result.shape).toEqual([2, 2]);
      expect(result.toArray()).toEqual([1, 3, 2, 4]);
    });

    it('add works with non-contiguous views', async () => {
      const arr = track(await NDArray.fromArray([[1, 2, 3], [4, 5, 6]]));
      const t1 = track(arr.T);
      const t2 = track(arr.T);
      const result = track(add(t1, t2));
      expect(result.shape).toEqual([3, 2]);
      expect(result.toArray()).toEqual([2, 8, 4, 10, 6, 12]);
    });
  });
});
