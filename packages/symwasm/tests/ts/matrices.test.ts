import { describe, it, expect, beforeAll } from 'vitest';
import * as sym from '../../src/ts/index';

describe('Matrix Construction', () => {
  beforeAll(async () => {
    await sym.core.loadWasmModule();
  });

  describe('Matrix constructor', () => {
    it('creates matrix from 2D array of numbers', () => {
      const m = new sym.matrices.Matrix([[1, 2], [3, 4]]);
      expect(m.rows).toBe(2);
      expect(m.cols).toBe(2);
      expect(m.shape).toEqual([2, 2]);
    });

    it('creates matrix from 2D array of Exprs', () => {
      const [x, y] = sym.core.symbols('x y');
      const m = new sym.matrices.Matrix([[x, y], [new sym.core.Integer(1), new sym.core.Integer(2)]]);
      expect(m.rows).toBe(2);
      expect(m.cols).toBe(2);
    });

    it('creates non-square matrix', () => {
      const m = new sym.matrices.Matrix([[1, 2, 3], [4, 5, 6]]);
      expect(m.rows).toBe(2);
      expect(m.cols).toBe(3);
      expect(m.shape).toEqual([2, 3]);
    });

    it('throws on empty data', () => {
      expect(() => new sym.matrices.Matrix([])).toThrow();
    });

    it('throws on empty row', () => {
      expect(() => new sym.matrices.Matrix([[]])).toThrow();
    });

    it('throws on non-rectangular data', () => {
      expect(() => new sym.matrices.Matrix([[1, 2], [3]])).toThrow();
    });
  });

  describe('Matrix.fromFlat', () => {
    it('creates matrix from flat array', () => {
      const m = sym.matrices.Matrix.fromFlat([1, 2, 3, 4, 5, 6], 2, 3);
      expect(m.rows).toBe(2);
      expect(m.cols).toBe(3);
      expect(m.get(0, 0).toString()).toBe('1');
      expect(m.get(0, 2).toString()).toBe('3');
      expect(m.get(1, 0).toString()).toBe('4');
      expect(m.get(1, 2).toString()).toBe('6');
    });

    it('throws on mismatched length', () => {
      expect(() => sym.matrices.Matrix.fromFlat([1, 2, 3], 2, 2)).toThrow();
    });
  });
});

describe('Matrix Properties', () => {
  beforeAll(async () => {
    await sym.core.loadWasmModule();
  });

  describe('get/set', () => {
    it('gets elements correctly', () => {
      const m = new sym.matrices.Matrix([[1, 2], [3, 4]]);
      expect(m.get(0, 0).toString()).toBe('1');
      expect(m.get(0, 1).toString()).toBe('2');
      expect(m.get(1, 0).toString()).toBe('3');
      expect(m.get(1, 1).toString()).toBe('4');
    });

    it('gets symbolic elements', () => {
      const [x, y] = sym.core.symbols('x y');
      const m = new sym.matrices.Matrix([[x, y]]);
      expect(m.get(0, 0).toString()).toBe('x');
      expect(m.get(0, 1).toString()).toBe('y');
    });

    it('sets elements correctly', () => {
      const m = sym.matrices.zeros(2, 2);
      m.set(0, 0, 5);
      expect(m.get(0, 0).toString()).toBe('5');
      expect(m.get(0, 1).toString()).toBe('0');
    });

    it('sets symbolic elements', () => {
      const m = sym.matrices.zeros(2, 2);
      const x = new sym.core.Symbol('x');
      m.set(1, 1, x);
      expect(m.get(1, 1).toString()).toBe('x');
    });

    it('throws RangeError on out-of-bounds get', () => {
      const m = new sym.matrices.Matrix([[1, 2]]);
      expect(() => m.get(1, 0)).toThrow(RangeError);
      expect(() => m.get(0, 2)).toThrow(RangeError);
      expect(() => m.get(-1, 0)).toThrow(RangeError);
    });

    it('throws RangeError on out-of-bounds set', () => {
      const m = new sym.matrices.Matrix([[1, 2]]);
      expect(() => m.set(1, 0, 5)).toThrow(RangeError);
      expect(() => m.set(0, 2, 5)).toThrow(RangeError);
    });
  });

  describe('toString', () => {
    it('returns string representation', () => {
      const m = new sym.matrices.Matrix([[1, 2], [3, 4]]);
      const str = m.toString();
      expect(str).toContain('1');
      expect(str).toContain('2');
      expect(str).toContain('3');
      expect(str).toContain('4');
    });

    it('includes symbolic elements in string', () => {
      const x = new sym.core.Symbol('x');
      const m = new sym.matrices.Matrix([[x, 1]]);
      expect(m.toString()).toContain('x');
    });
  });

  describe('equals', () => {
    it('returns true for equal matrices', () => {
      const m1 = new sym.matrices.Matrix([[1, 2], [3, 4]]);
      const m2 = new sym.matrices.Matrix([[1, 2], [3, 4]]);
      expect(m1.equals(m2)).toBe(true);
    });

    it('returns false for unequal matrices', () => {
      const m1 = new sym.matrices.Matrix([[1, 2], [3, 4]]);
      const m2 = new sym.matrices.Matrix([[1, 2], [3, 5]]);
      expect(m1.equals(m2)).toBe(false);
    });

    it('returns false for different dimensions', () => {
      const m1 = new sym.matrices.Matrix([[1, 2]]);
      const m2 = new sym.matrices.Matrix([[1], [2]]);
      expect(m1.equals(m2)).toBe(false);
    });

    it('returns true for equal symbolic matrices', () => {
      const x = new sym.core.Symbol('x');
      const m1 = new sym.matrices.Matrix([[x, 1]]);
      const m2 = new sym.matrices.Matrix([[x, 1]]);
      expect(m1.equals(m2)).toBe(true);
    });
  });

  describe('free', () => {
    it('frees matrix memory', () => {
      const m = new sym.matrices.Matrix([[1, 2]]);
      m.free();
      // After free, accessing methods should throw
      expect(() => m.rows).toThrow();
    });
  });
});

describe('Matrix Factory Functions', () => {
  beforeAll(async () => {
    await sym.core.loadWasmModule();
  });

  describe('eye', () => {
    it('creates identity matrix', () => {
      const m = sym.matrices.eye(3);
      expect(m.rows).toBe(3);
      expect(m.cols).toBe(3);
      expect(m.get(0, 0).toString()).toBe('1');
      expect(m.get(1, 1).toString()).toBe('1');
      expect(m.get(2, 2).toString()).toBe('1');
      expect(m.get(0, 1).toString()).toBe('0');
      expect(m.get(1, 0).toString()).toBe('0');
    });

    it('creates rectangular identity', () => {
      const m = sym.matrices.eye(2, 3);
      expect(m.rows).toBe(2);
      expect(m.cols).toBe(3);
      expect(m.get(0, 0).toString()).toBe('1');
      expect(m.get(1, 1).toString()).toBe('1');
      expect(m.get(0, 2).toString()).toBe('0');
    });

    it('supports positive diagonal offset', () => {
      const m = sym.matrices.eye(3, 3, 1);
      expect(m.get(0, 0).toString()).toBe('0');
      expect(m.get(0, 1).toString()).toBe('1');
      expect(m.get(1, 2).toString()).toBe('1');
      expect(m.get(1, 1).toString()).toBe('0');
    });

    it('supports negative diagonal offset', () => {
      const m = sym.matrices.eye(3, 3, -1);
      expect(m.get(0, 0).toString()).toBe('0');
      expect(m.get(1, 0).toString()).toBe('1');
      expect(m.get(2, 1).toString()).toBe('1');
    });
  });

  describe('zeros', () => {
    it('creates zero matrix', () => {
      const m = sym.matrices.zeros(2, 3);
      expect(m.rows).toBe(2);
      expect(m.cols).toBe(3);
      for (let i = 0; i < 2; i++) {
        for (let j = 0; j < 3; j++) {
          expect(m.get(i, j).toString()).toBe('0');
        }
      }
    });

    it('creates square zero matrix', () => {
      const m = sym.matrices.zeros(3, 3);
      expect(m.shape).toEqual([3, 3]);
      expect(m.get(1, 1).toString()).toBe('0');
    });
  });

  describe('ones', () => {
    it('creates ones matrix', () => {
      const m = sym.matrices.ones(2, 3);
      expect(m.rows).toBe(2);
      expect(m.cols).toBe(3);
      for (let i = 0; i < 2; i++) {
        for (let j = 0; j < 3; j++) {
          expect(m.get(i, j).toString()).toBe('1');
        }
      }
    });
  });

  describe('diag', () => {
    it('creates diagonal matrix from numbers', () => {
      const m = sym.matrices.diag([1, 2, 3]);
      expect(m.rows).toBe(3);
      expect(m.cols).toBe(3);
      expect(m.get(0, 0).toString()).toBe('1');
      expect(m.get(1, 1).toString()).toBe('2');
      expect(m.get(2, 2).toString()).toBe('3');
      expect(m.get(0, 1).toString()).toBe('0');
      expect(m.get(1, 0).toString()).toBe('0');
    });

    it('supports positive diagonal offset', () => {
      const m = sym.matrices.diag([1, 2], 1);
      expect(m.rows).toBe(3);
      expect(m.cols).toBe(3);
      expect(m.get(0, 1).toString()).toBe('1');
      expect(m.get(1, 2).toString()).toBe('2');
      expect(m.get(0, 0).toString()).toBe('0');
    });

    it('supports negative diagonal offset', () => {
      const m = sym.matrices.diag([1, 2], -1);
      expect(m.rows).toBe(3);
      expect(m.cols).toBe(3);
      expect(m.get(1, 0).toString()).toBe('1');
      expect(m.get(2, 1).toString()).toBe('2');
    });

    it('supports symbolic diagonals', () => {
      const [x, y] = sym.core.symbols('x y');
      const m = sym.matrices.diag([x, y]);
      expect(m.get(0, 0).toString()).toBe('x');
      expect(m.get(1, 1).toString()).toBe('y');
      expect(m.get(0, 1).toString()).toBe('0');
    });
  });
});

describe('Matrix stubs', () => {
  beforeAll(async () => {
    await sym.core.loadWasmModule();
  });

  it('det() throws NotImplementedError', () => {
    const m = new sym.matrices.Matrix([[1, 2], [3, 4]]);
    expect(() => m.det()).toThrow('not yet implemented');
  });

  it('inv() throws NotImplementedError', () => {
    const m = new sym.matrices.Matrix([[1, 2], [3, 4]]);
    expect(() => m.inv()).toThrow('not yet implemented');
  });

  it('transpose() throws NotImplementedError', () => {
    const m = new sym.matrices.Matrix([[1, 2], [3, 4]]);
    expect(() => m.transpose()).toThrow('not yet implemented');
  });
});
