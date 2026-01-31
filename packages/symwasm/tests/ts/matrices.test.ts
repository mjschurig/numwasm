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

describe('Matrix Submatrix Operations', () => {
  beforeAll(async () => {
    await sym.core.loadWasmModule();
  });

  describe('submatrix', () => {
    it('extracts submatrix with default step', () => {
      const m = new sym.matrices.Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]]);
      // Indices are inclusive: (0, 0) to (1, 1) gives 2x2
      const sub = m.submatrix(0, 0, 1, 1);
      expect(sub.rows).toBe(2);
      expect(sub.cols).toBe(2);
      expect(sub.get(0, 0).toString()).toBe('1');
      expect(sub.get(0, 1).toString()).toBe('2');
      expect(sub.get(1, 0).toString()).toBe('4');
      expect(sub.get(1, 1).toString()).toBe('5');
    });

    it('extracts single row', () => {
      const m = new sym.matrices.Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]]);
      // Row 1, columns 0 to 2 inclusive
      const sub = m.submatrix(1, 0, 1, 2);
      expect(sub.rows).toBe(1);
      expect(sub.cols).toBe(3);
      expect(sub.get(0, 0).toString()).toBe('4');
      expect(sub.get(0, 2).toString()).toBe('6');
    });

    it('extracts single column', () => {
      const m = new sym.matrices.Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]]);
      // Column 1, rows 0 to 2 inclusive
      const sub = m.submatrix(0, 1, 2, 1);
      expect(sub.rows).toBe(3);
      expect(sub.cols).toBe(1);
      expect(sub.get(0, 0).toString()).toBe('2');
      expect(sub.get(2, 0).toString()).toBe('8');
    });

    it('extracts with step', () => {
      const m = new sym.matrices.Matrix([[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12], [13, 14, 15, 16]]);
      // With step 2: only every other element is copied
      // Output still has dimensions (r2-r1+1) x (c2-c1+1) = 4x4
      // But only positions 0,2 x 0,2 get copied, others are uninitialized/zero
      const sub = m.submatrix(0, 0, 3, 3, 2, 2);
      expect(sub.rows).toBe(4);
      expect(sub.cols).toBe(4);
      expect(sub.get(0, 0).toString()).toBe('1');
      expect(sub.get(0, 2).toString()).toBe('3');
      expect(sub.get(2, 0).toString()).toBe('9');
      expect(sub.get(2, 2).toString()).toBe('11');
    });

    it('works with symbolic elements', () => {
      const [x, y] = sym.core.symbols('x y');
      const m = new sym.matrices.Matrix([[x, 1], [2, y]]);
      // Column 0, all rows
      const sub = m.submatrix(0, 0, 1, 0);
      expect(sub.get(0, 0).toString()).toBe('x');
      expect(sub.get(1, 0).toString()).toBe('2');
    });

    it('throws RangeError on out-of-bounds row', () => {
      const m = new sym.matrices.Matrix([[1, 2], [3, 4]]);
      expect(() => m.submatrix(0, 0, 2, 1)).toThrow(RangeError);
      expect(() => m.submatrix(-1, 0, 1, 1)).toThrow(RangeError);
    });

    it('throws RangeError on out-of-bounds column', () => {
      const m = new sym.matrices.Matrix([[1, 2], [3, 4]]);
      expect(() => m.submatrix(0, 0, 1, 2)).toThrow(RangeError);
      expect(() => m.submatrix(0, -1, 1, 1)).toThrow(RangeError);
    });

    it('throws on non-positive step', () => {
      const m = new sym.matrices.Matrix([[1, 2], [3, 4]]);
      expect(() => m.submatrix(0, 0, 1, 1, 0, 1)).toThrow();
      expect(() => m.submatrix(0, 0, 1, 1, 1, 0)).toThrow();
    });
  });

  describe('rowJoin', () => {
    it('joins matrices horizontally', () => {
      const a = new sym.matrices.Matrix([[1, 2], [3, 4]]);
      const b = new sym.matrices.Matrix([[5], [6]]);
      const result = a.rowJoin(b);
      expect(result.rows).toBe(2);
      expect(result.cols).toBe(3);
      expect(result.get(0, 0).toString()).toBe('1');
      expect(result.get(0, 2).toString()).toBe('5');
      expect(result.get(1, 2).toString()).toBe('6');
    });

    it('joins multiple columns', () => {
      const a = new sym.matrices.Matrix([[1], [2]]);
      const b = new sym.matrices.Matrix([[3, 4], [5, 6]]);
      const result = a.rowJoin(b);
      expect(result.cols).toBe(3);
      expect(result.get(0, 1).toString()).toBe('3');
      expect(result.get(1, 2).toString()).toBe('6');
    });

    it('works with symbolic elements', () => {
      const x = new sym.core.Symbol('x');
      const a = new sym.matrices.Matrix([[1, 2]]);
      const b = new sym.matrices.Matrix([[x]]);
      const result = a.rowJoin(b);
      expect(result.get(0, 2).toString()).toBe('x');
    });

    it('throws on mismatched row count', () => {
      const a = new sym.matrices.Matrix([[1, 2], [3, 4]]);
      const b = new sym.matrices.Matrix([[5, 6, 7]]);
      expect(() => a.rowJoin(b)).toThrow();
    });
  });

  describe('colJoin', () => {
    it('joins matrices vertically', () => {
      const a = new sym.matrices.Matrix([[1, 2], [3, 4]]);
      const b = new sym.matrices.Matrix([[5, 6]]);
      const result = a.colJoin(b);
      expect(result.rows).toBe(3);
      expect(result.cols).toBe(2);
      expect(result.get(2, 0).toString()).toBe('5');
      expect(result.get(2, 1).toString()).toBe('6');
    });

    it('joins multiple rows', () => {
      const a = new sym.matrices.Matrix([[1, 2]]);
      const b = new sym.matrices.Matrix([[3, 4], [5, 6]]);
      const result = a.colJoin(b);
      expect(result.rows).toBe(3);
      expect(result.get(1, 0).toString()).toBe('3');
      expect(result.get(2, 1).toString()).toBe('6');
    });

    it('works with symbolic elements', () => {
      const y = new sym.core.Symbol('y');
      const a = new sym.matrices.Matrix([[1, 2]]);
      const b = new sym.matrices.Matrix([[y, 3]]);
      const result = a.colJoin(b);
      expect(result.get(1, 0).toString()).toBe('y');
    });

    it('throws on mismatched column count', () => {
      const a = new sym.matrices.Matrix([[1, 2], [3, 4]]);
      const b = new sym.matrices.Matrix([[5, 6, 7]]);
      expect(() => a.colJoin(b)).toThrow();
    });
  });

  describe('rowDel', () => {
    it('deletes a row in-place', () => {
      const m = new sym.matrices.Matrix([[1, 2], [3, 4], [5, 6]]);
      m.rowDel(1);
      expect(m.rows).toBe(2);
      expect(m.cols).toBe(2);
      expect(m.get(0, 0).toString()).toBe('1');
      expect(m.get(1, 0).toString()).toBe('5');
    });

    it('deletes first row', () => {
      const m = new sym.matrices.Matrix([[1, 2], [3, 4]]);
      m.rowDel(0);
      expect(m.rows).toBe(1);
      expect(m.get(0, 0).toString()).toBe('3');
      expect(m.get(0, 1).toString()).toBe('4');
    });

    it('deletes last row', () => {
      const m = new sym.matrices.Matrix([[1, 2], [3, 4]]);
      m.rowDel(1);
      expect(m.rows).toBe(1);
      expect(m.get(0, 0).toString()).toBe('1');
      expect(m.get(0, 1).toString()).toBe('2');
    });

    it('throws RangeError on out-of-bounds index', () => {
      const m = new sym.matrices.Matrix([[1, 2], [3, 4]]);
      expect(() => m.rowDel(2)).toThrow(RangeError);
      expect(() => m.rowDel(-1)).toThrow(RangeError);
    });
  });

  describe('colDel', () => {
    it('deletes a column in-place', () => {
      const m = new sym.matrices.Matrix([[1, 2, 3], [4, 5, 6]]);
      m.colDel(1);
      expect(m.rows).toBe(2);
      expect(m.cols).toBe(2);
      expect(m.get(0, 0).toString()).toBe('1');
      expect(m.get(0, 1).toString()).toBe('3');
      expect(m.get(1, 0).toString()).toBe('4');
      expect(m.get(1, 1).toString()).toBe('6');
    });

    it('deletes first column', () => {
      const m = new sym.matrices.Matrix([[1, 2], [3, 4]]);
      m.colDel(0);
      expect(m.cols).toBe(1);
      expect(m.get(0, 0).toString()).toBe('2');
      expect(m.get(1, 0).toString()).toBe('4');
    });

    it('deletes last column', () => {
      const m = new sym.matrices.Matrix([[1, 2], [3, 4]]);
      m.colDel(1);
      expect(m.cols).toBe(1);
      expect(m.get(0, 0).toString()).toBe('1');
      expect(m.get(1, 0).toString()).toBe('3');
    });

    it('throws RangeError on out-of-bounds index', () => {
      const m = new sym.matrices.Matrix([[1, 2], [3, 4]]);
      expect(() => m.colDel(2)).toThrow(RangeError);
      expect(() => m.colDel(-1)).toThrow(RangeError);
    });
  });
});

describe('Matrix Basic Operations', () => {
  beforeAll(async () => {
    await sym.core.loadWasmModule();
  });

  describe('det()', () => {
    it('computes determinant of 2x2 matrix', () => {
      const m = new sym.matrices.Matrix([[1, 2], [3, 4]]);
      const d = m.det();
      // det = 1*4 - 2*3 = -2
      expect(d.toString()).toBe('-2');
    });

    it('computes determinant of 3x3 matrix', () => {
      const m = new sym.matrices.Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 10]]);
      const d = m.det();
      // det = 1*(5*10 - 6*8) - 2*(4*10 - 6*7) + 3*(4*8 - 5*7) = 1*(50-48) - 2*(40-42) + 3*(32-35)
      //     = 1*2 - 2*(-2) + 3*(-3) = 2 + 4 - 9 = -3
      expect(d.toString()).toBe('-3');
    });

    it('computes determinant with symbolic entries', () => {
      const x = new sym.core.Symbol('x');
      const m = new sym.matrices.Matrix([[x, 1], [0, x]]);
      const d = m.det();
      // det = x*x - 1*0 = x^2
      expect(d.toString()).toBe('x**2');
    });

    it('returns 0 for singular matrix', () => {
      const m = new sym.matrices.Matrix([[1, 2], [2, 4]]);
      const d = m.det();
      expect(d.toString()).toBe('0');
    });

    it('throws for non-square matrix', () => {
      const m = new sym.matrices.Matrix([[1, 2, 3], [4, 5, 6]]);
      expect(() => m.det()).toThrow('square matrix');
    });
  });

  describe('inv()', () => {
    it('computes inverse of 2x2 matrix', () => {
      const m = new sym.matrices.Matrix([[4, 7], [2, 6]]);
      const inv = m.inv();
      // det = 4*6 - 7*2 = 10
      // inv = 1/10 * [[6, -7], [-2, 4]] = [[3/5, -7/10], [-1/5, 2/5]]
      expect(inv.rows).toBe(2);
      expect(inv.cols).toBe(2);
      // Verify M * M^(-1) = I
      const product = m.mul(inv);
      expect(product.get(0, 0).toString()).toBe('1');
      expect(product.get(0, 1).toString()).toBe('0');
      expect(product.get(1, 0).toString()).toBe('0');
      expect(product.get(1, 1).toString()).toBe('1');
    });

    it('computes inverse of identity matrix', () => {
      const m = sym.matrices.eye(3);
      const inv = m.inv();
      expect(inv.equals(m)).toBe(true);
    });

    it('computes inverse with symbolic entries', () => {
      const x = new sym.core.Symbol('x');
      const m = new sym.matrices.Matrix([[x, 0], [0, 1]]);
      const inv = m.inv();
      // inv = [[1/x, 0], [0, 1]]
      // Note: SymEngine may represent 1/x as x**(-1)
      const invStr = inv.get(0, 0).toString();
      expect(invStr === '1/x' || invStr === 'x**(-1)').toBe(true);
      expect(inv.get(0, 1).toString()).toBe('0');
      expect(inv.get(1, 0).toString()).toBe('0');
      expect(inv.get(1, 1).toString()).toBe('1');
    });

    it('throws for non-square matrix', () => {
      const m = new sym.matrices.Matrix([[1, 2, 3]]);
      expect(() => m.inv()).toThrow('square matrix');
    });
  });

  describe('transpose()', () => {
    it('transposes square matrix', () => {
      const m = new sym.matrices.Matrix([[1, 2], [3, 4]]);
      const t = m.transpose();
      expect(t.rows).toBe(2);
      expect(t.cols).toBe(2);
      expect(t.get(0, 0).toString()).toBe('1');
      expect(t.get(0, 1).toString()).toBe('3');
      expect(t.get(1, 0).toString()).toBe('2');
      expect(t.get(1, 1).toString()).toBe('4');
    });

    it('transposes rectangular matrix', () => {
      const m = new sym.matrices.Matrix([[1, 2, 3], [4, 5, 6]]);
      const t = m.transpose();
      expect(t.rows).toBe(3);
      expect(t.cols).toBe(2);
      expect(t.get(0, 0).toString()).toBe('1');
      expect(t.get(0, 1).toString()).toBe('4');
      expect(t.get(1, 0).toString()).toBe('2');
      expect(t.get(1, 1).toString()).toBe('5');
      expect(t.get(2, 0).toString()).toBe('3');
      expect(t.get(2, 1).toString()).toBe('6');
    });

    it('double transpose equals original', () => {
      const m = new sym.matrices.Matrix([[1, 2], [3, 4]]);
      const t = m.transpose().transpose();
      expect(t.equals(m)).toBe(true);
    });
  });

  describe('add()', () => {
    it('adds two matrices', () => {
      const a = new sym.matrices.Matrix([[1, 2], [3, 4]]);
      const b = new sym.matrices.Matrix([[5, 6], [7, 8]]);
      const sum = a.add(b);
      expect(sum.get(0, 0).toString()).toBe('6');
      expect(sum.get(0, 1).toString()).toBe('8');
      expect(sum.get(1, 0).toString()).toBe('10');
      expect(sum.get(1, 1).toString()).toBe('12');
    });

    it('adds with symbolic entries', () => {
      const x = new sym.core.Symbol('x');
      const a = new sym.matrices.Matrix([[x, 1]]);
      const b = new sym.matrices.Matrix([[1, x]]);
      const sum = a.add(b);
      expect(sum.get(0, 0).toString()).toBe('1 + x');
      expect(sum.get(0, 1).toString()).toBe('1 + x');
    });

    it('throws for mismatched dimensions', () => {
      const a = new sym.matrices.Matrix([[1, 2]]);
      const b = new sym.matrices.Matrix([[1], [2]]);
      expect(() => a.add(b)).toThrow('dimensions must match');
    });
  });

  describe('mul()', () => {
    it('multiplies two 2x2 matrices', () => {
      const a = new sym.matrices.Matrix([[1, 2], [3, 4]]);
      const b = new sym.matrices.Matrix([[5, 6], [7, 8]]);
      const prod = a.mul(b);
      // [[1*5+2*7, 1*6+2*8], [3*5+4*7, 3*6+4*8]] = [[19, 22], [43, 50]]
      expect(prod.get(0, 0).toString()).toBe('19');
      expect(prod.get(0, 1).toString()).toBe('22');
      expect(prod.get(1, 0).toString()).toBe('43');
      expect(prod.get(1, 1).toString()).toBe('50');
    });

    it('multiplies rectangular matrices', () => {
      const a = new sym.matrices.Matrix([[1, 2, 3]]);  // 1x3
      const b = new sym.matrices.Matrix([[4], [5], [6]]);  // 3x1
      const prod = a.mul(b);
      // Result is 1x1: [[1*4 + 2*5 + 3*6]] = [[32]]
      expect(prod.rows).toBe(1);
      expect(prod.cols).toBe(1);
      expect(prod.get(0, 0).toString()).toBe('32');
    });

    it('multiplies with symbolic entries', () => {
      const x = new sym.core.Symbol('x');
      const a = new sym.matrices.Matrix([[x, 1], [0, x]]);
      const b = new sym.matrices.Matrix([[1, 0], [0, 1]]);
      const prod = a.mul(b);
      expect(prod.equals(a)).toBe(true);
    });

    it('throws for incompatible dimensions', () => {
      const a = new sym.matrices.Matrix([[1, 2]]);  // 1x2
      const b = new sym.matrices.Matrix([[1, 2]]);  // 1x2
      expect(() => a.mul(b)).toThrow('incompatible');
    });
  });

  describe('addScalar()', () => {
    it('adds scalar to all elements', () => {
      const m = new sym.matrices.Matrix([[1, 2], [3, 4]]);
      const result = m.addScalar(10);
      expect(result.get(0, 0).toString()).toBe('11');
      expect(result.get(0, 1).toString()).toBe('12');
      expect(result.get(1, 0).toString()).toBe('13');
      expect(result.get(1, 1).toString()).toBe('14');
    });

    it('adds symbolic scalar', () => {
      const x = new sym.core.Symbol('x');
      const m = new sym.matrices.Matrix([[1, 2]]);
      const result = m.addScalar(x);
      expect(result.get(0, 0).toString()).toBe('1 + x');
      expect(result.get(0, 1).toString()).toBe('2 + x');
    });

    it('preserves original matrix', () => {
      const m = new sym.matrices.Matrix([[1, 2]]);
      m.addScalar(10);
      expect(m.get(0, 0).toString()).toBe('1');
    });
  });

  describe('mulScalar()', () => {
    it('multiplies all elements by scalar', () => {
      const m = new sym.matrices.Matrix([[1, 2], [3, 4]]);
      const result = m.mulScalar(2);
      expect(result.get(0, 0).toString()).toBe('2');
      expect(result.get(0, 1).toString()).toBe('4');
      expect(result.get(1, 0).toString()).toBe('6');
      expect(result.get(1, 1).toString()).toBe('8');
    });

    it('multiplies by symbolic scalar', () => {
      const x = new sym.core.Symbol('x');
      const m = new sym.matrices.Matrix([[1, 2]]);
      const result = m.mulScalar(x);
      expect(result.get(0, 0).toString()).toBe('x');
      expect(result.get(0, 1).toString()).toBe('2*x');
    });

    it('multiplies by zero', () => {
      const m = new sym.matrices.Matrix([[1, 2], [3, 4]]);
      const result = m.mulScalar(0);
      expect(result.get(0, 0).toString()).toBe('0');
      expect(result.get(1, 1).toString()).toBe('0');
    });

    it('preserves original matrix', () => {
      const m = new sym.matrices.Matrix([[1, 2]]);
      m.mulScalar(10);
      expect(m.get(0, 0).toString()).toBe('1');
    });
  });
});
