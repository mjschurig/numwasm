/**
 * Tests for symwasm number types (Phase 1.4)
 *
 * Tests Integer, Rational, Float classes and S constants
 */
import { describe, it, expect, beforeAll } from 'vitest';
import * as sym from '../../src/ts/index';

describe('numbers: Integer', () => {
  beforeAll(async () => {
    await sym.core.loadWasmModule();
  });

  describe('constructor', () => {
    it('creates Integer with positive value', () => {
      const n = new sym.core.Integer(42);
      expect(n.value).toBe(42);
      expect(n.toString()).toBe('42');
    });

    it('creates Integer with zero', () => {
      const n = new sym.core.Integer(0);
      expect(n.value).toBe(0);
      expect(n.toString()).toBe('0');
    });

    it('creates Integer with negative value', () => {
      const n = new sym.core.Integer(-17);
      expect(n.value).toBe(-17);
      expect(n.toString()).toBe('-17');
    });

    it('truncates non-integer values', () => {
      const n = new sym.core.Integer(3.7);
      expect(n.value).toBe(3);
      expect(n.toString()).toBe('3');
    });
  });

  describe('equals()', () => {
    it('returns true for equal integers', () => {
      const a = new sym.core.Integer(42);
      const b = new sym.core.Integer(42);
      expect(a.equals(b)).toBe(true);
    });

    it('returns false for different integers', () => {
      const a = new sym.core.Integer(42);
      const b = new sym.core.Integer(43);
      expect(a.equals(b)).toBe(false);
    });
  });

  describe('hash()', () => {
    it('returns consistent hash for same value', () => {
      const a = new sym.core.Integer(42);
      const b = new sym.core.Integer(42);
      expect(a.hash()).toBe(b.hash());
    });
  });

  describe('get_type()', () => {
    it('returns SYMENGINE_INTEGER type', () => {
      const n = new sym.core.Integer(42);
      expect(n.get_type()).toBe(sym.core.SymEngineTypeID.SYMENGINE_INTEGER);
    });
  });

  describe('free_symbols()', () => {
    it('returns empty array', () => {
      const n = new sym.core.Integer(42);
      expect(n.free_symbols()).toEqual([]);
    });
  });
});

describe('numbers: Rational', () => {
  beforeAll(async () => {
    await sym.core.loadWasmModule();
  });

  describe('constructor', () => {
    it('creates Rational with positive values', () => {
      const r = new sym.core.Rational(1, 2);
      expect(r.p).toBe(1);
      expect(r.q).toBe(2);
      expect(r.toString()).toBe('1/2');
    });

    it('creates Rational with negative numerator', () => {
      const r = new sym.core.Rational(-3, 4);
      expect(r.p).toBe(-3);
      expect(r.q).toBe(4);
      expect(r.toString()).toBe('-3/4');
    });

    it('defaults denominator to 1', () => {
      const r = new sym.core.Rational(5);
      expect(r.p).toBe(5);
      expect(r.q).toBe(1);
      // SymEngine may return "5" instead of "5/1"
      expect(r.toString()).toBe('5');
    });

    it('throws on zero denominator', () => {
      expect(() => new sym.core.Rational(1, 0)).toThrow('Rational denominator cannot be zero');
    });

    it('simplifies rational automatically (SymEngine behavior)', () => {
      const r = new sym.core.Rational(4, 8);
      // SymEngine simplifies 4/8 to 1/2
      expect(r.toString()).toBe('1/2');
    });
  });

  describe('equals()', () => {
    it('returns true for equal rationals', () => {
      const a = new sym.core.Rational(1, 2);
      const b = new sym.core.Rational(1, 2);
      expect(a.equals(b)).toBe(true);
    });

    it('returns true for equivalent rationals (simplified)', () => {
      const a = new sym.core.Rational(1, 2);
      const b = new sym.core.Rational(2, 4);
      expect(a.equals(b)).toBe(true);
    });

    it('returns false for different rationals', () => {
      const a = new sym.core.Rational(1, 2);
      const b = new sym.core.Rational(1, 3);
      expect(a.equals(b)).toBe(false);
    });
  });

  describe('get_type()', () => {
    it('returns SYMENGINE_RATIONAL type', () => {
      const r = new sym.core.Rational(1, 2);
      expect(r.get_type()).toBe(sym.core.SymEngineTypeID.SYMENGINE_RATIONAL);
    });
  });

  describe('free_symbols()', () => {
    it('returns empty array', () => {
      const r = new sym.core.Rational(1, 2);
      expect(r.free_symbols()).toEqual([]);
    });
  });
});

describe('numbers: Float', () => {
  beforeAll(async () => {
    await sym.core.loadWasmModule();
  });

  describe('constructor', () => {
    it('creates Float with positive value', () => {
      const f = new sym.core.Float(3.14);
      expect(f.value).toBe(3.14);
      // SymEngine may format differently
      expect(f.toString()).toContain('3.14');
    });

    it('creates Float with zero', () => {
      const f = new sym.core.Float(0.0);
      expect(f.value).toBe(0);
      expect(f.toString()).toContain('0');
    });

    it('creates Float with negative value', () => {
      const f = new sym.core.Float(-2.5);
      expect(f.value).toBe(-2.5);
      expect(f.toString()).toContain('-2.5');
    });
  });

  describe('equals()', () => {
    it('returns true for equal floats', () => {
      const a = new sym.core.Float(3.14);
      const b = new sym.core.Float(3.14);
      expect(a.equals(b)).toBe(true);
    });

    it('returns false for different floats', () => {
      const a = new sym.core.Float(3.14);
      const b = new sym.core.Float(3.15);
      expect(a.equals(b)).toBe(false);
    });
  });

  describe('get_type()', () => {
    it('returns SYMENGINE_REAL_DOUBLE type', () => {
      const f = new sym.core.Float(3.14);
      expect(f.get_type()).toBe(sym.core.SymEngineTypeID.SYMENGINE_REAL_DOUBLE);
    });
  });

  describe('free_symbols()', () => {
    it('returns empty array', () => {
      const f = new sym.core.Float(3.14);
      expect(f.free_symbols()).toEqual([]);
    });
  });
});

describe('numbers: S constants', () => {
  beforeAll(async () => {
    await sym.core.loadWasmModule();
  });

  describe('S.Zero', () => {
    it('returns Integer(0)', () => {
      const zero = sym.core.S.Zero;
      expect(zero).toBeInstanceOf(sym.core.Integer);
      expect(zero.toString()).toBe('0');
    });

    it('returns the same instance on repeated access', () => {
      const a = sym.core.S.Zero;
      const b = sym.core.S.Zero;
      expect(a).toBe(b);
    });
  });

  describe('S.One', () => {
    it('returns Integer(1)', () => {
      const one = sym.core.S.One;
      expect(one).toBeInstanceOf(sym.core.Integer);
      expect(one.toString()).toBe('1');
    });

    it('returns the same instance on repeated access', () => {
      const a = sym.core.S.One;
      const b = sym.core.S.One;
      expect(a).toBe(b);
    });
  });

  describe('S.NegativeOne', () => {
    it('returns Integer(-1)', () => {
      const negOne = sym.core.S.NegativeOne;
      expect(negOne).toBeInstanceOf(sym.core.Integer);
      expect(negOne.toString()).toBe('-1');
    });

    it('returns the same instance on repeated access', () => {
      const a = sym.core.S.NegativeOne;
      const b = sym.core.S.NegativeOne;
      expect(a).toBe(b);
    });
  });

  describe('S.Half', () => {
    it('returns Rational(1, 2)', () => {
      const half = sym.core.S.Half;
      expect(half).toBeInstanceOf(sym.core.Rational);
      expect(half.toString()).toBe('1/2');
    });

    it('returns the same instance on repeated access', () => {
      const a = sym.core.S.Half;
      const b = sym.core.S.Half;
      expect(a).toBe(b);
    });
  });
});

describe('numbers: Complex', () => {
  beforeAll(async () => {
    await sym.core.loadWasmModule();
  });

  describe('constructor', () => {
    it('creates Complex with integer parts', () => {
      const c = new sym.core.Complex(1, 2);
      expect(c.toString()).toBe('1 + 2*I');
    });

    it('creates Complex with zero real part', () => {
      const c = new sym.core.Complex(0, 1);
      expect(c.toString()).toBe('I');
    });

    it('creates Complex with zero imaginary part', () => {
      const c = new sym.core.Complex(3, 0);
      // SymEngine may simplify to just "3"
      expect(c.toString()).toBe('3');
    });

    it('creates Complex with negative imaginary part', () => {
      const c = new sym.core.Complex(1, -2);
      expect(c.toString()).toBe('1 - 2*I');
    });

    it('creates Complex with Integer objects', () => {
      const re = new sym.core.Integer(3);
      const im = new sym.core.Integer(4);
      const c = new sym.core.Complex(re, im);
      expect(c.toString()).toBe('3 + 4*I');
    });

    it('creates Complex with Rational parts', () => {
      const re = new sym.core.Rational(1, 2);
      const im = new sym.core.Rational(3, 4);
      const c = new sym.core.Complex(re, im);
      expect(c.toString()).toBe('1/2 + 3/4*I');
    });
  });

  describe('equals()', () => {
    it('returns true for equal complex numbers', () => {
      const a = new sym.core.Complex(1, 2);
      const b = new sym.core.Complex(1, 2);
      expect(a.equals(b)).toBe(true);
    });

    it('returns false for different complex numbers', () => {
      const a = new sym.core.Complex(1, 2);
      const b = new sym.core.Complex(1, 3);
      expect(a.equals(b)).toBe(false);
    });
  });

  describe('get_type()', () => {
    it('returns SYMENGINE_COMPLEX type', () => {
      const c = new sym.core.Complex(1, 2);
      expect(c.get_type()).toBe(sym.core.SymEngineTypeID.SYMENGINE_COMPLEX);
    });
  });

  describe('free_symbols()', () => {
    it('returns empty array', () => {
      const c = new sym.core.Complex(1, 2);
      expect(c.free_symbols()).toEqual([]);
    });
  });
});

describe('numbers: cross-type comparisons', () => {
  beforeAll(async () => {
    await sym.core.loadWasmModule();
  });

  it('Integer and Rational(n, 1) are equal', () => {
    const int = new sym.core.Integer(5);
    const rat = new sym.core.Rational(5, 1);
    expect(int.equals(rat)).toBe(true);
  });

  it('Integer and Symbol are not equal', () => {
    const int = new sym.core.Integer(5);
    const sym_x = new sym.core.Symbol('x');
    expect(int.equals(sym_x)).toBe(false);
  });
});
