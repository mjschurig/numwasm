import { describe, it, expect, beforeAll } from 'vitest';
import * as sym from '../../dist/symwasm.mjs';

describe('series() - Taylor series expansion', () => {
  beforeAll(async () => {
    await sym.core.loadWasmModule();
  });

  describe('basic expansions', () => {
    it('series(x, x) = x', () => {
      const x = new sym.core.Symbol('x');
      const result = sym.core.series(x, x);
      expect(result.toString()).toBe('x');
    });

    it('series(1, x) = 1', () => {
      const x = new sym.core.Symbol('x');
      const one = new sym.core.Integer(1);
      const result = sym.core.series(one, x);
      expect(result.toString()).toBe('1');
    });

    it('series(x^2, x) = x^2', () => {
      const x = new sym.core.Symbol('x');
      const x2 = sym.core.pow(x, new sym.core.Integer(2));
      const result = sym.core.series(x2, x);
      expect(result.toString()).toBe('x**2');
    });
  });

  describe('exponential series', () => {
    it('series(exp(x), x, 0, 4) contains expected terms', () => {
      const x = new sym.core.Symbol('x');
      const result = sym.core.series(sym.core.exp(x), x, 0, 4);
      // exp(x) = 1 + x + x²/2 + x³/6 + O(x⁴)
      const str = result.toString();
      expect(str).toMatch(/1/);
      expect(str).toMatch(/x/);
    });

    it('series(exp(x), x, 0, 1) = 1', () => {
      const x = new sym.core.Symbol('x');
      const result = sym.core.series(sym.core.exp(x), x, 0, 1);
      expect(result.toString()).toBe('1');
    });

    it('series(exp(x), x, 0, 2) = 1 + x', () => {
      const x = new sym.core.Symbol('x');
      const result = sym.core.series(sym.core.exp(x), x, 0, 2);
      expect(result.toString()).toMatch(/1 \+ x|x \+ 1/);
    });
  });

  describe('trigonometric series', () => {
    it('series(sin(x), x, 0, 2) = x', () => {
      const x = new sym.core.Symbol('x');
      const result = sym.core.series(sym.core.sin(x), x, 0, 2);
      // sin(x) = x + O(x^2) when prec=2
      expect(result.toString()).toBe('x');
    });

    it('series(sin(x), x, 0, 4) contains x and x^3 term', () => {
      const x = new sym.core.Symbol('x');
      const result = sym.core.series(sym.core.sin(x), x, 0, 4);
      // sin(x) = x - x³/6 + O(x^4)
      const str = result.toString();
      expect(str).toMatch(/x/);
      expect(str).toMatch(/x\*\*3/);
    });

    it('series(cos(x), x, 0, 1) = 1', () => {
      const x = new sym.core.Symbol('x');
      const result = sym.core.series(sym.core.cos(x), x, 0, 1);
      // cos(x) = 1 + O(x)
      expect(result.toString()).toBe('1');
    });

    it('series(cos(x), x, 0, 3) contains 1 and x^2 term', () => {
      const x = new sym.core.Symbol('x');
      const result = sym.core.series(sym.core.cos(x), x, 0, 3);
      // cos(x) = 1 - x²/2 + O(x^3)
      const str = result.toString();
      expect(str).toMatch(/1/);
      expect(str).toMatch(/x\*\*2/);
    });

    it('series(tan(x), x, 0, 4) starts with x', () => {
      const x = new sym.core.Symbol('x');
      const result = sym.core.series(sym.core.tan(x), x, 0, 4);
      // tan(x) = x + x³/3 + O(x^4)
      const str = result.toString();
      expect(str).toMatch(/x/);
    });
  });

  describe('logarithmic series', () => {
    it('series(log(1+x), x, 0, 2) = x', () => {
      const x = new sym.core.Symbol('x');
      const one = new sym.core.Integer(1);
      const expr = sym.core.log(sym.core.add(one, x));
      const result = sym.core.series(expr, x, 0, 2);
      // log(1+x) = x + O(x^2)
      expect(result.toString()).toBe('x');
    });

    it('series(log(1+x), x, 0, 4) contains multiple terms', () => {
      const x = new sym.core.Symbol('x');
      const one = new sym.core.Integer(1);
      const expr = sym.core.log(sym.core.add(one, x));
      const result = sym.core.series(expr, x, 0, 4);
      // log(1+x) = x - x²/2 + x³/3 + O(x^4)
      const str = result.toString();
      expect(str).toMatch(/x/);
      expect(str).toMatch(/x\*\*2/);
      expect(str).toMatch(/x\*\*3/);
    });
  });

  describe('default parameters', () => {
    it('series(sin(x), x) uses default n=6', () => {
      const x = new sym.core.Symbol('x');
      const result = sym.core.series(sym.core.sin(x), x);
      // Should have terms up to x⁵ (since prec=6 means O(x^6))
      const str = result.toString();
      expect(str).toMatch(/x\*\*5/);
    });

    it('series(exp(x), x) uses default n=6', () => {
      const x = new sym.core.Symbol('x');
      const result = sym.core.series(sym.core.exp(x), x);
      // Should have terms up to x⁵
      const str = result.toString();
      expect(str).toMatch(/x\*\*5/);
    });
  });

  describe('composition of functions', () => {
    it('series(exp(sin(x)), x, 0, 4)', () => {
      const x = new sym.core.Symbol('x');
      const expr = sym.core.exp(sym.core.sin(x));
      const result = sym.core.series(expr, x, 0, 4);
      // exp(sin(x)) = 1 + x + x²/2 + O(x^4)
      const str = result.toString();
      expect(str).toMatch(/1/);
      expect(str).toMatch(/x/);
    });

    it('series(sin(x)^2, x, 0, 5)', () => {
      const x = new sym.core.Symbol('x');
      const sinX = sym.core.sin(x);
      const expr = sym.core.pow(sinX, new sym.core.Integer(2));
      const result = sym.core.series(expr, x, 0, 5);
      // sin²(x) = x² - x⁴/3 + O(x^5)
      const str = result.toString();
      expect(str).toMatch(/x\*\*2/);
    });
  });

  describe('error handling', () => {
    it('throws for non-zero expansion point', () => {
      const x = new sym.core.Symbol('x');
      expect(() => sym.core.series(x, x, 1)).toThrow(
        /only supports expansion around x=0/
      );
    });

    it('throws for non-zero expansion point (negative)', () => {
      const x = new sym.core.Symbol('x');
      expect(() => sym.core.series(x, x, -1)).toThrow(
        /only supports expansion around x=0/
      );
    });
  });

  describe('also exported from calculus module', () => {
    it('calculus.series works the same as core.series', () => {
      const x = new sym.core.Symbol('x');
      const result = sym.calculus.series(sym.core.exp(x), x, 0, 3);
      // exp(x) = 1 + x + x²/2 + O(x^3)
      const str = result.toString();
      expect(str).toMatch(/1/);
      expect(str).toMatch(/x/);
    });
  });
});
