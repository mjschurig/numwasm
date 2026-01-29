import { describe, it, expect, beforeAll } from 'vitest';
import * as sym from '../../dist/symwasm.mjs';

describe('diff() - symbolic differentiation', () => {
  beforeAll(async () => {
    await sym.core.loadWasmModule();
  });

  describe('polynomial derivatives', () => {
    it('diff(x, x) = 1', () => {
      const x = new sym.core.Symbol('x');
      const result = sym.core.diff(x, x);
      expect(result.toString()).toBe('1');
    });

    it('diff(x^2, x) = 2*x', () => {
      const x = new sym.core.Symbol('x');
      const x2 = sym.core.pow(x, new sym.core.Integer(2));
      const result = sym.core.diff(x2, x);
      expect(result.toString()).toBe('2*x');
    });

    it('diff(x^3, x) = 3*x^2', () => {
      const x = new sym.core.Symbol('x');
      const x3 = sym.core.pow(x, new sym.core.Integer(3));
      const result = sym.core.diff(x3, x);
      expect(result.toString()).toBe('3*x**2');
    });

    it('diff(5, x) = 0 (constant)', () => {
      const x = new sym.core.Symbol('x');
      const five = new sym.core.Integer(5);
      const result = sym.core.diff(five, x);
      expect(result.toString()).toBe('0');
    });

    it('diff(y, x) = 0 (different symbol)', () => {
      const x = new sym.core.Symbol('x');
      const y = new sym.core.Symbol('y');
      const result = sym.core.diff(y, x);
      expect(result.toString()).toBe('0');
    });

    it('diff(3*x + 2, x) = 3', () => {
      const x = new sym.core.Symbol('x');
      const three = new sym.core.Integer(3);
      const two = new sym.core.Integer(2);
      const expr = sym.core.add(sym.core.mul(three, x), two);
      const result = sym.core.diff(expr, x);
      expect(result.toString()).toBe('3');
    });
  });

  describe('trigonometric derivatives', () => {
    it('diff(sin(x), x) = cos(x)', () => {
      const x = new sym.core.Symbol('x');
      const sinX = sym.core.sin(x);
      const result = sym.core.diff(sinX, x);
      expect(result.toString()).toBe('cos(x)');
    });

    it('diff(cos(x), x) = -sin(x)', () => {
      const x = new sym.core.Symbol('x');
      const cosX = sym.core.cos(x);
      const result = sym.core.diff(cosX, x);
      expect(result.toString()).toBe('-sin(x)');
    });

    it('diff(tan(x), x) = sec(x)^2 or 1 + tan(x)^2', () => {
      const x = new sym.core.Symbol('x');
      const tanX = sym.core.tan(x);
      const result = sym.core.diff(tanX, x);
      // SymEngine may represent as 1 + tan(x)^2 or sec(x)^2
      expect(result.toString()).toMatch(/tan|sec|cos/);
    });
  });

  describe('exponential and logarithmic derivatives', () => {
    it('diff(exp(x), x) = exp(x)', () => {
      const x = new sym.core.Symbol('x');
      const expX = sym.core.exp(x);
      const result = sym.core.diff(expX, x);
      expect(result.toString()).toBe('exp(x)');
    });

    it('diff(log(x), x) = 1/x', () => {
      const x = new sym.core.Symbol('x');
      const logX = sym.core.log(x);
      const result = sym.core.diff(logX, x);
      expect(result.toString()).toMatch(/1\/x|x\*\*\(-1\)/);
    });
  });

  describe('chain rule', () => {
    it('diff(sin(x^2), x) = 2*x*cos(x^2)', () => {
      const x = new sym.core.Symbol('x');
      const x2 = sym.core.pow(x, new sym.core.Integer(2));
      const sinX2 = sym.core.sin(x2);
      const result = sym.core.diff(sinX2, x);
      expect(result.toString()).toMatch(/2\*x\*cos\(x\*\*2\)|cos\(x\*\*2\)\*2\*x/);
    });

    it('diff(exp(2*x), x) = 2*exp(2*x)', () => {
      const x = new sym.core.Symbol('x');
      const two = new sym.core.Integer(2);
      const twoX = sym.core.mul(two, x);
      const exp2x = sym.core.exp(twoX);
      const result = sym.core.diff(exp2x, x);
      expect(result.toString()).toMatch(/2\*exp\(2\*x\)|exp\(2\*x\)\*2/);
    });
  });

  describe('product rule', () => {
    it('diff(x*sin(x), x) = sin(x) + x*cos(x)', () => {
      const x = new sym.core.Symbol('x');
      const sinX = sym.core.sin(x);
      const xSinX = sym.core.mul(x, sinX);
      const result = sym.core.diff(xSinX, x);
      // Result should contain both sin(x) and x*cos(x) terms
      expect(result.toString()).toMatch(/sin\(x\)/);
      expect(result.toString()).toMatch(/cos\(x\)/);
    });
  });

  describe('quotient rule', () => {
    it('diff(1/x, x) = -1/x^2', () => {
      const x = new sym.core.Symbol('x');
      const one = new sym.core.Integer(1);
      const oneOverX = sym.core.div(one, x);
      const result = sym.core.diff(oneOverX, x);
      // SymEngine may represent as -x^(-2) or -1/x^2
      expect(result.toString()).toMatch(/-x\*\*\(-2\)|-1\/x\*\*2/);
    });
  });

  describe('higher-order derivatives', () => {
    it('diff(x^3, x, 2) = 6*x (second derivative)', () => {
      const x = new sym.core.Symbol('x');
      const x3 = sym.core.pow(x, new sym.core.Integer(3));
      const result = sym.core.diff(x3, x, 2);
      expect(result.toString()).toBe('6*x');
    });

    it('diff(x^4, x, 3) = 24*x (third derivative)', () => {
      const x = new sym.core.Symbol('x');
      const x4 = sym.core.pow(x, new sym.core.Integer(4));
      const result = sym.core.diff(x4, x, 3);
      expect(result.toString()).toBe('24*x');
    });

    it('diff(sin(x), x, 2) = -sin(x)', () => {
      const x = new sym.core.Symbol('x');
      const sinX = sym.core.sin(x);
      const result = sym.core.diff(sinX, x, 2);
      expect(result.toString()).toBe('-sin(x)');
    });

    it('diff(sin(x), x, 4) = sin(x)', () => {
      const x = new sym.core.Symbol('x');
      const sinX = sym.core.sin(x);
      const result = sym.core.diff(sinX, x, 4);
      expect(result.toString()).toBe('sin(x)');
    });

    it('diff(x^2, x, 0) = x^2 (zeroth derivative)', () => {
      const x = new sym.core.Symbol('x');
      const x2 = sym.core.pow(x, new sym.core.Integer(2));
      const result = sym.core.diff(x2, x, 0);
      expect(result.toString()).toBe('x**2');
    });

    it('diff(x^5, x, 5) = 120 (fifth derivative)', () => {
      const x = new sym.core.Symbol('x');
      const x5 = sym.core.pow(x, new sym.core.Integer(5));
      const result = sym.core.diff(x5, x, 5);
      expect(result.toString()).toBe('120');
    });
  });

  describe('error handling', () => {
    it('throws for negative derivative order', () => {
      const x = new sym.core.Symbol('x');
      expect(() => sym.core.diff(x, x, -1)).toThrow(/non-negative/);
    });
  });

  describe('multi-variable derivatives', () => {
    it('diff(x*y, x) = y', () => {
      const x = new sym.core.Symbol('x');
      const y = new sym.core.Symbol('y');
      const xy = sym.core.mul(x, y);
      const result = sym.core.diff(xy, x);
      expect(result.toString()).toBe('y');
    });

    it('diff(x*y, y) = x', () => {
      const x = new sym.core.Symbol('x');
      const y = new sym.core.Symbol('y');
      const xy = sym.core.mul(x, y);
      const result = sym.core.diff(xy, y);
      expect(result.toString()).toBe('x');
    });

    it('diff(x*y, x, y) = 1 (mixed partial)', () => {
      const x = new sym.core.Symbol('x');
      const y = new sym.core.Symbol('y');
      const xy = sym.core.mul(x, y);
      const result = sym.core.diff(xy, x, y);
      expect(result.toString()).toBe('1');
    });

    it('diff(x^2*y^3, x, 2, y, 3) = 12 (higher mixed)', () => {
      const x = new sym.core.Symbol('x');
      const y = new sym.core.Symbol('y');
      const two = new sym.core.Integer(2);
      const three = new sym.core.Integer(3);
      const expr = sym.core.mul(
        sym.core.pow(x, two),
        sym.core.pow(y, three)
      );
      const result = sym.core.diff(expr, x, 2, y, 3);
      expect(result.toString()).toBe('12');
    });

    it('diff(sin(x)*cos(y), x, y) = -cos(x)*sin(y)', () => {
      const x = new sym.core.Symbol('x');
      const y = new sym.core.Symbol('y');
      const expr = sym.core.mul(sym.core.sin(x), sym.core.cos(y));
      const result = sym.core.diff(expr, x, y);
      // d/dx[sin(x)*cos(y)] = cos(x)*cos(y)
      // d/dy[cos(x)*cos(y)] = -cos(x)*sin(y)
      expect(result.toString()).toMatch(/-sin\(y\)\*cos\(x\)|-cos\(x\)\*sin\(y\)/);
    });

    it('order independence: diff(f, x, y) = diff(f, y, x) for smooth functions', () => {
      const x = new sym.core.Symbol('x');
      const y = new sym.core.Symbol('y');
      // f = x*y^2
      const expr = sym.core.mul(x, sym.core.pow(y, new sym.core.Integer(2)));
      const dxdy = sym.core.diff(expr, x, y);
      const dydx = sym.core.diff(expr, y, x);
      // Both should equal 2*y
      expect(dxdy.equals(dydx)).toBe(true);
    });

    it('diff(x^2 + y^2 + z^2, x, y, z) = 0', () => {
      const x = new sym.core.Symbol('x');
      const y = new sym.core.Symbol('y');
      const z = new sym.core.Symbol('z');
      const two = new sym.core.Integer(2);
      const expr = sym.core.add(
        sym.core.add(sym.core.pow(x, two), sym.core.pow(y, two)),
        sym.core.pow(z, two)
      );
      const result = sym.core.diff(expr, x, y, z);
      expect(result.toString()).toBe('0');
    });

    it('diff(x*y*z, x, y, z) = 1', () => {
      const x = new sym.core.Symbol('x');
      const y = new sym.core.Symbol('y');
      const z = new sym.core.Symbol('z');
      const expr = sym.core.mul(sym.core.mul(x, y), z);
      const result = sym.core.diff(expr, x, y, z);
      expect(result.toString()).toBe('1');
    });
  });

  describe('error handling for multi-variable', () => {
    it('throws when no symbol argument provided', () => {
      const x = new sym.core.Symbol('x');
      expect(() => (sym.core.diff as any)(x)).toThrow(/at least one symbol/);
    });

    it('throws when number comes first without preceding symbol', () => {
      const x = new sym.core.Symbol('x');
      expect(() => (sym.core.diff as any)(x, 2)).toThrow(/must be a symbol/);
    });
  });

  describe('also exported from calculus module', () => {
    it('calculus.diff works the same as core.diff', () => {
      const x = new sym.core.Symbol('x');
      const x2 = sym.core.pow(x, new sym.core.Integer(2));
      const result = sym.calculus.diff(x2, x);
      expect(result.toString()).toBe('2*x');
    });
  });
});
