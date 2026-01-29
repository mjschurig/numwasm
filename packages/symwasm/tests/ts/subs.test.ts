/**
 * Tests for symwasm substitution (Phase 1.7)
 *
 * Tests Expr.subs() for single substitution
 */
import { describe, it, expect, beforeAll } from 'vitest';
import * as sym from '../../src/ts/index';

describe('subs: Expr.subs(old, new)', () => {
  beforeAll(async () => {
    await sym.core.loadWasmModule();
  });

  describe('basic substitution', () => {
    it('substitutes symbol with integer', () => {
      const x = new sym.core.Symbol('x');
      const two = new sym.core.Integer(2);
      const result = x.subs(x, two);
      expect(result.toString()).toBe('2');
    });

    it('substitutes symbol with symbol', () => {
      const x = new sym.core.Symbol('x');
      const y = new sym.core.Symbol('y');
      const result = x.subs(x, y);
      expect(result.toString()).toBe('y');
    });

    it('substitutes in polynomial x^2 + x', () => {
      const x = new sym.core.Symbol('x');
      const two = new sym.core.Integer(2);
      const x2 = sym.core.pow(x, two);
      const expr = sym.core.add(x2, x);

      const three = new sym.core.Integer(3);
      const result = expr.subs(x, three);
      // 3^2 + 3 = 9 + 3 = 12
      expect(result.toString()).toBe('12');
    });

    it('substitutes symbol with expression', () => {
      const x = new sym.core.Symbol('x');
      const y = new sym.core.Symbol('y');
      const one = new sym.core.Integer(1);

      // Substitute x -> y + 1
      const yPlusOne = sym.core.add(y, one);
      const result = x.subs(x, yPlusOne);
      expect(result.toString()).toBe('1 + y');
    });
  });

  describe('compound expression substitution', () => {
    it('substitutes in addition', () => {
      const x = new sym.core.Symbol('x');
      const y = new sym.core.Symbol('y');
      const expr = sym.core.add(x, y);

      const one = new sym.core.Integer(1);
      const result = expr.subs(x, one);

      expect(result.toString()).toBe('1 + y');
    });

    it('substitutes in multiplication', () => {
      const x = new sym.core.Symbol('x');
      const y = new sym.core.Symbol('y');
      const expr = sym.core.mul(x, y);

      const two = new sym.core.Integer(2);
      const result = expr.subs(x, two);

      expect(result.toString()).toBe('2*y');
    });

    it('substitutes in power base', () => {
      const x = new sym.core.Symbol('x');
      const two = new sym.core.Integer(2);
      const expr = sym.core.pow(x, two);

      const three = new sym.core.Integer(3);
      const result = expr.subs(x, three);

      expect(result.toString()).toBe('9');
    });

    it('substitutes in power exponent', () => {
      const x = new sym.core.Symbol('x');
      const n = new sym.core.Symbol('n');
      const expr = sym.core.pow(x, n);

      const three = new sym.core.Integer(3);
      const result = expr.subs(n, three);

      expect(result.toString()).toBe('x**3');
    });
  });

  describe('no-op substitution', () => {
    it('returns same expression when symbol not present', () => {
      const x = new sym.core.Symbol('x');
      const y = new sym.core.Symbol('y');
      const z = new sym.core.Symbol('z');

      const expr = sym.core.add(x, y);
      const one = new sym.core.Integer(1);

      // z is not in expr, so result should equal expr
      const result = expr.subs(z, one);
      expect(result.toString()).toBe('x + y');
    });

    it('integer remains unchanged', () => {
      const five = new sym.core.Integer(5);
      const x = new sym.core.Symbol('x');
      const one = new sym.core.Integer(1);

      const result = five.subs(x, one);
      expect(result.toString()).toBe('5');
    });
  });

  describe('nested substitution', () => {
    it('substitutes in nested expression', () => {
      const x = new sym.core.Symbol('x');
      const y = new sym.core.Symbol('y');
      const two = new sym.core.Integer(2);

      // Build (x + y)^2
      const sum = sym.core.add(x, y);
      const expr = sym.core.pow(sum, two);

      // Substitute x -> 1
      const one = new sym.core.Integer(1);
      const result = expr.subs(x, one);

      // (1 + y)^2
      expect(result.toString()).toBe('(1 + y)**2');
    });

    it('chained substitution', () => {
      const x = new sym.core.Symbol('x');
      const y = new sym.core.Symbol('y');

      // Build x + y
      const expr = sym.core.add(x, y);

      // Substitute x -> 1, then y -> 2 (chained)
      const one = new sym.core.Integer(1);
      const two = new sym.core.Integer(2);

      const step1 = expr.subs(x, one);
      const result = step1.subs(y, two);

      expect(result.toString()).toBe('3');
    });
  });

  describe('constants substitution', () => {
    it('substitutes with pi', () => {
      const x = new sym.core.Symbol('x');
      const result = x.subs(x, sym.core.pi);
      expect(result.toString()).toBe('pi');
    });

    it('substitutes with E', () => {
      const x = new sym.core.Symbol('x');
      const result = x.subs(x, sym.core.E);
      expect(result.toString()).toBe('E');
    });

    it('substitutes with I (imaginary unit)', () => {
      const x = new sym.core.Symbol('x');
      const result = x.subs(x, sym.core.I);
      expect(result.toString()).toBe('I');
    });
  });

  describe('result type preservation', () => {
    it('returns Integer when result is numeric', () => {
      const x = new sym.core.Symbol('x');
      const five = new sym.core.Integer(5);
      const result = x.subs(x, five);
      expect(result.get_type()).toBe(sym.core.SymEngineTypeID.SYMENGINE_INTEGER);
    });

    it('returns Symbol when substituting symbol for symbol', () => {
      const x = new sym.core.Symbol('x');
      const y = new sym.core.Symbol('y');
      const result = x.subs(x, y);
      expect(result.get_type()).toBe(sym.core.SymEngineTypeID.SYMENGINE_SYMBOL);
    });

    it('returns Add when result is symbolic addition', () => {
      const x = new sym.core.Symbol('x');
      const y = new sym.core.Symbol('y');
      const expr = sym.core.add(x, y);

      const z = new sym.core.Symbol('z');
      const result = expr.subs(x, z);

      expect(result.get_type()).toBe(sym.core.SymEngineTypeID.SYMENGINE_ADD);
    });
  });

  describe('substitution with rationals', () => {
    it('substitutes symbol with rational', () => {
      const x = new sym.core.Symbol('x');
      const half = new sym.core.Rational(1, 2);
      const result = x.subs(x, half);
      expect(result.toString()).toBe('1/2');
    });

    it('substitutes in expression with rational result', () => {
      const x = new sym.core.Symbol('x');
      const two = new sym.core.Integer(2);
      // x / 2
      const expr = sym.core.div(x, two);

      const one = new sym.core.Integer(1);
      const result = expr.subs(x, one);

      expect(result.toString()).toBe('1/2');
    });
  });
});

describe('subs: Expr.subs(map) - multiple substitutions', () => {
  beforeAll(async () => {
    await sym.core.loadWasmModule();
  });

  describe('Map-based substitution', () => {
    it('substitutes multiple symbols with Map', () => {
      const x = new sym.core.Symbol('x');
      const y = new sym.core.Symbol('y');
      const one = new sym.core.Integer(1);
      const two = new sym.core.Integer(2);

      // x + y with {x: 1, y: 2}
      const expr = sym.core.add(x, y);
      const result = expr.subs(new Map([[x, one], [y, two]]));
      expect(result.toString()).toBe('3');
    });

    it('substitutes in polynomial x^2 + y^2', () => {
      const x = new sym.core.Symbol('x');
      const y = new sym.core.Symbol('y');
      const two = new sym.core.Integer(2);
      const x2 = sym.core.pow(x, two);
      const y2 = sym.core.pow(y, two);
      const expr = sym.core.add(x2, y2);

      // x^2 + y^2 with {x: 3, y: 4} = 9 + 16 = 25
      const three = new sym.core.Integer(3);
      const four = new sym.core.Integer(4);
      const result = expr.subs(new Map([[x, three], [y, four]]));
      expect(result.toString()).toBe('25');
    });

    it('handles partial substitution', () => {
      const x = new sym.core.Symbol('x');
      const y = new sym.core.Symbol('y');
      const z = new sym.core.Symbol('z');
      const expr = sym.core.add(sym.core.add(x, y), z);

      // x + y + z with {x: 1} leaves y + z + 1
      const one = new sym.core.Integer(1);
      const result = expr.subs(new Map([[x, one]]));
      expect(result.toString()).toMatch(/1 \+ y \+ z|y \+ z \+ 1/);
    });

    it('handles empty map (returns same expression)', () => {
      const x = new sym.core.Symbol('x');
      const result = x.subs(new Map());
      expect(result.toString()).toBe('x');
    });

    it('substitutes symbol with expression', () => {
      const x = new sym.core.Symbol('x');
      const y = new sym.core.Symbol('y');
      const expr = sym.core.mul(x, x); // x^2

      // x^2 with {x: y+1} = (y+1)^2
      const yPlusOne = sym.core.add(y, new sym.core.Integer(1));
      const result = expr.subs(new Map([[x, yPlusOne]]));
      expect(result.toString()).toBe('(1 + y)**2');
    });

    it('performs atomic (simultaneous) substitution', () => {
      const x = new sym.core.Symbol('x');
      const y = new sym.core.Symbol('y');
      const expr = sym.core.add(x, y);

      // Swap x and y: {x: y, y: x}
      // Atomic: x+y -> y+x (same thing)
      // Sequential would be: x+y -> y+y -> y+x (via second sub)
      const result = expr.subs(new Map([[x, y], [y, x]]));
      // The result should still be x + y (swapped)
      expect(result.toString()).toBe('x + y');
    });
  });

  describe('Object-based substitution', () => {
    it('substitutes multiple symbols with object notation', () => {
      const [x, y] = sym.core.symbols('x y');
      const expr = sym.core.add(x, y);

      // x + y with { x: 1, y: 2 }
      const result = expr.subs({ x: 1, y: 2 });
      expect(result.toString()).toBe('3');
    });

    it('substitutes in polynomial with object notation', () => {
      const [x, y] = sym.core.symbols('x y');
      const two = new sym.core.Integer(2);
      const x2 = sym.core.pow(x, two);
      const y2 = sym.core.pow(y, two);
      const expr = sym.core.add(x2, y2);

      // x^2 + y^2 with {x: 3, y: 4} = 9 + 16 = 25
      const result = expr.subs({ x: 3, y: 4 });
      expect(result.toString()).toBe('25');
    });

    it('accepts Expr values in object', () => {
      const [x, y] = sym.core.symbols('x y');
      const expr = sym.core.add(x, y);

      // Use Expr values instead of numbers
      const result = expr.subs({ x: sym.core.pi, y: sym.core.E });
      expect(result.toString()).toBe('E + pi');
    });

    it('substitutes single symbol with object notation', () => {
      const x = new sym.core.Symbol('x');
      const two = new sym.core.Integer(2);
      const expr = sym.core.pow(x, two);

      // x^2 with {x: 5} = 25
      const result = expr.subs({ x: 5 });
      expect(result.toString()).toBe('25');
    });
  });

  describe('complex expressions', () => {
    it('substitutes in (x + y) * (x - y)', () => {
      const [x, y] = sym.core.symbols('x y');
      const sum = sym.core.add(x, y);
      const diff = sym.core.sub(x, y);
      const expr = sym.core.mul(sum, diff);

      // (x + y)(x - y) = x^2 - y^2 with {x: 5, y: 3} = 25 - 9 = 16
      const result = expr.subs({ x: 5, y: 3 });
      expect(result.toString()).toBe('16');
    });

    it('substitutes in nested power expression', () => {
      const [x, y] = sym.core.symbols('x y');
      const two = new sym.core.Integer(2);
      // (x * y)^2
      const xy = sym.core.mul(x, y);
      const expr = sym.core.pow(xy, two);

      // (x*y)^2 with {x: 2, y: 3} = (2*3)^2 = 36
      const result = expr.subs({ x: 2, y: 3 });
      expect(result.toString()).toBe('36');
    });
  });
});
