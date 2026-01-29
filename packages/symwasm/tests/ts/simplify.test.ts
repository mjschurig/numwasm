import { describe, it, expect, beforeAll } from 'vitest';
import * as sym from '../../dist/symwasm.mjs';

describe('simplify module', () => {
  beforeAll(async () => {
    await sym.core.loadWasmModule();
  });

  describe('expand()', () => {
    it('should expand (x+1)^2 to x^2 + 2*x + 1', () => {
      const x = new sym.core.Symbol('x');
      const one = new sym.core.Integer(1);
      const two = new sym.core.Integer(2);
      const expr = sym.core.pow(sym.core.add(x, one), two);
      const result = sym.simplify.expand(expr);
      // (x+1)^2 = x^2 + 2*x + 1
      const str = result.toString();
      expect(str).toContain('x**2');
      expect(str).toContain('2*x');
      expect(str).toContain('1');
    });

    it('should expand (a+b)*(c+d)', () => {
      const a = new sym.core.Symbol('a');
      const b = new sym.core.Symbol('b');
      const c = new sym.core.Symbol('c');
      const d = new sym.core.Symbol('d');
      const expr = sym.core.mul(sym.core.add(a, b), sym.core.add(c, d));
      const result = sym.simplify.expand(expr);
      const str = result.toString();
      // Should contain all four terms: a*c, a*d, b*c, b*d
      expect(str).toContain('a*c');
      expect(str).toContain('a*d');
      expect(str).toContain('b*c');
      expect(str).toContain('b*d');
    });

    it('should expand (x+y)^3', () => {
      const x = new sym.core.Symbol('x');
      const y = new sym.core.Symbol('y');
      const three = new sym.core.Integer(3);
      const expr = sym.core.pow(sym.core.add(x, y), three);
      const result = sym.simplify.expand(expr);
      const str = result.toString();
      // (x+y)^3 = x^3 + 3*x^2*y + 3*x*y^2 + y^3
      expect(str).toContain('x**3');
      expect(str).toContain('y**3');
    });

    it('should leave already expanded expressions unchanged', () => {
      const x = new sym.core.Symbol('x');
      const y = new sym.core.Symbol('y');
      const expr = sym.core.add(x, y);
      const result = sym.simplify.expand(expr);
      expect(result.toString()).toBe('x + y');
    });

    it('should expand 2*(x+1)', () => {
      const x = new sym.core.Symbol('x');
      const one = new sym.core.Integer(1);
      const two = new sym.core.Integer(2);
      const expr = sym.core.mul(two, sym.core.add(x, one));
      const result = sym.simplify.expand(expr);
      expect(result.toString()).toBe('2 + 2*x');
    });
  });

  describe('simplify()', () => {
    it('should simplify csc(x)^(-1) to sin(x)', () => {
      const x = new sym.core.Symbol('x');
      const negOne = new sym.core.Integer(-1);
      const expr = sym.core.pow(sym.core.csc(x), negOne);
      const result = sym.simplify.simplify(expr);
      expect(result.toString()).toBe('sin(x)');
    });

    it('should simplify sec(x)^(-1) to cos(x)', () => {
      const x = new sym.core.Symbol('x');
      const negOne = new sym.core.Integer(-1);
      const expr = sym.core.pow(sym.core.sec(x), negOne);
      const result = sym.simplify.simplify(expr);
      expect(result.toString()).toBe('cos(x)');
    });

    it('should simplify cot(x)^(-1) to tan(x)', () => {
      const x = new sym.core.Symbol('x');
      const negOne = new sym.core.Integer(-1);
      const expr = sym.core.pow(sym.core.cot(x), negOne);
      const result = sym.simplify.simplify(expr);
      expect(result.toString()).toBe('tan(x)');
    });

    it('should return symbols unchanged', () => {
      const x = new sym.core.Symbol('x');
      const result = sym.simplify.simplify(x);
      expect(result.toString()).toBe('x');
    });

    it('should return integers unchanged', () => {
      const n = new sym.core.Integer(42);
      const result = sym.simplify.simplify(n);
      expect(result.toString()).toBe('42');
    });
  });

  describe('numer()', () => {
    it('should extract numerator from x/y', () => {
      const x = new sym.core.Symbol('x');
      const y = new sym.core.Symbol('y');
      const expr = sym.core.div(x, y);
      const result = sym.simplify.numer(expr);
      expect(result.toString()).toBe('x');
    });

    it('should return integer as-is', () => {
      const n = new sym.core.Integer(5);
      const result = sym.simplify.numer(n);
      expect(result.toString()).toBe('5');
    });

    it('should extract numerator from Rational', () => {
      const r = new sym.core.Rational(3, 4);
      const result = sym.simplify.numer(r);
      expect(result.toString()).toBe('3');
    });

    it('should extract numerator from (x+1)/(x-1)', () => {
      const x = new sym.core.Symbol('x');
      const one = new sym.core.Integer(1);
      const expr = sym.core.div(sym.core.add(x, one), sym.core.sub(x, one));
      const result = sym.simplify.numer(expr);
      // SymEngine may order terms differently
      expect(result.toString()).toMatch(/^(x \+ 1|1 \+ x)$/);
    });
  });

  describe('denom()', () => {
    it('should extract denominator from x/y', () => {
      const x = new sym.core.Symbol('x');
      const y = new sym.core.Symbol('y');
      const expr = sym.core.div(x, y);
      const result = sym.simplify.denom(expr);
      expect(result.toString()).toBe('y');
    });

    it('should return 1 for integer', () => {
      const n = new sym.core.Integer(5);
      const result = sym.simplify.denom(n);
      expect(result.toString()).toBe('1');
    });

    it('should extract denominator from Rational', () => {
      const r = new sym.core.Rational(3, 4);
      const result = sym.simplify.denom(r);
      expect(result.toString()).toBe('4');
    });

    it('should extract denominator from (x+1)/(x-1)', () => {
      const x = new sym.core.Symbol('x');
      const one = new sym.core.Integer(1);
      const expr = sym.core.div(sym.core.add(x, one), sym.core.sub(x, one));
      const result = sym.simplify.denom(expr);
      // SymEngine may order terms differently
      expect(result.toString()).toMatch(/^(x - 1|-1 \+ x)$/);
    });

    it('should return 1 for symbol', () => {
      const x = new sym.core.Symbol('x');
      const result = sym.simplify.denom(x);
      expect(result.toString()).toBe('1');
    });
  });

  describe('trigsimp()', () => {
    it('should simplify csc(x)^(-1) to sin(x)', () => {
      const x = new sym.core.Symbol('x');
      const negOne = new sym.core.Integer(-1);
      const expr = sym.core.pow(sym.core.csc(x), negOne);
      const result = sym.simplify.trigsimp(expr);
      expect(result.toString()).toBe('sin(x)');
    });
  });

  describe('radsimp()', () => {
    it('should return expression (delegates to simplify)', () => {
      const x = new sym.core.Symbol('x');
      const result = sym.simplify.radsimp(x);
      expect(result.toString()).toBe('x');
    });
  });

  describe('powsimp()', () => {
    it('should return expression (delegates to simplify)', () => {
      const x = new sym.core.Symbol('x');
      const result = sym.simplify.powsimp(x);
      expect(result.toString()).toBe('x');
    });
  });

  describe('rewrite_as_exp()', () => {
    it('should rewrite sin(x) as exponentials', () => {
      const x = new sym.core.Symbol('x');
      const sinx = sym.core.sin(x);
      const result = sym.simplify.rewrite_as_exp(sinx);
      const str = result.toString();
      // sin(x) = (e^(ix) - e^(-ix)) / (2i)
      expect(str).toContain('exp');
      expect(str).toContain('I');
    });

    it('should rewrite cos(x) as exponentials', () => {
      const x = new sym.core.Symbol('x');
      const cosx = sym.core.cos(x);
      const result = sym.simplify.rewrite_as_exp(cosx);
      const str = result.toString();
      // cos(x) = (e^(ix) + e^(-ix)) / 2
      expect(str).toContain('exp');
    });

    it('should leave non-trig expressions unchanged', () => {
      const x = new sym.core.Symbol('x');
      const result = sym.simplify.rewrite_as_exp(x);
      expect(result.toString()).toBe('x');
    });
  });

  describe('rewrite_as_sin()', () => {
    it('should rewrite cos(x) in terms of sine', () => {
      const x = new sym.core.Symbol('x');
      const cosx = sym.core.cos(x);
      const result = sym.simplify.rewrite_as_sin(cosx);
      const str = result.toString();
      // cos(x) = sin(x + π/2)
      expect(str).toContain('sin');
    });

    it('should leave sin(x) as-is', () => {
      const x = new sym.core.Symbol('x');
      const sinx = sym.core.sin(x);
      const result = sym.simplify.rewrite_as_sin(sinx);
      expect(result.toString()).toBe('sin(x)');
    });
  });

  describe('rewrite_as_cos()', () => {
    it('should rewrite sin(x) in terms of cosine', () => {
      const x = new sym.core.Symbol('x');
      const sinx = sym.core.sin(x);
      const result = sym.simplify.rewrite_as_cos(sinx);
      const str = result.toString();
      // sin(x) = cos(π/2 - x)
      expect(str).toContain('cos');
    });

    it('should leave cos(x) as-is', () => {
      const x = new sym.core.Symbol('x');
      const cosx = sym.core.cos(x);
      const result = sym.simplify.rewrite_as_cos(cosx);
      expect(result.toString()).toBe('cos(x)');
    });
  });

  describe('as_real_imag()', () => {
    it('should extract real and imaginary parts of x + i*y', () => {
      const x = new sym.core.Symbol('x');
      const y = new sym.core.Symbol('y');
      const z = sym.core.add(x, sym.core.mul(sym.core.I, y)); // x + i*y
      const { real, imag } = sym.simplify.as_real_imag(z);
      expect(real.toString()).toBe('x');
      expect(imag.toString()).toBe('y');
    });

    it('should return 0 for imaginary part of a real symbol', () => {
      const x = new sym.core.Symbol('x');
      const { real, imag } = sym.simplify.as_real_imag(x);
      expect(real.toString()).toBe('x');
      expect(imag.toString()).toBe('0');
    });

    it('should extract parts from pure imaginary', () => {
      const y = new sym.core.Symbol('y');
      const iz = sym.core.mul(sym.core.I, y); // i*y
      const { real, imag } = sym.simplify.as_real_imag(iz);
      expect(real.toString()).toBe('0');
      expect(imag.toString()).toBe('y');
    });

    it('should extract parts from complex number', () => {
      const three = new sym.core.Integer(3);
      const four = new sym.core.Integer(4);
      const z = sym.core.add(three, sym.core.mul(sym.core.I, four)); // 3 + 4i
      const { real, imag } = sym.simplify.as_real_imag(z);
      expect(real.toString()).toBe('3');
      expect(imag.toString()).toBe('4');
    });
  });

  describe('expand_trig()', () => {
    it('should expand sin(x) using exponentials', () => {
      const x = new sym.core.Symbol('x');
      const sinx = sym.core.sin(x);
      const result = sym.simplify.expand_trig(sinx);
      const str = result.toString();
      // expand(rewrite_as_exp(sin(x)))
      expect(str).toContain('exp');
    });
  });

  describe('expand_complex()', () => {
    it('should be alias for as_real_imag', () => {
      const x = new sym.core.Symbol('x');
      const y = new sym.core.Symbol('y');
      const z = sym.core.add(x, sym.core.mul(sym.core.I, y));
      const { real, imag } = sym.simplify.expand_complex(z);
      expect(real.toString()).toBe('x');
      expect(imag.toString()).toBe('y');
    });
  });
});

describe('core module simplification exports', () => {
  beforeAll(async () => {
    await sym.core.loadWasmModule();
  });

  it('expand should be accessible from core', () => {
    const x = new sym.core.Symbol('x');
    const one = new sym.core.Integer(1);
    const two = new sym.core.Integer(2);
    const expr = sym.core.pow(sym.core.add(x, one), two);
    const result = sym.core.expand(expr);
    expect(result.toString()).toContain('x**2');
  });

  it('simplify should be accessible from core', () => {
    const x = new sym.core.Symbol('x');
    const negOne = new sym.core.Integer(-1);
    const expr = sym.core.pow(sym.core.csc(x), negOne);
    const result = sym.core.simplify(expr);
    expect(result.toString()).toBe('sin(x)');
  });

  it('numer should be accessible from core', () => {
    const r = new sym.core.Rational(3, 4);
    expect(sym.core.numer(r).toString()).toBe('3');
  });

  it('denom should be accessible from core', () => {
    const r = new sym.core.Rational(3, 4);
    expect(sym.core.denom(r).toString()).toBe('4');
  });
});
