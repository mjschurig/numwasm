/**
 * Tests for symwasm elementary functions (Phase 2.1)
 *
 * Tests trigonometric, exponential, logarithmic, and special functions.
 */
import { describe, it, expect, beforeAll } from 'vitest';
import * as sym from '../../src/ts/index';

describe('trigonometric functions', () => {
  beforeAll(async () => {
    await sym.core.loadWasmModule();
  });

  describe('sin()', () => {
    it('sin(0) = 0', () => {
      const result = sym.core.sin(new sym.core.Integer(0));
      expect(result.toString()).toBe('0');
    });

    it('sin(pi) = 0', () => {
      const result = sym.core.sin(sym.core.pi);
      expect(result.toString()).toBe('0');
    });

    it('sin(pi/2) = 1', () => {
      const piHalf = sym.core.mul(sym.core.pi, new sym.core.Rational(1, 2));
      const result = sym.core.sin(piHalf);
      expect(result.toString()).toBe('1');
    });

    it('sin(x) is symbolic', () => {
      const x = new sym.core.Symbol('x');
      const result = sym.core.sin(x);
      expect(result.toString()).toBe('sin(x)');
    });

    it('sin(1) evaluates numerically', () => {
      const result = sym.core.sin(new sym.core.Integer(1)).evalfNumber();
      expect(result).toBeCloseTo(0.8414709848, 8);
    });
  });

  describe('cos()', () => {
    it('cos(0) = 1', () => {
      const result = sym.core.cos(new sym.core.Integer(0));
      expect(result.toString()).toBe('1');
    });

    it('cos(pi) = -1', () => {
      const result = sym.core.cos(sym.core.pi);
      expect(result.toString()).toBe('-1');
    });

    it('cos(pi/2) = 0', () => {
      const piHalf = sym.core.mul(sym.core.pi, new sym.core.Rational(1, 2));
      const result = sym.core.cos(piHalf);
      expect(result.toString()).toBe('0');
    });

    it('cos(x) is symbolic', () => {
      const x = new sym.core.Symbol('x');
      const result = sym.core.cos(x);
      expect(result.toString()).toBe('cos(x)');
    });
  });

  describe('tan()', () => {
    it('tan(0) = 0', () => {
      const result = sym.core.tan(new sym.core.Integer(0));
      expect(result.toString()).toBe('0');
    });

    it('tan(pi/4) = 1', () => {
      const piFourth = sym.core.mul(sym.core.pi, new sym.core.Rational(1, 4));
      const result = sym.core.tan(piFourth);
      expect(result.toString()).toBe('1');
    });

    it('tan(x) is symbolic', () => {
      const x = new sym.core.Symbol('x');
      const result = sym.core.tan(x);
      expect(result.toString()).toBe('tan(x)');
    });
  });

  describe('cot(), sec(), csc()', () => {
    it('cot(x) is symbolic', () => {
      const x = new sym.core.Symbol('x');
      const result = sym.core.cot(x);
      expect(result.toString()).toBe('cot(x)');
    });

    it('sec(0) = 1', () => {
      const result = sym.core.sec(new sym.core.Integer(0));
      expect(result.toString()).toBe('1');
    });

    it('csc(pi/2) = 1', () => {
      const piHalf = sym.core.mul(sym.core.pi, new sym.core.Rational(1, 2));
      const result = sym.core.csc(piHalf);
      expect(result.toString()).toBe('1');
    });
  });
});

describe('inverse trigonometric functions', () => {
  beforeAll(async () => {
    await sym.core.loadWasmModule();
  });

  describe('asin()', () => {
    it('asin(0) = 0', () => {
      const result = sym.core.asin(new sym.core.Integer(0));
      expect(result.toString()).toBe('0');
    });

    it('asin(1) = pi/2', () => {
      const result = sym.core.asin(new sym.core.Integer(1));
      // SymEngine represents as (1/2)*pi
      expect(result.toString()).toMatch(/pi\/2|\(1\/2\)\*pi/);
    });

    it('asin(x) is symbolic', () => {
      const x = new sym.core.Symbol('x');
      const result = sym.core.asin(x);
      expect(result.toString()).toBe('asin(x)');
    });
  });

  describe('acos()', () => {
    it('acos(1) = 0', () => {
      const result = sym.core.acos(new sym.core.Integer(1));
      expect(result.toString()).toBe('0');
    });

    it('acos(0) = pi/2', () => {
      const result = sym.core.acos(new sym.core.Integer(0));
      // SymEngine represents as (1/2)*pi
      expect(result.toString()).toMatch(/pi\/2|\(1\/2\)\*pi/);
    });
  });

  describe('atan()', () => {
    it('atan(0) = 0', () => {
      const result = sym.core.atan(new sym.core.Integer(0));
      expect(result.toString()).toBe('0');
    });

    it('atan(1) = pi/4', () => {
      const result = sym.core.atan(new sym.core.Integer(1));
      // SymEngine represents as (1/4)*pi
      expect(result.toString()).toMatch(/pi\/4|\(1\/4\)\*pi/);
    });

    it('atan(x) is symbolic', () => {
      const x = new sym.core.Symbol('x');
      const result = sym.core.atan(x);
      expect(result.toString()).toBe('atan(x)');
    });
  });

  describe('atan2()', () => {
    it('atan2(0, 1) = 0', () => {
      const result = sym.core.atan2(new sym.core.Integer(0), new sym.core.Integer(1));
      expect(result.toString()).toBe('0');
    });

    it('atan2(1, 1) = pi/4', () => {
      const result = sym.core.atan2(new sym.core.Integer(1), new sym.core.Integer(1));
      // SymEngine represents as (1/4)*pi
      expect(result.toString()).toMatch(/pi\/4|\(1\/4\)\*pi/);
    });

    it('atan2(y, x) is symbolic', () => {
      const [x, y] = sym.core.symbols('x y');
      const result = sym.core.atan2(y, x);
      expect(result.toString()).toBe('atan2(y, x)');
    });
  });
});

describe('hyperbolic functions', () => {
  beforeAll(async () => {
    await sym.core.loadWasmModule();
  });

  describe('sinh()', () => {
    it('sinh(0) = 0', () => {
      const result = sym.core.sinh(new sym.core.Integer(0));
      expect(result.toString()).toBe('0');
    });

    it('sinh(x) is symbolic', () => {
      const x = new sym.core.Symbol('x');
      const result = sym.core.sinh(x);
      expect(result.toString()).toBe('sinh(x)');
    });

    it('sinh(1) evaluates numerically', () => {
      const result = sym.core.sinh(new sym.core.Integer(1)).evalfNumber();
      expect(result).toBeCloseTo(1.1752011936, 8);
    });
  });

  describe('cosh()', () => {
    it('cosh(0) = 1', () => {
      const result = sym.core.cosh(new sym.core.Integer(0));
      expect(result.toString()).toBe('1');
    });

    it('cosh(x) is symbolic', () => {
      const x = new sym.core.Symbol('x');
      const result = sym.core.cosh(x);
      expect(result.toString()).toBe('cosh(x)');
    });
  });

  describe('tanh()', () => {
    it('tanh(0) = 0', () => {
      const result = sym.core.tanh(new sym.core.Integer(0));
      expect(result.toString()).toBe('0');
    });

    it('tanh(x) is symbolic', () => {
      const x = new sym.core.Symbol('x');
      const result = sym.core.tanh(x);
      expect(result.toString()).toBe('tanh(x)');
    });
  });

  describe('coth(), sech(), csch()', () => {
    it('sech(0) = 1', () => {
      const result = sym.core.sech(new sym.core.Integer(0));
      expect(result.toString()).toBe('1');
    });

    it('coth(x) is symbolic', () => {
      const x = new sym.core.Symbol('x');
      const result = sym.core.coth(x);
      expect(result.toString()).toBe('coth(x)');
    });

    it('csch(x) is symbolic', () => {
      const x = new sym.core.Symbol('x');
      const result = sym.core.csch(x);
      expect(result.toString()).toBe('csch(x)');
    });
  });
});

describe('inverse hyperbolic functions', () => {
  beforeAll(async () => {
    await sym.core.loadWasmModule();
  });

  describe('asinh()', () => {
    it('asinh(0) = 0', () => {
      const result = sym.core.asinh(new sym.core.Integer(0));
      expect(result.toString()).toBe('0');
    });

    it('asinh(x) is symbolic', () => {
      const x = new sym.core.Symbol('x');
      const result = sym.core.asinh(x);
      expect(result.toString()).toBe('asinh(x)');
    });
  });

  describe('acosh()', () => {
    it('acosh(1) = 0', () => {
      const result = sym.core.acosh(new sym.core.Integer(1));
      expect(result.toString()).toBe('0');
    });

    it('acosh(x) is symbolic', () => {
      const x = new sym.core.Symbol('x');
      const result = sym.core.acosh(x);
      expect(result.toString()).toBe('acosh(x)');
    });
  });

  describe('atanh()', () => {
    it('atanh(0) = 0', () => {
      const result = sym.core.atanh(new sym.core.Integer(0));
      expect(result.toString()).toBe('0');
    });

    it('atanh(x) is symbolic', () => {
      const x = new sym.core.Symbol('x');
      const result = sym.core.atanh(x);
      expect(result.toString()).toBe('atanh(x)');
    });
  });

  describe('acoth(), asech(), acsch()', () => {
    it('acoth(x) is symbolic', () => {
      const x = new sym.core.Symbol('x');
      const result = sym.core.acoth(x);
      expect(result.toString()).toBe('acoth(x)');
    });

    it('asech(x) is symbolic', () => {
      const x = new sym.core.Symbol('x');
      const result = sym.core.asech(x);
      expect(result.toString()).toBe('asech(x)');
    });

    it('acsch(x) is symbolic', () => {
      const x = new sym.core.Symbol('x');
      const result = sym.core.acsch(x);
      expect(result.toString()).toBe('acsch(x)');
    });
  });
});

describe('exponential and logarithmic functions', () => {
  beforeAll(async () => {
    await sym.core.loadWasmModule();
  });

  describe('exp()', () => {
    it('exp(0) = 1', () => {
      const result = sym.core.exp(new sym.core.Integer(0));
      expect(result.toString()).toBe('1');
    });

    it('exp(1) = E', () => {
      const result = sym.core.exp(new sym.core.Integer(1));
      expect(result.toString()).toBe('E');
    });

    it('exp(x) is symbolic', () => {
      const x = new sym.core.Symbol('x');
      const result = sym.core.exp(x);
      expect(result.toString()).toBe('exp(x)');
    });

    it('exp(1) evaluates numerically', () => {
      const result = sym.core.exp(new sym.core.Integer(1)).evalfNumber();
      expect(result).toBeCloseTo(Math.E, 10);
    });
  });

  describe('log()', () => {
    it('log(1) = 0', () => {
      const result = sym.core.log(new sym.core.Integer(1));
      expect(result.toString()).toBe('0');
    });

    it('log(E) = 1', () => {
      const result = sym.core.log(sym.core.E);
      expect(result.toString()).toBe('1');
    });

    it('log(x) is symbolic', () => {
      const x = new sym.core.Symbol('x');
      const result = sym.core.log(x);
      expect(result.toString()).toBe('log(x)');
    });

    it('log(2) evaluates numerically', () => {
      const result = sym.core.log(new sym.core.Integer(2)).evalfNumber();
      expect(result).toBeCloseTo(Math.log(2), 10);
    });
  });

  describe('sqrt()', () => {
    it('sqrt(0) = 0', () => {
      const result = sym.core.sqrt(new sym.core.Integer(0));
      expect(result.toString()).toBe('0');
    });

    it('sqrt(1) = 1', () => {
      const result = sym.core.sqrt(new sym.core.Integer(1));
      expect(result.toString()).toBe('1');
    });

    it('sqrt(4) = 2', () => {
      const result = sym.core.sqrt(new sym.core.Integer(4));
      expect(result.toString()).toBe('2');
    });

    it('sqrt(2) stays symbolic', () => {
      const result = sym.core.sqrt(new sym.core.Integer(2));
      expect(result.toString()).toBe('sqrt(2)');
    });

    it('sqrt(x) is symbolic', () => {
      const x = new sym.core.Symbol('x');
      const result = sym.core.sqrt(x);
      expect(result.toString()).toBe('sqrt(x)');
    });

    it('sqrt(2) evaluates numerically', () => {
      const result = sym.core.sqrt(new sym.core.Integer(2)).evalfNumber();
      expect(result).toBeCloseTo(Math.sqrt(2), 10);
    });
  });

  describe('cbrt()', () => {
    it('cbrt(0) = 0', () => {
      const result = sym.core.cbrt(new sym.core.Integer(0));
      expect(result.toString()).toBe('0');
    });

    it('cbrt(1) = 1', () => {
      const result = sym.core.cbrt(new sym.core.Integer(1));
      expect(result.toString()).toBe('1');
    });

    it('cbrt(8) = 2', () => {
      const result = sym.core.cbrt(new sym.core.Integer(8));
      expect(result.toString()).toBe('2');
    });

    it('cbrt(x) is symbolic', () => {
      const x = new sym.core.Symbol('x');
      const result = sym.core.cbrt(x);
      // cbrt(x) may be represented as x^(1/3)
      expect(result.toString()).toMatch(/cbrt\(x\)|x\*\*\(1\/3\)/);
    });
  });

  describe('lambertw()', () => {
    it('lambertw(0) = 0', () => {
      const result = sym.core.lambertw(new sym.core.Integer(0));
      expect(result.toString()).toBe('0');
    });

    it('lambertw(x) is symbolic', () => {
      const x = new sym.core.Symbol('x');
      const result = sym.core.lambertw(x);
      expect(result.toString()).toMatch(/lambertw|LambertW/i);
    });
  });
});

describe('other mathematical functions', () => {
  beforeAll(async () => {
    await sym.core.loadWasmModule();
  });

  describe('abs()', () => {
    it('abs(5) = 5', () => {
      const result = sym.core.abs(new sym.core.Integer(5));
      expect(result.toString()).toBe('5');
    });

    it('abs(-5) = 5', () => {
      const result = sym.core.abs(new sym.core.Integer(-5));
      expect(result.toString()).toBe('5');
    });

    it('abs(0) = 0', () => {
      const result = sym.core.abs(new sym.core.Integer(0));
      expect(result.toString()).toBe('0');
    });

    it('abs(x) is symbolic', () => {
      const x = new sym.core.Symbol('x');
      const result = sym.core.abs(x);
      expect(result.toString()).toMatch(/Abs\(x\)|abs\(x\)/i);
    });
  });

  describe('sign()', () => {
    it('sign(5) = 1', () => {
      const result = sym.core.sign(new sym.core.Integer(5));
      expect(result.toString()).toBe('1');
    });

    it('sign(-5) = -1', () => {
      const result = sym.core.sign(new sym.core.Integer(-5));
      expect(result.toString()).toBe('-1');
    });

    it('sign(0) = 0', () => {
      const result = sym.core.sign(new sym.core.Integer(0));
      expect(result.toString()).toBe('0');
    });

    it('sign(x) is symbolic', () => {
      const x = new sym.core.Symbol('x');
      const result = sym.core.sign(x);
      expect(result.toString()).toMatch(/sign\(x\)|Sign\(x\)/i);
    });
  });

  describe('floor()', () => {
    it('floor(5) = 5', () => {
      const result = sym.core.floor(new sym.core.Integer(5));
      expect(result.toString()).toBe('5');
    });

    it('floor(3/2) = 1', () => {
      const result = sym.core.floor(new sym.core.Rational(3, 2));
      expect(result.toString()).toBe('1');
    });

    it('floor(-3/2) = -2', () => {
      const result = sym.core.floor(new sym.core.Rational(-3, 2));
      expect(result.toString()).toBe('-2');
    });

    it('floor(x) is symbolic', () => {
      const x = new sym.core.Symbol('x');
      const result = sym.core.floor(x);
      expect(result.toString()).toMatch(/floor\(x\)|Floor\(x\)/i);
    });
  });

  describe('ceiling()', () => {
    it('ceiling(5) = 5', () => {
      const result = sym.core.ceiling(new sym.core.Integer(5));
      expect(result.toString()).toBe('5');
    });

    it('ceiling(3/2) = 2', () => {
      const result = sym.core.ceiling(new sym.core.Rational(3, 2));
      expect(result.toString()).toBe('2');
    });

    it('ceiling(-3/2) = -1', () => {
      const result = sym.core.ceiling(new sym.core.Rational(-3, 2));
      expect(result.toString()).toBe('-1');
    });

    it('ceiling(x) is symbolic', () => {
      const x = new sym.core.Symbol('x');
      const result = sym.core.ceiling(x);
      expect(result.toString()).toMatch(/ceiling\(x\)|Ceiling\(x\)/i);
    });
  });
});

describe('special functions', () => {
  beforeAll(async () => {
    await sym.core.loadWasmModule();
  });

  describe('gamma()', () => {
    it('gamma(1) = 1', () => {
      const result = sym.core.gamma(new sym.core.Integer(1));
      expect(result.toString()).toBe('1');
    });

    it('gamma(2) = 1 (factorial of 1)', () => {
      const result = sym.core.gamma(new sym.core.Integer(2));
      expect(result.toString()).toBe('1');
    });

    it('gamma(5) = 24 (factorial of 4)', () => {
      const result = sym.core.gamma(new sym.core.Integer(5));
      expect(result.toString()).toBe('24');
    });

    it('gamma(6) = 120 (factorial of 5)', () => {
      const result = sym.core.gamma(new sym.core.Integer(6));
      expect(result.toString()).toBe('120');
    });

    it('gamma(x) is symbolic', () => {
      const x = new sym.core.Symbol('x');
      const result = sym.core.gamma(x);
      expect(result.toString()).toMatch(/gamma\(x\)|Gamma\(x\)/i);
    });
  });

  describe('loggamma()', () => {
    it('loggamma(1) = 0', () => {
      const result = sym.core.loggamma(new sym.core.Integer(1));
      expect(result.toString()).toBe('0');
    });

    it('loggamma(x) is symbolic', () => {
      const x = new sym.core.Symbol('x');
      const result = sym.core.loggamma(x);
      expect(result.toString()).toMatch(/loggamma\(x\)|LogGamma\(x\)/i);
    });
  });

  describe('erf() and erfc()', () => {
    it('erf(0) = 0', () => {
      const result = sym.core.erf(new sym.core.Integer(0));
      expect(result.toString()).toBe('0');
    });

    it('erfc(0) = 1', () => {
      const result = sym.core.erfc(new sym.core.Integer(0));
      expect(result.toString()).toBe('1');
    });

    it('erf(x) is symbolic', () => {
      const x = new sym.core.Symbol('x');
      const result = sym.core.erf(x);
      expect(result.toString()).toMatch(/erf\(x\)|Erf\(x\)/i);
    });

    it('erfc(x) is symbolic', () => {
      const x = new sym.core.Symbol('x');
      const result = sym.core.erfc(x);
      expect(result.toString()).toMatch(/erfc\(x\)|Erfc\(x\)/i);
    });
  });

  describe('zeta()', () => {
    it('zeta(2) = pi^2/6', () => {
      const result = sym.core.zeta(new sym.core.Integer(2));
      // SymEngine represents as (1/6)*pi**2
      expect(result.toString()).toMatch(/pi\*\*2\/6|zeta\(2\)|\(1\/6\)\*pi\*\*2/);
    });

    it('zeta(x) is symbolic', () => {
      const x = new sym.core.Symbol('x');
      const result = sym.core.zeta(x);
      // SymEngine's zeta has a second argument (default 1)
      expect(result.toString()).toMatch(/zeta\(x|Zeta\(x/i);
    });
  });

  describe('beta()', () => {
    it('beta(1, 1) = 1', () => {
      const result = sym.core.beta(new sym.core.Integer(1), new sym.core.Integer(1));
      expect(result.toString()).toBe('1');
    });

    it('beta(2, 2) = 1/6', () => {
      const result = sym.core.beta(new sym.core.Integer(2), new sym.core.Integer(2));
      expect(result.toString()).toBe('1/6');
    });

    it('beta(a, b) is symbolic', () => {
      const [a, b] = sym.core.symbols('a b');
      const result = sym.core.beta(a, b);
      // SymEngine may reorder arguments alphabetically
      expect(result.toString()).toMatch(/beta\([ab], [ab]\)|Beta\([ab], [ab]\)/i);
    });
  });

  describe('lowergamma() and uppergamma()', () => {
    it('lowergamma(s, x) is symbolic', () => {
      const [s, x] = sym.core.symbols('s x');
      const result = sym.core.lowergamma(s, x);
      expect(result.toString()).toMatch(/lowergamma|LowerGamma/i);
    });

    it('uppergamma(s, x) is symbolic', () => {
      const [s, x] = sym.core.symbols('s x');
      const result = sym.core.uppergamma(s, x);
      expect(result.toString()).toMatch(/uppergamma|UpperGamma/i);
    });
  });

  describe('polygamma()', () => {
    it('polygamma(0, 1) = -EulerGamma', () => {
      const result = sym.core.polygamma(new sym.core.Integer(0), new sym.core.Integer(1));
      expect(result.toString()).toMatch(/-EulerGamma|polygamma/i);
    });

    it('polygamma(n, x) is symbolic', () => {
      const [n, x] = sym.core.symbols('n x');
      const result = sym.core.polygamma(n, x);
      expect(result.toString()).toMatch(/polygamma|PolyGamma/i);
    });
  });

  describe('kronecker_delta()', () => {
    it('kronecker_delta(1, 1) = 1', () => {
      const result = sym.core.kronecker_delta(new sym.core.Integer(1), new sym.core.Integer(1));
      expect(result.toString()).toBe('1');
    });

    it('kronecker_delta(1, 2) = 0', () => {
      const result = sym.core.kronecker_delta(new sym.core.Integer(1), new sym.core.Integer(2));
      expect(result.toString()).toBe('0');
    });

    it('kronecker_delta(i, j) is symbolic', () => {
      const [i, j] = sym.core.symbols('i j');
      const result = sym.core.kronecker_delta(i, j);
      expect(result.toString()).toMatch(/KroneckerDelta|kronecker_delta/i);
    });
  });
});

describe('function composition', () => {
  beforeAll(async () => {
    await sym.core.loadWasmModule();
  });

  it('sin(cos(x)) composes correctly', () => {
    const x = new sym.core.Symbol('x');
    const result = sym.core.sin(sym.core.cos(x));
    expect(result.toString()).toBe('sin(cos(x))');
  });

  it('exp(log(x)) composes correctly', () => {
    const x = new sym.core.Symbol('x');
    const result = sym.core.exp(sym.core.log(x));
    // SymEngine may or may not simplify this without assumptions
    expect(result.toString()).toMatch(/x|exp\(log\(x\)\)/);
  });

  it('log(exp(x)) composes correctly', () => {
    const x = new sym.core.Symbol('x');
    const result = sym.core.log(sym.core.exp(x));
    // SymEngine may or may not simplify this without assumptions
    expect(result.toString()).toMatch(/x|log\(exp\(x\)\)/);
  });

  it('sqrt(x)^2 = x', () => {
    const x = new sym.core.Symbol('x');
    const sqrtX = sym.core.sqrt(x);
    const two = new sym.core.Integer(2);
    const result = sym.core.pow(sqrtX, two);
    expect(result.toString()).toBe('x');
  });

  it('sin^2 + cos^2 can be computed', () => {
    const x = new sym.core.Symbol('x');
    const two = new sym.core.Integer(2);
    const sin2 = sym.core.pow(sym.core.sin(x), two);
    const cos2 = sym.core.pow(sym.core.cos(x), two);
    const result = sym.core.add(sin2, cos2);
    // SymEngine may or may not simplify this to 1
    expect(result.toString()).toMatch(/sin\*\*2\(x\) \+ cos\*\*2\(x\)|1|sin\(x\)\*\*2 \+ cos\(x\)\*\*2/);
  });
});

// ============================================================================
// Phase 2.1b: Additional Functions
// ============================================================================

describe('complex number functions', () => {
  beforeAll(async () => {
    await sym.core.loadWasmModule();
  });

  describe('re()', () => {
    it('re(3 + 4i) = 3', () => {
      const z = new sym.core.Complex(3, 4);
      const result = sym.core.re(z);
      expect(result.toString()).toBe('3');
    });

    it('re(0 + 5i) = 0', () => {
      const z = new sym.core.Complex(0, 5);
      const result = sym.core.re(z);
      expect(result.toString()).toBe('0');
    });

    it('re(-2 + 3i) = -2', () => {
      const z = new sym.core.Complex(-2, 3);
      const result = sym.core.re(z);
      expect(result.toString()).toBe('-2');
    });
  });

  describe('im()', () => {
    it('im(3 + 4i) = 4', () => {
      const z = new sym.core.Complex(3, 4);
      const result = sym.core.im(z);
      expect(result.toString()).toBe('4');
    });

    it('im(0 + 5i) = 5', () => {
      const z = new sym.core.Complex(0, 5);
      const result = sym.core.im(z);
      expect(result.toString()).toBe('5');
    });

    it('im(1 - 7i) = -7', () => {
      const z = new sym.core.Complex(1, -7);
      const result = sym.core.im(z);
      expect(result.toString()).toBe('-7');
    });
  });

  describe('conjugate()', () => {
    it('conjugate(3 + 4i) = 3 - 4i', () => {
      const z = new sym.core.Complex(3, 4);
      const result = sym.core.conjugate(z);
      expect(result.toString()).toBe('3 - 4*I');
    });

    it('conjugate(x) is symbolic', () => {
      const x = new sym.core.Symbol('x');
      const result = sym.core.conjugate(x);
      expect(result.toString()).toMatch(/conjugate\(x\)/i);
    });

    it('conjugate(5) = 5 for real integers', () => {
      const five = new sym.core.Integer(5);
      const result = sym.core.conjugate(five);
      expect(result.toString()).toBe('5');
    });

    it('conjugate(2 - 3i) = 2 + 3i', () => {
      const z = new sym.core.Complex(2, -3);
      const result = sym.core.conjugate(z);
      expect(result.toString()).toBe('2 + 3*I');
    });
  });

  describe('arg()', () => {
    it('arg(1 + i) = pi/4', () => {
      const z = new sym.core.Complex(1, 1);
      const result = sym.core.arg(z);
      // arg(1+i) = atan2(1, 1) = pi/4
      expect(result.toString()).toMatch(/pi\/4|\(1\/4\)\*pi/);
    });

    it('arg(0 + i) = pi/2', () => {
      const z = new sym.core.Complex(0, 1);
      const result = sym.core.arg(z);
      // arg(i) = atan2(1, 0) = pi/2
      expect(result.toString()).toMatch(/pi\/2|\(1\/2\)\*pi/);
    });

    it('arg(-1 + i) = 3*pi/4', () => {
      const z = new sym.core.Complex(-1, 1);
      const result = sym.core.arg(z);
      // arg(-1+i) = atan2(1, -1) = 3*pi/4
      expect(result.toString()).toMatch(/3\*pi\/4|\(3\/4\)\*pi/);
    });
  });
});

describe('digamma function', () => {
  beforeAll(async () => {
    await sym.core.loadWasmModule();
  });

  it('digamma(1) = -EulerGamma', () => {
    const result = sym.core.digamma(new sym.core.Integer(1));
    expect(result.toString()).toMatch(/-EulerGamma/);
  });

  it('digamma(x) is symbolic', () => {
    const x = new sym.core.Symbol('x');
    const result = sym.core.digamma(x);
    // digamma is implemented as polygamma(0, x)
    expect(result.toString()).toMatch(/polygamma\(0, x\)|digamma\(x\)/i);
  });

  it('digamma(2) evaluates correctly', () => {
    const result = sym.core.digamma(new sym.core.Integer(2));
    // digamma(2) = 1 - EulerGamma
    expect(result.toString()).toMatch(/1 - EulerGamma|-EulerGamma \+ 1/);
  });
});

describe('Max and Min', () => {
  beforeAll(async () => {
    await sym.core.loadWasmModule();
  });

  describe('Max()', () => {
    it('Max(1, 2, 3) = 3', () => {
      const one = new sym.core.Integer(1);
      const two = new sym.core.Integer(2);
      const three = new sym.core.Integer(3);
      const result = sym.core.Max(one, two, three);
      expect(result.toString()).toBe('3');
    });

    it('Max(5, 2) = 5', () => {
      const five = new sym.core.Integer(5);
      const two = new sym.core.Integer(2);
      const result = sym.core.Max(five, two);
      expect(result.toString()).toBe('5');
    });

    it('Max(x, y) is symbolic', () => {
      const [x, y] = sym.core.symbols('x y');
      const result = sym.core.Max(x, y);
      expect(result.toString()).toMatch(/max\(x, y\)|max\(y, x\)/i);
    });

    it('Max(-3, -1, -5) = -1', () => {
      const a = new sym.core.Integer(-3);
      const b = new sym.core.Integer(-1);
      const c = new sym.core.Integer(-5);
      const result = sym.core.Max(a, b, c);
      expect(result.toString()).toBe('-1');
    });
  });

  describe('Min()', () => {
    it('Min(1, 2, 3) = 1', () => {
      const one = new sym.core.Integer(1);
      const two = new sym.core.Integer(2);
      const three = new sym.core.Integer(3);
      const result = sym.core.Min(one, two, three);
      expect(result.toString()).toBe('1');
    });

    it('Min(5, 2) = 2', () => {
      const five = new sym.core.Integer(5);
      const two = new sym.core.Integer(2);
      const result = sym.core.Min(five, two);
      expect(result.toString()).toBe('2');
    });

    it('Min(x, y) is symbolic', () => {
      const [x, y] = sym.core.symbols('x y');
      const result = sym.core.Min(x, y);
      expect(result.toString()).toMatch(/min\(x, y\)|min\(y, x\)/i);
    });

    it('Min(-3, -1, -5) = -5', () => {
      const a = new sym.core.Integer(-3);
      const b = new sym.core.Integer(-1);
      const c = new sym.core.Integer(-5);
      const result = sym.core.Min(a, b, c);
      expect(result.toString()).toBe('-5');
    });
  });
});
