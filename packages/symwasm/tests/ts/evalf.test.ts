/**
 * Tests for symwasm numerical evaluation (Phase 1.8)
 *
 * Tests Expr.evalf(), evalfNumber(), and evalfComplex()
 */
import { describe, it, expect, beforeAll } from 'vitest';
import * as sym from '../../src/ts/index';

describe('evalf: Numerical Evaluation', () => {
  beforeAll(async () => {
    await sym.core.loadWasmModule();
  });

  describe('Expr.evalf() basic usage', () => {
    it('evaluates pi to double precision', () => {
      const result = sym.core.pi.evalf();
      expect(result.toString()).toMatch(/^3\.14159/);
    });

    it('evaluates E to double precision', () => {
      const result = sym.core.E.evalf();
      expect(result.toString()).toMatch(/^2\.71828/);
    });

    it('evaluates rational to float', () => {
      const half = new sym.core.Rational(1, 2);
      const result = half.evalf();
      expect(result.toString()).toBe('0.5');
    });

    it('evaluates integer (returns same value as float)', () => {
      const five = new sym.core.Integer(5);
      const result = five.evalf();
      // SymEngine may return "5.0" or "5" depending on implementation
      expect(result.toString()).toMatch(/^5(\.0)?$/);
    });

    it('evaluates expression with constants', () => {
      // 2 * pi
      const two = new sym.core.Integer(2);
      const result = sym.core.mul(two, sym.core.pi).evalf();
      expect(result.toString()).toMatch(/^6\.28318/);
    });
  });

  describe('Expr.evalfNumber() convenience', () => {
    it('returns JavaScript number for pi', () => {
      const value = sym.core.pi.evalfNumber();
      expect(value).toBeCloseTo(Math.PI, 10);
    });

    it('returns JavaScript number for E', () => {
      const value = sym.core.E.evalfNumber();
      expect(value).toBeCloseTo(Math.E, 10);
    });

    it('returns JavaScript number for expression', () => {
      // pi + E
      const expr = sym.core.add(sym.core.pi, sym.core.E);
      const value = expr.evalfNumber();
      expect(value).toBeCloseTo(Math.PI + Math.E, 10);
    });

    it('returns JavaScript number for rational', () => {
      const half = new sym.core.Rational(1, 2);
      const value = half.evalfNumber();
      expect(value).toBeCloseTo(0.5, 10);
    });

    it('returns JavaScript number for integer', () => {
      const five = new sym.core.Integer(5);
      const value = five.evalfNumber();
      expect(value).toBe(5);
    });

    it('evaluates EulerGamma', () => {
      const value = sym.core.EulerGamma.evalfNumber();
      expect(value).toBeCloseTo(0.5772156649, 8);
    });

    it('evaluates GoldenRatio', () => {
      const value = sym.core.GoldenRatio.evalfNumber();
      expect(value).toBeCloseTo(1.6180339887, 8);
    });

    it('evaluates Catalan', () => {
      const value = sym.core.Catalan.evalfNumber();
      expect(value).toBeCloseTo(0.9159655941, 8);
    });
  });

  describe('complex evaluation', () => {
    it('evaluates I to complex', () => {
      const result = sym.core.I.evalfComplex();
      expect(result.real).toBeCloseTo(0, 10);
      expect(result.imag).toBeCloseTo(1, 10);
    });

    it('evaluates 1 + 2*I', () => {
      const one = new sym.core.Integer(1);
      const two = new sym.core.Integer(2);
      const twoI = sym.core.mul(two, sym.core.I);
      const expr = sym.core.add(one, twoI);

      const result = expr.evalfComplex();
      expect(result.real).toBeCloseTo(1, 10);
      expect(result.imag).toBeCloseTo(2, 10);
    });

    it('evaluates -3*I', () => {
      const three = new sym.core.Integer(3);
      const threeI = sym.core.mul(three, sym.core.I);
      const negThreeI = sym.core.neg(threeI);

      const result = negThreeI.evalfComplex();
      expect(result.real).toBeCloseTo(0, 10);
      expect(result.imag).toBeCloseTo(-3, 10);
    });

    it('evalfComplex returns real number with imag=0', () => {
      const result = sym.core.pi.evalfComplex();
      expect(result.real).toBeCloseTo(Math.PI, 10);
      expect(result.imag).toBe(0);
    });
  });

  describe('error handling', () => {
    it('evalfNumber throws for complex result', () => {
      expect(() => sym.core.I.evalfNumber()).toThrow(/complex/i);
    });
  });

  describe('precision parameter', () => {
    it('accepts default precision (53 bits)', () => {
      // Without MPFR, only 53-bit precision is supported
      const result = sym.core.pi.evalf();
      expect(result.toString()).toMatch(/^3\.14159/);
    });

    it('accepts explicit 53-bit precision', () => {
      const result = sym.core.pi.evalf(53);
      expect(result.toString()).toMatch(/^3\.14159/);
    });
  });

  describe('arithmetic with evaluated results', () => {
    it('can use evalf result in further calculations', () => {
      const piFloat = sym.core.pi.evalf();
      const two = new sym.core.Integer(2);
      const result = sym.core.mul(two, piFloat);
      // Result should be close to 2*pi
      const value = parseFloat(result.toString());
      expect(value).toBeCloseTo(2 * Math.PI, 8);
    });
  });
});
