/**
 * Tests for symwasm constants (Phase 1.6)
 *
 * Tests pi, E, I, oo, EulerGamma, Catalan, GoldenRatio
 * and S.Infinity, S.NegativeInfinity, S.ComplexInfinity, S.NaN
 */
import { describe, it, expect, beforeAll } from 'vitest';
import * as sym from '../../src/ts/index';

describe('constants: pi', () => {
  beforeAll(async () => {
    await sym.core.loadWasmModule();
  });

  it('has correct string representation', () => {
    expect(sym.core.pi.toString()).toBe('pi');
  });

  it('has correct type', () => {
    expect(sym.core.pi.get_type()).toBe(sym.core.SymEngineTypeID.SYMENGINE_CONSTANT);
  });

  it('can be used in arithmetic', () => {
    const two = new sym.core.Integer(2);
    const result = sym.core.mul(two, sym.core.pi);
    expect(result.toString()).toBe('2*pi');
  });

  it('has consistent hash', () => {
    const h1 = sym.core.pi.hash();
    const h2 = sym.core.pi.hash();
    expect(h1).toBe(h2);
  });

  it('equals itself', () => {
    expect(sym.core.pi.equals(sym.core.pi)).toBe(true);
  });

  it('has no free symbols', () => {
    const freeSyms = sym.core.pi.free_symbols();
    expect(freeSyms.length).toBe(0);
  });
});

describe('constants: E (Euler\'s number)', () => {
  beforeAll(async () => {
    await sym.core.loadWasmModule();
  });

  it('has correct string representation', () => {
    expect(sym.core.E.toString()).toBe('E');
  });

  it('has correct type', () => {
    expect(sym.core.E.get_type()).toBe(sym.core.SymEngineTypeID.SYMENGINE_CONSTANT);
  });

  it('can be used in arithmetic', () => {
    const x = new sym.core.Symbol('x');
    const result = sym.core.pow(sym.core.E, x);
    // SymEngine represents e^x as exp(x)
    expect(result.toString()).toBe('exp(x)');
  });
});

describe('constants: I (imaginary unit)', () => {
  beforeAll(async () => {
    await sym.core.loadWasmModule();
  });

  it('has correct string representation', () => {
    expect(sym.core.I.toString()).toBe('I');
  });

  it('I*I equals -1', () => {
    const result = sym.core.mul(sym.core.I, sym.core.I);
    expect(result.toString()).toBe('-1');
  });

  it('can create complex numbers', () => {
    const two = new sym.core.Integer(2);
    const three = new sym.core.Integer(3);
    // 2 + 3*I
    const threeI = sym.core.mul(three, sym.core.I);
    const result = sym.core.add(two, threeI);
    expect(result.toString()).toBe('2 + 3*I');
  });
});

describe('constants: oo (infinity)', () => {
  beforeAll(async () => {
    await sym.core.loadWasmModule();
  });

  it('has correct string representation', () => {
    expect(sym.core.oo.toString()).toBe('oo');
  });

  it('has correct type', () => {
    expect(sym.core.oo.get_type()).toBe(sym.core.SymEngineTypeID.SYMENGINE_INFTY);
  });

  it('oo + 1 equals oo', () => {
    const one = new sym.core.Integer(1);
    const result = sym.core.add(sym.core.oo, one);
    expect(result.toString()).toBe('oo');
  });
});

describe('constants: S.Infinity variants', () => {
  beforeAll(async () => {
    await sym.core.loadWasmModule();
  });

  it('S.Infinity has correct string representation', () => {
    expect(sym.core.S.Infinity.toString()).toBe('oo');
  });

  it('S.NegativeInfinity has correct string representation', () => {
    expect(sym.core.S.NegativeInfinity.toString()).toBe('-oo');
  });

  it('S.ComplexInfinity has correct string representation', () => {
    expect(sym.core.S.ComplexInfinity.toString()).toBe('zoo');
  });

  it('S.NaN has correct string representation', () => {
    expect(sym.core.S.NaN.toString()).toBe('nan');
  });

  it('S.NaN has correct type', () => {
    expect(sym.core.S.NaN.get_type()).toBe(sym.core.SymEngineTypeID.SYMENGINE_NOT_A_NUMBER);
  });
});

describe('constants: additional constants', () => {
  beforeAll(async () => {
    await sym.core.loadWasmModule();
  });

  it('EulerGamma has correct string representation', () => {
    expect(sym.core.EulerGamma.toString()).toBe('EulerGamma');
  });

  it('EulerGamma has correct type', () => {
    expect(sym.core.EulerGamma.get_type()).toBe(sym.core.SymEngineTypeID.SYMENGINE_CONSTANT);
  });

  it('Catalan has correct string representation', () => {
    expect(sym.core.Catalan.toString()).toBe('Catalan');
  });

  it('GoldenRatio has correct string representation', () => {
    expect(sym.core.GoldenRatio.toString()).toBe('GoldenRatio');
  });
});

describe('constants: arithmetic with multiple constants', () => {
  beforeAll(async () => {
    await sym.core.loadWasmModule();
  });

  it('can add pi and E', () => {
    const result = sym.core.add(sym.core.pi, sym.core.E);
    expect(result.toString()).toBe('E + pi');
  });

  it('can multiply constants by symbols', () => {
    const x = new sym.core.Symbol('x');
    const result = sym.core.mul(sym.core.pi, x);
    // SymEngine may order factors differently (alphabetically)
    expect(result.toString()).toBe('x*pi');
  });

  it('can create e^(i*pi)', () => {
    // e^(i*pi) = -1 (Euler's identity)
    const iPi = sym.core.mul(sym.core.I, sym.core.pi);
    const result = sym.core.pow(sym.core.E, iPi);
    // SymEngine uses exp() notation and may not simplify to -1
    expect(result.toString()).toMatch(/exp\(I\*pi\)|-1/);
  });
});
