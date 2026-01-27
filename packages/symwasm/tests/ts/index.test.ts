import { describe, it, expect } from 'vitest';
import * as sym from '../../src/ts/index';
import { NotImplementedError } from '../../src/ts/errors';

describe('symwasm module exports', () => {
  it('exports all 6 modules', () => {
    expect(sym.core).toBeDefined();
    expect(sym.simplify).toBeDefined();
    expect(sym.solvers).toBeDefined();
    expect(sym.calculus).toBeDefined();
    expect(sym.matrices).toBeDefined();
    expect(sym.printing).toBeDefined();
  });

  it('exports NotImplementedError', () => {
    expect(sym.NotImplementedError).toBeDefined();
  });

  it('core.Symbol constructor throws NotImplementedError', () => {
    expect(() => new sym.core.Symbol('x')).toThrow(NotImplementedError);
  });

  it('core.symbols throws NotImplementedError', () => {
    expect(() => sym.core.symbols('x y z')).toThrow(NotImplementedError);
  });

  it('core.Integer constructor throws NotImplementedError', () => {
    expect(() => new sym.core.Integer(42)).toThrow(NotImplementedError);
  });

  it('core sentinel constants exist', () => {
    expect(sym.core.pi).toBeDefined();
    expect(sym.core.E).toBeDefined();
    expect(sym.core.I).toBeDefined();
    expect(sym.core.oo).toBeDefined();
    expect(sym.core.pi.toString()).toBe('pi');
    expect(sym.core.E.toString()).toBe('E');
  });

  it('simplify.simplify throws NotImplementedError', () => {
    expect(() => sym.simplify.simplify({} as any)).toThrow(NotImplementedError);
  });

  it('solvers.solve throws NotImplementedError', () => {
    expect(() => sym.solvers.solve({} as any)).toThrow(NotImplementedError);
  });

  it('calculus.diff throws NotImplementedError', () => {
    expect(() => sym.calculus.diff({} as any, {} as any)).toThrow(NotImplementedError);
  });

  it('calculus.integrate throws NotImplementedError', () => {
    expect(() => sym.calculus.integrate({} as any, {} as any)).toThrow(NotImplementedError);
  });

  it('matrices.Matrix constructor throws NotImplementedError', () => {
    expect(() => new sym.matrices.Matrix([[1, 0], [0, 1]])).toThrow(NotImplementedError);
  });

  it('matrices.eye throws NotImplementedError', () => {
    expect(() => sym.matrices.eye(3)).toThrow(NotImplementedError);
  });

  it('printing.latex throws NotImplementedError', () => {
    expect(() => sym.printing.latex({} as any)).toThrow(NotImplementedError);
  });
});
