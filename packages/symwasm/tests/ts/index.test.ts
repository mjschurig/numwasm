/**
 * Tests for symwasm module exports and stub verification
 *
 * These tests verify that all modules are exported correctly and that
 * unimplemented functions throw NotImplementedError as expected.
 */
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
});

describe('symwasm stubs: core module', () => {
  // Note: Symbol and symbols() are now implemented (Phase 1.3)
  // Tests for them are in symbol.test.ts

  it('Integer constructor throws NotImplementedError', () => {
    expect(() => new sym.core.Integer(42)).toThrow(NotImplementedError);
  });

  it('Rational constructor throws NotImplementedError', () => {
    expect(() => new sym.core.Rational(1, 2)).toThrow(NotImplementedError);
  });

  it('Float constructor throws NotImplementedError', () => {
    expect(() => new sym.core.Float(3.14)).toThrow(NotImplementedError);
  });

  it('sentinel constants exist and have correct string representation', () => {
    expect(sym.core.pi).toBeDefined();
    expect(sym.core.E).toBeDefined();
    expect(sym.core.I).toBeDefined();
    expect(sym.core.oo).toBeDefined();
    expect(sym.core.pi.toString()).toBe('pi');
    expect(sym.core.E.toString()).toBe('E');
  });
});

describe('symwasm stubs: simplify module', () => {
  it('simplify() throws NotImplementedError', () => {
    expect(() => sym.simplify.simplify({} as any)).toThrow(NotImplementedError);
  });
});

describe('symwasm stubs: solvers module', () => {
  it('solve() throws NotImplementedError', () => {
    expect(() => sym.solvers.solve({} as any)).toThrow(NotImplementedError);
  });
});

describe('symwasm stubs: calculus module', () => {
  it('diff() throws NotImplementedError', () => {
    expect(() => sym.calculus.diff({} as any, {} as any)).toThrow(
      NotImplementedError
    );
  });

  it('integrate() throws NotImplementedError', () => {
    expect(() => sym.calculus.integrate({} as any, {} as any)).toThrow(
      NotImplementedError
    );
  });
});

describe('symwasm stubs: matrices module', () => {
  it('Matrix constructor throws NotImplementedError', () => {
    expect(() => new sym.matrices.Matrix([[1, 0], [0, 1]])).toThrow(
      NotImplementedError
    );
  });

  it('eye() throws NotImplementedError', () => {
    expect(() => sym.matrices.eye(3)).toThrow(NotImplementedError);
  });
});

describe('symwasm stubs: printing module', () => {
  it('latex() throws NotImplementedError', () => {
    expect(() => sym.printing.latex({} as any)).toThrow(NotImplementedError);
  });
});
