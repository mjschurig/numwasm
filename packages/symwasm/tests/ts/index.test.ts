/**
 * Tests for symwasm module exports and stub verification
 *
 * These tests verify that all modules are exported correctly and that
 * unimplemented functions throw NotImplementedError as expected.
 */
import { describe, it, expect, beforeAll } from 'vitest';
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
  // Note: Symbol, symbols(), Integer, Rational, Float are now implemented (Phase 1.3-1.4)
  // Tests for them are in symbol.test.ts and numbers.test.ts
  // Note: Constants (pi, E, I, oo) are now WASM-backed (Phase 1.6)
  // Tests for them are in constants.test.ts

  beforeAll(async () => {
    await sym.core.loadWasmModule();
  });

  it('WASM-backed constants exist and have correct string representation', () => {
    expect(sym.core.pi).toBeDefined();
    expect(sym.core.E).toBeDefined();
    expect(sym.core.I).toBeDefined();
    expect(sym.core.oo).toBeDefined();
    expect(sym.core.pi.toString()).toBe('pi');
    expect(sym.core.E.toString()).toBe('E');
  });
});

describe('symwasm: simplify module', () => {
  // Note: simplify(), expand(), numer(), denom(), and rewrite functions are now implemented (Phase 2.4)
  // factor(), collect(), cancel() were removed (no SymEngine support)
  // Tests for implemented functions are in simplify.test.ts
  it('simplify module is exported', () => {
    expect(sym.simplify).toBeDefined();
    expect(typeof sym.simplify.expand).toBe('function');
    expect(typeof sym.simplify.simplify).toBe('function');
  });
});

describe('symwasm stubs: solvers module', () => {
  it('solve() throws NotImplementedError', () => {
    expect(() => sym.solvers.solve({} as any)).toThrow(NotImplementedError);
  });
});

describe('symwasm stubs: calculus module', () => {
  // Note: diff() is now implemented (Phase 2.2)
  // Tests for it are in calculus.test.ts

  it('integrate() throws NotImplementedError', () => {
    expect(() => sym.calculus.integrate({} as any, {} as any)).toThrow(
      NotImplementedError
    );
  });
});

describe('symwasm: matrices module', () => {
  // Note: Matrix, eye, zeros, ones, diag are now implemented (Phase 3.1)
  // det(), inv(), transpose() still throw NotImplementedError
  // Tests for implemented functions are in matrices.test.ts
  it('matrices module is exported', () => {
    expect(sym.matrices).toBeDefined();
    expect(typeof sym.matrices.Matrix).toBe('function');
    expect(typeof sym.matrices.eye).toBe('function');
    expect(typeof sym.matrices.zeros).toBe('function');
  });
});

describe('symwasm stubs: printing module', () => {
  it('latex() throws NotImplementedError', () => {
    expect(() => sym.printing.latex({} as any)).toThrow(NotImplementedError);
  });
});
