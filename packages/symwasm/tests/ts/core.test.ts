/**
 * Tests for symwasm core module (Phase 1.2: Core Base Classes)
 *
 * Tests Expr base class methods: equals, hash, get_type, free_symbols, toString, free
 */
import { describe, it, expect, beforeAll } from 'vitest';
import * as sym from '../../src/ts/index';

describe('core: Expr base class', () => {
  beforeAll(async () => {
    await sym.core.loadWasmModule();
  });

  describe('WASM-backed constants (pi, E, I, oo)', () => {
    it('toString() returns correct string representation', () => {
      expect(sym.core.pi.toString()).toBe('pi');
      expect(sym.core.E.toString()).toBe('E');
      expect(sym.core.I.toString()).toBe('I');
      expect(sym.core.oo.toString()).toBe('oo');
    });

    it('equals() returns true for same constant', () => {
      expect(sym.core.pi.equals(sym.core.pi)).toBe(true);
      expect(sym.core.E.equals(sym.core.E)).toBe(true);
      expect(sym.core.I.equals(sym.core.I)).toBe(true);
      expect(sym.core.oo.equals(sym.core.oo)).toBe(true);
    });

    it('equals() returns false for different constants', () => {
      expect(sym.core.pi.equals(sym.core.E)).toBe(false);
      expect(sym.core.E.equals(sym.core.I)).toBe(false);
      expect(sym.core.I.equals(sym.core.oo)).toBe(false);
      expect(sym.core.oo.equals(sym.core.pi)).toBe(false);
    });

    it('free_symbols() returns empty array for constants', () => {
      expect(sym.core.pi.free_symbols()).toEqual([]);
      expect(sym.core.E.free_symbols()).toEqual([]);
      expect(sym.core.I.free_symbols()).toEqual([]);
      expect(sym.core.oo.free_symbols()).toEqual([]);
    });

    it('hash() returns consistent value for constants', () => {
      const h1 = sym.core.pi.hash();
      const h2 = sym.core.pi.hash();
      expect(h1).toBe(h2);
    });

    it('get_type() returns correct type for constants', () => {
      expect(sym.core.pi.get_type()).toBe(sym.core.SymEngineTypeID.SYMENGINE_CONSTANT);
      expect(sym.core.E.get_type()).toBe(sym.core.SymEngineTypeID.SYMENGINE_CONSTANT);
      expect(sym.core.oo.get_type()).toBe(sym.core.SymEngineTypeID.SYMENGINE_INFTY);
    });
  });

  describe('SymEngineTypeID enum', () => {
    it('exports SymEngineTypeID enum with correct values', () => {
      expect(sym.core.SymEngineTypeID).toBeDefined();
      // Note: Values depend on compile-time configuration
      // These are empirically determined from the WASM build
      expect(sym.core.SymEngineTypeID.SYMENGINE_INTEGER).toBe(0);
      expect(sym.core.SymEngineTypeID.SYMENGINE_RATIONAL).toBe(1);
      expect(sym.core.SymEngineTypeID.SYMENGINE_COMPLEX).toBe(2);
      expect(sym.core.SymEngineTypeID.SYMENGINE_REAL_DOUBLE).toBe(6);
      expect(sym.core.SymEngineTypeID.SYMENGINE_SYMBOL).toBe(13);
      expect(sym.core.SymEngineTypeID.SYMENGINE_ADD).toBe(16);
      expect(sym.core.SymEngineTypeID.SYMENGINE_MUL).toBe(15);
      expect(sym.core.SymEngineTypeID.SYMENGINE_POW).toBe(17);
    });
  });

  describe('Expr.free() method', () => {
    it('free() is callable on constants (no-op for singleton constants)', () => {
      // Constants should not throw when free() is called
      const piCopy = sym.core.pi;
      expect(() => piCopy.free()).not.toThrow();
    });
  });
});
