/**
 * Tests for symwasm Symbol and symbols() (Phase 1.3: Symbols and Variables)
 *
 * These tests require the WASM module to be loaded.
 */
import { describe, it, expect, beforeAll, afterEach } from 'vitest';
import * as sym from '../../src/ts/index';

describe('core: Symbol creation (Phase 1.3)', () => {
  beforeAll(async () => {
    // Load WASM module before tests
    await sym.core.loadWasmModule();
  });

  describe('Symbol constructor', () => {
    it('creates a symbol with the given name', () => {
      const x = new sym.core.Symbol('x');
      expect(x.name).toBe('x');
      x.free();
    });

    it('toString() returns the symbol name', () => {
      const y = new sym.core.Symbol('y');
      expect(y.toString()).toBe('y');
      y.free();
    });

    it('creates different symbols for different names', () => {
      const x = new sym.core.Symbol('x');
      const y = new sym.core.Symbol('y');
      expect(x.equals(y)).toBe(false);
      x.free();
      y.free();
    });

    it('creates equal symbols for the same name', () => {
      const x1 = new sym.core.Symbol('x');
      const x2 = new sym.core.Symbol('x');
      expect(x1.equals(x2)).toBe(true);
      x1.free();
      x2.free();
    });

    it('has get_type() returning SYMENGINE_SYMBOL', () => {
      const x = new sym.core.Symbol('x');
      expect(x.get_type()).toBe(sym.core.SymEngineTypeID.SYMENGINE_SYMBOL);
      x.free();
    });

    it('has hash() returning a number', () => {
      const x = new sym.core.Symbol('x');
      const hash = x.hash();
      expect(typeof hash).toBe('number');
      x.free();
    });

    it('same symbols have equal hashes', () => {
      const x1 = new sym.core.Symbol('x');
      const x2 = new sym.core.Symbol('x');
      expect(x1.hash()).toBe(x2.hash());
      x1.free();
      x2.free();
    });

    it('free_symbols() returns array containing itself', () => {
      const x = new sym.core.Symbol('x');
      const freeSyms = x.free_symbols();
      expect(freeSyms.length).toBe(1);
      expect(freeSyms[0].name).toBe('x');
      // Clean up
      freeSyms.forEach((s) => s.free());
      x.free();
    });

    it('handles multi-character symbol names', () => {
      const alpha = new sym.core.Symbol('alpha');
      expect(alpha.name).toBe('alpha');
      expect(alpha.toString()).toBe('alpha');
      alpha.free();
    });

    it('handles symbol names with underscores', () => {
      const x_1 = new sym.core.Symbol('x_1');
      expect(x_1.name).toBe('x_1');
      expect(x_1.toString()).toBe('x_1');
      x_1.free();
    });
  });

  describe('symbols() function', () => {
    it('creates multiple symbols from space-separated names', () => {
      const syms = sym.core.symbols('x y z');
      expect(syms.length).toBe(3);
      expect(syms[0].name).toBe('x');
      expect(syms[1].name).toBe('y');
      expect(syms[2].name).toBe('z');
      syms.forEach((s) => s.free());
    });

    it('creates multiple symbols from comma-separated names', () => {
      const syms = sym.core.symbols('a, b, c');
      expect(syms.length).toBe(3);
      expect(syms[0].name).toBe('a');
      expect(syms[1].name).toBe('b');
      expect(syms[2].name).toBe('c');
      syms.forEach((s) => s.free());
    });

    it('handles mixed separators', () => {
      const syms = sym.core.symbols('x,y z');
      expect(syms.length).toBe(3);
      expect(syms[0].name).toBe('x');
      expect(syms[1].name).toBe('y');
      expect(syms[2].name).toBe('z');
      syms.forEach((s) => s.free());
    });

    it('handles extra whitespace', () => {
      const syms = sym.core.symbols('  a   b   c  ');
      expect(syms.length).toBe(3);
      expect(syms[0].name).toBe('a');
      expect(syms[1].name).toBe('b');
      expect(syms[2].name).toBe('c');
      syms.forEach((s) => s.free());
    });

    it('returns empty array for empty string', () => {
      const syms = sym.core.symbols('');
      expect(syms.length).toBe(0);
    });

    it('returns empty array for whitespace-only string', () => {
      const syms = sym.core.symbols('   ');
      expect(syms.length).toBe(0);
    });

    it('creates a single symbol for a single name', () => {
      const syms = sym.core.symbols('x');
      expect(syms.length).toBe(1);
      expect(syms[0].name).toBe('x');
      syms[0].free();
    });

    it('allows destructuring assignment', () => {
      const [x, y] = sym.core.symbols('x y');
      expect(x.name).toBe('x');
      expect(y.name).toBe('y');
      x.free();
      y.free();
    });
  });

  describe('Symbol.free() method', () => {
    it('frees the underlying WASM memory', () => {
      const x = new sym.core.Symbol('x');
      expect(() => x.toString()).not.toThrow(); // Works before free
      x.free();
      // After free, the fallback string should be used
      expect(x.toString()).toBe('x');
    });

    it('can be called multiple times safely', () => {
      const x = new sym.core.Symbol('x');
      x.free();
      expect(() => x.free()).not.toThrow();
    });
  });
});
