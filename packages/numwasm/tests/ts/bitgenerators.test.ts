/**
 * Tests for Phase 27: Additional BitGenerators (MT19937, Philox, SFC64)
 */

import { describe, it, expect, beforeAll } from 'vitest';
import {
  loadWasmModule,
  MT19937,
  Philox,
  SFC64,
  PCG64,
  SeedSequence,
  Generator,
  default_rng,
  getBitGenerator,
  listBitGenerators,
} from '../../dist/numjs.mjs';

describe('Phase 27: Additional BitGenerators', () => {
  beforeAll(async () => {
    await loadWasmModule();
  });

  describe('SFC64 BitGenerator', () => {
    it('creates with integer seed', () => {
      const sfc = new SFC64(12345);
      expect(sfc).toBeInstanceOf(SFC64);
      sfc.dispose();
    });

    it('generates uint64 values', () => {
      const sfc = new SFC64(42);
      const val = sfc.next_uint64();
      expect(typeof val).toBe('bigint');
      expect(val).toBeGreaterThanOrEqual(0n);
      sfc.dispose();
    });

    it('generates double values in [0, 1)', () => {
      const sfc = new SFC64(42);
      for (let i = 0; i < 100; i++) {
        const val = sfc.next_double();
        expect(val).toBeGreaterThanOrEqual(0);
        expect(val).toBeLessThan(1);
      }
      sfc.dispose();
    });

    it('produces reproducible sequences', () => {
      const sfc1 = new SFC64(12345);
      const sfc2 = new SFC64(12345);
      for (let i = 0; i < 10; i++) {
        expect(sfc1.next_uint64()).toBe(sfc2.next_uint64());
      }
      sfc1.dispose();
      sfc2.dispose();
    });

    it('state serialization roundtrip', () => {
      const sfc = new SFC64(12345);
      for (let i = 0; i < 5; i++) sfc.next_uint64();
      const state = sfc.getState();
      expect(state.bit_generator).toBe('SFC64');
      const expected: bigint[] = [];
      for (let i = 0; i < 5; i++) expected.push(sfc.next_uint64());
      sfc.setState(state);
      for (let i = 0; i < 5; i++) expect(sfc.next_uint64()).toBe(expected[i]);
      sfc.dispose();
    });

    it('spawns independent generators', () => {
      const sfc = new SFC64(12345);
      const [child1, child2] = sfc.spawn(2);
      expect(child1.next_uint64()).not.toBe(child2.next_uint64());
      child1.dispose();
      child2.dispose();
      sfc.dispose();
    });
  });

  describe('MT19937 BitGenerator', () => {
    it('creates with integer seed', () => {
      const mt = new MT19937(12345);
      expect(mt).toBeInstanceOf(MT19937);
      mt.dispose();
    });

    it('generates uint32 values (native output)', () => {
      const mt = new MT19937(42);
      const val = mt.next_uint32();
      expect(typeof val).toBe('number');
      expect(val).toBeGreaterThanOrEqual(0);
      expect(val).toBeLessThanOrEqual(0xFFFFFFFF);
      mt.dispose();
    });

    it('generates double values in [0, 1)', () => {
      const mt = new MT19937(42);
      for (let i = 0; i < 100; i++) {
        const val = mt.next_double();
        expect(val).toBeGreaterThanOrEqual(0);
        expect(val).toBeLessThan(1);
      }
      mt.dispose();
    });

    it('produces reproducible sequences', () => {
      const mt1 = new MT19937(12345);
      const mt2 = new MT19937(12345);
      for (let i = 0; i < 10; i++) {
        expect(mt1.next_uint32()).toBe(mt2.next_uint32());
      }
      mt1.dispose();
      mt2.dispose();
    });

    it('state serialization roundtrip', () => {
      const mt = new MT19937(12345);
      for (let i = 0; i < 100; i++) mt.next_uint32();
      const state = mt.getState();
      expect(state.bit_generator).toBe('MT19937');
      expect(state.state.key.length).toBe(624);
      const expected: number[] = [];
      for (let i = 0; i < 10; i++) expected.push(mt.next_uint32());
      mt.setState(state);
      for (let i = 0; i < 10; i++) expect(mt.next_uint32()).toBe(expected[i]);
      mt.dispose();
    });
  });

  describe('Philox BitGenerator', () => {
    it('creates with integer seed', () => {
      const philox = new Philox(12345);
      expect(philox).toBeInstanceOf(Philox);
      philox.dispose();
    });

    it('generates uint64 values', () => {
      const philox = new Philox(42);
      const val = philox.next_uint64();
      expect(typeof val).toBe('bigint');
      expect(val).toBeGreaterThanOrEqual(0n);
      philox.dispose();
    });

    it('generates double values in [0, 1)', () => {
      const philox = new Philox(42);
      for (let i = 0; i < 100; i++) {
        const val = philox.next_double();
        expect(val).toBeGreaterThanOrEqual(0);
        expect(val).toBeLessThan(1);
      }
      philox.dispose();
    });

    it('produces reproducible sequences', () => {
      const p1 = new Philox(12345);
      const p2 = new Philox(12345);
      for (let i = 0; i < 10; i++) {
        expect(p1.next_uint64()).toBe(p2.next_uint64());
      }
      p1.dispose();
      p2.dispose();
    });

    it('jump advances generator', () => {
      const p1 = new Philox(12345);
      const p2 = new Philox(12345);
      p2.jump();
      expect(p1.next_uint64()).not.toBe(p2.next_uint64());
      p1.dispose();
      p2.dispose();
    });

    it('advance works correctly', () => {
      const p1 = new Philox(12345);
      const p2 = new Philox(12345);
      // Philox produces 4 values per counter increment
      // Advance counter by 25 = 100 random values
      p1.advance(25n);
      for (let i = 0; i < 100; i++) p2.next_uint64();
      expect(p1.next_uint64()).toBe(p2.next_uint64());
      p1.dispose();
      p2.dispose();
    });

    it('state serialization roundtrip', () => {
      const philox = new Philox(12345);
      for (let i = 0; i < 10; i++) philox.next_uint64();
      const state = philox.getState();
      expect(state.bit_generator).toBe('Philox');
      const expected: bigint[] = [];
      for (let i = 0; i < 5; i++) expected.push(philox.next_uint64());
      philox.setState(state);
      for (let i = 0; i < 5; i++) expect(philox.next_uint64()).toBe(expected[i]);
      philox.dispose();
    });
  });

  describe('BitGenerator Registry', () => {
    it('lists all available BitGenerators', () => {
      const generators = listBitGenerators();
      expect(generators).toContain('PCG64');
      expect(generators).toContain('MT19937');
      expect(generators).toContain('PHILOX');
      expect(generators).toContain('SFC64');
    });

    it('getBitGenerator returns correct class', () => {
      const PCG = getBitGenerator('PCG64');
      const MT = getBitGenerator('MT19937');
      expect(new PCG(42)).toBeInstanceOf(PCG64);
      expect(new MT(42)).toBeInstanceOf(MT19937);
    });

    it('getBitGenerator is case-insensitive', () => {
      expect(getBitGenerator('pcg64')).toBe(getBitGenerator('PCG64'));
      expect(getBitGenerator('mt19937')).toBe(getBitGenerator('MT19937'));
    });

    it('getBitGenerator throws for unknown generator', () => {
      expect(() => getBitGenerator('unknown')).toThrow(/Unknown BitGenerator/);
    });
  });

  describe('default_rng with BitGenerator option', () => {
    it('creates Generator with PCG64 by default', () => {
      const rng = default_rng(12345);
      expect(rng.bitGenerator).toBeInstanceOf(PCG64);
      rng.dispose();
    });

    it('creates Generator with MT19937', () => {
      const rng = default_rng(12345, { bitGenerator: 'MT19937' });
      expect(rng.bitGenerator).toBeInstanceOf(MT19937);
      rng.dispose();
    });

    it('creates Generator with Philox', () => {
      const rng = default_rng(12345, { bitGenerator: 'Philox' });
      expect(rng.bitGenerator).toBeInstanceOf(Philox);
      rng.dispose();
    });

    it('creates Generator with SFC64', () => {
      const rng = default_rng(12345, { bitGenerator: 'SFC64' });
      expect(rng.bitGenerator).toBeInstanceOf(SFC64);
      rng.dispose();
    });
  });

  describe('Generator integration', () => {
    it('works with MT19937', () => {
      const rng = new Generator(new MT19937(12345));
      const arr = rng.random(10);
      expect(arr.size).toBe(10);
      rng.dispose();
    });

    it('works with Philox', () => {
      const rng = new Generator(new Philox(12345));
      const arr = rng.random(10);
      expect(arr.size).toBe(10);
      rng.dispose();
    });

    it('works with SFC64', () => {
      const rng = new Generator(new SFC64(12345));
      const arr = rng.random(10);
      expect(arr.size).toBe(10);
      rng.dispose();
    });
  });
});
