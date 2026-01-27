/**
 * Constants Tests
 *
 * Tests for mathematical constants, special values, and type information classes.
 */

import { describe, it, expect } from 'vitest';
import {
  // Mathematical constants
  e,
  pi,
  euler_gamma,
  // Special floating-point values
  inf,
  PINF,
  NINF,
  nan,
  NAN,
  PZERO,
  NZERO,
  // Indexing constants
  newaxis,
  // Type information
  finfo,
  iinfo,
  DType,
} from '../../dist/numjs.mjs';

describe('Mathematical Constants', () => {
  it('e equals Euler\'s number', () => {
    expect(e).toBeCloseTo(2.718281828459045, 15);
    expect(e).toBe(Math.E);
  });

  it('pi equals π', () => {
    expect(pi).toBeCloseTo(3.141592653589793, 15);
    expect(pi).toBe(Math.PI);
  });

  it('euler_gamma equals Euler-Mascheroni constant', () => {
    expect(euler_gamma).toBeCloseTo(0.5772156649015329, 15);
  });

  it('constants can be used in calculations', () => {
    // e^1 = e
    expect(Math.exp(1)).toBeCloseTo(e, 15);

    // sin(pi) ≈ 0
    expect(Math.sin(pi)).toBeCloseTo(0, 10);

    // cos(pi) = -1
    expect(Math.cos(pi)).toBeCloseTo(-1, 15);

    // e^(i*pi) + 1 = 0 (Euler's identity, using real parts)
    expect(Math.cos(pi)).toBeCloseTo(-1, 15); // Real part of e^(i*pi)
  });
});

describe('Special Floating-Point Values', () => {
  describe('inf and PINF', () => {
    it('inf is positive infinity', () => {
      expect(inf).toBe(Infinity);
      expect(PINF).toBe(Infinity);
      expect(inf).toBe(PINF);
    });

    it('inf is not finite', () => {
      expect(Number.isFinite(inf)).toBe(false);
      expect(Number.isFinite(PINF)).toBe(false);
    });

    it('inf is greater than any finite number', () => {
      expect(inf > Number.MAX_VALUE).toBe(true);
      expect(inf > 1e308).toBe(true);
    });

    it('inf arithmetic works correctly', () => {
      expect(inf + 1).toBe(inf);
      expect(inf * 2).toBe(inf);
      expect(inf - 1000).toBe(inf);
      expect(1 / inf).toBe(0);
    });

    it('1/0 equals inf', () => {
      expect(1 / 0).toBe(inf);
    });
  });

  describe('NINF', () => {
    it('NINF is negative infinity', () => {
      expect(NINF).toBe(-Infinity);
    });

    it('NINF is less than any finite number', () => {
      expect(NINF < -Number.MAX_VALUE).toBe(true);
      expect(NINF < -1e308).toBe(true);
    });

    it('-1/0 equals NINF', () => {
      expect(-1 / 0).toBe(NINF);
    });
  });

  describe('nan and NAN', () => {
    it('nan is NaN', () => {
      expect(Number.isNaN(nan)).toBe(true);
      expect(Number.isNaN(NAN)).toBe(true);
    });

    it('nan is not equal to itself', () => {
      expect(nan === nan).toBe(false);
      expect(NAN === NAN).toBe(false);
      expect(nan === NAN).toBe(false);
    });

    it('nan comparisons are all false', () => {
      expect(nan > 0).toBe(false);
      expect(nan < 0).toBe(false);
      expect(nan >= 0).toBe(false);
      expect(nan <= 0).toBe(false);
      expect(nan === 0).toBe(false);
    });

    it('nan arithmetic produces nan', () => {
      expect(Number.isNaN(nan + 1)).toBe(true);
      expect(Number.isNaN(nan * 2)).toBe(true);
      expect(Number.isNaN(nan - nan)).toBe(true);
    });

    it('0/0 produces nan', () => {
      expect(Number.isNaN(0 / 0)).toBe(true);
    });
  });

  describe('PZERO and NZERO', () => {
    it('PZERO is positive zero', () => {
      expect(PZERO).toBe(0);
      expect(Object.is(PZERO, 0)).toBe(true);
      expect(Object.is(PZERO, +0)).toBe(true);
    });

    it('NZERO is negative zero', () => {
      expect(NZERO === 0).toBe(true); // Compares equal with ===
      expect(Object.is(NZERO, -0)).toBe(true);
    });

    it('PZERO and NZERO compare equal', () => {
      expect(PZERO === NZERO).toBe(true);
      expect(PZERO == NZERO).toBe(true);
    });

    it('PZERO and NZERO are distinguishable', () => {
      expect(Object.is(PZERO, NZERO)).toBe(false);
    });

    it('division by signed zeros produces different infinities', () => {
      expect(1 / PZERO).toBe(Infinity);
      expect(1 / NZERO).toBe(-Infinity);
    });
  });
});

describe('newaxis', () => {
  it('newaxis is a unique symbol', () => {
    expect(typeof newaxis).toBe('symbol');
  });

  it('newaxis has descriptive string', () => {
    expect(newaxis.toString()).toBe('Symbol(newaxis)');
  });

  it('newaxis is unique', () => {
    expect(newaxis).not.toBe(Symbol('newaxis'));
  });
});

describe('finfo', () => {
  describe('Float64', () => {
    it('provides correct Float64 limits', () => {
      const f64 = new finfo(DType.Float64);
      expect(f64.bits).toBe(64);
      expect(f64.dtype).toBe(DType.Float64);
    });

    it('eps is approximately 2.22e-16', () => {
      const f64 = new finfo(DType.Float64);
      expect(f64.eps).toBeCloseTo(2.220446049250313e-16, 25);
      expect(f64.eps).toBe(Number.EPSILON);
    });

    it('max is approximately 1.8e308', () => {
      const f64 = new finfo(DType.Float64);
      expect(f64.max).toBeCloseTo(1.7976931348623157e308, 290);
      expect(f64.max).toBe(Number.MAX_VALUE);
    });

    it('min is approximately 2.2e-308', () => {
      const f64 = new finfo(DType.Float64);
      expect(f64.min).toBeCloseTo(2.2250738585072014e-308, 320);
    });

    it('has correct exponent and mantissa bits', () => {
      const f64 = new finfo(DType.Float64);
      expect(f64.iexp).toBe(11);
      expect(f64.nmant).toBe(52);
      expect(f64.precision).toBe(15);
    });
  });

  describe('Float32', () => {
    it('provides correct Float32 limits', () => {
      const f32 = new finfo(DType.Float32);
      expect(f32.bits).toBe(32);
      expect(f32.dtype).toBe(DType.Float32);
    });

    it('eps is approximately 1.19e-7', () => {
      const f32 = new finfo(DType.Float32);
      expect(f32.eps).toBeCloseTo(1.1920928955078125e-7, 15);
    });

    it('max is approximately 3.4e38', () => {
      const f32 = new finfo(DType.Float32);
      expect(f32.max).toBeCloseTo(3.4028234663852886e38, 30);
    });

    it('has correct exponent and mantissa bits', () => {
      const f32 = new finfo(DType.Float32);
      expect(f32.iexp).toBe(8);
      expect(f32.nmant).toBe(23);
      expect(f32.precision).toBe(6);
    });
  });

  describe('Float16', () => {
    it('provides correct Float16 limits', () => {
      const f16 = new finfo(DType.Float16);
      expect(f16.bits).toBe(16);
      expect(f16.dtype).toBe(DType.Float16);
    });

    it('max is 65504', () => {
      const f16 = new finfo(DType.Float16);
      expect(f16.max).toBe(65504);
    });

    it('has correct exponent and mantissa bits', () => {
      const f16 = new finfo(DType.Float16);
      expect(f16.iexp).toBe(5);
      expect(f16.nmant).toBe(10);
      expect(f16.precision).toBe(3);
    });
  });

  describe('String dtype support', () => {
    it('accepts "float64" string', () => {
      const f64 = new finfo('float64');
      expect(f64.bits).toBe(64);
      expect(f64.dtype).toBe(DType.Float64);
    });

    it('accepts "float32" string', () => {
      const f32 = new finfo('float32');
      expect(f32.bits).toBe(32);
      expect(f32.dtype).toBe(DType.Float32);
    });

    it('accepts "float16" string', () => {
      const f16 = new finfo('float16');
      expect(f16.bits).toBe(16);
      expect(f16.dtype).toBe(DType.Float16);
    });
  });

  describe('Error handling', () => {
    it('throws for integer dtype', () => {
      expect(() => new finfo(DType.Int32)).toThrow('finfo not available');
    });

    it('throws for unknown string', () => {
      expect(() => new finfo('int32' as any)).toThrow('Unknown float dtype');
    });
  });

  describe('toString', () => {
    it('returns readable representation', () => {
      const f64 = new finfo(DType.Float64);
      const str = f64.toString();
      expect(str).toContain('finfo');
      expect(str).toContain('bits=64');
    });
  });
});

describe('iinfo', () => {
  describe('Int32', () => {
    it('provides correct Int32 limits', () => {
      const i32 = new iinfo(DType.Int32);
      expect(i32.bits).toBe(32);
      expect(i32.min).toBe(-2147483648);
      expect(i32.max).toBe(2147483647);
      expect(i32.dtype).toBe(DType.Int32);
    });
  });

  describe('Int8', () => {
    it('provides correct Int8 limits', () => {
      const i8 = new iinfo(DType.Int8);
      expect(i8.bits).toBe(8);
      expect(i8.min).toBe(-128);
      expect(i8.max).toBe(127);
    });
  });

  describe('Int16', () => {
    it('provides correct Int16 limits', () => {
      const i16 = new iinfo(DType.Int16);
      expect(i16.bits).toBe(16);
      expect(i16.min).toBe(-32768);
      expect(i16.max).toBe(32767);
    });
  });

  describe('Int64', () => {
    it('provides correct Int64 limits', () => {
      const i64 = new iinfo(DType.Int64);
      expect(i64.bits).toBe(64);
      expect(i64.min).toBe(-9223372036854775808);
      expect(i64.max).toBe(9223372036854775807);
    });
  });

  describe('Uint8', () => {
    it('provides correct Uint8 limits', () => {
      const u8 = new iinfo(DType.Uint8);
      expect(u8.bits).toBe(8);
      expect(u8.min).toBe(0);
      expect(u8.max).toBe(255);
    });
  });

  describe('Uint16', () => {
    it('provides correct Uint16 limits', () => {
      const u16 = new iinfo(DType.Uint16);
      expect(u16.bits).toBe(16);
      expect(u16.min).toBe(0);
      expect(u16.max).toBe(65535);
    });
  });

  describe('Uint32', () => {
    it('provides correct Uint32 limits', () => {
      const u32 = new iinfo(DType.Uint32);
      expect(u32.bits).toBe(32);
      expect(u32.min).toBe(0);
      expect(u32.max).toBe(4294967295);
    });
  });

  describe('Uint64', () => {
    it('provides correct Uint64 limits', () => {
      const u64 = new iinfo(DType.Uint64);
      expect(u64.bits).toBe(64);
      expect(u64.min).toBe(0);
      expect(u64.max).toBe(18446744073709551615);
    });
  });

  describe('Bool', () => {
    it('provides correct Bool limits', () => {
      const b = new iinfo(DType.Bool);
      expect(b.bits).toBe(8);
      expect(b.min).toBe(0);
      expect(b.max).toBe(1);
    });
  });

  describe('String dtype support', () => {
    it('accepts "int32" string', () => {
      const i32 = new iinfo('int32');
      expect(i32.bits).toBe(32);
      expect(i32.dtype).toBe(DType.Int32);
    });

    it('accepts "uint8" string', () => {
      const u8 = new iinfo('uint8');
      expect(u8.bits).toBe(8);
      expect(u8.min).toBe(0);
      expect(u8.max).toBe(255);
    });

    it('accepts all string dtypes', () => {
      expect(new iinfo('bool').bits).toBe(8);
      expect(new iinfo('int8').bits).toBe(8);
      expect(new iinfo('int16').bits).toBe(16);
      expect(new iinfo('int32').bits).toBe(32);
      expect(new iinfo('int64').bits).toBe(64);
      expect(new iinfo('uint8').bits).toBe(8);
      expect(new iinfo('uint16').bits).toBe(16);
      expect(new iinfo('uint32').bits).toBe(32);
      expect(new iinfo('uint64').bits).toBe(64);
    });
  });

  describe('Error handling', () => {
    it('throws for float dtype', () => {
      expect(() => new iinfo(DType.Float32)).toThrow('iinfo not available');
    });

    it('throws for unknown string', () => {
      expect(() => new iinfo('float32' as any)).toThrow('Unknown integer dtype');
    });
  });

  describe('toString', () => {
    it('returns readable representation', () => {
      const i32 = new iinfo(DType.Int32);
      const str = i32.toString();
      expect(str).toContain('iinfo');
      expect(str).toContain('bits=32');
      expect(str).toContain('min=-2147483648');
      expect(str).toContain('max=2147483647');
    });
  });
});
