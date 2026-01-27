/**
 * Tests for Phase 26: Miscellaneous Ufuncs
 */

import { describe, it, expect, beforeAll, afterEach } from 'vitest';
import {
  loadWasmModule,
  NDArray,
  DType,
  frexp,
  ldexp,
  nextafter,
  spacing,
  modf,
  gcd,
  lcm,
  sinc,
  heaviside,
  divmod,
  bitwise_count,
} from '../../dist/numjs.mjs';

const resources: NDArray[] = [];
function track<T extends NDArray>(arr: T): T {
  resources.push(arr);
  return arr;
}

describe('Phase 26: Miscellaneous Ufuncs', () => {
  beforeAll(async () => {
    await loadWasmModule();
  });

  afterEach(() => {
    for (const arr of resources) {
      if (!arr.isDisposed) {
        arr.dispose();
      }
    }
    resources.length = 0;
  });

  // ========================================================================
  // frexp - Floating Point Decomposition
  // ========================================================================
  describe('frexp', () => {
    it('decomposes floats into mantissa and exponent', async () => {
      const arr = track(await NDArray.fromArray([1.0, 2.0, 4.0, 8.0]));
      const [mantissa, exponent] = frexp(arr);
      track(mantissa);
      track(exponent);

      expect(mantissa.toArray()).toEqual([0.5, 0.5, 0.5, 0.5]);
      expect(exponent.toArray()).toEqual([1, 2, 3, 4]);
    });

    it('handles zero correctly', async () => {
      const arr = track(await NDArray.fromArray([0.0]));
      const [mantissa, exponent] = frexp(arr);
      track(mantissa);
      track(exponent);

      expect(mantissa.toArray()[0]).toBe(0);
      expect(exponent.toArray()[0]).toBe(0);
    });

    it('handles negative values', async () => {
      const arr = track(await NDArray.fromArray([-1.0, -2.0]));
      const [mantissa, exponent] = frexp(arr);
      track(mantissa);
      track(exponent);

      expect(mantissa.toArray()).toEqual([-0.5, -0.5]);
      expect(exponent.toArray()).toEqual([1, 2]);
    });

    it('roundtrips with ldexp', async () => {
      const original = track(await NDArray.fromArray([1.5, 3.14159, -2.5, 100.0]));
      const [mantissa, exponent] = frexp(original);
      track(mantissa);
      track(exponent);

      const reconstructed = track(ldexp(mantissa, exponent));
      const origVals = original.toArray() as number[];
      const reconVals = reconstructed.toArray() as number[];

      for (let i = 0; i < origVals.length; i++) {
        expect(reconVals[i]).toBeCloseTo(origVals[i], 10);
      }
    });

    it('handles infinity', async () => {
      const arr = track(await NDArray.fromArray([Infinity]));
      const [mantissa, exponent] = frexp(arr);
      track(mantissa);
      track(exponent);

      // mantissa should be Infinity for Infinity input
      expect(mantissa.toArray()[0]).toBe(Infinity);
      // exponent value for infinity is implementation-defined, just check it's a number
      expect(typeof exponent.toArray()[0]).toBe('number');
    });

    it('handles NaN', async () => {
      const arr = track(await NDArray.fromArray([NaN]));
      const [mantissa, exponent] = frexp(arr);
      track(mantissa);
      track(exponent);

      expect(Number.isNaN(mantissa.toArray()[0])).toBe(true);
    });
  });

  // ========================================================================
  // ldexp - Construct from mantissa and exponent
  // ========================================================================
  describe('ldexp', () => {
    it('computes x1 * 2^x2', async () => {
      const mantissa = track(await NDArray.fromArray([0.5, 0.5, 0.5]));
      const exponent = track(
        await NDArray.fromArray([1, 2, 3], undefined, { dtype: DType.Int32 })
      );
      const result = track(ldexp(mantissa, exponent));

      expect(result.toArray()).toEqual([1.0, 2.0, 4.0]);
    });

    it('handles negative exponents', async () => {
      const mantissa = track(await NDArray.fromArray([1.0, 1.0, 1.0]));
      const exponent = track(
        await NDArray.fromArray([-1, -2, -3], undefined, { dtype: DType.Int32 })
      );
      const result = track(ldexp(mantissa, exponent));

      expect(result.toArray()).toEqual([0.5, 0.25, 0.125]);
    });

    it('handles zero mantissa', async () => {
      const mantissa = track(await NDArray.fromArray([0.0]));
      const exponent = track(
        await NDArray.fromArray([10], undefined, { dtype: DType.Int32 })
      );
      const result = track(ldexp(mantissa, exponent));

      expect(result.toArray()[0]).toBe(0);
    });
  });

  // ========================================================================
  // nextafter - Next representable float
  // ========================================================================
  describe('nextafter', () => {
    it('returns next float toward direction (positive)', async () => {
      const x = track(await NDArray.fromArray([1.0]));
      const toward = track(await NDArray.fromArray([2.0]));
      const result = track(nextafter(x, toward));

      expect(result.toArray()[0]).toBeGreaterThan(1.0);
    });

    it('returns next float toward direction (negative)', async () => {
      const x = track(await NDArray.fromArray([1.0]));
      const toward = track(await NDArray.fromArray([0.0]));
      const result = track(nextafter(x, toward));

      expect(result.toArray()[0]).toBeLessThan(1.0);
    });

    it('returns same value when equal', async () => {
      const x = track(await NDArray.fromArray([1.0]));
      const same = track(await NDArray.fromArray([1.0]));
      const result = track(nextafter(x, same));

      expect(result.toArray()[0]).toBe(1.0);
    });

    it('handles zero', async () => {
      const x = track(await NDArray.fromArray([0.0]));
      const toward = track(await NDArray.fromArray([1.0]));
      const result = track(nextafter(x, toward));

      expect(result.toArray()[0]).toBeGreaterThan(0);
      expect(result.toArray()[0]).toBeLessThan(1e-300);
    });

    it('returns NaN if either input is NaN', async () => {
      const x = track(await NDArray.fromArray([NaN]));
      const toward = track(await NDArray.fromArray([1.0]));
      const result = track(nextafter(x, toward));

      expect(Number.isNaN(result.toArray()[0])).toBe(true);
    });
  });

  // ========================================================================
  // spacing - ULP distance
  // ========================================================================
  describe('spacing', () => {
    it('returns machine epsilon at 1.0', async () => {
      const arr = track(await NDArray.fromArray([1.0]));
      const result = track(spacing(arr));

      // Machine epsilon for float64 is approximately 2.22e-16
      expect(result.toArray()[0]).toBeCloseTo(2.220446049250313e-16, 25);
    });

    it('returns larger spacing for larger numbers', async () => {
      const arr = track(await NDArray.fromArray([1.0, 1e10]));
      const result = track(spacing(arr));
      const vals = result.toArray() as number[];

      expect(vals[1]).toBeGreaterThan(vals[0]);
    });

    it('handles zero', async () => {
      const arr = track(await NDArray.fromArray([0.0]));
      const result = track(spacing(arr));

      expect(result.toArray()[0]).toBeGreaterThan(0);
      expect(result.toArray()[0]).toBeLessThan(1e-300);
    });

    it('returns NaN for infinity', async () => {
      const arr = track(await NDArray.fromArray([Infinity]));
      const result = track(spacing(arr));

      expect(Number.isNaN(result.toArray()[0])).toBe(true);
    });
  });

  // ========================================================================
  // modf - Split into fractional and integral
  // ========================================================================
  describe('modf', () => {
    it('splits positive floats', async () => {
      const arr = track(await NDArray.fromArray([3.5]));
      const [frac, intg] = modf(arr);
      track(frac);
      track(intg);

      expect(frac.toArray()[0]).toBeCloseTo(0.5, 10);
      expect(intg.toArray()[0]).toBeCloseTo(3.0, 10);
    });

    it('splits negative floats correctly', async () => {
      const arr = track(await NDArray.fromArray([-2.7]));
      const [frac, intg] = modf(arr);
      track(frac);
      track(intg);

      expect(frac.toArray()[0]).toBeCloseTo(-0.7, 10);
      expect(intg.toArray()[0]).toBeCloseTo(-2.0, 10);
    });

    it('handles integer values', async () => {
      const arr = track(await NDArray.fromArray([5.0, -3.0]));
      const [frac, intg] = modf(arr);
      track(frac);
      track(intg);

      expect(frac.toArray()).toEqual([0, -0]);
      expect(intg.toArray()).toEqual([5, -3]);
    });

    it('handles zero', async () => {
      const arr = track(await NDArray.fromArray([0.0]));
      const [frac, intg] = modf(arr);
      track(frac);
      track(intg);

      expect(frac.toArray()[0]).toBe(0);
      expect(intg.toArray()[0]).toBe(0);
    });

    it('handles infinity', async () => {
      const arr = track(await NDArray.fromArray([Infinity]));
      const [frac, intg] = modf(arr);
      track(frac);
      track(intg);

      expect(frac.toArray()[0]).toBe(0);
      expect(intg.toArray()[0]).toBe(Infinity);
    });

    it('handles NaN', async () => {
      const arr = track(await NDArray.fromArray([NaN]));
      const [frac, intg] = modf(arr);
      track(frac);
      track(intg);

      expect(Number.isNaN(frac.toArray()[0])).toBe(true);
      expect(Number.isNaN(intg.toArray()[0])).toBe(true);
    });
  });

  // ========================================================================
  // gcd - Greatest Common Divisor
  // ========================================================================
  describe('gcd', () => {
    it('computes gcd correctly', async () => {
      const a = track(
        await NDArray.fromArray([12, 15, 20], undefined, { dtype: DType.Int32 })
      );
      const b = track(
        await NDArray.fromArray([8, 10, 15], undefined, { dtype: DType.Int32 })
      );
      const result = track(gcd(a, b));

      expect(result.toArray()).toEqual([4, 5, 5]);
    });

    it('handles zero correctly (gcd(0, n) = n)', async () => {
      const a = track(
        await NDArray.fromArray([0, 5], undefined, { dtype: DType.Int32 })
      );
      const b = track(
        await NDArray.fromArray([5, 0], undefined, { dtype: DType.Int32 })
      );
      const result = track(gcd(a, b));

      expect(result.toArray()).toEqual([5, 5]);
    });

    it('handles gcd(0, 0) = 0', async () => {
      const a = track(await NDArray.fromArray([0], undefined, { dtype: DType.Int32 }));
      const b = track(await NDArray.fromArray([0], undefined, { dtype: DType.Int32 }));
      const result = track(gcd(a, b));

      expect(result.toArray()[0]).toBe(0);
    });

    it('handles negative values (uses absolute)', async () => {
      const a = track(
        await NDArray.fromArray([-12, 12], undefined, { dtype: DType.Int32 })
      );
      const b = track(
        await NDArray.fromArray([8, -8], undefined, { dtype: DType.Int32 })
      );
      const result = track(gcd(a, b));

      expect(result.toArray()).toEqual([4, 4]);
    });

    it('throws on float input', async () => {
      const a = track(await NDArray.fromArray([1.5]));
      const b = track(await NDArray.fromArray([2.5]));

      expect(() => gcd(a, b)).toThrow(TypeError);
    });
  });

  // ========================================================================
  // lcm - Lowest Common Multiple
  // ========================================================================
  describe('lcm', () => {
    it('computes lcm correctly', async () => {
      const a = track(
        await NDArray.fromArray([12, 4, 6], undefined, { dtype: DType.Int32 })
      );
      const b = track(
        await NDArray.fromArray([8, 8, 9], undefined, { dtype: DType.Int32 })
      );
      const result = track(lcm(a, b));

      expect(result.toArray()).toEqual([24, 8, 18]);
    });

    it('returns zero when either input is zero', async () => {
      const a = track(
        await NDArray.fromArray([0, 5], undefined, { dtype: DType.Int32 })
      );
      const b = track(
        await NDArray.fromArray([5, 0], undefined, { dtype: DType.Int32 })
      );
      const result = track(lcm(a, b));

      expect(result.toArray()).toEqual([0, 0]);
    });

    it('handles negative values (uses absolute)', async () => {
      const a = track(
        await NDArray.fromArray([-4, 4], undefined, { dtype: DType.Int32 })
      );
      const b = track(
        await NDArray.fromArray([6, -6], undefined, { dtype: DType.Int32 })
      );
      const result = track(lcm(a, b));

      expect(result.toArray()).toEqual([12, 12]);
    });

    it('throws on float input', async () => {
      const a = track(await NDArray.fromArray([1.5]));
      const b = track(await NDArray.fromArray([2.5]));

      expect(() => lcm(a, b)).toThrow(TypeError);
    });
  });

  // ========================================================================
  // sinc - Normalized sinc function
  // ========================================================================
  describe('sinc', () => {
    it('returns 1 at zero', async () => {
      const arr = track(await NDArray.fromArray([0]));
      const result = track(sinc(arr));

      expect(result.toArray()[0]).toBe(1.0);
    });

    it('returns approximately 0 at integer values', async () => {
      const arr = track(await NDArray.fromArray([1, -1, 2, -2]));
      const result = track(sinc(arr));
      const vals = result.toArray() as number[];

      for (const val of vals) {
        expect(Math.abs(val)).toBeLessThan(1e-10);
      }
    });

    it('computes sin(pi*x)/(pi*x) correctly', async () => {
      const arr = track(await NDArray.fromArray([0.5]));
      const result = track(sinc(arr));

      // sinc(0.5) = sin(pi/2) / (pi/2) = 1 / (pi/2) = 2/pi
      expect(result.toArray()[0]).toBeCloseTo(2 / Math.PI, 10);
    });

    it('is symmetric', async () => {
      const arr = track(await NDArray.fromArray([0.3, -0.3]));
      const result = track(sinc(arr));
      const vals = result.toArray() as number[];

      expect(vals[0]).toBeCloseTo(vals[1], 10);
    });
  });

  // ========================================================================
  // heaviside - Step function
  // ========================================================================
  describe('heaviside', () => {
    it('returns 0 for negative, h0 for zero, 1 for positive', async () => {
      const x = track(await NDArray.fromArray([-1, 0, 1]));
      const h0 = track(await NDArray.fromArray([0.5, 0.5, 0.5]));
      const result = track(heaviside(x, h0));

      expect(result.toArray()).toEqual([0, 0.5, 1]);
    });

    it('handles scalar h0 via broadcasting', async () => {
      const x = track(await NDArray.fromArray([-2, -1, 0, 1, 2]));
      const h0 = track(await NDArray.fromArray([1.0]));
      const result = track(heaviside(x, h0));

      expect(result.toArray()).toEqual([0, 0, 1, 1, 1]);
    });

    it('handles NaN correctly', async () => {
      const x = track(await NDArray.fromArray([NaN]));
      const h0 = track(await NDArray.fromArray([0.5]));
      const result = track(heaviside(x, h0));

      expect(Number.isNaN(result.toArray()[0])).toBe(true);
    });

    it('works with different h0 values', async () => {
      const x = track(await NDArray.fromArray([0, 0, 0]));
      const h0 = track(await NDArray.fromArray([0, 0.5, 1]));
      const result = track(heaviside(x, h0));

      expect(result.toArray()).toEqual([0, 0.5, 1]);
    });
  });

  // ========================================================================
  // divmod - Quotient and Remainder
  // ========================================================================
  describe('divmod', () => {
    it('returns quotient and remainder for integers', async () => {
      const a = track(
        await NDArray.fromArray([10, 11, 12], undefined, { dtype: DType.Int32 })
      );
      const b = track(
        await NDArray.fromArray([3, 3, 3], undefined, { dtype: DType.Int32 })
      );
      const [quot, rem] = divmod(a, b);
      track(quot);
      track(rem);

      expect(quot.toArray()).toEqual([3, 3, 4]);
      expect(rem.toArray()).toEqual([1, 2, 0]);
    });

    it('returns quotient and remainder for floats', async () => {
      const a = track(await NDArray.fromArray([7.5]));
      const b = track(await NDArray.fromArray([2.5]));
      const [quot, rem] = divmod(a, b);
      track(quot);
      track(rem);

      expect(quot.toArray()[0]).toBe(3.0);
      expect(rem.toArray()[0]).toBeCloseTo(0.0, 10);
    });

    it('handles negative dividends', async () => {
      const a = track(
        await NDArray.fromArray([-7], undefined, { dtype: DType.Int32 })
      );
      const b = track(await NDArray.fromArray([3], undefined, { dtype: DType.Int32 }));
      const [quot, rem] = divmod(a, b);
      track(quot);
      track(rem);

      // Python-style floor division: -7 // 3 = -3, -7 % 3 = 2
      expect(quot.toArray()[0]).toBe(-3);
      expect(rem.toArray()[0]).toBe(2);
    });
  });

  // ========================================================================
  // bitwise_count - Population count
  // ========================================================================
  describe('bitwise_count', () => {
    it('counts set bits', async () => {
      const arr = track(
        await NDArray.fromArray([0, 1, 2, 3, 4, 5, 6, 7], undefined, {
          dtype: DType.Int32,
        })
      );
      const result = track(bitwise_count(arr));

      expect(result.toArray()).toEqual([0, 1, 1, 2, 1, 2, 2, 3]);
    });

    it('handles larger numbers', async () => {
      const arr = track(
        await NDArray.fromArray([255, 256], undefined, { dtype: DType.Int32 })
      );
      const result = track(bitwise_count(arr));

      // 255 = 0xFF = 8 bits set
      // 256 = 0x100 = 1 bit set
      expect(result.toArray()).toEqual([8, 1]);
    });

    it('throws on float input', async () => {
      const arr = track(await NDArray.fromArray([1.5]));

      expect(() => bitwise_count(arr)).toThrow(TypeError);
    });

    it('output dtype is uint8', async () => {
      const arr = track(
        await NDArray.fromArray([7], undefined, { dtype: DType.Int32 })
      );
      const result = track(bitwise_count(arr));

      expect(result.dtype).toBe(DType.Uint8);
    });
  });
});
