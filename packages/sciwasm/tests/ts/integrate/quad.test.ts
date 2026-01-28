/**
 * Tests for sciwasm.integrate.quad()
 *
 * Ported from scipy/integrate/tests/test_quadpack.py (TestQuad class).
 */

import { describe, it, expect } from 'vitest';
import { integrate } from '../../../src/ts/index.js';

const { quad } = integrate;

/**
 * Assert that a quad result matches an expected value within the reported error.
 * Matches scipy's assert_quad helper.
 */
function assertQuad(
  result: [number, number],
  expected: number,
  errorTol = 1.5e-8,
) {
  const [value, err] = result;
  expect(Math.abs(value - expected)).toBeLessThanOrEqual(Math.max(err, 1e-15));
  if (errorTol != null) {
    expect(err).toBeLessThan(errorTol);
  }
}

describe('integrate.quad', () => {
  // --- Basic polynomial integration ---

  it('should integrate x^2 from 0 to 1', async () => {
    const result = await quad((x: number) => x * x, 0, 1);
    assertQuad(result as [number, number], 1 / 3);
  });

  it('should integrate x^2 from 0 to 4', async () => {
    const result = await quad((x: number) => x * x, 0, 4);
    assertQuad(result as [number, number], 64 / 3);
  });

  it('should integrate constant function', async () => {
    const result = await quad(() => 5, 0, 3);
    assertQuad(result as [number, number], 15);
  });

  // --- scipy TestQuad.test_typical ---
  // Bessel function integrand with extra arguments
  it('should integrate Bessel-type integrand with extra args (scipy test_typical)', async () => {
    const myfunc = (x: number, n: number, z: number) =>
      Math.cos(n * x - z * Math.sin(x)) / Math.PI;

    const result = await quad(myfunc, 0, Math.PI, { args: [2, 1.8] });
    assertQuad(result as [number, number], 0.30614353532540296487);
  });

  // --- scipy TestQuad.test_indefinite ---
  // Euler's constant via infinite integral
  it('should integrate to Euler constant with infinite limit (scipy test_indefinite)', async () => {
    const myfunc = (x: number) => -Math.exp(-x) * Math.log(x);

    const result = await quad(myfunc, 0, Infinity);
    assertQuad(result as [number, number], 0.577215664901532860606512);
  });

  // --- Exponential decay over [0, +inf) ---
  it('should integrate exp(-x) from 0 to +inf', async () => {
    const result = await quad((x: number) => Math.exp(-x), 0, Infinity);
    assertQuad(result as [number, number], 1.0);
  });

  // --- Gaussian integral ---
  it('should integrate Gaussian over (-inf, +inf)', async () => {
    const result = await quad(
      (x: number) => Math.exp(-x * x),
      -Infinity,
      Infinity,
    );
    assertQuad(result as [number, number], Math.sqrt(Math.PI));
  });

  // --- scipy TestQuad.test_b_less_than_a ---
  it('should negate result when b < a (scipy test_b_less_than_a)', async () => {
    const f = (x: number, p: number, q: number) => p * Math.exp(-q * x);

    const [val1, err1] = (await quad(f, 0, Infinity, { args: [2, 3] })) as [number, number];
    const [val2, err2] = (await quad(f, Infinity, 0, { args: [2, 3] })) as [number, number];

    expect(Math.abs(val1 + val2)).toBeLessThanOrEqual(Math.max(err1, err2));
  });

  // --- scipy TestQuad.test_b_less_than_a_2 ---
  it('should negate result for double-infinite flipped (scipy test_b_less_than_a_2)', async () => {
    const f = (x: number, s: number) =>
      Math.exp((-x * x) / 2 / s) / Math.sqrt(2 * s);

    const [val1, err1] = (await quad(f, -Infinity, Infinity, { args: [2] })) as [number, number];
    const [val2, err2] = (await quad(f, Infinity, -Infinity, { args: [2] })) as [number, number];

    expect(Math.abs(val1 + val2)).toBeLessThanOrEqual(Math.max(err1, err2));
  });

  // --- scipy TestQuad.test_b_equals_a ---
  it('should return (0, 0) when a === b (scipy test_b_equals_a)', async () => {
    const f = (x: number) => 1 / x;

    const result = await quad(f, 0, 0);
    expect(result[0]).toBe(0);
    expect(result[1]).toBe(0);
  });

  it('should return infodict when a === b with fullOutput', async () => {
    const f = (x: number) => 1 / x;
    const limit = 50;

    const result = await quad(f, 0, 0, { fullOutput: true, limit });
    expect(result[0]).toBe(0);
    expect(result[1]).toBe(0);

    const infodict = result[2] as {
      neval: number;
      last: number;
      alist: Float64Array;
      blist: Float64Array;
      rlist: Float64Array;
      elist: Float64Array;
      iord: Int32Array;
    };
    expect(infodict.neval).toBe(0);
    expect(infodict.last).toBe(0);
    expect(infodict.alist.length).toBe(limit);
    expect(infodict.blist.length).toBe(limit);
    expect(infodict.rlist.length).toBe(limit);
    expect(infodict.elist.length).toBe(limit);
    expect(infodict.iord.length).toBe(limit);
  });

  // --- Error tolerance parameters ---
  it('should respect epsabs and epsrel', async () => {
    const f = (x: number) => Math.sin(x);

    // Default tolerances
    const [val1, err1] = (await quad(f, 0, Math.PI)) as [number, number];
    assertQuad([val1, err1], 2.0);

    // Tighter tolerance
    const [val2, err2] = (await quad(f, 0, Math.PI, {
      epsabs: 1e-12,
      epsrel: 1e-12,
    })) as [number, number];
    assertQuad([val2, err2], 2.0, 1e-11);
  });

  // --- fullOutput mode ---
  it('should return infodict with fullOutput', async () => {
    const result = await quad((x: number) => x * x, 0, 1, { fullOutput: true });
    expect(result.length).toBeGreaterThanOrEqual(3);

    const infodict = result[2] as {
      neval: number;
      last: number;
    };
    expect(infodict.neval).toBeGreaterThan(0);
    expect(infodict.last).toBeGreaterThan(0);
  });

  // --- Semi-infinite intervals ---
  it('should handle (-inf, b) intervals', async () => {
    const result = await quad(
      (x: number) => Math.exp(x),
      -Infinity,
      0,
    );
    assertQuad(result as [number, number], 1.0);
  });

  // --- Single arg (not array) ---
  it('should accept args as a single number', async () => {
    const f = (x: number, a: number) => Math.exp(a * (x - 1));
    const result = await quad(f, 0, 1, { args: -5 });
    // integral of exp(-5(x-1)) from 0 to 1 = (exp(5) - 1) / 5
    const expected = (Math.exp(5) - 1) / 5;
    assertQuad(result as [number, number], expected);
  });
});
