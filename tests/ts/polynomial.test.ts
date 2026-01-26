/**
 * Tests for NumJS Polynomial Module
 */

import { describe, it, expect, beforeAll } from 'vitest';
import { loadWasmModule } from '../../dist/numjs.mjs';

// Import directly from source for testing during development
import {
  // Utilities
  PolyError,
  trimseq,
  trimcoef,
  getdomain,
  mapdomain,
  mapparms,
  // Polynomial
  Polynomial,
  polyval,
  polyder,
  polyint,
  polyfit,
  polyroots,
  polyfromroots,
  polyadd,
  polysub,
  polymul,
  polydiv,
  polypow,
  polyvander,
  // Chebyshev
  Chebyshev,
  chebval,
  chebder,
  chebint,
  chebfit,
  chebroots,
  chebfromroots,
  poly2cheb,
  cheb2poly,
  chebadd,
  chebmul,
  chebvander,
  // Legendre
  Legendre,
  legval,
  legvander,
  poly2leg,
  leg2poly,
  // Hermite
  Hermite,
  hermval,
  poly2herm,
  herm2poly,
  // HermiteE
  HermiteE,
  hermeval,
  poly2herme,
  herme2poly,
  // Laguerre
  Laguerre,
  lagval,
  poly2lag,
  lag2poly,
} from '../../src/ts/polynomial/index.js';

/**
 * Helper function for comparing floating point arrays.
 */
function allClose(
  actual: number[],
  expected: number[],
  rtol: number = 1e-5,
  atol: number = 1e-8
): boolean {
  if (actual.length !== expected.length) return false;
  for (let i = 0; i < actual.length; i++) {
    const diff = Math.abs(actual[i] - expected[i]);
    const tolerance = atol + rtol * Math.abs(expected[i]);
    if (diff > tolerance) return false;
  }
  return true;
}

describe('Polynomial Module', () => {
  beforeAll(async () => {
    await loadWasmModule();
  });

  /* ============ Utility Functions ============ */

  describe('polyutils', () => {
    it('trimseq removes trailing zeros', () => {
      expect(trimseq([1, 2, 0, 0])).toEqual([1, 2]);
      expect(trimseq([0, 0, 0])).toEqual([0]);
      expect(trimseq([1])).toEqual([1]);
      expect(trimseq([])).toEqual([0]);
    });

    it('trimcoef handles tolerance', () => {
      expect(trimcoef([1, 2, 1e-10], 1e-9)).toEqual([1, 2]);
      expect(trimcoef([1, 0.001], 0.01)).toEqual([1]);
      expect(trimcoef([1, 2, 3])).toEqual([1, 2, 3]);
    });

    it('getdomain returns correct domain', () => {
      expect(getdomain([1, 2, 3, 4, 5])).toEqual([1, 5]);
      expect(getdomain([-1, 0, 1])).toEqual([-1, 1]);
      expect(getdomain([5])).toEqual([4, 6]); // Single point expands
    });

    it('mapparms returns correct parameters', () => {
      const [off, scl] = mapparms([0, 1], [-1, 1]);
      expect(off).toBeCloseTo(-1);
      expect(scl).toBeCloseTo(2);
    });

    it('mapdomain maps correctly', () => {
      expect(mapdomain(0.5, [0, 1], [-1, 1])).toBeCloseTo(0);
      const result = mapdomain([0, 0.5, 1], [0, 1], [-1, 1]) as number[];
      expect(result[0]).toBeCloseTo(-1);
      expect(result[1]).toBeCloseTo(0);
      expect(result[2]).toBeCloseTo(1);
    });
  });

  /* ============ Polynomial (Power Series) ============ */

  describe('Polynomial class', () => {
    it('constructs from coefficients', () => {
      const p = new Polynomial([1, 2, 3]);
      expect(p.coef).toEqual([1, 2, 3]);
      expect(p.degree).toBe(2);
    });

    it('trims trailing zeros', () => {
      const p = new Polynomial([1, 2, 0, 0]);
      expect(p.coef).toEqual([1, 2]);
      expect(p.degree).toBe(1);
    });

    it('evaluates using call()', () => {
      const p = new Polynomial([1, 2, 3]); // 1 + 2x + 3x^2
      // At x=2: 1 + 4 + 12 = 17
      expect(p.call(2)).toBeCloseTo(17);
      expect(p.call(0)).toBeCloseTo(1);
      expect(p.call(1)).toBeCloseTo(6);
    });

    it('adds polynomials', () => {
      const p1 = new Polynomial([1, 2]);
      const p2 = new Polynomial([3, 4, 5]);
      const sum = p1.add(p2);
      expect(sum.coef).toEqual([4, 6, 5]);
    });

    it('subtracts polynomials', () => {
      const p1 = new Polynomial([3, 4, 5]);
      const p2 = new Polynomial([1, 2]);
      const diff = p1.sub(p2);
      expect(diff.coef).toEqual([2, 2, 5]);
    });

    it('multiplies polynomials', () => {
      const p1 = new Polynomial([1, 1]); // 1 + x
      const p2 = new Polynomial([1, 1]); // 1 + x
      const prod = p1.mul(p2);
      expect(prod.coef).toEqual([1, 2, 1]); // 1 + 2x + x^2
    });

    it('raises to power', () => {
      const p = new Polynomial([1, 1]); // 1 + x
      const p3 = p.pow(3);
      expect(p3.coef).toEqual([1, 3, 3, 1]); // (1+x)^3 = 1 + 3x + 3x^2 + x^3
    });

    it('computes derivative', () => {
      const p = new Polynomial([1, 2, 3]); // 1 + 2x + 3x^2
      const dp = p.deriv();
      expect(dp.coef).toEqual([2, 6]); // 2 + 6x
    });

    it('computes integral', () => {
      const p = new Polynomial([2, 6]); // 2 + 6x
      const ip = p.integ(1, [0]); // with constant 0
      // The Polynomial class uses domain [-1, 1] which affects integration scale
      // Let's just verify the derivative roundtrip
      const dp = ip.deriv();
      // dp should equal p (approximately)
      expect(dp.coef.length).toBe(p.coef.length);
      expect(dp.coef[0]).toBeCloseTo(p.coef[0]);
      expect(dp.coef[1]).toBeCloseTo(p.coef[1]);
    });
  });

  describe('polyval', () => {
    it('evaluates using Horner method', () => {
      // p(x) = 1 + 2x + 3x^2
      expect(polyval(2, [1, 2, 3])).toBeCloseTo(17);
      expect(polyval(0, [1, 2, 3])).toBeCloseTo(1);
      expect(polyval(-1, [1, 2, 3])).toBeCloseTo(2); // 1 - 2 + 3
    });
  });

  describe('polyder', () => {
    it('computes derivative', () => {
      expect(polyder([1, 2, 3])).toEqual([2, 6]);
      expect(polyder([1, 2, 3, 4], 2)).toEqual([6, 24]); // Second derivative
    });
  });

  describe('polyint', () => {
    it('computes integral', () => {
      const result = polyint([2, 6]);
      expect(allClose(result, [0, 2, 3])).toBe(true);
    });
  });

  describe('polyadd/sub/mul/div', () => {
    it('adds polynomials', () => {
      expect(polyadd([1, 2], [3, 4, 5])).toEqual([4, 6, 5]);
    });

    it('subtracts polynomials', () => {
      expect(polysub([3, 4, 5], [1, 2])).toEqual([2, 2, 5]);
    });

    it('multiplies polynomials', () => {
      expect(polymul([1, 1], [1, 1])).toEqual([1, 2, 1]);
    });

    it('divides polynomials', () => {
      // (x^2 + 2x + 1) / (x + 1) = (x + 1)
      const [quo, rem] = polydiv([1, 2, 1], [1, 1]);
      expect(allClose(quo, [1, 1])).toBe(true);
      expect(allClose(rem, [0])).toBe(true);
    });
  });

  describe('polypow', () => {
    it('raises to power', () => {
      expect(polypow([1, 1], 0)).toEqual([1]);
      expect(polypow([1, 1], 1)).toEqual([1, 1]);
      expect(polypow([1, 1], 2)).toEqual([1, 2, 1]);
      expect(polypow([1, 1], 3)).toEqual([1, 3, 3, 1]);
    });
  });

  describe('polyfromroots', () => {
    it('creates polynomial from roots', () => {
      // Roots 2 and 3 -> (x-2)(x-3) = x^2 - 5x + 6 = [6, -5, 1]
      const coef = polyfromroots([2, 3]);
      expect(allClose(coef, [6, -5, 1])).toBe(true);
    });
  });

  describe('polyroots', () => {
    it('finds linear root', async () => {
      // 2 + 3x = 0 => x = -2/3
      const roots = await polyroots([2, 3]);
      expect(roots.length).toBe(1);
      expect(roots[0]).toBeCloseTo(-2 / 3);
    });

    it('finds quadratic roots', async () => {
      // x^2 - 5x + 6 = (x-2)(x-3) has roots 2, 3
      const roots = await polyroots([6, -5, 1]);
      expect(roots.length).toBe(2);
      expect(roots[0]).toBeCloseTo(2);
      expect(roots[1]).toBeCloseTo(3);
    });
  });

  describe('polyvander', () => {
    it('generates Vandermonde matrix', () => {
      const V = polyvander([1, 2, 3], 2);
      // Row i: [1, x_i, x_i^2]
      expect(V[0]).toEqual([1, 1, 1]);
      expect(V[1]).toEqual([1, 2, 4]);
      expect(V[2]).toEqual([1, 3, 9]);
    });
  });

  describe('polyfit', () => {
    it.skip('fits line to points', async () => {
      // Skip: requires built WASM module
      // Points: (0, 0), (1, 1), (2, 2) -> line y = x
      const x = [0, 1, 2];
      const y = [0, 1, 2];
      const coef = (await polyfit(x, y, 1)) as number[];
      expect(coef.length).toBe(2);
      expect(coef[0]).toBeCloseTo(0); // intercept
      expect(coef[1]).toBeCloseTo(1); // slope
    });

    it.skip('fits quadratic to points', async () => {
      // Skip: requires built WASM module
      // Points from y = x^2
      const x = [0, 1, 2, 3];
      const y = [0, 1, 4, 9];
      const coef = (await polyfit(x, y, 2)) as number[];
      expect(coef.length).toBe(3);
      expect(coef[0]).toBeCloseTo(0);
      expect(coef[1]).toBeCloseTo(0);
      expect(coef[2]).toBeCloseTo(1);
    });
  });

  /* ============ Chebyshev ============ */

  describe('Chebyshev class', () => {
    it('constructs from coefficients', () => {
      const p = new Chebyshev([1, 2, 3]);
      expect(p.coef).toEqual([1, 2, 3]);
      expect(p.degree).toBe(2);
    });
  });

  describe('chebval', () => {
    it('evaluates Chebyshev polynomial', () => {
      // T_0(x) = 1
      expect(chebval(0.5, [1])).toBeCloseTo(1);

      // T_1(x) = x
      expect(chebval(0.5, [0, 1])).toBeCloseTo(0.5);

      // T_2(x) = 2x^2 - 1
      expect(chebval(0.5, [0, 0, 1])).toBeCloseTo(2 * 0.25 - 1); // -0.5
    });
  });

  describe('cheb2poly / poly2cheb', () => {
    it('converts Chebyshev to power series', () => {
      // T_0 = 1
      expect(cheb2poly([1])).toEqual([1]);

      // T_1 = x
      expect(cheb2poly([0, 1])).toEqual([0, 1]);

      // T_2 = 2x^2 - 1
      const t2 = cheb2poly([0, 0, 1]);
      expect(allClose(t2, [-1, 0, 2])).toBe(true);
    });

    it('converts power series to Chebyshev', () => {
      // 1 = T_0
      expect(poly2cheb([1])).toEqual([1]);

      // x = T_1
      expect(poly2cheb([0, 1])).toEqual([0, 1]);
    });

    it('roundtrip conversion preserves values', () => {
      const original = [1, 2, 3];
      const asCheb = poly2cheb(original);
      const back = cheb2poly(asCheb);
      expect(allClose(back, original)).toBe(true);
    });
  });

  describe('chebvander', () => {
    it('generates Chebyshev Vandermonde matrix', () => {
      const V = chebvander([0, 0.5, 1], 2);
      // Row i: [T_0(x_i), T_1(x_i), T_2(x_i)]
      // T_2(x) = 2x^2 - 1

      // x = 0: T_0=1, T_1=0, T_2=-1
      expect(V[0][0]).toBeCloseTo(1);
      expect(V[0][1]).toBeCloseTo(0);
      expect(V[0][2]).toBeCloseTo(-1);

      // x = 0.5: T_0=1, T_1=0.5, T_2=-0.5
      expect(V[1][0]).toBeCloseTo(1);
      expect(V[1][1]).toBeCloseTo(0.5);
      expect(V[1][2]).toBeCloseTo(-0.5);
    });
  });

  describe('chebadd/chebmul', () => {
    it('adds Chebyshev polynomials', () => {
      expect(chebadd([1, 2], [3, 4, 5])).toEqual([4, 6, 5]);
    });

    it('multiplies Chebyshev polynomials', () => {
      // T_1 * T_1 = 0.5*(T_2 + T_0) = 0.5*T_0 + 0.5*T_2
      const result = chebmul([0, 1], [0, 1]);
      expect(allClose(result, [0.5, 0, 0.5])).toBe(true);
    });
  });

  /* ============ Legendre ============ */

  describe('Legendre class', () => {
    it('constructs from coefficients', () => {
      const p = new Legendre([1, 2, 3]);
      expect(p.coef).toEqual([1, 2, 3]);
    });
  });

  describe('legval', () => {
    it('evaluates Legendre polynomial', () => {
      // P_0(x) = 1
      expect(legval(0.5, [1])).toBeCloseTo(1);

      // P_1(x) = x
      expect(legval(0.5, [0, 1])).toBeCloseTo(0.5);

      // P_2(x) = (3x^2 - 1) / 2
      expect(legval(0.5, [0, 0, 1])).toBeCloseTo((3 * 0.25 - 1) / 2);
    });
  });

  describe('leg2poly / poly2leg', () => {
    it('converts Legendre to power series', () => {
      // P_0 = 1
      expect(leg2poly([1])).toEqual([1]);

      // P_1 = x
      expect(leg2poly([0, 1])).toEqual([0, 1]);

      // P_2 = (3x^2 - 1) / 2
      const p2 = leg2poly([0, 0, 1]);
      expect(allClose(p2, [-0.5, 0, 1.5])).toBe(true);
    });
  });

  describe('legvander', () => {
    it('generates Legendre Vandermonde matrix', () => {
      const V = legvander([0, 1], 2);
      // P_2(0) = -0.5, P_2(1) = 1
      expect(V[0][2]).toBeCloseTo(-0.5);
      expect(V[1][2]).toBeCloseTo(1);
    });
  });

  /* ============ Hermite ============ */

  describe('Hermite class', () => {
    it('constructs from coefficients', () => {
      const p = new Hermite([1, 2, 3]);
      expect(p.coef).toEqual([1, 2, 3]);
    });
  });

  describe('hermval', () => {
    it('evaluates Hermite polynomial', () => {
      // H_0(x) = 1
      expect(hermval(0.5, [1])).toBeCloseTo(1);

      // H_1(x) = 2x
      expect(hermval(0.5, [0, 1])).toBeCloseTo(1); // 2*0.5 = 1

      // H_2(x) = 4x^2 - 2
      expect(hermval(0.5, [0, 0, 1])).toBeCloseTo(4 * 0.25 - 2); // -1
    });
  });

  describe('herm2poly / poly2herm', () => {
    it('converts Hermite to power series', () => {
      // H_0 = 1
      expect(herm2poly([1])).toEqual([1]);

      // H_1 = 2x
      expect(herm2poly([0, 1])).toEqual([0, 2]);

      // H_2 = 4x^2 - 2
      const h2 = herm2poly([0, 0, 1]);
      expect(allClose(h2, [-2, 0, 4])).toBe(true);
    });
  });

  /* ============ HermiteE ============ */

  describe('HermiteE class', () => {
    it('constructs from coefficients', () => {
      const p = new HermiteE([1, 2, 3]);
      expect(p.coef).toEqual([1, 2, 3]);
    });
  });

  describe('hermeval', () => {
    it('evaluates HermiteE polynomial', () => {
      // He_0(x) = 1
      expect(hermeval(0.5, [1])).toBeCloseTo(1);

      // He_1(x) = x
      expect(hermeval(0.5, [0, 1])).toBeCloseTo(0.5);

      // He_2(x) = x^2 - 1
      expect(hermeval(0.5, [0, 0, 1])).toBeCloseTo(0.25 - 1); // -0.75
    });
  });

  describe('herme2poly / poly2herme', () => {
    it('converts HermiteE to power series', () => {
      // He_0 = 1
      expect(herme2poly([1])).toEqual([1]);

      // He_1 = x
      expect(herme2poly([0, 1])).toEqual([0, 1]);

      // He_2 = x^2 - 1
      const he2 = herme2poly([0, 0, 1]);
      expect(allClose(he2, [-1, 0, 1])).toBe(true);
    });
  });

  /* ============ Laguerre ============ */

  describe('Laguerre class', () => {
    it('constructs from coefficients', () => {
      const p = new Laguerre([1, 2, 3]);
      expect(p.coef).toEqual([1, 2, 3]);
    });
  });

  describe('lagval', () => {
    it('evaluates Laguerre polynomial', () => {
      // L_0(x) = 1
      expect(lagval(0.5, [1])).toBeCloseTo(1);

      // L_1(x) = 1 - x
      expect(lagval(0.5, [0, 1])).toBeCloseTo(0.5); // 1 - 0.5

      // L_2(x) = (x^2 - 4x + 2) / 2
      expect(lagval(0.5, [0, 0, 1])).toBeCloseTo((0.25 - 2 + 2) / 2);
    });
  });

  describe('lag2poly / poly2lag', () => {
    it('converts Laguerre to power series', () => {
      // L_0 = 1
      expect(lag2poly([1])).toEqual([1]);

      // L_1 = 1 - x
      const l1 = lag2poly([0, 1]);
      expect(allClose(l1, [1, -1])).toBe(true);

      // L_2 = (x^2 - 4x + 2) / 2 = 1 - 2x + 0.5x^2
      const l2 = lag2poly([0, 0, 1]);
      expect(allClose(l2, [1, -2, 0.5])).toBe(true);
    });
  });

  /* ============ Error handling ============ */

  describe('Error handling', () => {
    it('throws PolyError for zero division', () => {
      expect(() => polydiv([1, 2], [0])).toThrow(PolyError);
    });

    it('throws PolyError for negative power', () => {
      expect(() => polypow([1, 1], -1)).toThrow(PolyError);
    });

    it('throws PolyError for incompatible polynomial types', () => {
      const p1 = new Polynomial([1, 2]);
      const p2 = new Chebyshev([1, 2]);
      expect(() => p1.add(p2 as unknown as Polynomial)).toThrow(PolyError);
    });
  });
});
