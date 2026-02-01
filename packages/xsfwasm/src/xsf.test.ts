/**
 * Tests for XSF special functions
 */

import { describe, test, expect, beforeAll } from 'vitest';
import {
  loadXSFModule,
  getXSFModule,
  isXSFLoaded,
  resetXSFModule,
  // Gamma functions
  gamma,
  gammaln,
  rgamma,
  digamma,
  // Beta functions
  beta,
  betaln,
  // Error functions
  erf,
  erfc,
  erfcx,
  erfi,
  // Bessel functions
  j0,
  j1,
  jv,
  y0,
  y1,
  yv,
  i0,
  i1,
  iv,
  k0,
  k1,
  // Spherical Bessel
  spherical_jn,
  spherical_yn,
  spherical_in,
  spherical_kn,
  // Combinatorial
  binom,
  binomExact,
  poch,
  permExact,
  // Airy
  airy,
  // Elliptic
  ellipk,
  ellipe,
  ellipkinc,
  ellipeinc,
  ellipj,
  // Exponential integrals
  exp1,
  expi,
  // Fresnel
  fresnel,
  // Hypergeometric
  hyp2f1,
  // Kelvin
  ber,
  bei,
  ker,
  kei,
  berp,
  beip,
  kerp,
  keip,
  // Lambert W
  lambertw,
  // Sine/cosine integrals
  sici,
  shichi,
  // Struve
  struve_h,
  struve_l,
  // Zeta
  zeta,
  zetac,
  // Legendre
  legendre_p,
  // Statistics
  ndtr,
  ndtri,
  log_ndtr,
  chdtr,
  chdtri,
  fdtr,
  fdtri,
  gdtr,
  gdtrc,
  pdtr,
  pdtri,
  bdtr,
  bdtri,
  nbdtr,
  nbdtri,
  kolmogorov,
  kolmogi,
  smirnov,
  smirnovi,
  owens_t,
} from './index.js';

beforeAll(async () => {
  await loadXSFModule();
});

describe('Gamma functions', () => {
  test('gamma', () => {
    expect(gamma(5)).toBeCloseTo(24, 10);
    expect(gamma(0.5)).toBeCloseTo(Math.sqrt(Math.PI), 10);
  });

  test('gammaln', () => {
    expect(gammaln(5)).toBeCloseTo(Math.log(24), 10);
  });

  test('rgamma', () => {
    expect(rgamma(5)).toBeCloseTo(1 / 24, 10);
  });

  test('digamma', () => {
    expect(digamma(1)).toBeCloseTo(-0.5772156649, 8);
  });
});

describe('Beta functions', () => {
  test('beta', () => {
    expect(beta(1, 1)).toBeCloseTo(1, 10);
    expect(beta(2, 3)).toBeCloseTo(1 / 12, 10);
  });

  test('betaln', () => {
    expect(betaln(2, 3)).toBeCloseTo(Math.log(1 / 12), 10);
  });
});

describe('Error functions', () => {
  test('erf', () => {
    expect(erf(0)).toBeCloseTo(0, 10);
    expect(erf(1)).toBeCloseTo(0.8427007929, 8);
  });

  test('erfc', () => {
    expect(erfc(0)).toBeCloseTo(1, 10);
  });

  test('erfcx', () => {
    expect(erfcx(0)).toBeCloseTo(1, 10);
  });

  test('erfi', () => {
    expect(erfi(0)).toBeCloseTo(0, 10);
  });
});

describe('Bessel functions', () => {
  test('j0', () => {
    expect(j0(0)).toBeCloseTo(1, 10);
    expect(j0(2.4048)).toBeCloseTo(0, 3); // first zero
  });

  test('j1', () => {
    expect(j1(0)).toBeCloseTo(0, 10);
    expect(j1(1)).toBeCloseTo(0.4400505857, 8);
  });

  test('jv', () => {
    expect(jv(0, 0)).toBeCloseTo(1, 10);
    expect(jv(0.5, 1)).toBeCloseTo(0.6713967071, 8);
    expect(jv(2, 3)).toBeCloseTo(0.4860912606, 8);
  });

  test('y0', () => {
    expect(y0(1)).toBeCloseTo(0.0882569642, 8);
  });

  test('y1', () => {
    expect(y1(1)).toBeCloseTo(-0.7812128213, 8);
    expect(y1(2)).toBeCloseTo(-0.1070324315, 8);
  });

  test('yv', () => {
    expect(yv(0, 1)).toBeCloseTo(y0(1), 8);
    expect(yv(1, 1)).toBeCloseTo(y1(1), 8);
    expect(yv(0.5, 1)).toBeCloseTo(-0.4310988680, 8);
  });

  test('i0', () => {
    expect(i0(0)).toBeCloseTo(1, 10);
    expect(i0(1)).toBeCloseTo(1.2660658778, 8);
  });

  test('i1', () => {
    expect(i1(0)).toBeCloseTo(0, 10);
    expect(i1(1)).toBeCloseTo(0.5651591040, 8);
    expect(i1(2)).toBeCloseTo(1.5906368546, 8);
  });

  test('iv', () => {
    expect(iv(0, 1)).toBeCloseTo(i0(1), 8);
    expect(iv(1, 1)).toBeCloseTo(i1(1), 8);
    expect(iv(0.5, 1)).toBeCloseTo(0.9376748882, 8);
    expect(iv(2, 2)).toBeCloseTo(0.6889484476, 8);
  });

  test('k0', () => {
    expect(k0(1)).toBeCloseTo(0.4210244382, 8);
    expect(k0(2)).toBeCloseTo(0.1138938727, 8);
  });

  test('k1', () => {
    expect(k1(1)).toBeCloseTo(0.6019072302, 8);
    expect(k1(2)).toBeCloseTo(0.1398658819, 8);
    expect(k1(0.5)).toBeCloseTo(1.6564411200, 8);
  });
});

describe('Spherical Bessel functions', () => {
  test('spherical_jn', () => {
    // j_0(x) = sin(x)/x
    expect(spherical_jn(0, 1)).toBeCloseTo(Math.sin(1), 8);
  });

  test('spherical_yn', () => {
    // y_0(x) = -cos(x)/x
    expect(spherical_yn(0, 1)).toBeCloseTo(-Math.cos(1), 8);
  });

  test('spherical_in', () => {
    // i_0(x) = sinh(x)/x
    expect(spherical_in(0, 1)).toBeCloseTo(Math.sinh(1), 8);
  });

  test('spherical_kn', () => {
    // k_0(x) = (π/2) * e^(-x) / x
    expect(spherical_kn(0, 1)).toBeCloseTo((Math.PI / 2) * Math.exp(-1), 8);
  });
});

describe('Combinatorial functions', () => {
  test('binom', () => {
    expect(binom(5, 2)).toBeCloseTo(10, 10);
    expect(binom(10, 5)).toBeCloseTo(252, 10);
  });

  test('binomExact', () => {
    expect(binomExact(5, 2)).toBeCloseTo(10, 10);
  });

  test('poch', () => {
    expect(poch(1, 5)).toBeCloseTo(120, 10); // 5!
  });

  test('permExact', () => {
    expect(permExact(5, 2)).toBeCloseTo(20, 10);
  });
});

describe('Airy functions', () => {
  test('airy', () => {
    const result = airy(0);
    expect(result.ai).toBeCloseTo(0.3550280538, 8);
    expect(result.bi).toBeCloseTo(0.6149266274, 8);
  });
});

describe('Elliptic integrals', () => {
  test('ellipk', () => {
    expect(ellipk(0)).toBeCloseTo(Math.PI / 2, 10);
  });

  test('ellipe', () => {
    expect(ellipe(0)).toBeCloseTo(Math.PI / 2, 10);
  });

  test('ellipkinc', () => {
    expect(ellipkinc(0, 0.5)).toBeCloseTo(0, 10);
  });

  test('ellipeinc', () => {
    expect(ellipeinc(0, 0.5)).toBeCloseTo(0, 10);
  });

  test('ellipj', () => {
    const result = ellipj(0, 0.5);
    expect(result.sn).toBeCloseTo(0, 10);
    expect(result.cn).toBeCloseTo(1, 10);
    expect(result.dn).toBeCloseTo(1, 10);
  });
});

describe('Exponential integrals', () => {
  test('exp1', () => {
    expect(exp1(1)).toBeCloseTo(0.2193839343, 8);
  });

  test('expi', () => {
    expect(expi(1)).toBeCloseTo(1.8951178163, 8);
  });
});

describe('Fresnel integrals', () => {
  test('fresnel', () => {
    const result = fresnel(0);
    expect(result.s).toBeCloseTo(0, 10);
    expect(result.c).toBeCloseTo(0, 10);
  });
});

describe('Hypergeometric function', () => {
  test('hyp2f1', () => {
    // 2F1(1, 1, 2, x) = -ln(1-x)/x
    expect(hyp2f1(1, 1, 2, 0.5)).toBeCloseTo(-Math.log(0.5) / 0.5, 6);
  });
});

describe('Kelvin functions', () => {
  test('ber', () => {
    expect(ber(0)).toBeCloseTo(1, 10);
    expect(ber(1)).toBeCloseTo(0.9843817812, 8);
  });

  test('bei', () => {
    expect(bei(0)).toBeCloseTo(0, 10);
    expect(bei(1)).toBeCloseTo(0.2495660400, 6);
  });

  test('ker', () => {
    expect(ker(1)).toBeCloseTo(0.2867062088, 6);
  });

  test('kei', () => {
    expect(kei(1)).toBeCloseTo(-0.49499463, 5);
  });

  test('berp (derivative of ber)', () => {
    expect(berp(0)).toBeCloseTo(0, 10);
    expect(berp(1)).toBeCloseTo(-0.0624457521, 6);
  });

  test('beip (derivative of bei)', () => {
    expect(beip(0)).toBeCloseTo(0, 10);
    expect(beip(1)).toBeCloseTo(0.4973965114, 6);
  });

  test('kerp (derivative of ker)', () => {
    expect(kerp(1)).toBeCloseTo(-0.6946038911, 6);
  });

  test('keip (derivative of kei)', () => {
    expect(keip(1)).toBeCloseTo(0.3523699133, 6);
  });
});

describe('Lambert W function', () => {
  test('lambertw principal branch', () => {
    expect(lambertw(0, 0)).toBeCloseTo(0, 10);
    expect(lambertw(Math.E, 0)).toBeCloseTo(1, 10);
  });
});

describe('Sine/cosine integrals', () => {
  test('sici', () => {
    const result = sici(1);
    expect(result.si).toBeCloseTo(0.9460830703, 8);
    expect(result.ci).toBeCloseTo(0.3374039229, 8);
  });

  test('shichi', () => {
    const result = shichi(1);
    expect(result.shi).toBeCloseTo(1.0572508753, 8);
    expect(result.chi).toBeCloseTo(0.8378669409, 8);
  });
});

describe('Struve functions', () => {
  test('struve_h', () => {
    expect(struve_h(0, 0)).toBeCloseTo(0, 10);
  });

  test('struve_l', () => {
    expect(struve_l(0, 0)).toBeCloseTo(0, 10);
  });
});

describe('Zeta function', () => {
  test('zeta', () => {
    expect(zeta(2)).toBeCloseTo(Math.PI ** 2 / 6, 10);
    expect(zeta(4)).toBeCloseTo(Math.PI ** 4 / 90, 10);
  });

  test('zetac', () => {
    expect(zetac(2)).toBeCloseTo(Math.PI ** 2 / 6 - 1, 10);
  });
});

describe('Legendre polynomials', () => {
  test('legendre_p', () => {
    expect(legendre_p(0, 0.5)).toBeCloseTo(1, 10);
    expect(legendre_p(1, 0.5)).toBeCloseTo(0.5, 10);
    expect(legendre_p(2, 0.5)).toBeCloseTo(-0.125, 10);
  });
});

describe('Statistical distributions', () => {
  test('ndtr (normal CDF)', () => {
    expect(ndtr(0)).toBeCloseTo(0.5, 10);
    expect(ndtr(1)).toBeCloseTo(0.8413447460, 8);
    expect(ndtr(-1)).toBeCloseTo(0.1586552540, 8);
  });

  test('ndtri (normal quantile)', () => {
    expect(ndtri(0.5)).toBeCloseTo(0, 10);
    expect(ndtri(0.975)).toBeCloseTo(1.96, 2);
    expect(ndtri(0.025)).toBeCloseTo(-1.96, 2);
  });

  test('log_ndtr (log normal CDF)', () => {
    expect(log_ndtr(0)).toBeCloseTo(Math.log(0.5), 10);
    // For very negative x, log_ndtr is more accurate than log(ndtr(x))
    expect(log_ndtr(-5)).toBeCloseTo(-15.0649983939, 5);
  });

  test('chdtr (chi-square CDF)', () => {
    expect(chdtr(1, 0)).toBeCloseTo(0, 10);
    expect(chdtr(2, 2)).toBeCloseTo(0.6321205588, 8);
    expect(chdtr(5, 11.07)).toBeCloseTo(0.95, 2);
  });

  test('chdtri (chi-square inverse CDF)', () => {
    // chdtri returns x such that chdtrc(df, x) = p (it uses complementary CDF)
    const x = chdtri(2, 0.05); // pass 1-p since it uses survival function
    expect(chdtr(2, x)).toBeCloseTo(0.95, 5);
  });

  test('fdtr (F-distribution CDF)', () => {
    expect(fdtr(1, 1, 0)).toBeCloseTo(0, 10);
    expect(fdtr(5, 10, 3.33)).toBeCloseTo(0.95, 2);
  });

  test('fdtri (F-distribution inverse CDF)', () => {
    // fdtri returns x such that fdtr(dfn, dfd, x) = p
    const x = fdtri(5, 10, 0.95);
    expect(fdtr(5, 10, x)).toBeCloseTo(0.95, 5);
  });

  test('gdtr (gamma distribution CDF)', () => {
    // gdtr(a, b, x) = P(X <= x) for gamma distribution
    expect(gdtr(1, 1, 0)).toBeCloseTo(0, 10);
    expect(gdtr(1, 1, 1)).toBeCloseTo(0.6321205588, 8); // 1 - e^(-1)
  });

  test('gdtrc (gamma distribution survival)', () => {
    // gdtrc = 1 - gdtr (survival function)
    expect(gdtrc(1, 1, 1)).toBeCloseTo(1 - 0.6321205588, 8);
    expect(gdtr(1, 1, 1) + gdtrc(1, 1, 1)).toBeCloseTo(1, 10);
  });

  test('pdtr (Poisson CDF)', () => {
    expect(pdtr(0, 1)).toBeCloseTo(Math.exp(-1), 8);
    expect(pdtr(1, 1)).toBeCloseTo(2 * Math.exp(-1), 8);
  });

  test('pdtri (Poisson inverse)', () => {
    // Returns m such that pdtr(k, m) = p
    const m = pdtri(5, 0.5);
    expect(m).toBeGreaterThan(0);
  });

  test('bdtr (binomial CDF)', () => {
    expect(bdtr(0, 2, 0.5)).toBeCloseTo(0.25, 8);
    expect(bdtr(1, 2, 0.5)).toBeCloseTo(0.75, 8);
    expect(bdtr(2, 2, 0.5)).toBeCloseTo(1.0, 8);
  });

  test('bdtri (binomial inverse)', () => {
    // Returns p such that bdtr(k, n, p) = y
    const p = bdtri(1, 2, 0.75);
    expect(bdtr(1, 2, p)).toBeCloseTo(0.75, 5);
  });

  test('nbdtr (negative binomial CDF)', () => {
    // P(X <= k) where X is negative binomial
    expect(nbdtr(0, 1, 0.5)).toBeCloseTo(0.5, 8);
  });

  test('nbdtri (negative binomial inverse)', () => {
    // Returns p such that nbdtr(k, n, p) = y
    const p = nbdtri(5, 3, 0.5);
    expect(nbdtr(5, 3, p)).toBeCloseTo(0.5, 5);
  });

  test('kolmogorov', () => {
    expect(kolmogorov(0)).toBeCloseTo(1, 10);
    expect(kolmogorov(1)).toBeCloseTo(0.2699996716, 6);
  });

  test('kolmogi (Kolmogorov inverse)', () => {
    const x = kolmogi(0.5);
    expect(kolmogorov(x)).toBeCloseTo(0.5, 5);
  });

  test('smirnov', () => {
    // Smirnov one-sided distribution - returns survival function P(D_n > x)
    expect(smirnov(10, 0.5)).toBeCloseTo(0.003888705, 6);
  });

  test('smirnovi (Smirnov inverse)', () => {
    const x = smirnovi(10, 0.95);
    expect(smirnov(10, x)).toBeCloseTo(0.95, 5);
  });

  test('owens_t', () => {
    expect(owens_t(0, 0)).toBeCloseTo(0, 10);
    expect(owens_t(1, 1)).toBeCloseTo(0.0667418821, 6);
  });
});

describe('Vectorized operations', () => {
  test('gamma array', () => {
    const result = gamma([1, 2, 3, 4, 5]);
    expect(result).toBeInstanceOf(Float64Array);
    expect(result[0]).toBeCloseTo(1, 10);
    expect(result[1]).toBeCloseTo(1, 10);
    expect(result[2]).toBeCloseTo(2, 10);
    expect(result[3]).toBeCloseTo(6, 10);
    expect(result[4]).toBeCloseTo(24, 10);
  });

  test('erf array', () => {
    const result = erf([0, 1, -1]);
    expect(result).toBeInstanceOf(Float64Array);
    expect(result[0]).toBeCloseTo(0, 10);
    expect(result[1]).toBeCloseTo(0.8427007929, 8);
    expect(result[2]).toBeCloseTo(-0.8427007929, 8);
  });

  test('airy array', () => {
    const result = airy([0, 1]);
    expect(result.ai).toBeInstanceOf(Float64Array);
    expect(result.ai[0]).toBeCloseTo(0.3550280538, 8);
  });

  test('Float64Array input', () => {
    const input = new Float64Array([1, 2, 3, 4, 5]);
    const result = gamma(input);
    expect(result).toBeInstanceOf(Float64Array);
    expect(result[0]).toBeCloseTo(1, 10);
    expect(result[4]).toBeCloseTo(24, 10);
  });

  test('Bessel jv with array x', () => {
    const result = jv(0, [0, 1, 2]);
    expect(result).toBeInstanceOf(Float64Array);
    expect(result[0]).toBeCloseTo(1, 10);
    expect(result[1]).toBeCloseTo(j0(1), 8);
    expect(result[2]).toBeCloseTo(j0(2), 8);
  });

  test('Bessel jv with array v and x', () => {
    const result = jv([0, 1, 2], [1, 1, 1]);
    expect(result).toBeInstanceOf(Float64Array);
    expect(result[0]).toBeCloseTo(j0(1), 8);
    expect(result[1]).toBeCloseTo(j1(1), 8);
  });

  test('ndtr array', () => {
    const result = ndtr([-1, 0, 1]);
    expect(result).toBeInstanceOf(Float64Array);
    expect(result[0]).toBeCloseTo(0.1586552540, 8);
    expect(result[1]).toBeCloseTo(0.5, 10);
    expect(result[2]).toBeCloseTo(0.8413447460, 8);
  });

  test('ellipj array', () => {
    const result = ellipj([0, 1], 0.5);
    expect(result.sn).toBeInstanceOf(Float64Array);
    expect(result.cn).toBeInstanceOf(Float64Array);
    expect(result.dn).toBeInstanceOf(Float64Array);
    expect(result.sn[0]).toBeCloseTo(0, 10);
    expect(result.cn[0]).toBeCloseTo(1, 10);
  });

  test('fresnel array', () => {
    const result = fresnel([0, 1, 2]);
    expect(result.s).toBeInstanceOf(Float64Array);
    expect(result.c).toBeInstanceOf(Float64Array);
    expect(result.s[0]).toBeCloseTo(0, 10);
    expect(result.c[0]).toBeCloseTo(0, 10);
  });

  test('sici array', () => {
    const result = sici([1, 2]);
    expect(result.si).toBeInstanceOf(Float64Array);
    expect(result.ci).toBeInstanceOf(Float64Array);
    expect(result.si[0]).toBeCloseTo(0.9460830703, 8);
  });

  test('Kelvin functions array', () => {
    const result = ber([0, 1, 2]);
    expect(result).toBeInstanceOf(Float64Array);
    expect(result[0]).toBeCloseTo(1, 10);
    expect(result[1]).toBeCloseTo(0.9843817812, 6);
  });
});

describe('Module management', () => {
  test('isXSFLoaded returns true after load', () => {
    expect(isXSFLoaded()).toBe(true);
  });

  test('getXSFModule returns the module', () => {
    const module = getXSFModule();
    expect(module).toBeDefined();
    expect(typeof module._wasm_gamma).toBe('function');
  });

  test('loadXSFModule returns same module on subsequent calls', async () => {
    const module1 = await loadXSFModule();
    const module2 = await loadXSFModule();
    expect(module1).toBe(module2);
  });
});

describe('Edge cases and special values', () => {
  test('gamma at negative integers is NaN or Infinity', () => {
    // Gamma is undefined at non-positive integers
    expect(gamma(0)).toBe(Infinity);
    expect(Number.isNaN(gamma(-1))).toBe(true);
  });

  test('gammaln handles large values', () => {
    // Stirling approximation: gammaln(n) ≈ n*log(n) - n for large n
    const n = 100;
    expect(gammaln(n)).toBeCloseTo(359.1342053695, 5);
  });

  test('erf symmetry', () => {
    // erf(-x) = -erf(x)
    expect(erf(-2)).toBeCloseTo(-erf(2), 10);
  });

  test('erfc complement', () => {
    // erfc(x) = 1 - erf(x)
    expect(erfc(1)).toBeCloseTo(1 - erf(1), 10);
  });

  test('Bessel recurrence relation', () => {
    // J_{n-1}(x) + J_{n+1}(x) = (2n/x) * J_n(x)
    const x = 2;
    const n = 2;
    const lhs = jv(n - 1, x) + jv(n + 1, x);
    const rhs = ((2 * n) / x) * jv(n, x);
    expect(lhs).toBeCloseTo(rhs, 8);
  });

  test('spherical Bessel j_0 relation', () => {
    // j_0(x) = sin(x)/x
    const x = 2;
    expect(spherical_jn(0, x)).toBeCloseTo(Math.sin(x) / x, 10);
  });

  test('elliptic integrals at m=0', () => {
    // K(0) = E(0) = π/2
    expect(ellipk(0)).toBeCloseTo(Math.PI / 2, 10);
    expect(ellipe(0)).toBeCloseTo(Math.PI / 2, 10);
  });

  test('zeta at even positive integers', () => {
    // ζ(2) = π²/6, ζ(4) = π⁴/90
    expect(zeta(2)).toBeCloseTo(Math.PI ** 2 / 6, 10);
    expect(zeta(4)).toBeCloseTo(Math.PI ** 4 / 90, 10);
    expect(zeta(6)).toBeCloseTo(Math.PI ** 6 / 945, 8);
  });

  test('Lambert W at special points', () => {
    // W(0) = 0, W(e) = 1
    expect(lambertw(0, 0)).toBeCloseTo(0, 10);
    expect(lambertw(Math.E, 0)).toBeCloseTo(1, 10);
    // W(x) * e^W(x) = x
    const x = 5;
    const w = lambertw(x, 0);
    expect(w * Math.exp(w)).toBeCloseTo(x, 10);
  });

  test('beta relation to gamma', () => {
    // B(a, b) = Γ(a)Γ(b)/Γ(a+b)
    const a = 2,
      b = 3;
    expect(beta(a, b)).toBeCloseTo((gamma(a) * gamma(b)) / gamma(a + b), 10);
  });

  test('binomial coefficient symmetry', () => {
    // C(n, k) = C(n, n-k)
    expect(binom(10, 3)).toBeCloseTo(binom(10, 7), 10);
  });

  test('Pochhammer symbol relation', () => {
    // (x)_n = Γ(x+n)/Γ(x)
    const x = 3,
      n = 4;
    expect(poch(x, n)).toBeCloseTo(gamma(x + n) / gamma(x), 8);
  });
});
