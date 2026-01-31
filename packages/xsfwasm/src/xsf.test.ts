/**
 * Tests for XSF special functions
 */

import { describe, test, expect, beforeAll } from 'vitest';
import {
  loadXSFModule,
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
  chdtr,
  fdtr,
  pdtr,
  bdtr,
  kolmogorov,
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
  });

  test('jv', () => {
    expect(jv(0, 0)).toBeCloseTo(1, 10);
  });

  test('y0', () => {
    expect(y0(1)).toBeCloseTo(0.0882569642, 8);
  });

  test('i0', () => {
    expect(i0(0)).toBeCloseTo(1, 10);
  });

  test('k0', () => {
    expect(k0(1)).toBeCloseTo(0.4210244382, 8);
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
  });

  test('bei', () => {
    expect(bei(0)).toBeCloseTo(0, 10);
  });

  test('ker', () => {
    expect(ker(1)).toBeCloseTo(0.2867062088, 6);
  });

  test('kei', () => {
    expect(kei(1)).toBeCloseTo(-0.49499463, 5);
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
  });

  test('ndtri (normal quantile)', () => {
    expect(ndtri(0.5)).toBeCloseTo(0, 10);
    expect(ndtri(0.975)).toBeCloseTo(1.96, 2);
  });

  test('chdtr (chi-square CDF)', () => {
    expect(chdtr(1, 0)).toBeCloseTo(0, 10);
  });

  test('fdtr (F-distribution CDF)', () => {
    expect(fdtr(1, 1, 0)).toBeCloseTo(0, 10);
  });

  test('pdtr (Poisson CDF)', () => {
    expect(pdtr(0, 1)).toBeCloseTo(Math.exp(-1), 8);
  });

  test('bdtr (binomial CDF)', () => {
    expect(bdtr(0, 2, 0.5)).toBeCloseTo(0.25, 8);
  });

  test('kolmogorov', () => {
    expect(kolmogorov(0)).toBeCloseTo(1, 10);
  });

  test('owens_t', () => {
    expect(owens_t(0, 0)).toBeCloseTo(0, 10);
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
});
