/**
 * XSFWasm - XSF WebAssembly Module
 *
 * This package provides a WebAssembly build of XSF (Extended Special Functions),
 * a C++ library that implements special mathematical functions derived from
 * SciPy's Cephes, AMOS, and other sources.
 *
 * XSF provides routines for:
 * - Gamma and related functions (gamma, gammaln, rgamma, digamma, beta, betaln)
 * - Bessel functions (J, Y, I, K variants, spherical Bessel)
 * - Error functions (erf, erfc, erfcx, erfi)
 * - Combinatorial functions (binom, poch, permutations)
 * - Airy functions (Ai, Bi and derivatives)
 * - Elliptic integrals (complete and incomplete)
 * - Exponential integrals (E₁, Ei)
 * - Fresnel integrals (S, C)
 * - Hypergeometric function (₂F₁)
 * - Kelvin functions (ber, bei, ker, kei)
 * - Lambert W function
 * - Sine/cosine integrals (Si, Ci, Shi, Chi)
 * - Struve functions (H, L)
 * - Zeta function
 * - Legendre polynomials
 * - Statistical distributions (normal, chi-square, F, gamma, Poisson, binomial, etc.)
 *
 * @example
 * ```typescript
 * import { loadXSFModule, gamma, erf, j0, ndtr, airy } from 'xsfwasm';
 *
 * // Load the WASM module once
 * await loadXSFModule();
 *
 * // Use high-level API functions (synchronous after load)
 * gamma(5);           // 24 (= 4!)
 * erf(1);             // 0.8427...
 * j0(2.405);          // ~0 (first zero of J₀)
 * ndtr(0);            // 0.5 (normal CDF at 0)
 *
 * // Multi-value returns
 * const { ai, bi } = airy(0);  // Airy functions at x=0
 *
 * // Vectorized operations
 * gamma([1, 2, 3, 4, 5]);  // Float64Array([1, 1, 2, 6, 24])
 * ```
 *
 * @packageDocumentation
 */

// Module management
export {
  loadXSFModule,
  getXSFModule,
  isXSFLoaded,
  resetXSFModule,
  configureXSF,
  type XSFLoadConfig,
} from './ts/loader.js';

// Gamma functions
export { gamma, gammaln, rgamma } from './ts/gamma.js';

// Digamma function
export { digamma } from './ts/digamma.js';

// Beta functions
export { beta, betaln } from './ts/beta.js';

// Error functions
export { erf, erfc, erfcx, erfi } from './ts/error.js';

// Bessel functions
export {
  // First kind (J)
  j0,
  j1,
  jv,
  // Second kind (Y) - Neumann
  y0,
  y1,
  yv,
  // Modified first kind (I)
  i0,
  i1,
  iv,
  // Modified second kind (K) - MacDonald
  k0,
  k1,
} from './ts/bessel.js';

// Spherical Bessel functions
export {
  spherical_jn,
  spherical_yn,
  spherical_in,
  spherical_kn,
} from './ts/sph_bessel.js';

// Combinatorial functions
export { binom, binomExact, poch, permExact } from './ts/combinatorial.js';

// Airy functions
export { airy, type AiryResult } from './ts/airy.js';

// Elliptic integrals
export {
  ellipk,
  ellipe,
  ellipkinc,
  ellipeinc,
  ellipj,
  type EllipjResult,
} from './ts/elliptic.js';

// Exponential integrals
export { exp1, expi } from './ts/expint.js';

// Fresnel integrals
export { fresnel, type FresnelResult } from './ts/fresnel.js';

// Hypergeometric function
export { hyp2f1 } from './ts/hyper.js';

// Kelvin functions
export {
  ber,
  bei,
  ker,
  kei,
  berp,
  beip,
  kerp,
  keip,
} from './ts/kelvin.js';

// Lambert W function
export { lambertw } from './ts/lambertw.js';

// Sine/cosine integrals
export { sici, shichi, type SiCiResult, type ShiChiResult } from './ts/sici.js';

// Struve functions
export { struve_h, struve_l } from './ts/struve.js';

// Zeta function
export { zeta, zetac } from './ts/zeta.js';

// Legendre polynomials
export { legendre_p } from './ts/legendre.js';

// Statistical distributions
export {
  // Normal distribution
  ndtr,
  ndtri,
  log_ndtr,
  // Chi-square distribution
  chdtr,
  chdtri,
  // F-distribution
  fdtr,
  fdtri,
  // Gamma distribution
  gdtr,
  gdtrc,
  // Poisson distribution
  pdtr,
  pdtri,
  // Binomial distribution
  bdtr,
  bdtri,
  // Negative binomial distribution
  nbdtr,
  nbdtri,
  // Kolmogorov-Smirnov
  kolmogorov,
  kolmogi,
  smirnov,
  smirnovi,
  // Owen's T function
  owens_t,
} from './ts/stats.js';

// Types (re-export low-level types for advanced users)
export type { XSFModule, XSFModuleFactory } from './ts/types.js';
