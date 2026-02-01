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

// ============================================================
// MODULE MANAGEMENT
// ============================================================

export {
  loadXSFModule,
  getXSFModule,
  isXSFLoaded,
  resetXSFModule,
  configureXSF,
  type XSFLoadConfig,
} from './ts/core/loader.js';

// Types (re-export low-level types for advanced users)
export type { XSFModule, XSFModuleFactory } from './ts/core/types.js';

// ============================================================
// DIRECT FUNCTION EXPORTS (most common use case)
// ============================================================

// Gamma functions
export { gamma, gammaln, rgamma, beta, betaln, digamma } from './ts/gamma/index.js';

// Error functions
export { erf, erfc, erfcx, erfi } from './ts/error/index.js';

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
  // Spherical Bessel
  sphericalJn,
  sphericalYn,
  sphericalIn,
  sphericalKn,
  // Kelvin functions
  ber,
  bei,
  ker,
  kei,
  berp,
  beip,
  kerp,
  keip,
} from './ts/bessel/index.js';

// Special functions
export {
  airy,
  type AiryResult,
  type AiryArrayResult,
  struveH,
  struveL,
  legendreP,
  hyp2f1,
} from './ts/special/index.js';

// Integral functions
export {
  ellipk,
  ellipe,
  ellipkinc,
  ellipeinc,
  ellipj,
  type EllipjResult,
  type EllipjArrayResult,
  exp1,
  expi,
  fresnel,
  type FresnelResult,
  type FresnelArrayResult,
  sici,
  type SiCiResult,
  type SiCiArrayResult,
  shichi,
  type ShiChiResult,
  type ShiChiArrayResult,
} from './ts/integrals/index.js';

// Combinatorial functions
export { binom, binomExact, poch, permExact } from './ts/combinatorics/index.js';

// Analysis functions
export { zeta, zetac, lambertw } from './ts/analysis/index.js';

// Statistical distributions
export {
  // Normal distribution
  ndtr,
  ndtri,
  logNdtr,
  // Chi-square distribution
  chdtr,
  chdtrc,
  chdtri,
  // F-distribution
  fdtr,
  fdtrc,
  fdtri,
  // Gamma distribution
  gdtr,
  gdtrc,
  // Poisson distribution
  pdtr,
  pdtrc,
  pdtri,
  // Binomial distribution
  bdtr,
  bdtrc,
  bdtri,
  // Negative binomial distribution
  nbdtr,
  nbdtrc,
  nbdtri,
  // Kolmogorov-Smirnov
  kolmogorov,
  kolmogi,
  smirnov,
  smirnovi,
  // Owen's T function
  owensT,
} from './ts/stats/index.js';

// ============================================================
// NAMESPACE EXPORTS (for organized module access)
// ============================================================

/** Gamma function family (gamma, beta, digamma) */
export * as gammaFns from './ts/gamma/index.js';

/** Bessel function family (J, Y, I, K, spherical, Kelvin) */
export * as besselFns from './ts/bessel/index.js';

/** Error functions (erf, erfc, erfcx, erfi) */
export * as errorFns from './ts/error/index.js';

/** Special functions (Airy, Struve, Legendre, hypergeometric) */
export * as specialFns from './ts/special/index.js';

/** Integral functions (elliptic, exponential, Fresnel, sine/cosine) */
export * as integralFns from './ts/integrals/index.js';

/** Combinatorial functions (binomial, Pochhammer, permutations) */
export * as combinatoricsFns from './ts/combinatorics/index.js';

/** Analysis functions (zeta, Lambert W) */
export * as analysisFns from './ts/analysis/index.js';

/** Statistical distributions */
export * as statsFns from './ts/stats/index.js';
