/**
 * Statistical distribution functions - high-level API
 */

import { getXSFModule, isXSFLoaded } from './loader.js';
import type { XSFModule } from './types.js';

type ArrayInput = number[] | Float64Array;

function ensureLoaded(): XSFModule {
  if (!isXSFLoaded()) {
    throw new Error(
      'XSF WASM module not loaded. Call "await loadXSFModule()" before using special functions.'
    );
  }
  return getXSFModule();
}

// =============================================================================
// NORMAL DISTRIBUTION
// =============================================================================

/**
 * Normal distribution CDF (cumulative distribution function).
 *
 * Φ(x) = (1/√(2π)) ∫_{-∞}^x e^(-t²/2) dt
 *
 * @param x - Input value
 * @returns P(X ≤ x) for standard normal X
 *
 * @example
 * ```ts
 * ndtr(0);     // 0.5
 * ndtr(1.96);  // ≈ 0.975
 * ```
 */
export function ndtr(x: number): number;
export function ndtr(x: ArrayInput): Float64Array;
export function ndtr(x: number | ArrayInput): number | Float64Array {
  const xsf = ensureLoaded();

  if (typeof x === 'number') {
    return xsf._wasm_ndtr(x);
  }

  const arr = x instanceof Float64Array ? x : new Float64Array(x);
  const result = new Float64Array(arr.length);
  for (let i = 0; i < arr.length; i++) {
    result[i] = xsf._wasm_ndtr(arr[i]);
  }
  return result;
}

/**
 * Normal distribution inverse CDF (quantile function).
 *
 * Returns x such that Φ(x) = p.
 *
 * @param p - Probability (0 < p < 1)
 * @returns x such that P(X ≤ x) = p
 *
 * @example
 * ```ts
 * ndtri(0.5);    // 0
 * ndtri(0.975);  // ≈ 1.96
 * ```
 */
export function ndtri(p: number): number;
export function ndtri(p: ArrayInput): Float64Array;
export function ndtri(p: number | ArrayInput): number | Float64Array {
  const xsf = ensureLoaded();

  if (typeof p === 'number') {
    return xsf._wasm_ndtri(p);
  }

  const arr = p instanceof Float64Array ? p : new Float64Array(p);
  const result = new Float64Array(arr.length);
  for (let i = 0; i < arr.length; i++) {
    result[i] = xsf._wasm_ndtri(arr[i]);
  }
  return result;
}

/**
 * Log of normal distribution CDF.
 *
 * More accurate than log(ndtr(x)) for very negative x.
 *
 * @param x - Input value
 * @returns log(Φ(x))
 */
export function log_ndtr(x: number): number;
export function log_ndtr(x: ArrayInput): Float64Array;
export function log_ndtr(x: number | ArrayInput): number | Float64Array {
  const xsf = ensureLoaded();

  if (typeof x === 'number') {
    return xsf._wasm_log_ndtr(x);
  }

  const arr = x instanceof Float64Array ? x : new Float64Array(x);
  const result = new Float64Array(arr.length);
  for (let i = 0; i < arr.length; i++) {
    result[i] = xsf._wasm_log_ndtr(arr[i]);
  }
  return result;
}

// =============================================================================
// CHI-SQUARE DISTRIBUTION
// =============================================================================

/**
 * Chi-square distribution CDF.
 *
 * @param df - Degrees of freedom
 * @param x - Input value (x ≥ 0)
 * @returns P(X ≤ x) for chi-square with df degrees of freedom
 */
export function chdtr(df: number, x: number): number;
export function chdtr(df: number, x: ArrayInput): Float64Array;
export function chdtr(
  df: number,
  x: number | ArrayInput
): number | Float64Array {
  const xsf = ensureLoaded();

  if (typeof x === 'number') {
    return xsf._wasm_chdtr(df, x);
  }

  const arr = x instanceof Float64Array ? x : new Float64Array(x);
  const result = new Float64Array(arr.length);
  for (let i = 0; i < arr.length; i++) {
    result[i] = xsf._wasm_chdtr(df, arr[i]);
  }
  return result;
}

/**
 * Chi-square distribution survival function (complement of CDF).
 *
 * @param df - Degrees of freedom
 * @param x - Input value (x ≥ 0)
 * @returns P(X > x)
 */
export function chdtrc(df: number, x: number): number;
export function chdtrc(df: number, x: ArrayInput): Float64Array;
export function chdtrc(
  df: number,
  x: number | ArrayInput
): number | Float64Array {
  const xsf = ensureLoaded();

  if (typeof x === 'number') {
    return xsf._wasm_chdtrc(df, x);
  }

  const arr = x instanceof Float64Array ? x : new Float64Array(x);
  const result = new Float64Array(arr.length);
  for (let i = 0; i < arr.length; i++) {
    result[i] = xsf._wasm_chdtrc(df, arr[i]);
  }
  return result;
}

/**
 * Chi-square distribution inverse CDF.
 *
 * @param df - Degrees of freedom
 * @param p - Probability (0 < p < 1)
 * @returns x such that P(X ≤ x) = p
 */
export function chdtri(df: number, p: number): number;
export function chdtri(df: number, p: ArrayInput): Float64Array;
export function chdtri(
  df: number,
  p: number | ArrayInput
): number | Float64Array {
  const xsf = ensureLoaded();

  if (typeof p === 'number') {
    return xsf._wasm_chdtri(df, p);
  }

  const arr = p instanceof Float64Array ? p : new Float64Array(p);
  const result = new Float64Array(arr.length);
  for (let i = 0; i < arr.length; i++) {
    result[i] = xsf._wasm_chdtri(df, arr[i]);
  }
  return result;
}

// =============================================================================
// F DISTRIBUTION
// =============================================================================

/**
 * F distribution CDF.
 *
 * @param dfn - Numerator degrees of freedom
 * @param dfd - Denominator degrees of freedom
 * @param x - Input value (x ≥ 0)
 * @returns P(X ≤ x) for F(dfn, dfd)
 */
export function fdtr(dfn: number, dfd: number, x: number): number;
export function fdtr(dfn: number, dfd: number, x: ArrayInput): Float64Array;
export function fdtr(
  dfn: number,
  dfd: number,
  x: number | ArrayInput
): number | Float64Array {
  const xsf = ensureLoaded();

  if (typeof x === 'number') {
    return xsf._wasm_fdtr(dfn, dfd, x);
  }

  const arr = x instanceof Float64Array ? x : new Float64Array(x);
  const result = new Float64Array(arr.length);
  for (let i = 0; i < arr.length; i++) {
    result[i] = xsf._wasm_fdtr(dfn, dfd, arr[i]);
  }
  return result;
}

/**
 * F distribution survival function.
 */
export function fdtrc(dfn: number, dfd: number, x: number): number;
export function fdtrc(dfn: number, dfd: number, x: ArrayInput): Float64Array;
export function fdtrc(
  dfn: number,
  dfd: number,
  x: number | ArrayInput
): number | Float64Array {
  const xsf = ensureLoaded();

  if (typeof x === 'number') {
    return xsf._wasm_fdtrc(dfn, dfd, x);
  }

  const arr = x instanceof Float64Array ? x : new Float64Array(x);
  const result = new Float64Array(arr.length);
  for (let i = 0; i < arr.length; i++) {
    result[i] = xsf._wasm_fdtrc(dfn, dfd, arr[i]);
  }
  return result;
}

/**
 * F distribution inverse CDF.
 */
export function fdtri(dfn: number, dfd: number, p: number): number;
export function fdtri(dfn: number, dfd: number, p: ArrayInput): Float64Array;
export function fdtri(
  dfn: number,
  dfd: number,
  p: number | ArrayInput
): number | Float64Array {
  const xsf = ensureLoaded();

  if (typeof p === 'number') {
    return xsf._wasm_fdtri(dfn, dfd, p);
  }

  const arr = p instanceof Float64Array ? p : new Float64Array(p);
  const result = new Float64Array(arr.length);
  for (let i = 0; i < arr.length; i++) {
    result[i] = xsf._wasm_fdtri(dfn, dfd, arr[i]);
  }
  return result;
}

// =============================================================================
// GAMMA DISTRIBUTION
// =============================================================================

/**
 * Gamma distribution CDF.
 *
 * @param a - Shape parameter
 * @param b - Scale parameter
 * @param x - Input value (x ≥ 0)
 * @returns P(X ≤ x)
 */
export function gdtr(a: number, b: number, x: number): number;
export function gdtr(a: number, b: number, x: ArrayInput): Float64Array;
export function gdtr(
  a: number,
  b: number,
  x: number | ArrayInput
): number | Float64Array {
  const xsf = ensureLoaded();

  if (typeof x === 'number') {
    return xsf._wasm_gdtr(a, b, x);
  }

  const arr = x instanceof Float64Array ? x : new Float64Array(x);
  const result = new Float64Array(arr.length);
  for (let i = 0; i < arr.length; i++) {
    result[i] = xsf._wasm_gdtr(a, b, arr[i]);
  }
  return result;
}

/**
 * Gamma distribution survival function.
 */
export function gdtrc(a: number, b: number, x: number): number;
export function gdtrc(a: number, b: number, x: ArrayInput): Float64Array;
export function gdtrc(
  a: number,
  b: number,
  x: number | ArrayInput
): number | Float64Array {
  const xsf = ensureLoaded();

  if (typeof x === 'number') {
    return xsf._wasm_gdtrc(a, b, x);
  }

  const arr = x instanceof Float64Array ? x : new Float64Array(x);
  const result = new Float64Array(arr.length);
  for (let i = 0; i < arr.length; i++) {
    result[i] = xsf._wasm_gdtrc(a, b, arr[i]);
  }
  return result;
}

// =============================================================================
// POISSON DISTRIBUTION
// =============================================================================

/**
 * Poisson distribution CDF.
 *
 * @param k - Number of events
 * @param m - Expected number of events (mean)
 * @returns P(X ≤ k)
 */
export function pdtr(k: number, m: number): number;
export function pdtr(k: ArrayInput, m: number): Float64Array;
export function pdtr(
  k: number | ArrayInput,
  m: number
): number | Float64Array {
  const xsf = ensureLoaded();

  if (typeof k === 'number') {
    return xsf._wasm_pdtr(k, m);
  }

  const arr = k instanceof Float64Array ? k : new Float64Array(k);
  const result = new Float64Array(arr.length);
  for (let i = 0; i < arr.length; i++) {
    result[i] = xsf._wasm_pdtr(arr[i], m);
  }
  return result;
}

/**
 * Poisson distribution survival function.
 */
export function pdtrc(k: number, m: number): number;
export function pdtrc(k: ArrayInput, m: number): Float64Array;
export function pdtrc(
  k: number | ArrayInput,
  m: number
): number | Float64Array {
  const xsf = ensureLoaded();

  if (typeof k === 'number') {
    return xsf._wasm_pdtrc(k, m);
  }

  const arr = k instanceof Float64Array ? k : new Float64Array(k);
  const result = new Float64Array(arr.length);
  for (let i = 0; i < arr.length; i++) {
    result[i] = xsf._wasm_pdtrc(arr[i], m);
  }
  return result;
}

/**
 * Poisson distribution inverse CDF.
 */
export function pdtri(k: number, p: number): number {
  const xsf = ensureLoaded();
  return xsf._wasm_pdtri(k, p);
}

// =============================================================================
// BINOMIAL DISTRIBUTION
// =============================================================================

/**
 * Binomial distribution CDF.
 *
 * @param k - Number of successes
 * @param n - Number of trials
 * @param p - Probability of success
 * @returns P(X ≤ k)
 */
export function bdtr(k: number, n: number, p: number): number {
  const xsf = ensureLoaded();
  return xsf._wasm_bdtr(k, n, p);
}

/**
 * Binomial distribution survival function.
 */
export function bdtrc(k: number, n: number, p: number): number {
  const xsf = ensureLoaded();
  return xsf._wasm_bdtrc(k, n, p);
}

/**
 * Binomial distribution inverse CDF.
 */
export function bdtri(k: number, n: number, y: number): number {
  const xsf = ensureLoaded();
  return xsf._wasm_bdtri(k, n, y);
}

// =============================================================================
// NEGATIVE BINOMIAL DISTRIBUTION
// =============================================================================

/**
 * Negative binomial distribution CDF.
 *
 * @param k - Number of failures
 * @param n - Number of successes
 * @param p - Probability of success
 * @returns P(X ≤ k)
 */
export function nbdtr(k: number, n: number, p: number): number {
  const xsf = ensureLoaded();
  return xsf._wasm_nbdtr(k, n, p);
}

/**
 * Negative binomial distribution survival function.
 */
export function nbdtrc(k: number, n: number, p: number): number {
  const xsf = ensureLoaded();
  return xsf._wasm_nbdtrc(k, n, p);
}

/**
 * Negative binomial distribution inverse CDF.
 */
export function nbdtri(k: number, n: number, p: number): number {
  const xsf = ensureLoaded();
  return xsf._wasm_nbdtri(k, n, p);
}

// =============================================================================
// KOLMOGOROV-SMIRNOV DISTRIBUTION
// =============================================================================

/**
 * Kolmogorov distribution CDF.
 *
 * @param x - Input value
 * @returns P(D ≤ x) where D is the Kolmogorov statistic
 */
export function kolmogorov(x: number): number;
export function kolmogorov(x: ArrayInput): Float64Array;
export function kolmogorov(x: number | ArrayInput): number | Float64Array {
  const xsf = ensureLoaded();

  if (typeof x === 'number') {
    return xsf._wasm_kolmogorov(x);
  }

  const arr = x instanceof Float64Array ? x : new Float64Array(x);
  const result = new Float64Array(arr.length);
  for (let i = 0; i < arr.length; i++) {
    result[i] = xsf._wasm_kolmogorov(arr[i]);
  }
  return result;
}

/**
 * Kolmogorov distribution inverse CDF.
 */
export function kolmogi(p: number): number;
export function kolmogi(p: ArrayInput): Float64Array;
export function kolmogi(p: number | ArrayInput): number | Float64Array {
  const xsf = ensureLoaded();

  if (typeof p === 'number') {
    return xsf._wasm_kolmogi(p);
  }

  const arr = p instanceof Float64Array ? p : new Float64Array(p);
  const result = new Float64Array(arr.length);
  for (let i = 0; i < arr.length; i++) {
    result[i] = xsf._wasm_kolmogi(arr[i]);
  }
  return result;
}

/**
 * Smirnov distribution CDF.
 *
 * @param n - Sample size
 * @param x - Input value
 * @returns P(D_n ≤ x)
 */
export function smirnov(n: number, x: number): number;
export function smirnov(n: number, x: ArrayInput): Float64Array;
export function smirnov(
  n: number,
  x: number | ArrayInput
): number | Float64Array {
  const xsf = ensureLoaded();

  if (typeof x === 'number') {
    return xsf._wasm_smirnov(n, x);
  }

  const arr = x instanceof Float64Array ? x : new Float64Array(x);
  const result = new Float64Array(arr.length);
  for (let i = 0; i < arr.length; i++) {
    result[i] = xsf._wasm_smirnov(n, arr[i]);
  }
  return result;
}

/**
 * Smirnov distribution inverse CDF.
 */
export function smirnovi(n: number, p: number): number;
export function smirnovi(n: number, p: ArrayInput): Float64Array;
export function smirnovi(
  n: number,
  p: number | ArrayInput
): number | Float64Array {
  const xsf = ensureLoaded();

  if (typeof p === 'number') {
    return xsf._wasm_smirnovi(n, p);
  }

  const arr = p instanceof Float64Array ? p : new Float64Array(p);
  const result = new Float64Array(arr.length);
  for (let i = 0; i < arr.length; i++) {
    result[i] = xsf._wasm_smirnovi(n, arr[i]);
  }
  return result;
}

// =============================================================================
// OWEN'S T FUNCTION
// =============================================================================

/**
 * Owen's T function.
 *
 * T(h, a) = (1/2π) ∫₀^a exp(-h²(1+x²)/2) / (1+x²) dx
 *
 * Useful for computing bivariate normal probabilities.
 *
 * @param h - First parameter
 * @param a - Second parameter
 * @returns T(h, a)
 */
export function owens_t(h: number, a: number): number {
  const xsf = ensureLoaded();
  return xsf._wasm_owens_t(h, a);
}
