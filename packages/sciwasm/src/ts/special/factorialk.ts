/**
 * Multifactorial function: n!(k)
 * Mirrors scipy.special.factorialk
 */

import { gamma } from './gamma.js';
import { rangeProd, factorialArrayExact } from './_helpers.js';

export interface FactorialkOptions {
  /**
   * If true, calculate exactly using integer arithmetic.
   * If false (default), use gamma function approximation (faster, returns floats).
   */
  exact?: boolean;

  /**
   * How to handle negative values:
   * - 'zero': Return 0 for negative inputs (default)
   * - 'complex': Not yet implemented
   */
  extend?: 'zero';
}

/**
 * Multifactorial of n of order k, n(!!...!).
 *
 * This is the multifactorial of n skipping k values. For example:
 *
 *     factorialk(17, 4) = 17!!!! = 17 × 13 × 9 × 5 × 1 = 8721
 *
 * In particular, for any integer n:
 *
 *     factorialk(n, 1) = factorial(n)
 *     factorialk(n, 2) = factorial2(n)
 *
 * With `exact=false` (default), it can be approximated as (where r = n % k):
 *
 *     n!(k) = k^((n-r)/k) × Γ(n/k + 1) / Γ(r/k + 1) × max(r, 1)
 *
 * @param n - Input value(s). Can be a number or array of numbers.
 * @param k - Order of multifactorial (how many values to skip).
 * @param options - Optional settings for exact calculation and extend mode.
 * @returns Multifactorial of order k of n.
 *
 * @example
 * ```ts
 * import { special } from 'sciwasm';
 *
 * // Exact mode
 * const exact = await special.factorialk(17, 4, { exact: true });
 * console.log(exact); // 8721 (17 × 13 × 9 × 5 × 1)
 *
 * // Approximate mode
 * const approx = await special.factorialk(17, 4);
 * console.log(approx); // ~8721.0
 *
 * // Special cases
 * const fact = await special.factorialk(5, 1, { exact: true });
 * console.log(fact); // 120 (equivalent to factorial(5))
 *
 * const fact2 = await special.factorialk(7, 2, { exact: true });
 * console.log(fact2); // 105 (equivalent to factorial2(7))
 *
 * // Array input
 * const arr = await special.factorialk([5, 7, 9], 3, { exact: true });
 * console.log(arr); // [10, 28, 162]
 * ```
 */
export async function factorialk(
  n: number | number[],
  k: number,
  options?: FactorialkOptions
): Promise<number | number[] | bigint | (number | bigint)[]> {
  const { exact = false, extend = 'zero' } = options ?? {};

  // Validate k parameter
  if (k < 1 || !Number.isInteger(k)) {
    throw new ValueError('k must be a positive integer');
  }

  // Validate extend parameter
  if (extend !== 'zero') {
    throw new ValueError(`extend must be 'zero', got: ${extend}`);
  }

  // Scalar case
  if (typeof n === 'number') {
    return factorialkScalar(n, k, exact, extend);
  }

  // Array case
  return factorialkArray(n, k, exact, extend);
}

/**
 * Compute multifactorial for a scalar value.
 */
async function factorialkScalar(
  n: number,
  k: number,
  exact: boolean,
  extend: string
): Promise<number | bigint> {
  // Handle special cases
  if (isNaN(n)) {
    return exact ? 0 : NaN;
  }

  if (extend === 'zero' && n < 0) {
    return exact ? 0 : 0.0;
  }

  if (n <= 1) {
    return exact ? 1 : 1.0;
  }

  // Exact computation
  if (exact) {
    // Validate integer input
    if (!Number.isInteger(n)) {
      throw new ValueError('exact=true only supports integer inputs');
    }
    // Use rangeProd with step k
    // For n=17, k=4: rangeProd(1, 17, 4) = 17*13*9*5*1
    return rangeProd(1, n, k);
  }

  // Approximate via gamma function
  // Formula: n!(k) = k^((n-r)/k) * Γ(n/k + 1) / Γ(r/k + 1) * max(r, 1)
  // where r = n % k
  const r = n % k;
  const nDivK = n / k;
  const rDivK = r / k;

  // Compute k^((n-r)/k)
  const powerK = Math.pow(k, (n - r) / k);

  // Compute Γ(n/k + 1) / Γ(r/k + 1)
  const gammaNumer = await gamma(nDivK + 1);
  const gammaDenom = await gamma(rDivK + 1);
  const gammaNumerVal = Array.isArray(gammaNumer) ? gammaNumer[0] : gammaNumer;
  const gammaDenomVal = Array.isArray(gammaDenom) ? gammaDenom[0] : gammaDenom;
  const gammaRatio = gammaNumerVal / gammaDenomVal;

  // Multiply by max(r, 1)
  const maxR = Math.max(r, 1);

  return powerK * gammaRatio * maxR;
}

/**
 * Compute multifactorial for an array of values.
 */
async function factorialkArray(
  n: number[],
  k: number,
  exact: boolean,
  extend: string
): Promise<number[] | (number | bigint)[]> {
  if (n.length === 0) {
    return [];
  }

  if (exact) {
    return factorialkArrayExactImpl(n, k, extend);
  }

  return factorialkArrayApproxImpl(n, k, extend);
}

/**
 * Compute exact multifactorials for an array.
 */
function factorialkArrayExactImpl(
  n: number[],
  k: number,
  extend: string
): (number | bigint)[] {
  // Validate all inputs are integers
  for (const val of n) {
    if (!Number.isInteger(val) && val >= 0) {
      throw new ValueError('exact=true only supports integer inputs');
    }
  }

  // Use optimized array factorial computation with step k
  return factorialArrayExact(n, k);
}

/**
 * Compute approximate multifactorials for an array using gamma function.
 */
async function factorialkArrayApproxImpl(
  n: number[],
  k: number,
  extend: string
): Promise<number[]> {
  const results: number[] = [];

  for (const val of n) {
    if (val < 0 && extend === 'zero') {
      results.push(0.0);
      continue;
    }

    if (isNaN(val)) {
      results.push(NaN);
      continue;
    }

    if (val <= 1) {
      results.push(1.0);
      continue;
    }

    // Formula: n!(k) = k^((n-r)/k) * Γ(n/k + 1) / Γ(r/k + 1) * max(r, 1)
    const r = val % k;
    const nDivK = val / k;
    const rDivK = r / k;

    const powerK = Math.pow(k, (val - r) / k);

    const gammaNumer = await gamma(nDivK + 1);
    const gammaDenom = await gamma(rDivK + 1);
    const gammaNumerVal = Array.isArray(gammaNumer) ? gammaNumer[0] : gammaNumer;
    const gammaDenomVal = Array.isArray(gammaDenom) ? gammaDenom[0] : gammaDenom;
    const gammaRatio = gammaNumerVal / gammaDenomVal;

    const maxR = Math.max(r, 1);

    results.push(powerK * gammaRatio * maxR);
  }

  return results;
}

/**
 * ValueError for invalid input parameters.
 */
class ValueError extends Error {
  constructor(message: string) {
    super(message);
    this.name = 'ValueError';
  }
}
