/**
 * Double factorial function: n!!
 * Mirrors scipy.special.factorial2
 */

import { gamma } from './gamma.js';
import { rangeProd, factorialArrayExact } from './_helpers.js';

export interface Factorial2Options {
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
 * Double factorial.
 *
 * This is the factorial with every second value skipped. For example:
 *
 *     7!! = 7 × 5 × 3 × 1 = 105
 *     8!! = 8 × 6 × 4 × 2 = 384
 *
 * It can be approximated numerically as:
 *
 *     n!! = 2^(n/2) × Γ(n/2 + 1) × sqrt(2/π)    if n is odd
 *     n!! = 2^(n/2) × Γ(n/2 + 1)                if n is even
 *
 * @param n - Input value(s). Can be a number or array of numbers.
 * @param options - Optional settings for exact calculation and extend mode.
 * @returns Double factorial of n.
 *
 * @example
 * ```ts
 * import { special } from 'sciwasm';
 *
 * // Double factorial of odd number
 * const odd = await special.factorial2(7, { exact: true });
 * console.log(odd); // 105 (7 × 5 × 3 × 1)
 *
 * // Double factorial of even number
 * const even = await special.factorial2(8, { exact: true });
 * console.log(even); // 384 (8 × 6 × 4 × 2)
 *
 * // Approximate mode
 * const approx = await special.factorial2(7);
 * console.log(approx); // ~105.0
 *
 * // Array input
 * const arr = await special.factorial2([3, 5, 7], { exact: true });
 * console.log(arr); // [3, 15, 105]
 * ```
 */
export async function factorial2(
  n: number | number[],
  options?: Factorial2Options
): Promise<number | number[] | bigint | (number | bigint)[]> {
  const { exact = false, extend = 'zero' } = options ?? {};

  // Validate extend parameter
  if (extend !== 'zero') {
    throw new ValueError(`extend must be 'zero', got: ${extend}`);
  }

  // Scalar case
  if (typeof n === 'number') {
    return factorial2Scalar(n, exact, extend);
  }

  // Array case
  return factorial2Array(n, exact, extend);
}

/**
 * Compute double factorial for a scalar value.
 */
async function factorial2Scalar(
  n: number,
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
    // Use rangeProd with k=2 (every 2nd number)
    // For n=7: rangeProd(1, 7, 2) = 7*5*3*1
    // For n=8: rangeProd(2, 8, 2) = 8*6*4*2
    const start = n % 2 === 0 ? 2 : 1;
    return rangeProd(start, n, 2);
  }

  // Approximate via gamma function
  const isOdd = n % 2 === 1;
  const halfN = n / 2;

  // Compute 2^(n/2)
  const power2 = Math.pow(2, halfN);

  // Compute Γ(n/2 + 1)
  const gammaValue = await gamma(halfN + 1);
  const gammaResult = Array.isArray(gammaValue) ? gammaValue[0] : gammaValue;

  if (isOdd) {
    // n!! = 2^(n/2) × Γ(n/2 + 1) × sqrt(2/π)
    const sqrt2OverPi = Math.sqrt(2 / Math.PI);
    return power2 * gammaResult * sqrt2OverPi;
  } else {
    // n!! = 2^(n/2) × Γ(n/2 + 1)
    return power2 * gammaResult;
  }
}

/**
 * Compute double factorial for an array of values.
 */
async function factorial2Array(
  n: number[],
  exact: boolean,
  extend: string
): Promise<number[] | (number | bigint)[]> {
  if (n.length === 0) {
    return [];
  }

  if (exact) {
    return factorial2ArrayExactImpl(n, extend);
  }

  return factorial2ArrayApproxImpl(n, extend);
}

/**
 * Compute exact double factorials for an array.
 */
function factorial2ArrayExactImpl(
  n: number[],
  _extend: string
): (number | bigint)[] {
  // Validate all inputs are integers
  for (const val of n) {
    if (!Number.isInteger(val) && val >= 0) {
      throw new ValueError('exact=true only supports integer inputs');
    }
  }

  // Use optimized array factorial computation with k=2
  return factorialArrayExact(n, 2);
}

/**
 * Compute approximate double factorials for an array using gamma function.
 */
async function factorial2ArrayApproxImpl(
  n: number[],
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

    const isOdd = val % 2 === 1;
    const halfVal = val / 2;
    const power2 = Math.pow(2, halfVal);

    const gammaValue = await gamma(halfVal + 1);
    const gammaResult = Array.isArray(gammaValue) ? gammaValue[0] : gammaValue;

    if (isOdd) {
      const sqrt2OverPi = Math.sqrt(2 / Math.PI);
      results.push(power2 * gammaResult * sqrt2OverPi);
    } else {
      results.push(power2 * gammaResult);
    }
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
