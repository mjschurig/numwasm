/**
 * Factorial function: n!
 * Mirrors scipy.special.factorial
 */

import { gamma } from './gamma.js';
import { factorial64, factorialArrayExact } from './_helpers.js';

export interface FactorialOptions {
  /**
   * If true, calculate the factorial exactly using integer arithmetic.
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
 * The factorial of a number or array of numbers.
 *
 * The factorial of non-negative integer n is the product of all
 * positive integers less than or equal to n:
 *
 *     n! = n × (n-1) × (n-2) × ... × 2 × 1
 *
 * For example: 5! = 5 × 4 × 3 × 2 × 1 = 120
 *
 * With `exact=false` (default), the factorial is approximated using the
 * gamma function:
 *
 *     n! = Γ(n+1)
 *
 * @param n - Input value(s) for n!. Can be a number or array of numbers.
 * @param options - Optional settings for exact calculation and extend mode.
 * @returns Factorial of n. Type depends on options:
 *   - exact=false: returns number or number[]
 *   - exact=true: returns number, bigint, or array thereof
 *
 * @example
 * ```ts
 * import { special } from 'sciwasm';
 *
 * // Approximate mode (default)
 * const approx = await special.factorial(5);
 * console.log(approx); // 120.0
 *
 * // Exact integer mode
 * const exact = await special.factorial(5, { exact: true });
 * console.log(exact); // 120 (integer)
 *
 * // Large factorial with exact mode (returns bigint)
 * const large = await special.factorial(25, { exact: true });
 * console.log(large); // 15511210043330985984000000n
 *
 * // Array input
 * const arr = await special.factorial([3, 4, 5], { exact: false });
 * console.log(arr); // [6.0, 24.0, 120.0]
 *
 * // Negative values return 0 by default
 * const neg = await special.factorial(-1);
 * console.log(neg); // 0
 * ```
 */
export async function factorial(
  n: number | number[],
  options?: FactorialOptions
): Promise<number | number[] | bigint | (number | bigint)[]> {
  const { exact = false, extend = 'zero' } = options ?? {};

  // Validate extend parameter
  if (extend !== 'zero') {
    throw new ValueError(`extend must be 'zero', got: ${extend}`);
  }

  // Scalar case
  if (typeof n === 'number') {
    return factorialScalar(n, exact, extend);
  }

  // Array case
  return factorialArray(n, exact, extend);
}

/**
 * Compute factorial for a scalar value.
 */
async function factorialScalar(
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

  if (n === 0 || n === 1) {
    return exact ? 1 : 1.0;
  }

  // Exact computation
  if (exact) {
    // Validate integer input
    if (!Number.isInteger(n)) {
      throw new ValueError('exact=true only supports integer inputs');
    }
    return factorial64(n);
  }

  // Approximate via gamma function: n! = Γ(n+1)
  const result = await gamma(n + 1);
  return Array.isArray(result) ? result[0] : result;
}

/**
 * Compute factorial for an array of values.
 */
async function factorialArray(
  n: number[],
  exact: boolean,
  extend: string
): Promise<number[] | (number | bigint)[]> {
  if (n.length === 0) {
    return [];
  }

  if (exact) {
    return factorialArrayExactImpl(n, extend);
  }

  return factorialArrayApproxImpl(n, extend);
}

/**
 * Compute exact factorials for an array.
 */
function factorialArrayExactImpl(
  n: number[],
  _extend: string
): (number | bigint)[] {
  // Validate all inputs are integers
  for (const val of n) {
    if (!Number.isInteger(val) && val >= 0) {
      throw new ValueError('exact=true only supports integer inputs');
    }
  }

  // Use optimized array factorial computation
  return factorialArrayExact(n, 1);
}

/**
 * Compute approximate factorials for an array using gamma function.
 */
async function factorialArrayApproxImpl(
  n: number[],
  extend: string
): Promise<number[]> {
  // Handle extend='zero': negative values should return 0
  const adjusted = n.map(val => {
    if (val < 0) {
      return NaN; // Will be replaced with 0 after gamma computation
    }
    return val + 1; // n! = Γ(n+1)
  });

  // Compute gamma for all values
  const gammaResults = (await gamma(adjusted)) as number[];

  // Post-process: replace NaN with 0 for negative inputs in 'zero' mode
  if (extend === 'zero') {
    return gammaResults.map((val, i) => (n[i] < 0 ? 0.0 : val));
  }

  return gammaResults;
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
