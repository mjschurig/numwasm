/**
 * Internal helper functions for special function calculations.
 * Based on scipy.special._basic.py
 */

/**
 * Factorial lookup table for n=0 to n=20 (fits in float64 precisely).
 * Values beyond 20! may lose precision in JavaScript's number type.
 */
const FACTORIAL_TABLE: number[] = [
  1, // 0!
  1, // 1!
  2, // 2!
  6, // 3!
  24, // 4!
  120, // 5!
  720, // 6!
  5040, // 7!
  40320, // 8!
  362880, // 9!
  3628800, // 10!
  39916800, // 11!
  479001600, // 12!
  6227020800, // 13!
  87178291200, // 14!
  1307674368000, // 15!
  20922789888000, // 16!
  355687428096000, // 17!
  6402373705728000, // 18!
  121645100408832000, // 19!
  2432902008176640000, // 20!
];

/**
 * Compute factorial using BigInt for exact integer arithmetic.
 * For n <= 20, uses lookup table. For n > 20, computes using BigInt.
 */
export function factorial64(n: number): number | bigint {
  if (n < 0) {
    return 0;
  }
  if (n <= 20) {
    return FACTORIAL_TABLE[n];
  }
  // For n > 20, use BigInt to avoid overflow
  return factorialBigInt(BigInt(n));
}

/**
 * Compute factorial using BigInt arithmetic.
 */
function factorialBigInt(n: bigint): bigint {
  if (n <= 1n) {
    return 1n;
  }
  let result = 1n;
  for (let i = 2n; i <= n; i++) {
    result *= i;
  }
  return result;
}

/**
 * Product of a range of numbers spaced k apart (from hi).
 *
 * For k=1, this returns the product of:
 *   lo * (lo+1) * (lo+2) * ... * (hi-2) * (hi-1) * hi
 *   = hi! / (lo-1)!
 *
 * For k>1, it corresponds to taking only every k'th number when
 * counting down from hi - e.g. 18!!!! = rangeProd(1, 18, 4).
 *
 * Uses divide-and-conquer to break into smaller products for better performance.
 *
 * Based on scipy.special._basic._range_prod
 */
export function rangeProd(lo: number, hi: number, k: number = 1): number | bigint {
  // Fast path: factorial via lookup table or built-in
  if (lo === 1 && k === 1) {
    return factorial64(hi);
  }

  // Base cases
  if (lo + k < hi) {
    // Divide and conquer
    let mid = Math.floor((hi + lo) / 2);
    if (k > 1) {
      // Ensure mid is a multiple of k away from hi
      mid = mid - ((mid - hi) % k);
    }
    const left = rangeProd(lo, mid, k);
    const right = rangeProd(mid + k, hi, k);
    return multiplyBig(left, right);
  } else if (lo + k === hi) {
    // Two elements: lo and hi
    return multiplyBig(lo, hi);
  } else {
    // Single element
    return BigInt(hi);
  }
}

/**
 * Multiply two values that may be number or bigint.
 * Returns bigint if either input is bigint, otherwise number.
 */
function multiplyBig(a: number | bigint, b: number | bigint): number | bigint {
  // If either is bigint, convert both to bigint
  if (typeof a === 'bigint' || typeof b === 'bigint') {
    return BigInt(a) * BigInt(b);
  }
  // Both are numbers
  const product = a * b;
  // If product is too large for safe integer, convert to bigint
  if (!Number.isSafeInteger(product)) {
    return BigInt(a) * BigInt(b);
  }
  return product;
}

/**
 * Compute factorial for array of values using incremental approach.
 * This is more efficient than computing each factorial independently.
 *
 * Based on scipy.special._basic._factorialx_array_exact
 */
export function factorialArrayExact(n: number[], k: number = 1): (number | bigint)[] {
  if (n.length === 0) {
    return [];
  }

  // Get unique values and sort them
  const unique = Array.from(new Set(n)).sort((a, b) => a - b);

  // Map to store computed factorials
  const factorials = new Map<number, number | bigint>();

  // Handle negative and trivial values
  for (const val of unique) {
    if (val < 0) {
      factorials.set(val, 0);
    } else if (val < 2) {
      factorials.set(val, 1);
    }
  }

  // Filter to values > 1
  const validUnique = unique.filter(val => val > 1);

  if (validUnique.length === 0) {
    // Only trivial values
    return n.map(val => factorials.get(val)!);
  }

  // Process by lanes (residues modulo k)
  for (let lane = 0; lane < k; lane++) {
    // Get values in this lane
    const laneValues = validUnique.filter(val => (val % k) === lane);

    if (laneValues.length === 0) {
      continue;
    }

    // Compute factorials incrementally
    let val = rangeProd(1, laneValues[0], k);
    factorials.set(laneValues[0], val);

    for (let i = 0; i < laneValues.length - 1; i++) {
      const prev = laneValues[i];
      const current = laneValues[i + 1];
      // Multiply by the product of numbers between prev and current
      const increment = rangeProd(prev + k, current, k);
      val = multiplyBig(val, increment);
      factorials.set(current, val);
    }
  }

  // Map original input array to computed factorials
  return n.map(val => factorials.get(val)!);
}

/**
 * Check if a number is a safe integer (can be represented precisely as a number).
 */
function isSafeInteger(n: number): boolean {
  return Number.isSafeInteger(n);
}
