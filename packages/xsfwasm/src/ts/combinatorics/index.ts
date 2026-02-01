/**
 * Combinatorial functions
 *
 * Includes binomial coefficients, Pochhammer symbol, and permutations.
 */

import { ensureLoaded, toFloat64Array, type ArrayInput } from '../core/utils.js';

/**
 * Binomial coefficient C(n, k) = n! / (k! * (n-k)!).
 *
 * Also known as "n choose k", counts the number of ways to choose
 * k items from n items without regard to order.
 *
 * Uses the gamma function generalization:
 *   C(n, k) = Γ(n+1) / (Γ(k+1) * Γ(n-k+1))
 *
 * Properties:
 * - C(n, 0) = C(n, n) = 1
 * - C(n, k) = C(n, n-k) (symmetry)
 * - C(n, k) = C(n-1, k-1) + C(n-1, k) (Pascal's rule)
 *
 * @param n - Total number of items (can be any real for generalized binomial)
 * @param k - Number of items to choose (can be any real)
 * @returns C(n, k)
 *
 * @example
 * ```ts
 * binom(5, 2);     // 10
 * binom(10, 5);    // 252
 * binom(-0.5, 2);  // 0.375 (generalized binomial)
 *
 * // Vectorized
 * binom([5, 6, 7], 2);  // Float64Array([10, 15, 21])
 * ```
 */
export function binom(n: number, k: number): number;
export function binom(n: ArrayInput, k: number): Float64Array;
export function binom(n: number, k: ArrayInput): Float64Array;
export function binom(n: ArrayInput, k: ArrayInput): Float64Array;
export function binom(
  n: number | ArrayInput,
  k: number | ArrayInput
): number | Float64Array {
  const xsf = ensureLoaded();

  const nIsScalar = typeof n === 'number';
  const kIsScalar = typeof k === 'number';

  if (nIsScalar && kIsScalar) {
    return xsf._wasm_binom(n as number, k as number);
  }

  const nArr = nIsScalar ? null : toFloat64Array(n as ArrayInput);
  const kArr = kIsScalar ? null : toFloat64Array(k as ArrayInput);

  const length = nArr?.length ?? kArr!.length;
  if (nArr && kArr && nArr.length !== kArr.length) {
    throw new Error(
      `Array length mismatch: got arrays of length ${nArr.length} and ${kArr.length}`
    );
  }

  const result = new Float64Array(length);
  for (let i = 0; i < length; i++) {
    const nVal = nIsScalar ? (n as number) : nArr![i];
    const kVal = kIsScalar ? (k as number) : kArr![i];
    result[i] = xsf._wasm_binom(nVal, kVal);
  }
  return result;
}

/**
 * Exact binomial coefficient for integer arguments.
 *
 * Uses integer arithmetic where possible to maintain exact results
 * for small to moderate values of n and k.
 *
 * @param n - Total number of items (non-negative integer)
 * @param k - Number of items to choose (non-negative integer, k ≤ n)
 * @returns C(n, k) computed with maximum precision
 *
 * @example
 * ```ts
 * binomExact(5, 2);    // 10
 * binomExact(20, 10);  // 184756
 * binomExact(5, 6);    // 0 (k > n)
 * ```
 */
export function binomExact(n: number, k: number): number;
export function binomExact(n: ArrayInput, k: number): Float64Array;
export function binomExact(n: number, k: ArrayInput): Float64Array;
export function binomExact(n: ArrayInput, k: ArrayInput): Float64Array;
export function binomExact(
  n: number | ArrayInput,
  k: number | ArrayInput
): number | Float64Array {
  const xsf = ensureLoaded();

  const nIsScalar = typeof n === 'number';
  const kIsScalar = typeof k === 'number';

  if (nIsScalar && kIsScalar) {
    return xsf._wasm_binom_exact(n as number, k as number);
  }

  const nArr = nIsScalar ? null : toFloat64Array(n as ArrayInput);
  const kArr = kIsScalar ? null : toFloat64Array(k as ArrayInput);

  const length = nArr?.length ?? kArr!.length;
  if (nArr && kArr && nArr.length !== kArr.length) {
    throw new Error(
      `Array length mismatch: got arrays of length ${nArr.length} and ${kArr.length}`
    );
  }

  const result = new Float64Array(length);
  for (let i = 0; i < length; i++) {
    const nVal = nIsScalar ? (n as number) : nArr![i];
    const kVal = kIsScalar ? (k as number) : kArr![i];
    result[i] = xsf._wasm_binom_exact(nVal, kVal);
  }
  return result;
}

/**
 * Pochhammer symbol (rising factorial) (x)_m.
 *
 * Definition:
 *   (x)_m = x * (x+1) * (x+2) * ... * (x+m-1) = Γ(x+m) / Γ(x)
 *
 * Special cases:
 * - (x)_0 = 1
 * - (x)_1 = x
 * - (1)_m = m! (factorial)
 *
 * @param x - Base value (any real number)
 * @param m - Number of terms (typically a non-negative integer)
 * @returns (x)_m
 *
 * @example
 * ```ts
 * poch(1, 5);    // 120 (= 5!)
 * poch(2, 3);    // 24 (= 2 * 3 * 4)
 * poch(0.5, 3);  // 1.875 (= 0.5 * 1.5 * 2.5)
 * ```
 */
export function poch(x: number, m: number): number;
export function poch(x: ArrayInput, m: number): Float64Array;
export function poch(x: number, m: ArrayInput): Float64Array;
export function poch(x: ArrayInput, m: ArrayInput): Float64Array;
export function poch(
  x: number | ArrayInput,
  m: number | ArrayInput
): number | Float64Array {
  const xsf = ensureLoaded();

  const xIsScalar = typeof x === 'number';
  const mIsScalar = typeof m === 'number';

  if (xIsScalar && mIsScalar) {
    return xsf._wasm_poch(x as number, m as number);
  }

  const xArr = xIsScalar ? null : toFloat64Array(x as ArrayInput);
  const mArr = mIsScalar ? null : toFloat64Array(m as ArrayInput);

  const length = xArr?.length ?? mArr!.length;
  if (xArr && mArr && xArr.length !== mArr.length) {
    throw new Error(
      `Array length mismatch: got arrays of length ${xArr.length} and ${mArr.length}`
    );
  }

  const result = new Float64Array(length);
  for (let i = 0; i < length; i++) {
    const xVal = xIsScalar ? (x as number) : xArr![i];
    const mVal = mIsScalar ? (m as number) : mArr![i];
    result[i] = xsf._wasm_poch(xVal, mVal);
  }
  return result;
}

/**
 * Exact number of permutations P(n, k) = n! / (n-k)!.
 *
 * Counts the number of ways to arrange k items from n items
 * where order matters.
 *
 * Properties:
 * - P(n, 0) = 1
 * - P(n, 1) = n
 * - P(n, n) = n!
 * - P(n, k) = n * P(n-1, k-1)
 *
 * @param n - Total number of items (non-negative integer)
 * @param k - Number of items to arrange (non-negative integer, k ≤ n)
 * @returns P(n, k) = n!/(n-k)!
 *
 * @example
 * ```ts
 * permExact(5, 2);    // 20 (= 5 * 4)
 * permExact(5, 5);    // 120 (= 5!)
 * permExact(10, 3);   // 720 (= 10 * 9 * 8)
 * ```
 */
export function permExact(n: number, k: number): number;
export function permExact(n: ArrayInput, k: number): Float64Array;
export function permExact(n: number, k: ArrayInput): Float64Array;
export function permExact(n: ArrayInput, k: ArrayInput): Float64Array;
export function permExact(
  n: number | ArrayInput,
  k: number | ArrayInput
): number | Float64Array {
  const xsf = ensureLoaded();

  const nIsScalar = typeof n === 'number';
  const kIsScalar = typeof k === 'number';

  if (nIsScalar && kIsScalar) {
    return xsf._wasm_perm_exact(n as number, k as number);
  }

  const nArr = nIsScalar ? null : toFloat64Array(n as ArrayInput);
  const kArr = kIsScalar ? null : toFloat64Array(k as ArrayInput);

  const length = nArr?.length ?? kArr!.length;
  if (nArr && kArr && nArr.length !== kArr.length) {
    throw new Error(
      `Array length mismatch: got arrays of length ${nArr.length} and ${kArr.length}`
    );
  }

  const result = new Float64Array(length);
  for (let i = 0; i < length; i++) {
    const nVal = nIsScalar ? (n as number) : nArr![i];
    const kVal = kIsScalar ? (k as number) : kArr![i];
    result[i] = xsf._wasm_perm_exact(nVal, kVal);
  }
  return result;
}
