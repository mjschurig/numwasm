/**
 * Skewness computation.
 *
 * Computes the sample skewness of the data, a measure of asymmetry
 * of the probability distribution.
 */

import {
  NDArray,
  mean,
  multiply,
  divide,
  where,
  less_equal,
  full,
  finfo,
  absolute,
  power,
} from 'numwasm';
import { moment } from './moment.js';
import { scalarize } from './_utils.js';

/**
 * Compute the sample skewness of a data set.
 *
 * For normally distributed data, the skewness should be about zero.
 * For unimodal continuous distributions, a negative skew indicates
 * the tail on the left side is longer, while a positive skew
 * indicates the right side.
 *
 * Uses Fisher-Pearson coefficient of skewness: m3 / m2^(3/2).
 *
 * @param a - Input array
 * @param axis - Axis along which skewness is calculated. Default: 0
 * @param bias - If false, apply bias correction for statistical bias. Default: true
 * @returns Skewness values along the given axis
 */
export async function skew(
  a: NDArray,
  axis: number = 0,
  bias: boolean = true,
): Promise<NDArray | number> {
  const n = a.shape[axis];

  const m2 = await moment(a, 2, axis);
  const m3 = await moment(a, 3, axis);
  const mu = scalarize(await mean(a, axis));
  const eps = new finfo(a.dtype).eps;

  // Scalar path
  if (typeof m2 === 'number') {
    const muVal = mu as number;
    const threshold = (eps * Math.abs(muVal)) ** 2;
    if (m2 <= threshold) return NaN;
    const sk = (m3 as number) / m2 ** 1.5;
    if (!bias && n > 2) {
      return Math.sqrt(n * (n - 1)) / (n - 2) * sk;
    }
    return sk;
  }

  // Array path
  const m2Arr = m2 as NDArray;
  const m3Arr = m3 as NDArray;
  const muArr = mu as NDArray;

  // Build threshold: (eps * |mean|)^2
  const absMu = absolute(muArr);
  muArr.dispose();
  const epsScalar = await full(absMu.shape, eps);
  const epsMu = multiply(epsScalar, absMu);
  epsScalar.dispose();
  absMu.dispose();
  const threshold = multiply(epsMu, epsMu);
  epsMu.dispose();

  // Mask where m2 is effectively zero
  const zeroMask = less_equal(m2Arr, threshold);
  threshold.dispose();

  // m2^1.5
  const expArr = await full(m2Arr.shape, 1.5);
  const m2_15 = power(m2Arr, expArr);
  expArr.dispose();
  m2Arr.dispose();

  // Biased skewness: m3 / m2^1.5
  let sk = divide(m3Arr, m2_15);
  m3Arr.dispose();
  m2_15.dispose();

  // Apply bias correction if requested
  if (!bias && n > 2) {
    const correction = Math.sqrt(n * (n - 1)) / (n - 2);
    const corrArr = await full(sk.shape, correction);
    const corrected = multiply(corrArr, sk);
    corrArr.dispose();
    sk.dispose();
    sk = corrected;
  }

  // Replace near-zero-variance positions with NaN
  const nanArr = await full(sk.shape, NaN);
  const result = await where(zeroMask, nanArr, sk);
  nanArr.dispose();
  zeroMask.dispose();
  sk.dispose();

  return result;
}
