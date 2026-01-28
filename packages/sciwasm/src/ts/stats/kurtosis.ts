/**
 * Kurtosis computation.
 *
 * Computes the kurtosis (Fisher or Pearson) of the data.
 */

import {
  NDArray,
  mean,
  multiply,
  divide,
  subtract as sub,
  add,
  where,
  less_equal,
  full,
  finfo,
  absolute,
} from 'numwasm';
import { moment } from './moment.js';
import { scalarize } from './_utils.js';

/**
 * Compute the kurtosis (Fisher or Pearson) of a dataset.
 *
 * Kurtosis is the fourth central moment divided by the square of the
 * variance (second central moment). Fisher's definition subtracts 3
 * so that the normal distribution has kurtosis 0.
 *
 * @param a - Input array
 * @param axis - Axis along which kurtosis is calculated. Default: 0
 * @param fisher - If true (default), Fisher's definition (normal ==> 0.0).
 *                 If false, Pearson's definition (normal ==> 3.0).
 * @param bias - If false, apply bias correction. Default: true
 * @returns Kurtosis values along the given axis
 */
export async function kurtosis(
  a: NDArray,
  axis: number = 0,
  fisher: boolean = true,
  bias: boolean = true,
): Promise<NDArray | number> {
  const n = a.shape[axis];

  const m2 = await moment(a, 2, axis);
  const m4 = await moment(a, 4, axis);
  const mu = scalarize(await mean(a, axis));
  const eps = new finfo(a.dtype).eps;

  // Scalar path
  if (typeof m2 === 'number') {
    const muVal = mu as number;
    const threshold = (eps * Math.abs(muVal)) ** 2;
    if (m2 <= threshold) return NaN;

    let kurt = (m4 as number) / (m2 * m2);

    if (!bias && n > 3) {
      // Bias-corrected excess kurtosis:
      // (n^2-1)/((n-2)(n-3)) * (m4/m2^2) - 3*(n-1)^2/((n-2)(n-3))
      kurt = (1.0 / (n - 2) / (n - 3)) *
        ((n * n - 1.0) * (m4 as number) / (m2 * m2) - 3 * (n - 1) ** 2);
      if (!fisher) kurt += 3.0;
      return kurt;
    }

    return fisher ? kurt - 3.0 : kurt;
  }

  // Array path
  const m2Arr = m2 as NDArray;
  const m4Arr = m4 as NDArray;
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

  // m2^2
  const m2sq = multiply(m2Arr, m2Arr);

  // Raw kurtosis: m4 / m2^2
  let kurt = divide(m4Arr, m2sq);
  m4Arr.dispose();

  if (!bias && n > 3) {
    // Bias-corrected formula:
    // (1/(n-2)(n-3)) * ((n^2-1)*m4/m2^2 - 3*(n-1)^2)
    const nSq_1 = n * n - 1;
    const nSq1Arr = await full(kurt.shape, nSq_1);
    const term1 = multiply(nSq1Arr, kurt); // (n^2-1) * m4/m2^2
    nSq1Arr.dispose();
    kurt.dispose();

    const threeSq = 3 * (n - 1) ** 2;
    const threeSqArr = await full(term1.shape, threeSq);
    const diff = sub(term1, threeSqArr);
    term1.dispose();
    threeSqArr.dispose();

    const denom = (n - 2) * (n - 3);
    const denomArr = await full(diff.shape, denom);
    kurt = divide(diff, denomArr);
    diff.dispose();
    denomArr.dispose();

    if (!fisher) {
      const threeArr = await full(kurt.shape, 3.0);
      const adjusted = add(kurt, threeArr);
      threeArr.dispose();
      kurt.dispose();
      kurt = adjusted;
    }
  } else {
    // Fisher: subtract 3
    if (fisher) {
      const threeArr = await full(kurt.shape, 3.0);
      const adjusted = sub(kurt, threeArr);
      threeArr.dispose();
      kurt.dispose();
      kurt = adjusted;
    }
  }

  m2Arr.dispose();
  m2sq.dispose();

  // Replace near-zero-variance positions with NaN
  const nanArr = await full(kurt.shape, NaN);
  const result = await where(zeroMask, nanArr, kurt);
  nanArr.dispose();
  zeroMask.dispose();
  kurt.dispose();

  return result;
}
