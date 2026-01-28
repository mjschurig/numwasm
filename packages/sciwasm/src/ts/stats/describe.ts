/**
 * Descriptive statistics.
 *
 * Computes several descriptive statistics of the passed array.
 */

import {
  NDArray,
  mean as nwMean,
  min as nwMin,
  max as nwMax,
  variance as nwVariance,
  isnan,
  any,
} from 'numwasm';
import { skew } from './skew.js';
import { kurtosis } from './kurtosis.js';
import { scalarize } from './_utils.js';

/**
 * Result of the describe function.
 */
export interface DescribeResult {
  /** Number of observations along the given axis */
  nobs: number;
  /** Tuple of (minimum, maximum) along the given axis */
  minmax: [NDArray | number, NDArray | number];
  /** Arithmetic mean along the given axis */
  mean: NDArray | number;
  /** Variance with the given ddof along the given axis */
  variance: NDArray | number;
  /** Skewness along the given axis */
  skewness: NDArray | number;
  /** Fisher's kurtosis along the given axis */
  kurtosis: NDArray | number;
}

export interface DescribeOptions {
  /** Axis along which statistics are calculated. Default: 0.
   *  If null, compute over the whole flattened array. */
  axis?: number | null;
  /** Delta degrees of freedom for variance. Default: 1 */
  ddof?: number;
  /** If false, correct skewness and kurtosis for statistical bias. Default: true */
  bias?: boolean;
  /** How to handle NaN values: 'propagate' or 'raise'. Default: 'propagate' */
  nanPolicy?: 'propagate' | 'raise';
}

/**
 * Compute several descriptive statistics of the passed array.
 *
 * Returns a DescribeResult with: nobs, minmax, mean, variance,
 * skewness, and kurtosis.
 *
 * @param a - Input data array
 * @param options - Configuration options
 * @returns DescribeResult object
 *
 * @throws {Error} If the input is empty
 * @throws {Error} If nan_policy is 'raise' and input contains NaN
 */
export async function describe(
  a: NDArray,
  options: DescribeOptions = {},
): Promise<DescribeResult> {
  const {
    axis = 0,
    ddof = 1,
    bias = true,
    nanPolicy = 'propagate',
  } = options;

  // Flatten if axis is null
  let data: NDArray;
  let effectiveAxis: number;
  if (axis === null) {
    data = a.flatten();
    effectiveAxis = 0;
  } else {
    data = a;
    effectiveAxis = axis;
  }

  // Check for empty input
  if (data.size === 0) {
    if (data !== a) data.dispose();
    throw new Error('The input must not be empty.');
  }

  // NaN policy
  if (nanPolicy === 'raise') {
    const nanMask = await isnan(data);
    const hasNan = await any(nanMask) as boolean;
    nanMask.dispose();
    if (hasNan) {
      if (data !== a) data.dispose();
      throw new Error('The input contains nan values');
    }
  } else if (nanPolicy !== 'propagate') {
    if (data !== a) data.dispose();
    throw new Error(
      "nan_policy must be one of {'propagate', 'raise'}",
    );
  }

  const n = data.shape[effectiveAxis];

  // Check for NaN with propagate policy — numwasm min/max skip NaN,
  // but scipy propagates NaN, so we need to detect and override
  let containsNan = false;
  if (nanPolicy === 'propagate') {
    const nanMask = await isnan(data);
    containsNan = await any(nanMask) as boolean;
    nanMask.dispose();
  }

  if (containsNan) {
    if (data !== a) data.dispose();
    return {
      nobs: n,
      minmax: [NaN, NaN],
      mean: NaN,
      variance: NaN,
      skewness: NaN,
      kurtosis: NaN,
    };
  }

  // Compute all statistics
  const [minVal, maxVal, meanVal, varianceVal, skewnessVal, kurtosisVal] =
    await Promise.all([
      nwMin(data, effectiveAxis),
      nwMax(data, effectiveAxis),
      nwMean(data, effectiveAxis),
      nwVariance(data, effectiveAxis, ddof),
      skew(data, effectiveAxis, bias),
      kurtosis(data, effectiveAxis, true, bias),
    ]);

  if (data !== a) data.dispose();

  return {
    nobs: n,
    minmax: [scalarize(minVal), scalarize(maxVal)],
    mean: scalarize(meanVal),
    variance: scalarize(varianceVal),
    skewness: scalarize(skewnessVal),
    kurtosis: scalarize(kurtosisVal),
  };
}
