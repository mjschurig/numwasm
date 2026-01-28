/**
 * Central moment computation.
 *
 * Computes the nth central moment of the data along the given axis:
 *   moment(a, n) = mean((a - mean(a))^n)
 */

import { NDArray, mean, subtract, multiply } from 'numwasm';
import { scalarize } from './_utils.js';

/**
 * Compute the nth central moment of the data along the given axis.
 *
 * The nth central moment is E[(x - mean(x))^n].
 *
 * @param a - Input array
 * @param order - Order of the central moment (must be >= 0)
 * @param axis - Axis along which to compute the moment. Default: 0
 * @returns The nth central moment as an NDArray or number
 */
export async function moment(
  a: NDArray,
  order: number = 1,
  axis: number = 0,
): Promise<NDArray | number> {
  if (order === 0) {
    return 1.0;
  }

  // Compute mean with keepdims so we can broadcast-subtract
  const mu = await mean(a, axis, true) as NDArray;
  const centered = subtract(a, mu);
  mu.dispose();

  // Exponentiation by squaring for efficiency (matches scipy)
  // Build sequence of exponents via recursive halving
  const nList: number[] = [order];
  let current = order;
  while (current > 2) {
    if (current % 2 !== 0) {
      current = (current - 1) / 2;
    } else {
      current = current / 2;
    }
    nList.push(current);
  }

  // Start from the smallest exponent
  let s: NDArray;
  const smallest = nList[nList.length - 1];
  if (smallest === 1) {
    s = centered;
  } else {
    // smallest === 2
    s = multiply(centered, centered);
  }

  // Build up via squaring
  for (let i = nList.length - 2; i >= 0; i--) {
    const prev = s;
    s = multiply(s, s);
    if (prev !== centered) prev.dispose();

    if (nList[i] % 2 !== 0) {
      const prev2 = s;
      s = multiply(s, centered);
      prev2.dispose();
    }
  }

  if (s !== centered) {
    centered.dispose();
  }

  // Take the mean of the powered deviations along the axis
  const result = await mean(s, axis);
  s.dispose();

  return scalarize(result);
}
