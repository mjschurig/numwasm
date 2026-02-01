/**
 * Conjugate Dot Product (Pure TypeScript)
 *
 * Compute x^H * y for complex vectors
 *
 * Note: ZDOTC is for complex vectors. This implementation provides
 * a real-valued equivalent (same as regular dot product for real vectors).
 */

import { prepareVector } from '../helpers.js';
import type { Vector, DotcResult } from './types.js';

/**
 * Compute conjugate dot product x^H * y.
 *
 * For real vectors, this is equivalent to the regular dot product.
 * For complex vectors (stored as interleaved real/imag pairs),
 * this computes the conjugate dot product.
 *
 * @param x - First vector
 * @param y - Second vector (same length as x)
 * @returns Conjugate dot product result (real and imaginary parts)
 *
 * @example
 * ```typescript
 * // For real vectors, same as dot product
 * const x = [1, 2, 3];
 * const y = [4, 5, 6];
 * const { real } = dotc(x, y);
 * // real = 32, imag = 0
 *
 * // For complex vectors (interleaved format: [re1, im1, re2, im2, ...])
 * const cx = [1, 2, 3, 4]; // (1+2i), (3+4i)
 * const cy = [5, 6, 7, 8]; // (5+6i), (7+8i)
 * const { real: r, imag: i } = dotc(cx, cy);
 * // Computes conj(cx) . cy
 * ```
 */
export function dotc(x: Vector, y: Vector): DotcResult {
  const xData = prepareVector(x);
  const yData = prepareVector(y);

  // Check lengths match
  if (xData.length !== yData.length) {
    throw new Error(
      `Vector lengths must match: x has ${xData.length}, y has ${yData.length}`
    );
  }

  const n = xData.length;

  // Handle empty vectors
  if (n === 0) {
    return {
      real: 0,
      imag: 0,
      n: 0,
      success: true,
      message: 'Success',
    };
  }

  // Check if vectors are complex (even length for interleaved format)
  const isComplex = n % 2 === 0;

  if (isComplex) {
    // Treat as complex vectors in interleaved format
    // x = [re1, im1, re2, im2, ...], y = [re1, im1, re2, im2, ...]
    // dotc = sum of conj(x[k]) * y[k]
    // conj(a + bi) * (c + di) = (a - bi) * (c + di) = (ac + bd) + (ad - bc)i
    let realSum = 0;
    let imagSum = 0;

    for (let k = 0; k < n; k += 2) {
      const xRe = xData[k];
      const xIm = xData[k + 1];
      const yRe = yData[k];
      const yIm = yData[k + 1];

      // conj(x) * y = (xRe - i*xIm) * (yRe + i*yIm)
      // = xRe*yRe + xIm*yIm + i*(xRe*yIm - xIm*yRe)
      realSum += xRe * yRe + xIm * yIm;
      imagSum += xRe * yIm - xIm * yRe;
    }

    return {
      real: realSum,
      imag: imagSum,
      n: n / 2, // Number of complex elements
      success: true,
      message: 'Success',
    };
  } else {
    // Treat as real vectors (regular dot product)
    let sum = 0;
    for (let i = 0; i < n; i++) {
      sum += xData[i] * yData[i];
    }

    return {
      real: sum,
      imag: 0,
      n,
      success: true,
      message: 'Success',
    };
  }
}
