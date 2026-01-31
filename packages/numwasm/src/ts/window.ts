/**
 * Window functions for signal processing.
 *
 * Implements NumPy-compatible window functions:
 * - blackman(M) - Near-optimal window using three-term cosine sum
 * - bartlett(M) - Triangular window with zero endpoints
 * - hanning(M) - Cosine bell / Hann window
 * - hamming(M) - Similar to Hanning with raised endpoints
 * - kaiser(M, beta) - Bessel-based window with adjustable sidelobe level
 * - i0(x) - Modified Bessel function of the first kind, order 0
 *
 * Reference: numpy/lib/_function_base_impl.py lines 3037-3723
 */

import { NDArray } from "./_core/NDArray.js";
import { DType } from "./types.js";

// ============================================================================
// Bessel Function I₀ Implementation
// ============================================================================

/**
 * Chebyshev coefficients for I₀ computation (from Cephes library)
 * Used for |x| <= 8
 * @internal
 */
const _i0A: number[] = [
  -4.4153416464793393795e-18, 3.33079451882223809783e-17,
  -2.43127984654795469359e-16, 1.71539128555513303061e-15,
  -1.16853328779934516808e-14, 7.67618549860493561688e-14,
  -4.8564467831119294609e-13, 2.95505266312963983461e-12,
  -1.72682629144155570723e-11, 9.67580903537323691224e-11,
  -5.18979560163526290666e-10, 2.65982372468238665035e-9,
  -1.30002500998624804212e-8, 6.04699502254191894932e-8,
  -2.67079385394061173391e-7, 1.11738753912010371815e-6,
  -4.41673835845875056359e-6, 1.64484480707288970893e-5,
  -5.75419501008210370398e-5, 1.88502885095841655729e-4,
  -5.76375574538582365885e-4, 1.63947561694133579842e-3,
  -4.3243099950505759443e-3, 1.05464603945949983183e-2,
  -2.37374148058994688156e-2, 4.93052842396707084878e-2,
  -9.4901097048047644421e-2, 1.71620901522208775349e-1,
  -3.04682672343198398683e-1, 6.76795274409476084995e-1,
];

/**
 * Chebyshev coefficients for I₀ computation (from Cephes library)
 * Used for |x| > 8
 * @internal
 */
const _i0B: number[] = [
  -7.23318048787475395456e-18, -4.83050448594418207126e-18,
  4.46562142029675999901e-17, 3.4612228676974610931e-17,
  -2.82762398051658348494e-16, -3.42548561967721913462e-16,
  1.7725601330565263836e-15, 3.81168066935262242075e-15,
  -9.5548466988283076487e-15, -4.15056934728722208663e-14,
  1.54008621752140982691e-14, 3.85277838274214270114e-13,
  7.18012445138366623367e-13, -1.79417853150680611778e-12,
  -1.32158118404477131188e-11, -3.14991652796324136454e-11,
  1.18891471078464383424e-11, 4.9406023882249695891e-10,
  3.39623202570838634515e-9, 2.26666899049817806459e-8,
  2.04891858946906374183e-7, 2.89137052083475648297e-6,
  6.88975834691682398426e-5, 3.3691164782556940899e-3,
  8.04490411014108831608e-1,
];

/**
 * Evaluate Chebyshev polynomial using Clenshaw recurrence.
 * @internal
 */
function _chbevl(x: number, vals: number[]): number {
  let b0 = vals[0];
  let b1 = 0.0;
  let b2 = 0.0;

  for (let i = 1; i < vals.length; i++) {
    b2 = b1;
    b1 = b0;
    b0 = x * b1 - b2 + vals[i];
  }

  return 0.5 * (b0 - b2);
}

/**
 * Compute I₀(x) for |x| <= 8 using Chebyshev approximation.
 * @internal
 */
function _i0_1(x: number): number {
  return Math.exp(x) * _chbevl(x / 2.0 - 2, _i0A);
}

/**
 * Compute I₀(x) for |x| > 8 using asymptotic Chebyshev approximation.
 * @internal
 */
function _i0_2(x: number): number {
  return (Math.exp(x) * _chbevl(32.0 / x - 2.0, _i0B)) / Math.sqrt(x);
}

/**
 * Compute I₀(x) for a single scalar value.
 * @internal
 */
function _i0_scalar(x: number): number {
  const absX = Math.abs(x);
  return absX <= 8.0 ? _i0_1(absX) : _i0_2(absX);
}

/**
 * Modified Bessel function of the first kind, order 0.
 *
 * Usually denoted I₀. This implementation uses Chebyshev polynomial
 * approximations from the Cephes library, partitioning the domain
 * into [0, 8] and (8, ∞).
 *
 * @param x - Input array, number, or array of numbers
 * @returns Modified Bessel function I₀ evaluated at x
 *
 * @example
 * ```typescript
 * await i0(0)           // NDArray with value 1.0
 * await i0([0, 1, 2])   // NDArray([1.0, 1.266..., 2.279...])
 * ```
 *
 * @remarks
 * Complex values are not supported and will throw an error.
 * Relative error on [0, 30] is documented as having a peak of 5.8e-16.
 *
 * Reference: numpy/lib/_function_base_impl.py lines 3533-3590
 */
export async function i0(x: NDArray | number | number[]): Promise<NDArray> {
  // Handle scalar input
  if (typeof x === "number") {
    const result = await NDArray.fromArray([_i0_scalar(x)]);
    return result.reshape([]);
  }

  // Handle array input
  let arr: NDArray;
  if (Array.isArray(x)) {
    arr = await NDArray.fromArray(x);
  } else {
    arr = x;
  }

  // Validate dtype - complex not supported
  if (arr.dtype === DType.Complex64 || arr.dtype === DType.Complex128) {
    throw new TypeError("i0 not supported for complex values");
  }

  // Create output array
  const result = await NDArray.empty(arr.shape, { dtype: DType.Float64 });

  // Compute I₀ for each element
  const size = arr.size;
  for (let i = 0; i < size; i++) {
    const val = arr.getFlat(i);
    result.setFlat(i, _i0_scalar(val));
  }

  return result;
}

// ============================================================================
// Window Functions
// ============================================================================

/**
 * Return the Blackman window.
 *
 * The Blackman window is a taper formed by using the first three
 * terms of a summation of cosines. It was designed to have close to
 * the minimal leakage possible.
 *
 * Formula: w(n) = 0.42 + 0.5*cos(π*n/(M-1)) + 0.08*cos(2π*n/(M-1))
 *
 * @param M - Number of points in the output window. If zero or less,
 *            an empty array is returned.
 * @returns The window, with the maximum value normalized to one
 *          (appears only if M is odd)
 *
 * @example
 * ```typescript
 * await blackman(12)
 * // NDArray([-1.39e-17, 0.0326, 0.1599, 0.4144, 0.7360, 0.9670,
 * //          0.9670, 0.7360, 0.4144, 0.1599, 0.0326, -1.39e-17])
 * ```
 *
 * Reference: numpy/lib/_function_base_impl.py lines 3037-3133
 */
export async function blackman(M: number): Promise<NDArray> {
  M = Math.floor(M);

  if (M < 1) {
    return NDArray.fromArray([]);
  }
  if (M === 1) {
    return NDArray.fromArray([1.0]);
  }

  // n = arange(1 - M, M, 2) gives symmetric indices centered at 0
  // This creates values: [1-M, 3-M, 5-M, ..., M-3, M-1]
  const n = await NDArray.arange(1 - M, M, 2);
  const denom = M - 1;

  // w(n) = 0.42 + 0.5*cos(π*n/(M-1)) + 0.08*cos(2π*n/(M-1))
  const result = await NDArray.empty([M], { dtype: DType.Float64 });
  for (let i = 0; i < M; i++) {
    const ni = n.getFlat(i);
    const term1 = 0.42;
    const term2 = 0.5 * Math.cos((Math.PI * ni) / denom);
    const term3 = 0.08 * Math.cos((2.0 * Math.PI * ni) / denom);
    result.setFlat(i, term1 + term2 + term3);
  }

  n.dispose();
  return result;
}

/**
 * Return the Hanning window (also known as Hann window).
 *
 * The Hanning window is a taper formed by using a weighted cosine.
 * Named for Julius von Hann.
 *
 * Formula: w(n) = 0.5 + 0.5*cos(π*n/(M-1))
 *
 * @param M - Number of points in the output window. If zero or less,
 *            an empty array is returned.
 * @returns The window, with the maximum value normalized to one
 *          (appears only if M is odd)
 *
 * @example
 * ```typescript
 * await hanning(12)
 * // NDArray([0, 0.0794, 0.2923, 0.5712, 0.8274, 0.9797,
 * //          0.9797, 0.8274, 0.5712, 0.2923, 0.0794, 0])
 * ```
 *
 * Reference: numpy/lib/_function_base_impl.py lines 3244-3342
 */
export async function hanning(M: number): Promise<NDArray> {
  M = Math.floor(M);

  if (M < 1) {
    return NDArray.fromArray([]);
  }
  if (M === 1) {
    return NDArray.fromArray([1.0]);
  }

  const n = await NDArray.arange(1 - M, M, 2);
  const denom = M - 1;

  // w(n) = 0.5 + 0.5*cos(π*n/(M-1))
  const result = await NDArray.empty([M], { dtype: DType.Float64 });
  for (let i = 0; i < M; i++) {
    const ni = n.getFlat(i);
    result.setFlat(i, 0.5 + 0.5 * Math.cos((Math.PI * ni) / denom));
  }

  n.dispose();
  return result;
}

/**
 * Return the Hamming window.
 *
 * The Hamming window is a taper formed by using a weighted cosine,
 * similar to Hanning but with raised endpoints to reduce the first
 * side lobe.
 *
 * Formula: w(n) = 0.54 + 0.46*cos(π*n/(M-1))
 *
 * @param M - Number of points in the output window. If zero or less,
 *            an empty array is returned.
 * @returns The window, with the maximum value normalized to one
 *          (appears only if M is odd)
 *
 * @example
 * ```typescript
 * await hamming(12)
 * // NDArray([0.08, 0.153, 0.349, 0.605, 0.841, 0.981,
 * //          0.981, 0.841, 0.605, 0.349, 0.153, 0.08])
 * ```
 *
 * Reference: numpy/lib/_function_base_impl.py lines 3346-3441
 */
export async function hamming(M: number): Promise<NDArray> {
  M = Math.floor(M);

  if (M < 1) {
    return NDArray.fromArray([]);
  }
  if (M === 1) {
    return NDArray.fromArray([1.0]);
  }

  const n = await NDArray.arange(1 - M, M, 2);
  const denom = M - 1;

  // w(n) = 0.54 + 0.46*cos(π*n/(M-1))
  const result = await NDArray.empty([M], { dtype: DType.Float64 });
  for (let i = 0; i < M; i++) {
    const ni = n.getFlat(i);
    result.setFlat(i, 0.54 + 0.46 * Math.cos((Math.PI * ni) / denom));
  }

  n.dispose();
  return result;
}

/**
 * Return the Bartlett window.
 *
 * The Bartlett window is very similar to a triangular window, except
 * that the end points are at zero. Also known as the triangular or
 * Fejér window.
 *
 * Formula:
 *   w(n) = 1 + n/(M-1) for n <= 0
 *   w(n) = 1 - n/(M-1) for n > 0
 *
 * @param M - Number of points in the output window. If zero or less,
 *            an empty array is returned.
 * @returns The triangular window, with the maximum value normalized to one
 *          (appears only if M is odd), with first and last samples equal
 *          to zero.
 *
 * @example
 * ```typescript
 * await bartlett(12)
 * // NDArray([0, 0.182, 0.364, 0.545, 0.727, 0.909,
 * //          0.909, 0.727, 0.545, 0.364, 0.182, 0])
 * ```
 *
 * Reference: numpy/lib/_function_base_impl.py lines 3137-3240
 */
export async function bartlett(M: number): Promise<NDArray> {
  M = Math.floor(M);

  if (M < 1) {
    return NDArray.fromArray([]);
  }
  if (M === 1) {
    return NDArray.fromArray([1.0]);
  }

  const n = await NDArray.arange(1 - M, M, 2);
  const denom = M - 1;

  // w(n) = 1 + n/(M-1) for n <= 0
  // w(n) = 1 - n/(M-1) for n > 0
  const result = await NDArray.empty([M], { dtype: DType.Float64 });
  for (let i = 0; i < M; i++) {
    const ni = n.getFlat(i);
    if (ni <= 0) {
      result.setFlat(i, 1 + ni / denom);
    } else {
      result.setFlat(i, 1 - ni / denom);
    }
  }

  n.dispose();
  return result;
}

/**
 * Return the Kaiser window.
 *
 * The Kaiser window is a taper formed by using a Bessel function.
 * It provides a good approximation to the Digital Prolate Spheroidal
 * Sequence (DPSS/Slepian window).
 *
 * Formula: w(n) = I₀(β * sqrt(1 - (2n/(M-1))²)) / I₀(β)
 *
 * The beta parameter controls the window shape:
 * - β = 0: Rectangular window
 * - β ≈ 5: Similar to Hamming
 * - β ≈ 6: Similar to Hanning
 * - β ≈ 8.6: Similar to Blackman
 * - β = 14: Good starting point for most applications
 *
 * @param M - Number of points in the output window. If zero or less,
 *            an empty array is returned.
 * @param beta - Shape parameter for the window. Larger values produce
 *               narrower main lobe with lower side lobes.
 * @returns The window, with the maximum value normalized to one
 *          (appears only if M is odd)
 *
 * @example
 * ```typescript
 * await kaiser(12, 14)
 * // NDArray([7.73e-06, 0.00346, 0.0465, 0.230, 0.600, 0.946,
 * //          0.946, 0.600, 0.230, 0.0465, 0.00346, 7.73e-06])
 * ```
 *
 * @remarks
 * As beta increases, the window narrows. If beta is too large relative
 * to M, NaN values may be returned.
 *
 * Reference: numpy/lib/_function_base_impl.py lines 3596-3723
 */
export async function kaiser(M: number, beta: number): Promise<NDArray> {
  M = Math.floor(M);

  if (M < 1) {
    return NDArray.fromArray([]);
  }
  if (M === 1) {
    return NDArray.fromArray([1.0]);
  }

  const alpha = (M - 1) / 2.0;
  const i0Beta = _i0_scalar(beta);

  // w(n) = I₀(β * sqrt(1 - ((n - α) / α)²)) / I₀(β)
  const result = await NDArray.empty([M], { dtype: DType.Float64 });
  for (let n = 0; n < M; n++) {
    const x = (n - alpha) / alpha;
    const arg = 1 - x * x;
    // Handle numerical issues where arg might be slightly negative
    const sqrtArg = arg >= 0 ? Math.sqrt(arg) : 0;
    result.setFlat(n, _i0_scalar(beta * sqrtArg) / i0Beta);
  }

  return result;
}
