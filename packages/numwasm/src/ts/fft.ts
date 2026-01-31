/**
 * NumJS FFT (Fast Fourier Transform) Module
 *
 * Provides NumPy-compatible FFT operations:
 * - 1D transforms: fft, ifft, rfft, irfft, hfft, ihfft
 * - 2D transforms: fft2, ifft2, rfft2, irfft2
 * - N-D transforms: fftn, ifftn, rfftn, irfftn
 * - Helper functions: fftfreq, rfftfreq, fftshift, ifftshift
 */

import { NDArray } from "./_core/NDArray.js";
import { DType } from "./types.js";
import { getWasmModule } from "./wasm-loader.js";
import { roll } from "./manipulation.js";

/* ============ Types ============ */

/**
 * FFT normalization modes.
 * - "backward": No scaling on forward, 1/n on inverse (default)
 * - "ortho": 1/sqrt(n) on both forward and inverse
 * - "forward": 1/n on forward, no scaling on inverse
 */
export type FFTNorm = "backward" | "ortho" | "forward" | null;

/* ============ Internal Utilities ============ */

/**
 * Get normalization factor for FFT.
 */
function _getNormFactor(n: number, norm: FFTNorm, inverse: boolean): number {
  const effectiveNorm = norm ?? "backward";

  switch (effectiveNorm) {
    case "backward":
      return inverse ? 1 : 1; // Inverse scaling is done in WASM
    case "ortho":
      return 1 / Math.sqrt(n);
    case "forward":
      return inverse ? n : 1 / n; // Undo WASM scaling for inverse, apply for forward
    default:
      throw new Error(`Unknown normalization mode: ${norm}`);
  }
}

/**
 * Check if normalization requires post-processing.
 * WASM always uses "backward" convention (no scaling forward, 1/n inverse).
 */
function _needsNormAdjustment(norm: FFTNorm): boolean {
  const effectiveNorm = norm ?? "backward";
  if (effectiveNorm === "backward") return false;
  return true;
}

/**
 * Apply scalar multiplication to an array in-place.
 * Returns the same array (modified).
 */
function _scaleArrayInPlace(arr: NDArray, factor: number): void {
  const size = arr.size;
  const dtype = arr.dtype;

  if (dtype === DType.Complex64 || dtype === DType.Complex128) {
    // Scale complex array
    // Note: setComplex signature is (real, imag, ...indices)
    // For 1D array with flat indexing, we need to compute multi-dim indices
    const shape = arr.shape;

    for (let flatIdx = 0; flatIdx < size; flatIdx++) {
      const c = arr.getComplex(..._flatToMultiIndex(flatIdx, shape));
      const indices = _flatToMultiIndex(flatIdx, shape);
      arr.setComplex(c.real * factor, c.imag * factor, ...indices);
    }
  } else {
    // Scale real array
    for (let i = 0; i < size; i++) {
      arr.setFlat(i, arr.getFlat(i) * factor);
    }
  }
}

/**
 * Convert flat index to multi-dimensional indices.
 */
function _flatToMultiIndex(flatIdx: number, shape: number[]): number[] {
  const indices: number[] = new Array(shape.length);
  let remaining = flatIdx;

  for (let i = shape.length - 1; i >= 0; i--) {
    indices[i] = remaining % shape[i];
    remaining = Math.floor(remaining / shape[i]);
  }

  return indices;
}

/**
 * Compute complex conjugate of an array element-wise.
 * For non-complex types, returns a copy.
 */
function _conjugate(a: NDArray): NDArray {
  const dtype = a.dtype;

  // For non-complex types, just return a copy
  if (dtype !== DType.Complex64 && dtype !== DType.Complex128) {
    return a.copy();
  }

  // For complex types, negate the imaginary part
  // Create a new array and manually set conjugate values
  const result = a.copy();
  const size = a.size;
  const shape = a.shape;

  for (let flatIdx = 0; flatIdx < size; flatIdx++) {
    const indices = _flatToMultiIndex(flatIdx, shape);
    const complex = a.getComplex(...indices);
    // Note: setComplex signature is (real, imag, ...indices)
    result.setComplex(complex.real, -complex.imag, ...indices);
  }

  return result;
}

/**
 * Extract real part of a complex array.
 * Uses synchronous copy and modification.
 */
function _realPart(a: NDArray): NDArray {
  const dtype = a.dtype;

  // For non-complex types, just return a copy
  if (dtype !== DType.Complex64 && dtype !== DType.Complex128) {
    return a.copy();
  }

  // For complex types, we need to extract real parts
  // Use astype to convert to float64 (this will take real parts)
  return a.astype(DType.Float64);
}

/**
 * Normalize axis to positive index.
 */
function _normalizeAxis(axis: number, ndim: number): number {
  if (axis < 0) axis += ndim;
  if (axis < 0 || axis >= ndim) {
    throw new Error(
      `axis ${axis} is out of bounds for array of dimension ${ndim}`,
    );
  }
  return axis;
}

/* ============ 1D Complex FFT ============ */

/**
 * Compute the one-dimensional discrete Fourier Transform.
 *
 * @param a - Input array (can be complex or real)
 * @param n - Length of the transformed axis. If n < a.shape[axis], input is truncated.
 *            If n > a.shape[axis], input is zero-padded. Default: a.shape[axis]
 * @param axis - Axis over which to compute the FFT. Default: -1 (last axis)
 * @param norm - Normalization mode: "backward", "ortho", or "forward"
 * @returns Complex array containing the FFT result
 *
 * @example
 * // Simple FFT of a signal
 * const signal = NDArray.fromArray([1, 2, 3, 4]);
 * const spectrum = fft(signal);
 *
 * @example
 * // Zero-padded FFT for higher frequency resolution
 * const spectrum = fft(signal, 8);
 */
export function fft(
  a: NDArray,
  n: number | null = null,
  axis: number = -1,
  norm: FFTNorm = null,
): NDArray {
  if (a.ndim === 0) {
    throw new Error("FFT requires at least a 1-D array");
  }

  const module = getWasmModule();
  const normAxis = _normalizeAxis(axis, a.ndim);
  const fftSize = n ?? a.shape[normAxis];

  // Call WASM ndarray_fft (inverse=0 for forward)
  const resultPtr = module._ndarray_fft(a._wasmPtr, fftSize, normAxis, 0);
  if (resultPtr === 0) {
    throw new Error("FFT failed");
  }

  let result = NDArray._fromPtr(resultPtr, module);

  // Apply normalization adjustment if needed
  if (_needsNormAdjustment(norm)) {
    const factor = _getNormFactor(fftSize, norm, false);
    if (factor !== 1) {
      _scaleArrayInPlace(result, factor);
    }
  }

  return result;
}

/**
 * Compute the one-dimensional inverse discrete Fourier Transform.
 *
 * @param a - Input array
 * @param n - Length of the transformed axis
 * @param axis - Axis over which to compute the inverse FFT. Default: -1
 * @param norm - Normalization mode
 * @returns Complex array containing the inverse FFT result
 *
 * @example
 * // Round-trip: ifft(fft(x)) approximately equals x
 * const x = NDArray.fromArray([1, 2, 3, 4]);
 * const y = ifft(fft(x));
 */
export function ifft(
  a: NDArray,
  n: number | null = null,
  axis: number = -1,
  norm: FFTNorm = null,
): NDArray {
  if (a.ndim === 0) {
    throw new Error("IFFT requires at least a 1-D array");
  }

  const module = getWasmModule();
  const normAxis = _normalizeAxis(axis, a.ndim);
  const fftSize = n ?? a.shape[normAxis];

  // Call WASM ndarray_fft (inverse=1 for inverse)
  const resultPtr = module._ndarray_fft(a._wasmPtr, fftSize, normAxis, 1);
  if (resultPtr === 0) {
    throw new Error("IFFT failed");
  }

  let result = NDArray._fromPtr(resultPtr, module);

  // Apply normalization adjustment if needed
  if (_needsNormAdjustment(norm)) {
    const effectiveNorm = norm ?? "backward";
    let factor = 1;
    if (effectiveNorm === "ortho") {
      // WASM does 1/n, ortho needs 1/sqrt(n), so multiply by sqrt(n)
      factor = Math.sqrt(fftSize);
    } else if (effectiveNorm === "forward") {
      // WASM does 1/n, forward needs 1, so multiply by n
      factor = fftSize;
    }
    if (factor !== 1) {
      _scaleArrayInPlace(result, factor);
    }
  }

  return result;
}

/* ============ 1D Real FFT ============ */

/**
 * Compute the one-dimensional discrete Fourier Transform for real input.
 *
 * This function computes the one-dimensional n-point DFT of a real-valued
 * array. The output is Hermitian-symmetric, so only the positive frequencies
 * are returned: output has shape (..., n//2 + 1).
 *
 * @param a - Input array (real-valued)
 * @param n - Number of points in the FFT. If n < a.shape[axis], input is truncated.
 *            If n > a.shape[axis], input is zero-padded.
 * @param axis - Axis over which to compute the FFT. Default: -1
 * @param norm - Normalization mode
 * @returns Complex array with shape (..., n//2 + 1)
 *
 * @example
 * const signal = NDArray.fromArray([1, 2, 3, 4]);
 * const spectrum = rfft(signal);
 * // spectrum.shape = [3]  (4//2 + 1 = 3)
 */
export function rfft(
  a: NDArray,
  n: number | null = null,
  axis: number = -1,
  norm: FFTNorm = null,
): NDArray {
  if (a.ndim === 0) {
    throw new Error("RFFT requires at least a 1-D array");
  }

  const module = getWasmModule();
  const normAxis = _normalizeAxis(axis, a.ndim);
  const fftSize = n ?? a.shape[normAxis];

  // Call WASM ndarray_rfft
  const resultPtr = module._ndarray_rfft(a._wasmPtr, fftSize, normAxis);
  if (resultPtr === 0) {
    throw new Error("RFFT failed");
  }

  let result = NDArray._fromPtr(resultPtr, module);

  // Apply normalization adjustment if needed
  if (_needsNormAdjustment(norm)) {
    const factor = _getNormFactor(fftSize, norm, false);
    if (factor !== 1) {
      _scaleArrayInPlace(result, factor);
    }
  }

  return result;
}

/**
 * Compute the inverse of the one-dimensional discrete Fourier Transform for real input.
 *
 * This function computes the inverse of the one-dimensional n-point DFT of
 * Hermitian-symmetric input (from rfft). The input should have shape (..., n//2 + 1).
 *
 * @param a - Input array (complex, Hermitian-symmetric)
 * @param n - Length of the output. If n is not given, it is determined from the
 *            input shape: n = 2 * (a.shape[axis] - 1)
 * @param axis - Axis over which to compute the inverse FFT. Default: -1
 * @param norm - Normalization mode
 * @returns Real array with shape (..., n)
 *
 * @example
 * // Round-trip: irfft(rfft(x)) approximately equals x
 * const x = NDArray.fromArray([1, 2, 3, 4]);
 * const y = irfft(rfft(x));
 */
export function irfft(
  a: NDArray,
  n: number | null = null,
  axis: number = -1,
  norm: FFTNorm = null,
): NDArray {
  if (a.ndim === 0) {
    throw new Error("IRFFT requires at least a 1-D array");
  }

  const module = getWasmModule();
  const normAxis = _normalizeAxis(axis, a.ndim);

  // Determine output size from input if not specified
  const inputSize = a.shape[normAxis];
  const outputSize = n ?? 2 * (inputSize - 1);

  // Call WASM ndarray_irfft
  const resultPtr = module._ndarray_irfft(a._wasmPtr, outputSize, normAxis);
  if (resultPtr === 0) {
    throw new Error("IRFFT failed");
  }

  let result = NDArray._fromPtr(resultPtr, module);

  // Apply normalization adjustment if needed
  // WASM irfft applies 1/n scaling (backward convention)
  if (_needsNormAdjustment(norm)) {
    const effectiveNorm = norm ?? "backward";
    let factor = 1;
    if (effectiveNorm === "ortho") {
      factor = Math.sqrt(outputSize);
    } else if (effectiveNorm === "forward") {
      factor = outputSize;
    }
    if (factor !== 1) {
      _scaleArrayInPlace(result, factor);
    }
  }

  return result;
}

/* ============ Hermitian FFT ============ */

/**
 * Compute the FFT of a signal that has Hermitian symmetry (real spectrum).
 *
 * @param a - Input array with Hermitian symmetry
 * @param n - Length of the transformed axis
 * @param axis - Axis over which to compute the FFT. Default: -1
 * @param norm - Normalization mode
 * @returns Real array
 */
export function hfft(
  a: NDArray,
  n: number | null = null,
  axis: number = -1,
  norm: FFTNorm = null,
): NDArray {
  // hfft(a) = irfft(conj(a))
  // For Hermitian input, the FFT should give real output
  const conj_a = _conjugate(a);
  const result = irfft(conj_a, n, axis, norm);
  conj_a.dispose();
  return result;
}

/**
 * Compute the inverse FFT of a signal that has Hermitian symmetry.
 *
 * @param a - Input array (real-valued)
 * @param n - Length of the transformed axis
 * @param axis - Axis over which to compute the inverse FFT. Default: -1
 * @param norm - Normalization mode
 * @returns Complex array with Hermitian symmetry
 */
export function ihfft(
  a: NDArray,
  n: number | null = null,
  axis: number = -1,
  norm: FFTNorm = null,
): NDArray {
  // ihfft(a) = conj(rfft(a))
  const rfft_a = rfft(a, n, axis, norm);
  const result = _conjugate(rfft_a);
  rfft_a.dispose();
  return result;
}

/* ============ N-Dimensional FFT Implementation ============ */

/**
 * Internal function to apply FFT along multiple axes.
 */
function _fftn_impl(
  a: NDArray,
  s: number[] | null,
  axes: number[] | null,
  inverse: boolean,
  norm: FFTNorm,
): NDArray {
  const ndim = a.ndim;

  // Default axes: all axes
  let effectiveAxes = axes ?? [...Array(ndim).keys()];

  // Normalize negative axes
  effectiveAxes = effectiveAxes.map((ax) => _normalizeAxis(ax, ndim));

  // Default sizes: original shape along axes
  const effectiveS = s ?? effectiveAxes.map((ax) => a.shape[ax]);

  if (effectiveAxes.length !== effectiveS.length) {
    throw new Error("Shape and axes must have the same length");
  }

  if (effectiveAxes.length === 0) {
    return a.copy();
  }

  // Apply 1D FFT along each axis
  let result = a;
  let shouldDisposeResult = false;

  for (let i = 0; i < effectiveAxes.length; i++) {
    const ax = effectiveAxes[i];
    const n = effectiveS[i];

    const prevResult = result;
    if (inverse) {
      result = ifft(result, n, ax, norm);
    } else {
      result = fft(result, n, ax, norm);
    }

    if (shouldDisposeResult) {
      prevResult.dispose();
    }
    shouldDisposeResult = true;
  }

  return result;
}

/* ============ 2D FFT ============ */

/**
 * Compute the 2-dimensional discrete Fourier Transform.
 *
 * @param a - Input array
 * @param s - Shape (length of each axis) of the output. Default: shape of a
 * @param axes - Axes over which to compute the FFT. Default: (-2, -1)
 * @param norm - Normalization mode
 * @returns Complex array containing the 2D FFT result
 *
 * @example
 * const image = NDArray.zeros([64, 64]);
 * const spectrum = fft2(image);
 */
export function fft2(
  a: NDArray,
  s: [number, number] | null = null,
  axes: [number, number] = [-2, -1],
  norm: FFTNorm = null,
): NDArray {
  if (a.ndim < 2) {
    throw new Error("fft2 requires at least a 2-D array");
  }
  return _fftn_impl(a, s, axes, false, norm);
}

/**
 * Compute the 2-dimensional inverse discrete Fourier Transform.
 */
export function ifft2(
  a: NDArray,
  s: [number, number] | null = null,
  axes: [number, number] = [-2, -1],
  norm: FFTNorm = null,
): NDArray {
  if (a.ndim < 2) {
    throw new Error("ifft2 requires at least a 2-D array");
  }
  return _fftn_impl(a, s, axes, true, norm);
}

/**
 * Compute the 2-dimensional FFT of a real array.
 */
export function rfft2(
  a: NDArray,
  s: [number, number] | null = null,
  axes: [number, number] = [-2, -1],
  norm: FFTNorm = null,
): NDArray {
  if (a.ndim < 2) {
    throw new Error("rfft2 requires at least a 2-D array");
  }
  return rfftn(a, s, axes, norm);
}

/**
 * Compute the inverse of the 2-dimensional FFT of real input.
 */
export function irfft2(
  a: NDArray,
  s: [number, number] | null = null,
  axes: [number, number] = [-2, -1],
  norm: FFTNorm = null,
): NDArray {
  if (a.ndim < 2) {
    throw new Error("irfft2 requires at least a 2-D array");
  }
  return irfftn(a, s, axes, norm);
}

/* ============ N-Dimensional FFT ============ */

/**
 * Compute the N-dimensional discrete Fourier Transform.
 *
 * @param a - Input array
 * @param s - Shape of the output along transformed axes
 * @param axes - Axes over which to compute the FFT. Default: all axes
 * @param norm - Normalization mode
 * @returns Complex array containing the N-dimensional FFT result
 */
export function fftn(
  a: NDArray,
  s: number[] | null = null,
  axes: number[] | null = null,
  norm: FFTNorm = null,
): NDArray {
  return _fftn_impl(a, s, axes, false, norm);
}

/**
 * Compute the N-dimensional inverse discrete Fourier Transform.
 */
export function ifftn(
  a: NDArray,
  s: number[] | null = null,
  axes: number[] | null = null,
  norm: FFTNorm = null,
): NDArray {
  return _fftn_impl(a, s, axes, true, norm);
}

/**
 * Compute the N-dimensional discrete Fourier Transform for real input.
 *
 * The output has Hermitian symmetry along the last transformed axis.
 * Output shape along the last axis: s[-1]//2 + 1
 */
export function rfftn(
  a: NDArray,
  s: number[] | null = null,
  axes: number[] | null = null,
  norm: FFTNorm = null,
): NDArray {
  const ndim = a.ndim;

  // Default axes: all axes
  let effectiveAxes = axes ?? [...Array(ndim).keys()];
  effectiveAxes = effectiveAxes.map((ax) => _normalizeAxis(ax, ndim));

  // Default sizes
  const effectiveS = s ?? effectiveAxes.map((ax) => a.shape[ax]);

  if (effectiveAxes.length === 0) {
    // Convert to complex
    return a.astype(DType.Complex128);
  }

  // Apply real FFT on last axis first
  const lastAxis = effectiveAxes[effectiveAxes.length - 1];
  const lastN = effectiveS[effectiveS.length - 1];
  let result = rfft(a, lastN, lastAxis, norm);

  // Apply complex FFT on remaining axes
  for (let i = effectiveAxes.length - 2; i >= 0; i--) {
    const ax = effectiveAxes[i];
    const n = effectiveS[i];
    const prevResult = result;
    result = fft(result, n, ax, norm);
    prevResult.dispose();
  }

  return result;
}

/**
 * Compute the inverse of the N-dimensional FFT of real input.
 */
export function irfftn(
  a: NDArray,
  s: number[] | null = null,
  axes: number[] | null = null,
  norm: FFTNorm = null,
): NDArray {
  const ndim = a.ndim;

  // Default axes: all axes
  let effectiveAxes = axes ?? [...Array(ndim).keys()];
  effectiveAxes = effectiveAxes.map((ax) => _normalizeAxis(ax, ndim));

  // For irfftn, determine s from input if not provided
  let effectiveS: number[];
  if (s !== null) {
    effectiveS = s;
  } else {
    effectiveS = effectiveAxes.map((ax, i) => {
      if (i === effectiveAxes.length - 1) {
        // Last axis: n//2 + 1 -> n = 2 * (size - 1)
        return 2 * (a.shape[ax] - 1);
      }
      return a.shape[ax];
    });
  }

  if (effectiveAxes.length === 0) {
    // Return real part
    return _realPart(a);
  }

  // Apply complex inverse FFT on all but last axis
  let result: NDArray = a;
  let shouldDisposeResult = false;

  for (let i = 0; i < effectiveAxes.length - 1; i++) {
    const ax = effectiveAxes[i];
    const n = effectiveS[i];
    const prevResult = result;
    result = ifft(result, n, ax, norm);
    if (shouldDisposeResult) {
      prevResult.dispose();
    }
    shouldDisposeResult = true;
  }

  // Apply real inverse FFT on last axis
  const lastAxis = effectiveAxes[effectiveAxes.length - 1];
  const lastN = effectiveS[effectiveS.length - 1];
  const prevResult = result;
  result = irfft(result, lastN, lastAxis, norm);
  if (shouldDisposeResult) {
    prevResult.dispose();
  }

  return result;
}

/* ============ Helper Functions ============ */

/**
 * Return the Discrete Fourier Transform sample frequencies.
 *
 * Returns an array of length n with the sample frequencies:
 * f = [0, 1, ..., n/2-1, -n/2, ..., -1] / (d*n)   if n is even
 * f = [0, 1, ..., (n-1)/2, -(n-1)/2, ..., -1] / (d*n)   if n is odd
 *
 * @param n - Window length
 * @param d - Sample spacing (default: 1.0)
 * @returns Array of sample frequencies
 *
 * @example
 * const freqs = await fftfreq(8, 0.1);
 * // freqs = [0, 1.25, 2.5, 3.75, -5, -3.75, -2.5, -1.25]
 */
export async function fftfreq(n: number, d: number = 1.0): Promise<NDArray> {
  if (n <= 0) {
    throw new Error("n must be positive");
  }

  const freq = new Float64Array(n);
  const val = 1.0 / (n * d);

  // Number of positive frequencies (including zero)
  const N = Math.floor((n - 1) / 2) + 1;

  // Positive frequencies: 0, 1, ..., N-1
  for (let i = 0; i < N; i++) {
    freq[i] = i * val;
  }

  // Negative frequencies: -N, -(N-1), ..., -1 (if n is even, starts at -n/2)
  for (let i = N; i < n; i++) {
    freq[i] = (i - n) * val;
  }

  return NDArray.fromTypedArray(freq, [n], DType.Float64);
}

/**
 * Return the Discrete Fourier Transform sample frequencies for rfft.
 *
 * Returns an array of length n//2 + 1 with the positive sample frequencies:
 * f = [0, 1, ..., n//2] / (d*n)
 *
 * @param n - Window length
 * @param d - Sample spacing (default: 1.0)
 * @returns Array of positive sample frequencies
 *
 * @example
 * const freqs = await rfftfreq(8, 0.1);
 * // freqs = [0, 1.25, 2.5, 3.75, 5]
 */
export async function rfftfreq(n: number, d: number = 1.0): Promise<NDArray> {
  if (n <= 0) {
    throw new Error("n must be positive");
  }

  const outSize = Math.floor(n / 2) + 1;
  const freq = new Float64Array(outSize);
  const val = 1.0 / (n * d);

  for (let i = 0; i < outSize; i++) {
    freq[i] = i * val;
  }

  return NDArray.fromTypedArray(freq, [outSize], DType.Float64);
}

/**
 * Shift the zero-frequency component to the center of the spectrum.
 *
 * This function swaps half-spaces for all axes listed (or all axes by default).
 *
 * @param x - Input array
 * @param axes - Axes over which to shift. Default: all axes
 * @returns Shifted array
 *
 * @example
 * const freqs = fftfreq(4);
 * // freqs = [0, 0.25, -0.5, -0.25]
 * const shifted = fftshift(freqs);
 * // shifted = [-0.5, -0.25, 0, 0.25]
 */
export function fftshift(
  x: NDArray,
  axes: number | number[] | null = null,
): NDArray {
  const ndim = x.ndim;

  // Default: all axes
  let effectiveAxes: number[];
  if (axes === null) {
    effectiveAxes = [...Array(ndim).keys()];
  } else if (typeof axes === "number") {
    effectiveAxes = [axes];
  } else {
    effectiveAxes = axes;
  }

  // Normalize negative axes
  effectiveAxes = effectiveAxes.map((ax) => _normalizeAxis(ax, ndim));

  // Compute shift amounts: n // 2 for each axis
  const shift = effectiveAxes.map((ax) => Math.floor(x.shape[ax] / 2));

  // Use roll to perform the shift
  return roll(x, shift, effectiveAxes);
}

/**
 * The inverse of fftshift.
 *
 * Undoes the effect of fftshift.
 *
 * @param x - Input array
 * @param axes - Axes over which to shift. Default: all axes
 * @returns Shifted array
 *
 * @example
 * const x = fftfreq(4);
 * const y = fftshift(x);
 * const z = ifftshift(y);
 * // z approximately equals x
 */
export function ifftshift(
  x: NDArray,
  axes: number | number[] | null = null,
): NDArray {
  const ndim = x.ndim;

  // Default: all axes
  let effectiveAxes: number[];
  if (axes === null) {
    effectiveAxes = [...Array(ndim).keys()];
  } else if (typeof axes === "number") {
    effectiveAxes = [axes];
  } else {
    effectiveAxes = axes;
  }

  // Normalize negative axes
  effectiveAxes = effectiveAxes.map((ax) => _normalizeAxis(ax, ndim));

  // Compute shift amounts: -(n // 2) = (n - n//2) = ceil(n/2) for odd, n/2 for even
  // But we use negative shift which is equivalent to shifting by (n - n//2)
  const shift = effectiveAxes.map((ax) => -Math.floor(x.shape[ax] / 2));

  // Use roll to perform the shift
  return roll(x, shift, effectiveAxes);
}

/* ============ Module Export ============ */

/**
 * FFT module object containing all FFT functions.
 * Provides a namespace similar to numpy.fft.
 */
export const fftModule = {
  // 1D transforms
  fft,
  ifft,
  rfft,
  irfft,
  hfft,
  ihfft,

  // 2D transforms
  fft2,
  ifft2,
  rfft2,
  irfft2,

  // N-D transforms
  fftn,
  ifftn,
  rfftn,
  irfftn,

  // Helper functions
  fftfreq,
  rfftfreq,
  fftshift,
  ifftshift,
};
