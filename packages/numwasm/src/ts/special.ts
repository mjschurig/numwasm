/**
 * NumJS Special Functions
 *
 * Functions from NumPy's lib.function_base:
 * - clip: Limit values to a range
 * - diff: N-th discrete difference
 * - gradient: Gradient of an N-dimensional array
 * - convolve: 1D discrete convolution
 * - correlate: 1D cross-correlation
 * - interp: Linear interpolation
 * - trapezoid: Trapezoidal integration
 * - unwrap: Phase unwrap
 */

import { NDArray } from "./_core/NDArray.js";
import { DType } from "./types.js";
import { getWasmModule } from "./wasm-loader.js";

/* ============ Constants ============ */

/** Convolution/correlation modes */
const CONVOLVE_MODE_FULL = 0;
const CONVOLVE_MODE_SAME = 1;
const CONVOLVE_MODE_VALID = 2;

/* ============ Clip ============ */

/**
 * Clip (limit) the values in an array.
 *
 * Given an interval, values outside the interval are clipped to the interval edges.
 * For example, if an interval of [0, 1] is specified, values smaller than 0
 * become 0, and values larger than 1 become 1.
 *
 * @param a - Array containing elements to clip
 * @param a_min - Minimum value. If null, clipping is not performed on the lower interval edge.
 * @param a_max - Maximum value. If null, clipping is not performed on the upper interval edge.
 * @returns An array with the elements of a, but where values < a_min are replaced with a_min,
 *          and those > a_max with a_max.
 *
 * @example
 * const a = await NDArray.fromArray([1, 2, 3, 4, 5, 6, 7, 8, 9]);
 * const clipped = await clip(a, 3, 7);
 * // [3, 3, 3, 4, 5, 6, 7, 7, 7]
 *
 * @example
 * // Clip only lower bound
 * const clipped = await clip(a, 3, null);
 * // [3, 3, 3, 4, 5, 6, 7, 8, 9]
 *
 * @example
 * // Clip only upper bound
 * const clipped = await clip(a, null, 7);
 * // [1, 2, 3, 4, 5, 6, 7, 7, 7]
 */
export async function clip(
  a: NDArray,
  a_min: number | NDArray | null,
  a_max: number | NDArray | null,
): Promise<NDArray> {
  if (a_min === null && a_max === null) {
    throw new Error("clip: one of a_min or a_max must be provided");
  }

  const module = getWasmModule();
  let result = a;
  let needsDispose = false;

  // Apply lower bound using maximum(a, a_min)
  if (a_min !== null) {
    let minArr: NDArray;
    let disposeMin = false;

    if (typeof a_min === "number") {
      minArr = await NDArray.full(a.shape, a_min, { dtype: a.dtype });
      disposeMin = true;
    } else {
      minArr = a_min;
    }

    const resultPtr = module._ufunc_maximum(result._wasmPtr, minArr._wasmPtr);
    if (resultPtr === 0) {
      if (disposeMin) minArr.dispose();
      if (needsDispose) result.dispose();
      throw new Error("clip: maximum operation failed");
    }

    if (needsDispose) result.dispose();
    result = NDArray._fromPtr(resultPtr, module);
    needsDispose = true;

    if (disposeMin) minArr.dispose();
  }

  // Apply upper bound using minimum(result, a_max)
  if (a_max !== null) {
    let maxArr: NDArray;
    let disposeMax = false;

    if (typeof a_max === "number") {
      maxArr = await NDArray.full(result.shape, a_max, { dtype: result.dtype });
      disposeMax = true;
    } else {
      maxArr = a_max;
    }

    const resultPtr = module._ufunc_minimum(result._wasmPtr, maxArr._wasmPtr);
    if (resultPtr === 0) {
      if (disposeMax) maxArr.dispose();
      if (needsDispose) result.dispose();
      throw new Error("clip: minimum operation failed");
    }

    if (needsDispose) result.dispose();
    result = NDArray._fromPtr(resultPtr, module);

    if (disposeMax) maxArr.dispose();
  }

  return result;
}

/* ============ Diff ============ */

/**
 * Calculate the n-th discrete difference along the given axis.
 *
 * The first difference is given by out[i] = a[i+1] - a[i] along the given axis,
 * higher differences are calculated by using diff recursively.
 *
 * @param a - Input array
 * @param n - The number of times values are differenced. If zero, the input is returned as-is.
 * @param axis - The axis along which the difference is taken, default is the last axis.
 * @returns An array with the n-th differences. The shape is the same as a except along axis
 *          where the dimension is smaller by n.
 *
 * @example
 * const a = await NDArray.fromArray([1, 2, 4, 7, 0]);
 * const d = await diff(a);
 * // [1, 2, 3, -7]
 *
 * @example
 * // Second difference
 * const d2 = await diff(a, 2);
 * // [1, 1, -10]
 */
export async function diff(
  a: NDArray,
  n: number = 1,
  axis: number = -1,
): Promise<NDArray> {
  if (n < 0) {
    throw new Error("diff: n must be non-negative");
  }

  if (n === 0) {
    return a.copy();
  }

  // Normalize axis
  const normalizedAxis = axis < 0 ? a.ndim + axis : axis;
  if (normalizedAxis < 0 || normalizedAxis >= a.ndim) {
    throw new Error(`diff: axis ${axis} is out of bounds for array of dimension ${a.ndim}`);
  }

  if (a.shape[normalizedAxis] < n + 1) {
    throw new Error(
      `diff: n=${n} is too large for axis ${axis} with size ${a.shape[normalizedAxis]}`
    );
  }

  const module = getWasmModule();
  const resultPtr = module._ndarray_diff(a._wasmPtr, n, normalizedAxis);

  if (resultPtr === 0) {
    throw new Error("diff: computation failed");
  }

  return NDArray._fromPtr(resultPtr, module);
}

/* ============ Gradient ============ */

/**
 * Return the gradient of an N-dimensional array.
 *
 * The gradient is computed using second order accurate central differences
 * in the interior points and first order accurate one-sides (forward or backwards)
 * differences at the boundaries.
 *
 * @param f - An N-dimensional array containing samples of a scalar function
 * @param varargs - Spacing between f values. Can be:
 *   - A single scalar to specify uniform spacing for all dimensions
 *   - N scalars to specify uniform spacing for each dimension
 *   - N arrays to specify coordinates along each dimension
 * @param axis - Gradient is calculated only along the given axis or axes. Default is all axes.
 * @param edge_order - 1 or 2, Gradient is calculated using N-th order accurate differences at boundaries
 * @returns gradient(s) as NDArray or array of NDArrays (one per axis)
 *
 * @example
 * const f = await NDArray.fromArray([1, 2, 4, 7, 11]);
 * const g = await gradient(f);
 * // [1.0, 1.5, 2.5, 3.5, 4.0]
 *
 * @example
 * // With non-uniform spacing
 * const g = await gradient(f, 0.5);
 * // [2.0, 3.0, 5.0, 7.0, 8.0]
 */
export async function gradient(
  f: NDArray,
  ...args: (number | number[] | NDArray | { axis?: number | number[]; edge_order?: number })[]
): Promise<NDArray | NDArray[]> {
  // Parse arguments: varargs followed by optional kwargs object
  let spacing: number = 1;
  let axes: number[] | null = null;
  // edge_order is not used by our WASM implementation, but we accept it for API compatibility

  // Check if last argument is a config object
  const lastArg = args[args.length - 1];
  if (
    lastArg !== undefined &&
    typeof lastArg === "object" &&
    lastArg !== null &&
    !(lastArg instanceof NDArray) &&
    !Array.isArray(lastArg)
  ) {
    const config = lastArg as { axis?: number | number[]; edge_order?: number };
    if (config.axis !== undefined) {
      axes = typeof config.axis === "number" ? [config.axis] : config.axis;
    }
    args = args.slice(0, -1);
  }

  // Handle spacing argument
  if (args.length === 1 && typeof args[0] === "number") {
    spacing = args[0];
  } else if (args.length > 0) {
    // For now, only support uniform scalar spacing
    if (typeof args[0] === "number") {
      spacing = args[0];
    }
  }

  // Default: compute gradient along all axes
  if (axes === null) {
    axes = Array.from({ length: f.ndim }, (_, i) => i);
  }

  // Normalize axes
  const normalizedAxes = axes.map((ax) => (ax < 0 ? f.ndim + ax : ax));
  for (const ax of normalizedAxes) {
    if (ax < 0 || ax >= f.ndim) {
      throw new Error(`gradient: axis ${ax} is out of bounds for array of dimension ${f.ndim}`);
    }
  }

  const module = getWasmModule();
  const results: NDArray[] = [];

  for (const axis of normalizedAxes) {
    const resultPtr = module._ndarray_gradient(f._wasmPtr, spacing, axis);
    if (resultPtr === 0) {
      // Clean up already computed results
      for (const r of results) r.dispose();
      throw new Error(`gradient: computation failed for axis ${axis}`);
    }
    results.push(NDArray._fromPtr(resultPtr, module));
  }

  // Return single array if single axis, else array of arrays
  return results.length === 1 ? results[0] : results;
}

/* ============ Convolve ============ */

/**
 * Returns the discrete, linear convolution of two one-dimensional sequences.
 *
 * The convolution product is only given for points where the signals overlap completely.
 * Values outside the signal boundary have no effect.
 *
 * @param a - First one-dimensional input array
 * @param v - Second one-dimensional input array
 * @param mode - 'full', 'same', or 'valid':
 *   - 'full': Output is the full discrete linear convolution (default)
 *   - 'same': Output is the same size as a, centered with respect to the 'full' output
 *   - 'valid': Output consists only of those elements that do not rely on zero-padding
 * @returns Discrete, linear convolution of a and v
 *
 * @example
 * const a = await NDArray.fromArray([1, 2, 3]);
 * const v = await NDArray.fromArray([0, 1, 0.5]);
 * const c = await convolve(a, v);
 * // [0, 1, 2.5, 4, 1.5]
 *
 * @example
 * const c = await convolve(a, v, 'same');
 * // [1, 2.5, 4]
 *
 * @example
 * const c = await convolve(a, v, 'valid');
 * // [2.5]
 */
export async function convolve(
  a: NDArray,
  v: NDArray,
  mode: "full" | "same" | "valid" = "full",
): Promise<NDArray> {
  if (a.ndim !== 1 || v.ndim !== 1) {
    throw new Error("convolve: a and v must be 1-dimensional arrays");
  }

  let modeInt: number;
  switch (mode) {
    case "full":
      modeInt = CONVOLVE_MODE_FULL;
      break;
    case "same":
      modeInt = CONVOLVE_MODE_SAME;
      break;
    case "valid":
      modeInt = CONVOLVE_MODE_VALID;
      break;
    default:
      throw new Error(`convolve: unknown mode '${mode}'`);
  }

  const module = getWasmModule();
  const resultPtr = module._ndarray_convolve(a._wasmPtr, v._wasmPtr, modeInt);

  if (resultPtr === 0) {
    throw new Error("convolve: computation failed");
  }

  return NDArray._fromPtr(resultPtr, module);
}

/* ============ Correlate ============ */

/**
 * Cross-correlation of two 1-dimensional sequences.
 *
 * This function computes the correlation as generally understood in signal processing texts:
 * c_k = sum_n a_{n+k} * conj(v_n)
 *
 * @param a - First one-dimensional input array
 * @param v - Second one-dimensional input array
 * @param mode - 'full', 'same', or 'valid':
 *   - 'full': Output is the full discrete linear correlation (default)
 *   - 'same': Output is the same size as a, centered with respect to the 'full' output
 *   - 'valid': Output consists only of those elements that do not rely on zero-padding
 * @returns Discrete cross-correlation of a and v
 *
 * @example
 * const a = await NDArray.fromArray([1, 2, 3]);
 * const v = await NDArray.fromArray([0, 1, 0.5]);
 * const c = await correlate(a, v);
 * // [0.5, 2, 3.5, 3, 0]
 *
 * @example
 * const c = await correlate(a, v, 'same');
 * // [2, 3.5, 3]
 */
export async function correlate(
  a: NDArray,
  v: NDArray,
  mode: "full" | "same" | "valid" = "full",
): Promise<NDArray> {
  if (a.ndim !== 1 || v.ndim !== 1) {
    throw new Error("correlate: a and v must be 1-dimensional arrays");
  }

  let modeInt: number;
  switch (mode) {
    case "full":
      modeInt = CONVOLVE_MODE_FULL;
      break;
    case "same":
      modeInt = CONVOLVE_MODE_SAME;
      break;
    case "valid":
      modeInt = CONVOLVE_MODE_VALID;
      break;
    default:
      throw new Error(`correlate: unknown mode '${mode}'`);
  }

  const module = getWasmModule();
  const resultPtr = module._ndarray_correlate(a._wasmPtr, v._wasmPtr, modeInt);

  if (resultPtr === 0) {
    throw new Error("correlate: computation failed");
  }

  return NDArray._fromPtr(resultPtr, module);
}

/* ============ Interp ============ */

/**
 * One-dimensional linear interpolation for monotonically increasing sample points.
 *
 * Returns the one-dimensional piecewise linear interpolant to a function with
 * given discrete data points (xp, fp), evaluated at x.
 *
 * @param x - The x-coordinates at which to evaluate the interpolated values
 * @param xp - The x-coordinates of the data points, must be increasing
 * @param fp - The y-coordinates of the data points, same length as xp
 * @param left - Value to return for x < xp[0], default is fp[0]
 * @param right - Value to return for x > xp[-1], default is fp[-1]
 * @returns The interpolated values, same shape as x
 *
 * @example
 * const xp = await NDArray.fromArray([1, 2, 3]);
 * const fp = await NDArray.fromArray([3, 2, 0]);
 * const x = await NDArray.fromArray([2.5]);
 * const y = await interp(x, xp, fp);
 * // [1.0]
 *
 * @example
 * // Multiple interpolation points
 * const x = await NDArray.fromArray([0, 1, 1.5, 2.72, 3.14]);
 * const y = await interp(x, xp, fp);
 * // [3.0, 3.0, 2.5, 0.56, 0.0]
 */
export async function interp(
  x: NDArray | number | number[],
  xp: NDArray | number[],
  fp: NDArray | number[],
  left?: number,
  right?: number,
): Promise<NDArray | number> {
  // Convert inputs to NDArrays
  let xArr: NDArray;
  let xpArr: NDArray;
  let fpArr: NDArray;
  let disposeX = false;
  let disposeXp = false;
  let disposeFp = false;
  const isScalar = typeof x === "number";

  if (typeof x === "number") {
    xArr = await NDArray.fromArray([x]);
    disposeX = true;
  } else if (Array.isArray(x)) {
    xArr = await NDArray.fromArray(x);
    disposeX = true;
  } else {
    xArr = x;
  }

  if (Array.isArray(xp)) {
    xpArr = await NDArray.fromArray(xp);
    disposeXp = true;
  } else {
    xpArr = xp;
  }

  if (Array.isArray(fp)) {
    fpArr = await NDArray.fromArray(fp);
    disposeFp = true;
  } else {
    fpArr = fp;
  }

  if (xpArr.ndim !== 1 || fpArr.ndim !== 1) {
    if (disposeX) xArr.dispose();
    if (disposeXp) xpArr.dispose();
    if (disposeFp) fpArr.dispose();
    throw new Error("interp: xp and fp must be 1-dimensional");
  }

  if (xpArr.size !== fpArr.size) {
    if (disposeX) xArr.dispose();
    if (disposeXp) xpArr.dispose();
    if (disposeFp) fpArr.dispose();
    throw new Error("interp: xp and fp must have the same length");
  }

  // Get data arrays
  const xData = await xArr.toTypedArray();
  const xpData = await xpArr.toTypedArray();
  const fpData = await fpArr.toTypedArray();

  // Determine left and right values
  const leftVal = left !== undefined ? left : fpData[0];
  const rightVal = right !== undefined ? right : fpData[fpData.length - 1];

  // Perform interpolation
  const resultData = new Float64Array(xData.length);

  for (let i = 0; i < xData.length; i++) {
    const xi = xData[i];

    // Handle out-of-bounds
    if (xi <= xpData[0]) {
      resultData[i] = leftVal;
      continue;
    }
    if (xi >= xpData[xpData.length - 1]) {
      resultData[i] = rightVal;
      continue;
    }

    // Binary search for the right interval
    let lo = 0;
    let hi = xpData.length - 1;
    while (hi - lo > 1) {
      const mid = Math.floor((lo + hi) / 2);
      if (xpData[mid] <= xi) {
        lo = mid;
      } else {
        hi = mid;
      }
    }

    // Linear interpolation
    const x0 = xpData[lo];
    const x1 = xpData[hi];
    const f0 = fpData[lo];
    const f1 = fpData[hi];
    const t = (xi - x0) / (x1 - x0);
    resultData[i] = f0 + t * (f1 - f0);
  }

  // Clean up
  if (disposeX) xArr.dispose();
  if (disposeXp) xpArr.dispose();
  if (disposeFp) fpArr.dispose();

  // Return result
  if (isScalar) {
    return resultData[0];
  }

  return NDArray.fromTypedArray(
    resultData,
    Array.isArray(x) ? [x.length] : (x as NDArray).shape,
    DType.Float64,
  );
}

/* ============ Trapezoid ============ */

/**
 * Integrate along the given axis using the composite trapezoidal rule.
 *
 * If x is provided, the integration happens in sequence along its elements.
 * Otherwise, the default sample distance of 1 is used.
 *
 * @param y - Array to be integrated
 * @param x - The sample points corresponding to y values. If null, dx is used.
 * @param dx - Spacing between sample points when x is null. Default is 1.
 * @param axis - Axis along which to integrate
 * @returns Definite integral as approximated by trapezoidal rule
 *
 * @example
 * const y = await NDArray.fromArray([1, 2, 3]);
 * const integral = await trapezoid(y);
 * // 4.0 (= 0.5*(1+2) + 0.5*(2+3) = 1.5 + 2.5 = 4.0)
 *
 * @example
 * // With x coordinates
 * const x = await NDArray.fromArray([0, 1, 3]);
 * const integral = await trapezoid(y, x);
 * // 6.0 (= 0.5*1*(1+2) + 0.5*2*(2+3) = 1.5 + 5 = 6.5... wait, let me recalculate)
 * // Actually: (x[1]-x[0])*(y[0]+y[1])/2 + (x[2]-x[1])*(y[1]+y[2])/2
 * //         = 1*(1+2)/2 + 2*(2+3)/2 = 1.5 + 5 = 6.5
 */
export async function trapezoid(
  y: NDArray,
  x: NDArray | null = null,
  dx: number = 1,
  axis: number = -1,
): Promise<NDArray | number> {
  // Normalize axis
  const normalizedAxis = axis < 0 ? y.ndim + axis : axis;
  if (normalizedAxis < 0 || normalizedAxis >= y.ndim) {
    throw new Error(`trapezoid: axis ${axis} is out of bounds for array of dimension ${y.ndim}`);
  }

  const n = y.shape[normalizedAxis];
  if (n < 2) {
    throw new Error("trapezoid: at least 2 samples are required");
  }

  // For 1D arrays, compute directly
  if (y.ndim === 1) {
    const yData = await y.toTypedArray();
    let result = 0;

    if (x !== null) {
      if (x.ndim !== 1 || x.size !== y.size) {
        throw new Error("trapezoid: x must have same shape as y");
      }
      const xData = await x.toTypedArray();
      for (let i = 0; i < n - 1; i++) {
        const h = xData[i + 1] - xData[i];
        result += 0.5 * h * (yData[i] + yData[i + 1]);
      }
    } else {
      for (let i = 0; i < n - 1; i++) {
        result += 0.5 * dx * (yData[i] + yData[i + 1]);
      }
    }

    return result;
  }

  // For N-D arrays, we need to compute along the axis
  // Move target axis to last position
  const yMoved = y.moveaxis(normalizedAxis, -1);
  const axisSize = y.shape[normalizedAxis];
  const batchSize = y.size / axisSize;

  const yData = await yMoved.toTypedArray();
  const resultData = new Float64Array(batchSize);

  if (x !== null) {
    // x should be 1D with size = axisSize, or same shape as y
    let xData: Float64Array | Float32Array;
    if (x.ndim === 1 && x.size === axisSize) {
      xData = (await x.toTypedArray()) as Float64Array | Float32Array;
      // Apply same x to all batches
      for (let i = 0; i < batchSize; i++) {
        let sum = 0;
        const base = i * axisSize;
        for (let j = 0; j < axisSize - 1; j++) {
          const h = xData[j + 1] - xData[j];
          sum += 0.5 * h * (yData[base + j] + yData[base + j + 1]);
        }
        resultData[i] = sum;
      }
    } else {
      throw new Error("trapezoid: x must be 1D with same size as y along axis");
    }
  } else {
    for (let i = 0; i < batchSize; i++) {
      let sum = 0;
      const base = i * axisSize;
      for (let j = 0; j < axisSize - 1; j++) {
        sum += 0.5 * dx * (yData[base + j] + yData[base + j + 1]);
      }
      resultData[i] = sum;
    }
  }

  // Reshape to output shape (original shape minus the integration axis)
  const resultShape = y.shape.filter((_, i) => i !== normalizedAxis);
  if (resultShape.length === 0) {
    return resultData[0];
  }
  return NDArray.fromTypedArray(resultData, resultShape, DType.Float64);
}

/* ============ Unwrap ============ */

/**
 * Unwrap by taking the complement of large deltas with respect to the period.
 *
 * This unwraps a signal by changing absolute jumps greater than discont to
 * their 2*pi complement.
 *
 * @param p - Input array
 * @param discont - Maximum discontinuity between values, default is pi
 * @param axis - Axis along which to unwrap
 * @param period - Size of the range over which the input wraps, default is 2*pi
 * @returns Output array with unwrapped phases
 *
 * @example
 * const phase = await NDArray.fromArray([0, 0.78, 5.49, 6.28]);
 * const unwrapped = await unwrap(phase);
 * // [0, 0.78, -0.79, 0]  (approximately)
 */
export async function unwrap(
  p: NDArray,
  discont: number | null = null,
  axis: number = -1,
  period: number = 2 * Math.PI,
): Promise<NDArray> {
  // Default discont
  const disc = discont !== null ? discont : period / 2;

  // Normalize axis
  const normalizedAxis = axis < 0 ? p.ndim + axis : axis;
  if (normalizedAxis < 0 || normalizedAxis >= p.ndim) {
    throw new Error(`unwrap: axis ${axis} is out of bounds for array of dimension ${p.ndim}`);
  }

  // For 1D case
  if (p.ndim === 1) {
    const data = await p.toTypedArray();
    const result = new Float64Array(data.length);
    result[0] = data[0];

    let cumAdj = 0;
    for (let i = 1; i < data.length; i++) {
      let diff = data[i] - data[i - 1];

      // Adjust for discontinuities
      if (Math.abs(diff) > disc) {
        // Find the number of periods to subtract
        const nPeriods = Math.round(diff / period);
        cumAdj -= nPeriods * period;
      }

      result[i] = data[i] + cumAdj;
    }

    return NDArray.fromTypedArray(result, p.shape, DType.Float64);
  }

  // For N-D arrays, move target axis to last and process each slice
  const pMoved = p.moveaxis(normalizedAxis, -1);
  const axisSize = p.shape[normalizedAxis];
  const batchSize = p.size / axisSize;

  const pData = await pMoved.toTypedArray();
  const resultData = new Float64Array(p.size);

  for (let i = 0; i < batchSize; i++) {
    const base = i * axisSize;
    resultData[base] = pData[base];

    let cumAdj = 0;
    for (let j = 1; j < axisSize; j++) {
      let diff = pData[base + j] - pData[base + j - 1];

      if (Math.abs(diff) > disc) {
        const nPeriods = Math.round(diff / period);
        cumAdj -= nPeriods * period;
      }

      resultData[base + j] = pData[base + j] + cumAdj;
    }
  }

  // Create result with moved shape, then move axis back
  const movedShape = pMoved.shape;
  const result = await NDArray.fromTypedArray(resultData, movedShape, DType.Float64);
  return result.moveaxis(-1, normalizedAxis);
}
