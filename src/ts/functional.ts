/**
 * NumJS Functional Programming Utilities
 *
 * Level 10 implementation providing NumPy-style functional programming tools:
 * - applyAlongAxis: Apply a function to 1-D slices along an axis
 * - applyOverAxes: Apply a function repeatedly over multiple axes
 * - vectorize/Vectorize: Vectorize a scalar function with broadcasting
 * - frompyfunc: Create a ufunc-like wrapper from a function
 * - piecewise: Evaluate piecewise-defined functions
 */

import { NDArray } from './NDArray.js';
import { DType } from './types.js';
import { broadcastArrays } from './broadcast.js';
import { extract } from './indexing.js';
import { Slice, type IndexElement } from './slice.js';

/**
 * Normalize axis index, supporting negative indexing.
 *
 * @param axis - Axis index (can be negative)
 * @param ndim - Number of dimensions
 * @returns Normalized positive axis index
 * @throws Error if axis is out of bounds
 */
function normalizeAxisIndex(axis: number, ndim: number): number {
  if (axis < 0) {
    axis += ndim;
  }
  if (axis < 0 || axis >= ndim) {
    throw new Error(
      `axis ${axis} is out of bounds for array of dimension ${ndim}`
    );
  }
  return axis;
}

/**
 * Apply a function to 1-D slices along the given axis.
 *
 * Execute `func1d(a, ...args)` where `func1d` operates on 1-D arrays
 * and `a` is a 1-D slice of `arr` along `axis`.
 *
 * @param func1d - Function that accepts a 1-D array and returns a scalar or 1-D array
 * @param axis - Axis along which to apply `func1d`
 * @param arr - Input array
 * @param args - Additional arguments to pass to `func1d`
 * @returns Result array
 *
 * @example
 * ```typescript
 * // Sort each row
 * const sorted = await applyAlongAxis(
 *   async (x) => {
 *     const data = x.toArray().sort((a, b) => a - b);
 *     return NDArray.fromArray(data);
 *   },
 *   1,
 *   arr
 * );
 *
 * // Sum each column (returns scalar per column)
 * const sums = await applyAlongAxis(
 *   (x) => x.toArray().reduce((a, b) => a + b, 0),
 *   0,
 *   arr
 * );
 * ```
 */
export async function applyAlongAxis(
  func1d: (
    arr: NDArray,
    ...args: unknown[]
  ) => NDArray | number | Promise<NDArray | number>,
  axis: number,
  arr: NDArray,
  ...args: unknown[]
): Promise<NDArray> {
  const ndim = arr.ndim;

  if (ndim === 0) {
    throw new Error('applyAlongAxis requires at least 1-dimensional array');
  }

  axis = normalizeAxisIndex(axis, ndim);

  // Build shape for iteration (all dimensions except axis)
  const iterShape: number[] = [];
  for (let i = 0; i < ndim; i++) {
    if (i !== axis) {
      iterShape.push(arr.shape[i]);
    }
  }

  // Handle 1-D input: apply function directly
  if (iterShape.length === 0) {
    const result = await func1d(arr, ...args);
    if (typeof result === 'number') {
      return NDArray.fromArray([result]);
    }
    return result;
  }

  // Extract first slice to determine output shape
  const firstSlice = extractSliceAlongAxis(arr, axis, 0);
  const firstResult = await func1d(firstSlice, ...args);
  firstSlice.dispose();

  // Determine output shape based on first result
  let resultShape: number[];
  let resultIsScalar: boolean;
  let firstResultDtype: DType;
  let firstResultData: number[];

  if (typeof firstResult === 'number') {
    // Scalar result: output shape is iterShape
    resultShape = [...iterShape];
    resultIsScalar = true;
    firstResultDtype = DType.Float64;
    firstResultData = [firstResult];
  } else {
    // Array result: output shape is iterShape with result shape inserted at axis
    const resShape = firstResult.shape;
    resultShape = [];
    let iterIdx = 0;
    for (let i = 0; i < ndim; i++) {
      if (i === axis) {
        resultShape.push(...resShape);
      } else {
        resultShape.push(iterShape[iterIdx++]);
      }
    }
    // Handle case where result has different ndim than original axis size
    if (resShape.length === 0) {
      // 0-D result, remove the axis dimension
      resultShape = [...iterShape];
      resultIsScalar = true;
      firstResultDtype = firstResult.dtype;
      firstResultData = [firstResult.item()];
    } else {
      resultIsScalar = false;
      firstResultDtype = firstResult.dtype;
      firstResultData = firstResult.toArray();
    }
    firstResult.dispose();
  }

  // Allocate output buffer
  const output = await NDArray.zeros(resultShape, { dtype: firstResultDtype });

  // Store first result
  storeResult(output, iterShape, 0, firstResultData, resultIsScalar, axis);

  // Iterate over remaining indices
  let flatIdx = 1;
  const totalIter = iterShape.reduce((a, b) => a * b, 1);

  for (flatIdx = 1; flatIdx < totalIter; flatIdx++) {
    // Extract slice at this iteration index
    const slice = extractSliceAlongAxis(arr, axis, flatIdx, iterShape);
    const result = await func1d(slice, ...args);
    slice.dispose();

    // Store result
    let resultData: number[];
    if (typeof result === 'number') {
      resultData = [result];
    } else {
      resultData = result.toArray();
      result.dispose();
    }

    storeResult(output, iterShape, flatIdx, resultData, resultIsScalar, axis);
  }

  return output;
}

/**
 * Extract a 1-D slice along the specified axis at iteration index.
 */
function extractSliceAlongAxis(
  arr: NDArray,
  axis: number,
  iterFlatIdx: number,
  iterShape?: number[]
): NDArray {
  const ndim = arr.ndim;

  if (ndim === 1) {
    // For 1-D arrays, return the whole array
    return arr.view();
  }

  // Compute iteration indices from flat index
  const actualIterShape = iterShape ?? computeIterShape(arr.shape, axis);
  const iterIndices = unravelIndex(iterFlatIdx, actualIterShape);

  // Build slice indices
  const indices: IndexElement[] = [];
  let iterIdx = 0;

  for (let d = 0; d < ndim; d++) {
    if (d === axis) {
      indices.push(new Slice()); // Full slice along axis
    } else {
      indices.push(iterIndices[iterIdx++]);
    }
  }

  return arr.slice(indices);
}

/**
 * Compute iteration shape (all dims except axis).
 */
function computeIterShape(shape: number[], axis: number): number[] {
  const result: number[] = [];
  for (let i = 0; i < shape.length; i++) {
    if (i !== axis) {
      result.push(shape[i]);
    }
  }
  return result;
}

/**
 * Convert flat index to multi-dimensional indices.
 */
function unravelIndex(flatIdx: number, shape: number[]): number[] {
  const indices: number[] = new Array(shape.length);
  let remaining = flatIdx;

  for (let i = shape.length - 1; i >= 0; i--) {
    indices[i] = remaining % shape[i];
    remaining = Math.floor(remaining / shape[i]);
  }

  return indices;
}

/**
 * Store result data into output array at the correct position.
 */
function storeResult(
  output: NDArray,
  iterShape: number[],
  iterFlatIdx: number,
  data: number[],
  isScalar: boolean,
  _axis: number
): void {
  if (isScalar) {
    output.setFlat(iterFlatIdx, data[0]);
  } else {
    // For array results, we need to compute output indices
    const iterIndices = unravelIndex(iterFlatIdx, iterShape);
    const outputShape = output.shape;

    // Find where the result data should go
    // The output is laid out with iterShape dimensions surrounding the result dimensions
    for (let i = 0; i < data.length; i++) {
      // Compute flat output index
      // Result dimensions are inserted at the axis position
      let outFlatIdx = 0;
      let stride = 1;

      // Build output indices combining iterIndices and result index
      const resultNdim = output.ndim - iterShape.length;
      const resultIndices = unravelIndex(i, output.shape.slice(_axis, _axis + resultNdim));

      // Compute flat index for output
      let iterIdx = 0;
      let resIdx = 0;

      const outIndices: number[] = [];
      for (let d = 0; d < output.ndim; d++) {
        if (d >= _axis && d < _axis + resultNdim) {
          outIndices.push(resultIndices[resIdx++]);
        } else {
          outIndices.push(iterIndices[iterIdx++]);
        }
      }

      // Compute flat index from multi-index
      stride = 1;
      outFlatIdx = 0;
      for (let d = outputShape.length - 1; d >= 0; d--) {
        outFlatIdx += outIndices[d] * stride;
        stride *= outputShape[d];
      }

      output.setFlat(outFlatIdx, data[i]);
    }
  }
}

/**
 * Apply a function repeatedly over multiple axes.
 *
 * `func` is called as `func(a, axis)` where `a` is the current array
 * and `axis` is the current axis being processed. The result should
 * have the same number of dimensions as the input, or one less if
 * the function reduces along the axis (in which case `expandDims` is
 * automatically called to maintain dimensionality).
 *
 * @param func - Function to apply, signature: (a, axis) => NDArray
 * @param arr - Input array
 * @param axes - Axes over which to apply func
 * @returns Result array after applying func over all axes
 *
 * @example
 * ```typescript
 * // Sum over multiple axes, keeping dimensions
 * const result = await applyOverAxes(
 *   async (a, axis) => {
 *     // sum would return reduced array
 *     const summed = await sum(a, axis);
 *     return summed;
 *   },
 *   arr,
 *   [0, 2]
 * );
 * ```
 */
export async function applyOverAxes(
  func: (arr: NDArray, axis: number) => NDArray | Promise<NDArray>,
  arr: NDArray,
  axes: number | number[]
): Promise<NDArray> {
  const axesArr = Array.isArray(axes) ? axes : [axes];

  if (axesArr.length === 0) {
    return arr.view();
  }

  let result = arr;
  const isOriginal = (a: NDArray) => a === arr;

  for (const axis of axesArr) {
    const currentNdim = result.ndim;
    const normalizedAxis = normalizeAxisIndex(axis, currentNdim);

    // Apply function
    const newResult = await func(result, normalizedAxis);

    // Dispose previous intermediate (but not original input)
    if (!isOriginal(result)) {
      result.dispose();
    }

    // Check if dimensions were reduced
    if (newResult.ndim < currentNdim) {
      // Insert dimension at the axis position to maintain ndim
      // expandDims creates a view, so we need to make a copy to avoid
      // issues when newResult is disposed later
      const expanded = newResult.expandDims(normalizedAxis);
      // Copy the expanded view to own the data
      result = expanded.copy();
      expanded.dispose();
      newResult.dispose();
    } else {
      result = newResult;
    }
  }

  return result;
}

/**
 * Options for vectorize.
 */
export interface VectorizeOptions {
  /** Output data type(s). If not specified, determined from first call. */
  otypes?: DType | DType[] | string;
  /** Custom docstring for the vectorized function. */
  doc?: string;
  /** Arguments to exclude from vectorization (by index). */
  excluded?: Set<number>;
  /** Whether to cache the first function call result. */
  cache?: boolean;
}

/**
 * Generalized function class for vectorization.
 *
 * Takes scalars as inputs and returns a scalar or tuple of scalars.
 * The vectorized function evaluates the wrapped function over successive
 * elements of the input arrays, using NumPy's broadcasting rules.
 *
 * @example
 * ```typescript
 * // Basic vectorization
 * function myfunc(a: number, b: number): number {
 *   return a < b ? a : b;
 * }
 * const vfunc = new Vectorize(myfunc);
 * const result = await vfunc.call(
 *   await NDArray.fromArray([1, 2, 3, 4]),
 *   await NDArray.fromArray([2])
 * );
 * // result: [1, 2, 2, 2]
 *
 * // With output types
 * const vfunc2 = new Vectorize(
 *   (x: number) => x > 0,
 *   { otypes: [DType.Bool] }
 * );
 * ```
 */
export class Vectorize {
  private readonly pyfunc: (...args: unknown[]) => unknown;
  private readonly otypes: DType[] | null;
  private readonly excluded: Set<number>;
  private readonly cache: boolean;

  private _cachedResult: unknown[] | null = null;
  private _inferredOtypes: DType[] | null = null;
  private _nout: number = 0;

  constructor(
    pyfunc: (...args: unknown[]) => unknown,
    options: VectorizeOptions = {}
  ) {
    this.pyfunc = pyfunc;
    this.excluded = options.excluded ?? new Set();
    this.cache = options.cache ?? false;

    // Parse output types
    if (options.otypes !== undefined) {
      this.otypes = this.parseOtypes(options.otypes);
    } else {
      this.otypes = null;
    }
  }

  /**
   * Parse output types specification.
   */
  private parseOtypes(otypes: DType | DType[] | string): DType[] {
    if (typeof otypes === 'number') {
      return [otypes as DType];
    }
    if (Array.isArray(otypes)) {
      return otypes as DType[];
    }
    // String specification (e.g., 'ff' for two float32 outputs)
    const dtypeMap: Record<string, DType> = {
      b: DType.Bool,
      i: DType.Int32,
      l: DType.Int64,
      f: DType.Float32,
      d: DType.Float64,
    };
    return (otypes as string).split('').map((c) => {
      if (!(c in dtypeMap)) {
        throw new Error(`Unknown dtype character: ${c}`);
      }
      return dtypeMap[c];
    });
  }

  /**
   * Call the vectorized function.
   *
   * @param args - Input arrays or scalars
   * @returns Single array or array of arrays for multiple outputs
   */
  async call(...args: (NDArray | number)[]): Promise<NDArray | NDArray[]> {
    // Separate excluded and included arguments
    const includedArgs: (NDArray | number)[] = [];
    const excludedArgs: Map<number, unknown> = new Map();

    for (let i = 0; i < args.length; i++) {
      if (this.excluded.has(i)) {
        excludedArgs.set(i, args[i]);
      } else {
        includedArgs.push(args[i]);
      }
    }

    // Convert scalars to arrays
    const arrays: NDArray[] = [];
    const toDispose: NDArray[] = [];

    for (const arg of includedArgs) {
      if (typeof arg === 'number') {
        const arr = await NDArray.fromArray([arg]);
        arrays.push(arr);
        toDispose.push(arr);
      } else {
        arrays.push(arg);
      }
    }

    // Broadcast arrays to common shape
    const broadcastedArrays = await broadcastArrays(...arrays);
    const outputShape = broadcastedArrays[0].shape;
    const totalSize = broadcastedArrays[0].size;

    // Track which broadcast arrays need cleanup
    for (let i = 0; i < arrays.length; i++) {
      if (broadcastedArrays[i] !== arrays[i]) {
        toDispose.push(broadcastedArrays[i]);
      }
    }

    // Determine output types if not specified
    let outputTypes: DType[];
    let nout: number;

    if (this.otypes !== null) {
      outputTypes = this.otypes;
      nout = outputTypes.length;
    } else if (this._inferredOtypes !== null) {
      outputTypes = this._inferredOtypes;
      nout = this._nout;
    } else {
      // Infer from first call
      const firstArgs = this.buildArgsForIndex(
        0,
        broadcastedArrays,
        excludedArgs,
        args.length
      );
      const firstResult = this.pyfunc(...firstArgs);

      if (Array.isArray(firstResult)) {
        nout = firstResult.length;
        outputTypes = firstResult.map((r) => this.inferDtype(r));
        if (this.cache) {
          this._cachedResult = firstResult;
        }
      } else {
        nout = 1;
        outputTypes = [this.inferDtype(firstResult)];
        if (this.cache) {
          this._cachedResult = [firstResult];
        }
      }

      this._inferredOtypes = outputTypes;
      this._nout = nout;
    }

    // Allocate output arrays
    const outputs: NDArray[] = await Promise.all(
      outputTypes.map((dtype) => NDArray.zeros(outputShape, { dtype }))
    );

    // Iterate over all elements
    const startIdx = this.cache && this._cachedResult !== null ? 1 : 0;

    // Store cached result if applicable
    if (this.cache && this._cachedResult !== null && totalSize > 0) {
      for (let o = 0; o < nout; o++) {
        outputs[o].setFlat(0, this._cachedResult[o] as number);
      }
    }

    // Process remaining elements
    for (let i = startIdx; i < totalSize; i++) {
      const funcArgs = this.buildArgsForIndex(
        i,
        broadcastedArrays,
        excludedArgs,
        args.length
      );
      const result = this.pyfunc(...funcArgs);

      if (nout === 1) {
        const val = Array.isArray(result) ? result[0] : result;
        outputs[0].setFlat(i, val as number);
      } else {
        const results = result as unknown[];
        for (let o = 0; o < nout; o++) {
          outputs[o].setFlat(i, results[o] as number);
        }
      }
    }

    // Clean up
    for (const arr of toDispose) {
      arr.dispose();
    }

    // Return single array or tuple
    if (nout === 1) {
      return outputs[0];
    }
    return outputs;
  }

  /**
   * Build function arguments for a specific flat index.
   */
  private buildArgsForIndex(
    flatIdx: number,
    arrays: NDArray[],
    excludedArgs: Map<number, unknown>,
    totalArgs: number
  ): unknown[] {
    const result: unknown[] = [];
    let arrIdx = 0;

    for (let i = 0; i < totalArgs; i++) {
      if (excludedArgs.has(i)) {
        result.push(excludedArgs.get(i));
      } else {
        result.push(arrays[arrIdx].getFlat(flatIdx));
        arrIdx++;
      }
    }

    return result;
  }

  /**
   * Infer DType from a scalar value.
   */
  private inferDtype(value: unknown): DType {
    if (typeof value === 'boolean') {
      return DType.Bool;
    }
    if (typeof value === 'number') {
      if (Number.isInteger(value)) {
        return DType.Int32;
      }
      return DType.Float64;
    }
    return DType.Float64;
  }
}

/**
 * Factory function for creating a Vectorize instance.
 *
 * @param pyfunc - Function to vectorize
 * @param options - Vectorization options
 * @returns Vectorize instance
 *
 * @example
 * ```typescript
 * const vmin = vectorize((a: number, b: number) => Math.min(a, b));
 * const result = await vmin.call(arr1, arr2);
 * ```
 */
export function vectorize(
  pyfunc: (...args: unknown[]) => unknown,
  options?: VectorizeOptions
): Vectorize {
  return new Vectorize(pyfunc, options);
}

/**
 * A ufunc-like callable created from a JavaScript function.
 */
export interface UfuncLike {
  /** Call the ufunc on arrays. */
  (...args: (NDArray | number)[]): Promise<NDArray | NDArray[]>;

  /** Number of input arguments. */
  nin: number;

  /** Number of output values. */
  nout: number;

  /** Identity element for reductions (if any). */
  identity: unknown | null;
}

/**
 * Create a ufunc-like wrapper from a JavaScript function.
 *
 * Takes an arbitrary JavaScript function and returns a ufunc-like object
 * that can apply the function element-wise to arrays with broadcasting.
 *
 * @param func - A JavaScript function that operates on scalars
 * @param nin - The number of input arguments
 * @param nout - The number of objects returned by func
 * @param identity - Optional identity value for reduction operations
 * @returns A ufunc-like callable
 *
 * @example
 * ```typescript
 * // Simple function with 2 inputs, 1 output
 * const addone = frompyfunc((x, y) => x + y + 1, 2, 1);
 * const result = await addone(
 *   await NDArray.fromArray([1, 2, 3]),
 *   await NDArray.fromArray([10, 20, 30])
 * );
 * // result: [12, 23, 34]
 *
 * // Function with multiple outputs
 * const divmod = frompyfunc(
 *   (x, y) => [Math.floor(x / y), x % y],
 *   2, 2
 * );
 * const [quot, rem] = await divmod(
 *   await NDArray.fromArray([10, 11, 12]),
 *   await NDArray.fromArray([3])
 * );
 * ```
 */
export function frompyfunc(
  func: (...args: number[]) => number | number[],
  nin: number,
  nout: number,
  identity: unknown | null = null
): UfuncLike {
  if (nin < 1) {
    throw new Error('nin must be at least 1');
  }
  if (nout < 1) {
    throw new Error('nout must be at least 1');
  }

  const ufuncLike = async function (
    ...args: (NDArray | number)[]
  ): Promise<NDArray | NDArray[]> {
    if (args.length !== nin) {
      throw new Error(`Expected ${nin} arguments, got ${args.length}`);
    }

    // Convert scalars to arrays
    const arrays: NDArray[] = [];
    const toDispose: NDArray[] = [];

    for (const arg of args) {
      if (typeof arg === 'number') {
        const arr = await NDArray.fromArray([arg]);
        arrays.push(arr);
        toDispose.push(arr);
      } else {
        arrays.push(arg);
      }
    }

    // Broadcast to common shape
    const broadcastedArrays = await broadcastArrays(...arrays);
    const outputShape = broadcastedArrays[0].shape;
    const totalSize = broadcastedArrays[0].size;

    // Track which broadcast arrays need cleanup
    for (let i = 0; i < arrays.length; i++) {
      if (broadcastedArrays[i] !== arrays[i]) {
        toDispose.push(broadcastedArrays[i]);
      }
    }

    // Allocate output arrays (default to Float64)
    const outputs: NDArray[] = await Promise.all(
      Array(nout)
        .fill(null)
        .map(() => NDArray.zeros(outputShape, { dtype: DType.Float64 }))
    );

    // Apply function element-wise
    for (let i = 0; i < totalSize; i++) {
      const inputValues = broadcastedArrays.map((arr) => arr.getFlat(i));
      const result = func(...inputValues);

      if (nout === 1) {
        const val = Array.isArray(result) ? result[0] : result;
        outputs[0].setFlat(i, val as number);
      } else {
        const results = result as number[];
        for (let o = 0; o < nout; o++) {
          outputs[o].setFlat(i, results[o]);
        }
      }
    }

    // Clean up
    for (const arr of toDispose) {
      arr.dispose();
    }

    // Return single array or array of arrays
    if (nout === 1) {
      return outputs[0];
    }
    return outputs;
  } as UfuncLike;

  // Attach properties
  ufuncLike.nin = nin;
  ufuncLike.nout = nout;
  ufuncLike.identity = identity;

  return ufuncLike;
}

/**
 * Evaluate a piecewise-defined function.
 *
 * Given a set of conditions and corresponding functions, evaluate each
 * function on the input only where its condition is true.
 *
 * @param x - Input array (domain of the piecewise function)
 * @param condlist - List of boolean arrays defining conditions
 * @param funclist - List of functions or constants for each condition.
 *                   If funclist has one more element than condlist, the last
 *                   element is used as the default ("otherwise") case.
 * @param args - Additional arguments passed to all functions
 * @returns Array with values chosen based on conditions
 *
 * @example
 * ```typescript
 * // Absolute value implemented piecewise
 * const x = await NDArray.fromArray([-2, -1, 0, 1, 2]);
 *
 * // Create condition arrays (x < 0) and (x >= 0)
 * const cond1 = await less(x, await NDArray.fromArray([0]));
 * const cond2 = await greater_equal(x, await NDArray.fromArray([0]));
 *
 * const result = await piecewise(
 *   x,
 *   [cond1, cond2],
 *   [
 *     async (vals) => negative(vals),  // negate for x < 0
 *     async (vals) => vals.copy()      // identity for x >= 0
 *   ]
 * );
 * // result: [2, 1, 0, 1, 2]
 *
 * // With default value (one more function than conditions)
 * const y = await piecewise(
 *   x,
 *   [cond1],
 *   [
 *     -1,    // constant for x < 0
 *     1      // default for all other cases
 *   ]
 * );
 * // result: [-1, -1, 1, 1, 1]
 * ```
 */
export async function piecewise(
  x: NDArray,
  condlist: NDArray[],
  funclist: (
    | ((arr: NDArray, ...args: unknown[]) => NDArray | Promise<NDArray>)
    | number
  )[],
  ...args: unknown[]
): Promise<NDArray> {
  // Validate inputs
  if (condlist.length === 0) {
    throw new Error('condlist must contain at least one condition');
  }

  // Handle default case: len(funclist) == len(condlist) + 1
  const conditions: NDArray[] = [...condlist];
  const functions = [...funclist];
  let createdOtherwiseCondition = false;

  if (funclist.length === condlist.length + 1) {
    // Last function is the "otherwise" case
    // Create condition: ~(cond1 | cond2 | ... | condN)
    const otherwiseCondition = await computeOtherwiseCondition(condlist);
    conditions.push(otherwiseCondition);
    createdOtherwiseCondition = true;
  } else if (funclist.length !== condlist.length) {
    throw new Error(
      `funclist must have same length as condlist (${condlist.length}) ` +
        `or one more element`
    );
  }

  // Initialize output with zeros
  const output = await NDArray.zerosLike(x);

  // Apply each function where its condition is true
  for (let i = 0; i < conditions.length; i++) {
    const condition = conditions[i];
    const func = functions[i];

    // Check if any elements satisfy this condition
    const condData = condition.toArray();
    const hasTrue = condData.some((v) => v !== 0);

    if (!hasTrue) {
      continue;
    }

    if (typeof func === 'number') {
      // Constant value
      setWhereCondition(output, condData, func);
    } else {
      // Extract values where condition is true
      const vals = await extract(condition, x);

      if (vals.size > 0) {
        // Apply function
        const result = await func(vals, ...args);

        // Put results back into output at condition positions
        putWhereCondition(output, condData, result);

        result.dispose();
      }

      vals.dispose();
    }
  }

  // Clean up "otherwise" condition if we created it
  if (createdOtherwiseCondition) {
    conditions[conditions.length - 1].dispose();
  }

  return output;
}

/**
 * Compute the "otherwise" condition: ~(cond1 | cond2 | ... | condN)
 */
async function computeOtherwiseCondition(
  condlist: NDArray[]
): Promise<NDArray> {
  const size = condlist[0].size;
  const shape = condlist[0].shape;

  // Compute OR of all conditions
  const anyTrue = new Array(size).fill(false);

  for (const cond of condlist) {
    const condData = cond.toArray();
    for (let i = 0; i < size; i++) {
      if (condData[i] !== 0) {
        anyTrue[i] = true;
      }
    }
  }

  // Negate to get "otherwise" condition
  const otherwiseData = anyTrue.map((v) => (v ? 0 : 1));

  return NDArray.fromArray(otherwiseData, shape, { dtype: DType.Bool });
}

/**
 * Set values in array where condition is true.
 */
function setWhereCondition(
  arr: NDArray,
  condData: number[],
  value: number
): void {
  for (let i = 0; i < arr.size; i++) {
    if (condData[i] !== 0) {
      arr.setFlat(i, value);
    }
  }
}

/**
 * Put values into array at positions where condition is true.
 */
function putWhereCondition(
  arr: NDArray,
  condData: number[],
  values: NDArray
): void {
  const valData = values.toArray();
  let valIdx = 0;

  for (let i = 0; i < arr.size; i++) {
    if (condData[i] !== 0) {
      arr.setFlat(i, valData[valIdx++]);
    }
  }
}
