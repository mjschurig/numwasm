# Level 10: Functional Programming Implementation Plan

Level 10 implements NumPy's functional programming utilities that enable applying functions across array dimensions, vectorizing Python functions, and conditional evaluation. These tools bridge custom functions with NumJS's array infrastructure.

---

## Current State (Level 9 Complete)

```
src/wasm/
├── ndarray.h/c        # Core NDArray with views, slicing, shape ops
├── dtype.h/c          # DType utilities and promotion
├── broadcast.h/c      # Broadcasting infrastructure
├── indexing.h/c       # Index functions (take, put, where, etc.)
└── pairwise_sum.h/c   # Accurate summation algorithm

src/ts/
├── types.ts           # DType enum, WasmModule interface, flags
├── NDArray.ts         # Factory methods, element access, shape ops, slicing
├── iterators.ts       # FlatIterator, nditer, ndenumerate, ndindex
├── dtype.ts           # Type promotion utilities
├── broadcast.ts       # Broadcasting functions
├── indexing.ts        # Index operations
├── slice.ts           # Slice class, ellipsis, newaxis
├── wasm-loader.ts     # WASM module loading
└── index.ts           # Public exports
```

**Existing Infrastructure:**
- NDArray creation: `zeros()`, `ones()`, `empty()`, `full()`, `fromArray()`, `arange()`, `linspace()`
- Element access: `get()`, `set()`, `getFlat()`, `setFlat()`, `item()`
- Shape manipulation: `reshape()`, `transpose()`, `ravel()`, `flatten()`, `squeeze()`, `expandDims()`
- Slicing: `slice()`, `at()` methods
- Broadcasting: `broadcastTo()`, `broadcastArrays()`, `broadcastShapes()`
- Index functions: `take()`, `put()`, `nonzero()`, `where()`, `compress()`, `extract()`
- Iteration: `nditer()`, `ndenumerate()`, `ndindex()`, `FlatIterator`

---

## NumPy Reference

```
/numpy/numpy/lib/_shape_base_impl.py
├── apply_along_axis()    # Lines 277-418 (~141 lines)
└── apply_over_axes()     # Lines 422-495 (~73 lines)

/numpy/numpy/lib/_function_base_impl.py
├── vectorize class       # Lines 2268-2572 (~304 lines)
└── piecewise()           # Lines 686-802 (~116 lines)

/numpy/numpy/_core/_multiarray_umath (C extension)
└── frompyfunc()          # C implementation - wraps Python callable as ufunc
```

---

## Level 10 Implementation Tree

```
LEVEL 10: FUNCTIONAL PROGRAMMING
│
├── 10.1 Apply Along Axis (TypeScript)
│   ├── 10.1.1 applyAlongAxis() - core implementation
│   ├── 10.1.2 Axis normalization and validation
│   ├── 10.1.3 Multi-dimensional index iteration
│   ├── 10.1.4 Result buffer allocation
│   └── 10.1.5 Output shape computation
│
│   Dependencies: iterators.ts (ndindex), NDArray shape ops
│
├── 10.2 Apply Over Axes (TypeScript)
│   ├── 10.2.1 applyOverAxes() - sequential axis application
│   ├── 10.2.2 Dimension preservation logic
│   └── 10.2.3 Axis normalization for multiple axes
│
│   Dependencies: 10.1 (conceptually similar)
│
├── 10.3 Vectorize (TypeScript)
│   ├── 10.3.1 Vectorize class
│   ├── 10.3.2 Output type inference
│   ├── 10.3.3 Excluded argument handling
│   ├── 10.3.4 Caching support
│   ├── 10.3.5 Generalized ufunc signature parsing (optional)
│   └── 10.3.6 vectorize() factory function
│
│   Dependencies: broadcast.ts, NDArray iteration
│
├── 10.4 FromPyFunc (TypeScript)
│   ├── 10.4.1 frompyfunc() - create ufunc-like wrapper
│   ├── 10.4.2 Input/output count handling
│   ├── 10.4.3 Broadcasting application
│   └── 10.4.4 Identity element support
│
│   Dependencies: broadcast.ts
│
└── 10.5 Piecewise (TypeScript)
    ├── 10.5.1 piecewise() - conditional function application
    ├── 10.5.2 Condition list handling
    ├── 10.5.3 Default case ("otherwise") logic
    └── 10.5.4 Boolean indexing integration


    Dependencies: indexing.ts (where), NDArray boolean operations
```

---

## Detailed Implementation Specifications

### 10.1 Apply Along Axis

**File:** `src/ts/functional.ts` (new file)

```typescript
import { NDArray } from './NDArray.js';
import { ndindex } from './iterators.js';

/**
 * Normalize axis index, supporting negative indexing.
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
 * Execute `func1d(a, *args, **kwargs)` where `func1d` operates on 1-D arrays
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
 *   (x) => NDArray.fromArray(x.toArray().sort((a, b) => a - b)),
 *   1,
 *   arr
 * );
 *
 * // Custom function that returns scalar
 * const result = await applyAlongAxis(
 *   (x) => x.get(0) + x.get(-1),  // First + last element
 *   0,
 *   arr
 * );
 * ```
 */
export async function applyAlongAxis(
  func1d: (arr: NDArray, ...args: unknown[]) => NDArray | number | Promise<NDArray | number>,
  axis: number,
  arr: NDArray,
  ...args: unknown[]
): Promise<NDArray> {
  const ndim = arr.ndim;
  axis = normalizeAxisIndex(axis, ndim);

  // Build shape for iteration (all dimensions except axis)
  const iterShape: number[] = [];
  for (let i = 0; i < ndim; i++) {
    if (i !== axis) {
      iterShape.push(arr.shape[i]);
    }
  }

  // Handle scalar result case (single element along axis)
  if (iterShape.length === 0) {
    // arr is 1-D, apply function directly
    const result = await func1d(arr, ...args);
    if (typeof result === 'number') {
      return NDArray.fromArray([result]);
    }
    return result as NDArray;
  }

  // Build slice template for extracting 1-D slices
  // We'll use slice() with appropriate indices
  const sliceTemplate: (number | null)[] = new Array(ndim).fill(null);

  // Get first result to determine output shape
  const firstIndices = new Array(iterShape.length).fill(0);
  let sliceIdx = 0;
  for (let i = 0; i < ndim; i++) {
    if (i === axis) {
      sliceTemplate[i] = null; // Full slice along this axis
    } else {
      sliceTemplate[i] = firstIndices[sliceIdx++];
    }
  }

  // Extract first 1-D slice
  const firstSlice = extractSlice(arr, sliceTemplate, axis);
  const firstResult = await func1d(firstSlice, ...args);
  firstSlice.dispose();

  // Determine output shape based on first result
  let resultShape: number[];
  let resultIsScalar: boolean;
  let firstResultArr: NDArray;

  if (typeof firstResult === 'number') {
    // Scalar result: output shape is iterShape
    resultShape = [...iterShape];
    resultIsScalar = true;
    firstResultArr = await NDArray.fromArray([firstResult]);
  } else {
    // Array result: output shape includes result dimensions
    firstResultArr = firstResult;
    const resShape = firstResultArr.shape;

    // Insert result shape at the axis position
    resultShape = [];
    let resIdx = 0;
    for (let i = 0; i < ndim; i++) {
      if (i === axis) {
        resultShape.push(...resShape);
      } else {
        resultShape.push(arr.shape[i]);
      }
    }
    resultIsScalar = false;
  }

  // Allocate output buffer
  const output = await NDArray.zeros(resultShape, { dtype: firstResultArr.dtype });

  // Store first result
  if (resultIsScalar) {
    output.setFlat(0, firstResultArr.item());
  } else {
    copyIntoSlice(output, firstIndices, axis, firstResultArr);
  }
  firstResultArr.dispose();

  // Iterate over remaining indices
  let flatIdx = 1;
  for (const indices of ndindex(...iterShape)) {
    if (flatIdx === 0) {
      flatIdx++;
      continue; // Skip first (already processed)
    }

    // Build slice indices
    sliceIdx = 0;
    for (let i = 0; i < ndim; i++) {
      if (i !== axis) {
        sliceTemplate[i] = indices[sliceIdx++];
      }
    }

    // Extract slice and apply function
    const slice = extractSlice(arr, sliceTemplate, axis);
    const result = await func1d(slice, ...args);
    slice.dispose();

    // Store result
    if (resultIsScalar) {
      const val = typeof result === 'number' ? result : (result as NDArray).item();
      output.setFlat(flatIdx, val);
      if (typeof result !== 'number') {
        (result as NDArray).dispose();
      }
    } else {
      const resArr = result as NDArray;
      copyIntoSlice(output, Array.from(indices), axis, resArr);
      resArr.dispose();
    }

    flatIdx++;
  }

  return output;
}

/**
 * Extract a 1-D slice from array.
 * sliceTemplate has nulls for the axis dimension and integers for others.
 */
function extractSlice(arr: NDArray, template: (number | null)[], axis: number): NDArray {
  // Build index specifications
  const { slice, Slice } = require('./slice.js');

  const indices: IndexElement[] = [];
  for (let i = 0; i < template.length; i++) {
    if (template[i] === null) {
      indices.push(new Slice()); // Full slice
    } else {
      indices.push(template[i] as number);
    }
  }

  return arr.slice(indices);
}

/**
 * Copy result array into output at the given iteration indices.
 */
function copyIntoSlice(
  output: NDArray,
  iterIndices: number[],
  axis: number,
  result: NDArray
): void {
  // For simplicity, iterate element by element
  // Could be optimized with bulk copy for contiguous cases
  const resultData = result.toArray();
  const resultShape = result.shape;

  // Build output indices
  const outNdim = output.ndim;
  const resNdim = result.ndim;

  for (let i = 0; i < resultData.length; i++) {
    // Convert flat result index to multi-index
    const resIndices: number[] = [];
    let remaining = i;
    for (let d = resNdim - 1; d >= 0; d--) {
      resIndices.unshift(remaining % resultShape[d]);
      remaining = Math.floor(remaining / resultShape[d]);
    }

    // Build output indices
    const outIndices: number[] = [];
    let iterIdx = 0;
    let resIdx = 0;

    for (let d = 0; d < outNdim; d++) {
      // Check if this dimension is the expanded axis region
      if (d >= axis && d < axis + resNdim) {
        outIndices.push(resIndices[resIdx++]);
      } else {
        outIndices.push(iterIndices[iterIdx++]);
      }
    }

    output.set(resultData[i], ...outIndices);
  }
}
```

---

### 10.2 Apply Over Axes

**File:** `src/ts/functional.ts` (continued)

```typescript
/**
 * Apply a function repeatedly over multiple axes.
 *
 * `func` is called as `func(a, axis)` where `a` is the current array
 * and `axis` is the current axis being processed. The result should
 * have the same number of dimensions as the input.
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
 *     const summed = await sum(a, axis);
 *     return summed.expandDims(axis);
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
  // Normalize axes to array
  const axesArr = Array.isArray(axes) ? axes : [axes];

  let result = arr;
  const originalNdim = arr.ndim;

  for (const axis of axesArr) {
    // Normalize axis for current array dimensions
    const normalizedAxis = normalizeAxisIndex(axis, result.ndim);

    // Apply function
    const newResult = await func(result, normalizedAxis);

    // Check if dimensions were reduced
    if (newResult.ndim < result.ndim) {
      // Insert dimension at the axis position to maintain ndim
      const expanded = newResult.expandDims(normalizedAxis);
      newResult.dispose();

      // Dispose intermediate result (but not original input)
      if (result !== arr) {
        result.dispose();
      }
      result = expanded;
    } else {
      // Dispose intermediate result (but not original input)
      if (result !== arr) {
        result.dispose();
      }
      result = newResult;
    }
  }

  return result;
}
```

---

### 10.3 Vectorize Class

**File:** `src/ts/functional.ts` (continued)

```typescript
/**
 * Output type specification for vectorize.
 */
type OutputTypes = DType | DType[] | string;

/**
 * Options for vectorize.
 */
interface VectorizeOptions {
  /** Output data type(s). If not specified, determined from first call. */
  otypes?: OutputTypes;
  /** Custom docstring for the vectorized function. */
  doc?: string;
  /** Arguments to exclude from vectorization (indices or names). */
  excluded?: Set<number | string>;
  /** Whether to cache the first function call result. */
  cache?: boolean;
  /** Generalized ufunc signature (e.g., '(m,n),(n)->(m)'). */
  signature?: string;
}

/**
 * Generalized function class for vectorization.
 *
 * Takes a nested sequence of objects or scalars as inputs and returns a
 * single ndarray or tuple of ndarrays. The vectorized function evaluates
 * `pyfunc` over successive tuples of the input arrays like a Python map
 * function, except it uses the broadcasting rules of NumPy.
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
 *   (x) => x > 0,
 *   { otypes: [DType.Bool] }
 * );
 * ```
 */
export class Vectorize {
  private readonly pyfunc: (...args: unknown[]) => unknown;
  private readonly otypes: DType[] | null;
  private readonly excluded: Set<number | string>;
  private readonly cache: boolean;
  private readonly signature: string | null;
  private readonly doc: string | null;

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
    this.signature = options.signature ?? null;
    this.doc = options.doc ?? null;

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
  private parseOtypes(otypes: OutputTypes): DType[] {
    if (typeof otypes === 'number') {
      return [otypes as DType];
    }
    if (Array.isArray(otypes)) {
      return otypes as DType[];
    }
    // String specification (e.g., 'ff' for two float32 outputs)
    const dtypeMap: Record<string, DType> = {
      'b': DType.Bool,
      'i': DType.Int32,
      'l': DType.Int64,
      'f': DType.Float32,
      'd': DType.Float64,
    };
    return (otypes as string).split('').map(c => {
      if (!(c in dtypeMap)) {
        throw new Error(`Unknown dtype character: ${c}`);
      }
      return dtypeMap[c];
    });
  }

  /**
   * Call the vectorized function.
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
    const arrays: NDArray[] = await Promise.all(
      includedArgs.map(async (arg) => {
        if (typeof arg === 'number') {
          return NDArray.fromArray([arg]);
        }
        return arg;
      })
    );

    // Broadcast arrays to common shape
    const broadcastedArrays = broadcastArrays(...arrays);
    const outputShape = broadcastedArrays[0].shape;
    const totalSize = broadcastedArrays[0].size;

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
        0, broadcastedArrays, excludedArgs, args.length
      );
      const firstResult = this.pyfunc(...firstArgs);

      if (Array.isArray(firstResult)) {
        nout = firstResult.length;
        outputTypes = firstResult.map(r => this.inferDtype(r));
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
      outputTypes.map(dtype => NDArray.zeros(outputShape, { dtype }))
    );

    // Iterate over all elements
    const startIdx = (this.cache && this._cachedResult !== null) ? 1 : 0;

    // Store cached result if applicable
    if (this.cache && this._cachedResult !== null) {
      for (let o = 0; o < nout; o++) {
        outputs[o].setFlat(0, this._cachedResult[o] as number);
      }
    }

    // Process remaining elements
    for (let i = startIdx; i < totalSize; i++) {
      const funcArgs = this.buildArgsForIndex(
        i, broadcastedArrays, excludedArgs, args.length
      );
      const result = this.pyfunc(...funcArgs);

      if (nout === 1) {
        outputs[0].setFlat(i, result as number);
      } else {
        const results = result as unknown[];
        for (let o = 0; o < nout; o++) {
          outputs[o].setFlat(i, results[o] as number);
        }
      }
    }

    // Clean up temporary broadcast views
    for (let i = 0; i < arrays.length; i++) {
      if (broadcastedArrays[i] !== arrays[i]) {
        broadcastedArrays[i].dispose();
      }
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
 * Factory function for vectorize (can be used as decorator in JS).
 */
export function vectorize(
  pyfunc: (...args: unknown[]) => unknown,
  options?: VectorizeOptions
): Vectorize {
  return new Vectorize(pyfunc, options);
}
```

---

### 10.4 FromPyFunc

**File:** `src/ts/functional.ts` (continued)

```typescript
/**
 * A ufunc-like callable created from a Python/JS function.
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
 * @param func - A JavaScript function
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

  const ufuncLike = async function(...args: (NDArray | number)[]): Promise<NDArray | NDArray[]> {
    if (args.length !== nin) {
      throw new Error(`Expected ${nin} arguments, got ${args.length}`);
    }

    // Convert scalars to arrays
    const arrays: NDArray[] = await Promise.all(
      args.map(async (arg) => {
        if (typeof arg === 'number') {
          return NDArray.fromArray([arg]);
        }
        return arg;
      })
    );

    // Broadcast to common shape
    const broadcastedArrays = broadcastArrays(...arrays);
    const outputShape = broadcastedArrays[0].shape;
    const totalSize = broadcastedArrays[0].size;

    // Allocate output arrays (default to Float64)
    const outputs: NDArray[] = await Promise.all(
      Array(nout).fill(null).map(() =>
        NDArray.zeros(outputShape, { dtype: DType.Float64 })
      )
    );

    // Apply function element-wise
    for (let i = 0; i < totalSize; i++) {
      const inputValues = broadcastedArrays.map(arr => arr.getFlat(i));
      const result = func(...inputValues);

      if (nout === 1) {
        outputs[0].setFlat(i, result as number);
      } else {
        const results = result as number[];
        for (let o = 0; o < nout; o++) {
          outputs[o].setFlat(i, results[o]);
        }
      }
    }

    // Clean up temporary broadcast views
    for (let i = 0; i < arrays.length; i++) {
      if (broadcastedArrays[i] !== arrays[i]) {
        broadcastedArrays[i].dispose();
      }
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
```

---

### 10.5 Piecewise

**File:** `src/ts/functional.ts` (continued)

```typescript
/**
 * Evaluate a piecewise-defined function.
 *
 * Given a set of conditions and corresponding functions, evaluate each
 * function on the input only where its condition is true.
 *
 * @param x - Input array (domain of the piecewise function)
 * @param condlist - List of boolean arrays defining conditions
 * @param funclist - List of functions or constants for each condition
 * @param args - Additional arguments passed to all functions
 * @returns Array with values chosen based on conditions
 *
 * @example
 * ```typescript
 * // Absolute value implemented piecewise
 * const x = await NDArray.fromArray([-2, -1, 0, 1, 2]);
 * const condlist = [
 *   await lessThan(x, 0),  // x < 0
 *   await greaterEqual(x, 0)  // x >= 0
 * ];
 * const funclist = [
 *   (x) => x.mul(-1),  // negate for negative
 *   (x) => x           // identity for positive
 * ];
 * const result = await piecewise(x, condlist, funclist);
 * // result: [2, 1, 0, 1, 2]
 *
 * // With default value (one more function than conditions)
 * const y = await piecewise(
 *   x,
 *   [await lessThan(x, 0)],
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
  funclist: (((arr: NDArray, ...args: unknown[]) => NDArray | Promise<NDArray>) | number)[],
  ...args: unknown[]
): Promise<NDArray> {
  // Validate inputs
  if (condlist.length === 0) {
    throw new Error('condlist must contain at least one condition');
  }

  // Handle default case: len(funclist) == len(condlist) + 1
  let conditions: NDArray[] = [...condlist];
  let functions = [...funclist];

  if (funclist.length === condlist.length + 1) {
    // Last function is the "otherwise" case
    // Create condition: ~(cond1 | cond2 | ... | condN)

    // Compute union of all conditions
    let anyCondition = await NDArray.zerosLike(condlist[0], { dtype: DType.Bool });

    for (const cond of condlist) {
      // OR operation (will be implemented in ufuncs, for now use element-wise)
      const combined = await logicalOr(anyCondition, cond);
      anyCondition.dispose();
      anyCondition = combined;
    }

    // Negate to get "otherwise" condition
    const otherwiseCondition = await logicalNot(anyCondition);
    anyCondition.dispose();

    conditions.push(otherwiseCondition);
  } else if (funclist.length !== condlist.length) {
    throw new Error(
      `funclist must have same length as condlist (${condlist.length}) ` +
      `or one more element`
    );
  }

  // Validate condition shapes match x
  for (let i = 0; i < conditions.length; i++) {
    if (!arraysShapeEqual(conditions[i].shape, x.shape)) {
      // Broadcast condition to x's shape
      conditions[i] = broadcastTo(conditions[i], x.shape);
    }
  }

  // Initialize output with zeros
  const output = await NDArray.zerosLike(x);

  // Apply each function where its condition is true
  for (let i = 0; i < conditions.length; i++) {
    const condition = conditions[i];
    const func = functions[i];

    // Find indices where condition is true
    const indices = await nonzero(condition);

    if (indices[0].size === 0) {
      // No elements satisfy this condition
      continue;
    }

    if (typeof func === 'number') {
      // Constant value
      await setWhereCondition(output, condition, func);
    } else {
      // Extract values where condition is true
      const vals = await extract(condition, x);

      // Apply function
      const result = await func(vals, ...args);

      // Put results back into output at condition positions
      await putWhereCondition(output, condition, result);

      vals.dispose();
      result.dispose();
    }
  }

  // Clean up "otherwise" condition if we created it
  if (funclist.length === condlist.length + 1) {
    conditions[conditions.length - 1].dispose();
  }

  return output;
}

/**
 * Set values in array where condition is true.
 */
async function setWhereCondition(
  arr: NDArray,
  condition: NDArray,
  value: number
): Promise<void> {
  const condData = condition.toArray();
  for (let i = 0; i < arr.size; i++) {
    if (condData[i]) {
      arr.setFlat(i, value);
    }
  }
}

/**
 * Put values into array at positions where condition is true.
 */
async function putWhereCondition(
  arr: NDArray,
  condition: NDArray,
  values: NDArray
): Promise<void> {
  const condData = condition.toArray();
  const valData = values.toArray();
  let valIdx = 0;

  for (let i = 0; i < arr.size; i++) {
    if (condData[i]) {
      arr.setFlat(i, valData[valIdx++]);
    }
  }
}

/**
 * Check if two shapes are equal.
 */
function arraysShapeEqual(a: number[], b: number[]): boolean {
  if (a.length !== b.length) return false;
  return a.every((v, i) => v === b[i]);
}

// Placeholder for logical operations (to be implemented in ufuncs)
async function logicalOr(a: NDArray, b: NDArray): Promise<NDArray> {
  const result = await NDArray.zeros(a.shape, { dtype: DType.Bool });
  const aData = a.toArray();
  const bData = b.toArray();
  for (let i = 0; i < a.size; i++) {
    result.setFlat(i, aData[i] || bData[i] ? 1 : 0);
  }
  return result;
}

async function logicalNot(a: NDArray): Promise<NDArray> {
  const result = await NDArray.zeros(a.shape, { dtype: DType.Bool });
  const aData = a.toArray();
  for (let i = 0; i < a.size; i++) {
    result.setFlat(i, aData[i] ? 0 : 1);
  }
  return result;
}
```

---

## File Changes Summary

### New Files to Create

```
src/ts/
└── functional.ts      # All functional programming utilities
    ├── normalizeAxisIndex()
    ├── applyAlongAxis()
    ├── applyOverAxes()
    ├── Vectorize class
    ├── vectorize()
    ├── frompyfunc()
    ├── piecewise()
    └── Helper functions
```

### Files to Modify

```
src/ts/index.ts
├── Export applyAlongAxis
├── Export applyOverAxes
├── Export Vectorize, vectorize
├── Export frompyfunc
└── Export piecewise
```

---

## WasmModule Interface (No Changes Required)

Level 10 is primarily TypeScript-based and builds on existing WASM functions. No new C/WASM functions are required.

However, efficient implementations may benefit from future WASM optimizations:

```typescript
// Potential future WASM functions for optimization:
// _ndarray_apply_along_axis()  - C-level iteration
// _ndarray_piecewise()         - Optimized conditional assignment
```

---

## Implementation Order

```
Week 1: Core Functions
├── Day 1: normalizeAxisIndex(), helper functions
├── Day 2: applyAlongAxis() core implementation
├── Day 3: applyAlongAxis() edge cases and tests
├── Day 4: applyOverAxes() implementation
└── Day 5: applyOverAxes() tests

Week 2: Vectorize
├── Day 1: Vectorize class structure
├── Day 2: Output type inference
├── Day 3: Broadcasting integration
├── Day 4: Excluded arguments & caching
└── Day 5: Tests and documentation

Week 3: FromPyFunc and Piecewise
├── Day 1: frompyfunc() implementation
├── Day 2: frompyfunc() tests
├── Day 3: piecewise() core implementation
├── Day 4: piecewise() default case handling
└── Day 5: Integration tests

Week 4: Polish and Optimization
├── Day 1: Performance benchmarks
├── Day 2: Edge case handling
├── Day 3: Documentation and examples
├── Day 4: TypeScript type refinements
└── Day 5: Final review and cleanup
```

---

## Verification Plan

After Level 10 completion, verify:

```bash
# Build
npm run build

# Run all tests
npm test

# Specific Level 10 tests should pass:

# applyAlongAxis
✓ applyAlongAxis() applies function along axis 0
✓ applyAlongAxis() applies function along axis 1
✓ applyAlongAxis() handles negative axis
✓ applyAlongAxis() handles scalar results
✓ applyAlongAxis() handles array results (different shape)
✓ applyAlongAxis() passes additional arguments

# applyOverAxes
✓ applyOverAxes() applies function over single axis
✓ applyOverAxes() applies function over multiple axes
✓ applyOverAxes() maintains dimensions when result reduces

# Vectorize
✓ vectorize() creates callable from function
✓ Vectorize.call() broadcasts inputs
✓ Vectorize with otypes specifies output dtype
✓ Vectorize with excluded arguments works correctly
✓ Vectorize with cache option caches first result
✓ vectorize() handles multiple outputs

# frompyfunc
✓ frompyfunc() creates ufunc-like callable
✓ frompyfunc() with nin=2 broadcasts two arrays
✓ frompyfunc() with nout=2 returns tuple
✓ frompyfunc() identity property is set

# piecewise
✓ piecewise() evaluates constant functions
✓ piecewise() evaluates callable functions
✓ piecewise() handles default case (otherwise)
✓ piecewise() broadcasts conditions to input shape
✓ piecewise() passes additional arguments to functions
```

Generate NumPy comparison vectors:

```python
# tests/python/generate_level10_tests.py
import numpy as np
import json

def sorted_sum(x):
    """Sum after sorting."""
    return np.sum(np.sort(x))

tests = {
    "apply_along_axis": [
        {
            "name": "sort_rows",
            "data": [[3, 1, 2], [6, 4, 5]],
            "axis": 1,
            "func": "sort",
            "expected": [[1, 2, 3], [4, 5, 6]]
        },
        {
            "name": "sum_cols_scalar",
            "data": [[1, 2], [3, 4], [5, 6]],
            "axis": 0,
            "func": "sum",
            "expected": [9, 12]
        },
    ],
    "apply_over_axes": [
        {
            "name": "sum_axes_0_2",
            "shape": [2, 3, 4],
            "axes": [0, 2],
            "func": "sum",
            "expected_shape": [1, 3, 1]
        },
    ],
    "vectorize": [
        {
            "name": "min_func",
            "a": [1, 2, 3, 4],
            "b": [2, 2, 2, 2],
            "expected": [1, 2, 2, 2]
        },
    ],
    "frompyfunc": [
        {
            "name": "add_one",
            "nin": 2,
            "nout": 1,
            "a": [1, 2, 3],
            "b": [10, 20, 30],
            "expected": [12, 23, 34]  # a + b + 1
        },
    ],
    "piecewise": [
        {
            "name": "abs_piecewise",
            "x": [-2, -1, 0, 1, 2],
            "conditions": ["x < 0", "x >= 0"],
            "functions": ["-x", "x"],
            "expected": [2, 1, 0, 1, 2]
        },
        {
            "name": "with_default",
            "x": [-1, 0, 1, 2, 3],
            "conditions": ["x < 0"],
            "functions": [-1, 1],  # -1 for x<0, 1 otherwise
            "expected": [-1, 1, 1, 1, 1]
        },
    ],
}

with open("tests/fixtures/level10_vectors.json", "w") as f:
    json.dump(tests, f, indent=2)
```

---

## Dependencies

### From Existing Implementation
- `NDArray` class: element access, shape manipulation, slicing
- `iterators.ts`: `ndindex()` for multi-dimensional iteration
- `broadcast.ts`: `broadcastArrays()`, `broadcastTo()`
- `indexing.ts`: `nonzero()`, `extract()` for piecewise

### For Future Levels
Level 10 does NOT directly depend on ufuncs (Level 3), but will benefit from them:
- `piecewise()` placeholder logical operations will be replaced by proper ufuncs
- `vectorize` could leverage ufunc infrastructure for better performance

---

## Performance Considerations

### Current Approach (TypeScript)
- `applyAlongAxis()`: O(n) iterations where n = product of non-axis dimensions
- `vectorize`: O(n) iterations where n = total elements
- `piecewise`: O(k*n) where k = number of conditions

### Future Optimizations
1. **WASM acceleration**: Move inner loops to C for ~10x speedup
2. **Parallel iteration**: Use Web Workers for large arrays
3. **Lazy evaluation**: Defer computation until result is needed
4. **Signature parsing**: Full generalized ufunc support for `vectorize`

### Memory Management
- Functions should dispose intermediate arrays promptly
- Avoid unnecessary copies by using views where possible
- Consider memory pooling for frequently allocated temporary buffers

---

## API Summary

```typescript
// Apply function along single axis
function applyAlongAxis(
  func1d: (arr: NDArray, ...args: unknown[]) => NDArray | number | Promise<NDArray | number>,
  axis: number,
  arr: NDArray,
  ...args: unknown[]
): Promise<NDArray>;

// Apply function over multiple axes
function applyOverAxes(
  func: (arr: NDArray, axis: number) => NDArray | Promise<NDArray>,
  arr: NDArray,
  axes: number | number[]
): Promise<NDArray>;

// Vectorize a scalar function
class Vectorize {
  constructor(pyfunc: Function, options?: VectorizeOptions);
  call(...args: (NDArray | number)[]): Promise<NDArray | NDArray[]>;
}
function vectorize(pyfunc: Function, options?: VectorizeOptions): Vectorize;

// Create ufunc-like from function
function frompyfunc(
  func: (...args: number[]) => number | number[],
  nin: number,
  nout: number,
  identity?: unknown
): UfuncLike;

// Piecewise function evaluation
function piecewise(
  x: NDArray,
  condlist: NDArray[],
  funclist: (Function | number)[],
  ...args: unknown[]
): Promise<NDArray>;
```

---

## Critical Dependencies for Future Levels

Level 10 enables:
- **Custom reduction functions** via `applyAlongAxis()`
- **Element-wise custom operations** via `vectorize()` and `frompyfunc()`
- **Complex conditional logic** via `piecewise()`

These utilities are foundational for user-defined array operations and extend NumJS beyond the built-in functions.
