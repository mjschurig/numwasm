/**
 * NDArray - NumPy-inspired N-dimensional array for TypeScript/WebAssembly
 *
 * This class provides a JavaScript-friendly interface to the underlying
 * WASM implementation which uses NumPy's pairwise summation algorithm
 * for accurate floating-point reductions.
 *
 * @example
 * ```typescript
 * // Create an array and compute sum
 * const arr = await NDArray.fromArray([1, 2, 3, 4, 5]);
 * console.log(arr.sum()); // 15
 * arr.dispose(); // Free memory when done
 *
 * // Create zeros/ones
 * const zeros = await NDArray.zeros([3, 4]);
 * const ones = await NDArray.ones([5]);
 *
 * // Create range
 * const range = await NDArray.arange(0, 10);
 *
 * // Element access
 * const val = arr.get(0);
 * arr.set(42, 0);
 * ```
 */

import {
  DType,
  DTYPE_SIZES,
  DTYPE_NAMES,
  type NDArrayOptions,
  type WasmModule,
  type Complex,
  type NDArrayFlags,
  NDARRAY_OWNDATA,
  NDARRAY_WRITEABLE,
  NDARRAY_C_CONTIGUOUS,
  NDARRAY_F_CONTIGUOUS,
  NDARRAY_ALIGNED,
} from './types.js';
import { loadWasmModule } from './wasm-loader.js';
import { FlatIterator } from './iterators.js';
import {
  type IndexElement,
  expandEllipsis,
  buildIndexSpecs,
} from './slice.js';

/**
 * Type declaration for FileSystemFileHandle (File System Access API).
 * Not always available in all environments.
 */
interface FileSystemFileHandle {
  createWritable(): Promise<FileSystemWritableFileStream>;
}

interface FileSystemWritableFileStream {
  write(data: string | ArrayBuffer | Uint8Array): Promise<void>;
  close(): Promise<void>;
}

/**
 * Type for nested array input (supports arbitrary nesting depth).
 */
type NestedArray = number | NestedArray[];

/**
 * Infer shape from a nested array structure.
 */
function inferShape(data: NestedArray): number[] {
  const shape: number[] = [];
  let current: NestedArray = data;

  while (Array.isArray(current)) {
    shape.push(current.length);
    if (current.length === 0) break;
    current = current[0];
  }

  return shape;
}

/**
 * Flatten a nested array into a 1D array.
 */
function flattenArray(data: NestedArray): number[] {
  if (typeof data === 'number') {
    return [data];
  }

  const result: number[] = [];
  for (const item of data) {
    if (typeof item === 'number') {
      result.push(item);
    } else {
      result.push(...flattenArray(item));
    }
  }
  return result;
}

/**
 * Check if input is a nested array (at least 2D).
 */
function isNestedArray(data: unknown): data is NestedArray[] {
  return (
    Array.isArray(data) &&
    data.length > 0 &&
    Array.isArray(data[0])
  );
}

export class NDArray {
  private readonly _ptr: number;
  private readonly _module: WasmModule;
  private _disposed: boolean = false;

  /**
   * Private constructor - use static factory methods instead.
   */
  private constructor(ptr: number, module: WasmModule) {
    this._ptr = ptr;
    this._module = module;
  }

  /**
   * Create an NDArray from an existing WASM pointer.
   * Used internally for creating arrays from C functions.
   * @internal
   */
  static _fromPtr(ptr: number, module: WasmModule): NDArray {
    if (ptr === 0) {
      throw new Error('Cannot create NDArray from null pointer');
    }
    return new NDArray(ptr, module);
  }

  /**
   * Get the internal WASM pointer.
   * @internal
   */
  get _wasmPtr(): number {
    this.ensureNotDisposed();
    return this._ptr;
  }

  /**
   * Get the internal WASM module.
   * @internal
   */
  get _wasmModule(): WasmModule {
    return this._module;
  }

  /* ============ Static Factory Methods ============ */

  /**
   * Create a new NDArray filled with zeros.
   *
   * @param shape - Array dimensions (e.g., [3, 4] for 3x4 matrix)
   * @param options - Optional configuration (dtype)
   * @returns Promise resolving to the new NDArray
   *
   * @example
   * ```typescript
   * const arr = await NDArray.zeros([3, 4]); // 3x4 matrix of zeros
   * const vec = await NDArray.zeros([5], { dtype: DType.Float32 });
   * ```
   */
  static async zeros(
    shape: number[],
    options: NDArrayOptions = {}
  ): Promise<NDArray> {
    const module = await loadWasmModule();
    const dtype = options.dtype ?? DType.Float64;

    // Validate shape
    if (shape.length === 0) {
      throw new Error('Shape must have at least one dimension');
    }
    for (const dim of shape) {
      if (!Number.isInteger(dim) || dim < 0) {
        throw new Error(`Invalid dimension: ${dim}`);
      }
    }

    // Allocate shape array in WASM memory
    const shapePtr = module._malloc(shape.length * 4);
    for (let i = 0; i < shape.length; i++) {
      module.setValue(shapePtr + i * 4, shape[i], 'i32');
    }

    // Create the array
    const ptr = module._ndarray_create(shape.length, shapePtr, dtype);
    module._free(shapePtr);

    if (ptr === 0) {
      throw new Error('Failed to create NDArray: memory allocation failed');
    }

    return new NDArray(ptr, module);
  }

  /**
   * Create a new NDArray filled with ones.
   *
   * @param shape - Array dimensions
   * @param options - Optional configuration (dtype)
   * @returns Promise resolving to the new NDArray
   */
  static async ones(
    shape: number[],
    options: NDArrayOptions = {}
  ): Promise<NDArray> {
    const arr = await NDArray.zeros(shape, options);
    arr._module._ndarray_fill(arr._ptr, 1.0);
    return arr;
  }

  /**
   * Create a new NDArray without initializing values.
   * Faster than zeros() but contains arbitrary values.
   *
   * @param shape - Array dimensions
   * @param options - Optional configuration (dtype)
   * @returns Promise resolving to the new NDArray
   */
  static async empty(
    shape: number[],
    options: NDArrayOptions = {}
  ): Promise<NDArray> {
    const module = await loadWasmModule();
    const dtype = options.dtype ?? DType.Float64;

    if (shape.length === 0) {
      throw new Error('Shape must have at least one dimension');
    }
    for (const dim of shape) {
      if (!Number.isInteger(dim) || dim < 0) {
        throw new Error(`Invalid dimension: ${dim}`);
      }
    }

    const shapePtr = module._malloc(shape.length * 4);
    for (let i = 0; i < shape.length; i++) {
      module.setValue(shapePtr + i * 4, shape[i], 'i32');
    }

    const ptr = module._ndarray_empty(shape.length, shapePtr, dtype);
    module._free(shapePtr);

    if (ptr === 0) {
      throw new Error('Failed to create NDArray: memory allocation failed');
    }

    return new NDArray(ptr, module);
  }

  /**
   * Create a new NDArray filled with a constant value.
   *
   * @param shape - Array dimensions
   * @param fillValue - Value to fill with
   * @param options - Optional configuration (dtype)
   * @returns Promise resolving to the new NDArray
   */
  static async full(
    shape: number[],
    fillValue: number,
    options: NDArrayOptions = {}
  ): Promise<NDArray> {
    const module = await loadWasmModule();
    const dtype = options.dtype ?? DType.Float64;

    if (shape.length === 0) {
      throw new Error('Shape must have at least one dimension');
    }
    for (const dim of shape) {
      if (!Number.isInteger(dim) || dim < 0) {
        throw new Error(`Invalid dimension: ${dim}`);
      }
    }

    const shapePtr = module._malloc(shape.length * 4);
    for (let i = 0; i < shape.length; i++) {
      module.setValue(shapePtr + i * 4, shape[i], 'i32');
    }

    const ptr = module._ndarray_full(shape.length, shapePtr, dtype, fillValue);
    module._free(shapePtr);

    if (ptr === 0) {
      throw new Error('Failed to create NDArray: memory allocation failed');
    }

    return new NDArray(ptr, module);
  }

  /**
   * Create an NDArray from existing JavaScript data.
   *
   * @param data - Array of numbers, nested arrays, or TypedArray
   * @param shape - Optional shape (defaults to inferred from data)
   * @param options - Optional configuration (dtype)
   * @returns Promise resolving to the new NDArray
   *
   * @example
   * ```typescript
   * const arr = await NDArray.fromArray([1, 2, 3, 4, 5]);
   * const matrix = await NDArray.fromArray([[1, 2, 3], [4, 5, 6]]);
   * const reshaped = await NDArray.fromArray([1, 2, 3, 4, 5, 6], [2, 3]);
   * ```
   */
  static async fromArray(
    data: NestedArray | Float64Array | Float32Array | Int32Array,
    shape?: number[],
    options: NDArrayOptions = {}
  ): Promise<NDArray> {
    const module = await loadWasmModule();
    const dtype = options.dtype ?? DType.Float64;

    // Handle nested arrays
    let flatData: number[];
    let inferredShape: number[];

    if (data instanceof Float64Array || data instanceof Float32Array || data instanceof Int32Array) {
      flatData = Array.from(data as ArrayLike<number>);
      inferredShape = shape ?? [flatData.length];
    } else if (isNestedArray(data)) {
      // Nested array - infer shape and flatten
      inferredShape = shape ?? inferShape(data);
      flatData = flattenArray(data);
    } else if (Array.isArray(data)) {
      // Flat number array
      flatData = data as number[];
      inferredShape = shape ?? [flatData.length];
    } else {
      // Single number (scalar)
      flatData = [data as number];
      inferredShape = shape ?? [1]; // Default to 1-element 1D array
    }

    // Handle 0-dimensional arrays
    const expectedSize = inferredShape.length === 0 ? 1 : inferredShape.reduce((a, b) => a * b, 1);

    if (flatData.length !== expectedSize) {
      throw new Error(
        `Data length ${flatData.length} does not match shape [${inferredShape.join(', ')}] (expected ${expectedSize} elements)`
      );
    }

    // Allocate data buffer in WASM memory
    const bytesPerElement = DTYPE_SIZES[dtype];
    const dataPtr = module._malloc(flatData.length * bytesPerElement);

    // Copy data to WASM memory
    for (let i = 0; i < flatData.length; i++) {
      const offset = dataPtr + i * bytesPerElement;
      switch (dtype) {
        case DType.Float64:
          module.HEAPF64[offset / 8] = flatData[i];
          break;
        case DType.Float32:
          module.HEAPF32[offset / 4] = flatData[i];
          break;
        case DType.Int32:
          module.HEAP32[offset / 4] = flatData[i];
          break;
        case DType.Int16:
          module.HEAP16[offset / 2] = flatData[i];
          break;
        case DType.Int8:
          module.HEAP8[offset] = flatData[i];
          break;
        case DType.Uint32:
          module.HEAPU32[offset / 4] = flatData[i];
          break;
        case DType.Uint16:
          module.HEAPU16[offset / 2] = flatData[i];
          break;
        case DType.Uint8:
        case DType.Bool:
          module.HEAPU8[offset] = flatData[i];
          break;
        case DType.Int64:
          // Note: Int64 handling would need BigInt support
          module.HEAP32[offset / 4] = flatData[i];
          module.HEAP32[offset / 4 + 1] = flatData[i] < 0 ? -1 : 0;
          break;
        default:
          module.HEAPF64[offset / 8] = flatData[i];
      }
    }

    // Allocate shape array in WASM memory
    const shapePtr = module._malloc(inferredShape.length * 4);
    for (let i = 0; i < inferredShape.length; i++) {
      module.setValue(shapePtr + i * 4, inferredShape[i], 'i32');
    }

    // Create the array
    const ptr = module._ndarray_from_data(
      dataPtr,
      inferredShape.length,
      shapePtr,
      dtype
    );

    // Free temporary buffers
    module._free(dataPtr);
    module._free(shapePtr);

    if (ptr === 0) {
      throw new Error('Failed to create NDArray from data');
    }

    return new NDArray(ptr, module);
  }

  /**
   * Create an NDArray from a TypedArray.
   *
   * @param typedArray - Typed array containing data
   * @param shape - Array shape
   * @param dtype - Data type (inferred from typed array if not provided)
   * @returns Promise resolving to the new NDArray
   *
   * @example
   * ```typescript
   * const floats = new Float32Array([1, 2, 3, 4]);
   * const arr = await NDArray.fromTypedArray(floats, [2, 2], DType.Float32);
   * ```
   */
  static async fromTypedArray(
    typedArray: Float64Array | Float32Array | Int32Array | Int16Array | Int8Array |
      Uint32Array | Uint16Array | Uint8Array | BigInt64Array | BigUint64Array,
    shape: number[],
    dtype?: DType
  ): Promise<NDArray> {
    // Infer dtype from typed array if not provided
    const inferredDtype = dtype ?? NDArray._inferDtypeFromTypedArray(typedArray);

    // Convert BigInt arrays to regular numbers if needed
    let data: number[];
    if (typedArray instanceof BigInt64Array || typedArray instanceof BigUint64Array) {
      data = Array.from(typedArray, (v) => Number(v));
    } else {
      data = Array.from(typedArray as ArrayLike<number>);
    }

    return NDArray.fromArray(data, shape, { dtype: inferredDtype });
  }

  /**
   * Infer DType from a TypedArray constructor.
   */
  private static _inferDtypeFromTypedArray(
    typedArray: Float64Array | Float32Array | Int32Array | Int16Array | Int8Array |
      Uint32Array | Uint16Array | Uint8Array | BigInt64Array | BigUint64Array
  ): DType {
    if (typedArray instanceof Float64Array) return DType.Float64;
    if (typedArray instanceof Float32Array) return DType.Float32;
    if (typedArray instanceof Int32Array) return DType.Int32;
    if (typedArray instanceof Int16Array) return DType.Int16;
    if (typedArray instanceof Int8Array) return DType.Int8;
    if (typedArray instanceof Uint32Array) return DType.Uint32;
    if (typedArray instanceof Uint16Array) return DType.Uint16;
    if (typedArray instanceof Uint8Array) return DType.Uint8;
    if (typedArray instanceof BigInt64Array) return DType.Int64;
    if (typedArray instanceof BigUint64Array) return DType.Uint64;
    return DType.Float64;
  }

  /**
   * Create an NDArray with evenly spaced values.
   *
   * @param start - Start value (or end if end is not provided)
   * @param end - End value (exclusive)
   * @param step - Step between values (default: 1)
   * @returns Promise resolving to the new NDArray
   *
   * @example
   * ```typescript
   * const a = await NDArray.arange(5);      // [0, 1, 2, 3, 4]
   * const b = await NDArray.arange(2, 5);   // [2, 3, 4]
   * const c = await NDArray.arange(0, 10, 2); // [0, 2, 4, 6, 8]
   * ```
   */
  static async arange(
    start: number,
    end?: number,
    step: number = 1
  ): Promise<NDArray> {
    // Handle single argument case: arange(5) means arange(0, 5)
    if (end === undefined) {
      end = start;
      start = 0;
    }

    if (step === 0) {
      throw new Error('Step cannot be zero');
    }

    const data: number[] = [];
    if (step > 0) {
      for (let i = start; i < end; i += step) {
        data.push(i);
      }
    } else {
      for (let i = start; i > end; i += step) {
        data.push(i);
      }
    }

    return NDArray.fromArray(data);
  }

  /**
   * Create an NDArray with evenly spaced numbers over an interval.
   *
   * @param start - Start value
   * @param stop - End value
   * @param num - Number of samples (default: 50)
   * @param endpoint - Include stop value (default: true)
   * @param options - Optional configuration (dtype)
   * @returns Promise resolving to the new NDArray
   */
  static async linspace(
    start: number,
    stop: number,
    num: number = 50,
    endpoint: boolean = true,
    options: NDArrayOptions = {}
  ): Promise<NDArray> {
    if (num < 0) {
      throw new Error('Number of samples must be non-negative');
    }
    if (num === 0) {
      return NDArray.empty([0], options);
    }

    const div = endpoint ? num - 1 : num;
    const step = div > 0 ? (stop - start) / div : 0;

    const data: number[] = [];
    for (let i = 0; i < num; i++) {
      data.push(start + i * step);
    }

    return NDArray.fromArray(data, [num], options);
  }

  /**
   * Create an NDArray with numbers spaced evenly on a log scale.
   *
   * @param start - Start value (power of base)
   * @param stop - End value (power of base)
   * @param num - Number of samples (default: 50)
   * @param endpoint - Include stop value (default: true)
   * @param base - Base of the log scale (default: 10)
   * @param options - Optional configuration (dtype)
   * @returns Promise resolving to the new NDArray
   */
  static async logspace(
    start: number,
    stop: number,
    num: number = 50,
    endpoint: boolean = true,
    base: number = 10,
    options: NDArrayOptions = {}
  ): Promise<NDArray> {
    const linArr = await NDArray.linspace(start, stop, num, endpoint, options);
    const data = linArr.toArray().map((x) => Math.pow(base, x));
    linArr.dispose();
    return NDArray.fromArray(data, [num], options);
  }

  /**
   * Create an NDArray with numbers spaced evenly on a geometric scale.
   *
   * @param start - Start value (must be non-zero)
   * @param stop - End value (must be non-zero)
   * @param num - Number of samples (default: 50)
   * @param endpoint - Include stop value (default: true)
   * @param options - Optional configuration (dtype)
   * @returns Promise resolving to the new NDArray
   */
  static async geomspace(
    start: number,
    stop: number,
    num: number = 50,
    endpoint: boolean = true,
    options: NDArrayOptions = {}
  ): Promise<NDArray> {
    if (start === 0 || stop === 0) {
      throw new Error('Geometric sequence cannot include zero');
    }
    if ((start < 0) !== (stop < 0)) {
      throw new Error('Geometric sequence requires same sign for start/stop');
    }

    const logStart = Math.log10(Math.abs(start));
    const logStop = Math.log10(Math.abs(stop));
    const sign = start < 0 ? -1 : 1;

    const linArr = await NDArray.linspace(
      logStart,
      logStop,
      num,
      endpoint,
      options
    );
    const data = linArr.toArray().map((x) => sign * Math.pow(10, x));
    linArr.dispose();
    return NDArray.fromArray(data, [num], options);
  }

  /**
   * Create a 2D array with ones on the diagonal and zeros elsewhere.
   *
   * @param N - Number of rows
   * @param M - Number of columns (default: N)
   * @param k - Diagonal offset (default: 0, main diagonal)
   * @param options - Optional configuration (dtype)
   * @returns Promise resolving to the new NDArray
   */
  static async eye(
    N: number,
    M?: number,
    k: number = 0,
    options: NDArrayOptions = {}
  ): Promise<NDArray> {
    M = M ?? N;
    const arr = await NDArray.zeros([N, M], options);

    // Fill diagonal
    for (let i = 0; i < N; i++) {
      const j = i + k;
      if (j >= 0 && j < M) {
        arr.set(1, i, j);
      }
    }

    return arr;
  }

  /**
   * Create a square identity matrix.
   *
   * @param n - Size of the matrix
   * @param options - Optional configuration (dtype)
   * @returns Promise resolving to the new NDArray
   */
  static async identity(
    n: number,
    options: NDArrayOptions = {}
  ): Promise<NDArray> {
    return NDArray.eye(n, n, 0, options);
  }

  /**
   * Extract a diagonal or construct a diagonal array.
   *
   * @param v - Input array (1D to create diagonal matrix, 2D to extract diagonal)
   * @param k - Diagonal offset (default: 0)
   * @returns Promise resolving to the new NDArray
   */
  static async diag(v: NDArray, k: number = 0): Promise<NDArray> {
    if (v.ndim === 1) {
      // Create diagonal matrix from 1D array
      const n = v.size + Math.abs(k);
      const arr = await NDArray.zeros([n, n], { dtype: v.dtype });
      for (let i = 0; i < v.size; i++) {
        const row = k >= 0 ? i : i - k;
        const col = k >= 0 ? i + k : i;
        arr.set(v.get(i), row, col);
      }
      return arr;
    } else if (v.ndim === 2) {
      // Extract diagonal from 2D array
      const [rows, cols] = v.shape;
      const startRow = k >= 0 ? 0 : -k;
      const startCol = k >= 0 ? k : 0;
      const diagLen = Math.min(rows - startRow, cols - startCol);

      if (diagLen <= 0) {
        return NDArray.empty([0], { dtype: v.dtype });
      }

      const data: number[] = [];
      for (let i = 0; i < diagLen; i++) {
        data.push(v.get(startRow + i, startCol + i));
      }
      return NDArray.fromArray(data, [diagLen], { dtype: v.dtype });
    } else {
      throw new Error('Input must be 1D or 2D array');
    }
  }

  /**
   * Create array with ones at and below diagonal, zeros elsewhere.
   *
   * @param N - Number of rows
   * @param M - Number of columns (default: N)
   * @param k - Diagonal offset (default: 0)
   * @param options - Optional configuration (dtype)
   * @returns Promise resolving to the new NDArray
   */
  static async tri(
    N: number,
    M?: number,
    k: number = 0,
    options: NDArrayOptions = {}
  ): Promise<NDArray> {
    M = M ?? N;
    const arr = await NDArray.zeros([N, M], options);

    for (let i = 0; i < N; i++) {
      for (let j = 0; j <= i + k && j < M; j++) {
        arr.set(1, i, j);
      }
    }

    return arr;
  }

  /**
   * Lower triangle of an array.
   *
   * @param arr - Input array
   * @param k - Diagonal offset (default: 0)
   * @returns Promise resolving to the new NDArray
   */
  static async tril(arr: NDArray, k: number = 0): Promise<NDArray> {
    if (arr.ndim !== 2) {
      throw new Error('Input must be 2D array');
    }

    const copy = arr.copy();
    const [rows, cols] = copy.shape;

    for (let i = 0; i < rows; i++) {
      for (let j = i + k + 1; j < cols; j++) {
        copy.set(0, i, j);
      }
    }

    return copy;
  }

  /**
   * Upper triangle of an array.
   *
   * @param arr - Input array
   * @param k - Diagonal offset (default: 0)
   * @returns Promise resolving to the new NDArray
   */
  static async triu(arr: NDArray, k: number = 0): Promise<NDArray> {
    if (arr.ndim !== 2) {
      throw new Error('Input must be 2D array');
    }

    const copy = arr.copy();
    const [rows, cols] = copy.shape;

    for (let i = 0; i < rows; i++) {
      for (let j = 0; j < Math.min(i + k, cols); j++) {
        copy.set(0, i, j);
      }
    }

    return copy;
  }

  /**
   * Create array with same shape and dtype as another array, filled with zeros.
   */
  static async zerosLike(
    arr: NDArray,
    options: NDArrayOptions = {}
  ): Promise<NDArray> {
    return NDArray.zeros(arr.shape, { dtype: options.dtype ?? arr.dtype });
  }

  /**
   * Create array with same shape and dtype as another array, filled with ones.
   */
  static async onesLike(
    arr: NDArray,
    options: NDArrayOptions = {}
  ): Promise<NDArray> {
    return NDArray.ones(arr.shape, { dtype: options.dtype ?? arr.dtype });
  }

  /**
   * Create uninitialized array with same shape and dtype as another array.
   */
  static async emptyLike(
    arr: NDArray,
    options: NDArrayOptions = {}
  ): Promise<NDArray> {
    return NDArray.empty(arr.shape, { dtype: options.dtype ?? arr.dtype });
  }

  /**
   * Create array with same shape as another array, filled with a value.
   */
  static async fullLike(
    arr: NDArray,
    fillValue: number,
    options: NDArrayOptions = {}
  ): Promise<NDArray> {
    return NDArray.full(arr.shape, fillValue, {
      dtype: options.dtype ?? arr.dtype,
    });
  }

  /* ============ Properties ============ */

  /**
   * Get the shape of the array.
   *
   * @returns Array of dimension sizes
   */
  get shape(): number[] {
    this.ensureNotDisposed();
    const ndim = this._module._ndarray_get_ndim(this._ptr);
    const shapePtr = this._module._ndarray_get_shape(this._ptr);

    const shape: number[] = [];
    for (let i = 0; i < ndim; i++) {
      shape.push(this._module.getValue(shapePtr + i * 4, 'i32'));
    }
    return shape;
  }

  /**
   * Get the number of dimensions.
   */
  get ndim(): number {
    this.ensureNotDisposed();
    return this._module._ndarray_get_ndim(this._ptr);
  }

  /**
   * Get the total number of elements.
   */
  get size(): number {
    this.ensureNotDisposed();
    return this._module._ndarray_get_size(this._ptr);
  }

  /**
   * Get the data type.
   */
  get dtype(): DType {
    this.ensureNotDisposed();
    return this._module._ndarray_get_dtype(this._ptr) as DType;
  }

  /**
   * Get the data type name as a string.
   */
  get dtypeName(): string {
    return DTYPE_NAMES[this.dtype] ?? 'unknown';
  }

  /**
   * Get the strides array (bytes to step in each dimension).
   */
  get strides(): number[] {
    this.ensureNotDisposed();
    const ndim = this._module._ndarray_get_ndim(this._ptr);
    const stridesPtr = this._module._ndarray_get_strides(this._ptr);

    const strides: number[] = [];
    for (let i = 0; i < ndim; i++) {
      strides.push(this._module.getValue(stridesPtr + i * 4, 'i32'));
    }
    return strides;
  }

  /**
   * Get the total bytes consumed by the array data.
   */
  get nbytes(): number {
    this.ensureNotDisposed();
    return this.size * DTYPE_SIZES[this.dtype];
  }

  /**
   * Get the size of one element in bytes.
   */
  get itemsize(): number {
    return DTYPE_SIZES[this.dtype];
  }

  /**
   * Get array flags.
   */
  get flags(): NDArrayFlags {
    this.ensureNotDisposed();
    const flagsInt = this._module._ndarray_get_flags(this._ptr);
    return {
      owndata: (flagsInt & NDARRAY_OWNDATA) !== 0,
      writeable: (flagsInt & NDARRAY_WRITEABLE) !== 0,
      c_contiguous: (flagsInt & NDARRAY_C_CONTIGUOUS) !== 0,
      f_contiguous: (flagsInt & NDARRAY_F_CONTIGUOUS) !== 0,
      aligned: (flagsInt & NDARRAY_ALIGNED) !== 0,
    };
  }

  /**
   * Transpose of the array (reverses axes).
   * Returns a view that shares data with the original array.
   */
  get T(): NDArray {
    this.ensureNotDisposed();

    const resultPtr = this._module._ndarray_transpose(this._ptr, 0);
    if (resultPtr === 0) {
      throw new Error('Transpose failed');
    }

    return new NDArray(resultPtr, this._module);
  }

  /**
   * Whether this array is a view (shares data with another array).
   */
  get isView(): boolean {
    this.ensureNotDisposed();
    return this._module._ndarray_get_base(this._ptr) !== 0;
  }

  /**
   * Flat iterator over array elements in row-major order.
   */
  get flat(): FlatIterator {
    this.ensureNotDisposed();
    return new FlatIterator(this);
  }

  /* ============ Element Access ============ */

  /**
   * Get element at specified indices.
   *
   * @param indices - Multi-dimensional indices
   * @returns Element value
   */
  get(...indices: number[]): number {
    this.ensureNotDisposed();

    if (indices.length !== this.ndim) {
      throw new Error(`Expected ${this.ndim} indices, got ${indices.length}`);
    }

    // Allocate indices array in WASM memory
    const indicesPtr = this._module._malloc(indices.length * 4);
    for (let i = 0; i < indices.length; i++) {
      this._module.setValue(indicesPtr + i * 4, indices[i], 'i32');
    }

    // Check bounds first
    if (!this._module._ndarray_check_bounds(this._ptr, indicesPtr, indices.length)) {
      this._module._free(indicesPtr);
      throw new Error(`Index out of bounds: [${indices.join(', ')}]`);
    }

    const value = this._module._ndarray_get_item(
      this._ptr,
      indicesPtr,
      indices.length
    );
    this._module._free(indicesPtr);

    return value;
  }

  /**
   * Set element at specified indices.
   *
   * @param value - Value to set
   * @param indices - Multi-dimensional indices
   */
  set(value: number, ...indices: number[]): void {
    this.ensureNotDisposed();

    if (indices.length !== this.ndim) {
      throw new Error(`Expected ${this.ndim} indices, got ${indices.length}`);
    }

    const indicesPtr = this._module._malloc(indices.length * 4);
    for (let i = 0; i < indices.length; i++) {
      this._module.setValue(indicesPtr + i * 4, indices[i], 'i32');
    }

    // Check bounds first
    if (!this._module._ndarray_check_bounds(this._ptr, indicesPtr, indices.length)) {
      this._module._free(indicesPtr);
      throw new Error(`Index out of bounds: [${indices.join(', ')}]`);
    }

    this._module._ndarray_set_item(this._ptr, indicesPtr, indices.length, value);
    this._module._free(indicesPtr);
  }

  /**
   * Get element at flat index.
   */
  getFlat(index: number): number {
    this.ensureNotDisposed();
    if (index < 0 || index >= this.size) {
      throw new Error(`Flat index ${index} out of bounds for size ${this.size}`);
    }
    return this._module._ndarray_get_flat(this._ptr, index);
  }

  /**
   * Set element at flat index.
   */
  setFlat(index: number, value: number): void {
    this.ensureNotDisposed();
    if (index < 0 || index >= this.size) {
      throw new Error(`Flat index ${index} out of bounds for size ${this.size}`);
    }
    this._module._ndarray_set_flat(this._ptr, index, value);
  }

  /**
   * Get complex element at indices.
   */
  getComplex(...indices: number[]): Complex {
    this.ensureNotDisposed();

    if (this.dtype !== DType.Complex64 && this.dtype !== DType.Complex128) {
      throw new Error(`getComplex() requires a complex dtype, got ${this.dtypeName}`);
    }

    if (indices.length !== this.ndim) {
      throw new Error(`Expected ${this.ndim} indices, got ${indices.length}`);
    }

    const indicesPtr = this._module._malloc(indices.length * 4);
    for (let i = 0; i < indices.length; i++) {
      this._module.setValue(indicesPtr + i * 4, indices[i], 'i32');
    }

    const flatIdx = this._module._ndarray_flat_index(
      this._ptr,
      indicesPtr,
      indices.length
    );
    this._module._free(indicesPtr);

    if (flatIdx === 0xffffffff) {
      // SIZE_MAX
      throw new Error(`Index out of bounds: [${indices.join(', ')}]`);
    }

    return {
      real: this._module._ndarray_get_complex_real(this._ptr, flatIdx),
      imag: this._module._ndarray_get_complex_imag(this._ptr, flatIdx),
    };
  }

  /**
   * Set complex element at indices.
   */
  setComplex(real: number, imag: number, ...indices: number[]): void {
    this.ensureNotDisposed();

    if (this.dtype !== DType.Complex64 && this.dtype !== DType.Complex128) {
      throw new Error(`setComplex() requires a complex dtype, got ${this.dtypeName}`);
    }

    if (indices.length !== this.ndim) {
      throw new Error(`Expected ${this.ndim} indices, got ${indices.length}`);
    }

    const indicesPtr = this._module._malloc(indices.length * 4);
    for (let i = 0; i < indices.length; i++) {
      this._module.setValue(indicesPtr + i * 4, indices[i], 'i32');
    }

    const flatIdx = this._module._ndarray_flat_index(
      this._ptr,
      indicesPtr,
      indices.length
    );
    this._module._free(indicesPtr);

    if (flatIdx === 0xffffffff) {
      throw new Error(`Index out of bounds: [${indices.join(', ')}]`);
    }

    this._module._ndarray_set_complex(this._ptr, flatIdx, real, imag);
  }

  /**
   * Get single element for scalar or single-element array.
   */
  item(): number {
    this.ensureNotDisposed();
    if (this.size !== 1) {
      throw new Error('item() only works for single-element arrays');
    }
    return this._module._ndarray_get_flat(this._ptr, 0);
  }

  /**
   * Set single element for scalar or single-element array.
   */
  itemset(value: number): void {
    this.ensureNotDisposed();
    if (this.size !== 1) {
      throw new Error('itemset() only works for single-element arrays');
    }
    this._module._ndarray_set_flat(this._ptr, 0, value);
  }

  /* ============ Operations ============ */

  /**
   * Compute the sum of all elements.
   *
   * Uses NumPy's pairwise summation algorithm for O(lg n) rounding error
   * instead of O(n) for naive summation.
   *
   * @returns Sum of all elements
   */
  sum(): number {
    this.ensureNotDisposed();
    return this._module._ndarray_sum(this._ptr);
  }

  /**
   * Fill the array with a constant value.
   *
   * @param value - Value to fill with
   */
  fill(value: number): void {
    this.ensureNotDisposed();
    this._module._ndarray_fill(this._ptr, value);
  }

  /**
   * Get array data as a JavaScript array.
   *
   * @returns Flat array of values
   */
  toArray(): number[] {
    this.ensureNotDisposed();

    const size = this.size;
    const result: number[] = [];

    for (let i = 0; i < size; i++) {
      result.push(this._module._ndarray_get_flat(this._ptr, i));
    }

    return result;
  }

  /**
   * Get array data as a TypedArray.
   */
  toTypedArray():
    | Float64Array
    | Float32Array
    | Int32Array
    | Int16Array
    | Int8Array
    | Uint32Array
    | Uint16Array
    | Uint8Array {
    this.ensureNotDisposed();

    const data = this.toArray();
    const dtype = this.dtype;

    switch (dtype) {
      case DType.Float64:
        return new Float64Array(data);
      case DType.Float32:
        return new Float32Array(data);
      case DType.Int32:
        return new Int32Array(data);
      case DType.Int16:
        return new Int16Array(data);
      case DType.Int8:
        return new Int8Array(data);
      case DType.Uint32:
        return new Uint32Array(data);
      case DType.Uint16:
        return new Uint16Array(data);
      case DType.Uint8:
      case DType.Bool:
        return new Uint8Array(data);
      default:
        return new Float64Array(data);
    }
  }

  /**
   * Write array to a file as binary or text.
   *
   * @param file - File path (Node.js) or FileSystemFileHandle (browser)
   * @param options - Output options
   *
   * @example
   * // Binary format (default)
   * await arr.tofile('data.bin');
   *
   * // Text format with separator
   * await arr.tofile('data.txt', { sep: ',' });
   *
   * // With format string
   * await arr.tofile('data.txt', { sep: ' ', format: '%.4f' });
   */
  async tofile(
    file: string | FileSystemFileHandle,
    options: {
      /** Separator string (empty for binary) */
      sep?: string;
      /** Format string for text output */
      format?: string;
    } = {}
  ): Promise<void> {
    this.ensureNotDisposed();

    const { sep = '', format = '' } = options;

    // Text mode if separator is provided
    if (sep !== '') {
      await this._tofileText(file, sep, format);
    } else {
      await this._tofileBinary(file);
    }
  }

  /**
   * Write array as binary data.
   */
  private async _tofileBinary(file: string | FileSystemFileHandle): Promise<void> {
    const typedArray = this.toTypedArray();
    const bytes = new Uint8Array(typedArray.buffer, typedArray.byteOffset, typedArray.byteLength);

    if (typeof file === 'string') {
      // Node.js file path
      if (typeof process !== 'undefined' && process.versions?.node !== undefined) {
        const fs = await import('fs/promises');
        await fs.writeFile(file, bytes);
        return;
      }
      throw new Error('String file paths only supported in Node.js');
    }

    // FileSystemFileHandle (File System Access API)
    if ('createWritable' in file) {
      const writable = await (file as FileSystemFileHandle).createWritable();
      await writable.write(bytes);
      await writable.close();
      return;
    }

    throw new Error('Unsupported file target');
  }

  /**
   * Write array as text data with separator.
   */
  private async _tofileText(
    file: string | FileSystemFileHandle,
    sep: string,
    format: string
  ): Promise<void> {
    const data = this.toArray();

    // Format values
    let text: string;
    if (format) {
      text = data.map(v => this._formatValue(v, format)).join(sep);
    } else {
      text = data.join(sep);
    }

    if (typeof file === 'string') {
      // Node.js file path
      if (typeof process !== 'undefined' && process.versions?.node !== undefined) {
        const fs = await import('fs/promises');
        await fs.writeFile(file, text, { encoding: 'utf-8' });
        return;
      }
      throw new Error('String file paths only supported in Node.js');
    }

    // FileSystemFileHandle (File System Access API)
    if ('createWritable' in file) {
      const writable = await (file as FileSystemFileHandle).createWritable();
      await writable.write(text);
      await writable.close();
      return;
    }

    throw new Error('Unsupported file target');
  }

  /**
   * Format a value using printf-style format string.
   */
  private _formatValue(value: number, fmt: string): string {
    const match = fmt.match(/^%([+\- 0#]*)(\d*)(?:\.(\d+))?([diouxXeEfFgGaAcs%])$/);
    if (!match) {
      return value.toString();
    }

    const [, flags, width, precision, specifier] = match;
    const prec = precision ? parseInt(precision, 10) : 6;

    let result: string;
    switch (specifier) {
      case 'd':
      case 'i':
        result = Math.trunc(value).toString();
        break;
      case 'e':
        result = value.toExponential(prec);
        break;
      case 'E':
        result = value.toExponential(prec).toUpperCase();
        break;
      case 'f':
      case 'F':
        result = value.toFixed(prec);
        break;
      case 'g':
        result = value.toPrecision(prec || 1);
        break;
      case 'G':
        result = value.toPrecision(prec || 1).toUpperCase();
        break;
      default:
        result = value.toString();
    }

    // Apply sign
    if (flags.includes('+') && value >= 0 && !result.startsWith('-')) {
      result = '+' + result;
    } else if (flags.includes(' ') && value >= 0 && !result.startsWith('-')) {
      result = ' ' + result;
    }

    // Apply width padding
    if (width) {
      const w = parseInt(width, 10);
      const pad = flags.includes('-') ? 'end' : 'start';
      const char = flags.includes('0') && !flags.includes('-') ? '0' : ' ';
      result = pad === 'start' ? result.padStart(w, char) : result.padEnd(w, char);
    }

    return result;
  }

  /* ============ Shape Manipulation ============ */

  /**
   * Returns an array with the given shape.
   * Returns a view if the array is contiguous.
   *
   * @param newShape - New shape (can include one -1 for auto-calculation)
   * @returns View with new shape (or copy if non-contiguous)
   */
  reshape(newShape: number[]): NDArray {
    this.ensureNotDisposed();

    const shapePtr = this._module._malloc(newShape.length * 4);
    for (let i = 0; i < newShape.length; i++) {
      this._module.setValue(shapePtr + i * 4, newShape[i], 'i32');
    }

    const resultPtr = this._module._ndarray_reshape(
      this._ptr,
      newShape.length,
      shapePtr
    );
    this._module._free(shapePtr);

    if (resultPtr === 0) {
      throw new Error(
        'Cannot reshape: incompatible shape or non-contiguous array'
      );
    }

    return new NDArray(resultPtr, this._module);
  }

  /**
   * Transpose array with optional custom axes permutation.
   * Always returns a view.
   *
   * @param axes - Optional axes permutation (default: reverse)
   * @returns Transposed view
   */
  transpose(axes?: number[]): NDArray {
    this.ensureNotDisposed();

    let axesPtr = 0;
    if (axes) {
      if (axes.length !== this.ndim) {
        throw new Error(
          `Axes must have ${this.ndim} elements, got ${axes.length}`
        );
      }
      axesPtr = this._module._malloc(axes.length * 4);
      for (let i = 0; i < axes.length; i++) {
        this._module.setValue(axesPtr + i * 4, axes[i], 'i32');
      }
    }

    const resultPtr = this._module._ndarray_transpose(this._ptr, axesPtr);

    if (axesPtr) this._module._free(axesPtr);

    if (resultPtr === 0) {
      throw new Error('Transpose failed: invalid axes');
    }

    return new NDArray(resultPtr, this._module);
  }

  /**
   * Return a flattened array as a view (if contiguous) or copy.
   */
  ravel(): NDArray {
    this.ensureNotDisposed();

    const resultPtr = this._module._ndarray_ravel(this._ptr);
    if (resultPtr === 0) {
      // Non-contiguous: fall back to flatten (copy)
      return this.flatten();
    }

    return new NDArray(resultPtr, this._module);
  }

  /**
   * Return a flattened array (always a copy).
   */
  flatten(): NDArray {
    this.ensureNotDisposed();

    const resultPtr = this._module._ndarray_flatten(this._ptr);
    if (resultPtr === 0) {
      throw new Error('Flatten failed');
    }

    return new NDArray(resultPtr, this._module);
  }

  /**
   * Remove size-1 dimensions from the array.
   * Returns a view.
   *
   * @param axis - Specific axis to squeeze, or undefined for all size-1 axes
   */
  squeeze(axis?: number): NDArray {
    this.ensureNotDisposed();

    // Use INT32_MIN (0x80000000 = -2147483648) as sentinel for "squeeze all"
    const ax = axis ?? -2147483648;
    const resultPtr = this._module._ndarray_squeeze(this._ptr, ax);
    if (resultPtr === 0) {
      throw new Error(
        axis !== undefined
          ? `Cannot squeeze axis ${axis}: size is not 1`
          : 'Squeeze failed'
      );
    }

    return new NDArray(resultPtr, this._module);
  }

  /**
   * Add a size-1 dimension at the specified position.
   * Returns a view.
   *
   * @param axis - Position for new axis (supports negative indexing)
   */
  expandDims(axis: number): NDArray {
    this.ensureNotDisposed();

    const resultPtr = this._module._ndarray_expand_dims(this._ptr, axis);
    if (resultPtr === 0) {
      throw new Error(
        `Invalid axis ${axis} for array with ${this.ndim} dimensions`
      );
    }

    return new NDArray(resultPtr, this._module);
  }

  /**
   * Swap two axes of the array.
   * Returns a view.
   */
  swapaxes(axis1: number, axis2: number): NDArray {
    this.ensureNotDisposed();

    const resultPtr = this._module._ndarray_swapaxes(this._ptr, axis1, axis2);
    if (resultPtr === 0) {
      throw new Error(`Invalid axes: ${axis1}, ${axis2}`);
    }

    return new NDArray(resultPtr, this._module);
  }

  /**
   * Move axes of an array to new positions.
   * Other axes remain in their original order.
   * Returns a view.
   *
   * @param source - Original positions of the axes to move (can be negative)
   * @param destination - Destination positions for each axis (can be negative)
   * @returns View with reordered axes
   *
   * @example
   * ```typescript
   * const arr = await NDArray.zeros([3, 4, 5]);
   * const moved = arr.moveaxis(0, -1);  // shape: [4, 5, 3]
   * const moved2 = arr.moveaxis([0, 1], [-1, -2]);  // shape: [5, 4, 3]
   * ```
   */
  moveaxis(source: number | number[], destination: number | number[]): NDArray {
    this.ensureNotDisposed();

    const src = Array.isArray(source) ? source : [source];
    const dst = Array.isArray(destination) ? destination : [destination];

    if (src.length !== dst.length) {
      throw new Error('source and destination must have the same number of elements');
    }

    const ndim = this.ndim;

    // Normalize negative indices
    const normSrc = src.map(s => {
      const normalized = s < 0 ? s + ndim : s;
      if (normalized < 0 || normalized >= ndim) {
        throw new Error(`source axis ${s} is out of bounds for array with ${ndim} dimensions`);
      }
      return normalized;
    });

    const normDst = dst.map(d => {
      const normalized = d < 0 ? d + ndim : d;
      if (normalized < 0 || normalized >= ndim) {
        throw new Error(`destination axis ${d} is out of bounds for array with ${ndim} dimensions`);
      }
      return normalized;
    });

    // Check for duplicate source or destination axes
    if (new Set(normSrc).size !== normSrc.length) {
      throw new Error('repeated axis in source');
    }
    if (new Set(normDst).size !== normDst.length) {
      throw new Error('repeated axis in destination');
    }

    // Build the permutation order
    // 1. Get axes not in source (in original order)
    const remaining: number[] = [];
    for (let i = 0; i < ndim; i++) {
      if (!normSrc.includes(i)) {
        remaining.push(i);
      }
    }

    // 2. Create pairs of (dst, src) and sort by destination
    const pairs = normDst.map((d, i) => ({ dst: d, src: normSrc[i] }));
    pairs.sort((a, b) => a.dst - b.dst);

    // 3. Insert source axes at destination positions
    const order: number[] = new Array(ndim);
    let remIdx = 0;
    let pairIdx = 0;

    for (let i = 0; i < ndim; i++) {
      if (pairIdx < pairs.length && pairs[pairIdx].dst === i) {
        order[i] = pairs[pairIdx].src;
        pairIdx++;
      } else {
        order[i] = remaining[remIdx];
        remIdx++;
      }
    }

    return this.transpose(order);
  }

  /* ============ Slicing & Indexing (Level 2) ============ */

  /**
   * Get a sub-array by integer index along axis 0.
   * Returns a view with ndim-1 dimensions.
   *
   * @param index - Index along first axis (supports negative indices)
   * @returns View into this array
   *
   * @example
   * ```typescript
   * const arr = await NDArray.fromArray([[1, 2, 3], [4, 5, 6]]);
   * const row = arr.at(0);  // [1, 2, 3]
   * const lastRow = arr.at(-1);  // [4, 5, 6]
   * ```
   */
  at(index: number): NDArray {
    this.ensureNotDisposed();

    const resultPtr = this._module._ndarray_get_subarray(this._ptr, index);
    if (resultPtr === 0) {
      throw new Error(`Index ${index} out of bounds for axis 0 with size ${this.shape[0]}`);
    }

    return new NDArray(resultPtr, this._module);
  }

  /**
   * Slice the array using NumPy-style indexing.
   *
   * @param indices - Array of index elements (integers, Slice objects, ellipsis, newaxis)
   * @returns View into this array
   *
   * @example
   * ```typescript
   * import { slice, ellipsis, newaxis } from 'numjs';
   *
   * const arr = await NDArray.fromArray([[1, 2, 3], [4, 5, 6], [7, 8, 9]]);
   *
   * // Basic slicing
   * arr.slice([0]);                    // First row: [1, 2, 3]
   * arr.slice([slice(1, null)]);       // Rows 1 onwards: [[4,5,6], [7,8,9]]
   * arr.slice([slice(null, null, 2)]); // Every other row: [[1,2,3], [7,8,9]]
   *
   * // Multi-dimensional
   * arr.slice([0, 1]);                 // Element at (0,1): 2
   * arr.slice([slice(null), 0]);       // First column: [1, 4, 7]
   *
   * // With ellipsis and newaxis
   * arr.slice([ellipsis, 0]);          // Same as arr[:, 0]
   * arr.slice([newaxis]);              // Add dimension: shape (1, 3, 3)
   * ```
   */
  slice(indices: IndexElement[]): NDArray {
    this.ensureNotDisposed();

    // Expand ellipsis in indices
    const expanded = expandEllipsis(indices, this.ndim);

    // Build IndexSpec array for C interop
    const specs = buildIndexSpecs(expanded, this.shape);

    // Allocate IndexSpec array in WASM memory
    // Each IndexSpec is 5 int32s: type, value, start, stop, step (20 bytes)
    const specSize = 20;
    const specsPtr = this._module._malloc(specs.length * specSize);

    for (let i = 0; i < specs.length; i++) {
      const offset = specsPtr + i * specSize;
      this._module.setValue(offset, specs[i].type, 'i32');
      this._module.setValue(offset + 4, specs[i].value, 'i32');
      this._module.setValue(offset + 8, specs[i].start, 'i32');
      this._module.setValue(offset + 12, specs[i].stop, 'i32');
      this._module.setValue(offset + 16, specs[i].step, 'i32');
    }

    const resultPtr = this._module._ndarray_slice(this._ptr, specsPtr, specs.length);
    this._module._free(specsPtr);

    if (resultPtr === 0) {
      throw new Error('Slicing failed: invalid indices');
    }

    return new NDArray(resultPtr, this._module);
  }

  /**
   * Create a view of this array with different shape/strides.
   *
   * @param shape - New shape for the view
   * @param strides - New strides (in bytes) for the view
   * @returns View sharing data with this array
   */
  view(shape?: number[], strides?: number[]): NDArray {
    this.ensureNotDisposed();

    const newShape = shape ?? this.shape;
    const newStrides = strides ?? this.strides;

    if (newShape.length !== newStrides.length) {
      throw new Error('Shape and strides must have the same length');
    }

    const shapePtr = this._module._malloc(newShape.length * 4);
    const stridesPtr = this._module._malloc(newStrides.length * 4);

    for (let i = 0; i < newShape.length; i++) {
      this._module.setValue(shapePtr + i * 4, newShape[i], 'i32');
      this._module.setValue(stridesPtr + i * 4, newStrides[i], 'i32');
    }

    const resultPtr = this._module._ndarray_view(
      this._ptr,
      newShape.length,
      shapePtr,
      stridesPtr
    );

    this._module._free(shapePtr);
    this._module._free(stridesPtr);

    if (resultPtr === 0) {
      throw new Error('Failed to create view');
    }

    return new NDArray(resultPtr, this._module);
  }

  /**
   * Create a view with a different dtype interpretation.
   * The array must be C-contiguous.
   *
   * @param dtype - New dtype for the view
   * @returns View with different dtype
   */
  viewDtype(dtype: DType): NDArray {
    this.ensureNotDisposed();

    const resultPtr = this._module._ndarray_view_dtype(this._ptr, dtype);
    if (resultPtr === 0) {
      throw new Error('Failed to create dtype view: array must be C-contiguous and sizes must be compatible');
    }

    return new NDArray(resultPtr, this._module);
  }

  /**
   * Return a C-contiguous array.
   * Returns a view if already contiguous, otherwise a copy.
   */
  ascontiguousarray(): NDArray {
    this.ensureNotDisposed();

    const resultPtr = this._module._ndarray_ascontiguousarray(this._ptr);
    if (resultPtr === 0) {
      throw new Error('Failed to create contiguous array');
    }

    return new NDArray(resultPtr, this._module);
  }

  /**
   * Return a Fortran-contiguous (column-major) array.
   * Returns a view if already F-contiguous, otherwise a copy.
   */
  asfortranarray(): NDArray {
    this.ensureNotDisposed();

    const resultPtr = this._module._ndarray_asfortranarray(this._ptr);
    if (resultPtr === 0) {
      throw new Error('Failed to create F-contiguous array');
    }

    return new NDArray(resultPtr, this._module);
  }

  /* ============ Copy & Conversion ============ */

  /**
   * Create a deep copy of the array.
   */
  copy(): NDArray {
    this.ensureNotDisposed();
    const ptr = this._module._ndarray_copy(this._ptr);
    if (ptr === 0) {
      throw new Error('Failed to copy array');
    }
    return new NDArray(ptr, this._module);
  }

  /**
   * Create a copy with a different dtype.
   */
  astype(dtype: DType): NDArray {
    this.ensureNotDisposed();
    const ptr = this._module._ndarray_astype(this._ptr, dtype);
    if (ptr === 0) {
      throw new Error('Failed to convert array type');
    }
    return new NDArray(ptr, this._module);
  }

  /* ============ Memory Management ============ */

  /**
   * Free the underlying WASM memory.
   *
   * After calling dispose(), the array can no longer be used.
   * Always call dispose() when done with an array to prevent memory leaks.
   */
  dispose(): void {
    if (!this._disposed) {
      this._module._ndarray_free(this._ptr);
      this._disposed = true;
    }
  }

  /**
   * Check if this array has been disposed.
   */
  get isDisposed(): boolean {
    return this._disposed;
  }

  /* ============ Private Methods ============ */

  private ensureNotDisposed(): void {
    if (this._disposed) {
      throw new Error(
        'NDArray has been disposed. Cannot perform operations on a disposed array.'
      );
    }
  }
}
