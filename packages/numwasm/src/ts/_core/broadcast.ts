/**
 * NumJS Broadcasting Functions
 *
 * TypeScript wrappers for WASM broadcasting operations.
 * Adapted from NumPy's broadcasting implementation.
 */

import { NDArray } from "./NDArray.js";
import { loadWasmModule } from "./wasm-loader.js";

/**
 * Maximum number of dimensions supported
 */
const NPY_MAXDIMS = 32;

/* ============ Broadcast Class ============ */

/**
 * Broadcast class for iterating over multiple arrays as if they had been broadcast.
 *
 * This class encapsulates the result of broadcasting arrays together
 * and provides an iterable interface to iterate over them.
 *
 * @example
 * ```typescript
 * const x = await NDArray.fromArray([1, 2, 3]);
 * const y = await NDArray.fromArray([[4], [5]]);
 * const b = await broadcast.create(x, y);
 *
 * console.log(b.shape);  // [2, 3]
 * console.log(b.ndim);   // 2
 * console.log(b.size);   // 6
 *
 * // Iterate over broadcast indices
 * for (const idx of b) {
 *   console.log(idx);  // [0,0], [0,1], [0,2], [1,0], [1,1], [1,2]
 * }
 *
 * // Reset and iterate again
 * b.reset();
 * ```
 */
export class broadcast implements Iterable<number[]> {
  private _arrays: NDArray[];
  private _shape: number[];
  private _index: number[];
  private _position: number;
  private _numIterations: number;

  /**
   * Create a broadcast object.
   * Use broadcast.create() instead of constructor directly.
   */
  private constructor(arrays: NDArray[], shape: number[]) {
    this._arrays = arrays;
    this._shape = shape;
    this._index = new Array(shape.length).fill(0);
    this._position = 0;
    this._numIterations = shape.reduce((a, b) => a * b, 1);
  }

  /**
   * Create a broadcast object from multiple arrays.
   *
   * @param arrays - Arrays to broadcast together
   * @returns A new broadcast object
   * @throws Error if arrays cannot be broadcast together
   *
   * @example
   * ```typescript
   * const b = await broadcast.create(arr1, arr2, arr3);
   * ```
   */
  static async create(...arrays: NDArray[]): Promise<broadcast> {
    if (arrays.length === 0) {
      return new broadcast([], []);
    }

    const shapes = arrays.map((arr) => arr.shape);
    const broadcastShape = await broadcastShapesMulti(shapes);

    if (broadcastShape === null) {
      const shapeStrs = shapes.map((s) => `(${s.join(", ")})`).join(", ");
      throw new Error(
        `shape mismatch: objects cannot be broadcast to a single shape. Shapes: ${shapeStrs}`,
      );
    }

    return new broadcast(arrays, broadcastShape);
  }

  /**
   * The shape of the broadcast result.
   */
  get shape(): number[] {
    return [...this._shape];
  }

  /**
   * Number of dimensions of the broadcast result.
   */
  get ndim(): number {
    return this._shape.length;
  }

  /**
   * Total size of the broadcast result.
   */
  get size(): number {
    return this._numIterations;
  }

  /**
   * Number of arrays being broadcast.
   */
  get numiter(): number {
    return this._arrays.length;
  }

  /**
   * Current iteration index (flat index).
   */
  get index(): number {
    return this._position;
  }

  /**
   * The input arrays.
   */
  get iters(): NDArray[] {
    return this._arrays;
  }

  /**
   * Reset the broadcast iterator to the beginning.
   */
  reset(): void {
    this._index = new Array(this._shape.length).fill(0);
    this._position = 0;
  }

  /**
   * Iterate over multi-dimensional indices.
   */
  *[Symbol.iterator](): Iterator<number[]> {
    this.reset();

    while (this._position < this._numIterations) {
      yield [...this._index];

      // Increment multi-dimensional index
      this._position++;
      if (this._position < this._numIterations) {
        for (let d = this._shape.length - 1; d >= 0; d--) {
          this._index[d]++;
          if (this._index[d] < this._shape[d]) {
            break;
          }
          this._index[d] = 0;
        }
      }
    }
  }

  /**
   * Get values from all arrays at the current broadcast position.
   *
   * @param index - Multi-dimensional index in broadcast space
   * @returns Array of values, one from each input array
   */
  getValues(index: number[]): number[] {
    const values: number[] = [];

    for (const arr of this._arrays) {
      // Map broadcast index to array index
      const arrIndex: number[] = [];
      const offset = index.length - arr.ndim;

      for (let d = 0; d < arr.ndim; d++) {
        const broadcastDim = d + offset;
        if (arr.shape[d] === 1) {
          arrIndex.push(0);
        } else {
          arrIndex.push(index[broadcastDim]);
        }
      }

      // Get flat index and retrieve value
      let flatIdx = 0;
      for (let d = 0; d < arr.ndim; d++) {
        flatIdx = flatIdx * arr.shape[d] + arrIndex[d];
      }
      values.push(arr.getFlat(flatIdx));
    }

    return values;
  }
}

/**
 * Compute the broadcast shape for two shapes.
 *
 * @param shape1 - First shape
 * @param shape2 - Second shape
 * @returns Broadcast shape, or null if shapes are incompatible
 *
 * @example
 * ```typescript
 * broadcastShapes([3, 1], [1, 4]);  // [3, 4]
 * broadcastShapes([2, 3], [3]);      // [2, 3]
 * broadcastShapes([2, 3], [4]);      // null (incompatible)
 * ```
 */
export async function broadcastShapes(
  shape1: number[],
  shape2: number[],
): Promise<number[] | null> {
  const module = await loadWasmModule();

  // Allocate memory for inputs and outputs
  const shape1Ptr = module._malloc(shape1.length * 4);
  const shape2Ptr = module._malloc(shape2.length * 4);
  const outShapePtr = module._malloc(NPY_MAXDIMS * 4);
  const outNdimPtr = module._malloc(4);

  // Copy input shapes
  for (let i = 0; i < shape1.length; i++) {
    module.setValue(shape1Ptr + i * 4, shape1[i], "i32");
  }
  for (let i = 0; i < shape2.length; i++) {
    module.setValue(shape2Ptr + i * 4, shape2[i], "i32");
  }

  // Call WASM function
  const result = module._broadcast_shapes(
    shape1Ptr,
    shape1.length,
    shape2Ptr,
    shape2.length,
    outShapePtr,
    outNdimPtr,
  );

  if (result !== 0) {
    // Shapes are incompatible
    module._free(shape1Ptr);
    module._free(shape2Ptr);
    module._free(outShapePtr);
    module._free(outNdimPtr);
    return null;
  }

  // Read result
  const outNdim = module.getValue(outNdimPtr, "i32");
  const outShape: number[] = [];
  for (let i = 0; i < outNdim; i++) {
    outShape.push(module.getValue(outShapePtr + i * 4, "i32"));
  }

  // Free memory
  module._free(shape1Ptr);
  module._free(shape2Ptr);
  module._free(outShapePtr);
  module._free(outNdimPtr);

  return outShape;
}

/**
 * Compute the broadcast shape for multiple shapes.
 *
 * @param shapes - Array of shapes
 * @returns Broadcast shape, or null if shapes are incompatible
 *
 * @example
 * ```typescript
 * broadcastShapesMulti([[3, 1], [1, 4], [1, 1]]);  // [3, 4]
 * ```
 */
export async function broadcastShapesMulti(
  shapes: number[][],
): Promise<number[] | null> {
  if (shapes.length === 0) {
    return [];
  }
  if (shapes.length === 1) {
    return [...shapes[0]];
  }

  // Use pairwise broadcasting
  let result = shapes[0];
  for (let i = 1; i < shapes.length; i++) {
    const newResult = await broadcastShapes(result, shapes[i]);
    if (newResult === null) {
      return null;
    }
    result = newResult;
  }
  return result;
}

/**
 * Check if two shapes are broadcast-compatible.
 *
 * @param shape1 - First shape
 * @param shape2 - Second shape
 * @returns true if shapes can be broadcast together
 */
export async function shapesAreBroadcastable(
  shape1: number[],
  shape2: number[],
): Promise<boolean> {
  const result = await broadcastShapes(shape1, shape2);
  return result !== null;
}

/**
 * Broadcast an array to a target shape.
 *
 * Returns a view of the array with the target shape. The view shares
 * data with the original array but has zero strides for broadcast dimensions.
 *
 * @param arr - Array to broadcast
 * @param targetShape - Target shape
 * @returns Broadcast view (read-only)
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([1, 2, 3]);  // shape [3]
 * const broadcast = await broadcastTo(arr, [2, 3]); // shape [2, 3]
 * // broadcast looks like [[1,2,3], [1,2,3]] but shares data
 * ```
 */
export async function broadcastTo(
  arr: NDArray,
  targetShape: number[],
): Promise<NDArray> {
  const module = arr._wasmModule;

  // Allocate target shape in WASM memory
  const targetShapePtr = module._malloc(targetShape.length * 4);
  for (let i = 0; i < targetShape.length; i++) {
    module.setValue(targetShapePtr + i * 4, targetShape[i], "i32");
  }

  // Call WASM function
  const resultPtr = module._ndarray_broadcast_to(
    arr._wasmPtr,
    targetShapePtr,
    targetShape.length,
  );

  module._free(targetShapePtr);

  if (resultPtr === 0) {
    throw new Error(
      `Cannot broadcast array with shape [${arr.shape.join(", ")}] to [${targetShape.join(", ")}]`,
    );
  }

  return NDArray._fromPtr(resultPtr, module);
}

/**
 * Broadcast any number of arrays against each other.
 *
 * Returns views of the arrays with the same shape. All returned arrays
 * share data with their originals.
 *
 * @param arrays - Arrays to broadcast together
 * @returns Array of broadcast views with identical shapes
 *
 * @example
 * ```typescript
 * const a = await NDArray.fromArray([1, 2, 3]);       // shape [3]
 * const b = await NDArray.fromArray([[1], [2]]);     // shape [2, 1]
 * const [aBc, bBc] = await broadcastArrays(a, b);    // both shape [2, 3]
 * ```
 */
export async function broadcastArrays(
  ...arrays: NDArray[]
): Promise<NDArray[]> {
  if (arrays.length === 0) {
    return [];
  }
  if (arrays.length === 1) {
    // Return a view to maintain consistency
    return [arrays[0].view()];
  }

  // Compute broadcast shape
  const shapes = arrays.map((arr) => arr.shape);
  const targetShape = await broadcastShapesMulti(shapes);

  if (targetShape === null) {
    const shapeStrs = shapes.map((s) => `[${s.join(", ")}]`).join(", ");
    throw new Error(`Cannot broadcast arrays with shapes: ${shapeStrs}`);
  }

  // Broadcast each array to the target shape
  const results: NDArray[] = [];
  for (const arr of arrays) {
    // If already the right shape, just create a view
    if (
      arr.shape.length === targetShape.length &&
      arr.shape.every((d, i) => d === targetShape[i])
    ) {
      results.push(arr.view());
    } else {
      results.push(await broadcastTo(arr, targetShape));
    }
  }

  return results;
}

/**
 * Compute the strides needed to broadcast an array to a target shape.
 * This is mainly for internal use.
 *
 * @param arr - Array to broadcast
 * @param targetShape - Target shape
 * @returns Broadcast strides, or null if not broadcastable
 */
export async function computeBroadcastStrides(
  arr: NDArray,
  targetShape: number[],
): Promise<number[] | null> {
  const module = arr._wasmModule;

  // Allocate memory
  const targetShapePtr = module._malloc(targetShape.length * 4);
  const outStridesPtr = module._malloc(targetShape.length * 4);

  for (let i = 0; i < targetShape.length; i++) {
    module.setValue(targetShapePtr + i * 4, targetShape[i], "i32");
  }

  // Call WASM function
  const result = module._broadcast_strides(
    arr._wasmPtr,
    targetShapePtr,
    targetShape.length,
    outStridesPtr,
  );

  if (result !== 0) {
    module._free(targetShapePtr);
    module._free(outStridesPtr);
    return null;
  }

  // Read result
  const outStrides: number[] = [];
  for (let i = 0; i < targetShape.length; i++) {
    outStrides.push(module.getValue(outStridesPtr + i * 4, "i32"));
  }

  module._free(targetShapePtr);
  module._free(outStridesPtr);

  return outStrides;
}
