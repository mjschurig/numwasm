/**
 * Linear Algebra Matrix Products for NumJS-WASM
 *
 * Matrix multiplication, dot products, and related operations.
 */

import { NDArray } from "../_core/NDArray.js";
import { DType } from "../types.js";
import { getWasmModule } from "../wasm-loader.js";
import { LinAlgError, prepareMatrix, getLinalgDtype } from "./types.js";

/**
 * Matrix multiplication of two arrays.
 *
 * For 2-D arrays: regular matrix multiplication.
 * For N-D arrays: broadcast over batch dimensions, matmul on last two.
 *
 * @param a - First array (..., M, K)
 * @param b - Second array (..., K, N)
 * @returns Product array (..., M, N)
 *
 * @example
 * const a = await NDArray.fromArray([[1, 2], [3, 4]]);
 * const b = await NDArray.fromArray([[5, 6], [7, 8]]);
 * const c = await matmul(a, b);
 * // c = [[19, 22], [43, 50]]
 */
export async function matmul(a: NDArray, b: NDArray): Promise<NDArray> {
  if (a.ndim < 2 || b.ndim < 2) {
    throw new LinAlgError("matmul requires at least 2-dimensional arrays");
  }

  const k1 = a.shape[a.ndim - 1];
  const k2 = b.shape[b.ndim - 2];

  if (k1 !== k2) {
    throw new LinAlgError(
      `matmul: Input operand 1 has a mismatch in its core dimension 0, ` +
        `with gufunc signature (n?,k),(k,m?)->(n?,m?) (size ${k1} is different from ${k2})`,
    );
  }

  // Prepare arrays
  const aPrep = await prepareMatrix(a, "a");
  const bPrep = await prepareMatrix(b, "b");

  const module = getWasmModule();
  const resultPtr = module._linalg_matmul(aPrep._wasmPtr, bPrep._wasmPtr);

  if (resultPtr === 0) {
    throw new LinAlgError("matmul failed - check array shapes and types");
  }

  return NDArray._fromPtr(resultPtr, module);
}

/**
 * Dot product of two arrays.
 *
 * For 1-D arrays: inner product of vectors.
 * For 2-D arrays: matrix multiplication.
 *
 * @param a - First array
 * @param b - Second array
 * @returns Dot product result
 *
 * @example
 * // Vector dot product
 * const v1 = await NDArray.fromArray([1, 2, 3]);
 * const v2 = await NDArray.fromArray([4, 5, 6]);
 * const d = await dot(v1, v2);  // 32
 */
export async function dot(a: NDArray, b: NDArray): Promise<NDArray> {
  if (a.ndim === 0 || b.ndim === 0) {
    throw new LinAlgError("dot does not support 0-d arrays");
  }

  // Ensure appropriate dtype
  const targetDtype = getLinalgDtype(a);
  const aPrep = a.dtype !== targetDtype ? a.astype(targetDtype) : a;
  const bPrep = b.dtype !== targetDtype ? b.astype(targetDtype) : b;

  const module = getWasmModule();
  const resultPtr = module._linalg_dot(aPrep._wasmPtr, bPrep._wasmPtr);

  if (resultPtr === 0) {
    throw new LinAlgError("dot failed - check array shapes");
  }

  return NDArray._fromPtr(resultPtr, module);
}

/**
 * Compute the dot product of two vectors (flattened inputs).
 * For complex arrays, the first argument is conjugated.
 */
export async function vdot(a: NDArray, b: NDArray): Promise<number> {
  const aFlat = a.ravel();
  const bFlat = b.ravel();

  if (aFlat.size !== bFlat.size) {
    throw new LinAlgError("vdot requires arrays of equal size");
  }

  const result = await dot(aFlat, bFlat);
  return result.item() as number;
}

/**
 * Inner product of two arrays.
 * For 1-D arrays: dot product.
 * For N-D arrays: sum product over last axes.
 */
export async function inner(a: NDArray, b: NDArray): Promise<NDArray> {
  if (a.ndim === 1 && b.ndim === 1) {
    return dot(a, b);
  }

  // For N-D, reshape and use matmul
  // inner(a, b) = sum over last axis of a and last axis of b
  const aLast = a.shape[a.ndim - 1];
  const bLast = b.shape[b.ndim - 1];

  if (aLast !== bLast) {
    throw new LinAlgError(
      `shapes ${a.shape} and ${b.shape} not aligned: ` +
        `${aLast} (dim ${a.ndim - 1}) != ${bLast} (dim ${b.ndim - 1})`,
    );
  }

  // Implement via reshape + matmul + reshape
  throw new LinAlgError("inner for N-D arrays not yet implemented");
}

/**
 * Compute the outer product of two vectors.
 */
export async function outer(a: NDArray, b: NDArray): Promise<NDArray> {
  const aFlat = a.ravel();
  const bFlat = b.ravel();

  // Reshape to column and row vectors
  const aCol = aFlat.reshape([aFlat.size, 1]);
  const bRow = bFlat.reshape([1, bFlat.size]);

  return matmul(aCol, bRow);
}

/**
 * Compute tensor dot product along specified axes.
 *
 * For tensors a and b, tensordot(a, b, axes) computes the sum of products
 * over the specified axes.
 *
 * @param a - First tensor
 * @param b - Second tensor
 * @param axes - Axes to sum over:
 *   - number N: last N axes of a with first N axes of b
 *   - [axesA, axesB]: specific axes for each tensor
 * @returns Tensor contraction result
 *
 * @example
 * // Matrix multiplication
 * const c = await tensordot(a, b, 1);
 * // Same as: matmul(a, b)
 *
 * // Outer product
 * const c = await tensordot(a, b, 0);
 * // Shape: [...a.shape, ...b.shape]
 */
export async function tensordot(
  a: NDArray,
  b: NDArray,
  axes: number | [number[], number[]] = 2,
): Promise<NDArray> {
  // Parse axes specification
  let axesA: number[];
  let axesB: number[];

  if (typeof axes === "number") {
    if (axes < 0) {
      throw new LinAlgError("axes must be non-negative");
    }
    // Last N axes of a, first N axes of b
    axesA = [];
    axesB = [];
    for (let i = 0; i < axes; i++) {
      axesA.push(a.ndim - axes + i);
      axesB.push(i);
    }
  } else {
    [axesA, axesB] = axes;
    if (axesA.length !== axesB.length) {
      throw new LinAlgError("axes lists must have same length");
    }
  }

  // Check for duplicate axes
  if (new Set(axesA).size !== axesA.length) {
    throw new LinAlgError("duplicate axes not allowed in tensordot");
  }
  if (new Set(axesB).size !== axesB.length) {
    throw new LinAlgError("duplicate axes not allowed in tensordot");
  }

  // Normalize negative axes
  axesA = axesA.map((ax) => (ax < 0 ? a.ndim + ax : ax));
  axesB = axesB.map((ax) => (ax < 0 ? b.ndim + ax : ax));

  // Validate axes are in bounds
  for (const ax of axesA) {
    if (ax < 0 || ax >= a.ndim) {
      throw new LinAlgError(
        `axis ${ax} is out of bounds for array with ${a.ndim} dimensions`,
      );
    }
  }
  for (const ax of axesB) {
    if (ax < 0 || ax >= b.ndim) {
      throw new LinAlgError(
        `axis ${ax} is out of bounds for array with ${b.ndim} dimensions`,
      );
    }
  }

  // Validate axes match in size
  for (let i = 0; i < axesA.length; i++) {
    if (a.shape[axesA[i]] !== b.shape[axesB[i]]) {
      throw new LinAlgError(
        `shape mismatch for sum: ${a.shape[axesA[i]]} vs ${b.shape[axesB[i]]}`,
      );
    }
  }

  // Compute non-contracted axes (free axes)
  const freeA = a.shape.map((_, i) => i).filter((i) => !axesA.includes(i));
  const freeB = b.shape.map((_, i) => i).filter((i) => !axesB.includes(i));

  // Output shape from free axes
  const oldShapeA = freeA.map((i) => a.shape[i]);
  const oldShapeB = freeB.map((i) => b.shape[i]);

  // Transpose a: free axes first, then contracted
  const newAxesA = [...freeA, ...axesA];
  const aT = a.transpose(newAxesA);

  // Transpose b: contracted axes first, then free
  const newAxesB = [...axesB, ...freeB];
  const bT = b.transpose(newAxesB);

  // Compute reshape dimensions
  const N2 = axesA.reduce((p, i) => p * a.shape[i], 1);
  const freeASize = oldShapeA.reduce((p, x) => p * x, 1);
  const freeBSize = oldShapeB.reduce((p, x) => p * x, 1);

  // Handle edge case of no free axes
  const newShapeA: [number, number] = [freeASize || 1, N2 || 1];
  const newShapeB: [number, number] = [N2 || 1, freeBSize || 1];

  // Reshape to 2D
  const aReshaped = aT.reshape(newShapeA);
  const bReshaped = bT.reshape(newShapeB);

  // Matrix multiply
  const result = await dot(aReshaped, bReshaped);

  // Reshape to final output shape
  const finalShape = [...oldShapeA, ...oldShapeB];
  if (finalShape.length === 0) {
    // Scalar result
    return result.reshape([]);
  }
  return result.reshape(finalShape);
}

/**
 * Compute the dot product of two or more arrays in a single function call,
 * while automatically selecting the fastest evaluation order.
 *
 * Uses dynamic programming to find the optimal parenthesization that
 * minimizes the total number of scalar multiplications.
 *
 * @param arrays - List of arrays to multiply together
 * @returns Product of all arrays
 *
 * @example
 * // For arrays A (10x100), B (100x5), C (5x50)
 * // multi_dot([A, B, C]) is much faster than A @ B @ C
 * // because it computes (A @ B) @ C instead of A @ (B @ C)
 */
export async function multi_dot(arrays: NDArray[]): Promise<NDArray> {
  const n = arrays.length;

  if (n < 2) {
    throw new LinAlgError("multi_dot requires at least 2 arrays");
  }

  if (n === 2) {
    return dot(arrays[0], arrays[1]);
  }

  // Handle 1-D arrays: first as row vector, last as column vector
  const workArrays = [...arrays];
  let firstIs1D = false;
  let lastIs1D = false;

  if (workArrays[0].ndim === 1) {
    firstIs1D = true;
    workArrays[0] = workArrays[0].reshape([1, workArrays[0].shape[0]]);
  }
  if (workArrays[n - 1].ndim === 1) {
    lastIs1D = true;
    workArrays[n - 1] = workArrays[n - 1].reshape([
      workArrays[n - 1].shape[0],
      1,
    ]);
  }

  // Validate all arrays are 2D
  for (let i = 0; i < n; i++) {
    if (workArrays[i].ndim !== 2) {
      throw new LinAlgError(
        `multi_dot requires 1-D or 2-D arrays, got ${workArrays[i].ndim}-D at position ${i}`,
      );
    }
  }

  // Get dimensions for cost computation: p[i] = rows of matrix i, p[n] = cols of last matrix
  const p: number[] = [];
  p.push(workArrays[0].shape[0]);
  for (let i = 0; i < n; i++) {
    p.push(workArrays[i].shape[1]);
  }

  // Validate chain dimensions match
  for (let i = 0; i < n - 1; i++) {
    if (workArrays[i].shape[1] !== workArrays[i + 1].shape[0]) {
      throw new LinAlgError(
        `shapes (${workArrays[i].shape.join(",")}) and (${workArrays[i + 1].shape.join(",")}) not aligned`,
      );
    }
  }

  // For n=3, use simple cost comparison
  if (n === 3) {
    const cost1 = p[0] * p[1] * p[2] + p[0] * p[2] * p[3]; // (AB)C
    const cost2 = p[1] * p[2] * p[3] + p[0] * p[1] * p[3]; // A(BC)

    let result: NDArray;
    if (cost1 <= cost2) {
      const ab = await dot(workArrays[0], workArrays[1]);
      result = await dot(ab, workArrays[2]);
    } else {
      const bc = await dot(workArrays[1], workArrays[2]);
      result = await dot(workArrays[0], bc);
    }

    // Reshape result if input was 1-D
    if (firstIs1D && lastIs1D) {
      return result.reshape([]);
    } else if (firstIs1D) {
      return result.reshape([result.shape[1]]);
    } else if (lastIs1D) {
      return result.reshape([result.shape[0]]);
    }
    return result;
  }

  // Dynamic programming for n > 3: Matrix Chain Multiplication
  // m[i][j] = minimum cost to multiply matrices i..j
  // s[i][j] = optimal split point
  const m: number[][] = Array(n)
    .fill(null)
    .map(() => Array(n).fill(0));
  const s: number[][] = Array(n)
    .fill(null)
    .map(() => Array(n).fill(0));

  // Chain length from 2 to n
  for (let len = 2; len <= n; len++) {
    for (let i = 0; i <= n - len; i++) {
      const j = i + len - 1;
      m[i][j] = Infinity;

      for (let k = i; k < j; k++) {
        const cost = m[i][k] + m[k + 1][j] + p[i] * p[k + 1] * p[j + 1];
        if (cost < m[i][j]) {
          m[i][j] = cost;
          s[i][j] = k;
        }
      }
    }
  }

  // Execute multiplications in optimal order
  async function executeChain(i: number, j: number): Promise<NDArray> {
    if (i === j) {
      return workArrays[i];
    }
    const k = s[i][j];
    const left = await executeChain(i, k);
    const right = await executeChain(k + 1, j);
    return dot(left, right);
  }

  let result = await executeChain(0, n - 1);

  // Reshape result if input was 1-D
  if (firstIs1D && lastIs1D) {
    return result.reshape([]);
  } else if (firstIs1D) {
    return result.reshape([result.shape[1]]);
  } else if (lastIs1D) {
    return result.reshape([result.shape[0]]);
  }
  return result;
}

/**
 * Kronecker product of two arrays.
 *
 * Computes the Kronecker product, a composite array made of blocks of the
 * second array scaled by the first.
 *
 * @param a - First array
 * @param b - Second array
 * @returns Kronecker product
 *
 * @example
 * kron([1, 10, 100], [5, 6, 7])
 * // [5, 6, 7, 50, 60, 70, 500, 600, 700]
 *
 * kron([[1, 2], [3, 4]], [[1, 0], [0, 1]])
 * // [[1, 0, 2, 0],
 * //  [0, 1, 0, 2],
 * //  [3, 0, 4, 0],
 * //  [0, 3, 0, 4]]
 */
export async function kron(a: NDArray, b: NDArray): Promise<NDArray> {
  // Handle scalar case
  if (a.ndim === 0 && b.ndim === 0) {
    const val = (a.item() as number) * (b.item() as number);
    return NDArray.fromArray([val]).then((arr) => arr.reshape([]));
  }

  // Ensure at least 1D
  let aArr = a.ndim === 0 ? a.reshape([1]) : a;
  let bArr = b.ndim === 0 ? b.reshape([1]) : b;

  const nda = aArr.ndim;
  const ndb = bArr.ndim;
  const nd = Math.max(nda, ndb);

  // Equalize dimensions by prepending 1s
  if (nda < nd) {
    const newShape = [...Array(nd - nda).fill(1), ...aArr.shape];
    aArr = aArr.reshape(newShape);
  }
  if (ndb < nd) {
    const newShape = [...Array(nd - ndb).fill(1), ...bArr.shape];
    bArr = bArr.reshape(newShape);
  }

  // Build interleaved shape for expansion
  // a goes to positions 0, 2, 4, ... with b dimensions as 1
  // b goes to positions 1, 3, 5, ... with a dimensions as 1
  const aExpandShape: number[] = [];
  const bExpandShape: number[] = [];

  for (let i = 0; i < nd; i++) {
    aExpandShape.push(aArr.shape[i]);
    aExpandShape.push(1);
    bExpandShape.push(1);
    bExpandShape.push(bArr.shape[i]);
  }

  const aExpanded = aArr.reshape(aExpandShape);
  const bExpanded = bArr.reshape(bExpandShape);

  // Multiply (broadcasting handles it) using ufunc_multiply
  const module = getWasmModule();
  const productPtr = module._ufunc_multiply(
    aExpanded._wasmPtr,
    bExpanded._wasmPtr,
  );
  if (productPtr === 0) {
    throw new LinAlgError("kron multiply failed");
  }
  const product = NDArray._fromPtr(productPtr, module);

  // Compute final shape by multiplying corresponding dimensions
  const resultShape: number[] = [];
  for (let i = 0; i < nd; i++) {
    resultShape.push(aArr.shape[i] * bArr.shape[i]);
  }

  return product.reshape(resultShape);
}

/**
 * Return the cross product of two (arrays of) vectors.
 *
 * The cross product of a and b in R^3 is a vector perpendicular to both
 * a and b. The vectors are defined by the last axis by default.
 *
 * @param a - First vector(s)
 * @param b - Second vector(s)
 * @param axisa - Axis of a that defines the vector(s) (default: -1)
 * @param axisb - Axis of b that defines the vector(s) (default: -1)
 * @param axisc - Axis of c that contains the cross product vectors (default: -1)
 * @param axis - If specified, overrides axisa, axisb, and axisc
 * @returns Cross product vector(s)
 *
 * @example
 * // 3D cross product
 * cross([1, 0, 0], [0, 1, 0])  // [0, 0, 1]
 *
 * // Batch cross product
 * cross([[1, 0, 0], [0, 1, 0]], [[0, 1, 0], [0, 0, 1]])
 * // [[0, 0, 1], [1, 0, 0]]
 */
export async function cross(
  a: NDArray,
  b: NDArray,
  axisa: number = -1,
  axisb: number = -1,
  axisc: number = -1,
  axis: number | null = null,
): Promise<NDArray> {
  // Override with axis if provided
  if (axis !== null) {
    axisa = axis;
    axisb = axis;
    axisc = axis;
  }

  // Normalize axes
  axisa = axisa < 0 ? a.ndim + axisa : axisa;
  axisb = axisb < 0 ? b.ndim + axisb : axisb;

  // Get vector dimensions
  const dimA = a.shape[axisa];
  const dimB = b.shape[axisb];

  // Validate dimensions
  if (dimA !== 3 || dimB !== 3) {
    throw new LinAlgError(
      `cross product requires vectors of length 3, got ${dimA} and ${dimB}`,
    );
  }

  // Move vector axis to last position for easier computation
  const aMoved = a.moveaxis(axisa, -1);
  const bMoved = b.moveaxis(axisb, -1);

  // For simple 1D vectors, compute directly
  if (aMoved.ndim === 1 && bMoved.ndim === 1) {
    const aData = await aMoved.toTypedArray();
    const bData = await bMoved.toTypedArray();

    const c0 = aData[1] * bData[2] - aData[2] * bData[1];
    const c1 = aData[2] * bData[0] - aData[0] * bData[2];
    const c2 = aData[0] * bData[1] - aData[1] * bData[0];

    return NDArray.fromTypedArray(
      new Float64Array([c0, c1, c2]),
      [3],
      DType.Float64,
    );
  }

  // For batched case, ensure shapes match
  const batchShapeA = aMoved.shape.slice(0, -1);
  const batchShapeB = bMoved.shape.slice(0, -1);

  // Simple shape check (require exact match for now)
  if (
    batchShapeA.length !== batchShapeB.length ||
    !batchShapeA.every((s, i) => s === batchShapeB[i])
  ) {
    throw new LinAlgError(
      "cross: batch shapes must match for batched cross product",
    );
  }

  const batchShape = batchShapeA;
  const batchSize = batchShape.reduce((p, v) => p * v, 1);

  // Get typed arrays for computation
  const aData = await aMoved.toTypedArray();
  const bData = await bMoved.toTypedArray();

  // Compute cross product: c = a × b
  const resultData = new Float64Array(batchSize * 3);

  for (let i = 0; i < batchSize; i++) {
    const base = i * 3;
    const a0 = aData[base];
    const a1 = aData[base + 1];
    const a2 = aData[base + 2];
    const b0 = bData[base];
    const b1 = bData[base + 1];
    const b2 = bData[base + 2];

    resultData[base] = a1 * b2 - a2 * b1;
    resultData[base + 1] = a2 * b0 - a0 * b2;
    resultData[base + 2] = a0 * b1 - a1 * b0;
  }

  // Create result array
  const resultShape = [...batchShape, 3];
  let result = await NDArray.fromTypedArray(
    resultData,
    resultShape,
    DType.Float64,
  );

  // Move result axis if needed
  axisc = axisc < 0 ? result.ndim + axisc : axisc;
  if (axisc !== result.ndim - 1) {
    result = result.moveaxis(-1, axisc);
  }

  return result;
}

/* ============ NumPy 2.0 Functions ============ */

/**
 * Transposes the last two dimensions of an array.
 *
 * NumPy 2.0 addition: Equivalent to swapaxes(a, -2, -1).
 *
 * @param x - Input array, must have at least 2 dimensions
 * @returns Array with last two axes swapped
 *
 * @example
 * const a = await NDArray.fromArray([[1, 2], [3, 4]]);
 * const t = await matrix_transpose(a);
 * // [[1, 3], [2, 4]]
 */
export function matrix_transpose(x: NDArray): NDArray {
  if (x.ndim < 2) {
    throw new LinAlgError("matrix_transpose requires at least 2-dimensional input");
  }
  return x.swapaxes(-2, -1);
}

/**
 * Compute the vector dot product of two arrays along the last axis.
 *
 * NumPy 2.0 addition. For 1-D arrays, this is the same as dot.
 * For N-D arrays, sums the product of elements along the last axis.
 *
 * @param x1 - First input array
 * @param x2 - Second input array
 * @param axis - Axis along which to compute the dot product (default: -1)
 * @returns Vector dot product
 *
 * @example
 * const a = await NDArray.fromArray([1, 2, 3]);
 * const b = await NDArray.fromArray([4, 5, 6]);
 * const d = await vecdot(a, b);
 * // 32 (1*4 + 2*5 + 3*6)
 */
export async function vecdot(
  x1: NDArray,
  x2: NDArray,
  axis: number = -1,
): Promise<NDArray | number> {
  // Normalize axis
  const ax1 = axis < 0 ? x1.ndim + axis : axis;
  const ax2 = axis < 0 ? x2.ndim + axis : axis;

  if (ax1 < 0 || ax1 >= x1.ndim) {
    throw new LinAlgError(`axis ${axis} is out of bounds for x1 with ${x1.ndim} dimensions`);
  }
  if (ax2 < 0 || ax2 >= x2.ndim) {
    throw new LinAlgError(`axis ${axis} is out of bounds for x2 with ${x2.ndim} dimensions`);
  }

  // Check dimension sizes match
  if (x1.shape[ax1] !== x2.shape[ax2]) {
    throw new LinAlgError(
      `shape mismatch: x1[${ax1}]=${x1.shape[ax1]} != x2[${ax2}]=${x2.shape[ax2]}`
    );
  }

  // For 1D arrays, use dot directly
  if (x1.ndim === 1 && x2.ndim === 1) {
    return vdot(x1, x2);
  }

  // Move target axis to last position
  const x1Moved = ax1 !== x1.ndim - 1 ? x1.moveaxis(ax1, -1) : x1;
  const x2Moved = ax2 !== x2.ndim - 1 ? x2.moveaxis(ax2, -1) : x2;

  // Multiply and sum along last axis
  const module = getWasmModule();
  const productPtr = module._ufunc_multiply(x1Moved._wasmPtr, x2Moved._wasmPtr);
  if (productPtr === 0) {
    throw new LinAlgError("vecdot: multiply failed");
  }
  const product = NDArray._fromPtr(productPtr, module);

  // Sum along last axis
  const axisVal = -1;
  const sumPtr = module._ndarray_sum_axis(product._wasmPtr, axisVal, false, -1);
  product.dispose();

  if (sumPtr === 0) {
    throw new LinAlgError("vecdot: sum failed");
  }

  const result = NDArray._fromPtr(sumPtr, module);

  // Return scalar for 1D results
  if (result.ndim === 0) {
    const val = result.item();
    result.dispose();
    return val as number;
  }

  return result;
}

/**
 * Matrix-vector product.
 *
 * NumPy 2.0 addition. Computes the matrix-vector product of x1 and x2.
 * Equivalent to matmul with 2D @ 1D broadcasting.
 *
 * @param x1 - Matrix (..., M, N)
 * @param x2 - Vector (..., N)
 * @returns Result vector (..., M)
 *
 * @example
 * const A = await NDArray.fromArray([[1, 2], [3, 4]]);
 * const v = await NDArray.fromArray([1, 2]);
 * const result = await matvec(A, v);
 * // [5, 11] (1*1+2*2, 3*1+4*2)
 */
export async function matvec(x1: NDArray, x2: NDArray): Promise<NDArray> {
  if (x1.ndim < 2) {
    throw new LinAlgError("matvec: x1 must be at least 2-dimensional");
  }
  if (x2.ndim < 1) {
    throw new LinAlgError("matvec: x2 must be at least 1-dimensional");
  }

  // Check dimensions match
  const N1 = x1.shape[x1.ndim - 1];
  const N2 = x2.shape[x2.ndim - 1];
  if (N1 !== N2) {
    throw new LinAlgError(
      `matvec: x1.shape[-1]=${N1} does not match x2.shape[-1]=${N2}`
    );
  }

  // Reshape x2 to column vector for matmul
  const x2Col = x2.reshape([...x2.shape.slice(0, -1), N2, 1]);

  // Use matmul
  const result = await matmul(x1, x2Col);

  // Remove the trailing dimension of 1
  return result.reshape(result.shape.slice(0, -1));
}

/**
 * Vector-matrix product.
 *
 * NumPy 2.0 addition. Computes the vector-matrix product of x1 and x2.
 * Equivalent to matmul with 1D @ 2D broadcasting.
 *
 * @param x1 - Vector (..., M)
 * @param x2 - Matrix (..., M, N)
 * @returns Result vector (..., N)
 *
 * @example
 * const v = await NDArray.fromArray([1, 2]);
 * const A = await NDArray.fromArray([[1, 2], [3, 4]]);
 * const result = await vecmat(v, A);
 * // [7, 10] (1*1+2*3, 1*2+2*4)
 */
export async function vecmat(x1: NDArray, x2: NDArray): Promise<NDArray> {
  if (x1.ndim < 1) {
    throw new LinAlgError("vecmat: x1 must be at least 1-dimensional");
  }
  if (x2.ndim < 2) {
    throw new LinAlgError("vecmat: x2 must be at least 2-dimensional");
  }

  // Check dimensions match
  const M1 = x1.shape[x1.ndim - 1];
  const M2 = x2.shape[x2.ndim - 2];
  if (M1 !== M2) {
    throw new LinAlgError(
      `vecmat: x1.shape[-1]=${M1} does not match x2.shape[-2]=${M2}`
    );
  }

  // Reshape x1 to row vector for matmul
  const x1Row = x1.reshape([...x1.shape.slice(0, -1), 1, M1]);

  // Use matmul
  const result = await matmul(x1Row, x2);

  // Remove the second-to-last dimension of 1
  const newShape = [...result.shape.slice(0, -2), result.shape[result.ndim - 1]];
  return result.reshape(newShape);
}
