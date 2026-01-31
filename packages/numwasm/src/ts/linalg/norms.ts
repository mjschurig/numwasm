/**
 * Linear Algebra Norms and Numbers for NumJS-WASM
 *
 * Matrix/vector norms, determinants, and related computations.
 */

import { NDArray } from "../_core/NDArray.js";
import { DType } from "../types.js";
import { getWasmModule } from "../wasm-loader.js";
import { diagonal } from "../indexing.js";
import {
  LinAlgError,
  SlogdetResult,
  prepareMatrix,
  assertSquare,
} from "./types.js";
import { svdvals } from "./decomposition.js";
import { inv } from "./solve.js";

/**
 * Matrix or vector norm.
 *
 * @param x - Input array
 * @param ord - Order of the norm (default: 2 for vectors, 'fro' for matrices)
 *              Vector norms: 1, 2, inf, -inf, or any positive number
 *              Matrix norms: 1, 2, inf, -inf, 'fro', 'nuc'
 *
 * @example
 * const v = await NDArray.fromArray([1, 2, 3]);
 * const n = await norm(v);  // L2 norm
 */
export async function norm(
  x: NDArray,
  ord: number | "fro" | "nuc" | null = null,
): Promise<number> {
  const xPrep =
    x.dtype === DType.Float64 || x.dtype === DType.Float32
      ? x
      : x.astype(DType.Float64);

  // Map ord to integer for WASM
  let ordInt: number;
  if (ord === null || ord === "fro") {
    ordInt = 0; // Frobenius / L2
  } else if (ord === "nuc") {
    // Nuclear norm = sum of singular values
    const s = await svdvals(x);
    const sArr = s.toArray() as number[];
    return sArr.reduce((a, b) => a + b, 0);
  } else if (ord === Infinity) {
    ordInt = -1; // Inf norm
  } else if (ord === -Infinity) {
    ordInt = -2; // -Inf norm
  } else {
    ordInt = ord as number;
  }

  const module = getWasmModule();
  const resultPtr = module._linalg_norm(xPrep._wasmPtr, ordInt);

  if (resultPtr === 0) {
    throw new LinAlgError("norm computation failed");
  }

  const result = NDArray._fromPtr(resultPtr, module);
  return result.item() as number;
}

/**
 * Compute the determinant of an array.
 *
 * @param a - Square matrix (N x N)
 * @returns Determinant value
 *
 * @example
 * const a = await NDArray.fromArray([[1, 2], [3, 4]]);
 * const d = await det(a);  // -2
 */
export async function det(a: NDArray): Promise<number> {
  assertSquare(a, "a");

  const aPrep = await prepareMatrix(a, "a");

  const module = getWasmModule();
  const resultPtr = module._linalg_det(aPrep._wasmPtr);

  if (resultPtr === 0) {
    throw new LinAlgError("determinant computation failed");
  }

  const result = NDArray._fromPtr(resultPtr, module);
  return result.item() as number;
}

/**
 * Compute sign and (natural) logarithm of the determinant.
 * More numerically stable than computing det directly for large matrices.
 */
export async function slogdet(a: NDArray): Promise<SlogdetResult> {
  assertSquare(a, "a");

  const d = await det(a);
  const sign = Math.sign(d);
  const logabsdet = Math.log(Math.abs(d));

  return { sign, logabsdet };
}

/**
 * Return matrix rank using SVD method.
 *
 * @param A - Matrix to determine rank of
 * @param tol - Threshold for considering singular values as zero
 */
export async function matrix_rank(
  A: NDArray,
  tol: number | null = null,
): Promise<number> {
  if (A.ndim < 2) {
    throw new LinAlgError(
      "matrix_rank requires at least a 2-dimensional array",
    );
  }

  const s = await svdvals(A);
  const sArr = s.toArray() as number[];

  if (tol === null) {
    // Default tolerance: max(M, N) * eps * max(singular values)
    const m = A.shape[A.ndim - 2];
    const n = A.shape[A.ndim - 1];
    const eps = A.dtype === DType.Float32 ? 1.19e-7 : 2.22e-16;
    const sMax = Math.max(...sArr.map(Math.abs));
    tol = Math.max(m, n) * eps * sMax;
  }

  return sArr.filter((v) => Math.abs(v) > tol).length;
}

/**
 * Return the sum along diagonals of the array.
 */
export async function trace(
  a: NDArray,
  offset: number = 0,
  axis1: number = 0,
  axis2: number = 1,
): Promise<number> {
  const diag = await diagonal(a, offset, axis1, axis2);
  const sum = await diag.sum();
  return sum as number;
}

/**
 * Compute the condition number of a matrix.
 *
 * @param x - Matrix
 * @param p - Norm order (default: 2, using SVD)
 */
export async function cond(
  x: NDArray,
  p: number | null = null,
): Promise<number> {
  if (x.ndim < 2) {
    throw new LinAlgError("cond requires at least a 2-dimensional array");
  }

  if (p === null || p === 2 || p === -2) {
    // Use SVD
    const s = await svdvals(x);
    const sArr = s.toArray() as number[];

    if (p === -2) {
      return sArr[sArr.length - 1] / sArr[0];
    }
    return sArr[0] / sArr[sArr.length - 1];
  }

  // Use norm
  const normX = await norm(x, p);
  const xInv = await inv(x);
  const normXinv = await norm(xInv, p);
  return normX * normXinv;
}

/**
 * Compute the matrix norm.
 *
 * @param x - Input matrix (..., M, N)
 * @param ord - Order of the norm ('fro', 'nuc', 1, 2, -1, -2, Infinity, -Infinity)
 * @param keepdims - If true, axes are left with size one
 * @returns Matrix norm
 *
 * @example
 * matrix_norm([[1, 2], [3, 4]])  // Frobenius norm
 * matrix_norm([[1, 2], [3, 4]], 2)  // Spectral norm (largest singular value)
 */
export async function matrix_norm(
  x: NDArray,
  ord: number | "fro" | "nuc" = "fro",
  keepdims: boolean = false,
): Promise<NDArray | number> {
  if (x.ndim < 2) {
    throw new LinAlgError("matrix_norm requires at least 2-dimensional input");
  }

  // For 2D input, use existing norm function
  if (x.ndim === 2) {
    const result = await norm(x, ord);
    if (keepdims) {
      return (await NDArray.fromArray([[result]])) as NDArray;
    }
    return result;
  }

  // For N-D input, apply norm over last two axes
  // This requires iterating over batch dimensions
  const batchShape = x.shape.slice(0, -2);
  const batchSize = batchShape.reduce((p, v) => p * v, 1);
  const m = x.shape[x.ndim - 2];
  const n = x.shape[x.ndim - 1];
  const matrixSize = m * n;

  const results: number[] = [];

  // Flatten batch dimensions and get data
  const xReshaped = x.reshape([batchSize, m, n]);
  const xData = await xReshaped.toTypedArray();

  for (let i = 0; i < batchSize; i++) {
    // Extract matrix data at batch index i
    const matrixData = xData.slice(i * matrixSize, (i + 1) * matrixSize);
    const matrix = await NDArray.fromTypedArray(
      matrixData instanceof Float64Array
        ? matrixData
        : new Float64Array(matrixData),
      [m, n],
      DType.Float64,
    );
    const normValue = await norm(matrix, ord);
    results.push(normValue as number);
  }

  if (keepdims) {
    const resultShape = [...batchShape, 1, 1];
    return NDArray.fromTypedArray(
      new Float64Array(results),
      resultShape,
      DType.Float64,
    );
  }

  if (batchShape.length === 0) {
    return results[0];
  }
  return NDArray.fromTypedArray(
    new Float64Array(results),
    batchShape,
    DType.Float64,
  );
}

/**
 * Compute the vector norm.
 *
 * @param x - Input array
 * @param ord - Order of the norm (default: 2, Euclidean)
 * @param axis - Axis along which to compute. null = flatten
 * @param keepdims - If true, axes are left with size one
 * @returns Vector norm
 *
 * @example
 * vector_norm([3, 4])  // 5 (Euclidean)
 * vector_norm([3, 4], 1)  // 7 (Manhattan)
 * vector_norm([3, 4], Infinity)  // 4 (max)
 */
export async function vector_norm(
  x: NDArray,
  ord: number = 2,
  axis: number | null = null,
  keepdims: boolean = false,
): Promise<NDArray | number> {
  // Flatten case
  if (axis === null) {
    const flat = x.ravel();
    const data = await flat.toTypedArray();

    let result: number;

    if (ord === Infinity) {
      result = Math.max(...Array.from(data).map(Math.abs));
    } else if (ord === -Infinity) {
      result = Math.min(...Array.from(data).map(Math.abs));
    } else if (ord === 0) {
      result = Array.from(data).filter((v) => v !== 0).length;
    } else if (ord === 1) {
      result = Array.from(data).reduce((sum, v) => sum + Math.abs(v), 0);
    } else if (ord === 2) {
      result = Math.sqrt(Array.from(data).reduce((sum, v) => sum + v * v, 0));
    } else {
      // General p-norm
      const pSum = Array.from(data).reduce(
        (sum, v) => sum + Math.pow(Math.abs(v), ord),
        0,
      );
      result = Math.pow(pSum, 1 / ord);
    }

    if (keepdims) {
      const keepShape = x.shape.map(() => 1);
      return NDArray.fromTypedArray(
        new Float64Array([result]),
        keepShape,
        DType.Float64,
      );
    }
    return result;
  }

  // Axis-specific case
  const normalizedAxis = axis < 0 ? x.ndim + axis : axis;

  if (normalizedAxis < 0 || normalizedAxis >= x.ndim) {
    throw new LinAlgError(
      `axis ${axis} is out of bounds for array with ${x.ndim} dimensions`,
    );
  }

  // Compute norm along axis
  const resultShape = x.shape.filter((_, i) => i !== normalizedAxis);

  // Move target axis to last position for easier iteration
  const xMoved = x.moveaxis(normalizedAxis, -1);
  const axisSize = x.shape[normalizedAxis];
  const batchSize = x.size / axisSize;

  const data = await xMoved.toTypedArray();
  const results = new Float64Array(batchSize);

  for (let i = 0; i < batchSize; i++) {
    const start = i * axisSize;
    const slice = data.slice(start, start + axisSize);

    if (ord === Infinity) {
      results[i] = Math.max(...Array.from(slice).map(Math.abs));
    } else if (ord === -Infinity) {
      results[i] = Math.min(...Array.from(slice).map(Math.abs));
    } else if (ord === 0) {
      results[i] = Array.from(slice).filter((v) => v !== 0).length;
    } else if (ord === 1) {
      results[i] = Array.from(slice).reduce((sum, v) => sum + Math.abs(v), 0);
    } else if (ord === 2) {
      results[i] = Math.sqrt(
        Array.from(slice).reduce((sum, v) => sum + v * v, 0),
      );
    } else {
      const pSum = Array.from(slice).reduce(
        (sum, v) => sum + Math.pow(Math.abs(v), ord),
        0,
      );
      results[i] = Math.pow(pSum, 1 / ord);
    }
  }

  if (keepdims) {
    const keepShape = [...x.shape];
    keepShape[normalizedAxis] = 1;
    return NDArray.fromTypedArray(results, keepShape, DType.Float64);
  }

  if (resultShape.length === 0) {
    return results[0];
  }
  return NDArray.fromTypedArray(results, resultShape, DType.Float64);
}
