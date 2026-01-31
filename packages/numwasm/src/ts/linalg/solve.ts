/**
 * Linear Algebra Solving and Inverting for NumJS-WASM
 *
 * Linear system solving, matrix inversion, and related operations.
 */

import { NDArray } from "../_core/NDArray.js";
import { DType } from "../types.js";
import { getWasmModule } from "../wasm-loader.js";
import {
  LinAlgError,
  LstsqResult,
  prepareMatrix,
  assertSquare,
} from "./types.js";
import { svd, svdvals } from "./decomposition.js";
import { matmul } from "./products.js";

/**
 * Solve a linear matrix equation: a @ x = b
 *
 * @param a - Coefficient matrix (N x N)
 * @param b - Ordinate values (N,) or (N x NRHS)
 * @returns Solution x
 *
 * @example
 * const a = await NDArray.fromArray([[3, 1], [1, 2]]);
 * const b = await NDArray.fromArray([9, 8]);
 * const x = await solve(a, b);  // [2, 3]
 */
export async function solve(a: NDArray, b: NDArray): Promise<NDArray> {
  assertSquare(a, "a");

  const n = a.shape[a.ndim - 1];

  // Check b dimensions
  if (b.ndim >= 1) {
    const bFirstDim = b.ndim === 1 ? b.shape[0] : b.shape[b.ndim - 2];
    if (bFirstDim !== n) {
      throw new LinAlgError(
        `solve: last dimension of b (${bFirstDim}) must match matrix dimension (${n})`,
      );
    }
  }

  const aPrep = await prepareMatrix(a, "a");
  const bPrep = b.dtype !== aPrep.dtype ? b.astype(aPrep.dtype) : b;

  const module = getWasmModule();
  const resultPtr = module._linalg_solve(aPrep._wasmPtr, bPrep._wasmPtr);

  if (resultPtr === 0) {
    throw new LinAlgError("Singular matrix");
  }

  return NDArray._fromPtr(resultPtr, module);
}

/**
 * Compute the (multiplicative) inverse of a matrix.
 *
 * @param a - Square matrix (N x N)
 * @returns Inverse A^(-1)
 *
 * @example
 * const a = await NDArray.fromArray([[1, 2], [3, 4]]);
 * const aInv = await inv(a);
 */
export async function inv(a: NDArray): Promise<NDArray> {
  assertSquare(a, "a");

  const aPrep = await prepareMatrix(a, "a");

  const module = getWasmModule();
  const resultPtr = module._linalg_inv(aPrep._wasmPtr);

  if (resultPtr === 0) {
    throw new LinAlgError("Singular matrix");
  }

  return NDArray._fromPtr(resultPtr, module);
}

/**
 * Compute the (Moore-Penrose) pseudo-inverse of a matrix.
 *
 * @param a - Matrix to pseudo-invert
 * @param rcond - Cutoff for small singular values
 */
export async function pinv(
  a: NDArray,
  rcond: number = 1e-15,
): Promise<NDArray> {
  if (a.ndim < 2) {
    throw new LinAlgError("pinv requires at least a 2-dimensional array");
  }

  // For wide matrices (m < n), compute pinv(A) = pinv(A^T)^T
  // This works around an SVD bug with wide matrices
  const m0 = a.shape[0];
  const n0 = a.shape[1];
  if (m0 < n0) {
    const aT = a.T;
    const pinvAT = await pinv(aT, rcond);
    return pinvAT.T;
  }

  // SVD-based pseudo-inverse
  const { U, S, Vh } = await svd(a, false);

  // Get singular values as array
  const sArr = S.toArray() as number[];
  const sMax = Math.max(...sArr);
  const cutoff = rcond * sMax;

  // Compute reciprocal of singular values above threshold
  const sInv = sArr.map((s) => (s > cutoff ? 1 / s : 0));

  // pinv = Vh^H @ diag(s_inv) @ U^H
  // = V @ diag(s_inv) @ U^T  (for real matrices)

  // Create diagonal matrix from sInv
  const m = U.shape[0];
  const n = Vh.shape[1];
  const k = sArr.length;

  // V = Vh^T (transpose)
  const V = Vh.T;

  // U^T
  const Ut = U.T;

  // Scale columns of V by sInv
  // Then multiply by U^T
  // This is equivalent to V @ diag(sInv) @ U^T

  // For simplicity, construct explicitly
  // pinv[i,j] = sum_k V[i,k] * sInv[k] * Ut[k,j]

  // Use element access (.get) to correctly handle transposed views
  const resultData: number[] = new Array(n * m);

  for (let i = 0; i < n; i++) {
    for (let j = 0; j < m; j++) {
      let sum = 0;
      for (let l = 0; l < k; l++) {
        sum += (V.get(i, l) as number) * sInv[l] * (Ut.get(l, j) as number);
      }
      resultData[i * m + j] = sum;
    }
  }

  // Create result array from computed data
  return NDArray.fromArray(resultData, [n, m], { dtype: a.dtype });
}

/**
 * Return the least-squares solution to a linear matrix equation.
 *
 * Computes the vector x that approximately solves the equation a @ x = b.
 *
 * @param a - Coefficient matrix (M x N)
 * @param b - Ordinate values (M,) or (M x K)
 * @param rcond - Cutoff for small singular values
 */
export async function lstsq(
  a: NDArray,
  b: NDArray,
  rcond: number | null = null,
): Promise<LstsqResult> {
  if (a.ndim !== 2) {
    throw new LinAlgError("lstsq requires a 2-dimensional coefficient matrix");
  }

  const m = a.shape[0];
  const n = a.shape[1];

  if (rcond === null) {
    const eps = a.dtype === DType.Float32 ? 1.19e-7 : 2.22e-16;
    rcond = Math.max(m, n) * eps;
  }

  // Use pseudo-inverse approach: x = pinv(a) @ b
  const aPinv = await pinv(a, rcond);
  const x = await matmul(aPinv, b.ndim === 1 ? b.reshape([b.size, 1]) : b);

  // Get singular values for rank
  const s = await svdvals(a);
  const sArr = s.toArray() as number[];
  const sMax = Math.max(...sArr);
  const threshold = rcond * sMax;
  const rank = sArr.filter((v) => v > threshold).length;

  // Compute residuals if overdetermined
  let residuals: NDArray;
  if (m > n) {
    // residuals = sum((ax - b)^2, axis=0)
    // Simplified: return empty for now
    residuals = await NDArray.zeros([0], { dtype: a.dtype });
  } else {
    residuals = await NDArray.zeros([0], { dtype: a.dtype });
  }

  return {
    x: x.ndim === 2 && b.ndim === 1 ? x.ravel() : x,
    residuals,
    rank,
    s,
  };
}

/**
 * Raise a square matrix to the (integer) power n.
 *
 * @param a - Square matrix
 * @param n - Exponent (can be negative for inverse)
 */
export async function matrix_power(a: NDArray, n: number): Promise<NDArray> {
  assertSquare(a, "a");

  if (!Number.isInteger(n)) {
    throw new LinAlgError("exponent must be an integer");
  }

  const m = a.shape[a.ndim - 1];

  if (n === 0) {
    return NDArray.eye(m, m, 0, { dtype: a.dtype });
  }

  let base = a;
  if (n < 0) {
    base = await inv(a);
    n = -n;
  }

  // Binary exponentiation
  let result = await NDArray.eye(m, m, 0, { dtype: a.dtype });

  while (n > 0) {
    if (n & 1) {
      result = await matmul(result, base);
    }
    base = await matmul(base, base);
    n >>= 1;
  }

  return result;
}

/**
 * Solve the tensor equation a x = b for x.
 *
 * It is assumed that all indices of x are summed over in the product,
 * together with the rightmost indices of a.
 *
 * @param a - Coefficient tensor
 * @param b - Right-hand side tensor
 * @param axes - Axes in a to reorder to the right, for the solve
 * @returns Solution tensor x
 *
 * @example
 * const a = (await NDArray.eye(4)).reshape([2, 2, 2, 2]);
 * const b = (await NDArray.fromArray([1, 2, 3, 4])).reshape([2, 2]);
 * const x = await tensorsolve(a, b);
 * // tensordot(a, x, 2) ≈ b
 */
export async function tensorsolve(
  a: NDArray,
  b: NDArray,
  axes: number[] | null = null,
): Promise<NDArray> {
  let aWork = a;

  // Reorder axes if specified
  if (axes !== null) {
    const allAxes = Array.from({ length: a.ndim }, (_, i) => i);
    const remainingAxes = allAxes.filter((i) => !axes.includes(i));
    const newOrder = [...remainingAxes, ...axes];
    aWork = a.transpose(newOrder);
  }

  // The shape of x comes from the trailing dimensions of a
  const Q = aWork.ndim;
  const bNdim = b.ndim;

  // x.ndim = Q - bNdim
  const xShape = aWork.shape.slice(bNdim);
  const N = xShape.reduce((p, x) => p * x, 1);
  const M = b.size;

  // Validate dimensions
  if (N !== M) {
    throw new LinAlgError(
      `tensorsolve: incompatible dimensions. ` +
        `Product of trailing ${Q - bNdim} dims of a (${N}) must equal b.size (${M})`,
    );
  }

  // Reshape to 2D system: a becomes (M, N), b becomes (M,)
  const aMatrix = aWork.reshape([M, N]);
  const bVector = b.ravel();

  // Solve
  const xVector = await solve(aMatrix, bVector);

  // Reshape solution to tensor shape
  return xVector.reshape(xShape);
}

/**
 * Compute the 'inverse' of an N-dimensional array.
 *
 * The result is an inverse for a with respect to the tensordot operation
 * tensordot(a, b, ind).
 *
 * @param a - Tensor to pseudo-invert
 * @param ind - Number of first indices that are involved in the inverse sum
 * @returns Tensor inverse
 *
 * @example
 * const a = (await NDArray.eye(4)).reshape([2, 2, 2, 2]);
 * const ainv = await tensorinv(a);
 * // tensordot(ainv, a, ind=2) ≈ eye(4).reshape([2, 2, 2, 2])
 */
export async function tensorinv(a: NDArray, ind: number = 2): Promise<NDArray> {
  if (ind <= 0) {
    throw new LinAlgError("tensorinv: ind must be a positive integer");
  }

  const oldShape = a.shape;

  // Compute products of dimensions
  const prod1 = oldShape.slice(0, ind).reduce((p, x) => p * x, 1);
  const prod2 = oldShape.slice(ind).reduce((p, x) => p * x, 1);

  if (prod1 !== prod2) {
    throw new LinAlgError(
      `tensorinv: product of first ${ind} dimensions (${prod1}) ` +
        `must equal product of remaining dimensions (${prod2})`,
    );
  }

  // Reshape to square 2D matrix
  const aMatrix = a.reshape([prod1, prod2]);

  // Compute inverse
  const aInvMatrix = await inv(aMatrix);

  // Reshape back: new shape is shape[ind:] + shape[:ind]
  const newShape = [...oldShape.slice(ind), ...oldShape.slice(0, ind)];
  return aInvMatrix.reshape(newShape);
}
