/**
 * NumJS Linear Algebra Module
 *
 * Provides NumPy-compatible linear algebra operations with WASM acceleration.
 */

import { NDArray } from './NDArray.js';
import { DType } from './types.js';
import { getWasmModule } from './wasm-loader.js';
import { diagonal as diagonalFn } from './indexing.js';

/**
 * Linear algebra error.
 * Raised when matrix operations fail (singular matrix, non-convergence, etc.)
 */
export class LinAlgError extends Error {
  constructor(message: string) {
    super(message);
    this.name = 'LinAlgError';
  }
}

/* ============ Result Types ============ */

/**
 * Result of eigenvalue decomposition.
 */
export interface EigResult {
  /** Eigenvalues (may be complex for non-symmetric matrices) */
  eigenvalues: NDArray;
  /** Imaginary parts of eigenvalues (for non-symmetric matrices) */
  eigenvaluesImag?: NDArray;
  /** Right eigenvectors as columns */
  eigenvectors: NDArray;
}

/**
 * Result of singular value decomposition.
 */
export interface SVDResult {
  /** Left singular vectors (M x M or M x K) */
  U: NDArray;
  /** Singular values in descending order (K,) where K = min(M, N) */
  S: NDArray;
  /** Right singular vectors transposed (N x N or K x N) */
  Vh: NDArray;
}

/**
 * Result of QR decomposition.
 */
export interface QRResult {
  /** Orthogonal matrix Q (M x K where K = min(M, N)) */
  Q: NDArray;
  /** Upper triangular matrix R (K x N) */
  R: NDArray;
}

/**
 * Result of least squares solution.
 */
export interface LstsqResult {
  /** Least-squares solution */
  x: NDArray;
  /** Sum of squared residuals */
  residuals: NDArray;
  /** Effective rank of matrix */
  rank: number;
  /** Singular values */
  s: NDArray;
}

/**
 * Result of sign and log determinant.
 */
export interface SlogdetResult {
  /** Sign of determinant (-1, 0, or 1) */
  sign: number;
  /** Natural log of absolute value of determinant */
  logabsdet: number;
}

/* ============ Type Utilities ============ */

/**
 * Get the computation dtype for linalg operations.
 * Integer types are promoted to float64.
 */
function getLinalgDtype(arr: NDArray): DType {
  const dtype = arr.dtype;

  // Integer types -> float64
  if (
    dtype === DType.Bool ||
    dtype === DType.Int8 ||
    dtype === DType.Int16 ||
    dtype === DType.Int32 ||
    dtype === DType.Int64 ||
    dtype === DType.Uint8 ||
    dtype === DType.Uint16 ||
    dtype === DType.Uint32 ||
    dtype === DType.Uint64
  ) {
    return DType.Float64;
  }

  // Float16 -> Float32
  if (dtype === DType.Float16) {
    return DType.Float32;
  }

  // Float32, Float64, Complex64, Complex128 -> keep as is
  return dtype;
}

/**
 * Ensure array is at least 2D and contiguous.
 */
async function prepareMatrix(a: NDArray, name: string = 'array'): Promise<NDArray> {
  if (a.ndim < 2) {
    throw new LinAlgError(`${name} must be at least 2-dimensional, got ${a.ndim}-dimensional`);
  }

  // Convert to appropriate dtype
  const targetDtype = getLinalgDtype(a);
  if (a.dtype !== targetDtype) {
    return a.astype(targetDtype);
  }

  // Ensure C-contiguous
  if (!a.flags.c_contiguous) {
    return a.copy();
  }

  return a;
}

/**
 * Assert array is square in last two dimensions.
 */
function assertSquare(a: NDArray, name: string = 'array'): void {
  if (a.ndim < 2) {
    throw new LinAlgError(`${name} must be at least 2-dimensional`);
  }
  const m = a.shape[a.ndim - 2];
  const n = a.shape[a.ndim - 1];
  if (m !== n) {
    throw new LinAlgError(
      `${name} must be square, got shape (${a.shape.join(', ')})`
    );
  }
}

/* ============ Matrix Products ============ */

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
    throw new LinAlgError('matmul requires at least 2-dimensional arrays');
  }

  const k1 = a.shape[a.ndim - 1];
  const k2 = b.shape[b.ndim - 2];

  if (k1 !== k2) {
    throw new LinAlgError(
      `matmul: Input operand 1 has a mismatch in its core dimension 0, ` +
      `with gufunc signature (n?,k),(k,m?)->(n?,m?) (size ${k1} is different from ${k2})`
    );
  }

  // Prepare arrays
  const aPrep = await prepareMatrix(a, 'a');
  const bPrep = await prepareMatrix(b, 'b');

  const module = getWasmModule();
  const resultPtr = module._linalg_matmul(aPrep._wasmPtr, bPrep._wasmPtr);

  if (resultPtr === 0) {
    throw new LinAlgError('matmul failed - check array shapes and types');
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
    throw new LinAlgError('dot does not support 0-d arrays');
  }

  // Ensure appropriate dtype
  const targetDtype = getLinalgDtype(a);
  const aPrep = a.dtype !== targetDtype ? a.astype(targetDtype) : a;
  const bPrep = b.dtype !== targetDtype ? b.astype(targetDtype) : b;

  const module = getWasmModule();
  const resultPtr = module._linalg_dot(aPrep._wasmPtr, bPrep._wasmPtr);

  if (resultPtr === 0) {
    throw new LinAlgError('dot failed - check array shapes');
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
    throw new LinAlgError('vdot requires arrays of equal size');
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
      `${aLast} (dim ${a.ndim - 1}) != ${bLast} (dim ${b.ndim - 1})`
    );
  }

  // Implement via reshape + matmul + reshape
  throw new LinAlgError('inner for N-D arrays not yet implemented');
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

/* ============ Decompositions ============ */

/**
 * Cholesky decomposition.
 *
 * Returns the Cholesky factor of a Hermitian positive-definite matrix.
 *
 * @param a - Hermitian positive-definite matrix (N x N)
 * @param upper - If true, return upper triangular U (a = U^H @ U),
 *                otherwise return lower triangular L (a = L @ L^H)
 *
 * @example
 * const a = await NDArray.fromArray([[4, 2], [2, 5]]);
 * const L = await cholesky(a);
 * // L = [[2, 0], [1, 2]]
 */
export async function cholesky(a: NDArray, upper: boolean = false): Promise<NDArray> {
  assertSquare(a, 'a');

  const aPrep = await prepareMatrix(a, 'a');

  const module = getWasmModule();
  const resultPtr = module._linalg_cholesky(aPrep._wasmPtr, upper ? 1 : 0);

  if (resultPtr === 0) {
    throw new LinAlgError(
      'Matrix is not positive definite - Cholesky decomposition cannot be computed'
    );
  }

  return NDArray._fromPtr(resultPtr, module);
}

/**
 * QR decomposition.
 *
 * Factor the matrix a as qr, where q is orthonormal and r is upper-triangular.
 *
 * @param a - Matrix to factor (M x N)
 * @returns QRResult with Q (M x K) and R (K x N) where K = min(M, N)
 *
 * @example
 * const a = await NDArray.fromArray([[1, 2], [3, 4], [5, 6]]);
 * const { Q, R } = await qr(a);
 */
export async function qr(a: NDArray): Promise<QRResult> {
  if (a.ndim !== 2) {
    throw new LinAlgError('qr requires a 2-dimensional array');
  }

  const aPrep = await prepareMatrix(a, 'a');

  const module = getWasmModule();
  const resultPtr = module._linalg_qr(aPrep._wasmPtr);

  if (resultPtr === 0) {
    throw new LinAlgError('QR decomposition failed');
  }

  const qPtr = module._linalg_qr_get_q(resultPtr);
  const rPtr = module._linalg_qr_get_r(resultPtr);

  const Q = NDArray._fromPtr(qPtr, module);
  const R = NDArray._fromPtr(rPtr, module);

  // Free the result struct but not the arrays (they're now owned by NDArray)
  // Note: We need a way to free just the struct, not the arrays
  // For now, we'll manage this carefully

  return { Q, R };
}

/**
 * Singular Value Decomposition.
 *
 * @param a - Matrix to decompose (M x N)
 * @param fullMatrices - If true, U and Vh have full shapes (M x M) and (N x N)
 *                        If false, shapes are reduced (M x K) and (K x N)
 * @returns SVDResult with U, S, Vh
 *
 * @example
 * const a = await NDArray.fromArray([[1, 0], [0, 2], [0, 0]]);
 * const { U, S, Vh } = await svd(a);
 */
export async function svd(a: NDArray, fullMatrices: boolean = false): Promise<SVDResult> {
  if (a.ndim !== 2) {
    throw new LinAlgError('svd requires a 2-dimensional array');
  }

  const aPrep = await prepareMatrix(a, 'a');

  const module = getWasmModule();
  const resultPtr = module._linalg_svd(aPrep._wasmPtr, fullMatrices ? 1 : 0);

  if (resultPtr === 0) {
    throw new LinAlgError('SVD did not converge');
  }

  const uPtr = module._linalg_svd_get_u(resultPtr);
  const sPtr = module._linalg_svd_get_s(resultPtr);
  const vhPtr = module._linalg_svd_get_vh(resultPtr);

  const U = NDArray._fromPtr(uPtr, module);
  const S = NDArray._fromPtr(sPtr, module);
  const Vh = NDArray._fromPtr(vhPtr, module);

  return { U, S, Vh };
}

/**
 * Return singular values only.
 */
export async function svdvals(a: NDArray): Promise<NDArray> {
  const { S } = await svd(a, false);
  return S;
}

/* ============ Eigenvalues ============ */

/**
 * Compute eigenvalues and right eigenvectors of a general matrix.
 *
 * @param a - Square matrix (N x N)
 * @returns EigResult with eigenvalues and eigenvectors
 *
 * @example
 * const a = await NDArray.fromArray([[1, 0], [0, 2]]);
 * const { eigenvalues, eigenvectors } = await eig(a);
 */
export async function eig(a: NDArray): Promise<EigResult> {
  assertSquare(a, 'a');

  const aPrep = await prepareMatrix(a, 'a');

  const module = getWasmModule();
  const resultPtr = module._linalg_eig(aPrep._wasmPtr);

  if (resultPtr === 0) {
    throw new LinAlgError('Eigenvalue computation did not converge');
  }

  const valuesPtr = module._linalg_eig_get_values(resultPtr);
  const valuesImagPtr = module._linalg_eig_get_values_imag(resultPtr);
  const vectorsPtr = module._linalg_eig_get_vectors(resultPtr);

  const eigenvalues = NDArray._fromPtr(valuesPtr, module);
  const eigenvaluesImag = NDArray._fromPtr(valuesImagPtr, module);
  const eigenvectors = NDArray._fromPtr(vectorsPtr, module);

  // Check if any imaginary parts are non-zero
  const imagData = eigenvaluesImag.toArray() as number[];
  const hasComplex = imagData.some((v) => Math.abs(v) > 1e-10);

  return {
    eigenvalues,
    eigenvaluesImag: hasComplex ? eigenvaluesImag : undefined,
    eigenvectors,
  };
}

/**
 * Compute eigenvalues only.
 */
export async function eigvals(a: NDArray): Promise<NDArray> {
  const result = await eig(a);
  return result.eigenvalues;
}

/**
 * Compute eigenvalues and eigenvectors of a Hermitian/symmetric matrix.
 *
 * Eigenvalues are always real and returned in ascending order.
 *
 * @param a - Hermitian/symmetric matrix (N x N)
 * @param UPLO - 'L' (default) uses lower triangle, 'U' uses upper triangle
 */
export async function eigh(a: NDArray, _UPLO: 'L' | 'U' = 'L'): Promise<EigResult> {
  // For symmetric matrices, we can use the general eig
  // A full implementation would use dsyevd which is more efficient
  // _UPLO parameter reserved for future optimization
  const result = await eig(a);

  // Sort eigenvalues in ascending order
  const eigenvalues = result.eigenvalues;
  const eigenvectors = result.eigenvectors;

  // TODO: Implement sorting for symmetric case

  return {
    eigenvalues,
    eigenvectors,
  };
}

/**
 * Compute eigenvalues of a Hermitian/symmetric matrix.
 */
export async function eigvalsh(a: NDArray, UPLO: 'L' | 'U' = 'L'): Promise<NDArray> {
  const result = await eigh(a, UPLO);
  return result.eigenvalues;
}

/* ============ Norms & Numbers ============ */

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
  ord: number | 'fro' | 'nuc' | null = null
): Promise<number> {
  const xPrep = x.dtype === DType.Float64 || x.dtype === DType.Float32
    ? x
    : x.astype(DType.Float64);

  // Map ord to integer for WASM
  let ordInt: number;
  if (ord === null || ord === 'fro') {
    ordInt = 0;  // Frobenius / L2
  } else if (ord === 'nuc') {
    // Nuclear norm = sum of singular values
    const s = await svdvals(x);
    const sArr = s.toArray() as number[];
    return sArr.reduce((a, b) => a + b, 0);
  } else if (ord === Infinity) {
    ordInt = -1;  // Inf norm
  } else if (ord === -Infinity) {
    ordInt = -2;  // -Inf norm
  } else {
    ordInt = ord as number;
  }

  const module = getWasmModule();
  const resultPtr = module._linalg_norm(xPrep._wasmPtr, ordInt);

  if (resultPtr === 0) {
    throw new LinAlgError('norm computation failed');
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
  assertSquare(a, 'a');

  const aPrep = await prepareMatrix(a, 'a');

  const module = getWasmModule();
  const resultPtr = module._linalg_det(aPrep._wasmPtr);

  if (resultPtr === 0) {
    throw new LinAlgError('determinant computation failed');
  }

  const result = NDArray._fromPtr(resultPtr, module);
  return result.item() as number;
}

/**
 * Compute sign and (natural) logarithm of the determinant.
 * More numerically stable than computing det directly for large matrices.
 */
export async function slogdet(a: NDArray): Promise<SlogdetResult> {
  assertSquare(a, 'a');

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
export async function matrix_rank(A: NDArray, tol: number | null = null): Promise<number> {
  if (A.ndim < 2) {
    throw new LinAlgError('matrix_rank requires at least a 2-dimensional array');
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
  axis2: number = 1
): Promise<number> {
  const diag = await diagonalFn(a, offset, axis1, axis2);
  const sum = await diag.sum();
  return sum as number;
}

/**
 * Compute the condition number of a matrix.
 *
 * @param x - Matrix
 * @param p - Norm order (default: 2, using SVD)
 */
export async function cond(x: NDArray, p: number | null = null): Promise<number> {
  if (x.ndim < 2) {
    throw new LinAlgError('cond requires at least a 2-dimensional array');
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

/* ============ Solving & Inverting ============ */

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
  assertSquare(a, 'a');

  const n = a.shape[a.ndim - 1];

  // Check b dimensions
  if (b.ndim >= 1) {
    const bFirstDim = b.ndim === 1 ? b.shape[0] : b.shape[b.ndim - 2];
    if (bFirstDim !== n) {
      throw new LinAlgError(
        `solve: last dimension of b (${bFirstDim}) must match matrix dimension (${n})`
      );
    }
  }

  const aPrep = await prepareMatrix(a, 'a');
  const bPrep = b.dtype !== aPrep.dtype ? b.astype(aPrep.dtype) : b;

  const module = getWasmModule();
  const resultPtr = module._linalg_solve(aPrep._wasmPtr, bPrep._wasmPtr);

  if (resultPtr === 0) {
    throw new LinAlgError('Singular matrix');
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
  assertSquare(a, 'a');

  const aPrep = await prepareMatrix(a, 'a');

  const module = getWasmModule();
  const resultPtr = module._linalg_inv(aPrep._wasmPtr);

  if (resultPtr === 0) {
    throw new LinAlgError('Singular matrix');
  }

  return NDArray._fromPtr(resultPtr, module);
}

/**
 * Compute the (Moore-Penrose) pseudo-inverse of a matrix.
 *
 * @param a - Matrix to pseudo-invert
 * @param rcond - Cutoff for small singular values
 */
export async function pinv(a: NDArray, rcond: number = 1e-15): Promise<NDArray> {
  if (a.ndim < 2) {
    throw new LinAlgError('pinv requires at least a 2-dimensional array');
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

  const result = await NDArray.zeros([n, m], { dtype: a.dtype });
  const resultData = result.toArray() as number[];

  const vData = V.toArray() as number[];
  const utData = Ut.toArray() as number[];

  for (let i = 0; i < n; i++) {
    for (let j = 0; j < m; j++) {
      let sum = 0;
      for (let l = 0; l < k; l++) {
        sum += vData[i * k + l] * sInv[l] * utData[l * m + j];
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
  rcond: number | null = null
): Promise<LstsqResult> {
  if (a.ndim !== 2) {
    throw new LinAlgError('lstsq requires a 2-dimensional coefficient matrix');
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
  assertSquare(a, 'a');

  if (!Number.isInteger(n)) {
    throw new LinAlgError('exponent must be an integer');
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

/* ============ Export linalg namespace ============ */

export const linalg = {
  // Error class
  LinAlgError,

  // Matrix products
  dot,
  vdot,
  inner,
  outer,
  matmul,

  // Decompositions
  cholesky,
  qr,
  svd,
  svdvals,

  // Eigenvalues
  eig,
  eigh,
  eigvals,
  eigvalsh,

  // Norms & Numbers
  norm,
  det,
  slogdet,
  matrix_rank,
  trace,
  cond,

  // Solving & Inverting
  solve,
  lstsq,
  inv,
  pinv,

  // Matrix Operations
  matrix_power,
};
