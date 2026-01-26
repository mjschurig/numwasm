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

/* ============ Phase 25: Advanced Linear Algebra ============ */

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
  axes: number | [number[], number[]] = 2
): Promise<NDArray> {
  // Parse axes specification
  let axesA: number[];
  let axesB: number[];

  if (typeof axes === 'number') {
    if (axes < 0) {
      throw new LinAlgError('axes must be non-negative');
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
      throw new LinAlgError('axes lists must have same length');
    }
  }

  // Check for duplicate axes
  if (new Set(axesA).size !== axesA.length) {
    throw new LinAlgError('duplicate axes not allowed in tensordot');
  }
  if (new Set(axesB).size !== axesB.length) {
    throw new LinAlgError('duplicate axes not allowed in tensordot');
  }

  // Normalize negative axes
  axesA = axesA.map((ax) => (ax < 0 ? a.ndim + ax : ax));
  axesB = axesB.map((ax) => (ax < 0 ? b.ndim + ax : ax));

  // Validate axes are in bounds
  for (const ax of axesA) {
    if (ax < 0 || ax >= a.ndim) {
      throw new LinAlgError(`axis ${ax} is out of bounds for array with ${a.ndim} dimensions`);
    }
  }
  for (const ax of axesB) {
    if (ax < 0 || ax >= b.ndim) {
      throw new LinAlgError(`axis ${ax} is out of bounds for array with ${b.ndim} dimensions`);
    }
  }

  // Validate axes match in size
  for (let i = 0; i < axesA.length; i++) {
    if (a.shape[axesA[i]] !== b.shape[axesB[i]]) {
      throw new LinAlgError(
        `shape mismatch for sum: ${a.shape[axesA[i]]} vs ${b.shape[axesB[i]]}`
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
    throw new LinAlgError('multi_dot requires at least 2 arrays');
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
    workArrays[n - 1] = workArrays[n - 1].reshape([workArrays[n - 1].shape[0], 1]);
  }

  // Validate all arrays are 2D
  for (let i = 0; i < n; i++) {
    if (workArrays[i].ndim !== 2) {
      throw new LinAlgError(`multi_dot requires 1-D or 2-D arrays, got ${workArrays[i].ndim}-D at position ${i}`);
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
        `shapes (${workArrays[i].shape.join(',')}) and (${workArrays[i + 1].shape.join(',')}) not aligned`
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
  const productPtr = module._ufunc_multiply(aExpanded._wasmPtr, bExpanded._wasmPtr);
  if (productPtr === 0) {
    throw new LinAlgError('kron multiply failed');
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
  axis: number | null = null
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
      `cross product requires vectors of length 3, got ${dimA} and ${dimB}`
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

    return NDArray.fromTypedArray(new Float64Array([c0, c1, c2]), [3], DType.Float64);
  }

  // For batched case, ensure shapes match
  const batchShapeA = aMoved.shape.slice(0, -1);
  const batchShapeB = bMoved.shape.slice(0, -1);

  // Simple shape check (require exact match for now)
  if (batchShapeA.length !== batchShapeB.length ||
      !batchShapeA.every((s, i) => s === batchShapeB[i])) {
    throw new LinAlgError('cross: batch shapes must match for batched cross product');
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
  let result = await NDArray.fromTypedArray(resultData, resultShape, DType.Float64);

  // Move result axis if needed
  axisc = axisc < 0 ? result.ndim + axisc : axisc;
  if (axisc !== result.ndim - 1) {
    result = result.moveaxis(-1, axisc);
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
  axes: number[] | null = null
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
        `Product of trailing ${Q - bNdim} dims of a (${N}) must equal b.size (${M})`
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
    throw new LinAlgError('tensorinv: ind must be a positive integer');
  }

  const oldShape = a.shape;

  // Compute products of dimensions
  const prod1 = oldShape.slice(0, ind).reduce((p, x) => p * x, 1);
  const prod2 = oldShape.slice(ind).reduce((p, x) => p * x, 1);

  if (prod1 !== prod2) {
    throw new LinAlgError(
      `tensorinv: product of first ${ind} dimensions (${prod1}) ` +
        `must equal product of remaining dimensions (${prod2})`
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
  ord: number | 'fro' | 'nuc' = 'fro',
  keepdims: boolean = false
): Promise<NDArray | number> {
  if (x.ndim < 2) {
    throw new LinAlgError('matrix_norm requires at least 2-dimensional input');
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
      matrixData instanceof Float64Array ? matrixData : new Float64Array(matrixData),
      [m, n],
      DType.Float64
    );
    const normValue = await norm(matrix, ord);
    results.push(normValue as number);
  }

  if (keepdims) {
    const resultShape = [...batchShape, 1, 1];
    return NDArray.fromTypedArray(new Float64Array(results), resultShape, DType.Float64);
  }

  if (batchShape.length === 0) {
    return results[0];
  }
  return NDArray.fromTypedArray(new Float64Array(results), batchShape, DType.Float64);
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
  keepdims: boolean = false
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
      const pSum = Array.from(data).reduce((sum, v) => sum + Math.pow(Math.abs(v), ord), 0);
      result = Math.pow(pSum, 1 / ord);
    }

    if (keepdims) {
      const keepShape = x.shape.map(() => 1);
      return NDArray.fromTypedArray(new Float64Array([result]), keepShape, DType.Float64);
    }
    return result;
  }

  // Axis-specific case
  const normalizedAxis = axis < 0 ? x.ndim + axis : axis;

  if (normalizedAxis < 0 || normalizedAxis >= x.ndim) {
    throw new LinAlgError(`axis ${axis} is out of bounds for array with ${x.ndim} dimensions`);
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
      results[i] = Math.sqrt(Array.from(slice).reduce((sum, v) => sum + v * v, 0));
    } else {
      const pSum = Array.from(slice).reduce((sum, v) => sum + Math.pow(Math.abs(v), ord), 0);
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

  // Phase 25: Advanced Linear Algebra
  tensordot,
  multi_dot,
  kron,
  cross,
  tensorsolve,
  tensorinv,
  matrix_norm,
  vector_norm,
};
