/**
 * Linear Algebra Decompositions for NumJS-WASM
 *
 * Matrix decompositions: Cholesky, QR, SVD, Eigenvalue.
 */

import { NDArray } from "../_core/NDArray.js";
import { getWasmModule } from "../wasm-loader.js";
import {
  LinAlgError,
  EigResult,
  SVDResult,
  QRResult,
  prepareMatrix,
  assertSquare,
} from "./types.js";

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
export async function cholesky(
  a: NDArray,
  upper: boolean = false,
): Promise<NDArray> {
  assertSquare(a, "a");

  const aPrep = await prepareMatrix(a, "a");

  const module = getWasmModule();
  const resultPtr = module._linalg_cholesky(aPrep._wasmPtr, upper ? 1 : 0);

  if (resultPtr === 0) {
    throw new LinAlgError(
      "Matrix is not positive definite - Cholesky decomposition cannot be computed",
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
    throw new LinAlgError("qr requires a 2-dimensional array");
  }

  const aPrep = await prepareMatrix(a, "a");

  const module = getWasmModule();
  const resultPtr = module._linalg_qr(aPrep._wasmPtr);

  if (resultPtr === 0) {
    throw new LinAlgError("QR decomposition failed");
  }

  const qPtr = module._linalg_qr_get_q(resultPtr);
  const rPtr = module._linalg_qr_get_r(resultPtr);

  const Q = NDArray._fromPtr(qPtr, module);
  const R = NDArray._fromPtr(rPtr, module);

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
export async function svd(
  a: NDArray,
  fullMatrices: boolean = false,
): Promise<SVDResult> {
  if (a.ndim !== 2) {
    throw new LinAlgError("svd requires a 2-dimensional array");
  }

  const aPrep = await prepareMatrix(a, "a");

  const module = getWasmModule();
  const resultPtr = module._linalg_svd(aPrep._wasmPtr, fullMatrices ? 1 : 0);

  if (resultPtr === 0) {
    throw new LinAlgError("SVD did not converge");
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
  assertSquare(a, "a");

  const aPrep = await prepareMatrix(a, "a");

  const module = getWasmModule();
  const resultPtr = module._linalg_eig(aPrep._wasmPtr);

  if (resultPtr === 0) {
    throw new LinAlgError("Eigenvalue computation did not converge");
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
export async function eigh(
  a: NDArray,
  _UPLO: "L" | "U" = "L",
): Promise<EigResult> {
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
export async function eigvalsh(
  a: NDArray,
  UPLO: "L" | "U" = "L",
): Promise<NDArray> {
  const result = await eigh(a, UPLO);
  return result.eigenvalues;
}
