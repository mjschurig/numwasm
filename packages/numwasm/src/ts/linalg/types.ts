/**
 * Linear Algebra Types and Utilities for NumJS-WASM
 *
 * Error classes, result types, and helper functions.
 */

import { NDArray } from "../_core/NDArray.js";
import { DType } from "../types.js";

/**
 * Linear algebra error.
 * Raised when matrix operations fail (singular matrix, non-convergence, etc.)
 */
export class LinAlgError extends Error {
  constructor(message: string) {
    super(message);
    this.name = "LinAlgError";
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
export function getLinalgDtype(arr: NDArray): DType {
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
export async function prepareMatrix(
  a: NDArray,
  name: string = "array",
): Promise<NDArray> {
  if (a.ndim < 2) {
    throw new LinAlgError(
      `${name} must be at least 2-dimensional, got ${a.ndim}-dimensional`,
    );
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
export function assertSquare(a: NDArray, name: string = "array"): void {
  if (a.ndim < 2) {
    throw new LinAlgError(`${name} must be at least 2-dimensional`);
  }
  const m = a.shape[a.ndim - 2];
  const n = a.shape[a.ndim - 1];
  if (m !== n) {
    throw new LinAlgError(
      `${name} must be square, got shape (${a.shape.join(", ")})`,
    );
  }
}
