/**
 * BLAS Types
 *
 * Type definitions for BLAS (Basic Linear Algebra Subprograms) operations.
 */

import type { RealArray } from '../helpers.js';

// ============================================================
// INPUT TYPES
// ============================================================

/**
 * Matrix input type - supports both 2D arrays (row-major) and 1D typed arrays (column-major).
 */
export type Matrix = number[][] | RealArray;

/**
 * Vector input type.
 */
export type Vector = number[] | RealArray;

/**
 * Transpose operation type.
 */
export type TransposeOp = 'N' | 'T' | 'C';

// ============================================================
// BLAS LEVEL 3 (Matrix-Matrix)
// ============================================================

/**
 * Options for general matrix multiplication.
 */
export interface MatmulOptions {
  /**
   * Transpose operation on A: 'N' (none), 'T' (transpose), 'C' (conjugate transpose).
   * @default 'N'
   */
  transA?: TransposeOp;

  /**
   * Transpose operation on B: 'N' (none), 'T' (transpose), 'C' (conjugate transpose).
   * @default 'N'
   */
  transB?: TransposeOp;

  /**
   * Scalar multiplier for A*B.
   * @default 1.0
   */
  alpha?: number;

  /**
   * Scalar multiplier for C.
   * @default 0.0
   */
  beta?: number;

  /**
   * Existing matrix C to accumulate into. If not provided, a new matrix is created.
   */
  C?: Matrix;
}

/**
 * Result from matrix multiplication.
 */
export interface MatmulResult {
  /**
   * Result matrix C = alpha*op(A)*op(B) + beta*C.
   */
  C: Float64Array;

  /**
   * Number of rows of result.
   */
  m: number;

  /**
   * Number of columns of result.
   */
  n: number;

  /**
   * Whether the operation was successful.
   */
  success: boolean;

  /**
   * Status message.
   */
  message: string;
}

/**
 * Options for triangular matrix multiplication.
 */
export interface MatmulTriangularOptions {
  /**
   * Side: 'L' for A*B, 'R' for B*A.
   * @default 'L'
   */
  side?: 'L' | 'R';

  /**
   * Whether A is upper triangular.
   * @default true
   */
  upper?: boolean;

  /**
   * Transpose operation on A.
   * @default 'N'
   */
  transA?: TransposeOp;

  /**
   * Whether A has unit diagonal.
   * @default false
   */
  unitDiagonal?: boolean;

  /**
   * Scalar multiplier.
   * @default 1.0
   */
  alpha?: number;
}

/**
 * Result from triangular matrix multiplication.
 */
export interface MatmulTriangularResult {
  /**
   * Result matrix B = alpha*op(A)*B or B = alpha*B*op(A).
   */
  B: Float64Array;

  /**
   * Number of rows.
   */
  m: number;

  /**
   * Number of columns.
   */
  n: number;

  /**
   * Whether the operation was successful.
   */
  success: boolean;

  /**
   * Status message.
   */
  message: string;
}

/**
 * Options for solving triangular matrix system.
 */
export interface SolveMatrixTriangularOptions {
  /**
   * Side: 'L' for op(A)*X=B, 'R' for X*op(A)=B.
   * @default 'L'
   */
  side?: 'L' | 'R';

  /**
   * Whether A is upper triangular.
   * @default true
   */
  upper?: boolean;

  /**
   * Transpose operation on A.
   * @default 'N'
   */
  transA?: TransposeOp;

  /**
   * Whether A has unit diagonal.
   * @default false
   */
  unitDiagonal?: boolean;

  /**
   * Scalar multiplier for B.
   * @default 1.0
   */
  alpha?: number;
}

/**
 * Result from triangular matrix solve.
 */
export interface SolveMatrixTriangularResult {
  /**
   * Solution matrix X.
   */
  X: Float64Array;

  /**
   * Number of rows.
   */
  m: number;

  /**
   * Number of columns.
   */
  n: number;

  /**
   * Whether the operation was successful.
   */
  success: boolean;

  /**
   * Status message.
   */
  message: string;
}

/**
 * Options for symmetric rank-k update.
 */
export interface SyrkOptions {
  /**
   * Whether to update upper or lower triangle.
   * @default 'lower'
   */
  uplo?: 'upper' | 'lower';

  /**
   * Transpose: 'N' for C = alpha*A*A^T, 'T' for C = alpha*A^T*A.
   * @default 'N'
   */
  trans?: 'N' | 'T';

  /**
   * Scalar multiplier for A*A^T or A^T*A.
   * @default 1.0
   */
  alpha?: number;

  /**
   * Scalar multiplier for C.
   * @default 0.0
   */
  beta?: number;

  /**
   * Existing symmetric matrix C to accumulate into.
   */
  C?: Matrix;
}

/**
 * Result from symmetric rank-k update.
 */
export interface SyrkResult {
  /**
   * Result symmetric matrix C.
   */
  C: Float64Array;

  /**
   * Matrix dimension.
   */
  n: number;

  /**
   * Whether the operation was successful.
   */
  success: boolean;

  /**
   * Status message.
   */
  message: string;
}

/**
 * Options for symmetric rank-2k update.
 */
export interface Syr2kOptions {
  /**
   * Whether to update upper or lower triangle.
   * @default 'lower'
   */
  uplo?: 'upper' | 'lower';

  /**
   * Transpose: 'N' for C = alpha*A*B^T + alpha*B*A^T, 'T' for C = alpha*A^T*B + alpha*B^T*A.
   * @default 'N'
   */
  trans?: 'N' | 'T';

  /**
   * Scalar multiplier.
   * @default 1.0
   */
  alpha?: number;

  /**
   * Scalar multiplier for C.
   * @default 0.0
   */
  beta?: number;

  /**
   * Existing symmetric matrix C to accumulate into.
   */
  C?: Matrix;
}

/**
 * Result from symmetric rank-2k update.
 */
export interface Syr2kResult {
  /**
   * Result symmetric matrix C.
   */
  C: Float64Array;

  /**
   * Matrix dimension.
   */
  n: number;

  /**
   * Whether the operation was successful.
   */
  success: boolean;

  /**
   * Status message.
   */
  message: string;
}

/**
 * Options for Hermitian rank-k update.
 */
export interface HerkOptions {
  /**
   * Whether to update upper or lower triangle.
   * @default 'lower'
   */
  uplo?: 'upper' | 'lower';

  /**
   * Transpose: 'N' for C = alpha*A*A^H, 'C' for C = alpha*A^H*A.
   * @default 'N'
   */
  trans?: 'N' | 'C';

  /**
   * Real scalar multiplier.
   * @default 1.0
   */
  alpha?: number;

  /**
   * Real scalar multiplier for C.
   * @default 0.0
   */
  beta?: number;

  /**
   * Existing Hermitian matrix C to accumulate into.
   */
  C?: Matrix;
}

/**
 * Result from Hermitian rank-k update.
 */
export interface HerkResult {
  /**
   * Result Hermitian matrix C.
   */
  C: Float64Array;

  /**
   * Matrix dimension.
   */
  n: number;

  /**
   * Whether the operation was successful.
   */
  success: boolean;

  /**
   * Status message.
   */
  message: string;
}

/**
 * Options for Hermitian rank-2k update.
 */
export interface Her2kOptions {
  /**
   * Whether to update upper or lower triangle.
   * @default 'lower'
   */
  uplo?: 'upper' | 'lower';

  /**
   * Transpose: 'N' for C = alpha*A*B^H + conj(alpha)*B*A^H.
   * @default 'N'
   */
  trans?: 'N' | 'C';

  /**
   * Complex scalar multiplier.
   * @default 1.0
   */
  alpha?: number;

  /**
   * Real scalar multiplier for C.
   * @default 0.0
   */
  beta?: number;

  /**
   * Existing Hermitian matrix C to accumulate into.
   */
  C?: Matrix;
}

/**
 * Result from Hermitian rank-2k update.
 */
export interface Her2kResult {
  /**
   * Result Hermitian matrix C.
   */
  C: Float64Array;

  /**
   * Matrix dimension.
   */
  n: number;

  /**
   * Whether the operation was successful.
   */
  success: boolean;

  /**
   * Status message.
   */
  message: string;
}

// ============================================================
// BLAS LEVEL 2 (Matrix-Vector)
// ============================================================

/**
 * Options for matrix-vector multiplication.
 */
export interface MatvecOptions {
  /**
   * Transpose operation on A.
   * @default 'N'
   */
  trans?: TransposeOp;

  /**
   * Scalar multiplier for A*x.
   * @default 1.0
   */
  alpha?: number;

  /**
   * Scalar multiplier for y.
   * @default 0.0
   */
  beta?: number;

  /**
   * Existing vector y to accumulate into.
   */
  y?: Vector;
}

/**
 * Result from matrix-vector multiplication.
 */
export interface MatvecResult {
  /**
   * Result vector y = alpha*op(A)*x + beta*y.
   */
  y: Float64Array;

  /**
   * Length of result vector.
   */
  n: number;

  /**
   * Whether the operation was successful.
   */
  success: boolean;

  /**
   * Status message.
   */
  message: string;
}

/**
 * Options for triangular matrix-vector multiplication.
 */
export interface MatvecTriangularOptions {
  /**
   * Whether A is upper triangular.
   * @default true
   */
  upper?: boolean;

  /**
   * Transpose operation on A.
   * @default 'N'
   */
  trans?: TransposeOp;

  /**
   * Whether A has unit diagonal.
   * @default false
   */
  unitDiagonal?: boolean;
}

/**
 * Result from triangular matrix-vector multiplication.
 */
export interface MatvecTriangularResult {
  /**
   * Result vector x = op(A)*x.
   */
  x: Float64Array;

  /**
   * Length of result vector.
   */
  n: number;

  /**
   * Whether the operation was successful.
   */
  success: boolean;

  /**
   * Status message.
   */
  message: string;
}

/**
 * Options for triangular vector solve.
 */
export interface SolveVectorTriangularOptions {
  /**
   * Whether A is upper triangular.
   * @default true
   */
  upper?: boolean;

  /**
   * Transpose operation on A.
   * @default 'N'
   */
  trans?: TransposeOp;

  /**
   * Whether A has unit diagonal.
   * @default false
   */
  unitDiagonal?: boolean;
}

/**
 * Result from triangular vector solve.
 */
export interface SolveVectorTriangularResult {
  /**
   * Solution vector x such that op(A)*x = b.
   */
  x: Float64Array;

  /**
   * Length of result vector.
   */
  n: number;

  /**
   * Whether the operation was successful.
   */
  success: boolean;

  /**
   * Status message.
   */
  message: string;
}

/**
 * Options for rank-1 update.
 */
export interface GerOptions {
  /**
   * Scalar multiplier.
   * @default 1.0
   */
  alpha?: number;

  /**
   * Existing matrix A to accumulate into.
   */
  A?: Matrix;
}

/**
 * Result from rank-1 update.
 */
export interface GerResult {
  /**
   * Result matrix A = alpha*x*y^T + A.
   */
  A: Float64Array;

  /**
   * Number of rows.
   */
  m: number;

  /**
   * Number of columns.
   */
  n: number;

  /**
   * Whether the operation was successful.
   */
  success: boolean;

  /**
   * Status message.
   */
  message: string;
}

/**
 * Options for symmetric rank-1 update.
 */
export interface SyrOptions {
  /**
   * Whether to update upper or lower triangle.
   * @default 'lower'
   */
  uplo?: 'upper' | 'lower';

  /**
   * Scalar multiplier.
   * @default 1.0
   */
  alpha?: number;

  /**
   * Existing symmetric matrix A to accumulate into.
   */
  A?: Matrix;
}

/**
 * Result from symmetric rank-1 update.
 */
export interface SyrResult {
  /**
   * Result symmetric matrix A = alpha*x*x^T + A.
   */
  A: Float64Array;

  /**
   * Matrix dimension.
   */
  n: number;

  /**
   * Whether the operation was successful.
   */
  success: boolean;

  /**
   * Status message.
   */
  message: string;
}

/**
 * Options for Hermitian rank-1 update.
 */
export interface HerOptions {
  /**
   * Whether to update upper or lower triangle.
   * @default 'lower'
   */
  uplo?: 'upper' | 'lower';

  /**
   * Real scalar multiplier.
   * @default 1.0
   */
  alpha?: number;

  /**
   * Existing Hermitian matrix A to accumulate into.
   */
  A?: Matrix;
}

/**
 * Result from Hermitian rank-1 update.
 */
export interface HerResult {
  /**
   * Result Hermitian matrix A = alpha*x*x^H + A.
   */
  A: Float64Array;

  /**
   * Matrix dimension.
   */
  n: number;

  /**
   * Whether the operation was successful.
   */
  success: boolean;

  /**
   * Status message.
   */
  message: string;
}

// ============================================================
// BLAS LEVEL 1 (Vector)
// ============================================================

/**
 * Result from dot product.
 */
export interface DotResult {
  /**
   * The dot product x^T * y.
   */
  dot: number;

  /**
   * Length of vectors.
   */
  n: number;

  /**
   * Whether the operation was successful.
   */
  success: boolean;

  /**
   * Status message.
   */
  message: string;
}

/**
 * Result from conjugate dot product.
 */
export interface DotcResult {
  /**
   * Real part of conjugate dot product x^H * y.
   */
  real: number;

  /**
   * Imaginary part of conjugate dot product.
   */
  imag: number;

  /**
   * Length of vectors.
   */
  n: number;

  /**
   * Whether the operation was successful.
   */
  success: boolean;

  /**
   * Status message.
   */
  message: string;
}

/**
 * Result from axpy operation.
 */
export interface AxpyResult {
  /**
   * Result vector y = alpha*x + y.
   */
  y: Float64Array;

  /**
   * Length of vectors.
   */
  n: number;

  /**
   * Whether the operation was successful.
   */
  success: boolean;

  /**
   * Status message.
   */
  message: string;
}

/**
 * Result from scal operation.
 */
export interface ScalResult {
  /**
   * Result vector x = alpha*x.
   */
  x: Float64Array;

  /**
   * Length of vector.
   */
  n: number;

  /**
   * Whether the operation was successful.
   */
  success: boolean;

  /**
   * Status message.
   */
  message: string;
}

/**
 * Result from copy operation.
 */
export interface CopyResult {
  /**
   * Destination vector y = x.
   */
  y: Float64Array;

  /**
   * Length of vectors.
   */
  n: number;

  /**
   * Whether the operation was successful.
   */
  success: boolean;

  /**
   * Status message.
   */
  message: string;
}

/**
 * Result from swap operation.
 */
export interface SwapResult {
  /**
   * First vector (swapped).
   */
  x: Float64Array;

  /**
   * Second vector (swapped).
   */
  y: Float64Array;

  /**
   * Length of vectors.
   */
  n: number;

  /**
   * Whether the operation was successful.
   */
  success: boolean;

  /**
   * Status message.
   */
  message: string;
}

/**
 * Result from nrm2 operation.
 */
export interface Nrm2Result {
  /**
   * Euclidean norm ||x||_2.
   */
  norm: number;

  /**
   * Length of vector.
   */
  n: number;

  /**
   * Whether the operation was successful.
   */
  success: boolean;

  /**
   * Status message.
   */
  message: string;
}

/**
 * Result from asum operation.
 */
export interface AsumResult {
  /**
   * Sum of absolute values.
   */
  asum: number;

  /**
   * Length of vector.
   */
  n: number;

  /**
   * Whether the operation was successful.
   */
  success: boolean;

  /**
   * Status message.
   */
  message: string;
}

/**
 * Result from iamax operation.
 */
export interface IamaxResult {
  /**
   * 0-based index of element with maximum absolute value.
   */
  index: number;

  /**
   * The maximum absolute value.
   */
  value: number;

  /**
   * Length of vector.
   */
  n: number;

  /**
   * Whether the operation was successful.
   */
  success: boolean;

  /**
   * Status message.
   */
  message: string;
}
