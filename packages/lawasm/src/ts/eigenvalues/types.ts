/**
 * Eigenvalue Problem Types
 *
 * Type definitions for eigenvalue computation functions.
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
 * Complex number representation.
 */
export interface Complex {
  re: number;
  im: number;
}

// ============================================================
// GENERAL EIGENVALUE PROBLEM
// ============================================================

/**
 * Options for general eigenvalue computation (eig).
 */
export interface EigOptions {
  /**
   * Whether to compute eigenvectors.
   * @default true
   */
  computeVectors?: boolean;

  /**
   * Whether to compute left eigenvectors (in addition to right).
   * @default false
   */
  computeLeftVectors?: boolean;

  /**
   * Whether to overwrite the input matrix A.
   * @default false
   */
  overwriteA?: boolean;
}

/**
 * Options for computing eigenvalues only (eigvals).
 */
export interface EigvalsOptions {
  /**
   * Whether to overwrite the input matrix A.
   * @default false
   */
  overwriteA?: boolean;
}

/**
 * Result from general eigenvalue computation.
 */
export interface EigResult {
  /**
   * Eigenvalues as complex numbers.
   * For a real matrix, complex eigenvalues come in conjugate pairs.
   */
  values: Complex[];

  /**
   * Real parts of eigenvalues.
   */
  wr: Float64Array;

  /**
   * Imaginary parts of eigenvalues.
   */
  wi: Float64Array;

  /**
   * Right eigenvectors (if computed).
   * For complex eigenvalue λ with eigenvector v = vr + i*vi:
   * - If wi[j] = 0, the j-th eigenvector is real: v[:,j]
   * - If wi[j] > 0, the j-th and (j+1)-th eigenvalues form a complex pair,
   *   with eigenvectors v[:,j] ± i*v[:,j+1]
   *
   * Stored in column-major order.
   */
  vectors?: Float64Array;

  /**
   * Left eigenvectors (if computed).
   */
  leftVectors?: Float64Array;

  /**
   * Matrix dimension.
   */
  n: number;

  /**
   * LAPACK info code.
   */
  info: number;

  /**
   * Whether the computation was successful.
   */
  success: boolean;

  /**
   * Status message.
   */
  message: string;
}

/**
 * Result from eigenvalues-only computation.
 */
export interface EigvalsResult {
  /**
   * Eigenvalues as complex numbers.
   */
  values: Complex[];

  /**
   * Real parts of eigenvalues.
   */
  wr: Float64Array;

  /**
   * Imaginary parts of eigenvalues.
   */
  wi: Float64Array;

  /**
   * Matrix dimension.
   */
  n: number;

  /**
   * LAPACK info code.
   */
  info: number;

  /**
   * Whether the computation was successful.
   */
  success: boolean;

  /**
   * Status message.
   */
  message: string;
}

// ============================================================
// SYMMETRIC EIGENVALUE PROBLEM
// ============================================================

/**
 * Options for symmetric eigenvalue computation.
 */
export interface EigSymmetricOptions {
  /**
   * Whether to compute eigenvectors.
   * @default true
   */
  computeVectors?: boolean;

  /**
   * Which algorithm to use:
   * - 'syev': Standard algorithm (DSYEV)
   * - 'syevd': Divide and conquer (DSYEVD) - faster for large matrices
   * - 'syevr': RRR algorithm (DSYEVR) - fastest, can compute subset
   * @default 'syevd'
   */
  algorithm?: 'syev' | 'syevd' | 'syevr';

  /**
   * Which triangle of A contains the data.
   * @default 'lower'
   */
  uplo?: 'upper' | 'lower';

  /**
   * Whether to overwrite the input matrix A.
   * @default false
   */
  overwriteA?: boolean;
}

/**
 * Result from symmetric eigenvalue computation.
 * All eigenvalues of a symmetric matrix are real.
 */
export interface EigSymmetricResult {
  /**
   * Eigenvalues in ascending order.
   */
  values: Float64Array;

  /**
   * Eigenvectors (if computed).
   * The j-th column is the eigenvector corresponding to values[j].
   * Stored in column-major order.
   */
  vectors?: Float64Array;

  /**
   * Matrix dimension.
   */
  n: number;

  /**
   * LAPACK info code.
   */
  info: number;

  /**
   * Whether the computation was successful.
   */
  success: boolean;

  /**
   * Status message.
   */
  message: string;
}

// ============================================================
// HERMITIAN EIGENVALUE PROBLEM
// ============================================================

/**
 * Options for Hermitian eigenvalue computation.
 */
export interface EigHermitianOptions {
  /**
   * Whether to compute eigenvectors.
   * @default true
   */
  computeVectors?: boolean;

  /**
   * Which algorithm to use:
   * - 'heev': Standard algorithm (ZHEEV)
   * - 'heevd': Divide and conquer (ZHEEVD) - faster for large matrices
   * - 'heevr': RRR algorithm (ZHEEVR) - fastest, can compute subset
   * @default 'heevd'
   */
  algorithm?: 'heev' | 'heevd' | 'heevr';

  /**
   * Which triangle of A contains the data.
   * @default 'lower'
   */
  uplo?: 'upper' | 'lower';

  /**
   * Whether to overwrite the input matrix A.
   * @default false
   */
  overwriteA?: boolean;
}

/**
 * Result from Hermitian eigenvalue computation.
 * All eigenvalues of a Hermitian matrix are real.
 */
export interface EigHermitianResult {
  /**
   * Eigenvalues in ascending order (always real for Hermitian matrices).
   */
  values: Float64Array;

  /**
   * Complex eigenvectors (if computed).
   * Stored as interleaved real/imaginary pairs in column-major order.
   * Length is 2*n*n.
   */
  vectors?: Float64Array;

  /**
   * Matrix dimension.
   */
  n: number;

  /**
   * LAPACK info code.
   */
  info: number;

  /**
   * Whether the computation was successful.
   */
  success: boolean;

  /**
   * Status message.
   */
  message: string;
}

// ============================================================
// GENERALIZED EIGENVALUE PROBLEM
// ============================================================

/**
 * Options for generalized eigenvalue problem Ax = λBx.
 */
export interface EigGeneralizedOptions {
  /**
   * Whether to compute eigenvectors.
   * @default true
   */
  computeVectors?: boolean;

  /**
   * Whether to compute left eigenvectors (in addition to right).
   * @default false
   */
  computeLeftVectors?: boolean;

  /**
   * Whether to overwrite input matrices.
   * @default false
   */
  overwriteInputs?: boolean;
}

/**
 * Result from generalized eigenvalue problem.
 * Generalized eigenvalues are given by alphar[j]/beta[j] + i*alphai[j]/beta[j].
 */
export interface EigGeneralizedResult {
  /**
   * Real parts of numerators of eigenvalues.
   */
  alphar: Float64Array;

  /**
   * Imaginary parts of numerators of eigenvalues.
   */
  alphai: Float64Array;

  /**
   * Denominators of eigenvalues.
   * The j-th eigenvalue is (alphar[j] + i*alphai[j]) / beta[j].
   * If beta[j] = 0, the eigenvalue is infinite.
   */
  beta: Float64Array;

  /**
   * Eigenvalues as complex numbers.
   * Computed from alphar, alphai, and beta.
   */
  values: Complex[];

  /**
   * Right eigenvectors (if computed).
   */
  vectors?: Float64Array;

  /**
   * Left eigenvectors (if computed).
   */
  leftVectors?: Float64Array;

  /**
   * Matrix dimension.
   */
  n: number;

  /**
   * LAPACK info code.
   */
  info: number;

  /**
   * Whether the computation was successful.
   */
  success: boolean;

  /**
   * Status message.
   */
  message: string;
}

// ============================================================
// SYMMETRIC GENERALIZED EIGENVALUE PROBLEM
// ============================================================

/**
 * Options for symmetric generalized eigenvalue problem.
 */
export interface EigGeneralizedSymmetricOptions {
  /**
   * Which problem to solve:
   * - 1: A*x = λ*B*x
   * - 2: A*B*x = λ*x
   * - 3: B*A*x = λ*x
   * @default 1
   */
  itype?: 1 | 2 | 3;

  /**
   * Whether to compute eigenvectors.
   * @default true
   */
  computeVectors?: boolean;

  /**
   * Which triangle of A and B contains the data.
   * @default 'lower'
   */
  uplo?: 'upper' | 'lower';

  /**
   * Whether to overwrite input matrices.
   * @default false
   */
  overwriteInputs?: boolean;
}

/**
 * Result from symmetric generalized eigenvalue problem.
 * All eigenvalues are real for symmetric/symmetric positive definite pair.
 */
export interface EigGeneralizedSymmetricResult {
  /**
   * Eigenvalues in ascending order.
   */
  values: Float64Array;

  /**
   * Eigenvectors (if computed).
   */
  vectors?: Float64Array;

  /**
   * Matrix dimension.
   */
  n: number;

  /**
   * LAPACK info code.
   */
  info: number;

  /**
   * Whether the computation was successful.
   */
  success: boolean;

  /**
   * Status message.
   */
  message: string;
}

// ============================================================
// BANDED EIGENVALUE PROBLEM
// ============================================================

/**
 * Options for banded matrix eigenvalue computation.
 */
export interface EigBandedOptions {
  /**
   * Whether to compute eigenvectors.
   * @default true
   */
  computeVectors?: boolean;

  /**
   * Which triangle of the banded matrix contains the data.
   * @default 'lower'
   */
  uplo?: 'upper' | 'lower';
}

/**
 * Result from banded eigenvalue computation.
 * For symmetric banded matrices, all eigenvalues are real.
 */
export interface EigBandedResult {
  /**
   * Eigenvalues in ascending order.
   */
  values: Float64Array;

  /**
   * Eigenvectors (if computed).
   */
  vectors?: Float64Array;

  /**
   * Matrix dimension.
   */
  n: number;

  /**
   * LAPACK info code.
   */
  info: number;

  /**
   * Whether the computation was successful.
   */
  success: boolean;

  /**
   * Status message.
   */
  message: string;
}

// ============================================================
// TRIDIAGONAL EIGENVALUE PROBLEM
// ============================================================

/**
 * Options for tridiagonal eigenvalue computation.
 */
export interface EigTridiagonalOptions {
  /**
   * Whether to compute eigenvectors.
   * @default true
   */
  computeVectors?: boolean;

  /**
   * Which algorithm to use:
   * - 'stev': Standard algorithm (DSTEV)
   * - 'stevd': Divide and conquer (DSTEVD) - faster
   * - 'stevr': RRR algorithm (DSTEVR) - fastest
   * @default 'stevd'
   */
  algorithm?: 'stev' | 'stevd' | 'stevr';
}

/**
 * Result from tridiagonal eigenvalue computation.
 * All eigenvalues of a symmetric tridiagonal matrix are real.
 */
export interface EigTridiagonalResult {
  /**
   * Eigenvalues in ascending order.
   */
  values: Float64Array;

  /**
   * Eigenvectors (if computed).
   */
  vectors?: Float64Array;

  /**
   * Matrix dimension.
   */
  n: number;

  /**
   * LAPACK info code.
   */
  info: number;

  /**
   * Whether the computation was successful.
   */
  success: boolean;

  /**
   * Status message.
   */
  message: string;
}

// ============================================================
// SELECTED EIGENVALUES
// ============================================================

/**
 * Options for computing selected eigenvalues.
 */
export interface EigSelectOptions {
  /**
   * Which eigenvalues to compute:
   * - 'all': All eigenvalues
   * - 'value': Eigenvalues in interval (vl, vu]
   * - 'index': Eigenvalues with indices il to iu (1-based)
   * @default 'all'
   */
  range?: 'all' | 'value' | 'index';

  /**
   * Lower bound of interval (for range='value').
   */
  vl?: number;

  /**
   * Upper bound of interval (for range='value').
   */
  vu?: number;

  /**
   * Index of smallest eigenvalue to compute (for range='index', 1-based).
   */
  il?: number;

  /**
   * Index of largest eigenvalue to compute (for range='index', 1-based).
   */
  iu?: number;

  /**
   * Whether to compute eigenvectors.
   * @default true
   */
  computeVectors?: boolean;

  /**
   * Which triangle of the matrix contains the data.
   * @default 'lower'
   */
  uplo?: 'upper' | 'lower';

  /**
   * Absolute tolerance for eigenvalue accuracy.
   * @default 0 (use default based on machine precision)
   */
  abstol?: number;
}

/**
 * Result from selected eigenvalue computation.
 */
export interface EigSelectResult {
  /**
   * Eigenvalues found (in ascending order).
   */
  values: Float64Array;

  /**
   * Number of eigenvalues found.
   */
  m: number;

  /**
   * Eigenvectors (if computed).
   * Only m columns are valid.
   */
  vectors?: Float64Array;

  /**
   * Matrix dimension.
   */
  n: number;

  /**
   * LAPACK info code.
   */
  info: number;

  /**
   * Whether the computation was successful.
   */
  success: boolean;

  /**
   * Status message.
   */
  message: string;
}
