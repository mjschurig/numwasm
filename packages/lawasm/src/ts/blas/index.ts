/**
 * BLAS (Basic Linear Algebra Subprograms)
 *
 * This module provides functions for basic linear algebra operations:
 * - Level 3: Matrix-matrix operations
 * - Level 2: Matrix-vector operations
 * - Level 1: Vector operations
 *
 * @packageDocumentation
 */

// ============================================================
// BLAS Level 3 (Matrix-Matrix)
// ============================================================

export { matmul } from './matmul.js';
export { matmulTriangular } from './matmul-triangular.js';
export { solveMatrixTriangular } from './solve-matrix-triangular.js';
export { syrk } from './syrk.js';
export { syr2k } from './syr2k.js';
export { herk } from './herk.js';
export { her2k } from './her2k.js';

// ============================================================
// BLAS Level 2 (Matrix-Vector)
// ============================================================

export { matvec } from './matvec.js';
export { matvecTriangular } from './matvec-triangular.js';
export { solveVectorTriangular } from './solve-vector-triangular.js';
export { ger } from './ger.js';
export { syr } from './syr.js';
export { her } from './her.js';

// ============================================================
// BLAS Level 1 (Vector)
// ============================================================

export { dot } from './dot.js';
export { dotc } from './dotc.js';
export { axpy } from './axpy.js';
export { scal } from './scal.js';
export { copy } from './copy.js';
export { swap } from './swap.js';
export { nrm2 } from './nrm2.js';
export { asum } from './asum.js';
export { iamax } from './iamax.js';

// ============================================================
// Types
// ============================================================

export type {
  Matrix,
  Vector,
  TransposeOp,
  // Level 3
  MatmulOptions,
  MatmulResult,
  MatmulTriangularOptions,
  MatmulTriangularResult,
  SolveMatrixTriangularOptions,
  SolveMatrixTriangularResult,
  SyrkOptions,
  SyrkResult,
  Syr2kOptions,
  Syr2kResult,
  HerkOptions,
  HerkResult,
  Her2kOptions,
  Her2kResult,
  // Level 2
  MatvecOptions,
  MatvecResult,
  MatvecTriangularOptions,
  MatvecTriangularResult,
  SolveVectorTriangularOptions,
  SolveVectorTriangularResult,
  GerOptions,
  GerResult,
  SyrOptions,
  SyrResult,
  HerOptions,
  HerResult,
  // Level 1
  DotResult,
  DotcResult,
  AxpyResult,
  ScalResult,
  CopyResult,
  SwapResult,
  Nrm2Result,
  AsumResult,
  IamaxResult,
} from './types.js';
