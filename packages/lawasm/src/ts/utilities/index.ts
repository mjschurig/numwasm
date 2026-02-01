/**
 * Matrix Utilities & Properties
 *
 * This module provides utility functions for matrix operations
 * and property tests.
 *
 * @packageDocumentation
 */

// ============================================================
// Matrix Utilities
// ============================================================

export { transpose } from './transpose.js';
export { conjugate } from './conjugate.js';
export { hermitian } from './hermitian.js';
export { triu } from './triu.js';
export { tril } from './tril.js';
export { diag } from './diag.js';
export { trace } from './trace.js';
export { balance } from './balance.js';

// ============================================================
// Matrix Properties & Tests
// ============================================================

export { isSymmetric } from './is-symmetric.js';
export { isHermitian } from './is-hermitian.js';
export { isPositiveDefinite } from './is-positive-definite.js';
export { isOrthogonal } from './is-orthogonal.js';
export { isUnitary } from './is-unitary.js';
export { isSingular } from './is-singular.js';

// ============================================================
// Types
// ============================================================

export type {
  Matrix,
  ComplexMatrix,
  // Utility results
  TransposeResult,
  ConjugateResult,
  HermitianResult,
  TriuOptions,
  TriuResult,
  TrilOptions,
  TrilResult,
  DiagResult,
  TraceResult,
  BalanceResult,
  // Property test options and results
  IsSymmetricOptions,
  IsSymmetricResult,
  IsHermitianOptions,
  IsHermitianResult,
  IsPositiveDefiniteResult,
  IsOrthogonalOptions,
  IsOrthogonalResult,
  IsUnitaryOptions,
  IsUnitaryResult,
  IsSingularOptions,
  IsSingularResult,
} from './types.js';
