/**
 * Matrix Functions & Specialized Decompositions
 *
 * This module provides matrix functions (expm, logm, sqrtm, etc.)
 * and specialized decompositions (polar, RRQR, CSD).
 *
 * @packageDocumentation
 */

// ============================================================
// Matrix Functions
// ============================================================

export { expm } from './expm.js';
export { logm } from './logm.js';
export { sqrtm } from './sqrtm.js';
export { powm } from './powm.js';
export { funm } from './funm.js';

// ============================================================
// Specialized Decompositions
// ============================================================

export { polarDecomposition } from './polar-decomposition.js';
export { rrqr } from './rrqr.js';
export { csd } from './csd.js';

// ============================================================
// Complex-Specific Functions
// ============================================================

export { choleskyHermitian } from './cholesky-hermitian.js';
export type { CholeskyHermitianOptions } from './cholesky-hermitian.js';

// Re-export existing complex functions for convenience
// (eigHermitian is in eigenvalues, solveHermitian is in linear-solvers)

// ============================================================
// Types
// ============================================================

export type {
  Matrix,
  // Matrix Functions
  ExpmResult,
  LogmResult,
  SqrtmResult,
  PowmResult,
  FunmResult,
  // Specialized Decompositions
  PolarDecompositionResult,
  RRQROptions,
  RRQRResult,
  CSDResult,
  // Complex-Specific
  Complex,
  EigHermitianResult,
  CholeskyHermitianResult,
  SolveHermitianResult,
} from './types.js';
