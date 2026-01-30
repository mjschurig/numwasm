/**
 * LinWasm - LINPACK WebAssembly Module
 *
 * This package provides a WebAssembly build of LINPACK, a classic
 * Fortran library for solving linear algebra problems.
 *
 * LINPACK provides routines for:
 * - LU factorization (general matrices)
 * - Cholesky factorization (positive definite matrices)
 * - QR decomposition
 * - Singular value decomposition
 * - Linear system solving
 * - Matrix condition estimation
 *
 * @example
 * ```typescript
 * import { loadLINPACKModule, getLINPACKModule } from 'linwasm';
 *
 * // Load the WASM module
 * await loadLINPACKModule();
 *
 * // Get the module instance
 * const linpack = getLINPACKModule();
 *
 * // Use LINPACK routines via the module interface
 * ```
 *
 * @packageDocumentation
 */

// Re-export loader functions
export {
  loadLINPACKModule,
  getLINPACKModule,
  isLINPACKLoaded,
  resetLINPACKModule,
  configureLINPACK,
  type LINPACKLoadConfig,
} from './ts/loader.js';

// Re-export types
export type {
  LINPACKModule,
  LINPACKModuleFactory,
} from './ts/types.js';
