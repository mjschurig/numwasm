/**
 * ARWasm - ARPACK WebAssembly Module
 *
 * This package provides a WebAssembly build of ARPACK (ARnoldi PACKage),
 * a collection of Fortran77 subroutines designed to solve large scale
 * eigenvalue problems.
 *
 * ARPACK provides routines for:
 * - Real symmetric eigenvalue problems (dsaupd/dseupd, ssaupd/sseupd)
 * - Real nonsymmetric eigenvalue problems (dnaupd/dneupd, snaupd/sneupd)
 * - Complex eigenvalue problems (cnaupd/cneupd, znaupd/zneupd)
 *
 * @example
 * ```typescript
 * import { loadARPACKModule, getARPACKModule } from 'arwasm';
 *
 * // Load the WASM module
 * await loadARPACKModule();
 *
 * // Get the module instance
 * const arpack = getARPACKModule();
 *
 * // Use ARPACK routines via the module interface
 * ```
 *
 * @packageDocumentation
 */

// Re-export loader functions
export {
  loadARPACKModule,
  getARPACKModule,
  isARPACKLoaded,
  resetARPACKModule,
  configureARPACK,
  type ARPACKLoadConfig,
} from './ts/loader.js';

// Re-export types
export type {
  ARPACKModule,
  ARPACKModuleFactory,
} from './ts/types.js';

// Re-export error code mappings
export {
  DSAUPD_ERRORS,
  DSEUPD_ERRORS,
  DNAUPD_ERRORS,
  DNEUPD_ERRORS,
} from './ts/types.js';
