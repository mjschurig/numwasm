/**
 * QuadWasm - QUADPACK WebAssembly Module
 *
 * This package provides a WebAssembly build of QUADPACK (Quadrature PACKage),
 * a Fortran library for numerical integration of one-dimensional functions.
 *
 * QUADPACK provides routines for:
 * - Automatic integration with adaptive subdivision (QAGS, QAG)
 * - Integration over infinite intervals (QAGI)
 * - Integration with singularities (QAGP, QAWS, QAWC)
 * - Oscillatory integration (QAWO, QAWF)
 *
 * @example
 * ```typescript
 * import { loadQUADPACKModule, getQUADPACKModule } from 'quadwasm';
 *
 * // Load the WASM module
 * await loadQUADPACKModule();
 *
 * // Get the module instance
 * const quadpack = getQUADPACKModule();
 *
 * // Use QUADPACK routines via the module interface
 * ```
 *
 * @packageDocumentation
 */

// Re-export loader functions
export {
  loadQUADPACKModule,
  getQUADPACKModule,
  isQUADPACKLoaded,
  resetQUADPACKModule,
  configureQUADPACK,
  type QUADPACKLoadConfig,
} from './ts/loader.js';

// Re-export types
export type {
  QUADPACKModule,
  QUADPACKModuleFactory,
} from './ts/types.js';
