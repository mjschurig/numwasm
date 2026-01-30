/**
 * XSFWasm - XSF WebAssembly Module
 *
 * This package provides a WebAssembly build of XSF (Extended Special Functions),
 * a C++ library that implements special mathematical functions derived from
 * SciPy's Cephes, AMOS, and other sources.
 *
 * XSF provides routines for:
 * - Gamma and related functions (gamma, lgamma, digamma, beta)
 * - Bessel functions (J, Y, I, K variants)
 * - Airy functions (Ai, Bi and derivatives)
 * - Error functions (erf, erfc)
 * - Elliptic integrals
 * - Hypergeometric functions
 * - Statistical distributions
 *
 * @example
 * ```typescript
 * import { loadXSFModule, getXSFModule } from 'xsfwasm';
 *
 * // Load the WASM module
 * await loadXSFModule();
 *
 * // Get the module instance
 * const xsf = getXSFModule();
 *
 * // Use XSF routines via the module interface
 * const gammaValue = xsf._gamma(5.0); // Returns 24.0
 * ```
 *
 * @packageDocumentation
 */

// Re-export loader functions
export {
  loadXSFModule,
  getXSFModule,
  isXSFLoaded,
  resetXSFModule,
  configureXSF,
  type XSFLoadConfig,
} from './ts/loader.js';

// Re-export types
export type {
  XSFModule,
  XSFModuleFactory,
} from './ts/types.js';
