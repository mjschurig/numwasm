/**
 * QuadWasm - QUADPACK WebAssembly Module
 *
 * This package provides a WebAssembly build of QUADPACK (Quadrature PACKage),
 * a Fortran library for numerical integration of one-dimensional functions.
 *
 * ## High-Level API
 *
 * For most use cases, use the high-level functions that handle all the
 * complexity of memory management and callback registration:
 *
 * @example
 * ```typescript
 * import { quad, quad_inf, integrate } from 'quadwasm';
 *
 * // Simple definite integral
 * const result = await quad(x => x * x, 0, 1);
 * console.log(result.result);  // ~0.333...
 *
 * // Integral over infinite interval
 * const result2 = await quad_inf(x => Math.exp(-x * x), 0, Infinity);
 * console.log(result2.result);  // ~0.886... (√π/2)
 *
 * // Unified API with auto-detection
 * const result3 = await integrate(x => Math.sin(x), 0, Math.PI);
 * console.log(result3.result);  // ~2.0
 * ```
 *
 * ## Low-Level API
 *
 * For advanced use cases, you can use the raw WASM module interface:
 *
 * @example
 * ```typescript
 * import { loadQUADPACKModule, getQUADPACKModule } from 'quadwasm';
 *
 * await loadQUADPACKModule();
 * const quadpack = getQUADPACKModule();
 * // Use QUADPACK routines directly via the module interface
 * ```
 *
 * @packageDocumentation
 */

// ============================================================
// HIGH-LEVEL API (recommended for most users)
// ============================================================

// High-level integration functions
export { quad } from './ts/quad.js';
export { quad_inf } from './ts/quad_inf.js';
export { quad_osc, quad_fourier } from './ts/quad_osc.js';
export { quad_singular, quad_cauchy } from './ts/quad_singular.js';
export { quad_break } from './ts/quad_break.js';
export { quad_rule } from './ts/quad_rule.js';
export { quad_ng } from './ts/quad_ng.js';
export { quad_fixed, type GKRule } from './ts/quad_fixed.js';
export { integrate } from './ts/integrate.js';

// High-level type definitions
export type {
  IntegrandFunction,
  QuadOptionsBase,
  QuadOptions,
  QuadRuleOptions,
  QuadInfOptions,
  QuadOscOptions,
  QuadFourierOptions,
  QuadSingularOptions,
  QuadCauchyOptions,
  QuadBreakOptions,
  QuadNgOptions,
  QuadResult,
  QuadFixedResult,
  IntegrationType,
  IntegrateOptions,
} from './ts/high-level-types.js';

// Helper utilities
export {
  quadWorkSize,
  quadIworkSize,
  quadOscWorkSize,
  quadOscIworkSize,
  quadFourierWorkSize,
  quadFourierIworkSize,
  quadBreakWorkSize,
  quadBreakIworkSize,
  getQuadpackMessage,
} from './ts/helpers.js';

// ============================================================
// LOW-LEVEL API (for advanced usage)
// ============================================================

// Module loader functions
export {
  loadQUADPACKModule,
  getQUADPACKModule,
  isQUADPACKLoaded,
  resetQUADPACKModule,
  configureQUADPACK,
  type QUADPACKLoadConfig,
} from './ts/loader.js';

// Low-level module types
export type {
  QUADPACKModule,
  QUADPACKModuleFactory,
  QUADPACKModuleOptions,
} from './ts/types.js';

// Error codes and constants
export {
  QUADPACK_ERRORS,
  QUADPACKRule,
  QUADPACKInf,
  QUADPACKWeight,
  QUADPACKOscillatory,
} from './ts/types.js';
