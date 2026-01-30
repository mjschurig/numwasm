/**
 * SciWasm - SciPy-inspired scientific computing in TypeScript/WebAssembly
 */

export * as integrate from './integrate/index.js';
export * as optimize from './optimize/index.js';
export * as sparse from './sparse/index.js';
export * as spatial from './spatial/index.js';
export * as special from './special/index.js';
export * as stats from './stats/index.js';

// WASM module loaders
export { loadWasmModule, configureWasm, isWasmLoaded, resetWasmModule } from './wasm-loader.js';
export { loadSuperLUModule, configureSuperLU, isSuperLULoaded, resetSuperLUModule, getSuperLUModule } from './superlu-loader.js';
export { loadARPACKModule, configureARPACK, isARPACKLoaded, resetARPACKModule, getARPACKModule } from './arpack-loader.js';

// WASM module types
export type { SuperLUModule, SuperLUModuleFactory, SuperLUModuleOptions } from './superlu-types.js';
export type { ARPACKModule, ARPACKModuleFactory, ARPACKModuleOptions } from './arpack-types.js';

// Errors
export {
  NotImplementedError,
  IntegrationWarning,
  SingularMatrixError,
  DimensionMismatchError,
  ARPACKError,
  ARPACKNoConvergence,
} from './errors.js';
