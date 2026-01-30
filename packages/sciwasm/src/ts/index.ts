/**
 * SciWasm - SciPy-inspired scientific computing in TypeScript/WebAssembly
 */

export * as integrate from './integrate/index.js';
export * as optimize from './optimize/index.js';
export * as sparse from './sparse/index.js';
export * as spatial from './spatial/index.js';
export * as special from './special/index.js';
export * as stats from './stats/index.js';

// WASM module loaders - re-exported from separate packages
export {
  loadSuperLUModule,
  configureSuperLU,
  isSuperLULoaded,
  resetSuperLUModule,
  getSuperLUModule,
  type SuperLUModule,
  type SuperLUModuleFactory,
  type SuperLULoadConfig,
} from 'superluwasm';

export {
  loadARPACKModule,
  configureARPACK,
  isARPACKLoaded,
  resetARPACKModule,
  getARPACKModule,
  type ARPACKModule,
  type ARPACKModuleFactory,
  type ARPACKLoadConfig,
} from 'arwasm';

export {
  loadXSFModule,
  configureXSF,
  isXSFLoaded,
  resetXSFModule,
  getXSFModule,
  type XSFModule,
  type XSFModuleFactory,
  type XSFLoadConfig,
} from 'xsfwasm';

export {
  loadQUADPACKModule,
  configureQUADPACK,
  isQUADPACKLoaded,
  resetQUADPACKModule,
  getQUADPACKModule,
  type QUADPACKModule,
  type QUADPACKModuleFactory,
  type QUADPACKLoadConfig,
} from 'quadwasm';

export {
  loadLAPACKModule,
  configureLAPACK,
  isLAPACKLoaded,
  resetLAPACKModule,
  getLAPACKModule,
  type LAPACKModule,
  type LAPACKModuleFactory,
  type LAPACKLoadConfig,
} from 'lawasm';

export {
  loadLINPACKModule,
  configureLINPACK,
  isLINPACKLoaded,
  resetLINPACKModule,
  getLINPACKModule,
  type LINPACKModule,
  type LINPACKModuleFactory,
  type LINPACKLoadConfig,
} from 'linwasm';

// Internal WASM module loader (sparsetools, spatial, optimize)
export { loadWasmModule, configureWasm, isWasmLoaded, resetWasmModule } from './wasm-loader.js';

// Errors
export {
  NotImplementedError,
  IntegrationWarning,
  SingularMatrixError,
  DimensionMismatchError,
  ARPACKError,
  ARPACKNoConvergence,
} from './errors.js';
