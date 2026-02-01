/**
 * MFEM WebAssembly Module
 *
 * A finite element library for solving partial differential equations,
 * compiled to WebAssembly for use in browsers and Node.js.
 *
 * @packageDocumentation
 */

// Types
export type {
  MFEMModule,
  MFEMModuleFactory,
  MFEMModuleOptions,
} from './ts/types.js';

// Loader functions
export {
  loadMFEMModule,
  getMFEMModule,
  isMFEMLoaded,
  resetMFEMModule,
  configureMFEM,
} from './ts/loader.js';

export type { MFEMLoadConfig } from './ts/loader.js';

// High-level classes
export { Mesh, BoundingBox } from './ts/mesh/index.js';
export { FiniteElementSpace } from './ts/fespace/index.js';
export { GridFunction } from './ts/gridfunction/index.js';

// Mesh operations
export * as mesh from './ts/mesh/index.js';

// Finite element space operations
export * as fespace from './ts/fespace/index.js';

// Grid function operations
export * as gridfunction from './ts/gridfunction/index.js';
