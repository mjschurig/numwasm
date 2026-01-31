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
export { Mesh, FiniteElementSpace, GridFunction } from './ts/mesh.js';

export type { BoundingBox } from './ts/mesh.js';
