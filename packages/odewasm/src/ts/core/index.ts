/**
 * Core Module
 *
 * Provides module loading, type definitions, and shared utilities.
 */

export {
  loadODEModule,
  getODEModule,
  isODELoaded,
  resetODEModule,
  configureODE,
} from "./loader.js";

export {
  STATUS_MESSAGES,
  ExplicitMethod,
  ImplicitMethod,
  StabilizedMethod,
} from "./types.js";

export type {
  // Module types
  ODEModule,
  ODEModuleFactory,
  ODEModuleOptions,
  // Function types
  ODEFunction,
  JacobianFunction,
  EventFunction,
  SoloutCallback,
  DenseOutput,
  // Method enums
  ODEMethod,
  // Solver-specific options
  Dopri5Options,
  Dop853Options,
  OdexOptions,
  Radau5Options,
  ODESolverResult,
  // Unified options
  SolveIVPOptions,
  SolveIVPResult,
} from "./types.js";
