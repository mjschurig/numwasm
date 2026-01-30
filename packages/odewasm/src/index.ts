/**
 * ODE WASM - ODE solvers compiled to WebAssembly
 *
 * Provides TypeScript bindings to classic Fortran ODE solvers:
 * - DOPRI5: Dormand-Prince 5(4) explicit Runge-Kutta (RK45)
 * - DOP853: Dormand-Prince 8(5,3) high-order explicit Runge-Kutta
 * - RADAU5: Implicit Runge-Kutta for stiff systems
 */

// Dedicated solver APIs
export { dopri5 } from "./ts/dopri5.js";
export { dop853 } from "./ts/dop853.js";
export { radau5 } from "./ts/radau5.js";

// Unified API (SciPy-compatible)
export { solve_ivp } from "./ts/solve_ivp.js";

// Constants and enums
export { STATUS_MESSAGES, ExplicitMethod, ImplicitMethod } from "./ts/types.js";

// Module management
export {
  loadODEModule,
  getODEModule,
  isODELoaded,
  resetODEModule,
  configureODE,
} from "./ts/loader.js";

// Types
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
  Radau5Options,
  ODESolverResult,
  // Legacy unified options
  SolveIVPOptions,
  SolveIVPResult,
} from "./ts/types.js";
