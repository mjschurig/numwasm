/**
 * ODE WASM - ODE solvers compiled to WebAssembly
 *
 * Provides TypeScript bindings to classic Fortran ODE solvers:
 *
 * **Hairer Solvers:**
 * - DOPRI5: Dormand-Prince 5(4) explicit Runge-Kutta (RK45)
 * - DOP853: Dormand-Prince 8(5,3) high-order explicit Runge-Kutta
 * - RADAU5: Implicit Runge-Kutta for stiff systems
 *
 * **Netlib Solvers:**
 * - RKF45: Runge-Kutta-Fehlberg 4(5)
 * - DVERK: Verner 6(5) Runge-Kutta
 * - ODE: Adams-Bashforth-Moulton multistep
 * - VODE: Variable-coefficient ODE solver (BDF/Adams)
 * - ZVODE: VODE for complex-valued ODEs
 * - VODPK: VODE with Krylov methods for large systems
 * - RKSUITE: RK suite with (2,3), (4,5), and (7,8) pairs
 * - RKC: Runge-Kutta-Chebyshev for mildly stiff problems
 */

// Hairer solver APIs
export { dopri5 } from "./ts/dopri5.js";
export { dop853 } from "./ts/dop853.js";
export { radau5 } from "./ts/radau5.js";

// Netlib solver APIs
export { rkf45 } from "./ts/rkf45.js";
export { dverk } from "./ts/dverk.js";
export { ode_abm } from "./ts/ode.js";
export { vode } from "./ts/vode.js";
export { zvode } from "./ts/zvode.js";
export { vodpk } from "./ts/vodpk.js";
export { rksuite } from "./ts/rksuite.js";
export { rkc } from "./ts/rkc.js";

// Unified API (SciPy-compatible)
export { solve_ivp } from "./ts/solve_ivp.js";

// Constants and enums
export {
  STATUS_MESSAGES,
  ExplicitMethod,
  ImplicitMethod,
  StabilizedMethod,
} from "./ts/types.js";

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

// Solver-specific option types
export type { Rkf45Options } from "./ts/rkf45.js";
export type { DverkOptions } from "./ts/dverk.js";
export type { OdeOptions } from "./ts/ode.js";
export type { VodeOptions, VodeMethodFlag } from "./ts/vode.js";
export type {
  ZvodeOptions,
  ZvodeMethodFlag,
  ZvodeSolverResult,
  Complex,
  ComplexODEFunction,
} from "./ts/zvode.js";
export type { VodpkOptions, VodpkMethodFlag, PsolFunction } from "./ts/vodpk.js";
export type { RksuiteOptions, RksuiteMethod } from "./ts/rksuite.js";
export type { RkcOptions, SpcradFunction } from "./ts/rkc.js";
