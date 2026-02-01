/**
 * ODE WASM - ODE solvers compiled to WebAssembly
 *
 * Provides TypeScript bindings to classic Fortran ODE solvers:
 *
 * **Hairer Solvers:**
 * - DOPRI5: Dormand-Prince 5(4) explicit Runge-Kutta (RK45)
 * - DOP853: Dormand-Prince 8(5,3) high-order explicit Runge-Kutta
 * - ODEX: GBS extrapolation based on explicit midpoint rule
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

// Explicit solver APIs (non-stiff)
export { dopri5 } from "./ts/explicit/dopri5.js";
export { dop853 } from "./ts/explicit/dop853.js";
export { odex } from "./ts/explicit/odex.js";
export { rkf45 } from "./ts/explicit/rkf45.js";
export { dverk } from "./ts/explicit/dverk.js";
export { ode_abm } from "./ts/explicit/ode.js";
export { rksuite } from "./ts/explicit/rksuite.js";

// Implicit solver APIs (stiff)
export { radau5 } from "./ts/implicit/radau5.js";
export { vode } from "./ts/implicit/vode.js";
export { zvode } from "./ts/implicit/zvode.js";
export { vodpk } from "./ts/implicit/vodpk.js";

// Stabilized solver APIs (mildly stiff)
export { rkc } from "./ts/stabilized/rkc.js";

// Unified API (SciPy-compatible)
export { solve_ivp } from "./ts/unified/solve_ivp.js";

// Constants and enums
export {
  STATUS_MESSAGES,
  ExplicitMethod,
  ImplicitMethod,
  StabilizedMethod,
} from "./ts/core/types.js";

// Module management
export {
  loadODEModule,
  getODEModule,
  isODELoaded,
  resetODEModule,
  configureODE,
} from "./ts/core/loader.js";

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
  OdexOptions,
  Radau5Options,
  ODESolverResult,
  // Legacy unified options
  SolveIVPOptions,
  SolveIVPResult,
} from "./ts/core/types.js";

// Explicit solver option types
export type { Rkf45Options } from "./ts/explicit/rkf45.js";
export type { DverkOptions } from "./ts/explicit/dverk.js";
export type { OdeOptions } from "./ts/explicit/ode.js";
export type { RksuiteOptions, RksuiteMethod } from "./ts/explicit/rksuite.js";

// Implicit solver option types
export type { VodeOptions, VodeMethodFlag } from "./ts/implicit/vode.js";
export type {
  ZvodeOptions,
  ZvodeMethodFlag,
  ZvodeSolverResult,
  Complex,
  ComplexODEFunction,
} from "./ts/implicit/zvode.js";
export type { VodpkOptions, VodpkMethodFlag, PsolFunction } from "./ts/implicit/vodpk.js";

// Stabilized solver option types
export type { RkcOptions, SpcradFunction } from "./ts/stabilized/rkc.js";
