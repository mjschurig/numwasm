/**
 * Implicit ODE Solvers (Stiff)
 *
 * These solvers are designed for stiff problems where explicit methods
 * would require impractically small step sizes. They use implicit methods
 * with Newton iteration and Jacobian computation.
 */

export { radau5 } from "./radau5.js";
export { vode } from "./vode.js";
export { zvode } from "./zvode.js";
export { vodpk } from "./vodpk.js";

// Re-export option types
export type { VodeOptions, VodeMethodFlag } from "./vode.js";
export type {
  ZvodeOptions,
  ZvodeMethodFlag,
  ZvodeSolverResult,
  Complex,
  ComplexODEFunction,
} from "./zvode.js";
export type { VodpkOptions, VodpkMethodFlag, PsolFunction } from "./vodpk.js";
