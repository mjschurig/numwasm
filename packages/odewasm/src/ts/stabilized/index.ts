/**
 * Stabilized ODE Solvers (Mildly Stiff)
 *
 * These solvers are explicit but use stabilization techniques to handle
 * mildly stiff problems, particularly those arising from parabolic PDEs.
 */

export { rkc } from "./rkc.js";

// Re-export option types
export type { RkcOptions, SpcradFunction } from "./rkc.js";
