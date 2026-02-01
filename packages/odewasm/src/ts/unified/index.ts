/**
 * Unified ODE Solver API
 *
 * Provides a SciPy-compatible interface for solving initial value problems.
 * The solve_ivp function automatically dispatches to the appropriate solver
 * based on the specified method.
 */

export { solve_ivp } from "./solve_ivp.js";
