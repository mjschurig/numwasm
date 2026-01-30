/**
 * Unified solve_ivp interface for ODE solving
 *
 * Provides a SciPy-like API that dispatches to dedicated solver implementations.
 */

import { dopri5 } from "./dopri5.js";
import { dop853 } from "./dop853.js";
import { radau5 } from "./radau5.js";
import {
  ExplicitMethod,
  ODEFunction,
  SolveIVPOptions,
  SolveIVPResult,
} from "./types.js";

/**
 * Solve an initial value problem for a system of ODEs.
 *
 * Integrates a system of ordinary differential equations:
 *   dy/dt = f(t, y)
 *   y(t0) = y0
 *
 * @param fun - ODE function: f(t, y) => dy/dt
 * @param t_span - Time span [t0, tf]
 * @param y0 - Initial state vector
 * @param options - Solver options
 * @returns Solution result
 *
 * @example
 * ```ts
 * import { solve_ivp } from 'odewasm';
 *
 * // Simple exponential decay: dy/dt = -y
 * const result = await solve_ivp(
 *   (t, y) => [-y[0]],
 *   [0, 5],
 *   [1]
 * );
 *
 * // Harmonic oscillator: y'' = -y => [y, y'] => [y', -y]
 * const harmonic = await solve_ivp(
 *   (t, [y, v]) => [v, -y],
 *   [0, 10],
 *   [1, 0]  // y(0) = 1, y'(0) = 0
 * );
 * ```
 */
export async function solve_ivp(
  fun: ODEFunction,
  t_span: [number, number],
  y0: number[],
  options?: SolveIVPOptions,
): Promise<SolveIVPResult> {
  const {
    method = ExplicitMethod.RK45,
    rtol = 1e-3,
    atol = 1e-6,
    max_step = Infinity,
    first_step = 0,
    dense_output = false,
    t_eval,
    jac,
  } = options ?? {};

  if (method === "RK45") {
    return dopri5(fun, t_span, y0, {
      rtol,
      atol,
      max_step,
      first_step,
      dense_output,
      t_eval,
    });
  } else if (method === "DOP853") {
    return dop853(fun, t_span, y0, {
      rtol,
      atol,
      max_step,
      first_step,
      dense_output,
      t_eval,
    });
  } else if (method === "Radau") {
    return radau5(fun, t_span, y0, {
      rtol,
      atol,
      first_step,
      dense_output,
      t_eval,
      jac,
    });
  } else {
    throw new Error(`Unknown method: ${method}`);
  }
}
