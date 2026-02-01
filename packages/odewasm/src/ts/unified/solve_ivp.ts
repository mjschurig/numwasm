/**
 * Unified solve_ivp interface for ODE solving
 *
 * Provides a SciPy-like API that dispatches to dedicated solver implementations.
 */

import { dopri5 } from "../explicit/dopri5.js";
import { dop853 } from "../explicit/dop853.js";
import { odex } from "../explicit/odex.js";
import { radau5 } from "../implicit/radau5.js";
import { rkf45 } from "../explicit/rkf45.js";
import { dverk } from "../explicit/dverk.js";
import { ode_abm } from "../explicit/ode.js";
import { vode } from "../implicit/vode.js";
import { rksuite } from "../explicit/rksuite.js";
import { rkc } from "../stabilized/rkc.js";
import {
  ExplicitMethod,
  ODEFunction,
  SolveIVPOptions,
  SolveIVPResult,
} from "../core/types.js";

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
 *
 * // Using a different method
 * const rkcResult = await solve_ivp(
 *   (t, y) => y.map((yi, i) => -i * yi),
 *   [0, 1],
 *   [1, 1, 1],
 *   { method: 'RKC' }
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

  // Convert method to string for comparison
  const methodStr = String(method);

  switch (methodStr) {
    // Hairer explicit methods
    case "RK45":
      return dopri5(fun, t_span, y0, {
        rtol,
        atol,
        max_step,
        first_step,
        dense_output,
        t_eval,
      });

    case "DOP853":
      return dop853(fun, t_span, y0, {
        rtol,
        atol,
        max_step,
        first_step,
        dense_output,
        t_eval,
      });

    case "ODEX":
      return odex(fun, t_span, y0, {
        rtol,
        atol,
        max_step,
        first_step,
        dense_output,
        t_eval,
      });

    // Hairer implicit method
    case "Radau":
      return radau5(fun, t_span, y0, {
        rtol,
        atol,
        first_step,
        dense_output,
        t_eval,
        jac,
      });

    // Netlib explicit methods
    case "RKF45":
      return rkf45(fun, t_span, y0, {
        rtol,
        atol: typeof atol === "number" ? atol : atol[0],
      });

    case "DVERK":
      return dverk(fun, t_span, y0, {
        tol: rtol,
      });

    case "ODE":
      return ode_abm(fun, t_span, y0, {
        rtol,
        atol: typeof atol === "number" ? atol : atol[0],
      });

    // RKSUITE methods
    case "RKSUITE23":
      return rksuite(fun, t_span, y0, {
        tol: rtol,
        method: 1,
      });

    case "RKSUITE45":
      return rksuite(fun, t_span, y0, {
        tol: rtol,
        method: 2,
      });

    case "RKSUITE78":
      return rksuite(fun, t_span, y0, {
        tol: rtol,
        method: 3,
      });

    // Netlib implicit methods
    case "VODE":
      return vode(fun, t_span, y0, {
        rtol,
        atol,
        mf: 10, // Adams method
        jac,
      });

    case "BDF":
      return vode(fun, t_span, y0, {
        rtol,
        atol,
        mf: 22, // BDF with internal Jacobian
        jac,
      });

    // Stabilized explicit method
    case "RKC":
      return rkc(fun, t_span, y0, {
        rtol,
        atol,
      });

    default:
      throw new Error(`Unknown method: ${method}`);
  }
}
