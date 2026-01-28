/**
 * Main minimize() dispatcher — unified interface for optimization.
 *
 * Mirrors scipy.optimize.minimize.
 */

import type {
  ObjectiveFunction,
  MinimizeOptions,
  OptimizeResult,
  Bounds,
  NelderMeadOptions,
  BFGSOptions,
  LBFGSBOptions,
} from './types.js';
import { minimizeNelderMead } from './nelder_mead.js';
import { minimizeBFGS } from './bfgs.js';
import { minimizeLBFGSB } from './lbfgsb.js';

/**
 * Normalize bounds from tuple-array format to Bounds object.
 * Accepts either a Bounds object or an array of [lower, upper] pairs.
 */
function normalizeBounds(
  bounds?: Bounds | Array<[number | null, number | null]>
): Bounds | undefined {
  if (!bounds) return undefined;
  if ('lb' in bounds && 'ub' in bounds) return bounds;

  const tuples = bounds as Array<[number | null, number | null]>;
  const lb: number[] = new Array(tuples.length);
  const ub: number[] = new Array(tuples.length);
  for (let i = 0; i < tuples.length; i++) {
    lb[i] = tuples[i][0] ?? -Infinity;
    ub[i] = tuples[i][1] ?? Infinity;
  }
  return { lb, ub };
}

/**
 * Minimization of scalar function of one or more variables.
 *
 * @param fun - Objective function: (x: number[]) => number
 * @param x0 - Initial guess
 * @param options - Method, bounds, tolerances, etc.
 * @returns OptimizeResult with solution, status, and diagnostics
 *
 * @example
 * ```ts
 * import { optimize } from 'sciwasm';
 *
 * // Rosenbrock function
 * const rosen = (x: number[]) =>
 *   100 * (x[1] - x[0] ** 2) ** 2 + (1 - x[0]) ** 2;
 *
 * const result = optimize.minimize(rosen, [-1, -1], {
 *   method: 'Nelder-Mead',
 * });
 * console.log(result.x); // ~[1, 1]
 * ```
 */
export async function minimize(
  fun: ObjectiveFunction,
  x0: number[],
  options?: MinimizeOptions
): Promise<OptimizeResult> {
  const bounds = normalizeBounds(options?.bounds);

  // Determine method — mirrors scipy defaults:
  //   no bounds, no constraints -> BFGS
  //   bounds, no constraints -> L-BFGS-B
  //   constraints -> SLSQP
  let method = options?.method;
  if (!method) {
    if (bounds) {
      method = 'L-BFGS-B';
    } else {
      method = 'BFGS';
    }
  }

  const methodLower = method.toLowerCase();

  switch (methodLower) {
    case 'nelder-mead':
      return await minimizeNelderMead(
        fun,
        x0,
        options?.options as NelderMeadOptions | undefined,
        bounds
      );

    case 'bfgs':
      return await minimizeBFGS(
        fun,
        x0,
        options?.options as BFGSOptions | undefined,
        typeof options?.jac === 'function' ? options.jac : undefined
      );

    case 'l-bfgs-b':
      return await minimizeLBFGSB(
        fun,
        x0,
        options?.options as LBFGSBOptions | undefined,
        bounds,
        typeof options?.jac === 'function' ? options.jac : undefined
      );

    case 'powell':
    case 'cg':
    case 'newton-cg':
    case 'tnc':
    case 'slsqp':
      throw new Error(
        `Method '${method}' is not yet implemented. Available methods: Nelder-Mead, BFGS, L-BFGS-B`
      );

    default:
      throw new Error(
        `Unknown method '${method}'. Available methods: Nelder-Mead, BFGS, L-BFGS-B, Powell, CG, Newton-CG, TNC, SLSQP`
      );
  }
}
