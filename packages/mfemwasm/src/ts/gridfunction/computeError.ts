/**
 * Compute error compared to exact solution.
 *
 * @module gridfunction
 */

import { GridFunction } from './GridFunction.js';

/**
 * Computes the L2 error between a grid function and a constant exact solution.
 *
 * The L2 error is defined as:
 * ```
 * ||u - u_exact||_L2 = sqrt( ∫_Ω |u(x) - u_exact|² dx )
 * ```
 *
 * This is useful for verifying solutions against known manufactured solutions
 * or for convergence studies.
 *
 * @param gf - GridFunction (numerical solution)
 * @param exactValue - The constant exact solution value
 * @returns The L2 error value
 * @throws Error if the grid function has been destroyed
 * @throws Error if error computation fails
 *
 * @category GridFunction
 *
 * @example Verify manufactured solution
 * ```typescript
 * // If exact solution is constant 1.0
 * const error = computeError(numerical, 1.0);
 * console.log(`L2 error: ${error}`);
 * ```
 *
 * @example Convergence study
 * ```typescript
 * const errors: number[] = [];
 * for (const h of meshSizes) {
 *   const solution = await solve(h);
 *   errors.push(computeError(solution, exactValue));
 *   solution.destroy();
 * }
 * // Compute convergence rate
 * for (let i = 1; i < errors.length; i++) {
 *   const rate = Math.log(errors[i-1] / errors[i]) / Math.log(2);
 *   console.log(`Convergence rate: ${rate.toFixed(2)}`);
 * }
 * ```
 *
 * @example Check if solution is close to zero
 * ```typescript
 * const error = computeError(residual, 0.0);
 * if (error < tolerance) {
 *   console.log('Converged!');
 * }
 * ```
 */
export function computeError(gf: GridFunction, exactValue: number): number {
  const module = gf.getModule();
  const ptr = gf.getPointer();

  const error = module._mfem_gridfunc_compute_l2_error_const(ptr, exactValue);

  if (error < 0) {
    throw new Error('Failed to compute L2 error');
  }

  return error;
}
