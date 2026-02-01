/**
 * Compute L-infinity (maximum) norm of grid function.
 *
 * @module gridfunction
 */

import { GridFunction } from './GridFunction.js';

/**
 * Computes the L-infinity (maximum) norm of a grid function.
 *
 * The L-infinity norm is defined as:
 * ```
 * ||u||_∞ = max |u(x)| over all x in Ω
 * ```
 *
 * In practice, this is computed as the maximum absolute value of the DOF coefficients.
 *
 * @param gf - GridFunction to compute norm of
 * @returns The L-infinity norm value
 * @throws Error if the grid function has been destroyed
 * @throws Error if norm computation fails
 *
 * @category GridFunction
 *
 * @example Check maximum value
 * ```typescript
 * const maxVal = normLinf(solution);
 * console.log(`Maximum absolute value: ${maxVal}`);
 * ```
 *
 * @example Check for blow-up in time-dependent problems
 * ```typescript
 * const maxNorm = normLinf(solution);
 * if (maxNorm > 1e10) {
 *   console.error('Solution is blowing up!');
 * }
 * ```
 *
 * @see {@link normL2} for integral norm
 * @see {@link normH1} for H1 norm including gradients
 */
export function normLinf(gf: GridFunction): number {
  const module = gf.getModule();
  const ptr = gf.getPointer();

  const norm = module._mfem_gridfunc_norm_linf(ptr);

  if (norm < 0) {
    throw new Error('Failed to compute L-infinity norm');
  }

  return norm;
}
