/**
 * Compute H1 norm of grid function.
 *
 * @module gridfunction
 */

import { GridFunction } from './GridFunction.js';

/**
 * Computes the H1 norm of a grid function.
 *
 * The H1 norm is defined as:
 * ```
 * ||u||_H1 = sqrt( ||u||²_L2 + ||∇u||²_L2 )
 *          = sqrt( ∫_Ω |u(x)|² dx + ∫_Ω |∇u(x)|² dx )
 * ```
 *
 * This norm measures both the function values and their derivatives,
 * making it suitable for error analysis in elliptic PDEs.
 *
 * **Note:** This function returns the H1 seminorm (gradient norm only).
 * To get the full H1 norm, combine with the L2 norm:
 * `sqrt(normL2(gf)**2 + normH1(gf)**2)`
 *
 * @param gf - GridFunction to compute norm of
 * @returns The H1 seminorm value (gradient norm)
 * @throws Error if the grid function has been destroyed
 * @throws Error if norm computation fails
 *
 * @category GridFunction
 *
 * @example Compute full H1 norm
 * ```typescript
 * const l2 = normL2(solution);
 * const h1Semi = normH1(solution);
 * const h1Full = Math.sqrt(l2 * l2 + h1Semi * h1Semi);
 * console.log(`Full H1 norm: ${h1Full}`);
 * ```
 *
 * @example Check gradient magnitude
 * ```typescript
 * const gradNorm = normH1(temperature);
 * console.log(`Temperature gradient norm: ${gradNorm}`);
 * ```
 *
 * @see {@link normL2} for L2 norm (function values only)
 * @see {@link normLinf} for maximum norm
 */
export function normH1(gf: GridFunction): number {
  const module = gf.getModule();
  const ptr = gf.getPointer();

  const norm = module._mfem_gridfunc_norm_h1_semi(ptr);

  if (norm < 0) {
    throw new Error('Failed to compute H1 seminorm');
  }

  return norm;
}
