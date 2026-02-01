/**
 * Compute L2 norm of grid function.
 *
 * @module gridfunction
 */

import { GridFunction } from './GridFunction.js';

/**
 * Computes the L2 norm of a grid function.
 *
 * The L2 norm is defined as:
 * ```
 * ||u||_L2 = sqrt( ∫_Ω |u(x)|² dx )
 * ```
 *
 * This is computed using numerical quadrature over all elements in the mesh.
 *
 * @param gf - GridFunction to compute norm of
 * @returns The L2 norm value
 * @throws Error if the grid function has been destroyed
 * @throws Error if norm computation fails
 *
 * @category GridFunction
 *
 * @example Check solution magnitude
 * ```typescript
 * const norm = normL2(solution);
 * console.log(`Solution L2 norm: ${norm}`);
 * ```
 *
 * @example Normalize a function
 * ```typescript
 * const norm = normL2(gf);
 * const data = gf.getData();
 * for (let i = 0; i < data.length; i++) {
 *   data[i] /= norm;
 * }
 * gf.setData(data);
 * ```
 *
 * @see {@link normLinf} for maximum norm
 * @see {@link normH1} for H1 norm including gradients
 */
export function normL2(gf: GridFunction): number {
  const module = gf.getModule();
  const ptr = gf.getPointer();

  const norm = module._mfem_gridfunc_norm_l2(ptr);

  if (norm < 0) {
    throw new Error('Failed to compute L2 norm');
  }

  return norm;
}
