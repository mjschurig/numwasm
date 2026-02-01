/**
 * Discontinuous Galerkin finite element space creation.
 *
 * @module fespace
 */

import type { Mesh } from '../mesh/Mesh.js';
import { FiniteElementSpace } from './FiniteElementSpace.js';
import { createL2 } from './createL2.js';

/**
 * Creates a Discontinuous Galerkin (DG) finite element space.
 *
 * This is an alias for {@link createL2}. DG spaces are L2 spaces used in
 * Discontinuous Galerkin methods, where numerical fluxes connect elements.
 *
 * DG methods are popular for:
 * - Hyperbolic conservation laws
 * - Advection-dominated transport
 * - Problems with sharp gradients or shocks
 * - hp-adaptivity
 *
 * @param mesh - The mesh on which to define the space
 * @param order - Polynomial order (0 = piecewise constant, 1 = linear, etc.)
 * @param vdim - Vector dimension (default 1 for scalar fields)
 * @returns A promise that resolves to the new FiniteElementSpace
 * @throws Error if order is negative
 * @throws Error if space creation fails
 *
 * @category FiniteElementSpace
 *
 * @example DG space for transport equation
 * ```typescript
 * const mesh = await makeCartesian2D(20, 20);
 * const fespace = await createDG(mesh, 2); // P2 DG
 * ```
 *
 * @see {@link createL2} - The underlying L2 space creation function
 */
export async function createDG(
  mesh: Mesh,
  order: number,
  vdim: number = 1
): Promise<FiniteElementSpace> {
  return createL2(mesh, order, vdim);
}
