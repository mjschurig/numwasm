/**
 * H(div) finite element space creation.
 *
 * @module fespace
 */

import { loadMFEMModule } from '../loader.js';
import type { Mesh } from '../mesh/Mesh.js';
import { FiniteElementSpace } from './FiniteElementSpace.js';

/**
 * Creates an H(div) (Raviart-Thomas) finite element space.
 *
 * H(div) spaces are designed for vector fields where the divergence is
 * well-defined and square-integrable. They are essential for:
 * - Fluid flow (velocity fields)
 * - Mixed finite element methods
 * - Darcy flow and porous media
 * - Magnetic flux density in electromagnetics
 *
 * The normal component is continuous across element faces while the
 * tangential component may be discontinuous.
 *
 * @param mesh - The mesh on which to define the space (must be 2D or 3D)
 * @param order - Polynomial order (0 = lowest order RT, etc.)
 * @returns A promise that resolves to the new FiniteElementSpace
 * @throws Error if order is negative
 * @throws Error if space creation fails
 *
 * @category FiniteElementSpace
 *
 * @example H(div) space for mixed Poisson
 * ```typescript
 * const mesh = await makeCartesian2D(10, 10);
 * const fespace = await createHdiv(mesh, 0); // Lowest order RT
 * // DOFs are associated with faces/edges
 * ```
 *
 * @see {@link createRT} - Alias using Raviart-Thomas (RT) naming convention
 */
export async function createHdiv(
  mesh: Mesh,
  order: number
): Promise<FiniteElementSpace> {
  if (order < 0) {
    throw new Error('H(div) order must be at least 0');
  }

  const module = await loadMFEMModule();
  const ptr = module._mfem_fespace_create_rt(mesh.getPointer(), order);

  if (ptr === 0) {
    throw new Error('Failed to create H(div) finite element space');
  }

  return new FiniteElementSpace(module, ptr);
}
