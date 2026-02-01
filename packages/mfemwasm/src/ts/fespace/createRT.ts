/**
 * Raviart-Thomas finite element space creation.
 *
 * @module fespace
 */

import type { Mesh } from '../mesh/Mesh.js';
import { FiniteElementSpace } from './FiniteElementSpace.js';
import { createHdiv } from './createHdiv.js';

/**
 * Creates a Raviart-Thomas (RT) finite element space.
 *
 * This is an alias for {@link createHdiv}. Raviart-Thomas elements are face
 * elements that conform to H(div), named after Pierre-Arnaud Raviart and
 * Jean-Marie Thomas who introduced them.
 *
 * Common applications:
 * - Velocity in mixed formulations
 * - Flux variables in conservation laws
 * - Magnetic flux density
 *
 * @param mesh - The mesh on which to define the space (must be 2D or 3D)
 * @param order - Polynomial order (0 = lowest order, etc.)
 * @returns A promise that resolves to the new FiniteElementSpace
 * @throws Error if order is negative
 * @throws Error if space creation fails
 *
 * @category FiniteElementSpace
 *
 * @example RT space for Darcy flow
 * ```typescript
 * const mesh = await makeCartesian2D(10, 10);
 * const fespace = await createRT(mesh, 1); // First-order RT
 * ```
 *
 * @see {@link createHdiv} - The underlying H(div) space creation function
 */
export async function createRT(
  mesh: Mesh,
  order: number
): Promise<FiniteElementSpace> {
  return createHdiv(mesh, order);
}
