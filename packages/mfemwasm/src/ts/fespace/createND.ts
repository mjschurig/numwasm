/**
 * Nedelec finite element space creation.
 *
 * @module fespace
 */

import type { Mesh } from '../mesh/Mesh.js';
import { FiniteElementSpace } from './FiniteElementSpace.js';
import { createHcurl } from './createHcurl.js';

/**
 * Creates a Nedelec (ND) finite element space.
 *
 * This is an alias for {@link createHcurl}. Nedelec elements are edge elements
 * that conform to H(curl), named after Jean-Claude Nédélec who introduced them.
 *
 * Common applications:
 * - Electric field in electromagnetics
 * - Magnetic vector potential
 * - Any vector field where curl is the primary differential operator
 *
 * @param mesh - The mesh on which to define the space (must be 2D or 3D)
 * @param order - Polynomial order (1 = lowest order, etc.)
 * @returns A promise that resolves to the new FiniteElementSpace
 * @throws Error if order is less than 1
 * @throws Error if space creation fails
 *
 * @category FiniteElementSpace
 *
 * @example Nedelec space for curl-curl problem
 * ```typescript
 * const mesh = await makeCartesian3D(5, 5, 5);
 * const fespace = await createND(mesh, 2); // Second-order Nedelec
 * ```
 *
 * @see {@link createHcurl} - The underlying H(curl) space creation function
 */
export async function createND(
  mesh: Mesh,
  order: number
): Promise<FiniteElementSpace> {
  return createHcurl(mesh, order);
}
