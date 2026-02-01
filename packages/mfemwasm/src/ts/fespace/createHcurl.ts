/**
 * H(curl) finite element space creation.
 *
 * @module fespace
 */

import { loadMFEMModule } from '../loader.js';
import type { Mesh } from '../mesh/Mesh.js';
import { FiniteElementSpace } from './FiniteElementSpace.js';

/**
 * Creates an H(curl) (Nedelec) finite element space.
 *
 * H(curl) spaces are designed for vector fields where the curl is well-defined
 * and square-integrable. They are essential for:
 * - Electromagnetics (electric field, magnetic vector potential)
 * - Maxwell's equations
 * - Edge element discretizations
 *
 * The tangential component is continuous across element faces while the
 * normal component may be discontinuous.
 *
 * @param mesh - The mesh on which to define the space (must be 2D or 3D)
 * @param order - Polynomial order (1 = lowest order Nedelec, etc.)
 * @returns A promise that resolves to the new FiniteElementSpace
 * @throws Error if order is less than 1
 * @throws Error if space creation fails
 *
 * @category FiniteElementSpace
 *
 * @example H(curl) space for electromagnetics
 * ```typescript
 * const mesh = await makeCartesian3D(5, 5, 5);
 * const fespace = await createHcurl(mesh, 1);
 * // DOFs are associated with edges
 * ```
 *
 * @see {@link createND} - Alias using Nedelec (ND) naming convention
 */
export async function createHcurl(
  mesh: Mesh,
  order: number
): Promise<FiniteElementSpace> {
  if (order < 1) {
    throw new Error('H(curl) order must be at least 1');
  }

  const module = await loadMFEMModule();
  const ptr = module._mfem_fespace_create_nd(mesh.getPointer(), order);

  if (ptr === 0) {
    throw new Error('Failed to create H(curl) finite element space');
  }

  return new FiniteElementSpace(module, ptr);
}
