/**
 * L2 (discontinuous) finite element space creation.
 *
 * @module fespace
 */

import { loadMFEMModule } from '../loader.js';
import type { Mesh } from '../mesh/Mesh.js';
import { FiniteElementSpace } from './FiniteElementSpace.js';

/**
 * Creates an L2 (discontinuous) finite element space.
 *
 * L2 spaces use piecewise polynomial functions with no continuity constraints
 * across element boundaries. This allows for:
 * - Discontinuous Galerkin (DG) methods
 * - Higher accuracy for advection-dominated problems
 * - Local conservation properties
 * - Representation of discontinuous fields
 *
 * L2 spaces have more DOFs than H1 spaces for the same mesh and order,
 * as each element has its own independent polynomial representation.
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
 * @example Piecewise constant space for DG
 * ```typescript
 * const mesh = await makeCartesian2D(10, 10);
 * const fespace = await createL2(mesh, 0); // P0 (constant per element)
 * console.log(`DOFs: ${fespace.ndofs}`); // Equals number of elements
 * ```
 *
 * @example High-order DG for advection
 * ```typescript
 * const mesh = await makeCartesian2D(10, 10);
 * const fespace = await createL2(mesh, 3); // Cubic polynomials per element
 * ```
 */
export async function createL2(
  mesh: Mesh,
  order: number,
  vdim: number = 1
): Promise<FiniteElementSpace> {
  if (order < 0) {
    throw new Error('L2 order must be at least 0');
  }

  const module = await loadMFEMModule();
  const ptr = module._mfem_fespace_create_l2(mesh.getPointer(), order, vdim);

  if (ptr === 0) {
    throw new Error('Failed to create L2 finite element space');
  }

  return new FiniteElementSpace(module, ptr);
}
