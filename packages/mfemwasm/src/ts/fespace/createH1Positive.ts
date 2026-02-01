/**
 * H1 positive (Bernstein) finite element space creation.
 *
 * @module fespace
 */

import { loadMFEMModule } from '../loader.js';
import type { Mesh } from '../mesh/Mesh.js';
import { FiniteElementSpace } from './FiniteElementSpace.js';

/**
 * Creates an H1 finite element space with positive (Bernstein) basis.
 *
 * Bernstein basis functions are non-negative and form a partition of unity,
 * making them ideal for applications requiring positivity preservation:
 * - Concentration fields that must remain positive
 * - Probability distributions
 * - NURBS-like representations
 *
 * The space has the same approximation properties as standard H1 but with
 * different basis function properties.
 *
 * @param mesh - The mesh on which to define the space
 * @param order - Polynomial order (1 = linear, 2 = quadratic, etc.)
 * @param vdim - Vector dimension (default 1 for scalar fields)
 * @returns A promise that resolves to the new FiniteElementSpace
 * @throws Error if order is less than 1
 * @throws Error if space creation fails
 *
 * @category FiniteElementSpace
 *
 * @example Positive basis for concentration field
 * ```typescript
 * const mesh = await makeCartesian2D(10, 10);
 * const fespace = await createH1Positive(mesh, 2);
 * // Coefficients directly control positivity
 * ```
 */
export async function createH1Positive(
  mesh: Mesh,
  order: number,
  vdim: number = 1
): Promise<FiniteElementSpace> {
  if (order < 1) {
    throw new Error('H1 order must be at least 1');
  }

  const module = await loadMFEMModule();
  const ptr = module._mfem_fespace_create_h1_positive(mesh.getPointer(), order, vdim);

  if (ptr === 0) {
    throw new Error('Failed to create H1 positive finite element space');
  }

  return new FiniteElementSpace(module, ptr);
}
