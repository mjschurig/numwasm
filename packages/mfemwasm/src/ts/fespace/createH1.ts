/**
 * H1 finite element space creation.
 *
 * @module fespace
 */

import { loadMFEMModule } from '../loader.js';
import type { Mesh } from '../mesh/Mesh.js';
import { FiniteElementSpace } from './FiniteElementSpace.js';

/**
 * Creates an H1 (continuous Lagrange) finite element space.
 *
 * H1 spaces use continuous piecewise polynomial functions that are globally
 * conforming in H1 (functions with square-integrable first derivatives).
 * They are the most common choice for:
 * - Scalar problems (heat transfer, Poisson equation)
 * - Vector problems (elasticity, Stokes flow)
 *
 * @param mesh - The mesh on which to define the space
 * @param order - Polynomial order (1 = linear, 2 = quadratic, 3 = cubic, etc.)
 * @param vdim - Vector dimension (default 1 for scalar fields)
 * @returns A promise that resolves to the new FiniteElementSpace
 * @throws Error if order is less than 1
 * @throws Error if space creation fails
 *
 * @category FiniteElementSpace
 *
 * @example Scalar H1 space for Poisson problem
 * ```typescript
 * const mesh = await makeCartesian2D(10, 10);
 * const fespace = await createH1(mesh, 2); // Quadratic
 * console.log(`DOFs: ${fespace.ndofs}`);
 * ```
 *
 * @example Vector H1 space for 2D elasticity
 * ```typescript
 * const mesh = await makeCartesian2D(10, 10);
 * const fespace = await createH1(mesh, 1, 2); // Linear, 2 components
 * console.log(`Total DOFs: ${fespace.ndofs}`);
 * ```
 */
export async function createH1(
  mesh: Mesh,
  order: number,
  vdim: number = 1
): Promise<FiniteElementSpace> {
  if (order < 1) {
    throw new Error('H1 order must be at least 1');
  }

  const module = await loadMFEMModule();
  const ptr = module._mfem_fespace_create_h1(mesh.getPointer(), order, vdim);

  if (ptr === 0) {
    throw new Error('Failed to create H1 finite element space');
  }

  return new FiniteElementSpace(module, ptr);
}
