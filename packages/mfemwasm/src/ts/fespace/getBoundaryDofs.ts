/**
 * Boundary DOF extraction.
 *
 * @module fespace
 */

import { FiniteElementSpace } from './FiniteElementSpace.js';

/**
 * Gets the DOF indices on a boundary with a specific attribute.
 *
 * Returns the true (global) DOF indices that lie on boundary elements
 * with the specified attribute. This is useful for applying Dirichlet
 * boundary conditions on specific parts of the domain boundary.
 *
 * @param fespace - The finite element space
 * @param bdrAttr - The boundary attribute (1-based, as defined in the mesh)
 * @returns Array of boundary DOF indices
 * @throws Error if the space has been destroyed
 * @throws Error if memory allocation fails
 *
 * @category FiniteElementSpace
 *
 * @example Get DOFs on inlet boundary (attribute 1)
 * ```typescript
 * const fespace = await createH1(mesh, 2);
 * const inletDofs = getBoundaryDofs(fespace, 1);
 * console.log(`${inletDofs.length} DOFs on inlet`);
 * ```
 *
 * @example Apply Dirichlet BC on specific boundary
 * ```typescript
 * const wallDofs = getBoundaryDofs(fespace, 2);
 * for (const dof of wallDofs) {
 *   // Set boundary value at this DOF
 *   solution.setData(dof, 0.0);
 * }
 * ```
 */
export function getBoundaryDofs(fespace: FiniteElementSpace, bdrAttr: number): number[] {
  const module = fespace.getModule();
  const ptr = fespace.getPointer();

  // Allocate buffer for DOF indices
  const maxDofs = fespace.ndofs;
  const dofsPtr = module._malloc(maxDofs * 4);

  if (dofsPtr === 0) {
    throw new Error('Failed to allocate memory for DOF array');
  }

  try {
    const numDofs = module._mfem_fespace_get_boundary_dofs(ptr, bdrAttr, dofsPtr, maxDofs);

    if (numDofs < 0) {
      throw new Error(`Failed to get boundary DOFs for attribute ${bdrAttr}`);
    }

    const dofs: number[] = [];
    const count = Math.min(numDofs, maxDofs);
    for (let i = 0; i < count; i++) {
      dofs.push(module.HEAP32[(dofsPtr >> 2) + i]);
    }
    return dofs;
  } finally {
    module._free(dofsPtr);
  }
}
