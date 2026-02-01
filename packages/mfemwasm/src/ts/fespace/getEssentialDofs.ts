/**
 * Essential (Dirichlet) DOF extraction.
 *
 * @module fespace
 */

import { FiniteElementSpace } from './FiniteElementSpace.js';

/**
 * Gets the essential (Dirichlet) DOF indices for specified boundary attributes.
 *
 * Returns the true (global) DOF indices that should be constrained for
 * Dirichlet boundary conditions on the specified boundary attributes.
 * Multiple boundary attributes can be specified to get DOFs from several
 * boundary regions at once.
 *
 * @param fespace - The finite element space
 * @param bdrAttrs - Array of boundary attributes (1-based, as defined in the mesh)
 * @returns Array of essential DOF indices
 * @throws Error if the space has been destroyed
 * @throws Error if memory allocation fails
 *
 * @category FiniteElementSpace
 *
 * @example Get all Dirichlet DOFs
 * ```typescript
 * const fespace = await createH1(mesh, 2);
 * // DOFs on boundaries 1 and 3
 * const essDofs = getEssentialDofs(fespace, [1, 3]);
 * console.log(`${essDofs.length} constrained DOFs`);
 * ```
 *
 * @example Eliminate DOFs from linear system
 * ```typescript
 * const essDofs = getEssentialDofs(fespace, [1, 2, 3, 4]);
 * // Use essDofs to eliminate rows/columns from stiffness matrix
 * ```
 */
export function getEssentialDofs(fespace: FiniteElementSpace, bdrAttrs: number[]): number[] {
  const module = fespace.getModule();
  const ptr = fespace.getPointer();

  // Allocate buffer for boundary attributes
  const attrsPtr = module._malloc(bdrAttrs.length * 4);
  if (attrsPtr === 0) {
    throw new Error('Failed to allocate memory for attributes array');
  }

  // Allocate buffer for DOF indices
  const maxDofs = fespace.ndofs;
  const dofsPtr = module._malloc(maxDofs * 4);
  if (dofsPtr === 0) {
    module._free(attrsPtr);
    throw new Error('Failed to allocate memory for DOF array');
  }

  try {
    // Copy boundary attributes to WASM memory
    for (let i = 0; i < bdrAttrs.length; i++) {
      module.HEAP32[(attrsPtr >> 2) + i] = bdrAttrs[i];
    }

    const numDofs = module._mfem_fespace_get_essential_dofs(
      ptr,
      attrsPtr,
      bdrAttrs.length,
      dofsPtr,
      maxDofs
    );

    if (numDofs < 0) {
      throw new Error('Failed to get essential DOFs');
    }

    const dofs: number[] = [];
    const count = Math.min(numDofs, maxDofs);
    for (let i = 0; i < count; i++) {
      dofs.push(module.HEAP32[(dofsPtr >> 2) + i]);
    }
    return dofs;
  } finally {
    module._free(attrsPtr);
    module._free(dofsPtr);
  }
}
