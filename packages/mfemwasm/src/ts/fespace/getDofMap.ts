/**
 * Finite element space DOF mapping.
 *
 * @module fespace
 */

import { FiniteElementSpace } from './FiniteElementSpace.js';

/**
 * Gets the DOF indices for a specific element.
 *
 * Returns the global DOF indices associated with the specified element.
 * This is useful for assembling local element contributions into global
 * matrices and vectors.
 *
 * @param fespace - The finite element space
 * @param elemIndex - The element index (0-based)
 * @returns Array of global DOF indices for this element
 * @throws Error if the space has been destroyed
 * @throws Error if the element index is out of bounds
 *
 * @category FiniteElementSpace
 *
 * @example Get DOFs for element assembly
 * ```typescript
 * const fespace = await createH1(mesh, 2);
 * const dofs = getDofMap(fespace, 0); // DOFs for first element
 * console.log(`Element 0 has ${dofs.length} DOFs`);
 * ```
 */
export function getDofMap(fespace: FiniteElementSpace, elemIndex: number): Int32Array {
  const module = fespace.getModule();
  const ptr = fespace.getPointer();

  // Allocate enough space for DOFs (worst case: high-order hex has many DOFs)
  const maxDofs = 512;
  const dofsPtr = module._malloc(maxDofs * 4);

  if (dofsPtr === 0) {
    throw new Error('Failed to allocate memory for DOF array');
  }

  try {
    const numDofs = module._mfem_fespace_get_element_dofs(ptr, elemIndex, dofsPtr);

    if (numDofs < 0) {
      throw new Error(`Failed to get DOFs for element ${elemIndex}`);
    }

    const dofs = new Int32Array(numDofs);
    for (let i = 0; i < numDofs; i++) {
      dofs[i] = module.HEAP32[(dofsPtr >> 2) + i];
    }
    return dofs;
  } finally {
    module._free(dofsPtr);
  }
}
