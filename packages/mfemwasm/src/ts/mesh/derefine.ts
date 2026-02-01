/**
 * Mesh derefinement (coarsening) operations.
 *
 * @module mesh
 */

import { Mesh } from './Mesh.js';

/**
 * Derefines (coarsens) specified elements in a non-conforming mesh.
 *
 * Derefinement is the reverse of refinement - it merges previously refined
 * elements back into their parent elements. This is useful for adaptive
 * mesh refinement (AMR) where elements may need to be coarsened after
 * solution features have moved or diminished.
 *
 * **Requirements:**
 * - The mesh must have NCMesh (non-conforming mesh) support enabled
 * - The elements must have been previously refined
 * - The derefinement indices refer to rows in the NCMesh derefinement table
 *
 * To enable NCMesh support before refinement, use {@link enableNCMesh}.
 *
 * @param mesh - Mesh instance with NCMesh support
 * @param derefinements - Array of derefinement table row indices
 * @throws Error if the mesh does not have NCMesh support
 * @throws Error if derefinement fails
 *
 * @category Mesh
 *
 * @example Adaptive mesh refinement/derefinement cycle
 * ```typescript
 * // Enable NCMesh before any refinement
 * enableNCMesh(mesh);
 *
 * // Refine some elements
 * refineLocal(mesh, [0, 1, 2, 3]);
 *
 * // Later, coarsen some elements (using derefinement table indices)
 * const numDerefs = getDerefTableSize(mesh);
 * if (numDerefs > 0) {
 *   derefine(mesh, [0, 1]); // Derefine first two possible coarsenings
 * }
 * ```
 *
 * @see {@link enableNCMesh} - Enable non-conforming mesh support
 * @see {@link isNCMesh} - Check if mesh has NCMesh support
 * @see {@link getDerefTableSize} - Get number of possible derefinements
 * @see {@link refineLocal} - Local mesh refinement
 */
export function derefine(mesh: Mesh, derefinements: number[]): void {
  const module = mesh.getModule();
  const ptr = mesh.getPointer();

  // Verify this is an NCMesh
  if (!module._mfem_mesh_is_ncmesh(ptr)) {
    throw new Error(
      'Cannot derefine: mesh does not have NCMesh support. ' +
      'Enable NCMesh with enableNCMesh() before refinement.'
    );
  }

  if (derefinements.length === 0) {
    return; // Nothing to do
  }

  // Allocate buffer for derefinement indices
  const derefsPtr = module._malloc(derefinements.length * 4);
  if (derefsPtr === 0) {
    throw new Error('Failed to allocate memory for derefinement array');
  }

  try {
    // Copy derefinement indices to WASM memory
    for (let i = 0; i < derefinements.length; i++) {
      module.HEAP32[(derefsPtr >> 2) + i] = derefinements[i];
    }

    const result = module._mfem_mesh_derefine(ptr, derefsPtr, derefinements.length);

    if (result < 0) {
      throw new Error('Derefinement failed');
    }
  } finally {
    module._free(derefsPtr);
  }
}

/**
 * Enables non-conforming mesh (NCMesh) support.
 *
 * NCMesh support is required for local refinement with hanging nodes
 * and for derefinement (coarsening). This function must be called
 * before any local refinement operations if derefinement is desired later.
 *
 * @param mesh - Mesh instance
 * @throws Error if enabling NCMesh fails
 *
 * @category Mesh
 *
 * @example Enable NCMesh for adaptive refinement
 * ```typescript
 * const mesh = makeCartesian2D(10, 10);
 * enableNCMesh(mesh);
 *
 * // Now local refinement will create hanging nodes
 * refineLocal(mesh, [0, 1, 2]);
 *
 * // And derefinement is possible
 * derefine(mesh, [0]);
 * ```
 */
export function enableNCMesh(mesh: Mesh): void {
  const module = mesh.getModule();
  const ptr = mesh.getPointer();

  const result = module._mfem_mesh_enable_ncmesh(ptr);

  if (result < 0) {
    throw new Error('Failed to enable NCMesh support');
  }
}

/**
 * Checks if a mesh has NCMesh (non-conforming mesh) support.
 *
 * @param mesh - Mesh instance
 * @returns True if the mesh has NCMesh support enabled
 *
 * @category Mesh
 *
 * @example Check before derefinement
 * ```typescript
 * if (isNCMesh(mesh)) {
 *   const numDerefs = getDerefTableSize(mesh);
 *   console.log(`${numDerefs} possible derefinements`);
 * }
 * ```
 */
export function isNCMesh(mesh: Mesh): boolean {
  const module = mesh.getModule();
  return module._mfem_mesh_is_ncmesh(mesh.getPointer()) === 1;
}

/**
 * Gets the size of the derefinement table.
 *
 * The derefinement table contains information about which element groups
 * can be coarsened. Each row in the table represents one possible
 * derefinement operation.
 *
 * @param mesh - Mesh instance with NCMesh support
 * @returns Number of possible derefinements, or 0 if not an NCMesh
 * @throws Error if query fails
 *
 * @category Mesh
 *
 * @example Get available derefinements
 * ```typescript
 * enableNCMesh(mesh);
 * refineLocal(mesh, [0, 1, 2, 3]);
 *
 * const numDerefs = getDerefTableSize(mesh);
 * console.log(`Can derefine ${numDerefs} element groups`);
 *
 * // Derefine all possible elements
 * if (numDerefs > 0) {
 *   const allDerefs = Array.from({ length: numDerefs }, (_, i) => i);
 *   derefine(mesh, allDerefs);
 * }
 * ```
 */
export function getDerefTableSize(mesh: Mesh): number {
  const module = mesh.getModule();
  const ptr = mesh.getPointer();

  if (!module._mfem_mesh_is_ncmesh(ptr)) {
    return 0;
  }

  const size = module._mfem_mesh_get_deref_table_size(ptr);

  if (size < 0) {
    throw new Error('Failed to get derefinement table size');
  }

  return size;
}
