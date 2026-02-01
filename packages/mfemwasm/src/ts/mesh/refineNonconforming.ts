/**
 * Non-conforming mesh refinement.
 *
 * @module mesh
 */

import { getMFEMModule } from '../loader.js';
import { Mesh } from './Mesh.js';

/**
 * Performs non-conforming refinement of specified elements.
 *
 * Non-conforming refinement allows hanging nodes (vertices that lie on
 * edges/faces of neighboring elements) which enables more flexible local
 * refinement without requiring transition elements.
 *
 * **Note:** The current implementation uses the same mechanism as
 * {@link refineLocal}. True non-conforming refinement with hanging node
 * support requires full NCMesh integration which is not yet available.
 *
 * @param mesh - The mesh to refine (modified in place)
 * @param elements - Array of element indices to refine (0-based)
 * @throws Error if the MFEM module is not loaded
 * @throws Error if memory allocation fails
 * @throws Error if refinement fails
 *
 * @category Mesh
 *
 * @example Basic non-conforming refinement
 * ```typescript
 * const mesh = await makeCartesian2D(4, 4);
 *
 * // Refine just the corner elements
 * refineNonconforming(mesh, [0, 3, 12, 15]);
 * mesh.destroy();
 * ```
 *
 * @example AMR with non-conforming refinement
 * ```typescript
 * const mesh = await makeCartesian2D(8, 8);
 *
 * // Iterative refinement loop
 * for (let iter = 0; iter < 3; iter++) {
 *   const centroids = getElementCentroids(mesh);
 *
 *   // Mark elements near singularity at origin
 *   const toRefine = centroids
 *     .map((c, i) => ({ i, r: Math.hypot(c[0], c[1]) }))
 *     .filter(e => e.r < 0.1 * (iter + 1))
 *     .map(e => e.i);
 *
 *   if (toRefine.length > 0) {
 *     refineNonconforming(mesh, toRefine);
 *   }
 * }
 * mesh.destroy();
 * ```
 *
 * @see {@link refineLocal} - Conforming local refinement
 * @see {@link refineUniform} - Uniform global refinement
 * @see {@link derefine} - Mesh coarsening (not yet implemented)
 */
export function refineNonconforming(mesh: Mesh, elements: number[]): void {
  if (elements.length === 0) {
    return;
  }

  const module = getMFEMModule();
  if (!module) {
    throw new Error('MFEM module not loaded');
  }

  const elementsPtr = module._malloc(elements.length * 4);
  if (elementsPtr === 0) {
    throw new Error('Failed to allocate memory for elements array');
  }

  try {
    // Copy elements to WASM memory
    for (let i = 0; i < elements.length; i++) {
      module.HEAP32[(elementsPtr >> 2) + i] = elements[i];
    }

    const result = module._mfem_mesh_refine_local(
      mesh.getPointer(),
      elementsPtr,
      elements.length
    );

    if (result !== 0) {
      throw new Error('Failed to refine mesh');
    }
  } finally {
    module._free(elementsPtr);
  }
}
