/**
 * Local mesh refinement.
 *
 * @module mesh
 */

import { getMFEMModule } from '../loader.js';
import { Mesh } from './Mesh.js';

/**
 * Locally refines specified elements in the mesh.
 *
 * Performs adaptive refinement by subdividing only the specified elements,
 * leaving other elements unchanged. This is useful for adaptive mesh refinement
 * (AMR) strategies where refinement is focused on regions with high error.
 *
 * The mesh maintains conformity by adding necessary transition elements between
 * refined and unrefined regions.
 *
 * @param mesh - The mesh to refine (modified in place)
 * @param elements - Array of element indices to refine (0-based)
 * @throws Error if the MFEM module is not loaded
 * @throws Error if memory allocation fails
 * @throws Error if refinement fails
 *
 * @category Mesh
 *
 * @example Refine specific elements
 * ```typescript
 * const mesh = await makeCartesian2D(4, 4);
 * console.log(mesh.numElements); // 16
 *
 * // Refine elements 0, 1, and 5
 * refineLocal(mesh, [0, 1, 5]);
 * console.log(mesh.numElements); // More than 16 (refined elements split)
 * mesh.destroy();
 * ```
 *
 * @example Adaptive refinement based on error indicator
 * ```typescript
 * const mesh = await makeCartesian2D(8, 8);
 *
 * // Compute error indicator for each element (placeholder)
 * const errors = getElementVolumes(mesh).map((_, i) => Math.random());
 * const threshold = 0.8;
 *
 * // Refine elements with high error
 * const elementsToRefine = errors
 *   .map((err, i) => ({ i, err }))
 *   .filter(e => e.err > threshold)
 *   .map(e => e.i);
 *
 * refineLocal(mesh, elementsToRefine);
 * mesh.destroy();
 * ```
 *
 * @example Refine elements near a point
 * ```typescript
 * const mesh = await makeCartesian2D(10, 10);
 * const centroids = getElementCentroids(mesh);
 * const target = [0.5, 0.5];
 *
 * // Find elements within distance 0.2 of target
 * const nearTarget = centroids
 *   .map((c, i) => ({
 *     index: i,
 *     dist: Math.hypot(c[0] - target[0], c[1] - target[1])
 *   }))
 *   .filter(e => e.dist < 0.2)
 *   .map(e => e.index);
 *
 * refineLocal(mesh, nearTarget);
 * mesh.destroy();
 * ```
 *
 * @see {@link refineUniform} - Refine all elements uniformly
 * @see {@link refineNonconforming} - Non-conforming local refinement
 * @see {@link getElementCentroids} - Find element locations
 */
export function refineLocal(mesh: Mesh, elements: number[]): void {
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
      throw new Error('Failed to refine mesh locally');
    }
  } finally {
    module._free(elementsPtr);
  }
}
