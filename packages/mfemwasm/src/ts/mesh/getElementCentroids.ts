/**
 * Element centroid computation.
 *
 * @module mesh
 */

import { getMFEMModule } from '../loader.js';
import { Mesh } from './Mesh.js';

/**
 * Computes the centroid (geometric center) of each element in the mesh.
 *
 * The centroid is the average of the element's vertex coordinates and
 * represents the center of mass assuming uniform density.
 *
 * @param mesh - The mesh to compute element centroids for
 * @returns Array of coordinate arrays, where centroids[i] contains the
 *          coordinates [x, y] (2D) or [x, y, z] (3D) of element i's centroid
 * @throws Error if the MFEM module is not loaded
 * @throws Error if memory allocation fails
 *
 * @category Mesh
 *
 * @example Get centroids of all elements
 * ```typescript
 * const mesh = await makeCartesian2D(2, 2);
 * const centroids = getElementCentroids(mesh);
 * console.log(centroids.length); // 4
 * console.log(centroids[0]); // [0.25, 0.25] (center of first element)
 * mesh.destroy();
 * ```
 *
 * @example Compute average location of elements
 * ```typescript
 * const mesh = await makeCartesian3D(4, 4, 4);
 * const centroids = getElementCentroids(mesh);
 * const avgX = centroids.reduce((sum, c) => sum + c[0], 0) / centroids.length;
 * const avgY = centroids.reduce((sum, c) => sum + c[1], 0) / centroids.length;
 * const avgZ = centroids.reduce((sum, c) => sum + c[2], 0) / centroids.length;
 * console.log(`Center of mass: (${avgX}, ${avgY}, ${avgZ})`);
 * mesh.destroy();
 * ```
 *
 * @example Find elements in a region
 * ```typescript
 * const mesh = await makeCartesian2D(10, 10);
 * const centroids = getElementCentroids(mesh);
 *
 * // Find elements whose centers are in the circle r < 0.3
 * const elementsInCircle = centroids
 *   .map((c, i) => ({ index: i, x: c[0], y: c[1] }))
 *   .filter(e => {
 *     const dx = e.x - 0.5, dy = e.y - 0.5;
 *     return Math.sqrt(dx*dx + dy*dy) < 0.3;
 *   })
 *   .map(e => e.index);
 *
 * console.log(`${elementsInCircle.length} elements in circle`);
 * mesh.destroy();
 * ```
 *
 * @example Select elements for local refinement
 * ```typescript
 * const mesh = await makeCartesian2D(8, 8);
 * const centroids = getElementCentroids(mesh);
 *
 * // Refine elements near the origin
 * const nearOrigin = centroids
 *   .map((c, i) => ({ i, dist: Math.hypot(c[0], c[1]) }))
 *   .filter(e => e.dist < 0.2)
 *   .map(e => e.i);
 *
 * refineLocal(mesh, nearOrigin);
 * mesh.destroy();
 * ```
 *
 * @see {@link getElementVolumes} - Get element volumes/areas
 */
export function getElementCentroids(mesh: Mesh): number[][] {
  const module = getMFEMModule();
  if (!module) {
    throw new Error('MFEM module not loaded');
  }

  const ne = mesh.numElements;
  const sdim = mesh.spaceDimension;
  const centroidsPtr = module._malloc(ne * sdim * 8);

  if (centroidsPtr === 0) {
    throw new Error('Failed to allocate memory for centroids');
  }

  try {
    const result = module._mfem_mesh_get_element_centroids(mesh.getPointer(), centroidsPtr);
    if (result !== 0) {
      throw new Error('Failed to get element centroids');
    }

    const centroids: number[][] = [];
    for (let i = 0; i < ne; i++) {
      const centroid: number[] = [];
      for (let d = 0; d < sdim; d++) {
        centroid.push(module.HEAPF64[(centroidsPtr >> 3) + i * sdim + d]);
      }
      centroids.push(centroid);
    }
    return centroids;
  } finally {
    module._free(centroidsPtr);
  }
}
