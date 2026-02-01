/**
 * Circular mesh generation.
 *
 * @module mesh
 */

import { loadMFEMModule } from '../loader.js';
import { Mesh } from './Mesh.js';

/**
 * Creates a 2D circular disk mesh.
 *
 * Generates a triangulated mesh approximating a circular disk by creating a
 * square mesh and applying an elliptic conformal mapping. This mapping
 * preserves mesh quality while smoothly transforming the square to a circle.
 *
 * The mesh is centered at the origin with the specified radius.
 *
 * @param npts - Controls mesh density. Higher values create finer meshes (minimum 3)
 * @param radius - Radius of the circular disk (default 1.0)
 * @returns A promise that resolves to a 2D triangular Mesh
 * @throws Error if npts is less than 3
 * @throws Error if mesh creation fails
 *
 * @category Mesh
 *
 * @example Basic unit circle mesh
 * ```typescript
 * const mesh = await makeCircle(20);
 * console.log(mesh.dimension); // 2
 * // Mesh is centered at origin with radius 1
 * mesh.destroy();
 * ```
 *
 * @example Circular domain with custom radius
 * ```typescript
 * // Create a disk of radius 0.5
 * const mesh = await makeCircle(30, 0.5);
 * mesh.destroy();
 * ```
 *
 * @example Solving Laplace equation on disk
 * ```typescript
 * const mesh = await makeCircle(40, 1.0);
 * // Apply high-order curvature for better boundary approximation
 * setCurvature(mesh, 3);
 * // Now use for FEM solve...
 * mesh.destroy();
 * ```
 *
 * @see {@link makeDisk} - Alternative disk mesh generator
 * @see {@link setCurvature} - Add high-order curvature for better circle approximation
 */
export async function makeCircle(npts: number, radius: number = 1.0): Promise<Mesh> {
  if (npts < 3) {
    throw new Error('Number of points must be at least 3');
  }

  const module = await loadMFEMModule();

  // Create a square mesh and transform to disk
  const n = Math.max(2, Math.floor(npts / 4));
  const ptr = module._mfem_mesh_make_cartesian_2d(n, n, 1, 2.0, 2.0, 0); // triangles

  if (ptr === 0) {
    throw new Error('Failed to create mesh');
  }

  const mesh = new Mesh(module, ptr);

  // Transform to circular shape using elliptic mapping
  const sdim = mesh.dimension;
  const nv = mesh.numVertices;

  for (let i = 0; i < nv; i++) {
    const coordsPtr = module._malloc(sdim * 8);
    if (coordsPtr === 0) {
      throw new Error('Failed to allocate memory');
    }

    try {
      module._mfem_mesh_get_vertex(mesh.getPointer(), i, coordsPtr);

      const u = module.HEAPF64[coordsPtr >> 3] - 1.0;
      const v = module.HEAPF64[(coordsPtr >> 3) + 1] - 1.0;

      // Elliptic mapping from square to disk
      const x = u * Math.sqrt(1 - 0.5 * v * v) * radius;
      const y = v * Math.sqrt(1 - 0.5 * u * u) * radius;

      module.HEAPF64[coordsPtr >> 3] = x;
      module.HEAPF64[(coordsPtr >> 3) + 1] = y;

      module._mfem_mesh_set_vertex(mesh.getPointer(), i, coordsPtr);
    } finally {
      module._free(coordsPtr);
    }
  }

  return mesh;
}
