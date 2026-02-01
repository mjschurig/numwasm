/**
 * Sphere/ball mesh generation.
 *
 * @module mesh
 */

import { loadMFEMModule } from '../loader.js';
import { Mesh } from './Mesh.js';

/**
 * Creates a 3D solid sphere (ball) mesh.
 *
 * Generates a tetrahedral mesh of a solid ball by creating a cube mesh
 * and applying a radial projection transformation. Interior points are
 * mapped proportionally while boundary points are projected to the sphere surface.
 *
 * The mesh is centered at the origin with the specified radius.
 *
 * @param npts - Controls mesh density via ∛npts elements per side (minimum 4)
 * @param radius - Radius of the sphere (default 1.0)
 * @returns A promise that resolves to a 3D tetrahedral Mesh
 * @throws Error if npts is less than 4
 * @throws Error if mesh creation fails
 *
 * @category Mesh
 *
 * @example Basic unit sphere mesh
 * ```typescript
 * const mesh = await makeSphere(125); // ~5×5×5 base mesh
 * console.log(mesh.dimension); // 3
 * mesh.destroy();
 * ```
 *
 * @example Sphere for heat equation
 * ```typescript
 * const mesh = await makeSphere(1000, 1.0);
 * // Refine for better accuracy
 * refineUniform(mesh);
 * // Solve heat equation on ball domain
 * mesh.destroy();
 * ```
 *
 * @example High-order curved sphere
 * ```typescript
 * const mesh = await makeSphere(500);
 * // Apply high-order curvature for accurate sphere boundary
 * setCurvature(mesh, 4);
 * // Now the boundary better approximates the true sphere
 * mesh.destroy();
 * ```
 *
 * @see {@link makeDisk} - 2D disk mesh
 * @see {@link makeCylinder} - 3D cylinder mesh
 * @see {@link setCurvature} - Add high-order curvature
 */
export async function makeSphere(npts: number, radius: number = 1.0): Promise<Mesh> {
  if (npts < 4) {
    throw new Error('Number of points must be at least 4');
  }

  const module = await loadMFEMModule();

  // Create a cube mesh and transform to ball
  const n = Math.max(2, Math.floor(Math.cbrt(npts)));
  const ptr = module._mfem_mesh_make_cartesian_3d(n, n, n, 1, 2.0, 2.0, 2.0, 0); // tetrahedra

  if (ptr === 0) {
    throw new Error('Failed to create mesh');
  }

  const mesh = new Mesh(module, ptr);

  // Get space dimension
  const sdim = module._mfem_mesh_get_space_dimension(mesh.getPointer());
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
      const w = module.HEAPF64[(coordsPtr >> 3) + 2] - 1.0;

      // Normalize to sphere surface, but scale by original distance for interior points
      const r = Math.sqrt(u * u + v * v + w * w);
      if (r > 1e-10) {
        const scale = radius * Math.min(1.0, r) / r;
        module.HEAPF64[coordsPtr >> 3] = u * scale;
        module.HEAPF64[(coordsPtr >> 3) + 1] = v * scale;
        module.HEAPF64[(coordsPtr >> 3) + 2] = w * scale;
      }

      module._mfem_mesh_set_vertex(mesh.getPointer(), i, coordsPtr);
    } finally {
      module._free(coordsPtr);
    }
  }

  return mesh;
}
