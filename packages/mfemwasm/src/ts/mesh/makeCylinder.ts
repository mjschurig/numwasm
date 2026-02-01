/**
 * Cylinder mesh generation.
 *
 * @module mesh
 */

import { loadMFEMModule } from '../loader.js';
import { Mesh } from './Mesh.js';

/**
 * Creates a 3D solid cylinder mesh.
 *
 * Generates a hexahedral mesh of a solid cylinder by creating a box mesh
 * and applying an elliptic conformal mapping to the x-y cross-section while
 * preserving the z-coordinate. The cylinder is centered on the z-axis with
 * z ∈ [0, length].
 *
 * @param npts - Controls mesh density via ∛npts elements per side (minimum 4)
 * @param length - Length of the cylinder in the z-direction (default 1.0)
 * @param radius - Radius of the cylinder in the x-y plane (default 1.0)
 * @returns A promise that resolves to a 3D hexahedral Mesh
 * @throws Error if npts is less than 4
 * @throws Error if mesh creation fails
 *
 * @category Mesh
 *
 * @example Basic unit cylinder mesh
 * ```typescript
 * const mesh = await makeCylinder(125); // ~5×5×5 base mesh
 * // Cylinder along z-axis, radius=1, length=1
 * console.log(mesh.dimension); // 3
 * mesh.destroy();
 * ```
 *
 * @example Long thin cylinder (pipe)
 * ```typescript
 * // Pipe with length 10, radius 0.5
 * const mesh = await makeCylinder(500, 10.0, 0.5);
 * mesh.destroy();
 * ```
 *
 * @example Cylinder for axisymmetric problems
 * ```typescript
 * const mesh = await makeCylinder(300, 2.0, 1.0);
 * // Refine near boundaries for better resolution
 * refineUniform(mesh);
 * // Use for solving PDEs on cylindrical domains
 * mesh.destroy();
 * ```
 *
 * @example High-order curved cylinder
 * ```typescript
 * const mesh = await makeCylinder(200);
 * // Apply curvature for better circular boundary approximation
 * setCurvature(mesh, 3);
 * mesh.destroy();
 * ```
 *
 * @see {@link makeSphere} - 3D sphere mesh
 * @see {@link makeDisk} - 2D disk mesh (cylinder cross-section)
 */
export async function makeCylinder(
  npts: number,
  length: number = 1.0,
  radius: number = 1.0
): Promise<Mesh> {
  if (npts < 4) {
    throw new Error('Number of points must be at least 4');
  }

  const module = await loadMFEMModule();

  // Create a box mesh and transform to cylinder
  const n = Math.max(2, Math.floor(Math.cbrt(npts)));
  const ptr = module._mfem_mesh_make_cartesian_3d(n, n, n, 0, 2.0, 2.0, length, 0); // hexahedra

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
      const z = module.HEAPF64[(coordsPtr >> 3) + 2]; // Keep z as-is

      // Elliptic mapping from square to disk for x-y cross section
      const x = u * Math.sqrt(1 - 0.5 * v * v) * radius;
      const y = v * Math.sqrt(1 - 0.5 * u * u) * radius;

      module.HEAPF64[coordsPtr >> 3] = x;
      module.HEAPF64[(coordsPtr >> 3) + 1] = y;
      module.HEAPF64[(coordsPtr >> 3) + 2] = z;

      module._mfem_mesh_set_vertex(mesh.getPointer(), i, coordsPtr);
    } finally {
      module._free(coordsPtr);
    }
  }

  return mesh;
}
