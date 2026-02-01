/**
 * Disk mesh generation.
 *
 * @module mesh
 */

import { loadMFEMModule } from '../loader.js';
import { Mesh } from './Mesh.js';

/**
 * Creates a 2D disk mesh using elliptic mapping.
 *
 * Generates a triangulated mesh of a circular disk by mapping a square mesh
 * to a disk using the elliptic conformal mapping:
 * ```
 * x = u × √(1 - v²/2) × radius
 * y = v × √(1 - u²/2) × radius
 * ```
 *
 * This mapping produces high-quality elements with smooth variation across
 * the domain, making it suitable for finite element analysis.
 *
 * @param npts - Controls mesh density via √npts elements per side (minimum 3)
 * @param radius - Radius of the disk (default 1.0)
 * @returns A promise that resolves to a 2D triangular Mesh
 * @throws Error if npts is less than 3
 * @throws Error if mesh creation fails
 *
 * @category Mesh
 *
 * @example Basic unit disk mesh
 * ```typescript
 * const mesh = await makeDisk(100); // ~10×10 elements
 * console.log(mesh.dimension); // 2
 * mesh.destroy();
 * ```
 *
 * @example Disk for Poisson equation
 * ```typescript
 * const mesh = await makeDisk(400, 1.0); // ~20×20 mesh
 * // Refine for better accuracy
 * refineUniform(mesh, 2);
 * // Use for solving -Δu = f on disk
 * mesh.destroy();
 * ```
 *
 * @example Annular domain (disk with hole)
 * ```typescript
 * const mesh = await makeDisk(200);
 * // Transform to create annulus: keep only r > 0.3
 * // (requires custom filtering or use mesh attributes)
 * mesh.destroy();
 * ```
 *
 * @see {@link makeCircle} - Alternative disk mesh generator
 * @see {@link makeSphere} - 3D sphere mesh generator
 */
export async function makeDisk(npts: number, radius: number = 1.0): Promise<Mesh> {
  if (npts < 3) {
    throw new Error('Number of points must be at least 3');
  }

  const module = await loadMFEMModule();

  // Create a square mesh and transform to disk
  const n = Math.max(2, Math.floor(Math.sqrt(npts)));
  const ptr = module._mfem_mesh_make_cartesian_2d(n, n, 1, 2.0, 2.0, 0); // triangles

  if (ptr === 0) {
    throw new Error('Failed to create mesh');
  }

  const mesh = new Mesh(module, ptr);

  // Transform square [-1,1]^2 to unit disk using elliptic mapping
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
