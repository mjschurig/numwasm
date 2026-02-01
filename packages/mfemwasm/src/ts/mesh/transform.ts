/**
 * Mesh coordinate transformation.
 *
 * @module mesh
 */

import { Mesh } from './Mesh.js';

/**
 * Transforms all mesh vertex coordinates using a custom function.
 *
 * This function applies an arbitrary coordinate transformation to every vertex
 * in the mesh. It's useful for creating curved domains, deforming meshes, or
 * applying custom mappings.
 *
 * The transformation function receives the current coordinates as an array
 * and must return an array of the same length with the new coordinates.
 *
 * @param mesh - The mesh to transform (modified in place)
 * @param fn - Transformation function that maps old coordinates to new coordinates
 * @throws Error if the mesh has been destroyed
 * @throws Error if memory allocation fails
 * @throws Error if the transformation function returns wrong number of coordinates
 *
 * @category Mesh
 *
 * @example Shift mesh by (1, 2)
 * ```typescript
 * transform(mesh, ([x, y]) => [x + 1, y + 2]);
 * ```
 *
 * @example Scale x-axis by 2
 * ```typescript
 * transform(mesh, ([x, y]) => [2 * x, y]);
 * ```
 *
 * @example Create circular domain from square
 * ```typescript
 * const mesh = await makeCartesian2D(10, 10, 'quad', 2, 2);
 * // Shift to [-1, 1] × [-1, 1]
 * transform(mesh, ([x, y]) => [x - 1, y - 1]);
 * // Apply conformal map to create disk
 * transform(mesh, ([x, y]) => {
 *   const r = Math.sqrt(x*x + y*y);
 *   if (r < 1e-10) return [0, 0];
 *   const scale = Math.min(1, r) / r;
 *   return [x * scale, y * scale];
 * });
 * ```
 *
 * @example Apply shear transformation
 * ```typescript
 * const shearFactor = 0.5;
 * transform(mesh, ([x, y]) => [x + shearFactor * y, y]);
 * ```
 *
 * @example 3D rotation around z-axis
 * ```typescript
 * const angle = Math.PI / 4; // 45 degrees
 * const cos = Math.cos(angle);
 * const sin = Math.sin(angle);
 * transform(mesh, ([x, y, z]) => [
 *   cos * x - sin * y,
 *   sin * x + cos * y,
 *   z
 * ]);
 * ```
 */
export function transform(mesh: Mesh, fn: (coords: number[]) => number[]): void {
  const module = mesh.getModule();
  const ptr = mesh.getPointer();

  const sdim = module._mfem_mesh_get_space_dimension(ptr);
  const nv = mesh.numVertices;

  const coordsPtr = module._malloc(sdim * 8);
  if (coordsPtr === 0) {
    throw new Error('Failed to allocate memory');
  }

  try {
    for (let i = 0; i < nv; i++) {
      module._mfem_mesh_get_vertex(ptr, i, coordsPtr);

      const coords: number[] = [];
      for (let j = 0; j < sdim; j++) {
        coords.push(module.HEAPF64[(coordsPtr >> 3) + j]);
      }

      const newCoords = fn(coords);

      if (newCoords.length !== sdim) {
        throw new Error(`Transformation function must return ${sdim} coordinates`);
      }

      for (let j = 0; j < sdim; j++) {
        module.HEAPF64[(coordsPtr >> 3) + j] = newCoords[j];
      }

      module._mfem_mesh_set_vertex(ptr, i, coordsPtr);
    }
  } finally {
    module._free(coordsPtr);
  }
}
