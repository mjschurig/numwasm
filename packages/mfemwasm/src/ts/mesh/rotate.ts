/**
 * Mesh rotation transformation.
 *
 * @module mesh
 */

import { Mesh } from './Mesh.js';

/**
 * Rotates a mesh around an axis.
 *
 * For 2D meshes, rotates around the origin (z-axis).
 * For 3D meshes, rotates around the specified coordinate axis passing
 * through the origin.
 *
 * Rotation follows the right-hand rule: positive angles rotate
 * counterclockwise when looking along the axis toward the origin.
 *
 * @param mesh - The mesh to rotate (modified in place)
 * @param angle - Rotation angle in radians (positive = counterclockwise)
 * @param axis - Rotation axis for 3D: 0=x, 1=y, 2=z (default). Ignored for 2D.
 * @throws Error if rotation fails
 * @throws Error if mesh is 1D (rotation not supported)
 *
 * @category Mesh
 *
 * @example Rotate 2D mesh by 45 degrees
 * ```typescript
 * const mesh = await makeCartesian2D(4, 4);
 * // Center the mesh first
 * translate(mesh, [-0.5, -0.5]);
 * // Rotate 45 degrees (π/4 radians)
 * rotate(mesh, Math.PI / 4);
 * mesh.destroy();
 * ```
 *
 * @example Rotate 3D mesh around z-axis (default)
 * ```typescript
 * const mesh = await makeCartesian3D(3, 3, 3);
 * translate(mesh, [-0.5, -0.5, -0.5]); // Center at origin
 * rotate(mesh, Math.PI / 6); // Rotate 30 degrees around z
 * mesh.destroy();
 * ```
 *
 * @example Rotate around different axes
 * ```typescript
 * const mesh = await makeCylinder(100, 2.0, 0.5);
 * // Cylinder is along z-axis by default
 *
 * // Rotate 90 degrees around x-axis to make it along y
 * rotate(mesh, Math.PI / 2, 0);
 *
 * // Or rotate around y-axis to make it along x
 * // rotate(mesh, Math.PI / 2, 1);
 * mesh.destroy();
 * ```
 *
 * @example Create rotated copies of a mesh
 * ```typescript
 * async function createRotatedMesh(angleDegrees: number) {
 *   const mesh = await makeDisk(50, 1.0);
 *   translate(mesh, [2.0, 0]); // Move away from origin
 *   rotate(mesh, angleDegrees * Math.PI / 180);
 *   return mesh;
 * }
 *
 * // Create meshes at 0°, 90°, 180°, 270°
 * const meshes = await Promise.all([0, 90, 180, 270].map(createRotatedMesh));
 * meshes.forEach(m => m.destroy());
 * ```
 *
 * @see {@link translate} - Translate mesh
 * @see {@link scale} - Scale mesh
 * @see {@link transform} - Apply arbitrary transformation
 */
export function rotate(mesh: Mesh, angle: number, axis: number = 2): void {
  const module = mesh.getModule();
  const ptr = mesh.getPointer();
  const dim = mesh.dimension;

  if (dim === 2) {
    // 2D rotation
    const result = module._mfem_mesh_rotate_2d(ptr, angle);
    if (result !== 0) {
      throw new Error('Failed to rotate 2D mesh');
    }
  } else if (dim === 3) {
    // 3D rotation - apply via vertex transformation
    const cos = Math.cos(angle);
    const sin = Math.sin(angle);

    const sdim = module._mfem_mesh_get_space_dimension(ptr);
    const nv = mesh.numVertices;

    const coordsPtr = module._malloc(sdim * 8);
    if (coordsPtr === 0) {
      throw new Error('Failed to allocate memory');
    }

    try {
      for (let i = 0; i < nv; i++) {
        module._mfem_mesh_get_vertex(ptr, i, coordsPtr);

        const x = module.HEAPF64[(coordsPtr >> 3)];
        const y = module.HEAPF64[(coordsPtr >> 3) + 1];
        const z = module.HEAPF64[(coordsPtr >> 3) + 2];

        let nx: number, ny: number, nz: number;

        if (axis === 0) {
          // Rotate around x-axis
          nx = x;
          ny = cos * y - sin * z;
          nz = sin * y + cos * z;
        } else if (axis === 1) {
          // Rotate around y-axis
          nx = cos * x + sin * z;
          ny = y;
          nz = -sin * x + cos * z;
        } else {
          // Rotate around z-axis (default)
          nx = cos * x - sin * y;
          ny = sin * x + cos * y;
          nz = z;
        }

        module.HEAPF64[(coordsPtr >> 3)] = nx;
        module.HEAPF64[(coordsPtr >> 3) + 1] = ny;
        module.HEAPF64[(coordsPtr >> 3) + 2] = nz;

        module._mfem_mesh_set_vertex(ptr, i, coordsPtr);
      }
    } finally {
      module._free(coordsPtr);
    }
  } else {
    throw new Error('Rotation only supported for 2D and 3D meshes');
  }
}
