/**
 * Mesh translation transformation.
 *
 * @module mesh
 */

import { Mesh } from './Mesh.js';

/**
 * Translates a mesh by an offset vector.
 *
 * Adds the offset to all vertex coordinates, effectively moving the entire
 * mesh in space without changing its shape or size.
 *
 * @param mesh - The mesh to translate (modified in place)
 * @param offset - Translation vector [dx, dy] for 2D or [dx, dy, dz] for 3D
 * @throws Error if translation fails
 *
 * @category Mesh
 *
 * @example Center mesh at origin (2D)
 * ```typescript
 * const mesh = await makeCartesian2D(4, 4); // Domain [0,1] × [0,1]
 * translate(mesh, [-0.5, -0.5]);
 * // Domain is now [-0.5,0.5] × [-0.5,0.5], centered at origin
 * mesh.destroy();
 * ```
 *
 * @example Move mesh to specific location
 * ```typescript
 * const mesh = await makeDisk(50, 1.0); // Unit disk centered at origin
 * // Move disk to center at (5, 3)
 * translate(mesh, [5.0, 3.0]);
 * mesh.destroy();
 * ```
 *
 * @example Combine multiple meshes (conceptually)
 * ```typescript
 * // Create two spheres at different locations
 * const sphere1 = await makeSphere(200, 0.5);
 * translate(sphere1, [-1.0, 0.0, 0.0]);
 *
 * const sphere2 = await makeSphere(200, 0.5);
 * translate(sphere2, [1.0, 0.0, 0.0]);
 *
 * // Now sphere1 is centered at (-1, 0, 0)
 * // and sphere2 is centered at (1, 0, 0)
 * sphere1.destroy();
 * sphere2.destroy();
 * ```
 *
 * @example 3D translation
 * ```typescript
 * const mesh = await makeCartesian3D(3, 3, 3);
 * // Move the cube up by 2 units in z
 * translate(mesh, [0, 0, 2.0]);
 * mesh.destroy();
 * ```
 *
 * @see {@link scale} - Scale mesh dimensions
 * @see {@link rotate} - Rotate mesh
 * @see {@link transform} - Apply arbitrary coordinate transformation
 */
export function translate(mesh: Mesh, offset: number[]): void {
  const module = mesh.getModule();
  const ptr = mesh.getPointer();

  const dx = offset[0] ?? 0.0;
  const dy = offset[1] ?? 0.0;
  const dz = offset[2] ?? 0.0;

  const result = module._mfem_mesh_translate(ptr, dx, dy, dz);
  if (result !== 0) {
    throw new Error('Failed to translate mesh');
  }
}
