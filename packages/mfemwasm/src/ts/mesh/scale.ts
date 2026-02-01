/**
 * Mesh scaling transformation.
 *
 * @module mesh
 */

import { Mesh } from './Mesh.js';

/**
 * Scales a mesh uniformly or along each axis.
 *
 * Multiplies all vertex coordinates by the scale factor(s), effectively
 * resizing the mesh. Scaling is performed relative to the origin.
 *
 * @param mesh - The mesh to scale (modified in place)
 * @param factor - Either a single uniform scale factor, or an array of
 *                 per-axis factors [sx, sy] (2D) or [sx, sy, sz] (3D)
 * @throws Error if factor is not a number or array
 * @throws Error if scaling fails
 *
 * @category Mesh
 *
 * @example Uniform scaling
 * ```typescript
 * const mesh = await makeCartesian2D(4, 4); // Domain [0,1] × [0,1]
 * scale(mesh, 2.0);
 * // Domain is now [0,2] × [0,2]
 * mesh.destroy();
 * ```
 *
 * @example Scale to match physical dimensions
 * ```typescript
 * // Create unit square mesh
 * const mesh = await makeCartesian2D(10, 10);
 * // Scale to 5m × 3m physical domain
 * scale(mesh, [5.0, 3.0]);
 * mesh.destroy();
 * ```
 *
 * @example Non-uniform 3D scaling
 * ```typescript
 * const mesh = await makeCartesian3D(4, 4, 4);
 * // Stretch in z-direction (thin plate geometry)
 * scale(mesh, [1.0, 1.0, 0.1]);
 * mesh.destroy();
 * ```
 *
 * @example Scale and translate to specific location
 * ```typescript
 * const mesh = await makeCartesian2D(5, 5);
 * // Scale to [0,2] × [0,2]
 * scale(mesh, 2.0);
 * // Then translate to center at origin
 * translate(mesh, [-1.0, -1.0]);
 * // Domain is now [-1,1] × [-1,1]
 * mesh.destroy();
 * ```
 *
 * @see {@link translate} - Translate mesh by offset
 * @see {@link rotate} - Rotate mesh
 * @see {@link transform} - Apply arbitrary coordinate transformation
 */
export function scale(mesh: Mesh, factor: number | number[]): void {
  const module = mesh.getModule();
  const ptr = mesh.getPointer();

  let sx: number, sy: number, sz: number;

  if (typeof factor === 'number') {
    sx = sy = sz = factor;
  } else if (Array.isArray(factor)) {
    sx = factor[0] ?? 1.0;
    sy = factor[1] ?? 1.0;
    sz = factor[2] ?? 1.0;
  } else {
    throw new Error('Factor must be a number or array');
  }

  const result = module._mfem_mesh_scale(ptr, sx, sy, sz);
  if (result !== 0) {
    throw new Error('Failed to scale mesh');
  }
}
