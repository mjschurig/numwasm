/**
 * High-order mesh curvature.
 *
 * @module mesh
 */

import { Mesh } from './Mesh.js';

/**
 * Sets the curvature order for high-order geometry representation.
 *
 * Converts the mesh to use high-order polynomial nodes for representing
 * curved element boundaries. This is essential for accurately representing
 * curved domains (circles, spheres, etc.) and achieving optimal convergence
 * rates with high-order finite elements.
 *
 * Higher orders provide better approximation of curved boundaries but
 * increase the number of geometry nodes per element.
 *
 * @param mesh - The mesh to set curvature on (modified in place)
 * @param order - Polynomial order for geometry representation (≥ 1)
 * @throws Error if order is less than 1
 * @throws Error if setting curvature fails
 *
 * @category Mesh
 *
 * @example Basic curvature for disk mesh
 * ```typescript
 * const mesh = await makeDisk(50);
 * // Linear mesh (order 1) has piecewise linear boundary
 * // Set quadratic curvature for better circle approximation
 * setCurvature(mesh, 2);
 * mesh.destroy();
 * ```
 *
 * @example High-order sphere mesh
 * ```typescript
 * const mesh = await makeSphere(200);
 * // Use cubic geometry for smooth sphere boundary
 * setCurvature(mesh, 3);
 * // Now the mesh boundary closely approximates the true sphere
 * mesh.destroy();
 * ```
 *
 * @example Match geometry order to FE order
 * ```typescript
 * const mesh = await makeCircle(40);
 * // For p=4 finite elements, use order 4 geometry
 * // (isoparametric or super-parametric elements)
 * setCurvature(mesh, 4);
 * mesh.destroy();
 * ```
 *
 * @example Curved cylinder
 * ```typescript
 * const mesh = await makeCylinder(300, 2.0, 1.0);
 * // Add curvature for smooth circular cross-section
 * setCurvature(mesh, 2);
 * // Export for visualization
 * const vtk = toVTK(mesh);
 * mesh.destroy();
 * ```
 *
 * @see {@link makeCircle} - Create circular mesh
 * @see {@link makeDisk} - Create disk mesh
 * @see {@link makeSphere} - Create sphere mesh
 */
export function setCurvature(mesh: Mesh, order: number): void {
  if (order < 1) {
    throw new Error('Curvature order must be at least 1');
  }

  const module = mesh.getModule();
  const ptr = mesh.getPointer();

  const result = module._mfem_mesh_set_curvature(ptr, order);
  if (result !== 0) {
    throw new Error('Failed to set mesh curvature');
  }
}
