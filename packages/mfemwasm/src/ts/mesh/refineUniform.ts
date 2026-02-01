/**
 * Uniform mesh refinement.
 *
 * @module mesh
 */

import { Mesh } from './Mesh.js';

/**
 * Uniformly refines a mesh one or more times.
 *
 * Each refinement level subdivides every element into smaller elements:
 * - 2D triangles → 4 triangles
 * - 2D quadrilaterals → 4 quadrilaterals
 * - 3D tetrahedra → 8 tetrahedra
 * - 3D hexahedra → 8 hexahedra
 *
 * The number of elements grows exponentially: after `times` refinements,
 * a 2D mesh will have 4^times × original elements.
 *
 * @param mesh - The mesh to refine (modified in place)
 * @param times - Number of refinement levels to apply (default 1)
 * @throws Error if the mesh has been destroyed
 * @throws Error if refinement fails
 *
 * @category Mesh
 *
 * @example Single refinement
 * ```typescript
 * const mesh = await makeCartesian2D(4, 4); // 16 elements
 * refineUniform(mesh);
 * console.log(mesh.numElements); // 64 elements
 * ```
 *
 * @example Multiple refinement levels
 * ```typescript
 * const mesh = await makeCartesian2D(2, 2); // 4 elements
 * refineUniform(mesh, 3); // Refine 3 times
 * console.log(mesh.numElements); // 4 × 4³ = 256 elements
 * ```
 *
 * @example Convergence study
 * ```typescript
 * const mesh = await makeCartesian2D(2, 2);
 * for (let level = 0; level < 5; level++) {
 *   console.log(`Level ${level}: ${mesh.numElements} elements`);
 *   // Solve problem and compute error...
 *   if (level < 4) refineUniform(mesh);
 * }
 * ```
 */
export function refineUniform(mesh: Mesh, times: number = 1): void {
  for (let i = 0; i < times; i++) {
    mesh.refineUniform();
  }
}
