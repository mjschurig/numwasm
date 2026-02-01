/**
 * 1D Cartesian mesh generation.
 *
 * @module mesh
 */

import { loadMFEMModule } from '../loader.js';
import { Mesh } from './Mesh.js';

/**
 * Creates a 1D Cartesian (line segment) mesh.
 *
 * Generates a uniform mesh on the interval [0, sx] with `n` line segment elements.
 * This is useful for 1D finite element problems like ODEs, 1D heat transfer,
 * or testing numerical methods.
 *
 * @param n - Number of elements (must be ≥ 1)
 * @param sx - Domain size / length of the interval (default 1.0)
 * @returns A promise that resolves to the created 1D Mesh
 * @throws Error if n is less than 1
 * @throws Error if mesh creation fails
 *
 * @category Mesh
 *
 * @example Basic 10-element mesh on [0, 1]
 * ```typescript
 * const mesh = await makeCartesian1D(10);
 * console.log(mesh.numElements); // 10
 * console.log(mesh.numVertices); // 11
 * console.log(mesh.dimension); // 1
 * mesh.destroy();
 * ```
 *
 * @example Custom domain length
 * ```typescript
 * // Mesh on [0, 2π] for periodic problems
 * const mesh = await makeCartesian1D(100, 2 * Math.PI);
 * mesh.destroy();
 * ```
 *
 * @example Convergence study
 * ```typescript
 * for (const n of [10, 20, 40, 80]) {
 *   const mesh = await makeCartesian1D(n);
 *   const h = 1.0 / n; // Element size
 *   console.log(`n=${n}, h=${h}`);
 *   // Solve problem and compute error...
 *   mesh.destroy();
 * }
 * ```
 */
export async function makeCartesian1D(n: number, sx: number = 1.0): Promise<Mesh> {
  if (n < 1) {
    throw new Error('Number of elements must be at least 1');
  }

  const module = await loadMFEMModule();
  const ptr = module._mfem_mesh_make_cartesian_1d(n, sx);

  if (ptr === 0) {
    throw new Error('Failed to create 1D Cartesian mesh');
  }

  return new Mesh(module, ptr);
}
