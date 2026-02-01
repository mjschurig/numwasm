/**
 * 2D Cartesian mesh generation.
 *
 * @module mesh
 */

import { loadMFEMModule } from '../loader.js';
import { Mesh } from './Mesh.js';

/**
 * Element types available for 2D meshes.
 *
 * - `'quad'` - Quadrilateral elements (4-sided)
 * - `'tri'` - Triangular elements (3-sided)
 *
 * @category Mesh
 */
export type ElementType2D = 'quad' | 'tri';

/**
 * Creates a 2D Cartesian (rectangular) mesh.
 *
 * Generates a structured mesh on the rectangular domain [0, sx] × [0, sy]
 * with the specified number of elements in each direction.
 *
 * @param nx - Number of elements in the x-direction (must be ≥ 1)
 * @param ny - Number of elements in the y-direction (must be ≥ 1)
 * @param type - Element type: `'quad'` for quadrilaterals (default), `'tri'` for triangles
 * @param sx - Domain size in x-direction (default 1.0)
 * @param sy - Domain size in y-direction (default 1.0)
 * @returns A promise that resolves to the created Mesh
 * @throws Error if nx or ny is less than 1
 * @throws Error if mesh creation fails
 *
 * @category Mesh
 *
 * @example Basic 10×10 quadrilateral mesh
 * ```typescript
 * const mesh = await makeCartesian2D(10, 10);
 * console.log(mesh.numElements); // 100
 * console.log(mesh.numVertices); // 121
 * mesh.destroy();
 * ```
 *
 * @example Triangular mesh on rectangular domain
 * ```typescript
 * const mesh = await makeCartesian2D(5, 10, 'tri', 2.0, 4.0);
 * // Creates triangular mesh on [0,2] × [0,4]
 * // 5×10 = 50 quads → 100 triangles
 * console.log(mesh.numElements); // 100
 * mesh.destroy();
 * ```
 *
 * @example High-resolution mesh for accuracy study
 * ```typescript
 * for (const n of [4, 8, 16, 32, 64]) {
 *   const mesh = await makeCartesian2D(n, n);
 *   console.log(`${n}×${n}: ${mesh.numElements} elements`);
 *   mesh.destroy();
 * }
 * ```
 */
export async function makeCartesian2D(
  nx: number,
  ny: number,
  type: ElementType2D = 'quad',
  sx: number = 1.0,
  sy: number = 1.0
): Promise<Mesh> {
  if (nx < 1 || ny < 1) {
    throw new Error('Number of elements must be at least 1 in each direction');
  }

  const module = await loadMFEMModule();
  const typeCode = type === 'quad' ? 0 : 1;
  const ptr = module._mfem_mesh_make_cartesian_2d(nx, ny, typeCode, sx, sy, 0);

  if (ptr === 0) {
    throw new Error('Failed to create 2D Cartesian mesh');
  }

  return new Mesh(module, ptr);
}
