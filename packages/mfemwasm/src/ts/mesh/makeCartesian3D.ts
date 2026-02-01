/**
 * 3D Cartesian mesh generation.
 *
 * @module mesh
 */

import { loadMFEMModule } from '../loader.js';
import { Mesh } from './Mesh.js';

/**
 * Element types available for 3D meshes.
 *
 * - `'hex'` - Hexahedral elements (8 vertices, default)
 * - `'tet'` - Tetrahedral elements (4 vertices)
 * - `'wedge'` - Wedge/prism elements (6 vertices)
 *
 * @category Mesh
 */
export type ElementType3D = 'hex' | 'tet' | 'wedge';

/**
 * Creates a 3D Cartesian (box) mesh.
 *
 * Generates a structured mesh on the box domain [0, sx] × [0, sy] × [0, sz]
 * with the specified number of elements in each direction.
 *
 * @param nx - Number of elements in the x-direction (must be ≥ 1)
 * @param ny - Number of elements in the y-direction (must be ≥ 1)
 * @param nz - Number of elements in the z-direction (must be ≥ 1)
 * @param type - Element type: `'hex'` for hexahedra (default), `'tet'` for tetrahedra, `'wedge'` for wedges
 * @param sx - Domain size in x-direction (default 1.0)
 * @param sy - Domain size in y-direction (default 1.0)
 * @param sz - Domain size in z-direction (default 1.0)
 * @returns A promise that resolves to the created 3D Mesh
 * @throws Error if nx, ny, or nz is less than 1
 * @throws Error if mesh creation fails
 *
 * @category Mesh
 *
 * @example Basic 4×4×4 hexahedral mesh
 * ```typescript
 * const mesh = await makeCartesian3D(4, 4, 4);
 * console.log(mesh.numElements); // 64
 * console.log(mesh.numVertices); // 125
 * console.log(mesh.dimension); // 3
 * mesh.destroy();
 * ```
 *
 * @example Tetrahedral mesh on custom domain
 * ```typescript
 * // Create tetrahedral mesh on [0,2] × [0,1] × [0,0.5]
 * const mesh = await makeCartesian3D(8, 4, 2, 'tet', 2.0, 1.0, 0.5);
 * // Each hex cell is split into 6 tetrahedra
 * console.log(mesh.numElements); // 8 × 4 × 2 × 6 = 384
 * mesh.destroy();
 * ```
 *
 * @example Wedge mesh for layered problems
 * ```typescript
 * // Create wedge mesh for a thin layer problem
 * const mesh = await makeCartesian3D(10, 10, 2, 'wedge', 1.0, 1.0, 0.1);
 * mesh.destroy();
 * ```
 *
 * @example 3D convergence study
 * ```typescript
 * for (const n of [2, 4, 8, 16]) {
 *   const mesh = await makeCartesian3D(n, n, n);
 *   console.log(`${n}³: ${mesh.numElements} elements`);
 *   mesh.destroy();
 * }
 * ```
 */
export async function makeCartesian3D(
  nx: number,
  ny: number,
  nz: number,
  type: ElementType3D = 'hex',
  sx: number = 1.0,
  sy: number = 1.0,
  sz: number = 1.0
): Promise<Mesh> {
  if (nx < 1 || ny < 1 || nz < 1) {
    throw new Error('Number of elements must be at least 1 in each direction');
  }

  const module = await loadMFEMModule();
  let typeCode: number;
  switch (type) {
    case 'tet':
      typeCode = 1;
      break;
    case 'wedge':
      typeCode = 2;
      break;
    default:
      typeCode = 0; // hex
      break;
  }

  const ptr = module._mfem_mesh_make_cartesian_3d(nx, ny, nz, typeCode, sx, sy, sz, 0);

  if (ptr === 0) {
    throw new Error('Failed to create 3D Cartesian mesh');
  }

  return new Mesh(module, ptr);
}
