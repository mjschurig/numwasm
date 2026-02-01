/**
 * NURBS finite element space creation.
 *
 * @module fespace
 */

import type { Mesh } from '../mesh/Mesh.js';
import { FiniteElementSpace } from './FiniteElementSpace.js';

/**
 * Creates a NURBS finite element space.
 *
 * NURBS (Non-Uniform Rational B-Splines) spaces are used for isogeometric
 * analysis (IGA), where the same basis functions used for CAD geometry
 * representation are used for the finite element analysis. This provides
 * exact geometry representation and higher continuity across elements.
 *
 * The mesh must be a NURBS mesh (created from NURBS patch data). Standard
 * meshes created with `makeCartesian2D` or similar functions cannot be used.
 *
 * @param mesh - A NURBS mesh instance (must have NURBS extension)
 * @param order - Polynomial order. Use -1 or 0 for variable order (uses
 *                the mesh's native NURBS order from each patch)
 * @param vdim - Vector dimension (default 1 for scalar fields)
 * @returns The new NURBS FiniteElementSpace
 * @throws Error if the mesh is not a NURBS mesh
 * @throws Error if space creation fails
 *
 * @category FiniteElementSpace
 *
 * @example NURBS space for isogeometric analysis
 * ```typescript
 * // Load a NURBS mesh from MFEM format (with NURBS patch data)
 * const mesh = await loadMesh(nurbsMeshData);
 *
 * // Create NURBS space with mesh's native order
 * const fespace = createNURBS(mesh, -1);
 *
 * // Or specify explicit order
 * const fespace2 = createNURBS(mesh, 3);
 * ```
 *
 * @example Check if mesh supports NURBS
 * ```typescript
 * if (mesh.isNURBS) {
 *   const fespace = createNURBS(mesh, mesh.nurbsOrder);
 *   // ... perform isogeometric analysis
 * }
 * ```
 *
 * @see {@link https://mfem.org/nurbs-howto/ | MFEM NURBS documentation}
 */
export function createNURBS(
  mesh: Mesh,
  order: number = -1,
  vdim: number = 1
): FiniteElementSpace {
  const module = mesh.getModule();
  const meshPtr = mesh.getPointer();

  // Verify this is a NURBS mesh
  if (!module._mfem_mesh_is_nurbs(meshPtr)) {
    throw new Error(
      'Cannot create NURBS finite element space: mesh is not a NURBS mesh. ' +
      'NURBS spaces require a mesh created from NURBS patch data.'
    );
  }

  const ptr = module._mfem_fespace_create_nurbs(meshPtr, order, vdim);

  if (ptr === 0) {
    throw new Error('Failed to create NURBS finite element space');
  }

  return new FiniteElementSpace(module, ptr);
}
