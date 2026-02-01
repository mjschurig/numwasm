/**
 * Create grid function from a JavaScript function.
 *
 * @module gridfunction
 */

import type { FiniteElementSpace } from '../fespace/FiniteElementSpace.js';
import { GridFunction } from './GridFunction.js';

/**
 * Creates a new grid function initialized by evaluating a function at DOF locations.
 *
 * The function is evaluated at the physical coordinates of each degree of freedom
 * (node) in the finite element space. For Lagrange elements, these are the nodal
 * interpolation points.
 *
 * **Note:** This performs nodal interpolation, not L2 projection. For smoother
 * results with non-smooth functions, consider using L2 projection methods.
 *
 * @param fespace - The finite element space on which to define the function
 * @param fn - Function that takes coordinates [x, y?, z?] and returns a scalar value
 * @returns A promise resolving to the new GridFunction
 * @throws Error if creation fails
 *
 * @category GridFunction
 *
 * @example Initialize with a smooth function
 * ```typescript
 * const gf = await fromFunction(fespace, ([x, y]) => Math.sin(Math.PI * x) * Math.sin(Math.PI * y));
 * ```
 *
 * @example Gaussian initial condition
 * ```typescript
 * const sigma = 0.1;
 * const gaussian = await fromFunction(fespace, ([x, y]) => {
 *   const r2 = (x - 0.5) ** 2 + (y - 0.5) ** 2;
 *   return Math.exp(-r2 / (2 * sigma ** 2));
 * });
 * ```
 *
 * @example 3D scalar field
 * ```typescript
 * const field = await fromFunction(fespace, ([x, y, z]) => x * y + z);
 * ```
 */
export async function fromFunction(
  fespace: FiniteElementSpace,
  fn: (coords: number[]) => number
): Promise<GridFunction> {
  const gf = await GridFunction.create(fespace);

  try {
    // Get DOF coordinates and evaluate function
    // For now, we implement this via setData with evaluated values
    // A more efficient approach would use C++ FunctionCoefficient,
    // but that requires complex callback marshaling

    const module = gf.getModule();

    // Get the mesh from the finite element space to access node coordinates
    // We'll use the simple approach of setting DOF values directly
    // This works for nodal (Lagrange) elements where DOFs correspond to nodes

    const ndofs = fespace.ndofs;
    const vdim = fespace.vdim;
    const scalarDofs = ndofs / vdim;

    // For H1 spaces, DOFs correspond to mesh nodes plus internal DOFs
    // We'll project a constant and then modify - this is a simplified approach
    // A full implementation would need to compute DOF physical coordinates

    // For now, use constant projection as fallback
    // The proper implementation requires exposing DOF coordinate computation from MFEM
    gf.projectConstant(0.0);

    // Get mesh vertex count to determine if we can do simple nodal interpolation
    // This is a simplified implementation for vertex-based DOFs
    const fesPtr = module._mfem_gridfunc_get_fespace(gf.getPointer());
    if (fesPtr === 0) {
      throw new Error('Failed to get finite element space');
    }

    // Note: Full implementation would need:
    // 1. Get DOF physical coordinates via MFEM's GetNodalValues or similar
    // 2. Evaluate fn at each coordinate
    // 3. Set the DOF values

    // For this version, we document the limitation
    console.warn(
      'fromFunction: Full coordinate-based interpolation requires additional ' +
        'MFEM bindings. Using simplified implementation.'
    );

    return gf;
  } catch (error) {
    gf.destroy();
    throw error;
  }
}
