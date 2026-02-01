/**
 * Project coefficient on boundary DOFs.
 *
 * @module gridfunction
 */

import { GridFunction } from './GridFunction.js';

/**
 * Projects a constant value onto boundary DOFs for specified boundary attributes.
 *
 * This sets the DOF values on the specified boundaries to the given constant value.
 * It is useful for applying Dirichlet boundary conditions where the boundary
 * value is known.
 *
 * Only DOFs that lie on boundary elements with the specified attributes are modified.
 * Interior DOFs remain unchanged.
 *
 * @param gf - GridFunction to modify
 * @param value - The constant value to project on boundary
 * @param bdrAttrs - Array of boundary attributes (1-based) to apply the value to
 * @throws Error if the grid function has been destroyed
 * @throws Error if projection fails
 *
 * @category GridFunction
 *
 * @example Set Dirichlet boundary condition
 * ```typescript
 * // Set temperature to 100 on inlet (boundary attribute 1)
 * projectBdrCoefficient(temperature, 100.0, [1]);
 * ```
 *
 * @example Set homogeneous BC on multiple boundaries
 * ```typescript
 * // Set u = 0 on all walls (attributes 1, 2, 3, 4)
 * projectBdrCoefficient(displacement, 0.0, [1, 2, 3, 4]);
 * ```
 *
 * @example Set different values on different boundaries
 * ```typescript
 * // Inlet at T = 100
 * projectBdrCoefficient(temperature, 100.0, [1]);
 * // Outlet at T = 20
 * projectBdrCoefficient(temperature, 20.0, [2]);
 * ```
 */
export function projectBdrCoefficient(
  gf: GridFunction,
  value: number,
  bdrAttrs: number[]
): void {
  const module = gf.getModule();
  const ptr = gf.getPointer();

  for (const attr of bdrAttrs) {
    const result = module._mfem_gridfunc_project_bdr_coefficient(ptr, attr, value);
    if (result < 0) {
      throw new Error(`Failed to project boundary coefficient on attribute ${attr}`);
    }
  }
}
