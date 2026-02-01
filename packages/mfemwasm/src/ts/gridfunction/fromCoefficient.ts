/**
 * Create grid function from a constant coefficient.
 *
 * @module gridfunction
 */

import type { FiniteElementSpace } from '../fespace/FiniteElementSpace.js';
import { GridFunction } from './GridFunction.js';

/**
 * Creates a new grid function initialized with a constant value.
 *
 * This is a convenience function that creates a grid function and projects
 * a constant coefficient onto it. For more complex coefficients, create a
 * grid function manually and use projection functions.
 *
 * @param fespace - The finite element space on which to define the function
 * @param value - The constant value to project
 * @returns A promise resolving to the new GridFunction
 * @throws Error if creation or projection fails
 *
 * @category GridFunction
 *
 * @example Initialize temperature field to constant value
 * ```typescript
 * const temperature = await fromCoefficient(fespace, 20.0);
 * console.log(`Initialized ${temperature.size} DOFs to 20.0`);
 * ```
 *
 * @example Initialize vector field components
 * ```typescript
 * const fespace = await createH1(mesh, 2, 3); // 3D vector field
 * const velocity = await fromCoefficient(fespace, 0.0);
 * ```
 */
export async function fromCoefficient(
  fespace: FiniteElementSpace,
  value: number
): Promise<GridFunction> {
  const gf = await GridFunction.create(fespace);

  try {
    gf.projectConstant(value);
    return gf;
  } catch (error) {
    gf.destroy();
    throw error;
  }
}
