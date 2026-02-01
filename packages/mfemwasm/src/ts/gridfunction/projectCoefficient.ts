/**
 * Project coefficient onto grid function.
 *
 * @module gridfunction
 */

import { GridFunction } from './GridFunction.js';

/**
 * Projects a constant coefficient onto a grid function.
 *
 * This sets all DOF values such that the grid function represents
 * the constant function f(x) = value everywhere in the domain.
 * This is equivalent to calling `gf.projectConstant(value)`.
 *
 * @param gf - GridFunction to modify
 * @param value - The constant value to project
 * @throws Error if the grid function has been destroyed
 * @throws Error if projection fails
 *
 * @category GridFunction
 *
 * @example Set initial condition
 * ```typescript
 * const temperature = await GridFunction.create(fespace);
 * projectCoefficient(temperature, 20.0); // Room temperature
 * ```
 *
 * @example Reset solution to zero
 * ```typescript
 * projectCoefficient(solution, 0.0);
 * ```
 */
export function projectCoefficient(gf: GridFunction, value: number): void {
  gf.projectConstant(value);
}
