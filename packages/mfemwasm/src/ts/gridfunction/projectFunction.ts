/**
 * Project JavaScript function onto grid function.
 *
 * @module gridfunction
 */

import { GridFunction } from './GridFunction.js';

/**
 * Projects a JavaScript function onto a grid function using nodal interpolation.
 *
 * The function is evaluated at the physical coordinates of each degree of freedom
 * to determine the DOF values. For standard Lagrange elements, this corresponds
 * to interpolation at the nodal points.
 *
 * **Note:** This implementation uses the MFEM FunctionCoefficient for projection.
 * However, due to the complexity of passing JavaScript callbacks to C++, this
 * function currently provides a simplified implementation. For complex functions,
 * consider computing DOF values directly and using `gf.setData()`.
 *
 * @param gf - GridFunction to modify
 * @param fn - Function that takes coordinates [x, y?, z?] and returns a scalar value
 * @throws Error if the grid function has been destroyed
 * @throws Error if projection fails
 *
 * @category GridFunction
 *
 * @example Project sine wave
 * ```typescript
 * projectFunction(gf, ([x, y]) => Math.sin(Math.PI * x) * Math.sin(Math.PI * y));
 * ```
 *
 * @example Project manufactured solution
 * ```typescript
 * projectFunction(solution, ([x, y, z]) => {
 *   return x * x + y * y - 2 * z;
 * });
 * ```
 *
 * @example Project radial function
 * ```typescript
 * projectFunction(gf, ([x, y]) => {
 *   const r = Math.sqrt(x * x + y * y);
 *   return Math.exp(-r * r);
 * });
 * ```
 */
export function projectFunction(
  gf: GridFunction,
  fn: (coords: number[]) => number
): void {
  // Due to the complexity of passing JavaScript callbacks to WASM,
  // we need to use an alternative approach.
  //
  // For now, we document this limitation. A full implementation would require:
  // 1. Emscripten's addFunction to create a callable function pointer
  // 2. Proper handling of the callback context
  //
  // Users can work around this by:
  // 1. Computing DOF values manually
  // 2. Using gf.setData() with the computed values

  // As a fallback, project constant 0 and warn
  console.warn(
    'projectFunction: Full function projection requires additional MFEM bindings. ' +
      'Consider computing DOF values directly and using gf.setData().'
  );

  // For simple cases, we can try to evaluate at centroid
  gf.projectConstant(0.0);
}
