/**
 * Apply essential (Dirichlet) boundary conditions
 * @param fespace - FiniteElementSpace
 * @param bdrAttrs - Boundary attributes
 * @param value - Boundary value (constant or coefficient)
 * @param x - GridFunction to modify
 */
export function applyEssentialBC(
  fespace: unknown,
  bdrAttrs: number[],
  value: number | unknown,
  x: unknown
): void {
  throw new Error('Not implemented');
}
