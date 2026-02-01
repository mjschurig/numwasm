/**
 * Add Robin boundary condition integrator (α*u + β*∂u/∂n = g)
 * @param blf - BilinearForm instance
 * @param alpha - Robin alpha coefficient
 * @param beta - Robin beta coefficient
 * @param bdrAttrs - Boundary attributes
 */
export function addBoundaryRobinIntegrator(
  blf: unknown,
  alpha: unknown,
  beta: unknown,
  bdrAttrs?: number[]
): void {
  throw new Error('Not implemented');
}
