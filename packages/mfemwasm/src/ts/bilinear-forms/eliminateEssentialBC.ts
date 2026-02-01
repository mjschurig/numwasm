/**
 * Eliminate essential BCs from system
 * @param blf - BilinearForm instance
 * @param essentialDofs - Essential DOF indices
 * @param sol - Solution grid function
 * @param rhs - Right-hand side linear form
 */
export function eliminateEssentialBC(
  blf: unknown,
  essentialDofs: number[],
  sol: unknown,
  rhs: unknown
): void {
  throw new Error('Not implemented');
}
