/**
 * Get integration point at index
 * @param rule - IntegrationRule instance
 * @param i - Point index
 * @returns IntegrationPoint {x, y?, z?, weight}
 */
export function getPoint(
  rule: unknown,
  i: number
): { x: number; y?: number; z?: number; weight: number } {
  throw new Error('Not implemented');
}
