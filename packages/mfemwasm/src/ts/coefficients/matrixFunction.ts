/**
 * Create matrix coefficient from JavaScript function
 * @param fn - Function (x: number[], t?: number) => number[][]
 * @returns MatrixCoefficient instance
 */
export function matrixFunction(
  fn: (x: number[], t?: number) => number[][]
): unknown {
  throw new Error('Not implemented');
}
