/**
 * Utility functions for special functions module
 */

import { NDArray } from 'numwasm';

/**
 * Check if a number is effectively an integer
 */
export function isInteger(x: number): boolean {
  return Number.isInteger(x) && Number.isFinite(x);
}

/**
 * Check if value is a scalar (number)
 */
export function isScalar(value: unknown): value is number {
  return typeof value === 'number';
}

/**
 * Validate that inputs are scalars for exact mode
 */
export function validateScalarInputs(
  N: number | NDArray,
  k: number | NDArray,
  functionName: string
): void {
  if (!isScalar(N) || !isScalar(k)) {
    throw new Error(
      `\`N\` and \`k\` must be scalar integers with \`exact=True\` in ${functionName}.`
    );
  }
}

/**
 * Validate that scalar inputs are integers for exact mode
 */
export function validateIntegerInputs(
  N: number,
  k: number,
  functionName: string
): void {
  if (!isInteger(N) || !isInteger(k)) {
    throw new Error(
      `Non-integer \`N\` and \`k\` with \`exact=True\` is not supported in ${functionName}.`
    );
  }
}

/**
 * Check if result indicates overflow from WASM function
 */
export function checkOverflow(result: number, functionName: string): void {
  if (result === -1.0) {
    throw new Error(
      `Result exceeds maximum representable value in ${functionName} exact mode. ` +
      `Values must fit within unsigned 64-bit integer range (up to 2^64-1).`
    );
  }
}

/**
 * Apply a function element-wise to arrays
 */
export async function elementWise(
  N: NDArray,
  k: NDArray,
  fn: (n: number, k: number) => number
): Promise<NDArray> {
  // Ensure both are NDArrays
  const NArray = N instanceof NDArray ? N : await NDArray.array(N);
  const kArray = k instanceof NDArray ? k : await NDArray.array(k);

  // Get data as typed arrays
  const nData = await NArray.getData();
  const kData = await kArray.getData();

  // Broadcast shapes if needed (for now, assume compatible shapes)
  const size = Math.max(nData.length, kData.length);
  const result = new Float64Array(size);

  for (let i = 0; i < size; i++) {
    const nVal = nData[Math.min(i, nData.length - 1)];
    const kVal = kData[Math.min(i, kData.length - 1)];
    result[i] = fn(nVal, kVal);
  }

  // Create result array with same shape as input (or broadcast shape)
  return NDArray.array(Array.from(result));
}
