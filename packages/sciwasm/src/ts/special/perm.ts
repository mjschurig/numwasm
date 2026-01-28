/**
 * Permutations: P(N, k) = N! / (N-k)!
 *
 * Computes the number of ways to arrange k items from N items,
 * where order matters. Also known as k-permutations of N or partial permutations.
 */

import { NDArray } from 'numwasm';
import { getWasmModule } from '../wasm-loader.js';
import {
  isScalar,
  validateScalarInputs,
  validateIntegerInputs,
  checkOverflow,
  elementWise,
} from './_utils.js';

export interface PermOptions {
  /**
   * If true, compute exactly using integer arithmetic.
   * Only supports scalar integer inputs.
   * Range limited to values that fit in unsigned 64-bit integer (up to ~2^64-1).
   * Default: false
   */
  exact?: boolean;
}

/**
 * Permutations (k-permutations of N).
 *
 * Computes P(N, k) = N! / (N-k)!, the number of ways to arrange k items
 * from N items where order matters.
 *
 * @param N - Number of items (integer or array)
 * @param k - Number of items to arrange (integer or array)
 * @param options - Configuration options
 * @returns Number of k-permutations
 *
 * @example
 * ```typescript
 * // Non-exact mode (floating point)
 * await perm(10, 3);  // 720
 * await perm([10, 10], [3, 4]);  // [720, 5040]
 *
 * // Exact mode (integer arithmetic)
 * await perm(10, 3, { exact: true });  // 720
 * ```
 *
 * Notes:
 * - Returns 0 if k > N, k < 0, or N < 0
 * - Returns 1 if k = 0
 * - In exact mode, only scalar integer inputs are supported
 * - Exact mode range limited to unsigned 64-bit integers
 * - Non-exact mode uses Pochhammer symbol (rising factorial) via lgamma
 */
export async function perm(
  N: number | NDArray,
  k: number | NDArray,
  options: PermOptions = {}
): Promise<number | NDArray> {
  const { exact = false } = options;

  // Exact mode
  if (exact) {
    validateScalarInputs(N, k, 'perm');
    const nVal = N as number;
    const kVal = k as number;
    validateIntegerInputs(nVal, kVal, 'perm');

    const wasm = getWasmModule();
    const result = wasm._wasm_perm_exact(nVal, kVal);

    // Check for overflow
    checkOverflow(result, 'perm');

    return result;
  }

  // Non-exact mode using Pochhammer symbol
  // P(N, k) = poch(N-k+1, k) = (N-k+1) * (N-k+2) * ... * N
  if (isScalar(N) && isScalar(k)) {
    const wasm = getWasmModule();
    return wasm._wasm_poch((N as number) - (k as number) + 1, k as number);
  }

  // Array mode: apply element-wise
  const NArray = N instanceof NDArray ? N : await NDArray.array([N]);
  const kArray = k instanceof NDArray ? k : await NDArray.array([k]);

  return elementWise(NArray, kArray, (n, k_val) => {
    const wasm = getWasmModule();
    return wasm._wasm_poch(n - k_val + 1, k_val);
  });
}
