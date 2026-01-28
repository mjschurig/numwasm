/**
 * Binomial coefficient: C(N, k) = N! / (k! * (N-k)!)
 *
 * Computes the number of ways to choose k items from N items without replacement
 * and without order.
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

export interface CombOptions {
  /**
   * If true, compute exactly using integer arithmetic.
   * Only supports scalar integer inputs.
   * Range limited to values that fit in unsigned 64-bit integer (up to ~2^64-1).
   * Default: false
   */
  exact?: boolean;

  /**
   * If true, compute combinations with repetition (multicombinations).
   * Formula: C(N+k-1, k)
   * Default: false
   */
  repetition?: boolean;
}

/**
 * Binomial coefficient (N choose k).
 *
 * @param N - Number of items (integer or array)
 * @param k - Number of items to choose (integer or array)
 * @param options - Configuration options
 * @returns Binomial coefficient C(N, k)
 *
 * @example
 * ```typescript
 * // Non-exact mode (floating point)
 * await comb(10, 3);  // 120
 * await comb([10, 10], [3, 4]);  // [120, 210]
 *
 * // Exact mode (integer arithmetic)
 * await comb(10, 3, { exact: true });  // 120
 *
 * // Combinations with repetition
 * await comb(10, 3, { exact: true, repetition: true });  // 220
 * ```
 *
 * Notes:
 * - Returns 0 if k > N (when repetition=false), k < 0, or N < 0
 * - In exact mode, only scalar integer inputs are supported
 * - Exact mode range limited to unsigned 64-bit integers
 * - Non-exact mode uses lgamma for computation
 */
export async function comb(
  N: number | NDArray,
  k: number | NDArray,
  options: CombOptions = {}
): Promise<number | NDArray> {
  const { exact = false, repetition = false } = options;

  // Handle repetition mode: transform to comb(N+k-1, k)
  if (repetition) {
    if (exact) {
      // Validate scalar inputs
      validateScalarInputs(N, k, 'comb');
      const nVal = N as number;
      const kVal = k as number;
      validateIntegerInputs(nVal, kVal, 'comb');

      // Special case: C(n, 0, repetition=True) = 1 for n >= 0
      if (kVal === 0 && nVal >= 0) {
        return 1;
      }

      // Transform: C(n, k, repetition=True) = C(n+k-1, k)
      return comb(nVal + kVal - 1, kVal, { exact: true });
    } else {
      // Non-exact array mode
      if (isScalar(N) && isScalar(k)) {
        const nVal = N as number;
        const kVal = k as number;

        // Special case for k=0
        if (kVal === 0 && nVal >= 0) {
          return 1.0;
        }

        const wasm = getWasmModule();
        const cond = kVal === 0 && nVal >= 0;
        if (cond) {
          return 1.0;
        }

        // Use non-exact binom on transformed inputs
        return wasm._wasm_binom(nVal + kVal - 1, kVal);
      } else {
        // Array case: apply element-wise
        const NArray = N instanceof NDArray ? N : await NDArray.array([N]);
        const kArray = k instanceof NDArray ? k : await NDArray.array([k]);

        return elementWise(NArray, kArray, (n, k_val) => {
          if (k_val === 0 && n >= 0) return 1.0;
          const wasm = getWasmModule();
          return wasm._wasm_binom(n + k_val - 1, k_val);
        });
      }
    }
  }

  // Exact mode
  if (exact) {
    validateScalarInputs(N, k, 'comb');
    const nVal = N as number;
    const kVal = k as number;
    validateIntegerInputs(nVal, kVal, 'comb');

    const wasm = getWasmModule();
    const result = wasm._wasm_binom_exact(nVal, kVal);

    // Check for overflow
    checkOverflow(result, 'comb');

    return result;
  }

  // Non-exact mode
  if (isScalar(N) && isScalar(k)) {
    const wasm = getWasmModule();
    return wasm._wasm_binom(N as number, k as number);
  }

  // Array mode: apply element-wise
  const NArray = N instanceof NDArray ? N : await NDArray.array([N]);
  const kArray = k instanceof NDArray ? k : await NDArray.array([k]);

  return elementWise(NArray, kArray, (n, k_val) => {
    const wasm = getWasmModule();
    return wasm._wasm_binom(n, k_val);
  });
}
