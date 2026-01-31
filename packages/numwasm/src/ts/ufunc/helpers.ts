/**
 * Ufunc helper functions for NumJS-WASM
 *
 * Internal helpers for applying unary and binary ufuncs.
 */

import { NDArray } from "../_core/NDArray.js";
import { getWasmModule } from "../wasm-loader.js";

/**
 * Helper to apply a unary ufunc.
 * @internal
 */
export function applyUnary(
  input: NDArray,
  wasmFunc: (ptr: number) => number,
): NDArray {
  const module = getWasmModule();
  const resultPtr = wasmFunc(input._wasmPtr);
  if (resultPtr === 0) {
    throw new Error("Ufunc operation failed");
  }
  return NDArray._fromPtr(resultPtr, module);
}

/**
 * Helper to apply a binary ufunc.
 * @internal
 */
export function applyBinary(
  x1: NDArray,
  x2: NDArray,
  wasmFunc: (ptr1: number, ptr2: number) => number,
): NDArray {
  const module = getWasmModule();
  const resultPtr = wasmFunc(x1._wasmPtr, x2._wasmPtr);
  if (resultPtr === 0) {
    throw new Error(
      "Ufunc operation failed (possibly due to incompatible shapes for broadcasting)",
    );
  }
  return NDArray._fromPtr(resultPtr, module);
}
