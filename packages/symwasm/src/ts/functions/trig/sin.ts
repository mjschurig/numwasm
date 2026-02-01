/**
 * Sine function.
 * @module functions/trig/sin
 */

import { getWasmModule } from '../../wasm-loader.js';
import { createBasic, checkException } from '../../wasm-memory.js';
import type { Expr } from '../../core/expr.js';
import { exprFromWasm } from '../../core/expr-factory.js';

/**
 * Sine: sin(x)
 * @param x The angle in radians
 * @returns The sine of x
 */
export function sin(x: Expr): Expr {
  const wasm = getWasmModule();
  const obj = createBasic();
  try {
    const code = wasm._basic_sin(obj.getPtr(), x.getWasmPtr());
    checkException(code);
    return exprFromWasm(obj);
  } catch (e) {
    obj.free();
    throw e;
  }
}
