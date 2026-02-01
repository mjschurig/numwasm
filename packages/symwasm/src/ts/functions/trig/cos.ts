/**
 * Cosine function.
 * @module functions/trig/cos
 */

import { getWasmModule } from '../../wasm-loader.js';
import { createBasic, checkException } from '../../wasm-memory.js';
import type { Expr } from '../../core/expr.js';
import { exprFromWasm } from '../../core/expr-factory.js';

/**
 * Cosine: cos(x)
 * @param x The angle in radians
 * @returns The cosine of x
 */
export function cos(x: Expr): Expr {
  const wasm = getWasmModule();
  const obj = createBasic();
  try {
    const code = wasm._basic_cos(obj.getPtr(), x.getWasmPtr());
    checkException(code);
    return exprFromWasm(obj);
  } catch (e) {
    obj.free();
    throw e;
  }
}
