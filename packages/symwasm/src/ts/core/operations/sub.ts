/**
 * Subtraction operation.
 * @module core/operations/sub
 */

import { getWasmModule } from '../../wasm-loader.js';
import { createBasic, checkException } from '../../wasm-memory.js';
import { Expr } from '../expr.js';
import { exprFromWasm } from '../expr-factory.js';

/**
 * Subtract two expressions: a - b
 * @param a First operand
 * @param b Second operand
 * @returns The difference (may be simplified by SymEngine)
 */
export function sub(a: Expr, b: Expr): Expr {
  const wasm = getWasmModule();
  const obj = createBasic();
  const code = wasm._basic_sub(obj.getPtr(), a.getWasmPtr(), b.getWasmPtr());
  checkException(code);
  return exprFromWasm(obj);
}
