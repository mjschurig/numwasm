/**
 * Negation operation.
 * @module core/operations/neg
 */

import { getWasmModule } from '../../wasm-loader.js';
import { createBasic, checkException } from '../../wasm-memory.js';
import { Expr } from '../expr.js';
import { exprFromWasm } from '../expr-factory.js';

/**
 * Negate an expression: -a
 * @param a The operand
 * @returns The negation
 */
export function neg(a: Expr): Expr {
  const wasm = getWasmModule();
  const obj = createBasic();
  const code = wasm._basic_neg(obj.getPtr(), a.getWasmPtr());
  checkException(code);
  return exprFromWasm(obj);
}
