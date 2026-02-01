/**
 * Multiplication operation.
 * @module core/operations/mul
 */

import { getWasmModule } from '../../wasm-loader.js';
import { createBasic, checkException } from '../../wasm-memory.js';
import { Expr } from '../expr.js';
import { exprFromWasm } from '../expr-factory.js';

/**
 * Multiply two expressions: a * b
 * @param a First operand
 * @param b Second operand
 * @returns The product (may be simplified by SymEngine)
 */
export function mul(a: Expr, b: Expr): Expr {
  const wasm = getWasmModule();
  const obj = createBasic();
  const code = wasm._basic_mul(obj.getPtr(), a.getWasmPtr(), b.getWasmPtr());
  checkException(code);
  return exprFromWasm(obj);
}
