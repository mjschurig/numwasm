/**
 * Division operation.
 * @module core/operations/div
 */

import { getWasmModule } from '../../wasm-loader.js';
import { createBasic, checkException } from '../../wasm-memory.js';
import { Expr } from '../expr.js';
import { exprFromWasm } from '../expr-factory.js';

/**
 * Divide two expressions: a / b
 * @param a Numerator
 * @param b Denominator
 * @returns The quotient (may be simplified by SymEngine)
 */
export function div(a: Expr, b: Expr): Expr {
  const wasm = getWasmModule();
  const obj = createBasic();
  const code = wasm._basic_div(obj.getPtr(), a.getWasmPtr(), b.getWasmPtr());
  checkException(code);
  return exprFromWasm(obj);
}
