/**
 * Power operation.
 * @module core/operations/pow
 */

import { getWasmModule } from '../../wasm-loader.js';
import { createBasic, checkException } from '../../wasm-memory.js';
import { Expr } from '../expr.js';
import { exprFromWasm } from '../expr-factory.js';

/**
 * Raise to a power: base ** exp
 * @param base The base
 * @param exp The exponent
 * @returns The power (may be simplified by SymEngine)
 */
export function pow(base: Expr, exp: Expr): Expr {
  const wasm = getWasmModule();
  const obj = createBasic();
  const code = wasm._basic_pow(obj.getPtr(), base.getWasmPtr(), exp.getWasmPtr());
  checkException(code);
  return exprFromWasm(obj);
}
