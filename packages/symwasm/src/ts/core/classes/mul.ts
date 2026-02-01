/**
 * Mul class for symbolic multiplication.
 * @module core/classes/mul
 */

import { getWasmModule } from '../../wasm-loader.js';
import { SymEngineObject, SymEngineVec, checkException } from '../../wasm-memory.js';
import { Expr } from '../expr.js';

// Lazy import to avoid circular dependency
let _exprFromWasm: ((obj: SymEngineObject) => Expr) | null = null;

function getExprFromWasm() {
  if (!_exprFromWasm) {
    _exprFromWasm = require('../expr-factory.js').exprFromWasm;
  }
  return _exprFromWasm;
}

/**
 * Symbolic multiplication. Mirrors sympy.Mul.
 * Note: Use the mul() function to create Mul expressions.
 */
export class Mul extends Expr {
  private _args: Expr[] | null = null;

  /** Get the factors of this multiplication (lazily extracted from WASM) */
  get args(): Expr[] {
    if (!this._args) {
      this._args = this._extractArgs();
    }
    return this._args;
  }

  private _extractArgs(): Expr[] {
    if (!this._obj) return [];
    const wasm = getWasmModule();
    const vec = new SymEngineVec();
    const exprFromWasm = getExprFromWasm();
    try {
      const code = wasm._basic_get_args(this._obj.getPtr(), vec.getPtr());
      checkException(code);
      const result: Expr[] = [];
      const count = vec.size();
      for (let i = 0; i < count; i++) {
        result.push(exprFromWasm(vec.get(i)));
      }
      return result;
    } finally {
      vec.free();
    }
  }

  protected _fallbackString(): string {
    return this.args.map((a) => a.toString()).join('*');
  }

  /**
   * @internal Create a Mul from an existing WASM object.
   */
  static _fromWasm(obj: SymEngineObject): Mul {
    const mul = Object.create(Mul.prototype) as Mul;
    (mul as any)._obj = obj;
    (mul as any)._args = null;
    return mul;
  }
}
