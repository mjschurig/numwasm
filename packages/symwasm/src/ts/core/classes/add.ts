/**
 * Add class for symbolic addition.
 * @module core/classes/add
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
 * Symbolic addition. Mirrors sympy.Add.
 * Note: Use the add() function to create Add expressions.
 */
export class Add extends Expr {
  private _args: Expr[] | null = null;

  /** Get the terms of this addition (lazily extracted from WASM) */
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
    return this.args.map((a) => a.toString()).join(' + ');
  }

  /**
   * @internal Create an Add from an existing WASM object.
   */
  static _fromWasm(obj: SymEngineObject): Add {
    const add = Object.create(Add.prototype) as Add;
    (add as any)._obj = obj;
    (add as any)._args = null;
    return add;
  }
}
