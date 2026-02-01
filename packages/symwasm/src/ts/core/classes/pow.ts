/**
 * Pow class for symbolic exponentiation.
 * @module core/classes/pow
 */

import { getWasmModule } from '../../wasm-loader.js';
import { SymEngineObject, SymEngineVec, checkException } from '../../wasm-memory.js';
import { Expr } from '../expr.js';
import { Integer } from '../numbers/integer.js';

// Lazy import to avoid circular dependency
let _exprFromWasm: ((obj: SymEngineObject) => Expr) | null = null;

function getExprFromWasm() {
  if (!_exprFromWasm) {
    _exprFromWasm = require('../expr-factory.js').exprFromWasm;
  }
  return _exprFromWasm;
}

/**
 * Symbolic exponentiation. Mirrors sympy.Pow.
 * Note: Use the pow() function to create Pow expressions.
 */
export class Pow extends Expr {
  private _base: Expr | null = null;
  private _exp: Expr | null = null;

  /** Get the base of this power */
  get base(): Expr {
    if (!this._base) {
      this._extractBaseExp();
    }
    return this._base!;
  }

  /** Get the exponent of this power */
  get exp(): Expr {
    if (!this._exp) {
      this._extractBaseExp();
    }
    return this._exp!;
  }

  private _extractBaseExp(): void {
    if (!this._obj) {
      this._base = new Integer(0);
      this._exp = new Integer(0);
      return;
    }
    const wasm = getWasmModule();
    const vec = new SymEngineVec();
    const exprFromWasm = getExprFromWasm();
    try {
      const code = wasm._basic_get_args(this._obj.getPtr(), vec.getPtr());
      checkException(code);
      if (vec.size() >= 2) {
        this._base = exprFromWasm(vec.get(0));
        this._exp = exprFromWasm(vec.get(1));
      } else {
        this._base = new Integer(0);
        this._exp = new Integer(0);
      }
    } finally {
      vec.free();
    }
  }

  protected _fallbackString(): string {
    return this.base + '**' + this.exp;
  }

  /**
   * @internal Create a Pow from an existing WASM object.
   */
  static _fromWasm(obj: SymEngineObject): Pow {
    const p = Object.create(Pow.prototype) as Pow;
    (p as any)._obj = obj;
    (p as any)._base = null;
    (p as any)._exp = null;
    return p;
  }
}
