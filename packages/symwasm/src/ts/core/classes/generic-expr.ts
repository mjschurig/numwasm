/**
 * GenericExpr class for unrecognized expression types.
 * @module core/classes/generic-expr
 * @internal
 */

import { SymEngineObject } from '../../wasm-memory.js';
import { Expr } from '../expr.js';

/**
 * Generic expression wrapper for types not yet fully implemented.
 * @internal
 */
export class GenericExpr extends Expr {
  protected _fallbackString(): string {
    return this._obj ? this._obj.toString() : '<unknown>';
  }

  static _fromWasm(obj: SymEngineObject): GenericExpr {
    const e = Object.create(GenericExpr.prototype) as GenericExpr;
    (e as any)._obj = obj;
    return e;
  }
}
