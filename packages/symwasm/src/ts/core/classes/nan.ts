/**
 * NaN class for not-a-number values.
 * @module core/classes/nan
 */

import { getWasmModule } from '../../wasm-loader.js';
import { SymEngineObject, createBasic } from '../../wasm-memory.js';
import { Expr } from '../expr.js';

/**
 * Not a Number (undefined/indeterminate result).
 * WASM-backed version.
 */
export class NaN_ extends Expr {
  private constructor() {
    super();
  }

  protected _fallbackString(): string {
    return 'nan';
  }

  /** Create NaN */
  static create(): NaN_ {
    const n = new NaN_();
    const obj = createBasic();
    getWasmModule()._basic_const_nan(obj.getPtr());
    (n as any)._obj = obj;
    return n;
  }

  /** @internal Create from WASM object */
  static _fromWasm(obj: SymEngineObject): NaN_ {
    const n = Object.create(NaN_.prototype) as NaN_;
    (n as any)._obj = obj;
    return n;
  }
}
