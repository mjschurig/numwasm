/**
 * Integer number type.
 * @module core/numbers/integer
 */

import { getWasmModule } from '../../wasm-loader.js';
import { SymEngineObject, createBasic, checkException } from '../../wasm-memory.js';
import { Expr } from '../expr.js';

/**
 * Exact integer value. Mirrors sympy.Integer.
 */
export class Integer extends Expr {
  readonly value: number;

  /**
   * Create an exact integer value.
   * @param value The integer value (must fit in a signed 32-bit integer)
   */
  constructor(value: number) {
    super();
    this.value = Math.trunc(value); // Ensure it's an integer

    // Create WASM-backed integer
    const wasm = getWasmModule();
    const obj = createBasic();
    const code = wasm._integer_set_si(obj.getPtr(), this.value);
    checkException(code);
    this._obj = obj;
  }

  protected _fallbackString(): string {
    return String(this.value);
  }

  /**
   * @internal Create an Integer from an existing WASM object.
   */
  static _fromWasm(obj: SymEngineObject): Integer {
    const int = Object.create(Integer.prototype) as Integer;
    (int as any)._obj = obj;
    (int as any).value = parseInt(obj.toString(), 10);
    return int;
  }
}
