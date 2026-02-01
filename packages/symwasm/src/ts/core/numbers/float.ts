/**
 * Float number type.
 * @module core/numbers/float
 */

import { getWasmModule } from '../../wasm-loader.js';
import { SymEngineObject, createBasic, checkException } from '../../wasm-memory.js';
import { Expr } from '../expr.js';

/**
 * Floating-point number. Mirrors sympy.Float.
 * Note: SymEngine uses double precision internally (no arbitrary precision without MPFR).
 */
export class Float extends Expr {
  readonly value: number;

  /**
   * Create a floating-point number.
   * @param value The numeric value
   * @param _precision Ignored (SymEngine uses double precision without MPFR)
   */
  constructor(value: number, _precision?: number) {
    super();
    this.value = value;

    // Create WASM-backed real double
    const wasm = getWasmModule();
    const obj = createBasic();
    const code = wasm._real_double_set_d(obj.getPtr(), this.value);
    checkException(code);
    this._obj = obj;
  }

  protected _fallbackString(): string {
    return String(this.value);
  }

  /**
   * @internal Create a Float from an existing WASM object.
   */
  static _fromWasm(obj: SymEngineObject): Float {
    const f = Object.create(Float.prototype) as Float;
    (f as any)._obj = obj;
    (f as any).value = parseFloat(obj.toString());
    return f;
  }
}
