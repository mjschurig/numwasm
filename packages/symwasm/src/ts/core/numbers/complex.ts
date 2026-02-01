/**
 * Complex number type.
 * @module core/numbers/complex
 */

import { getWasmModule } from '../../wasm-loader.js';
import { SymEngineObject, createBasic, checkException } from '../../wasm-memory.js';
import { Expr } from '../expr.js';
import { Integer } from './integer.js';
import { Rational } from './rational.js';

/**
 * Exact complex number re + im*I.
 * Uses SymEngine's Complex class which stores rational real and imaginary parts.
 */
export class Complex extends Expr {
  readonly re: Integer | Rational;
  readonly im: Integer | Rational;

  /**
   * Create an exact complex number re + im*I.
   * @param re Real part (number, Integer, or Rational)
   * @param im Imaginary part (number, Integer, or Rational)
   */
  constructor(re: number | Integer | Rational, im: number | Integer | Rational) {
    super();

    // Convert numbers to Integer
    const reExpr = typeof re === 'number' ? new Integer(re) : re;
    const imExpr = typeof im === 'number' ? new Integer(im) : im;

    this.re = reExpr;
    this.im = imExpr;

    // Create WASM-backed complex
    const wasm = getWasmModule();
    const obj = createBasic();
    const code = wasm._complex_set(
      obj.getPtr(),
      (reExpr as any)._obj.getPtr(),
      (imExpr as any)._obj.getPtr()
    );
    checkException(code);
    this._obj = obj;
  }

  protected _fallbackString(): string {
    return `${this.re} + ${this.im}*I`;
  }

  /**
   * @internal Create a Complex from an existing WASM object.
   * Note: re/im fields are placeholders - full parsing not implemented.
   */
  static _fromWasm(obj: SymEngineObject): Complex {
    const c = Object.create(Complex.prototype) as Complex;
    (c as any)._obj = obj;
    // Placeholder values - full parsing of "a + b*I" format not implemented
    (c as any).re = new Integer(0);
    (c as any).im = new Integer(0);
    return c;
  }
}
