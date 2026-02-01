/**
 * ImaginaryUnit class for the imaginary unit i.
 * @module core/classes/imaginary-unit
 */

import { getWasmModule } from '../../wasm-loader.js';
import { SymEngineObject, createBasic } from '../../wasm-memory.js';
import { Expr } from '../expr.js';

/**
 * Imaginary unit (i = sqrt(-1)).
 * Special constant backed by SymEngine's I constant.
 */
export class ImaginaryUnit extends Expr {
  private constructor() {
    super();
  }

  protected _fallbackString(): string {
    return 'I';
  }

  /** Create the imaginary unit I */
  static create(): ImaginaryUnit {
    const i = new ImaginaryUnit();
    const obj = createBasic();
    getWasmModule()._basic_const_I(obj.getPtr());
    (i as any)._obj = obj;
    return i;
  }

  /** @internal Create from WASM object */
  static _fromWasm(obj: SymEngineObject): ImaginaryUnit {
    const i = Object.create(ImaginaryUnit.prototype) as ImaginaryUnit;
    (i as any)._obj = obj;
    return i;
  }
}
