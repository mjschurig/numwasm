/**
 * Factory function for creating Expr instances from WASM objects.
 * @module core/expr-factory
 */

import { SymEngineObject } from '../wasm-memory.js';
import { SymEngineTypeID } from '../wasm-types.js';
import { Expr } from './expr.js';
import { Integer } from './numbers/integer.js';
import { Rational } from './numbers/rational.js';
import { Float } from './numbers/float.js';
import { Complex } from './numbers/complex.js';
import { Symbol } from './classes/symbol.js';
import { Add } from './classes/add.js';
import { Mul } from './classes/mul.js';
import { Pow } from './classes/pow.js';
import { Constant } from './classes/constant.js';
import { Infinity_ } from './classes/infinity.js';
import { NaN_ } from './classes/nan.js';
import { GenericExpr } from './classes/generic-expr.js';

/**
 * Create an Expr subclass instance from a WASM object based on its type.
 * Used internally and by the matrices module to convert WASM basic pointers to Expr.
 */
export function exprFromWasm(obj: SymEngineObject): Expr {
  const typeId = obj.getType();
  switch (typeId) {
    case SymEngineTypeID.SYMENGINE_INTEGER:
      return Integer._fromWasm(obj);
    case SymEngineTypeID.SYMENGINE_RATIONAL:
      return Rational._fromWasm(obj);
    case SymEngineTypeID.SYMENGINE_REAL_DOUBLE:
      return Float._fromWasm(obj);
    case SymEngineTypeID.SYMENGINE_COMPLEX_DOUBLE:
      return GenericExpr._fromWasm(obj); // ComplexDouble - handled via evalfComplex
    case SymEngineTypeID.SYMENGINE_COMPLEX:
      return Complex._fromWasm(obj);
    case SymEngineTypeID.SYMENGINE_SYMBOL:
      return Symbol._fromWasm(obj);
    case SymEngineTypeID.SYMENGINE_ADD:
      return Add._fromWasm(obj);
    case SymEngineTypeID.SYMENGINE_MUL:
      return Mul._fromWasm(obj);
    case SymEngineTypeID.SYMENGINE_POW:
      return Pow._fromWasm(obj);
    case SymEngineTypeID.SYMENGINE_CONSTANT:
      return Constant._fromWasm(obj);
    case SymEngineTypeID.SYMENGINE_INFTY:
      return Infinity_._fromWasm(obj);
    case SymEngineTypeID.SYMENGINE_NOT_A_NUMBER:
      return NaN_._fromWasm(obj);
    default:
      return GenericExpr._fromWasm(obj);
  }
}
