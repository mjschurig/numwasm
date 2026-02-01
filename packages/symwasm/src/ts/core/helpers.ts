/**
 * Internal helper functions for creating WASM function wrappers.
 * @module core/helpers
 * @internal
 */

import { getWasmModule } from '../wasm-loader.js';
import { SymEngineObject, createBasic, checkException } from '../wasm-memory.js';
import type { Expr } from './expr.js';

/**
 * Internal helper to create a 1-argument function wrapper.
 * @internal
 */
export function makeOneArgFunc(
  wasmFnName: string
): (x: Expr) => Expr {
  // Lazy import to avoid circular dependency
  let exprFromWasm: ((obj: SymEngineObject) => Expr) | null = null;

  return (x: Expr): Expr => {
    const wasm = getWasmModule();
    const obj = createBasic();
    try {
      const fn = (wasm as any)[wasmFnName] as (result: number, arg: number) => number;
      const code = fn.call(wasm, obj.getPtr(), x.getWasmPtr());
      checkException(code);

      // Lazy load exprFromWasm
      if (!exprFromWasm) {
        exprFromWasm = require('./expr-factory.js').exprFromWasm;
      }
      return exprFromWasm(obj);
    } catch (e) {
      obj.free();
      throw e;
    }
  };
}

/**
 * Internal helper to create a 2-argument function wrapper.
 * @internal
 */
export function makeTwoArgFunc(
  wasmFnName: string
): (a: Expr, b: Expr) => Expr {
  // Lazy import to avoid circular dependency
  let exprFromWasm: ((obj: SymEngineObject) => Expr) | null = null;

  return (a: Expr, b: Expr): Expr => {
    const wasm = getWasmModule();
    const obj = createBasic();
    try {
      const fn = (wasm as any)[wasmFnName] as (
        result: number,
        arg1: number,
        arg2: number
      ) => number;
      const code = fn.call(wasm, obj.getPtr(), a.getWasmPtr(), b.getWasmPtr());
      checkException(code);

      // Lazy load exprFromWasm
      if (!exprFromWasm) {
        exprFromWasm = require('./expr-factory.js').exprFromWasm;
      }
      return exprFromWasm(obj);
    } catch (e) {
      obj.free();
      throw e;
    }
  };
}

/**
 * Parse a complex number string representation from SymEngine.
 * Handles formats like: "1.0 + 2.0*I", "1.5 - 2.3*I", "2.3*I", "-2.3*I", "1.5"
 * @internal
 */
export function parseComplexString(str: string): { real: number; imag: number } {
  const s = str.trim();

  // Pure imaginary: just "I" or "-I"
  if (s === 'I') {
    return { real: 0, imag: 1 };
  }
  if (s === '-I') {
    return { real: 0, imag: -1 };
  }

  // Check for *I pattern (pure imaginary): "2.3*I" or "-2.3*I"
  const pureImagMatch = s.match(/^(-?[\d.eE+-]+)\*I$/);
  if (pureImagMatch) {
    return { real: 0, imag: parseFloat(pureImagMatch[1]) };
  }

  // Complex with both parts with spaces: "real + imag*I" or "real - imag*I"
  const complexWithSpacesMatch = s.match(/^(-?[\d.eE+-]+)\s*([+-])\s*([\d.eE+-]+)\*I$/);
  if (complexWithSpacesMatch) {
    const real = parseFloat(complexWithSpacesMatch[1]);
    const sign = complexWithSpacesMatch[2] === '-' ? -1 : 1;
    const imag = sign * parseFloat(complexWithSpacesMatch[3]);
    return { real, imag };
  }

  // Handle case where imag coefficient is 1: "real + I" or "real - I"
  const complexOneMatch = s.match(/^(-?[\d.eE+-]+)\s*([+-])\s*I$/);
  if (complexOneMatch) {
    const real = parseFloat(complexOneMatch[1]);
    const imag = complexOneMatch[2] === '-' ? -1 : 1;
    return { real, imag };
  }

  // Pure real number
  return { real: parseFloat(s), imag: 0 };
}
