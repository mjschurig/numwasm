/**
 * WASM module type definitions for sciwasm
 */

export interface SciWasmModule {
  // Nelder-Mead optimization
  _nelder_mead_minimize(
    n: number,
    x0Ptr: number,
    xOutPtr: number,
    simOutPtr: number,
    fsimOutPtr: number,
    resultPtr: number,
    xatol: number,
    fatol: number,
    maxiter: number,
    maxfev: number,
    adaptive: number,
    hasBounds: number,
    lowerPtr: number,
    upperPtr: number,
    initialSimplexPtr: number,
    funcPtr: number
  ): number;

  // QUADPACK adaptive quadrature
  _wasm_dqagse(
    fcn: number, a: number, b: number,
    epsabs: number, epsrel: number, limit: number,
    result: number, abserr: number, neval: number, ier: number,
    alist: number, blist: number, rlist: number, elist: number,
    iord: number, last: number,
  ): void;

  _wasm_dqagie(
    fcn: number, bound: number, inf: number,
    epsabs: number, epsrel: number, limit: number,
    result: number, abserr: number, neval: number, ier: number,
    alist: number, blist: number, rlist: number, elist: number,
    iord: number, last: number,
  ): void;

  // Memory management
  _malloc(size: number): number;
  _free(ptr: number): void;

  // Emscripten runtime methods
  getValue(ptr: number, type: string): number;
  setValue(ptr: number, value: number, type: string): void;

  // Heap views
  HEAPF64: Float64Array;
  HEAPF32: Float32Array;
  HEAP32: Int32Array;
  HEAP8: Int8Array;
  HEAPU8: Uint8Array;
  HEAPU32: Uint32Array;

  // Function table management (for JS→C callbacks)
  addFunction(func: Function, signature: string): number;
  removeFunction(ptr: number): void;
}

export interface WasmModuleOptions {
  locateFile?: (path: string, scriptDirectory: string) => string;
}

export type WasmModuleFactory = (
  options?: WasmModuleOptions
) => Promise<SciWasmModule>;
