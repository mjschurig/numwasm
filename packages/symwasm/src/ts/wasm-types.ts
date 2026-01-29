/**
 * Type definitions for the SymWASM module
 * These interfaces define the WASM module's API and memory access
 */

export interface EmscriptenModule {
  // Memory views
  HEAPF64: Float64Array;
  HEAP32: Int32Array;
  HEAPU8: Uint8Array;
  HEAPF32: Float32Array;
  HEAP16: Int16Array;
  HEAPU16: Uint16Array;
  HEAP8: Int8Array;
  HEAPU32: Uint32Array;

  // Memory management
  _malloc(size: number): number;
  _free(ptr: number): void;

  // Utility functions
  ccall(
    name: string,
    returnType: string | null,
    argTypes: string[],
    args: any[]
  ): any;
  cwrap(
    name: string,
    returnType: string | null,
    argTypes: string[]
  ): (...args: any[]) => any;
  UTF8ToString(ptr: number, maxBytesToRead?: number): string;
  stringToUTF8(
    str: string,
    outPtr: number,
    maxBytesToWrite: number
  ): void;
  lengthBytesUTF8(str: string): number;
  getValue(ptr: number, type: string): number;
  setValue(ptr: number, value: number, type: string): void;
}

/**
 * SymWASM module interface
 * Extends Emscripten module with SymEngine C API functions
 */
export interface SymwasmModule extends EmscriptenModule {
  // === Memory Management ===
  _basic_new_stack(): number; // For stack-allocated basic (requires pre-allocated memory)
  _basic_free_stack(ptr: number): void;
  _basic_new_heap(): number; // For heap-allocated basic (allocates and returns pointer)
  _basic_free_heap(ptr: number): void;

  // === Symbol Creation ===
  _symbol_set(ptr: number, namePtr: number): number; // namePtr is a C string pointer

  // === Number Creation ===
  _integer_set_si(ptr: number, value: number): number;
  _integer_set_ui(ptr: number, value: number): number;
  _rational_set(ptr: number, num: number, den: number): number;
  _rational_set_si(ptr: number, num: number, den: number): number;
  _rational_set_ui(ptr: number, num: number, den: number): number;
  _real_double_set_d(ptr: number, value: number): number;

  // === Basic Arithmetic ===
  _basic_add(result: number, a: number, b: number): number;
  _basic_sub(result: number, a: number, b: number): number;
  _basic_mul(result: number, a: number, b: number): number;
  _basic_div(result: number, a: number, b: number): number;
  _basic_pow(result: number, base: number, exp: number): number;

  // === Constants ===
  _basic_const_zero(ptr: number): number;
  _basic_const_one(ptr: number): number;
  _basic_const_minus_one(ptr: number): number;
  _basic_const_I(ptr: number): number;
  _basic_const_pi(ptr: number): number;
  _basic_const_E(ptr: number): number;

  // === Type Information ===
  _basic_get_type(ptr: number): number;

  // === Comparison ===
  _basic_eq(a: number, b: number): number;
  _basic_hash(ptr: number): number;

  // === String Conversion ===
  _basic_str(ptr: number): number; // Returns char* (must be freed)

  // === Free Symbols Extraction ===
  _basic_free_symbols(self: number, symbols: number): number; // Returns exception code

  // === Set Container Operations (CSetBasic) ===
  _setbasic_new(): number; // Returns CSetBasic*
  _setbasic_free(self: number): void;
  _setbasic_get(self: number, n: number, result: number): void; // Gets nth element into result
  _setbasic_size(self: number): number; // Returns size_t
}

/**
 * SymWASM module factory function type
 */
export interface SymwasmModuleFactory {
  (options?: Partial<EmscriptenModule>): Promise<SymwasmModule>;
}

/**
 * SymEngine type IDs (from type_codes.inc)
 * These correspond to the TypeID enum in SymEngine
 * Note: Values depend on compile-time configuration (MPFR, MPC, PIRANHA, FLINT)
 * These values are empirically determined from the WASM build
 */
export enum SymEngineTypeID {
  SYMENGINE_INTEGER = 0,
  SYMENGINE_RATIONAL = 1,
  SYMENGINE_COMPLEX = 2,
  SYMENGINE_COMPLEX_DOUBLE = 3,
  SYMENGINE_REAL_MPFR = 4,
  SYMENGINE_COMPLEX_MPC = 5,
  SYMENGINE_REAL_DOUBLE = 6,
  SYMENGINE_INFTY = 7,
  SYMENGINE_NOT_A_NUMBER = 8,
  SYMENGINE_URATPSERIESPIRANHA = 9,
  SYMENGINE_UPSERIESPIRANHA = 10,
  SYMENGINE_URATPSERIESFLINT = 11,
  SYMENGINE_NUMBER_WRAPPER = 12,
  // NUMBER_WRAPPER marks the end of Number subclasses
  SYMENGINE_SYMBOL = 13,
  SYMENGINE_DUMMY = 14,
  SYMENGINE_MUL = 15,
  SYMENGINE_ADD = 16,
  SYMENGINE_POW = 17,
  // ... more types follow
  SYMENGINE_CONSTANT = 31,
  // Add more as needed
}

/**
 * SymEngine exception codes (from symengine_exception.h)
 */
export enum SymEngineException {
  NO_EXCEPTION = 0,
  RUNTIME_ERROR = 1,
  DIV_BY_ZERO = 2,
  NOT_IMPLEMENTED = 3,
  DOMAIN_ERROR = 4,
  PARSE_ERROR = 5,
}
