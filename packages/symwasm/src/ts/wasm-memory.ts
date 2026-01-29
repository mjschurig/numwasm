/**
 * WASM Memory Management Utilities
 * Provides wrapper classes and helpers for managing SymEngine WASM objects
 */

import { getWasmModule } from './wasm-loader.js';

/**
 * SymEngine exception codes
 * Maps to symengine_exceptions_t enum in cwrapper.h
 */
export enum SymEngineException {
  NO_EXCEPTION = 0,
  RUNTIME_ERROR = 1,
  DIV_BY_ZERO = 2,
  NOT_IMPLEMENTED = 3,
  DOMAIN_ERROR = 4,
  PARSE_ERROR = 5,
}

/**
 * Wrapper class for SymEngine Basic objects
 * Manages lifecycle and memory cleanup
 */
export class SymEngineObject {
  private ptr: number;

  constructor(ptr: number) {
    this.ptr = ptr;
  }

  /**
   * Get the WASM pointer to this object
   */
  getPtr(): number {
    return this.ptr;
  }

  /**
   * Check if this object has been freed
   */
  isValid(): boolean {
    return this.ptr !== 0;
  }

  /**
   * Free the WASM memory for this object
   * WARNING: After calling this, the object is invalid
   */
  free(): void {
    if (this.ptr) {
      const wasm = getWasmModule();
      wasm._basic_free_stack(this.ptr);
      this.ptr = 0;
    }
  }

  /**
   * Convert to string representation
   */
  toString(): string {
    if (!this.isValid()) {
      throw new Error('Cannot convert freed SymEngineObject to string');
    }

    const wasm = getWasmModule();
    const strPtr = wasm._basic_str(this.ptr);
    const str = wasm.UTF8ToString(strPtr);
    wasm._free(strPtr);
    return str;
  }

  /**
   * Get the type code of this object
   */
  getType(): number {
    if (!this.isValid()) {
      throw new Error('Cannot get type of freed SymEngineObject');
    }

    const wasm = getWasmModule();
    return wasm._basic_get_type(this.ptr);
  }

  /**
   * Check equality with another SymEngineObject
   */
  equals(other: SymEngineObject): boolean {
    if (!this.isValid() || !other.isValid()) {
      throw new Error('Cannot compare freed SymEngineObject');
    }

    const wasm = getWasmModule();
    return wasm._basic_eq(this.ptr, other.getPtr()) !== 0;
  }

  /**
   * Get hash code for this object
   */
  hash(): number {
    if (!this.isValid()) {
      throw new Error('Cannot hash freed SymEngineObject');
    }

    const wasm = getWasmModule();
    return wasm._basic_hash(this.ptr);
  }
}

/**
 * Helper to allocate and manage SymEngine Basic objects
 * @returns A new SymEngineObject wrapper
 */
export function createBasic(): SymEngineObject {
  const wasm = getWasmModule();
  const ptr = wasm._basic_new_stack();
  return new SymEngineObject(ptr);
}

/**
 * Check SymEngine exception code and throw JavaScript error if non-zero
 * @param code Exception code returned from SymEngine C API
 * @throws Error if code is not NO_EXCEPTION
 */
export function checkException(code: number): void {
  if (code !== SymEngineException.NO_EXCEPTION) {
    const errorNames = [
      'No error',
      'Runtime error',
      'Division by zero',
      'Not implemented',
      'Domain error',
      'Parse error',
    ];
    throw new Error(`SymEngine error: ${errorNames[code] || 'Unknown error'}`);
  }
}

/**
 * Execute a WASM operation that returns an exception code
 * Automatically checks the exception and throws if non-zero
 * @param operation Function that returns exception code
 * @returns The exception code (always 0 if successful)
 */
export function withExceptionCheck(operation: () => number): void {
  const code = operation();
  checkException(code);
}

/**
 * Create a string in WASM memory from JavaScript string
 * @param str JavaScript string
 * @returns Pointer to C string in WASM memory (caller must free)
 */
export function allocateString(str: string): number {
  const wasm = getWasmModule();
  const length = wasm.lengthBytesUTF8(str) + 1; // +1 for null terminator
  const ptr = wasm._malloc(length);
  wasm.stringToUTF8(str, ptr, length);
  return ptr;
}

/**
 * Free a string allocated in WASM memory
 * @param ptr Pointer to C string
 */
export function freeString(ptr: number): void {
  const wasm = getWasmModule();
  wasm._free(ptr);
}

/**
 * Helper to execute operation with temporary string allocation
 * Automatically frees the string after operation completes
 * @param str JavaScript string
 * @param operation Function that uses the allocated string pointer
 * @returns Result from operation
 */
export function withTempString<T>(str: string, operation: (ptr: number) => T): T {
  const ptr = allocateString(str);
  try {
    return operation(ptr);
  } finally {
    freeString(ptr);
  }
}
