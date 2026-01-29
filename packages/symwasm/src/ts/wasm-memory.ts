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
      wasm._basic_free_heap(this.ptr);
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
 * Uses heap allocation for proper WASM memory management
 * @returns A new SymEngineObject wrapper
 */
export function createBasic(): SymEngineObject {
  const wasm = getWasmModule();
  const ptr = wasm._basic_new_heap();
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

/**
 * Wrapper class for SymEngine CSetBasic container
 * Used for operations that return sets of Basic objects (e.g., free_symbols)
 */
export class SymEngineSet {
  private ptr: number;

  constructor() {
    const wasm = getWasmModule();
    this.ptr = wasm._setbasic_new();
  }

  /**
   * Get the WASM pointer to this set
   */
  getPtr(): number {
    return this.ptr;
  }

  /**
   * Check if this set has been freed
   */
  isValid(): boolean {
    return this.ptr !== 0;
  }

  /**
   * Free the WASM memory for this set
   */
  free(): void {
    if (this.ptr) {
      const wasm = getWasmModule();
      wasm._setbasic_free(this.ptr);
      this.ptr = 0;
    }
  }

  /**
   * Get the number of elements in this set
   */
  size(): number {
    if (!this.isValid()) {
      throw new Error('Cannot get size of freed SymEngineSet');
    }
    const wasm = getWasmModule();
    return wasm._setbasic_size(this.ptr);
  }

  /**
   * Get the nth element from the set
   * @param n Index (0-based)
   * @returns A new SymEngineObject containing the element
   */
  get(n: number): SymEngineObject {
    if (!this.isValid()) {
      throw new Error('Cannot get element from freed SymEngineSet');
    }
    const wasm = getWasmModule();
    const resultPtr = wasm._basic_new_heap();
    wasm._setbasic_get(this.ptr, n, resultPtr);
    return new SymEngineObject(resultPtr);
  }

  /**
   * Convert the set to an array of SymEngineObjects
   * Note: The caller is responsible for freeing the returned objects
   */
  toArray(): SymEngineObject[] {
    const count = this.size();
    const result: SymEngineObject[] = [];
    for (let i = 0; i < count; i++) {
      result.push(this.get(i));
    }
    return result;
  }
}

/**
 * Wrapper class for SymEngine CVecBasic container
 * Used for operations that return vectors of Basic objects (e.g., get_args)
 */
export class SymEngineVec {
  private ptr: number;

  constructor() {
    const wasm = getWasmModule();
    this.ptr = wasm._vecbasic_new();
  }

  /**
   * Get the WASM pointer to this vector
   */
  getPtr(): number {
    return this.ptr;
  }

  /**
   * Check if this vector has been freed
   */
  isValid(): boolean {
    return this.ptr !== 0;
  }

  /**
   * Free the WASM memory for this vector
   */
  free(): void {
    if (this.ptr) {
      const wasm = getWasmModule();
      wasm._vecbasic_free(this.ptr);
      this.ptr = 0;
    }
  }

  /**
   * Get the number of elements in this vector
   */
  size(): number {
    if (!this.isValid()) {
      throw new Error('Cannot get size of freed SymEngineVec');
    }
    const wasm = getWasmModule();
    return wasm._vecbasic_size(this.ptr);
  }

  /**
   * Get the nth element from the vector
   * @param n Index (0-based)
   * @returns A new SymEngineObject containing the element
   */
  get(n: number): SymEngineObject {
    if (!this.isValid()) {
      throw new Error('Cannot get element from freed SymEngineVec');
    }
    const wasm = getWasmModule();
    const resultPtr = wasm._basic_new_heap();
    wasm._vecbasic_get(this.ptr, n, resultPtr);
    return new SymEngineObject(resultPtr);
  }

  /**
   * Append an element to the vector
   * @param obj The SymEngineObject to append
   */
  push(obj: SymEngineObject): void {
    if (!this.isValid()) {
      throw new Error('Cannot push to freed SymEngineVec');
    }
    const wasm = getWasmModule();
    wasm._vecbasic_push_back(this.ptr, obj.getPtr());
  }

  /**
   * Convert the vector to an array of SymEngineObjects
   * Note: The caller is responsible for freeing the returned objects
   */
  toArray(): SymEngineObject[] {
    const count = this.size();
    const result: SymEngineObject[] = [];
    for (let i = 0; i < count; i++) {
      result.push(this.get(i));
    }
    return result;
  }
}

/**
 * Wrapper class for SymEngine CMapBasicBasic container
 * Used for map-based substitution operations
 */
export class SymEngineMap {
  private ptr: number;

  constructor() {
    const wasm = getWasmModule();
    this.ptr = wasm._mapbasicbasic_new();
  }

  /**
   * Get the WASM pointer to this map
   */
  getPtr(): number {
    return this.ptr;
  }

  /**
   * Check if this map has been freed
   */
  isValid(): boolean {
    return this.ptr !== 0;
  }

  /**
   * Free the WASM memory for this map
   */
  free(): void {
    if (this.ptr) {
      const wasm = getWasmModule();
      wasm._mapbasicbasic_free(this.ptr);
      this.ptr = 0;
    }
  }

  /**
   * Insert a key-value pair into the map
   * @param keyPtr WASM pointer to the key Basic
   * @param valuePtr WASM pointer to the value Basic
   */
  insert(keyPtr: number, valuePtr: number): void {
    if (!this.isValid()) {
      throw new Error('Cannot insert into freed SymEngineMap');
    }
    const wasm = getWasmModule();
    wasm._mapbasicbasic_insert(this.ptr, keyPtr, valuePtr);
  }

  /**
   * Get the number of entries in this map
   */
  size(): number {
    if (!this.isValid()) {
      throw new Error('Cannot get size of freed SymEngineMap');
    }
    const wasm = getWasmModule();
    return wasm._mapbasicbasic_size(this.ptr);
  }
}

/**
 * Wrapper class for SymEngine CDenseMatrix objects
 * Manages lifecycle and memory cleanup for dense matrices
 */
export class DenseMatrixObject {
  private ptr: number;

  constructor(ptr: number) {
    this.ptr = ptr;
  }

  /**
   * Get the WASM pointer to this matrix
   */
  getPtr(): number {
    return this.ptr;
  }

  /**
   * Check if this matrix has been freed
   */
  isValid(): boolean {
    return this.ptr !== 0;
  }

  /**
   * Free the WASM memory for this matrix
   */
  free(): void {
    if (this.ptr) {
      const wasm = getWasmModule();
      wasm._dense_matrix_free(this.ptr);
      this.ptr = 0;
    }
  }

  /**
   * Get the string representation of this matrix
   */
  toString(): string {
    if (!this.isValid()) {
      throw new Error('Cannot convert freed DenseMatrixObject to string');
    }
    const wasm = getWasmModule();
    const strPtr = wasm._dense_matrix_str(this.ptr);
    const str = wasm.UTF8ToString(strPtr);
    wasm._free(strPtr);
    return str;
  }

  /**
   * Get the number of rows in this matrix
   */
  rows(): number {
    if (!this.isValid()) {
      throw new Error('Cannot get rows of freed DenseMatrixObject');
    }
    const wasm = getWasmModule();
    return wasm._dense_matrix_rows(this.ptr);
  }

  /**
   * Get the number of columns in this matrix
   */
  cols(): number {
    if (!this.isValid()) {
      throw new Error('Cannot get cols of freed DenseMatrixObject');
    }
    const wasm = getWasmModule();
    return wasm._dense_matrix_cols(this.ptr);
  }
}

/**
 * Create a new empty DenseMatrixObject
 * @returns A new DenseMatrixObject wrapper
 */
export function createDenseMatrix(): DenseMatrixObject {
  const wasm = getWasmModule();
  const ptr = wasm._dense_matrix_new();
  return new DenseMatrixObject(ptr);
}

/**
 * Create a DenseMatrixObject with specified dimensions
 * @param rows Number of rows
 * @param cols Number of columns
 * @returns A new DenseMatrixObject wrapper
 */
export function createDenseMatrixWithSize(rows: number, cols: number): DenseMatrixObject {
  const wasm = getWasmModule();
  const ptr = wasm._dense_matrix_new_rows_cols(rows, cols);
  return new DenseMatrixObject(ptr);
}
