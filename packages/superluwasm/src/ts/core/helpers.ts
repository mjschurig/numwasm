/**
 * Memory and Array Helper Functions for SuperLU WASM
 *
 * Provides utilities for allocating, reading, and writing data
 * to/from the WebAssembly heap.
 *
 * @module core/helpers
 */

import type { SuperLUModule } from '../types.js';
import type { RealArray, ComplexArray, IntArray, DataType } from '../high-level-types.js';

// ============================================================
// Constants
// ============================================================

/** Size of a double (float64) in bytes */
export const DOUBLE_SIZE = 8;

/** Size of a float (float32) in bytes */
export const FLOAT_SIZE = 4;

/** Size of an int32 in bytes */
export const INT_SIZE = 4;

/** Size of superlu_options_t structure (approximate, platform-dependent) */
export const OPTIONS_STRUCT_SIZE = 256;

/** Size of SuperLUStat_t structure (approximate) */
export const STAT_STRUCT_SIZE = 256;

/** Size of SuperMatrix structure */
export const SUPERMATRIX_SIZE = 48;

/** Size of GlobalLU_t structure */
export const GLOBALLU_SIZE = 128;

/** Size of mem_usage_t structure */
export const MEM_USAGE_SIZE = 32;

// ============================================================
// Memory Allocation Helpers
// ============================================================

/**
 * Allocate a Float64Array in WASM memory and optionally initialize with values.
 *
 * @param module - SuperLU WASM module
 * @param values - Initial values, or null for uninitialized allocation
 * @param size - Size to allocate (defaults to values.length)
 * @returns Pointer to allocated memory
 * @throws Error if allocation fails
 */
export function allocateDoubles(
  module: SuperLUModule,
  values: RealArray | null,
  size?: number
): number {
  const len = size ?? (values?.length ?? 0);
  if (len === 0) return 0;

  const ptr = module._malloc(len * DOUBLE_SIZE);
  if (ptr === 0) {
    throw new Error('Failed to allocate WASM memory for doubles');
  }

  if (values) {
    const baseIdx = ptr >> 3;
    for (let i = 0; i < values.length; i++) {
      module.HEAPF64[baseIdx + i] = values[i];
    }
  }

  return ptr;
}

/**
 * Allocate a Float32Array in WASM memory and optionally initialize with values.
 *
 * @param module - SuperLU WASM module
 * @param values - Initial values, or null for uninitialized allocation
 * @param size - Size to allocate (defaults to values.length)
 * @returns Pointer to allocated memory
 * @throws Error if allocation fails
 */
export function allocateFloats(
  module: SuperLUModule,
  values: RealArray | null,
  size?: number
): number {
  const len = size ?? (values?.length ?? 0);
  if (len === 0) return 0;

  const ptr = module._malloc(len * FLOAT_SIZE);
  if (ptr === 0) {
    throw new Error('Failed to allocate WASM memory for floats');
  }

  if (values) {
    const baseIdx = ptr >> 2;
    for (let i = 0; i < values.length; i++) {
      module.HEAPF32[baseIdx + i] = values[i];
    }
  }

  return ptr;
}

/**
 * Allocate an Int32Array in WASM memory and optionally initialize with values.
 *
 * @param module - SuperLU WASM module
 * @param values - Initial values, or null for uninitialized allocation
 * @param size - Size to allocate (defaults to values.length)
 * @returns Pointer to allocated memory
 * @throws Error if allocation fails
 */
export function allocateInts(
  module: SuperLUModule,
  values: IntArray | null,
  size?: number
): number {
  const len = size ?? (values?.length ?? 0);
  if (len === 0) return 0;

  const ptr = module._malloc(len * INT_SIZE);
  if (ptr === 0) {
    throw new Error('Failed to allocate WASM memory for integers');
  }

  if (values) {
    const baseIdx = ptr >> 2;
    for (let i = 0; i < values.length; i++) {
      module.HEAP32[baseIdx + i] = values[i];
    }
  }

  return ptr;
}

/**
 * Allocate values based on data type.
 *
 * @param module - SuperLU WASM module
 * @param values - Initial values
 * @param dtype - Data type
 * @returns Pointer to allocated memory
 */
export function allocateValues(
  module: SuperLUModule,
  values: RealArray | ComplexArray | null,
  dtype: DataType,
  size?: number
): number {
  switch (dtype) {
    case 'float32':
    case 'complex64':
      return allocateFloats(module, values as RealArray, size);
    case 'float64':
    case 'complex128':
    default:
      return allocateDoubles(module, values as RealArray, size);
  }
}

// ============================================================
// Memory Reading Helpers
// ============================================================

/**
 * Read doubles from WASM memory.
 *
 * @param module - SuperLU WASM module
 * @param ptr - Pointer to read from
 * @param length - Number of elements to read
 * @returns New Float64Array with the values
 */
export function readDoubles(
  module: SuperLUModule,
  ptr: number,
  length: number
): Float64Array {
  const result = new Float64Array(length);
  const baseIdx = ptr >> 3;
  for (let i = 0; i < length; i++) {
    result[i] = module.HEAPF64[baseIdx + i];
  }
  return result;
}

/**
 * Read floats from WASM memory.
 *
 * @param module - SuperLU WASM module
 * @param ptr - Pointer to read from
 * @param length - Number of elements to read
 * @returns New Float32Array with the values
 */
export function readFloats(
  module: SuperLUModule,
  ptr: number,
  length: number
): Float32Array {
  const result = new Float32Array(length);
  const baseIdx = ptr >> 2;
  for (let i = 0; i < length; i++) {
    result[i] = module.HEAPF32[baseIdx + i];
  }
  return result;
}

/**
 * Read integers from WASM memory.
 *
 * @param module - SuperLU WASM module
 * @param ptr - Pointer to read from
 * @param length - Number of elements to read
 * @returns New Int32Array with the values
 */
export function readInts(
  module: SuperLUModule,
  ptr: number,
  length: number
): Int32Array {
  const result = new Int32Array(length);
  const baseIdx = ptr >> 2;
  for (let i = 0; i < length; i++) {
    result[i] = module.HEAP32[baseIdx + i];
  }
  return result;
}

/**
 * Read a single integer from WASM memory.
 *
 * @param module - SuperLU WASM module
 * @param ptr - Pointer to read from
 * @returns The integer value
 */
export function readInt(module: SuperLUModule, ptr: number): number {
  return module.HEAP32[ptr >> 2];
}

/**
 * Read a single double from WASM memory.
 *
 * @param module - SuperLU WASM module
 * @param ptr - Pointer to read from
 * @returns The double value
 */
export function readDouble(module: SuperLUModule, ptr: number): number {
  return module.HEAPF64[ptr >> 3];
}

/**
 * Read a single byte/char from WASM memory.
 *
 * @param module - SuperLU WASM module
 * @param ptr - Pointer to read from
 * @returns The byte value
 */
export function readByte(module: SuperLUModule, ptr: number): number {
  return module.HEAP8[ptr];
}

/**
 * Read values based on data type.
 *
 * @param module - SuperLU WASM module
 * @param ptr - Pointer to read from
 * @param length - Number of elements
 * @param dtype - Data type
 * @returns Typed array with values
 */
export function readValues(
  module: SuperLUModule,
  ptr: number,
  length: number,
  dtype: DataType
): Float64Array | Float32Array {
  switch (dtype) {
    case 'float32':
    case 'complex64':
      return readFloats(module, ptr, length);
    case 'float64':
    case 'complex128':
    default:
      return readDoubles(module, ptr, length);
  }
}

// ============================================================
// Memory Writing Helpers
// ============================================================

/**
 * Write doubles to WASM memory.
 *
 * @param module - SuperLU WASM module
 * @param ptr - Pointer to write to
 * @param values - Values to write
 */
export function writeDoubles(
  module: SuperLUModule,
  ptr: number,
  values: RealArray
): void {
  const baseIdx = ptr >> 3;
  for (let i = 0; i < values.length; i++) {
    module.HEAPF64[baseIdx + i] = values[i];
  }
}

/**
 * Write floats to WASM memory.
 *
 * @param module - SuperLU WASM module
 * @param ptr - Pointer to write to
 * @param values - Values to write
 */
export function writeFloats(
  module: SuperLUModule,
  ptr: number,
  values: RealArray
): void {
  const baseIdx = ptr >> 2;
  for (let i = 0; i < values.length; i++) {
    module.HEAPF32[baseIdx + i] = values[i];
  }
}

/**
 * Write integers to WASM memory.
 *
 * @param module - SuperLU WASM module
 * @param ptr - Pointer to write to
 * @param values - Values to write
 */
export function writeInts(
  module: SuperLUModule,
  ptr: number,
  values: IntArray
): void {
  const baseIdx = ptr >> 2;
  for (let i = 0; i < values.length; i++) {
    module.HEAP32[baseIdx + i] = values[i];
  }
}

/**
 * Write a single integer to WASM memory.
 *
 * @param module - SuperLU WASM module
 * @param ptr - Pointer to write to
 * @param value - Value to write
 */
export function writeInt(module: SuperLUModule, ptr: number, value: number): void {
  module.HEAP32[ptr >> 2] = value;
}

/**
 * Write a single double to WASM memory.
 *
 * @param module - SuperLU WASM module
 * @param ptr - Pointer to write to
 * @param value - Value to write
 */
export function writeDouble(module: SuperLUModule, ptr: number, value: number): void {
  module.HEAPF64[ptr >> 3] = value;
}

/**
 * Write a single byte to WASM memory.
 *
 * @param module - SuperLU WASM module
 * @param ptr - Pointer to write to
 * @param value - Value to write
 */
export function writeByte(module: SuperLUModule, ptr: number, value: number): void {
  module.HEAP8[ptr] = value;
}

// ============================================================
// Memory Management
// ============================================================

/**
 * Free multiple pointers at once.
 * Safely handles zero/null pointers (skips them).
 *
 * @param module - SuperLU WASM module
 * @param ptrs - Array of pointers to free
 */
export function freeAll(module: SuperLUModule, ptrs: number[]): void {
  for (const ptr of ptrs) {
    if (ptr !== 0) {
      module._free(ptr);
    }
  }
}

/**
 * RAII-style resource management.
 * Executes a function with allocated resources and ensures cleanup.
 *
 * @param module - SuperLU WASM module
 * @param allocations - Object mapping names to allocation functions
 * @param fn - Function to execute with allocated pointers
 * @returns Result of fn
 */
export async function withAllocations<T>(
  module: SuperLUModule,
  allocations: Record<string, () => number>,
  fn: (ptrs: Record<string, number>) => T | Promise<T>
): Promise<T> {
  const ptrs: Record<string, number> = {};
  const allocated: number[] = [];

  try {
    // Allocate all resources
    for (const [name, allocFn] of Object.entries(allocations)) {
      const ptr = allocFn();
      ptrs[name] = ptr;
      if (ptr !== 0) {
        allocated.push(ptr);
      }
    }

    // Execute function
    return await fn(ptrs);
  } finally {
    // Free all allocated resources
    freeAll(module, allocated);
  }
}

// ============================================================
// Conversion Helpers
// ============================================================

/**
 * Convert a JavaScript array to Float64Array.
 */
export function toFloat64Array(arr: RealArray): Float64Array {
  if (arr instanceof Float64Array) return arr;
  return new Float64Array(arr);
}

/**
 * Convert a JavaScript array to Float32Array.
 */
export function toFloat32Array(arr: RealArray): Float32Array {
  if (arr instanceof Float32Array) return arr;
  return new Float32Array(arr);
}

/**
 * Convert a JavaScript array to Int32Array.
 */
export function toInt32Array(arr: IntArray): Int32Array {
  if (arr instanceof Int32Array) return arr;
  return new Int32Array(arr);
}

/**
 * Get element size in bytes for a data type.
 */
export function getElementSize(dtype: DataType): number {
  switch (dtype) {
    case 'float32':
      return FLOAT_SIZE;
    case 'float64':
      return DOUBLE_SIZE;
    case 'complex64':
      return 2 * FLOAT_SIZE;
    case 'complex128':
      return 2 * DOUBLE_SIZE;
    default:
      return DOUBLE_SIZE;
  }
}

/**
 * Check if a data type is complex.
 */
export function isComplexType(dtype: DataType): boolean {
  return dtype === 'complex64' || dtype === 'complex128';
}

// ============================================================
// Error Handling
// ============================================================

/**
 * SuperLU error codes and messages.
 */
export const SUPERLU_ERRORS: Record<number, string> = {
  0: 'Success',
  [-1]: 'Illegal value in argument 1',
  [-2]: 'Illegal value in argument 2',
  [-3]: 'Illegal value in argument 3',
  [-4]: 'Illegal value in argument 4',
  [-5]: 'Illegal value in argument 5',
  [-6]: 'Illegal value in argument 6',
  [-7]: 'Illegal value in argument 7',
  [-8]: 'Illegal value in argument 8',
  [-9]: 'Illegal value in argument 9',
  [-10]: 'Illegal value in argument 10',
};

/**
 * Get error message for a SuperLU info code.
 *
 * @param info - Info code from SuperLU routine
 * @param routine - Name of the routine (for context)
 * @returns Human-readable error message
 */
export function getSuperLUErrorMessage(info: number, routine?: string): string {
  if (info === 0) return 'Success';

  if (info < 0) {
    return SUPERLU_ERRORS[info] ?? `Illegal value in argument ${-info}`;
  }

  // Positive info codes indicate singularity
  if (routine?.includes('gstrf') || routine?.includes('gssv')) {
    return `Matrix is singular: U(${info},${info}) is exactly zero. ` +
           `The factorization has been completed, but U is exactly singular.`;
  }

  return `SuperLU${routine ? ` ${routine}` : ''} returned info=${info}`;
}

/**
 * Check SuperLU return code and throw if error.
 *
 * @param info - Info code from SuperLU routine
 * @param routine - Name of the routine
 * @throws Error if info indicates failure
 */
export function checkSuperLUError(info: number, routine: string): void {
  if (info !== 0) {
    throw new Error(getSuperLUErrorMessage(info, routine));
  }
}
