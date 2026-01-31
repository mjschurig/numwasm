/**
 * QUADPACK Helper Functions
 *
 * Utility functions for workspace sizing, memory management, and error handling.
 */

import type { QUADPACKModule } from './types.js';
import { QUADPACK_ERRORS } from './types.js';

const DOUBLE_SIZE = 8;
const INT_SIZE = 4;

// ============================================================
// WORKSPACE SIZING
// ============================================================

/**
 * Calculate workspace size for basic adaptive routines (DQAGS, DQAG, etc.).
 * @param limit - Maximum number of subdivisions
 * @returns Required size for work array (in doubles)
 */
export function quadWorkSize(limit: number): number {
  return 4 * limit;
}

/**
 * Calculate integer workspace size for basic adaptive routines.
 * @param limit - Maximum number of subdivisions
 * @returns Required size for iwork array (in integers)
 */
export function quadIworkSize(limit: number): number {
  return limit;
}

/**
 * Calculate workspace size for oscillatory routines (DQAWO).
 * @param limit - Maximum number of subdivisions
 * @param maxp1 - Maximum number of Chebyshev moments
 * @returns Required size for work array (in doubles)
 */
export function quadOscWorkSize(limit: number, maxp1: number): number {
  // From QUADPACK docs: lenw >= leniw*2 + maxp1*25
  // where leniw = limit
  return limit * 2 + maxp1 * 25;
}

/**
 * Calculate integer workspace size for oscillatory routines (DQAWO).
 * @param limit - Maximum number of subdivisions
 * @returns Required size for iwork array
 */
export function quadOscIworkSize(limit: number): number {
  return limit;
}

/**
 * Calculate workspace size for Fourier routines (DQAWF).
 * @param limlst - Maximum number of cycles
 * @param limit - Maximum subdivisions per cycle
 * @param maxp1 - Maximum Chebyshev moments
 * @returns Required size for work array
 */
export function quadFourierWorkSize(
  limlst: number,
  limit: number,
  maxp1: number
): number {
  // From QUADPACK: lenw >= leniw*2 + maxp1*25
  // where leniw >= limlst + limit
  return (limlst + limit) * 2 + maxp1 * 25;
}

/**
 * Calculate integer workspace size for Fourier routines (DQAWF).
 * @param limlst - Maximum number of cycles
 * @param limit - Maximum subdivisions per cycle
 * @returns Required size for iwork array
 */
export function quadFourierIworkSize(limlst: number, limit: number): number {
  return limlst + limit;
}

/**
 * Calculate workspace size for break point routines (DQAGP).
 * @param limit - Maximum number of subdivisions
 * @param npts - Number of break points (not including a and b)
 * @returns Required size for work array (in doubles)
 */
export function quadBreakWorkSize(limit: number, npts: number): number {
  // From QUADPACK: lenw >= limit*4 + npts2*4 where npts2 = npts + 2
  const npts2 = npts + 2;
  return limit * 4 + npts2 * 4;
}

/**
 * Calculate integer workspace size for break point routines (DQAGP).
 * @param limit - Maximum number of subdivisions
 * @param npts - Number of break points (not including a and b)
 * @returns Required size for iwork array
 */
export function quadBreakIworkSize(limit: number, npts: number): number {
  // From QUADPACK: leniw >= limit*2 + npts2 where npts2 = npts + 2
  const npts2 = npts + 2;
  return limit * 2 + npts2;
}

// ============================================================
// ERROR MESSAGE UTILITIES
// ============================================================

/**
 * Get human-readable error message for QUADPACK error code.
 * @param ier - QUADPACK error code
 * @returns Human-readable description
 */
export function getQuadpackMessage(ier: number): string {
  return QUADPACK_ERRORS[ier] ?? `Unknown error code: ${ier}`;
}

// ============================================================
// MEMORY HELPERS
// ============================================================

/**
 * Allocate a single double in WASM memory and initialize it.
 * @param Module - QUADPACK WASM module
 * @param value - Value to store
 * @returns Pointer to allocated memory
 */
export function allocateDouble(
  Module: QUADPACKModule,
  value: number
): number {
  const ptr = Module._malloc(DOUBLE_SIZE);
  if (ptr === 0) {
    throw new Error('Failed to allocate WASM memory for double');
  }
  Module.HEAPF64[ptr >> 3] = value;
  return ptr;
}

/**
 * Allocate a double array in WASM memory and optionally initialize with values.
 * @param Module - QUADPACK WASM module
 * @param values - Initial values, or null for uninitialized
 * @param size - Size to allocate (defaults to values.length if values provided)
 * @returns Pointer to allocated memory
 */
export function allocateDoubles(
  Module: QUADPACKModule,
  values: number[] | Float64Array | null,
  size?: number
): number {
  const len = size ?? (values?.length ?? 0);
  if (len === 0) return 0;

  const ptr = Module._malloc(len * DOUBLE_SIZE);
  if (ptr === 0) {
    throw new Error('Failed to allocate WASM memory for doubles');
  }

  if (values) {
    const baseIdx = ptr >> 3;
    for (let i = 0; i < values.length; i++) {
      Module.HEAPF64[baseIdx + i] = values[i];
    }
  }

  return ptr;
}

/**
 * Allocate a single integer in WASM memory and initialize it.
 * @param Module - QUADPACK WASM module
 * @param value - Value to store
 * @returns Pointer to allocated memory
 */
export function allocateInt(Module: QUADPACKModule, value: number): number {
  const ptr = Module._malloc(INT_SIZE);
  if (ptr === 0) {
    throw new Error('Failed to allocate WASM memory for integer');
  }
  Module.HEAP32[ptr >> 2] = value;
  return ptr;
}

/**
 * Allocate an integer array in WASM memory and optionally initialize with values.
 * @param Module - QUADPACK WASM module
 * @param values - Initial values, or null for uninitialized
 * @param size - Size to allocate (defaults to values.length if values provided)
 * @returns Pointer to allocated memory
 */
export function allocateInts(
  Module: QUADPACKModule,
  values: number[] | Int32Array | null,
  size?: number
): number {
  const len = size ?? (values?.length ?? 0);
  if (len === 0) return 0;

  const ptr = Module._malloc(len * INT_SIZE);
  if (ptr === 0) {
    throw new Error('Failed to allocate WASM memory for integers');
  }

  if (values) {
    const baseIdx = ptr >> 2;
    for (let i = 0; i < values.length; i++) {
      Module.HEAP32[baseIdx + i] = values[i];
    }
  }

  return ptr;
}

/**
 * Read a double from WASM memory.
 * @param Module - QUADPACK WASM module
 * @param ptr - Pointer to the double
 * @returns The value
 */
export function readDouble(Module: QUADPACKModule, ptr: number): number {
  return Module.HEAPF64[ptr >> 3];
}

/**
 * Read an integer from WASM memory.
 * @param Module - QUADPACK WASM module
 * @param ptr - Pointer to the integer
 * @returns The value
 */
export function readInt(Module: QUADPACKModule, ptr: number): number {
  return Module.HEAP32[ptr >> 2];
}

/**
 * Read a double array from WASM memory.
 * @param Module - QUADPACK WASM module
 * @param ptr - Pointer to the array
 * @param length - Number of elements to read
 * @returns New Float64Array containing the values
 */
export function readDoubles(
  Module: QUADPACKModule,
  ptr: number,
  length: number
): Float64Array {
  const result = new Float64Array(length);
  const baseIdx = ptr >> 3;
  for (let i = 0; i < length; i++) {
    result[i] = Module.HEAPF64[baseIdx + i];
  }
  return result;
}

/**
 * Free multiple pointers at once.
 * Safely handles zero pointers (skips them).
 * @param Module - QUADPACK WASM module
 * @param ptrs - Array of pointers to free
 */
export function freeAll(Module: QUADPACKModule, ptrs: number[]): void {
  for (const ptr of ptrs) {
    if (ptr !== 0) {
      Module._free(ptr);
    }
  }
}

// ============================================================
// CALLBACK MANAGEMENT
// ============================================================

/**
 * Register an integrand function as a WASM callback.
 *
 * The integrand signature is `double f(double x)`, which maps to
 * Emscripten signature 'dd' (double return, double parameter).
 *
 * @param Module - QUADPACK WASM module
 * @param f - JavaScript integrand function
 * @returns Function pointer for use with QUADPACK routines
 */
export function registerIntegrand(
  Module: QUADPACKModule,
  f: (x: number) => number
): number {
  // QUADPACK expects f(x: *double) -> double via pointer
  // The f2c convention passes x by pointer, so we need a wrapper
  const wrapper = (xPtr: number): number => {
    const x = Module.HEAPF64[xPtr >> 3];
    return f(x);
  };
  return Module.addFunction(wrapper, 'di');
}

/**
 * Unregister a previously registered integrand function.
 * @param Module - QUADPACK WASM module
 * @param ptr - Function pointer from registerIntegrand
 */
export function unregisterIntegrand(
  Module: QUADPACKModule,
  ptr: number
): void {
  Module.removeFunction(ptr);
}
