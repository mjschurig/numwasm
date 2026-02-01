/**
 * ARPACK Helper Functions
 *
 * Utility functions for work array sizing and memory management.
 */

import type { ARPACKModule } from './types.js';
import {
  DSAUPD_ERRORS,
  DNAUPD_ERRORS,
  DSEUPD_ERRORS,
  DNEUPD_ERRORS,
  ZNAUPD_ERRORS,
  ZNEUPD_ERRORS,
} from './types.js';

const DOUBLE_SIZE = 8;
const INT_SIZE = 4;

/** Array type that accepts both number[] and Float64Array */
export type RealArray = number[] | Float64Array;

// ============================================================
// WORK ARRAY SIZING
// ============================================================

/**
 * Calculate required workl size for symmetric problems (dsaupd).
 * @param ncv - Number of Lanczos vectors
 * @returns Required size in doubles
 */
export function dsaupdWorklSize(ncv: number): number {
  // From ARPACK docs: lworkl >= ncv*(ncv+8)
  return ncv * (ncv + 8);
}

/**
 * Calculate required workl size for non-symmetric problems (dnaupd).
 * @param ncv - Number of Arnoldi vectors
 * @returns Required size in doubles
 */
export function dnaupdWorklSize(ncv: number): number {
  // From ARPACK docs: lworkl >= 3*ncv² + 6*ncv
  return 3 * ncv * ncv + 6 * ncv;
}

/**
 * Calculate required workev size for non-symmetric extraction (dneupd).
 * @param ncv - Number of Arnoldi vectors
 * @returns Required size in doubles
 */
export function dneupdWorkevSize(ncv: number): number {
  return 3 * ncv;
}

/**
 * Calculate required workl size for complex problems (znaupd).
 * @param ncv - Number of Arnoldi vectors
 * @returns Required size in complex numbers (multiply by 2 for doubles)
 */
export function znaupdWorklSize(ncv: number): number {
  // From ARPACK docs: lworkl >= 3*ncv² + 5*ncv (in complex numbers)
  return 3 * ncv * ncv + 5 * ncv;
}

/**
 * Calculate required workev size for complex extraction (zneupd).
 * @param ncv - Number of Arnoldi vectors
 * @returns Required size in complex numbers (multiply by 2 for doubles)
 */
export function zneupdWorkevSize(ncv: number): number {
  return 2 * ncv;
}

/**
 * Calculate required rwork size for complex problems (znaupd).
 * @param ncv - Number of Arnoldi vectors
 * @returns Required size in doubles (real array)
 */
export function znaupdRworkSize(ncv: number): number {
  return ncv;
}

/**
 * Calculate default ncv value.
 * @param n - Matrix dimension
 * @param nev - Number of eigenvalues requested
 * @param symmetric - Whether the problem is symmetric
 * @returns Recommended ncv value
 */
export function defaultNcv(
  n: number,
  nev: number,
  symmetric: boolean
): number {
  if (symmetric) {
    // ARPACK recommends ncv >= 2*nev for symmetric
    return Math.min(n, Math.max(2 * nev, 20));
  } else {
    // ARPACK recommends ncv >= 2*nev+1 for non-symmetric
    return Math.min(n, Math.max(2 * nev + 1, 20));
  }
}

// ============================================================
// ERROR MESSAGE UTILITIES
// ============================================================

/**
 * Get error message for dsaupd info code.
 */
export function getDsaupdMessage(info: number): string {
  return DSAUPD_ERRORS[info] ?? `Unknown error code: ${info}`;
}

/**
 * Get error message for dnaupd info code.
 */
export function getDnaupdMessage(info: number): string {
  return DNAUPD_ERRORS[info] ?? `Unknown error code: ${info}`;
}

/**
 * Get error message for dseupd info code.
 */
export function getDseupdMessage(info: number): string {
  return DSEUPD_ERRORS[info] ?? `Unknown error code: ${info}`;
}

/**
 * Get error message for dneupd info code.
 */
export function getDneupdMessage(info: number): string {
  return DNEUPD_ERRORS[info] ?? `Unknown error code: ${info}`;
}

/**
 * Get error message for znaupd info code.
 */
export function getZnaupdMessage(info: number): string {
  return ZNAUPD_ERRORS[info] ?? `Unknown error code: ${info}`;
}

/**
 * Get error message for zneupd info code.
 */
export function getZneupdMessage(info: number): string {
  return ZNEUPD_ERRORS[info] ?? `Unknown error code: ${info}`;
}

// ============================================================
// MEMORY HELPERS
// ============================================================

/**
 * Allocate a double array in WASM memory and optionally initialize with values.
 * @param Module - ARPACK WASM module
 * @param values - Initial values (number[] or Float64Array), or null for uninitialized
 * @param size - Size to allocate (defaults to values.length if values provided)
 * @returns Pointer to allocated memory
 */
export function allocateDoubles(
  Module: ARPACKModule,
  values: RealArray | null,
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
 * Allocate an integer array in WASM memory and optionally initialize with values.
 * @param Module - ARPACK WASM module
 * @param values - Initial values, or null for uninitialized
 * @param size - Size to allocate (defaults to values.length if values provided)
 * @returns Pointer to allocated memory
 */
export function allocateInts(
  Module: ARPACKModule,
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
 * Write a string to WASM memory.
 * Allocates memory for the string plus null terminator.
 * @param Module - ARPACK WASM module
 * @param str - String to write
 * @returns Pointer to allocated memory
 */
export function allocateString(Module: ARPACKModule, str: string): number {
  const ptr = Module._malloc(str.length + 1);
  if (ptr === 0) {
    throw new Error('Failed to allocate WASM memory for string');
  }

  for (let i = 0; i < str.length; i++) {
    Module.HEAPU8[ptr + i] = str.charCodeAt(i);
  }
  Module.HEAPU8[ptr + str.length] = 0; // null terminator

  return ptr;
}

/**
 * Read a double array from WASM memory.
 * @param Module - ARPACK WASM module
 * @param ptr - Pointer to the array
 * @param length - Number of elements to read
 * @returns New Float64Array containing the values
 */
export function readDoubles(
  Module: ARPACKModule,
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
 * Read an integer from WASM memory.
 * @param Module - ARPACK WASM module
 * @param ptr - Pointer to the integer
 * @returns The integer value
 */
export function readInt(Module: ARPACKModule, ptr: number): number {
  return Module.HEAP32[ptr >> 2];
}

/**
 * Write doubles to WASM memory.
 * @param Module - ARPACK WASM module
 * @param ptr - Pointer to write to
 * @param values - Values to write
 */
export function writeDoubles(
  Module: ARPACKModule,
  ptr: number,
  values: RealArray
): void {
  const baseIdx = ptr >> 3;
  for (let i = 0; i < values.length; i++) {
    Module.HEAPF64[baseIdx + i] = values[i];
  }
}

/**
 * Write an integer to WASM memory.
 * @param Module - ARPACK WASM module
 * @param ptr - Pointer to write to
 * @param value - Value to write
 */
export function writeInt(
  Module: ARPACKModule,
  ptr: number,
  value: number
): void {
  Module.HEAP32[ptr >> 2] = value;
}

/**
 * Free multiple pointers at once.
 * Safely handles zero pointers (skips them).
 * @param Module - ARPACK WASM module
 * @param ptrs - Array of pointers to free
 */
export function freeAll(Module: ARPACKModule, ptrs: number[]): void {
  for (const ptr of ptrs) {
    if (ptr !== 0) {
      Module._free(ptr);
    }
  }
}

// ============================================================
// COMPLEX ARRAY HELPERS
// ============================================================

/** Complex array type with separate real and imaginary parts */
export interface ComplexArray {
  re: Float64Array;
  im: Float64Array;
}

/**
 * Allocate a complex array in WASM memory (interleaved format).
 * Complex numbers are stored as [re0, im0, re1, im1, ...].
 * @param Module - ARPACK WASM module
 * @param values - Initial values as ComplexArray (re/im arrays), or null for uninitialized
 * @param size - Number of complex elements to allocate
 * @returns Pointer to allocated memory
 */
export function allocateComplex(
  Module: ARPACKModule,
  values: ComplexArray | null,
  size?: number
): number {
  const len = size ?? (values?.re.length ?? 0);
  if (len === 0) return 0;

  // Each complex number is 2 doubles (16 bytes)
  const ptr = Module._malloc(len * 2 * DOUBLE_SIZE);
  if (ptr === 0) {
    throw new Error('Failed to allocate WASM memory for complex array');
  }

  if (values) {
    const baseIdx = ptr >> 3;
    for (let i = 0; i < values.re.length; i++) {
      Module.HEAPF64[baseIdx + 2 * i] = values.re[i];
      Module.HEAPF64[baseIdx + 2 * i + 1] = values.im[i];
    }
  }

  return ptr;
}

/**
 * Read a complex array from WASM memory (interleaved format).
 * @param Module - ARPACK WASM module
 * @param ptr - Pointer to the array
 * @param length - Number of complex elements to read
 * @returns ComplexArray with separate re and im Float64Arrays
 */
export function readComplex(
  Module: ARPACKModule,
  ptr: number,
  length: number
): ComplexArray {
  const re = new Float64Array(length);
  const im = new Float64Array(length);
  const baseIdx = ptr >> 3;

  for (let i = 0; i < length; i++) {
    re[i] = Module.HEAPF64[baseIdx + 2 * i];
    im[i] = Module.HEAPF64[baseIdx + 2 * i + 1];
  }

  return { re, im };
}

/**
 * Write complex values to WASM memory (interleaved format).
 * @param Module - ARPACK WASM module
 * @param ptr - Pointer to write to
 * @param values - ComplexArray to write
 */
export function writeComplex(
  Module: ARPACKModule,
  ptr: number,
  values: ComplexArray
): void {
  const baseIdx = ptr >> 3;
  for (let i = 0; i < values.re.length; i++) {
    Module.HEAPF64[baseIdx + 2 * i] = values.re[i];
    Module.HEAPF64[baseIdx + 2 * i + 1] = values.im[i];
  }
}
