/**
 * LAPACK Helper Functions
 *
 * Utility functions for memory management and array conversion.
 */

import type { LAPACKModule } from './types.js';

const DOUBLE_SIZE = 8;
const INT_SIZE = 4;

/** Array type that accepts both number[] and Float64Array */
export type RealArray = number[] | Float64Array;

/** Complex array type - interleaved real/imaginary pairs */
export type ComplexArray = number[] | Float64Array;

// ============================================================
// CHARACTER CODE CONSTANTS
// ============================================================

/** Character codes for LAPACK routines (Fortran expects single chars as integers) */
export const CHAR = {
  // Transpose options
  N: 78, // 'N' - No transpose
  T: 84, // 'T' - Transpose
  C: 67, // 'C' - Conjugate transpose

  // Triangle selection
  U: 85, // 'U' - Upper triangular
  L: 76, // 'L' - Lower triangular

  // Diagonal type
  // N: 78, // 'N' - Non-unit diagonal (already defined)
  UNIT: 85, // 'U' - Unit diagonal

  // Norm types
  ONE: 49, // '1' - One norm
  INF: 73, // 'I' - Infinity norm
  FRO: 70, // 'F' - Frobenius norm
  MAX: 77, // 'M' - Max absolute value

  // Side
  LEFT: 76, // 'L' - Left side
  RIGHT: 82, // 'R' - Right side

  // Job options
  V: 86, // 'V' - Compute vectors
  A: 65, // 'A' - All
  S: 83, // 'S' - Some/Subset
  O: 79, // 'O' - Overwrite
} as const;

// ============================================================
// ERROR MESSAGES
// ============================================================

/** Error messages for LAPACK info codes */
export const LAPACK_ERRORS: Record<string, Record<number, string>> = {
  DGESV: {
    0: 'Success',
    [-1]: 'Illegal value for N',
    [-2]: 'Illegal value for NRHS',
    [-3]: 'Illegal value for A',
    [-4]: 'Illegal value for LDA',
    [-5]: 'Illegal value for IPIV',
    [-6]: 'Illegal value for B',
    [-7]: 'Illegal value for LDB',
  },
  DPOSV: {
    0: 'Success',
    [-1]: 'Illegal value for UPLO',
    [-2]: 'Illegal value for N',
    [-3]: 'Illegal value for NRHS',
    [-4]: 'Illegal value for A',
    [-5]: 'Illegal value for LDA',
    [-6]: 'Illegal value for B',
    [-7]: 'Illegal value for LDB',
  },
  DTRTRS: {
    0: 'Success',
    [-1]: 'Illegal value for UPLO',
    [-2]: 'Illegal value for TRANS',
    [-3]: 'Illegal value for DIAG',
    [-4]: 'Illegal value for N',
    [-5]: 'Illegal value for NRHS',
    [-6]: 'Illegal value for A',
    [-7]: 'Illegal value for LDA',
    [-8]: 'Illegal value for B',
    [-9]: 'Illegal value for LDB',
  },
};

/**
 * Get error message for LAPACK info code.
 */
export function getLapackErrorMessage(routine: string, info: number): string {
  if (info === 0) return 'Success';
  if (info < 0) {
    const routineErrors = LAPACK_ERRORS[routine];
    if (routineErrors && routineErrors[info]) {
      return routineErrors[info];
    }
    return `Illegal value for argument ${-info}`;
  }
  // Positive info codes are routine-specific
  if (routine === 'DGESV' || routine === 'ZGESV') {
    return `U(${info},${info}) is exactly zero. The factorization has been completed, but the factor U is exactly singular, so the solution could not be computed.`;
  }
  if (routine === 'DPOSV' || routine === 'ZPOSV') {
    return `The leading minor of order ${info} is not positive definite, so the factorization could not be completed.`;
  }
  if (routine === 'DTRTRS' || routine === 'ZTRTRS') {
    return `The diagonal element ${info} of A is zero, indicating a singular matrix.`;
  }
  return `LAPACK routine ${routine} returned info=${info}`;
}

// ============================================================
// MEMORY ALLOCATION HELPERS
// ============================================================

/**
 * Allocate a double array in WASM memory and optionally initialize with values.
 * @param Module - LAPACK WASM module
 * @param values - Initial values (number[] or Float64Array), or null for uninitialized
 * @param size - Size to allocate (defaults to values.length if values provided)
 * @returns Pointer to allocated memory
 */
export function allocateDoubles(
  Module: LAPACKModule,
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
 * @param Module - LAPACK WASM module
 * @param values - Initial values, or null for uninitialized
 * @param size - Size to allocate (defaults to values.length if values provided)
 * @returns Pointer to allocated memory
 */
export function allocateInts(
  Module: LAPACKModule,
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
 * Read a double array from WASM memory.
 * @param Module - LAPACK WASM module
 * @param ptr - Pointer to the array
 * @param length - Number of elements to read
 * @returns New Float64Array containing the values
 */
export function readDoubles(
  Module: LAPACKModule,
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
 * @param Module - LAPACK WASM module
 * @param ptr - Pointer to the integer
 * @returns The integer value
 */
export function readInt(Module: LAPACKModule, ptr: number): number {
  return Module.HEAP32[ptr >> 2];
}

/**
 * Write doubles to WASM memory.
 * @param Module - LAPACK WASM module
 * @param ptr - Pointer to write to
 * @param values - Values to write
 */
export function writeDoubles(
  Module: LAPACKModule,
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
 * @param Module - LAPACK WASM module
 * @param ptr - Pointer to write to
 * @param value - Value to write
 */
export function writeInt(Module: LAPACKModule, ptr: number, value: number): void {
  Module.HEAP32[ptr >> 2] = value;
}

/**
 * Free multiple pointers at once.
 * Safely handles zero pointers (skips them).
 * @param Module - LAPACK WASM module
 * @param ptrs - Array of pointers to free
 */
export function freeAll(Module: LAPACKModule, ptrs: number[]): void {
  for (const ptr of ptrs) {
    if (ptr !== 0) {
      Module._free(ptr);
    }
  }
}

// ============================================================
// MATRIX CONVERSION HELPERS
// ============================================================

/**
 * Convert a row-major 2D array to column-major 1D array (Fortran order).
 * LAPACK expects column-major ordering.
 *
 * @param matrix - Row-major 2D array [row][col]
 * @returns Column-major 1D Float64Array
 */
export function toColumnMajor(matrix: number[][] | RealArray[]): Float64Array {
  const m = matrix.length;
  if (m === 0) return new Float64Array(0);
  const n = matrix[0].length;
  const result = new Float64Array(m * n);

  for (let j = 0; j < n; j++) {
    for (let i = 0; i < m; i++) {
      result[j * m + i] = matrix[i][j];
    }
  }
  return result;
}

/**
 * Convert a column-major 1D array to row-major 2D array.
 *
 * @param data - Column-major 1D array
 * @param m - Number of rows
 * @param n - Number of columns
 * @returns Row-major 2D array
 */
export function fromColumnMajor(data: Float64Array, m: number, n: number): number[][] {
  const result: number[][] = [];
  for (let i = 0; i < m; i++) {
    const row: number[] = [];
    for (let j = 0; j < n; j++) {
      row.push(data[j * m + i]);
    }
    result.push(row);
  }
  return result;
}

/**
 * Flatten a 2D array to 1D, assuming it's already in the desired order.
 *
 * @param matrix - 2D array
 * @returns Flattened 1D Float64Array
 */
export function flatten(matrix: number[][] | RealArray[]): Float64Array {
  const m = matrix.length;
  if (m === 0) return new Float64Array(0);
  const n = matrix[0].length;
  const result = new Float64Array(m * n);
  let idx = 0;
  for (let i = 0; i < m; i++) {
    for (let j = 0; j < n; j++) {
      result[idx++] = matrix[i][j];
    }
  }
  return result;
}

/**
 * Check if input is a 2D array (array of arrays).
 */
export function is2DArray(arr: unknown): arr is number[][] {
  return Array.isArray(arr) && arr.length > 0 && Array.isArray(arr[0]);
}

/**
 * Get matrix dimensions from various input formats.
 * Supports:
 * - 2D arrays: number[][]
 * - 1D arrays with explicit dimensions
 *
 * @param matrix - Input matrix
 * @param providedM - Provided number of rows (for 1D input)
 * @param providedN - Provided number of columns (for 1D input)
 * @returns [rows, cols]
 */
export function getMatrixDimensions(
  matrix: number[][] | RealArray,
  providedM?: number,
  providedN?: number
): [number, number] {
  if (is2DArray(matrix)) {
    return [matrix.length, matrix[0]?.length ?? 0];
  }

  // 1D array - need provided dimensions
  if (providedM !== undefined && providedN !== undefined) {
    return [providedM, providedN];
  }

  // Assume square matrix if only one dimension or none provided
  const len = matrix.length;
  const n = Math.floor(Math.sqrt(len));
  if (n * n !== len) {
    throw new Error(
      'Cannot determine matrix dimensions. Provide m and n for non-square matrices.'
    );
  }
  return [n, n];
}

/**
 * Prepare a matrix for LAPACK: convert to column-major Float64Array.
 * Handles both 2D arrays (row-major) and 1D arrays (assumed column-major).
 *
 * @param matrix - Input matrix (2D row-major or 1D column-major)
 * @returns Column-major Float64Array
 */
export function prepareMatrix(matrix: number[][] | RealArray): Float64Array {
  if (is2DArray(matrix)) {
    return toColumnMajor(matrix);
  }
  // Assume 1D is already column-major
  return matrix instanceof Float64Array ? matrix : new Float64Array(matrix);
}

/**
 * Prepare a vector for LAPACK.
 *
 * @param vector - Input vector
 * @returns Float64Array
 */
export function prepareVector(vector: RealArray): Float64Array {
  return vector instanceof Float64Array ? vector : new Float64Array(vector);
}
