/**
 * Test Setup for LAWasm
 *
 * Loads the LAPACK module before all tests.
 */

import { loadLAPACKModule, isLAPACKLoaded } from '../src/index.js';
import { beforeAll, expect } from 'vitest';

// Load module once before all tests
beforeAll(async () => {
  if (!isLAPACKLoaded()) {
    await loadLAPACKModule();
  }
});

/**
 * Helper to check if two arrays are approximately equal.
 */
export function expectArrayClose(
  actual: ArrayLike<number>,
  expected: ArrayLike<number>,
  tol = 1e-10
): void {
  expect(actual.length).toBe(expected.length);
  for (let i = 0; i < actual.length; i++) {
    expect(actual[i]).toBeCloseTo(expected[i], -Math.log10(tol));
  }
}

/**
 * Helper to check if two matrices are approximately equal.
 */
export function expectMatrixClose(
  actual: Float64Array,
  expected: number[][],
  m: number,
  n: number,
  tol = 1e-10
): void {
  for (let i = 0; i < m; i++) {
    for (let j = 0; j < n; j++) {
      // Column-major indexing
      expect(actual[i + j * m]).toBeCloseTo(expected[i][j], -Math.log10(tol));
    }
  }
}

/**
 * Helper to multiply two matrices (for verification).
 */
export function matmulRef(
  A: number[][],
  B: number[][]
): number[][] {
  const m = A.length;
  const n = B[0].length;
  const k = B.length;
  const C: number[][] = [];
  for (let i = 0; i < m; i++) {
    C.push([]);
    for (let j = 0; j < n; j++) {
      let sum = 0;
      for (let l = 0; l < k; l++) {
        sum += A[i][l] * B[l][j];
      }
      C[i].push(sum);
    }
  }
  return C;
}

/**
 * Helper to transpose a matrix.
 */
export function transposeRef(A: number[][]): number[][] {
  const m = A.length;
  const n = A[0].length;
  const At: number[][] = [];
  for (let j = 0; j < n; j++) {
    At.push([]);
    for (let i = 0; i < m; i++) {
      At[j].push(A[i][j]);
    }
  }
  return At;
}

/**
 * Helper to create an identity matrix.
 */
export function eye(n: number): number[][] {
  const I: number[][] = [];
  for (let i = 0; i < n; i++) {
    I.push([]);
    for (let j = 0; j < n; j++) {
      I[i].push(i === j ? 1 : 0);
    }
  }
  return I;
}

/**
 * Helper to create a zero matrix.
 */
export function zeros(m: number, n: number): number[][] {
  const Z: number[][] = [];
  for (let i = 0; i < m; i++) {
    Z.push(new Array(n).fill(0));
  }
  return Z;
}

/**
 * Convert column-major Float64Array to 2D array.
 */
export function toArray2D(data: Float64Array, m: number, n: number): number[][] {
  const result: number[][] = [];
  for (let i = 0; i < m; i++) {
    result.push([]);
    for (let j = 0; j < n; j++) {
      result[i].push(data[i + j * m]);
    }
  }
  return result;
}
