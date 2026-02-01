/**
 * Tests for Linear System Solvers
 */

import { describe, it, expect } from 'vitest';
import {
  solve,
  solveTriangular,
  solveSymmetric,
  solveTridiagonal,
} from '../src/index.js';
import { expectArrayClose } from './setup.js';
import './setup.js';

describe('solve - General Linear System', () => {
  it('solves a simple 2x2 system', () => {
    const A = [
      [2, 1],
      [1, 3],
    ];
    const b = [1, 2];

    const result = solve(A, b);

    expect(result.success).toBe(true);
    // Verify: A * x = b
    const Ax = [
      A[0][0] * result.x[0] + A[0][1] * result.x[1],
      A[1][0] * result.x[0] + A[1][1] * result.x[1],
    ];
    expectArrayClose(Ax, b);
  });

  it('solves a 3x3 system', () => {
    const A = [
      [1, 2, 3],
      [4, 5, 6],
      [7, 8, 10], // Not 9 to avoid singularity
    ];
    const b = [1, 2, 3];

    const result = solve(A, b);

    expect(result.success).toBe(true);
    // Verify A * x = b
    for (let i = 0; i < 3; i++) {
      let sum = 0;
      for (let j = 0; j < 3; j++) {
        sum += A[i][j] * result.x[j];
      }
      expect(sum).toBeCloseTo(b[i], 10);
    }
  });

  it('solves with multiple right-hand sides', () => {
    const A = [
      [2, 1],
      [1, 3],
    ];
    const B = [
      [1, 4],
      [2, 5],
    ];

    const result = solve(A, B);

    expect(result.success).toBe(true);
    expect(result.nrhs).toBe(2);
  });

  it('handles singular matrix gracefully', () => {
    const A = [
      [1, 2],
      [2, 4], // Row 2 is 2 * Row 1
    ];
    const b = [1, 2];

    const result = solve(A, b);

    // Should fail or return info > 0
    expect(result.info).toBeGreaterThan(0);
  });
});

describe('solveTriangular - Triangular System', () => {
  it('solves upper triangular system', () => {
    const U = [
      [2, 1],
      [0, 3],
    ];
    const b = [4, 3];

    const result = solveTriangular(U, b, { uplo: 'upper' });

    expect(result.success).toBe(true);
    // Verify: U * x = b
    // x[1] = 3/3 = 1
    // x[0] = (4 - 1*1)/2 = 1.5
    expect(result.x[1]).toBeCloseTo(1, 10);
    expect(result.x[0]).toBeCloseTo(1.5, 10);
  });

  it('solves lower triangular system', () => {
    const L = [
      [2, 0],
      [1, 3],
    ];
    const b = [2, 5];

    const result = solveTriangular(L, b, { uplo: 'lower' });

    expect(result.success).toBe(true);
    // Just verify it returns a solution
    expect(result.x.length).toBe(2);
  });
});

describe('solveSymmetric - Symmetric Positive Definite', () => {
  it('solves a symmetric positive definite system', () => {
    // SPD matrix
    const A = [
      [4, 2],
      [2, 5],
    ];
    const b = [6, 9];

    const result = solveSymmetric(A, b);

    expect(result.success).toBe(true);
    // Verify A * x = b
    const Ax = [
      A[0][0] * result.x[0] + A[0][1] * result.x[1],
      A[1][0] * result.x[0] + A[1][1] * result.x[1],
    ];
    expectArrayClose(Ax, b);
  });

  it('solves a 3x3 SPD system', () => {
    // Create SPD: A = L * L^T
    const A = [
      [4, 2, 2],
      [2, 5, 3],
      [2, 3, 6],
    ];
    const b = [8, 10, 11];

    const result = solveSymmetric(A, b);

    expect(result.success).toBe(true);
    // Verify
    for (let i = 0; i < 3; i++) {
      let sum = 0;
      for (let j = 0; j < 3; j++) {
        sum += A[i][j] * result.x[j];
      }
      expect(sum).toBeCloseTo(b[i], 10);
    }
  });

  it('fails on non-positive definite matrix', () => {
    const A = [
      [1, 2],
      [2, 1], // Not positive definite
    ];
    const b = [1, 1];

    const result = solveSymmetric(A, b);

    expect(result.success).toBe(false);
  });
});

describe('solveTridiagonal - Tridiagonal System', () => {
  it('solves a simple tridiagonal system', () => {
    // Matrix:
    // [2  -1   0]
    // [-1  2  -1]
    // [0  -1   2]
    const dl = [-1, -1]; // subdiagonal
    const d = [2, 2, 2]; // diagonal
    const du = [-1, -1]; // superdiagonal
    const b = [1, 0, 1];

    const result = solveTridiagonal(dl, d, du, b);

    expect(result.success).toBe(true);
    // This is a discretized Laplacian, solution should be symmetric
    expect(result.x[0]).toBeCloseTo(result.x[2], 10);
  });
});
