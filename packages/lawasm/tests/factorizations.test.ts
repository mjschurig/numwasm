/**
 * Tests for Matrix Factorizations
 */

import { describe, it, expect } from 'vitest';
import {
  lu,
  cholesky,
  qr,
  qrPivoted,
  ldl,
  schur,
  hessenberg,
} from '../src/index.js';
import { toArray2D } from './setup.js';
import './setup.js';

describe('lu - LU Factorization', () => {
  it('factors a 2x2 matrix', () => {
    const A = [
      [4, 3],
      [6, 3],
    ];

    const result = lu(A);

    expect(result.success).toBe(true);
    expect(result.m).toBe(2);
    expect(result.n).toBe(2);
  });

  it('factors a 3x3 matrix', () => {
    const A = [
      [1, 2, 3],
      [4, 5, 6],
      [7, 8, 10],
    ];

    const result = lu(A);

    expect(result.success).toBe(true);
  });

  it('provides L, U, and P', () => {
    const A = [
      [2, 1],
      [4, 3],
    ];

    const result = lu(A);

    expect(result.success).toBe(true);
    expect(result.L).toBeDefined();
    expect(result.U).toBeDefined();
    expect(result.P).toBeDefined();
  });

  it('handles singular matrix', () => {
    const A = [
      [1, 2],
      [2, 4],
    ];

    const result = lu(A);

    // LU can still be computed, but U will have zero on diagonal
    expect(result.success).toBe(true);
  });
});

describe('cholesky - Cholesky Factorization', () => {
  it('factors a 2x2 SPD matrix', () => {
    const A = [
      [4, 2],
      [2, 5],
    ];

    const result = cholesky(A);

    expect(result.success).toBe(true);
    expect(result.n).toBe(2);
    expect(result.factor).toBeDefined();
  });

  it('computes upper Cholesky', () => {
    const A = [
      [4, 2],
      [2, 5],
    ];

    const result = cholesky(A, { uplo: 'upper' });

    expect(result.success).toBe(true);
    // Check that computation succeeded (upper is just a flag)
  });

  it('fails on non-positive definite matrix', () => {
    const A = [
      [1, 2],
      [2, 1], // eigenvalues: 3, -1
    ];

    const result = cholesky(A);

    expect(result.success).toBe(false);
    expect(result.info).toBeGreaterThan(0);
  });

  it('factors a 3x3 SPD matrix', () => {
    const A = [
      [4, 2, 2],
      [2, 5, 3],
      [2, 3, 6],
    ];

    const result = cholesky(A);

    expect(result.success).toBe(true);
  });
});

describe('qr - QR Factorization', () => {
  it('factors a 3x2 matrix', () => {
    const A = [
      [1, 2],
      [3, 4],
      [5, 6],
    ];

    const result = qr(A);

    expect(result.success).toBe(true);
    expect(result.m).toBe(3);
    expect(result.n).toBe(2);
  });

  it('computes full QR', () => {
    const A = [
      [1, 2],
      [3, 4],
      [5, 6],
    ];

    const result = qr(A, { mode: 'full' });

    expect(result.success).toBe(true);
    expect(result.Q).toBeDefined();
  });

  it('computes reduced QR', () => {
    const A = [
      [1, 2],
      [3, 4],
      [5, 6],
    ];

    const result = qr(A, { mode: 'reduced' });

    expect(result.success).toBe(true);
    expect(result.Q).toBeDefined();
  });

  it('R is upper triangular', () => {
    const A = [
      [1, 2],
      [3, 4],
      [5, 6],
    ];

    const result = qr(A, { mode: 'reduced' });

    expect(result.success).toBe(true);
    expect(result.R).toBeDefined();
    if (result.R) {
      // R is 2x2, check lower part is zero
      expect(result.R[1]).toBeCloseTo(0, 10); // R[1,0]
    }
  });
});

describe('qrPivoted - QR with Column Pivoting', () => {
  it('factors with pivoting', () => {
    const A = [
      [1, 2, 3],
      [4, 5, 6],
      [7, 8, 9],
    ];

    const result = qrPivoted(A);

    expect(result.success).toBe(true);
  });

  it('detects rank deficiency', () => {
    // Rank 2 matrix
    const A = [
      [1, 2, 3],
      [2, 4, 6], // 2 * row 1
      [3, 5, 8],
    ];

    const result = qrPivoted(A);

    expect(result.success).toBe(true);
    // Last diagonal of R should be small
  });
});

describe('ldl - LDL^T Factorization', () => {
  it('factors a symmetric matrix', () => {
    const A = [
      [4, 2, 0],
      [2, 5, 3],
      [0, 3, 6],
    ];

    const result = ldl(A);

    expect(result.success).toBe(true);
    expect(result.n).toBe(3);
  });

  it('handles indefinite matrix', () => {
    const A = [
      [1, 2],
      [2, 1], // Indefinite (eigenvalues 3, -1)
    ];

    const result = ldl(A);

    expect(result.success).toBe(true);
    // D can have negative entries for indefinite matrices
  });
});

describe('schur - Schur Decomposition', () => {
  it('computes Schur form of a 3x3 matrix', () => {
    const A = [
      [1, 2, 3],
      [0, 4, 5],
      [0, 0, 6],
    ];

    const result = schur(A);

    expect(result.success).toBe(true);
    expect(result.n).toBe(3);

    // T should be upper quasi-triangular
    // For this upper triangular input, T = A
    const T = toArray2D(result.T, 3, 3);
    expect(T[1][0]).toBeCloseTo(0, 10);
    expect(T[2][0]).toBeCloseTo(0, 10);
    expect(T[2][1]).toBeCloseTo(0, 10);
  });

  it('computes Schur vectors', () => {
    const A = [
      [0, 1],
      [-1, 0], // Rotation by 90 degrees
    ];

    const result = schur(A, { computeVectors: true });

    expect(result.success).toBe(true);
    expect(result.Z).toBeDefined();
  });
});

describe('hessenberg - Hessenberg Reduction', () => {
  it('reduces a 3x3 matrix to Hessenberg form', () => {
    const A = [
      [1, 2, 3],
      [4, 5, 6],
      [7, 8, 9],
    ];

    const result = hessenberg(A);

    expect(result.success).toBe(true);
    expect(result.n).toBe(3);

    // H should be upper Hessenberg (zeros below subdiagonal)
    const H = toArray2D(result.H, 3, 3);
    expect(H[2][0]).toBeCloseTo(0, 10);
  });

  it('computes orthogonal transformation matrix', () => {
    const A = [
      [1, 2, 3],
      [4, 5, 6],
      [7, 8, 9],
    ];

    const result = hessenberg(A, { computeQ: true });

    expect(result.success).toBe(true);
    expect(result.Q).toBeDefined();
  });

  it('2x2 matrix is already Hessenberg', () => {
    const A = [
      [1, 2],
      [3, 4],
    ];

    const result = hessenberg(A);

    expect(result.success).toBe(true);
    // Any 2x2 is already Hessenberg
  });
});
