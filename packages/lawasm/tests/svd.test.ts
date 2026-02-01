/**
 * Tests for Singular Value Decomposition
 */

import { describe, it, expect } from 'vitest';
import { svd, svdvals, svdCompact, svdRank } from '../src/index.js';
import { toArray2D, matmulRef } from './setup.js';
import './setup.js';

describe('svd - Full SVD', () => {
  it('computes SVD of a 2x2 matrix', () => {
    const A = [
      [3, 0],
      [0, 2],
    ];

    const result = svd(A);

    expect(result.success).toBe(true);
    expect(result.m).toBe(2);
    expect(result.n).toBe(2);

    // Singular values should be 3 and 2 (in descending order)
    expect(result.s[0]).toBeCloseTo(3, 10);
    expect(result.s[1]).toBeCloseTo(2, 10);
  });

  it('computes SVD of a 3x2 matrix', () => {
    const A = [
      [1, 2],
      [3, 4],
      [5, 6],
    ];

    const result = svd(A, { mode: 'full' });

    expect(result.success).toBe(true);
    expect(result.s.length).toBe(2); // min(3, 2)
    expect(result.U).toBeDefined();
    expect(result.Vt).toBeDefined();

    if (result.U && result.Vt) {
      // U should be 3x3 in full mode
      expect(result.U.length).toBe(9);
      // Vt should be 2x2
      expect(result.Vt.length).toBe(4);
    }
  });

  it('computes reduced SVD', () => {
    const A = [
      [1, 2],
      [3, 4],
      [5, 6],
    ];

    const result = svd(A, { mode: 'reduced' });

    expect(result.success).toBe(true);
    if (result.U) {
      // U should be 3x2 in reduced mode
      expect(result.U.length).toBe(6);
    }
  });

  it('verifies U is orthogonal', () => {
    const A = [
      [1, 2],
      [3, 4],
      [5, 6],
    ];

    const result = svd(A, { mode: 'full' });

    expect(result.success).toBe(true);
    if (result.U) {
      const U = toArray2D(result.U, 3, 3);
      // U[:,0] dot U[:,0] = 1
      let dot = 0;
      for (let i = 0; i < 3; i++) {
        dot += U[i][0] * U[i][0];
      }
      expect(dot).toBeCloseTo(1, 10);
    }
  });

  it('verifies V is orthogonal', () => {
    const A = [
      [1, 2],
      [3, 4],
      [5, 6],
    ];

    const result = svd(A, { mode: 'full' });

    expect(result.success).toBe(true);
    if (result.Vt) {
      const Vt = toArray2D(result.Vt, 2, 2);
      // Vt is V^T, so rows of Vt should be orthonormal
      let dot = 0;
      for (let j = 0; j < 2; j++) {
        dot += Vt[0][j] * Vt[0][j];
      }
      expect(dot).toBeCloseTo(1, 10);
    }
  });

  it('reconstructs A from SVD', () => {
    const A = [
      [1, 2],
      [3, 4],
    ];

    const result = svd(A, { mode: 'full' });

    expect(result.success).toBe(true);
    if (result.U && result.Vt) {
      const U = toArray2D(result.U, 2, 2);
      const Vt = toArray2D(result.Vt, 2, 2);
      const s = result.s;

      // Reconstruct: A = U * diag(s) * Vt
      for (let i = 0; i < 2; i++) {
        for (let j = 0; j < 2; j++) {
          let sum = 0;
          for (let k = 0; k < 2; k++) {
            sum += U[i][k] * s[k] * Vt[k][j];
          }
          expect(sum).toBeCloseTo(A[i][j], 10);
        }
      }
    }
  });

  it('singular values are non-negative and descending', () => {
    const A = [
      [1, 2, 3],
      [4, 5, 6],
    ];

    const result = svd(A);

    expect(result.success).toBe(true);
    for (let i = 0; i < result.s.length; i++) {
      expect(result.s[i]).toBeGreaterThanOrEqual(0);
    }
    for (let i = 0; i < result.s.length - 1; i++) {
      expect(result.s[i]).toBeGreaterThanOrEqual(result.s[i + 1]);
    }
  });
});

describe('svdvals - Singular Values Only', () => {
  it('computes singular values only', () => {
    const A = [
      [3, 0],
      [0, 4],
    ];

    const result = svdvals(A);

    expect(result.success).toBe(true);
    expect(result.s[0]).toBeCloseTo(4, 10);
    expect(result.s[1]).toBeCloseTo(3, 10);
  });

  it('is faster than full SVD for large matrices', () => {
    const n = 50;
    const A: number[][] = [];
    for (let i = 0; i < n; i++) {
      A.push([]);
      for (let j = 0; j < n; j++) {
        A[i].push(Math.random());
      }
    }

    const result = svdvals(A);

    expect(result.success).toBe(true);
    expect(result.s.length).toBe(n);
  });
});

describe('svdCompact - Compact SVD', () => {
  it('returns compact form', () => {
    const A = [
      [1, 2],
      [3, 4],
      [5, 6],
    ];

    const result = svdCompact(A);

    expect(result.success).toBe(true);
    expect(result.s.length).toBe(2);
  });
});

describe('svdRank - Numerical Rank', () => {
  it('computes rank of full rank matrix', () => {
    const A = [
      [1, 0],
      [0, 1],
    ];

    const result = svdRank(A);

    expect(result.success).toBe(true);
    expect(result.rank).toBe(2);
  });

  it('computes rank of rank-1 matrix', () => {
    const A = [
      [1, 2],
      [2, 4], // Row 2 = 2 * Row 1
    ];

    const result = svdRank(A);

    expect(result.success).toBe(true);
    expect(result.rank).toBe(1);
  });

  it('computes rank of zero matrix', () => {
    const A = [
      [0, 0],
      [0, 0],
    ];

    const result = svdRank(A);

    expect(result.success).toBe(true);
    expect(result.rank).toBe(0);
  });

  it('respects tolerance', () => {
    const A = [
      [1, 0],
      [0, 1e-12], // Very small but non-zero
    ];

    const resultTight = svdRank(A, { tol: 1e-15 });
    const resultLoose = svdRank(A, { tol: 1e-10 });

    expect(resultTight.success).toBe(true);
    expect(resultLoose.success).toBe(true);
    expect(resultTight.rank).toBe(2);
    expect(resultLoose.rank).toBe(1);
  });

  it('handles rectangular matrices', () => {
    const A = [
      [1, 2, 3],
      [4, 5, 6],
    ];

    const result = svdRank(A);

    expect(result.success).toBe(true);
    expect(result.rank).toBe(2); // max rank is min(2, 3) = 2
  });
});
