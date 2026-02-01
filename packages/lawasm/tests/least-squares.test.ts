/**
 * Tests for Least Squares & Minimum Norm
 */

import { describe, it, expect } from 'vitest';
import {
  lstsq,
  lstsqSVD,
  lstsqGelsy,
  generalizedLstSq,
} from '../src/index.js';
import './setup.js';

describe('lstsq - Least Squares via QR', () => {
  it('solves overdetermined system', () => {
    // 3 equations, 2 unknowns
    const A = [
      [1, 1],
      [1, 2],
      [1, 3],
    ];
    const b = [1, 2, 2];

    const result = lstsq(A, b);

    expect(result.success).toBe(true);
    expect(result.x.length).toBe(2);
  });

  it('solves underdetermined system (minimum norm)', () => {
    // 2 equations, 3 unknowns
    const A = [
      [1, 0, 1],
      [0, 1, 1],
    ];
    const b = [1, 1];

    const result = lstsq(A, b);

    expect(result.success).toBe(true);
    expect(result.x.length).toBe(3);
  });

  it('exact solution for square system', () => {
    const A = [
      [2, 1],
      [1, 3],
    ];
    const b = [4, 5];

    const result = lstsq(A, b);

    expect(result.success).toBe(true);
    // Verify A * x = b
    const Ax0 = A[0][0] * result.x[0] + A[0][1] * result.x[1];
    const Ax1 = A[1][0] * result.x[0] + A[1][1] * result.x[1];
    expect(Ax0).toBeCloseTo(b[0], 10);
    expect(Ax1).toBeCloseTo(b[1], 10);
  });

  it('minimizes residual for overdetermined system', () => {
    const A = [
      [1, 0],
      [0, 1],
      [1, 1],
    ];
    const b = [1, 1, 3]; // Inconsistent

    const result = lstsq(A, b);

    expect(result.success).toBe(true);
    // Residual should be non-zero
    const r0 = b[0] - (A[0][0] * result.x[0] + A[0][1] * result.x[1]);
    const r1 = b[1] - (A[1][0] * result.x[0] + A[1][1] * result.x[1]);
    const r2 = b[2] - (A[2][0] * result.x[0] + A[2][1] * result.x[1]);
    const residualNorm = Math.sqrt(r0 * r0 + r1 * r1 + r2 * r2);
    expect(residualNorm).toBeGreaterThan(0);
  });
});

describe('lstsqSVD - Least Squares via SVD', () => {
  it('handles rank-deficient matrix', () => {
    // Rank 1 matrix
    const A = [
      [1, 2],
      [2, 4], // 2 * row 1
      [3, 6], // 3 * row 1
    ];
    const b = [1, 2, 3];

    const result = lstsqSVD(A, b);

    expect(result.success).toBe(true);
    expect(result.rank).toBe(1);
  });

  it('computes rank correctly', () => {
    const A = [
      [1, 0],
      [0, 1],
      [0, 0],
    ];
    const b = [1, 1, 0];

    const result = lstsqSVD(A, b);

    expect(result.success).toBe(true);
    expect(result.rank).toBe(2);
  });

  it('respects rcond parameter', () => {
    const A = [
      [1, 0],
      [0, 1e-10],
    ];
    const b = [1, 1e-10];

    const resultTight = lstsqSVD(A, b, { rcond: 1e-15 });
    const resultLoose = lstsqSVD(A, b, { rcond: 1e-5 });

    expect(resultTight.success).toBe(true);
    expect(resultLoose.success).toBe(true);
    expect(resultTight.rank).toBe(2);
    expect(resultLoose.rank).toBe(1);
  });

  it('returns singular values', () => {
    const A = [
      [3, 0],
      [0, 4],
    ];
    const b = [3, 4];

    const result = lstsqSVD(A, b);

    expect(result.success).toBe(true);
    expect(result.s).toBeDefined();
    if (result.s) {
      expect(result.s[0]).toBeCloseTo(4, 10);
      expect(result.s[1]).toBeCloseTo(3, 10);
    }
  });
});

describe('lstsqGelsy - Least Squares via QR with Pivoting', () => {
  it('handles rank-deficient matrix', () => {
    const A = [
      [1, 2, 3],
      [2, 4, 6], // 2 * row 1
    ];
    const b = [6, 12];

    const result = lstsqGelsy(A, b);

    expect(result.success).toBe(true);
    expect(result.rank).toBeLessThanOrEqual(2);
  });

  it('provides solution for full rank', () => {
    const A = [
      [1, 0],
      [0, 1],
    ];
    const b = [1, 1];

    const result = lstsqGelsy(A, b);

    expect(result.success).toBe(true);
    expect(result.x[0]).toBeCloseTo(1, 10);
    expect(result.x[1]).toBeCloseTo(1, 10);
  });
});

describe('generalizedLstSq - Generalized Least Squares', () => {
  it('solves generalized problem', () => {
    // min ||y||² subject to Ax + By = d
    const A = [
      [1, 0],
      [0, 1],
    ];
    const B = [
      [1, 0],
      [0, 1],
    ];
    const d = [1, 1];

    const result = generalizedLstSq(A, B, d);

    expect(result.success).toBe(true);
    expect(result.x.length).toBe(2);
    expect(result.y.length).toBe(2);

    // Check Ax + By = d
    for (let i = 0; i < 2; i++) {
      let sum = 0;
      for (let j = 0; j < 2; j++) {
        sum += A[i][j] * result.x[j] + B[i][j] * result.y[j];
      }
      expect(sum).toBeCloseTo(d[i], 10);
    }
  });
});
