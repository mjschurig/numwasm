/**
 * Tests for Matrix Inverses, Norms, and Related Functions
 */

import { describe, it, expect } from 'vitest';
import {
  inv,
  invTriangular,
  invSymmetric,
  pinv,
  norm,
  cond,
  condEst,
  rcond,
  det,
  logdet,
  slogdet,
  rank,
} from '../src/index.js';
import { toArray2D, matmulRef, eye } from './setup.js';
import './setup.js';

describe('inv - General Matrix Inverse', () => {
  it('inverts a 2x2 matrix', () => {
    const A = [
      [4, 7],
      [2, 6],
    ];
    // det = 24 - 14 = 10
    // inv = [6, -7; -2, 4] / 10

    const result = inv(A);

    expect(result.success).toBe(true);

    // Verify A * A^{-1} = I
    const Ainv = toArray2D(result.inv, 2, 2);
    const prod = matmulRef(A, Ainv);
    expect(prod[0][0]).toBeCloseTo(1, 10);
    expect(prod[0][1]).toBeCloseTo(0, 10);
    expect(prod[1][0]).toBeCloseTo(0, 10);
    expect(prod[1][1]).toBeCloseTo(1, 10);
  });

  it('inverts a 3x3 matrix', () => {
    const A = [
      [1, 2, 3],
      [0, 1, 4],
      [5, 6, 0],
    ];

    const result = inv(A);

    expect(result.success).toBe(true);

    // Verify
    const Ainv = toArray2D(result.inv, 3, 3);
    const prod = matmulRef(A, Ainv);
    for (let i = 0; i < 3; i++) {
      for (let j = 0; j < 3; j++) {
        expect(prod[i][j]).toBeCloseTo(i === j ? 1 : 0, 10);
      }
    }
  });

  it('fails on singular matrix', () => {
    const A = [
      [1, 2],
      [2, 4],
    ];

    const result = inv(A);

    expect(result.success).toBe(false);
  });
});

describe('invTriangular - Triangular Inverse', () => {
  it('inverts upper triangular matrix', () => {
    const U = [
      [2, 3],
      [0, 4],
    ];

    const result = invTriangular(U, { uplo: 'upper' });

    expect(result.success).toBe(true);

    // Verify
    const Uinv = toArray2D(result.inv, 2, 2);
    expect(Uinv[0][0]).toBeCloseTo(0.5, 10);
    expect(Uinv[1][1]).toBeCloseTo(0.25, 10);
    expect(Uinv[1][0]).toBeCloseTo(0, 10); // Should remain triangular
  });

  it('inverts lower triangular matrix', () => {
    const L = [
      [2, 0],
      [3, 4],
    ];

    const result = invTriangular(L, { uplo: 'lower' });

    expect(result.success).toBe(true);

    const Linv = toArray2D(result.inv, 2, 2);
    expect(Linv[0][1]).toBeCloseTo(0, 10); // Should remain triangular
  });
});

describe('invSymmetric - Symmetric Positive Definite Inverse', () => {
  it('inverts SPD matrix', () => {
    const A = [
      [4, 2],
      [2, 5],
    ];

    const result = invSymmetric(A);

    expect(result.success).toBe(true);

    // Verify
    const Ainv = toArray2D(result.inv, 2, 2);
    const prod = matmulRef(A, Ainv);
    expect(prod[0][0]).toBeCloseTo(1, 10);
    expect(prod[1][1]).toBeCloseTo(1, 10);
  });

  it('fails on non-positive definite', () => {
    const A = [
      [1, 2],
      [2, 1],
    ];

    const result = invSymmetric(A);

    expect(result.success).toBe(false);
  });
});

describe('pinv - Moore-Penrose Pseudoinverse', () => {
  it('computes pseudoinverse of full rank square matrix', () => {
    const A = [
      [1, 0],
      [0, 2],
    ];

    const result = pinv(A);

    expect(result.success).toBe(true);

    // For full rank, pinv(A) = inv(A)
    const Apinv = toArray2D(result.pinv, 2, 2);
    expect(Apinv[0][0]).toBeCloseTo(1, 10);
    expect(Apinv[1][1]).toBeCloseTo(0.5, 10);
  });

  it('computes pseudoinverse of tall matrix', () => {
    const A = [
      [1, 0],
      [0, 1],
      [0, 0],
    ];

    const result = pinv(A);

    expect(result.success).toBe(true);
    // Result should have pinv data - m and n refer to input dimensions
    expect(result.m).toBe(3);
    expect(result.n).toBe(2);
  });

  it('computes pseudoinverse of wide matrix', () => {
    const A = [
      [1, 0, 0],
      [0, 1, 0],
    ];

    const result = pinv(A);

    expect(result.success).toBe(true);
    // Result m and n refer to input dimensions
    expect(result.m).toBe(2);
    expect(result.n).toBe(3);
  });

  it('handles rank deficient matrix', () => {
    const A = [
      [1, 2],
      [2, 4], // rank 1
    ];

    const result = pinv(A);

    expect(result.success).toBe(true);
    expect(result.rank).toBe(1);
  });
});

describe('norm - Matrix Norms', () => {
  it('computes Frobenius norm', () => {
    const A = [
      [1, 2],
      [3, 4],
    ];
    // ||A||_F = sqrt(1 + 4 + 9 + 16) = sqrt(30)

    const result = norm(A, { ord: 'fro' });

    expect(result.success).toBe(true);
    expect(result.norm).toBeCloseTo(Math.sqrt(30), 10);
  });

  it('computes 1-norm (max column sum)', () => {
    const A = [
      [1, -2],
      [3, 4],
    ];
    // Column sums: |1|+|3|=4, |-2|+|4|=6

    const result = norm(A, { ord: '1' });

    expect(result.success).toBe(true);
    expect(result.norm).toBeCloseTo(6, 10);
  });

  it('computes infinity-norm (max row sum)', () => {
    const A = [
      [1, -2],
      [3, 4],
    ];
    // Row sums: |1|+|-2|=3, |3|+|4|=7

    const result = norm(A, { ord: 'inf' });

    expect(result.success).toBe(true);
    expect(result.norm).toBeCloseTo(7, 10);
  });

  it('computes 2-norm (largest singular value)', () => {
    const A = [
      [3, 0],
      [0, 4],
    ];

    const result = norm(A, { ord: '2' });

    expect(result.success).toBe(true);
    expect(result.norm).toBeCloseTo(4, 10);
  });
});

describe('cond - Condition Number', () => {
  it('computes condition number of well-conditioned matrix', () => {
    const A = [
      [1, 0],
      [0, 1],
    ];

    const result = cond(A);

    expect(result.success).toBe(true);
    expect(result.cond).toBeCloseTo(1, 10);
  });

  it('computes condition number of ill-conditioned matrix', () => {
    const A = [
      [1, 0],
      [0, 1e-10],
    ];

    const result = cond(A);

    expect(result.success).toBe(true);
    expect(result.cond).toBeCloseTo(1e10, 5);
  });
});

describe('condEst - Condition Number Estimate', () => {
  it('estimates condition number', () => {
    const A = [
      [1, 0],
      [0, 2],
    ];

    const result = condEst(A);

    expect(result.success).toBe(true);
    // Estimate should be close to actual
    expect(result.cond).toBeDefined();
  });
});

describe('rcond - Reciprocal Condition Number', () => {
  it('computes rcond for well-conditioned matrix', () => {
    const A = [
      [1, 0],
      [0, 1],
    ];

    const result = rcond(A);

    expect(result.success).toBe(true);
    expect(result.rcond).toBeCloseTo(1, 5);
  });

  it('rcond is small for ill-conditioned matrix', () => {
    const A = [
      [1, 0],
      [0, 1e-10],
    ];

    const result = rcond(A);

    expect(result.success).toBe(true);
    expect(result.rcond).toBeLessThan(1e-8);
  });
});

describe('det - Determinant', () => {
  it('computes determinant of 2x2 matrix', () => {
    const A = [
      [4, 7],
      [2, 6],
    ];
    // det = 4*6 - 7*2 = 24 - 14 = 10

    const result = det(A);

    expect(result.success).toBe(true);
    expect(result.det).toBeCloseTo(10, 10);
  });

  it('computes determinant of 3x3 matrix', () => {
    const A = [
      [1, 2, 3],
      [4, 5, 6],
      [7, 8, 10], // det = -3
    ];

    const result = det(A);

    expect(result.success).toBe(true);
    expect(result.det).toBeCloseTo(-3, 10);
  });

  it('determinant of singular matrix is zero', () => {
    const A = [
      [1, 2],
      [2, 4],
    ];

    const result = det(A);

    expect(result.success).toBe(true);
    expect(result.det).toBeCloseTo(0, 10);
  });

  it('determinant of identity is 1', () => {
    const I = eye(3);

    const result = det(I);

    expect(result.success).toBe(true);
    expect(result.det).toBeCloseTo(1, 10);
  });
});

describe('logdet - Log Determinant', () => {
  it('computes log determinant of positive definite matrix', () => {
    const A = [
      [Math.E, 0],
      [0, Math.E * Math.E],
    ];
    // det = e * e^2 = e^3
    // log(det) = 3

    const result = logdet(A);

    expect(result.success).toBe(true);
    expect(result.logdet).toBeCloseTo(3, 10);
  });
});

describe('slogdet - Sign and Log Determinant', () => {
  it('computes sign and log for positive determinant', () => {
    const A = [
      [2, 0],
      [0, 3],
    ];

    const result = slogdet(A);

    expect(result.success).toBe(true);
    expect(result.sign).toBe(1);
    expect(result.logabsdet).toBeCloseTo(Math.log(6), 10);
  });

  it('computes sign and log for negative determinant', () => {
    const A = [
      [0, 1],
      [1, 0],
    ];
    // Permutation matrix, det = -1

    const result = slogdet(A);

    expect(result.success).toBe(true);
    expect(result.sign).toBe(-1);
    expect(result.logabsdet).toBeCloseTo(0, 10); // log(1) = 0
  });
});

describe('rank - Matrix Rank', () => {
  it('computes rank of full rank matrix', () => {
    const A = [
      [1, 0],
      [0, 1],
    ];

    const result = rank(A);

    expect(result.success).toBe(true);
    expect(result.rank).toBe(2);
  });

  it('computes rank of rank-deficient matrix', () => {
    const A = [
      [1, 2],
      [2, 4],
    ];

    const result = rank(A);

    expect(result.success).toBe(true);
    expect(result.rank).toBe(1);
  });

  it('computes rank of rectangular matrix', () => {
    const A = [
      [1, 0, 0],
      [0, 1, 0],
    ];

    const result = rank(A);

    expect(result.success).toBe(true);
    expect(result.rank).toBe(2);
  });

  it('rank of zero matrix is 0', () => {
    const A = [
      [0, 0],
      [0, 0],
    ];

    const result = rank(A);

    expect(result.success).toBe(true);
    expect(result.rank).toBe(0);
  });
});
