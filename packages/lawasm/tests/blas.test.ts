/**
 * Tests for BLAS Operations
 */

import { describe, it, expect } from 'vitest';
import {
  // Level 3
  matmul,
  matmulTriangular,
  solveMatrixTriangular,
  syrk,
  syr2k,
  // Level 2
  matvec,
  matvecTriangular,
  solveVectorTriangular,
  ger,
  syr,
  // Level 1
  dot,
  axpy,
  scal,
  copy,
  swap,
  nrm2,
  asum,
  iamax,
} from '../src/index.js';
import { toArray2D, matmulRef } from './setup.js';
import './setup.js';

describe('BLAS Level 3 - Matrix-Matrix Operations', () => {
  describe('matmul - General Matrix Multiplication', () => {
    it('multiplies two 2x2 matrices', () => {
      const A = [
        [1, 2],
        [3, 4],
      ];
      const B = [
        [5, 6],
        [7, 8],
      ];

      const result = matmul(A, B);

      expect(result.success).toBe(true);
      const C = toArray2D(result.C, 2, 2);
      // C = A * B
      expect(C[0][0]).toBeCloseTo(19, 10); // 1*5 + 2*7
      expect(C[0][1]).toBeCloseTo(22, 10); // 1*6 + 2*8
      expect(C[1][0]).toBeCloseTo(43, 10); // 3*5 + 4*7
      expect(C[1][1]).toBeCloseTo(50, 10); // 3*6 + 4*8
    });

    it('multiplies with alpha and beta', () => {
      const A = [
        [1, 0],
        [0, 1],
      ];
      const B = [
        [2, 0],
        [0, 2],
      ];
      const C = [
        [1, 1],
        [1, 1],
      ];

      const result = matmul(A, B, { alpha: 2, beta: 3, C });

      expect(result.success).toBe(true);
      // C = 2*A*B + 3*C = 2*[[2,0],[0,2]] + 3*[[1,1],[1,1]]
      //   = [[4,0],[0,4]] + [[3,3],[3,3]] = [[7,3],[3,7]]
      const Cout = toArray2D(result.C, 2, 2);
      expect(Cout[0][0]).toBeCloseTo(7, 10);
      expect(Cout[0][1]).toBeCloseTo(3, 10);
    });

    it('handles transpose operations', () => {
      const A = [
        [1, 2],
        [3, 4],
      ];
      const B = [
        [1, 0],
        [0, 1],
      ];

      const result = matmul(A, B, { transA: 'T' });

      expect(result.success).toBe(true);
      // C = A^T * B = [[1,3],[2,4]] * I = [[1,3],[2,4]]
      const C = toArray2D(result.C, 2, 2);
      expect(C[0][0]).toBeCloseTo(1, 10);
      expect(C[1][0]).toBeCloseTo(2, 10);
    });

    it('multiplies non-square matrices', () => {
      const A = [
        [1, 2, 3],
        [4, 5, 6],
      ]; // 2x3
      const B = [
        [1, 2],
        [3, 4],
        [5, 6],
      ]; // 3x2

      const result = matmul(A, B);

      expect(result.success).toBe(true);
      expect(result.m).toBe(2);
      expect(result.n).toBe(2);
    });
  });

  describe('matmulTriangular - Triangular Matrix Multiply', () => {
    it('multiplies with upper triangular matrix', () => {
      const A = [
        [2, 1],
        [0, 3],
      ];
      const B = [
        [1, 2],
        [3, 4],
      ];

      const result = matmulTriangular(A, B, { side: 'left', uplo: 'upper' });

      expect(result.success).toBe(true);
      // B = A * B (left multiply)
    });
  });

  describe('solveMatrixTriangular - Triangular Matrix Solve', () => {
    it('solves AX = B with upper triangular A', () => {
      const A = [
        [2, 1],
        [0, 3],
      ];
      const B = [
        [4, 5],
        [3, 6],
      ];

      const result = solveMatrixTriangular(A, B, { side: 'left', uplo: 'upper' });

      expect(result.success).toBe(true);
      // Verify A * X = B
    });
  });

  describe('syrk - Symmetric Rank-k Update', () => {
    it('computes C = A * A^T', () => {
      const A = [
        [1, 2],
        [3, 4],
      ];

      const result = syrk(A);

      expect(result.success).toBe(true);
      // C = A * A^T
      // [[1,2],[3,4]] * [[1,3],[2,4]] = [[5,11],[11,25]]
      const C = toArray2D(result.C, 2, 2);
      expect(C[0][0]).toBeCloseTo(5, 10);
      expect(C[1][1]).toBeCloseTo(25, 10);
    });
  });

  describe('syr2k - Symmetric Rank-2k Update', () => {
    it('computes C = A * B^T + B * A^T', () => {
      const A = [
        [1, 0],
        [0, 1],
      ];
      const B = [
        [1, 0],
        [0, 1],
      ];

      const result = syr2k(A, B);

      expect(result.success).toBe(true);
      // C = I * I^T + I * I^T = 2*I
      const C = toArray2D(result.C, 2, 2);
      expect(C[0][0]).toBeCloseTo(2, 10);
      expect(C[1][1]).toBeCloseTo(2, 10);
    });
  });
});

describe('BLAS Level 2 - Matrix-Vector Operations', () => {
  describe('matvec - General Matrix-Vector Multiply', () => {
    it('computes y = A * x', () => {
      const A = [
        [1, 2],
        [3, 4],
      ];
      const x = [1, 2];

      const result = matvec(A, x);

      expect(result.success).toBe(true);
      expect(result.y[0]).toBeCloseTo(5, 10); // 1*1 + 2*2
      expect(result.y[1]).toBeCloseTo(11, 10); // 3*1 + 4*2
    });

    it('computes y = alpha * A * x + beta * y', () => {
      const A = [
        [1, 0],
        [0, 1],
      ];
      const x = [2, 3];
      const y = [1, 1];

      const result = matvec(A, x, { alpha: 2, beta: 1, y });

      expect(result.success).toBe(true);
      // y = 2 * I * [2,3] + 1 * [1,1] = [4,6] + [1,1] = [5,7]
      expect(result.y[0]).toBeCloseTo(5, 10);
      expect(result.y[1]).toBeCloseTo(7, 10);
    });

    it('handles transpose', () => {
      const A = [
        [1, 2],
        [3, 4],
      ];
      const x = [1, 1];

      const result = matvec(A, x, { trans: 'T' });

      expect(result.success).toBe(true);
      // y = A^T * x = [[1,3],[2,4]] * [1,1] = [4, 6]
      expect(result.y[0]).toBeCloseTo(4, 10);
      expect(result.y[1]).toBeCloseTo(6, 10);
    });
  });

  describe('matvecTriangular - Triangular Matrix-Vector Multiply', () => {
    it('computes x = A * x for upper triangular A', () => {
      const A = [
        [2, 1],
        [0, 3],
      ];
      const x = [1, 2];

      const result = matvecTriangular(A, x, { uplo: 'upper' });

      expect(result.success).toBe(true);
      // x = [[2,1],[0,3]] * [1,2] = [4, 6]
      expect(result.x[0]).toBeCloseTo(4, 10);
      expect(result.x[1]).toBeCloseTo(6, 10);
    });
  });

  describe('solveVectorTriangular - Triangular Vector Solve', () => {
    it('solves A * x = b for upper triangular A', () => {
      const A = [
        [2, 1],
        [0, 3],
      ];
      const b = [4, 6];

      const result = solveVectorTriangular(A, b, { uplo: 'upper' });

      expect(result.success).toBe(true);
      // x[1] = 6/3 = 2
      // x[0] = (4 - 1*2)/2 = 1
      expect(result.x[1]).toBeCloseTo(2, 10);
      expect(result.x[0]).toBeCloseTo(1, 10);
    });
  });

  describe('ger - Rank-1 Update', () => {
    it('computes A = alpha * x * y^T + A', () => {
      const x = [1, 2];
      const y = [3, 4];

      const result = ger(x, y);

      expect(result.success).toBe(true);
      // A = x * y^T = [[3,4],[6,8]]
      const A = toArray2D(result.A, 2, 2);
      expect(A[0][0]).toBeCloseTo(3, 10);
      expect(A[0][1]).toBeCloseTo(4, 10);
      expect(A[1][0]).toBeCloseTo(6, 10);
      expect(A[1][1]).toBeCloseTo(8, 10);
    });
  });

  describe('syr - Symmetric Rank-1 Update', () => {
    it('computes A = alpha * x * x^T + A', () => {
      const x = [1, 2];

      const result = syr(x);

      expect(result.success).toBe(true);
      // A = x * x^T = [[1,2],[2,4]]
      const A = toArray2D(result.A, 2, 2);
      expect(A[0][0]).toBeCloseTo(1, 10);
      expect(A[1][1]).toBeCloseTo(4, 10);
    });
  });
});

describe('BLAS Level 1 - Vector Operations', () => {
  describe('dot - Dot Product', () => {
    it('computes dot product', () => {
      const x = [1, 2, 3];
      const y = [4, 5, 6];

      const result = dot(x, y);

      expect(result.success).toBe(true);
      expect(result.dot).toBeCloseTo(32, 10); // 1*4 + 2*5 + 3*6
    });

    it('handles orthogonal vectors', () => {
      const x = [1, 0];
      const y = [0, 1];

      const result = dot(x, y);

      expect(result.success).toBe(true);
      expect(result.dot).toBeCloseTo(0, 10);
    });
  });

  describe('axpy - y = alpha*x + y', () => {
    it('computes y = 2*x + y', () => {
      const x = [1, 2, 3];
      const y = [4, 5, 6];

      const result = axpy(2, x, y);

      expect(result.success).toBe(true);
      expect(result.y[0]).toBeCloseTo(6, 10); // 2*1 + 4
      expect(result.y[1]).toBeCloseTo(9, 10); // 2*2 + 5
      expect(result.y[2]).toBeCloseTo(12, 10); // 2*3 + 6
    });
  });

  describe('scal - x = alpha*x', () => {
    it('scales vector', () => {
      const x = [1, 2, 3];

      const result = scal(3, x);

      expect(result.success).toBe(true);
      expect(result.x[0]).toBeCloseTo(3, 10);
      expect(result.x[1]).toBeCloseTo(6, 10);
      expect(result.x[2]).toBeCloseTo(9, 10);
    });
  });

  describe('copy - y = x', () => {
    it('copies vector', () => {
      const x = [1, 2, 3];
      const y = [0, 0, 0];

      const result = copy(x, y);

      expect(result.success).toBe(true);
      expect(result.y[0]).toBeCloseTo(1, 10);
      expect(result.y[1]).toBeCloseTo(2, 10);
      expect(result.y[2]).toBeCloseTo(3, 10);
    });
  });

  describe('swap - swap x and y', () => {
    it('swaps vectors', () => {
      const x = [1, 2, 3];
      const y = [4, 5, 6];

      const result = swap(x, y);

      expect(result.success).toBe(true);
      expect(result.x[0]).toBeCloseTo(4, 10);
      expect(result.y[0]).toBeCloseTo(1, 10);
    });
  });

  describe('nrm2 - Euclidean Norm', () => {
    it('computes 2-norm', () => {
      const x = [3, 4];

      const result = nrm2(x);

      expect(result.success).toBe(true);
      expect(result.norm).toBeCloseTo(5, 10);
    });

    it('handles zero vector', () => {
      const x = [0, 0, 0];

      const result = nrm2(x);

      expect(result.success).toBe(true);
      expect(result.norm).toBeCloseTo(0, 10);
    });
  });

  describe('asum - Sum of Absolute Values', () => {
    it('computes absolute sum', () => {
      const x = [1, -2, 3, -4];

      const result = asum(x);

      expect(result.success).toBe(true);
      expect(result.asum).toBeCloseTo(10, 10);
    });
  });

  describe('iamax - Index of Max Absolute Value', () => {
    it('finds index of max absolute value', () => {
      const x = [1, -5, 3, -2];

      const result = iamax(x);

      expect(result.success).toBe(true);
      expect(result.index).toBe(1); // |-5| is max
    });

    it('returns first index on tie', () => {
      const x = [3, -3, 2];

      const result = iamax(x);

      expect(result.success).toBe(true);
      expect(result.index).toBe(0); // First 3
    });
  });
});
