/**
 * Tests for Matrix Functions and Specialized Decompositions
 */

import { describe, it, expect } from 'vitest';
import {
  expm,
  logm,
  sqrtm,
  powm,
  funm,
  polarDecomposition,
  rrqr,
  csd,
} from '../src/index.js';
import { toArray2D, matmulRef } from './setup.js';
import './setup.js';

describe('Matrix Functions', () => {
  describe('expm - Matrix Exponential', () => {
    it('computes exp(0) = I', () => {
      const A = [
        [0, 0],
        [0, 0],
      ];

      const result = expm(A);

      expect(result.success).toBe(true);
      const E = toArray2D(result.E, 2, 2);
      expect(E[0][0]).toBeCloseTo(1, 10);
      expect(E[0][1]).toBeCloseTo(0, 10);
      expect(E[1][0]).toBeCloseTo(0, 10);
      expect(E[1][1]).toBeCloseTo(1, 10);
    });

    it('computes exp of diagonal matrix', () => {
      const A = [
        [1, 0],
        [0, 2],
      ];

      const result = expm(A);

      expect(result.success).toBe(true);
      const E = toArray2D(result.E, 2, 2);
      expect(E[0][0]).toBeCloseTo(Math.E, 8);
      expect(E[1][1]).toBeCloseTo(Math.E * Math.E, 8);
      expect(E[0][1]).toBeCloseTo(0, 10);
    });

    it('computes exp of nilpotent matrix', () => {
      const A = [
        [0, 1],
        [0, 0],
      ];
      // A^2 = 0, so exp(A) = I + A

      const result = expm(A);

      expect(result.success).toBe(true);
      const E = toArray2D(result.E, 2, 2);
      expect(E[0][0]).toBeCloseTo(1, 10);
      expect(E[0][1]).toBeCloseTo(1, 10);
      expect(E[1][0]).toBeCloseTo(0, 10);
      expect(E[1][1]).toBeCloseTo(1, 10);
    });
  });

  describe('logm - Matrix Logarithm', () => {
    it('computes log(I) = 0', () => {
      const A = [
        [1, 0],
        [0, 1],
      ];

      const result = logm(A);

      expect(result.success).toBe(true);
      const L = toArray2D(result.L, 2, 2);
      expect(L[0][0]).toBeCloseTo(0, 8);
      expect(L[1][1]).toBeCloseTo(0, 8);
    });

    it('computes log of positive diagonal matrix', () => {
      const A = [
        [Math.E, 0],
        [0, Math.E * Math.E],
      ];

      const result = logm(A);

      expect(result.success).toBe(true);
      const L = toArray2D(result.L, 2, 2);
      expect(L[0][0]).toBeCloseTo(1, 6);
      expect(L[1][1]).toBeCloseTo(2, 6);
    });
  });

  describe('sqrtm - Matrix Square Root', () => {
    it('computes sqrt(I) = I', () => {
      const A = [
        [1, 0],
        [0, 1],
      ];

      const result = sqrtm(A);

      expect(result.success).toBe(true);
      const S = toArray2D(result.S, 2, 2);
      expect(S[0][0]).toBeCloseTo(1, 8);
      expect(S[1][1]).toBeCloseTo(1, 8);
    });

    it('computes sqrt of diagonal matrix', () => {
      const A = [
        [4, 0],
        [0, 9],
      ];

      const result = sqrtm(A);

      expect(result.success).toBe(true);
      const S = toArray2D(result.S, 2, 2);
      expect(S[0][0]).toBeCloseTo(2, 8);
      expect(S[1][1]).toBeCloseTo(3, 8);
    });

    it('verifies S*S = A', () => {
      const A = [
        [4, 0],
        [0, 9],
      ];

      const result = sqrtm(A);

      expect(result.success).toBe(true);
      const S = toArray2D(result.S, 2, 2);
      const SS = matmulRef(S, S);
      expect(SS[0][0]).toBeCloseTo(4, 6);
      expect(SS[1][1]).toBeCloseTo(9, 6);
    });
  });

  describe('powm - Matrix Power', () => {
    it('computes A^0 = I', () => {
      const A = [
        [2, 1],
        [1, 2],
      ];

      const result = powm(A, 0);

      expect(result.success).toBe(true);
      const P = toArray2D(result.P, 2, 2);
      expect(P[0][0]).toBeCloseTo(1, 10);
      expect(P[0][1]).toBeCloseTo(0, 10);
    });

    it('computes A^1 = A', () => {
      const A = [
        [2, 1],
        [1, 2],
      ];

      const result = powm(A, 1);

      expect(result.success).toBe(true);
      const P = toArray2D(result.P, 2, 2);
      expect(P[0][0]).toBeCloseTo(2, 10);
      expect(P[0][1]).toBeCloseTo(1, 10);
    });

    it('computes A^2', () => {
      const A = [
        [1, 1],
        [0, 1],
      ];

      const result = powm(A, 2);

      expect(result.success).toBe(true);
      const P = toArray2D(result.P, 2, 2);
      // A^2 = [[1,2],[0,1]]
      expect(P[0][0]).toBeCloseTo(1, 10);
      expect(P[0][1]).toBeCloseTo(2, 10);
    });

    it('computes negative integer power', () => {
      const A = [
        [2, 0],
        [0, 4],
      ];

      const result = powm(A, -1);

      expect(result.success).toBe(true);
      const P = toArray2D(result.P, 2, 2);
      expect(P[0][0]).toBeCloseTo(0.5, 10);
      expect(P[1][1]).toBeCloseTo(0.25, 10);
    });
  });

  describe('funm - General Matrix Function', () => {
    it('applies exp function', () => {
      const A = [
        [0, 0],
        [0, 0],
      ];

      const result = funm(A, Math.exp);

      expect(result.success).toBe(true);
      const F = toArray2D(result.F, 2, 2);
      expect(F[0][0]).toBeCloseTo(1, 8);
      expect(F[1][1]).toBeCloseTo(1, 8);
    });

    it('applies sqrt function to diagonal', () => {
      const A = [
        [4, 0],
        [0, 9],
      ];

      const result = funm(A, Math.sqrt);

      expect(result.success).toBe(true);
      const F = toArray2D(result.F, 2, 2);
      expect(F[0][0]).toBeCloseTo(2, 8);
      expect(F[1][1]).toBeCloseTo(3, 8);
    });

    it('applies custom function', () => {
      const A = [
        [1, 0],
        [0, 2],
      ];

      const result = funm(A, (x) => x * x + 1);

      expect(result.success).toBe(true);
      const F = toArray2D(result.F, 2, 2);
      expect(F[0][0]).toBeCloseTo(2, 8); // 1^2 + 1
      expect(F[1][1]).toBeCloseTo(5, 8); // 2^2 + 1
    });
  });
});

describe('Specialized Decompositions', () => {
  describe('polarDecomposition', () => {
    it('decomposes identity matrix', () => {
      const A = [
        [1, 0],
        [0, 1],
      ];

      const result = polarDecomposition(A);

      expect(result.success).toBe(true);
      // U = I, P = I
      const U = toArray2D(result.U, 2, 2);
      const P = toArray2D(result.P, 2, 2);
      expect(U[0][0]).toBeCloseTo(1, 10);
      expect(P[0][0]).toBeCloseTo(1, 10);
    });

    it('decomposes rotation matrix', () => {
      const theta = Math.PI / 4;
      const c = Math.cos(theta);
      const s = Math.sin(theta);
      const Q = [
        [c, -s],
        [s, c],
      ];

      const result = polarDecomposition(Q);

      expect(result.success).toBe(true);
      // For orthogonal Q: U = Q, P = I
      const P = toArray2D(result.P, 2, 2);
      expect(P[0][0]).toBeCloseTo(1, 8);
      expect(P[1][1]).toBeCloseTo(1, 8);
    });

    it('decomposes scaling matrix', () => {
      const A = [
        [3, 0],
        [0, 2],
      ];

      const result = polarDecomposition(A);

      expect(result.success).toBe(true);
      // For positive diagonal: U = I, P = A
      const U = toArray2D(result.U, 2, 2);
      const P = toArray2D(result.P, 2, 2);
      expect(U[0][0]).toBeCloseTo(1, 8);
      expect(P[0][0]).toBeCloseTo(3, 8);
      expect(P[1][1]).toBeCloseTo(2, 8);
    });

    it('verifies U is orthogonal', () => {
      const A = [
        [1, 2],
        [3, 4],
      ];

      const result = polarDecomposition(A);

      expect(result.success).toBe(true);
      const U = toArray2D(result.U, 2, 2);
      // U^T * U should be close to I
      const UtU = matmulRef(
        [
          [U[0][0], U[1][0]],
          [U[0][1], U[1][1]],
        ],
        U
      );
      expect(UtU[0][0]).toBeCloseTo(1, 8);
      expect(UtU[1][1]).toBeCloseTo(1, 8);
    });
  });

  describe('rrqr - Rank-Revealing QR', () => {
    it('computes RRQR of full rank matrix', () => {
      const A = [
        [1, 2],
        [3, 4],
        [5, 6],
      ];

      const result = rrqr(A);

      expect(result.success).toBe(true);
      expect(result.rank).toBe(2);
      expect(result.P).toBeDefined();
    });

    it('detects rank deficiency', () => {
      const A = [
        [1, 2, 3],
        [2, 4, 6], // 2 * row 1
        [3, 6, 9], // 3 * row 1
      ];

      const result = rrqr(A);

      expect(result.success).toBe(true);
      expect(result.rank).toBe(1);
    });

    it('respects tolerance', () => {
      const A = [
        [1, 0, 0],
        [0, 1e-12, 0],
        [0, 0, 1],
      ];

      const resultTight = rrqr(A, { tol: 1e-15 });
      const resultLoose = rrqr(A, { tol: 1e-10 });

      expect(resultTight.rank).toBe(3);
      expect(resultLoose.rank).toBe(2);
    });
  });

  describe('csd - Cosine-Sine Decomposition', () => {
    it('decomposes rotation matrix with p=q=1', () => {
      const theta = Math.PI / 4;
      const c = Math.cos(theta);
      const s = Math.sin(theta);
      const X = [
        [c, -s],
        [s, c],
      ];

      const result = csd(X, 1, 1);

      expect(result.success).toBe(true);
      expect(result.theta).toBeDefined();
      expect(result.theta.length).toBe(1);
      // theta[0] should be close to π/4
      expect(result.theta[0]).toBeCloseTo(Math.PI / 4, 6);
    });

    it('handles identity matrix', () => {
      const X = [
        [1, 0],
        [0, 1],
      ];

      const result = csd(X, 1, 1);

      expect(result.success).toBe(true);
      // For identity, theta = 0 (cos(0) = 1)
      expect(result.theta[0]).toBeCloseTo(0, 8);
    });
  });
});
