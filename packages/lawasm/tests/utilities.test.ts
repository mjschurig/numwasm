/**
 * Tests for Matrix Utilities and Property Tests
 */

import { describe, it, expect } from 'vitest';
import {
  transpose,
  triu,
  tril,
  diag,
  trace,
  balance,
  isSymmetric,
  isPositiveDefinite,
  isOrthogonal,
  isSingular,
} from '../src/index.js';
import { toArray2D } from './setup.js';
import './setup.js';

describe('Matrix Utilities', () => {
  describe('transpose', () => {
    it('transposes a square matrix', () => {
      const A = [
        [1, 2],
        [3, 4],
      ];

      const result = transpose(A);

      expect(result.success).toBe(true);
      const At = toArray2D(result.T, 2, 2);
      expect(At[0][0]).toBeCloseTo(1, 10);
      expect(At[0][1]).toBeCloseTo(3, 10);
      expect(At[1][0]).toBeCloseTo(2, 10);
      expect(At[1][1]).toBeCloseTo(4, 10);
    });

    it('transposes a rectangular matrix', () => {
      const A = [
        [1, 2, 3],
        [4, 5, 6],
      ];

      const result = transpose(A);

      expect(result.success).toBe(true);
      expect(result.m).toBe(3);
      expect(result.n).toBe(2);
    });
  });

  describe('triu', () => {
    it('extracts upper triangular part', () => {
      const A = [
        [1, 2, 3],
        [4, 5, 6],
        [7, 8, 9],
      ];

      const result = triu(A);

      expect(result.success).toBe(true);
      const U = toArray2D(result.U, 3, 3);
      expect(U[0][0]).toBeCloseTo(1, 10);
      expect(U[1][0]).toBeCloseTo(0, 10);
      expect(U[2][0]).toBeCloseTo(0, 10);
      expect(U[2][1]).toBeCloseTo(0, 10);
      expect(U[0][2]).toBeCloseTo(3, 10);
    });

    it('handles k parameter', () => {
      const A = [
        [1, 2, 3],
        [4, 5, 6],
        [7, 8, 9],
      ];

      const result = triu(A, { k: 1 });

      expect(result.success).toBe(true);
      const U = toArray2D(result.U, 3, 3);
      // k=1 means start from first superdiagonal
      expect(U[0][0]).toBeCloseTo(0, 10);
      expect(U[1][1]).toBeCloseTo(0, 10);
      expect(U[0][1]).toBeCloseTo(2, 10);
    });
  });

  describe('tril', () => {
    it('extracts lower triangular part', () => {
      const A = [
        [1, 2, 3],
        [4, 5, 6],
        [7, 8, 9],
      ];

      const result = tril(A);

      expect(result.success).toBe(true);
      const L = toArray2D(result.L, 3, 3);
      expect(L[0][0]).toBeCloseTo(1, 10);
      expect(L[0][1]).toBeCloseTo(0, 10);
      expect(L[0][2]).toBeCloseTo(0, 10);
      expect(L[2][0]).toBeCloseTo(7, 10);
    });

    it('handles k parameter', () => {
      const A = [
        [1, 2, 3],
        [4, 5, 6],
        [7, 8, 9],
      ];

      const result = tril(A, { k: -1 });

      expect(result.success).toBe(true);
      const L = toArray2D(result.L, 3, 3);
      // k=-1 means only below main diagonal
      expect(L[0][0]).toBeCloseTo(0, 10);
      expect(L[1][0]).toBeCloseTo(4, 10);
    });
  });

  describe('diag', () => {
    it('extracts diagonal from matrix', () => {
      const A = [
        [1, 2, 3],
        [4, 5, 6],
        [7, 8, 9],
      ];

      const result = diag(A);

      expect(result.success).toBe(true);
      expect(result.diag[0]).toBeCloseTo(1, 10);
      expect(result.diag[1]).toBeCloseTo(5, 10);
      expect(result.diag[2]).toBeCloseTo(9, 10);
    });
  });

  describe('trace', () => {
    it('computes trace of matrix', () => {
      const A = [
        [1, 2, 3],
        [4, 5, 6],
        [7, 8, 9],
      ];

      const result = trace(A);

      expect(result.success).toBe(true);
      expect(result.trace).toBeCloseTo(15, 10); // 1 + 5 + 9
    });

    it('trace of identity is n', () => {
      const I = [
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1],
      ];

      const result = trace(I);

      expect(result.success).toBe(true);
      expect(result.trace).toBeCloseTo(3, 10);
    });
  });

  describe('balance', () => {
    it('balances a matrix for eigenvalue computation', () => {
      const A = [
        [1, 1e10],
        [1e-10, 1],
      ];

      const result = balance(A);

      expect(result.success).toBe(true);
      expect(result.B).toBeDefined();
      expect(result.scale).toBeDefined();
    });
  });
});

describe('Matrix Property Tests', () => {
  describe('isSymmetric', () => {
    it('returns true for symmetric matrix', () => {
      const A = [
        [1, 2, 3],
        [2, 4, 5],
        [3, 5, 6],
      ];

      const result = isSymmetric(A);

      expect(result.success).toBe(true);
      expect(result.isSymmetric).toBe(true);
    });

    it('returns false for non-symmetric matrix', () => {
      const A = [
        [1, 2],
        [3, 4],
      ];

      const result = isSymmetric(A);

      expect(result.success).toBe(true);
      expect(result.isSymmetric).toBe(false);
    });

    it('handles tolerance', () => {
      const A = [
        [1, 2 + 1e-12],
        [2, 4],
      ];

      const result = isSymmetric(A, { tol: 1e-10 });

      expect(result.success).toBe(true);
      expect(result.isSymmetric).toBe(true);
    });
  });

  describe('isPositiveDefinite', () => {
    it('returns true for SPD matrix', () => {
      const A = [
        [4, 2],
        [2, 5],
      ];

      const result = isPositiveDefinite(A);

      expect(result.success).toBe(true);
      expect(result.isPositiveDefinite).toBe(true);
    });

    it('returns false for indefinite matrix', () => {
      const A = [
        [1, 2],
        [2, 1],
      ];

      const result = isPositiveDefinite(A);

      expect(result.success).toBe(true);
      expect(result.isPositiveDefinite).toBe(false);
    });

    it('returns false for negative definite', () => {
      const A = [
        [-2, 0],
        [0, -3],
      ];

      const result = isPositiveDefinite(A);

      expect(result.success).toBe(true);
      expect(result.isPositiveDefinite).toBe(false);
    });
  });

  describe('isOrthogonal', () => {
    it('returns true for rotation matrix', () => {
      const theta = Math.PI / 4;
      const c = Math.cos(theta);
      const s = Math.sin(theta);
      const Q = [
        [c, -s],
        [s, c],
      ];

      const result = isOrthogonal(Q);

      expect(result.success).toBe(true);
      expect(result.isOrthogonal).toBe(true);
    });

    it('returns true for identity', () => {
      const I = [
        [1, 0],
        [0, 1],
      ];

      const result = isOrthogonal(I);

      expect(result.success).toBe(true);
      expect(result.isOrthogonal).toBe(true);
    });

    it('returns true for permutation matrix', () => {
      const P = [
        [0, 1],
        [1, 0],
      ];

      const result = isOrthogonal(P);

      expect(result.success).toBe(true);
      expect(result.isOrthogonal).toBe(true);
    });

    it('returns false for non-orthogonal', () => {
      const A = [
        [1, 1],
        [0, 1],
      ];

      const result = isOrthogonal(A);

      expect(result.success).toBe(true);
      expect(result.isOrthogonal).toBe(false);
    });
  });

  describe('isSingular', () => {
    it('returns true for singular matrix', () => {
      const A = [
        [1, 2],
        [2, 4],
      ];

      const result = isSingular(A);

      expect(result.success).toBe(true);
      expect(result.isSingular).toBe(true);
    });

    it('returns false for invertible matrix', () => {
      const A = [
        [1, 0],
        [0, 1],
      ];

      const result = isSingular(A);

      expect(result.success).toBe(true);
      expect(result.isSingular).toBe(false);
    });

    it('handles near-singular with tolerance', () => {
      const A = [
        [1, 0],
        [0, 1e-15],
      ];

      const resultTight = isSingular(A, { tol: 1e-20 });
      const resultLoose = isSingular(A, { tol: 1e-10 });

      expect(resultTight.isSingular).toBe(false);
      expect(resultLoose.isSingular).toBe(true);
    });
  });
});
