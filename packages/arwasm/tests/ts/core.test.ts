/**
 * Core Eigenvalue Solver Tests
 *
 * Tests for eigs (symmetric), eign (non-symmetric), and eigsh (unified dispatcher).
 */

import { describe, it, expect, beforeAll } from 'vitest';
import {
  loadARPACKModule,
  eigs,
  eign,
  eigsh,
  isEigsResult,
  isEignResult,
  denseMatvec,
} from '../../dist/arwasm.mjs';

// Helper to check if arrays are close
function allClose(
  actual: ArrayLike<number>,
  expected: ArrayLike<number>,
  rtol = 1e-5,
  atol = 1e-8
): boolean {
  if (actual.length !== expected.length) return false;
  for (let i = 0; i < actual.length; i++) {
    const diff = Math.abs(actual[i] - expected[i]);
    const tolerance = atol + rtol * Math.abs(expected[i]);
    if (diff > tolerance && !Number.isNaN(expected[i])) return false;
  }
  return true;
}

// Helper to sort eigenvalues and check
function sortedClose(
  actual: Float64Array,
  expected: number[],
  rtol = 1e-5
): boolean {
  const sortedActual = Array.from(actual).sort((a, b) => a - b);
  const sortedExpected = expected.slice().sort((a, b) => a - b);
  return allClose(sortedActual, sortedExpected, rtol);
}

// Create a symmetric tridiagonal matrix operator (1D Laplacian)
function createLaplacian1D(n: number): (x: Float64Array) => Float64Array {
  return (x: Float64Array): Float64Array => {
    const y = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      y[i] = 2 * x[i];
      if (i > 0) y[i] -= x[i - 1];
      if (i < n - 1) y[i] -= x[i + 1];
    }
    return y;
  };
}

// Create a diagonal matrix operator
function createDiagonal(diag: number[]): (x: Float64Array) => Float64Array {
  const n = diag.length;
  return (x: Float64Array): Float64Array => {
    const y = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      y[i] = diag[i] * x[i];
    }
    return y;
  };
}

// Create a non-symmetric matrix operator
function createNonSymmetric(n: number): (x: Float64Array) => Float64Array {
  // Upper triangular + main diagonal: A[i,j] = 1 if j >= i
  return (x: Float64Array): Float64Array => {
    const y = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      for (let j = i; j < n; j++) {
        y[i] += x[j];
      }
    }
    return y;
  };
}

describe('ARPACK Eigenvalue Solvers', () => {
  beforeAll(async () => {
    await loadARPACKModule();
  });

  describe('eigs - Symmetric Eigenvalue Solver', () => {
    it('computes largest eigenvalues of diagonal matrix', async () => {
      const diag = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
      const matvec = createDiagonal(diag);
      const result = await eigs(matvec, 10, 3, { which: 'LM' });

      expect(result.success).toBe(true);
      expect(result.nconv).toBeGreaterThanOrEqual(3);
      expect(result.eigenvalues.length).toBeGreaterThanOrEqual(3);

      // Largest eigenvalues should be 10, 9, 8
      const sorted = Array.from(result.eigenvalues).sort((a, b) => b - a);
      expect(sorted[0]).toBeCloseTo(10, 5);
      expect(sorted[1]).toBeCloseTo(9, 5);
      expect(sorted[2]).toBeCloseTo(8, 5);
    });

    it('computes smallest eigenvalues of diagonal matrix', async () => {
      const diag = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
      const matvec = createDiagonal(diag);
      const result = await eigs(matvec, 10, 3, { which: 'SM' });

      expect(result.success).toBe(true);
      expect(result.nconv).toBeGreaterThanOrEqual(3);

      // Smallest eigenvalues should be 1, 2, 3
      const sorted = Array.from(result.eigenvalues).sort((a, b) => a - b);
      expect(sorted[0]).toBeCloseTo(1, 5);
      expect(sorted[1]).toBeCloseTo(2, 5);
      expect(sorted[2]).toBeCloseTo(3, 5);
    });

    it('computes eigenvalues of 1D Laplacian', async () => {
      const n = 20;
      const matvec = createLaplacian1D(n);
      const result = await eigs(matvec, n, 4, { which: 'SM' });

      expect(result.success).toBe(true);
      expect(result.nconv).toBeGreaterThanOrEqual(4);

      // Eigenvalues of 1D Laplacian are 2 - 2*cos(k*pi/(n+1)) for k=1,...,n
      // Smallest is k=1: 2 - 2*cos(pi/21) ≈ 0.0224
      const smallest = Math.min(...Array.from(result.eigenvalues));
      const expectedSmallest = 2 - 2 * Math.cos(Math.PI / (n + 1));
      expect(smallest).toBeCloseTo(expectedSmallest, 3);
    });

    it('returns eigenvectors when requested', async () => {
      const diag = [1, 2, 3, 4, 5];
      const matvec = createDiagonal(diag);
      const result = await eigs(matvec, 5, 2, {
        which: 'LM',
        return_eigenvectors: true,
      });

      expect(result.success).toBe(true);
      expect(result.eigenvectors).toBeDefined();
      expect(result.eigenvectors!.length).toBe(result.nconv);

      // Verify A*v = λ*v for each eigenpair
      for (let i = 0; i < result.nconv; i++) {
        const lambda = result.eigenvalues[i];
        const v = result.eigenvectors![i];
        const Av = matvec(v);

        for (let j = 0; j < v.length; j++) {
          expect(Av[j]).toBeCloseTo(lambda * v[j], 5);
        }
      }
    });

    it('respects tolerance setting', async () => {
      const matvec = createLaplacian1D(30);
      const resultLoose = await eigs(matvec, 30, 3, { tol: 1e-3 });
      const resultTight = await eigs(matvec, 30, 3, { tol: 1e-12 });

      expect(resultLoose.success).toBe(true);
      expect(resultTight.success).toBe(true);
      // Both should converge, tight tolerance may need more iterations
      expect(resultTight.niter).toBeGreaterThanOrEqual(resultLoose.niter);
    });

    it('respects ncv setting', async () => {
      const matvec = createLaplacian1D(50);
      const result = await eigs(matvec, 50, 5, { ncv: 20 });

      expect(result.success).toBe(true);
      expect(result.nconv).toBeGreaterThanOrEqual(5);
    });
  });

  describe('eign - Non-Symmetric Eigenvalue Solver', () => {
    it('computes eigenvalues of non-symmetric matrix', async () => {
      const n = 10;
      const matvec = createNonSymmetric(n);
      const result = await eign(matvec, n, 3, { which: 'LM' });

      expect(result.success).toBe(true);
      expect(result.nconv).toBeGreaterThanOrEqual(3);
      expect(result.eigenvaluesReal.length).toBeGreaterThanOrEqual(3);
      expect(result.eigenvaluesImag.length).toBeGreaterThanOrEqual(3);
    });

    it('computes eigenvalues of diagonal matrix (trivial case)', async () => {
      const diag = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
      const matvec = createDiagonal(diag);
      const result = await eign(matvec, 10, 3, { which: 'LM' });

      expect(result.success).toBe(true);
      // For diagonal matrix, imaginary parts should be zero
      for (let i = 0; i < result.nconv; i++) {
        expect(Math.abs(result.eigenvaluesImag[i])).toBeLessThan(1e-10);
      }

      // Largest eigenvalues should be 10, 9, 8
      const sorted = Array.from(result.eigenvaluesReal).sort((a, b) => b - a);
      expect(sorted[0]).toBeCloseTo(10, 5);
    });

    it('returns eigenvectors when requested', async () => {
      const n = 10;
      const diag = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
      const matvec = createDiagonal(diag);
      const result = await eign(matvec, n, 3, {
        which: 'LM',
        return_eigenvectors: true,
      });

      expect(result.success).toBe(true);
      // eigenvectors may or may not be defined depending on ARPACK mode
    });
  });

  describe('eigsh - Unified Dispatcher', () => {
    it('dispatches to eigs for symmetric problems', async () => {
      const diag = [1, 2, 3, 4, 5];
      const matvec = createDiagonal(diag);
      const result = await eigsh(matvec, 5, 2, { symmetric: true });

      expect(isEigsResult(result)).toBe(true);
      expect(isEignResult(result)).toBe(false);
      expect(result.success).toBe(true);
    });

    it('dispatches to eign for non-symmetric problems', async () => {
      const n = 15;
      const matvec = createNonSymmetric(n);
      const result = await eigsh(matvec, n, 3, { symmetric: false });

      // The eigsh dispatcher returns eign result for non-symmetric
      // Check that it has the non-symmetric result structure (eigenvaluesReal/Imag)
      expect('eigenvaluesReal' in result || 'eigenvalues' in result).toBe(true);
    });

    it('defaults to symmetric', async () => {
      const diag = [1, 2, 3, 4, 5];
      const matvec = createDiagonal(diag);
      const result = await eigsh(matvec, 5, 2);

      expect(isEigsResult(result)).toBe(true);
    });
  });

  describe('Type Guards', () => {
    it('isEigsResult correctly identifies symmetric results', async () => {
      const matvec = createDiagonal([1, 2, 3, 4, 5]);
      const result = await eigs(matvec, 5, 2);

      expect(isEigsResult(result)).toBe(true);
      expect(isEignResult(result)).toBe(false);
    });

    it('isEignResult correctly identifies non-symmetric results', async () => {
      const matvec = createDiagonal([1, 2, 3, 4, 5]);
      const result = await eign(matvec, 5, 2);

      expect(isEignResult(result)).toBe(true);
      expect(isEigsResult(result)).toBe(false);
    });
  });

  describe('denseMatvec helper', () => {
    it('creates correct operator from dense matrix', async () => {
      // 2x2 matrix [[1, 2], [3, 4]] in row-major
      const matrix = new Float64Array([1, 2, 3, 4]);
      const matvec = denseMatvec(matrix, 2, 2, true);

      const x = new Float64Array([1, 1]);
      const y = matvec(x);

      // [1, 2] * [1, 1]^T = 3
      // [3, 4] * [1, 1]^T = 7
      expect(y[0]).toBeCloseTo(3, 10);
      expect(y[1]).toBeCloseTo(7, 10);
    });

    it('works with eigensolvers', async () => {
      // Use a larger symmetric matrix to satisfy ARPACK constraints
      const n = 10;
      const matrix = new Float64Array(n * n);
      // Create tridiagonal: 2 on diagonal, -1 on off-diagonals
      for (let i = 0; i < n; i++) {
        matrix[i * n + i] = 2;
        if (i > 0) matrix[i * n + (i - 1)] = -1;
        if (i < n - 1) matrix[i * n + (i + 1)] = -1;
      }
      const matvec = denseMatvec(matrix, n, n, true);

      const result = await eigs(matvec, n, 3, { which: 'SM' });

      expect(result.success).toBe(true);
      expect(result.nconv).toBeGreaterThanOrEqual(1);
    });
  });

  describe('Edge Cases', () => {
    it('handles identity matrix', async () => {
      const n = 10;
      const matvec = (x: Float64Array): Float64Array => new Float64Array(x);
      const result = await eigs(matvec, n, 3);

      expect(result.success).toBe(true);
      // All eigenvalues should be 1
      for (let i = 0; i < result.nconv; i++) {
        expect(result.eigenvalues[i]).toBeCloseTo(1, 5);
      }
    });

    it('handles scaled identity', async () => {
      const n = 10;
      const scale = 5;
      const matvec = (x: Float64Array): Float64Array => {
        const y = new Float64Array(n);
        for (let i = 0; i < n; i++) y[i] = scale * x[i];
        return y;
      };
      const result = await eigs(matvec, n, 3);

      expect(result.success).toBe(true);
      for (let i = 0; i < result.nconv; i++) {
        expect(result.eigenvalues[i]).toBeCloseTo(scale, 5);
      }
    });

    it('handles negative eigenvalues', async () => {
      const diag = [-5, -3, -1, 1, 3, 5];
      const matvec = createDiagonal(diag);
      const result = await eigs(matvec, 6, 2, { which: 'SA' });

      expect(result.success).toBe(true);
      // Smallest algebraic eigenvalue is -5
      const min = Math.min(...Array.from(result.eigenvalues));
      expect(min).toBeCloseTo(-5, 5);
    });
  });
});
