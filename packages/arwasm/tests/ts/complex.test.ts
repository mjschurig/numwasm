/**
 * Complex Eigenvalue Solver Tests
 *
 * Tests for zeigs (complex non-symmetric) and zeigsh (complex Hermitian).
 */

import { describe, it, expect, beforeAll } from 'vitest';
import {
  loadARPACKModule,
  zeigs,
  zeigsh,
  type ComplexArray,
} from '../../dist/arwasm.mjs';

// Helper to create a complex diagonal matrix operator
function createComplexDiagonal(
  diagReal: number[],
  diagImag: number[]
): (x: ComplexArray) => ComplexArray {
  const n = diagReal.length;
  return (x: ComplexArray): ComplexArray => {
    const re = new Float64Array(n);
    const im = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      // (a + bi) * (c + di) = (ac - bd) + (ad + bc)i
      const a = diagReal[i];
      const b = diagImag[i];
      const c = x.re[i];
      const d = x.im[i];
      re[i] = a * c - b * d;
      im[i] = a * d + b * c;
    }
    return { re, im };
  };
}

// Helper to create a Hermitian matrix operator (A = A*)
function createHermitian(n: number): (x: ComplexArray) => ComplexArray {
  // Create a simple Hermitian matrix: tridiagonal with real diagonal
  // and complex off-diagonal elements that are conjugates
  return (x: ComplexArray): ComplexArray => {
    const re = new Float64Array(n);
    const im = new Float64Array(n);

    for (let i = 0; i < n; i++) {
      // Diagonal: 2 (real)
      re[i] = 2 * x.re[i];
      im[i] = 2 * x.im[i];

      // Off-diagonal: -1 + 0.5i on upper, -1 - 0.5i on lower (conjugate)
      if (i > 0) {
        // Lower: (-1 - 0.5i) * x[i-1]
        re[i] += -1 * x.re[i - 1] - (-0.5) * x.im[i - 1];
        im[i] += -1 * x.im[i - 1] + (-0.5) * x.re[i - 1];
      }
      if (i < n - 1) {
        // Upper: (-1 + 0.5i) * x[i+1]
        re[i] += -1 * x.re[i + 1] - 0.5 * x.im[i + 1];
        im[i] += -1 * x.im[i + 1] + 0.5 * x.re[i + 1];
      }
    }
    return { re, im };
  };
}

// Helper to create a non-Hermitian complex matrix
function createNonHermitian(n: number): (x: ComplexArray) => ComplexArray {
  // Upper triangular with complex entries
  return (x: ComplexArray): ComplexArray => {
    const re = new Float64Array(n);
    const im = new Float64Array(n);

    for (let i = 0; i < n; i++) {
      for (let j = i; j < n; j++) {
        // A[i,j] = (j - i + 1) + 0.5i
        const aReal = j - i + 1;
        const aImag = 0.5;
        // (a + bi) * (c + di)
        re[i] += aReal * x.re[j] - aImag * x.im[j];
        im[i] += aReal * x.im[j] + aImag * x.re[j];
      }
    }
    return { re, im };
  };
}

describe('Complex Eigenvalue Solvers', () => {
  beforeAll(async () => {
    await loadARPACKModule();
  });

  describe('zeigs - Complex Non-Symmetric Solver', () => {
    it('computes eigenvalues of complex diagonal matrix', async () => {
      // Diagonal matrix with complex eigenvalues
      const diagReal = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
      const diagImag = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0];
      const matvec = createComplexDiagonal(diagReal, diagImag);

      const result = await zeigs(matvec, 10, 3, { which: 'LM' });

      expect(result.success).toBe(true);
      expect(result.nconv).toBeGreaterThanOrEqual(3);
      expect(result.eigenvaluesReal.length).toBeGreaterThanOrEqual(3);
      expect(result.eigenvaluesImag.length).toBeGreaterThanOrEqual(3);

      // Largest magnitude eigenvalues should be 10+i, 9+0.9i, 8+0.8i
      // |10+i| = sqrt(101) ≈ 10.05
      const magnitudes = [];
      for (let i = 0; i < result.nconv; i++) {
        const mag = Math.sqrt(
          result.eigenvaluesReal[i] ** 2 + result.eigenvaluesImag[i] ** 2
        );
        magnitudes.push(mag);
      }
      magnitudes.sort((a, b) => b - a);
      expect(magnitudes[0]).toBeGreaterThan(9);
    });

    it('computes eigenvalues of purely real diagonal', async () => {
      const diagReal = [1, 2, 3, 4, 5];
      const diagImag = [0, 0, 0, 0, 0];
      const matvec = createComplexDiagonal(diagReal, diagImag);

      const result = await zeigs(matvec, 5, 2, { which: 'LM' });

      expect(result.success).toBe(true);
      // Imaginary parts should be essentially zero
      for (let i = 0; i < result.nconv; i++) {
        expect(Math.abs(result.eigenvaluesImag[i])).toBeLessThan(1e-8);
      }
    });

    it('handles non-Hermitian complex matrix', async () => {
      const n = 10;
      const matvec = createNonHermitian(n);

      const result = await zeigs(matvec, n, 3, { which: 'LM' });

      expect(result.success).toBe(true);
      expect(result.nconv).toBeGreaterThanOrEqual(1);
    });

    it('returns eigenvectors when requested', async () => {
      const diagReal = [1, 2, 3, 4, 5];
      const diagImag = [0, 0, 0, 0, 0];
      const matvec = createComplexDiagonal(diagReal, diagImag);

      const result = await zeigs(matvec, 5, 2, {
        which: 'LM',
        return_eigenvectors: true,
      });

      expect(result.success).toBe(true);
      expect(result.eigenvectors).toBeDefined();
      expect(result.eigenvectors!.length).toBe(result.nconv);
    });
  });

  describe('zeigsh - Complex Hermitian Solver', () => {
    it('computes real eigenvalues of Hermitian matrix', async () => {
      const n = 10;
      const matvec = createHermitian(n);

      const result = await zeigsh(matvec, n, 3, { which: 'SM' });

      expect(result.success).toBe(true);
      expect(result.nconv).toBeGreaterThanOrEqual(1);
      // Hermitian matrices have real eigenvalues
      expect(result.eigenvalues).toBeDefined();
    });

    it('computes eigenvalues of real symmetric via complex', async () => {
      // Real symmetric matrix as complex Hermitian
      const n = 8;
      const matvec = (x: ComplexArray): ComplexArray => {
        const re = new Float64Array(n);
        const im = new Float64Array(n);
        // 1D Laplacian (real symmetric)
        for (let i = 0; i < n; i++) {
          re[i] = 2 * x.re[i];
          im[i] = 2 * x.im[i];
          if (i > 0) {
            re[i] -= x.re[i - 1];
            im[i] -= x.im[i - 1];
          }
          if (i < n - 1) {
            re[i] -= x.re[i + 1];
            im[i] -= x.im[i + 1];
          }
        }
        return { re, im };
      };

      const result = await zeigsh(matvec, n, 2, { which: 'SM' });

      expect(result.success).toBe(true);
      // Smallest eigenvalue of 1D Laplacian: 2 - 2*cos(pi/(n+1))
      const expectedSmallest = 2 - 2 * Math.cos(Math.PI / (n + 1));
      const minEig = Math.min(...Array.from(result.eigenvalues));
      expect(minEig).toBeCloseTo(expectedSmallest, 2);
    });

    it('returns eigenvectors when requested', async () => {
      const n = 8;
      const matvec = createHermitian(n);

      const result = await zeigsh(matvec, n, 2, {
        which: 'LM',
        return_eigenvectors: true,
      });

      expect(result.success).toBe(true);
      expect(result.eigenvectors).toBeDefined();
    });
  });

  describe('Which parameter options', () => {
    it('LM finds largest magnitude eigenvalues', async () => {
      const diagReal = [1, 2, 3, 4, 5];
      const diagImag = [0, 0, 0, 0, 0];
      const matvec = createComplexDiagonal(diagReal, diagImag);

      const result = await zeigs(matvec, 5, 2, { which: 'LM' });

      expect(result.success).toBe(true);
      // Should get eigenvalues 5 and 4
      const sorted = Array.from(result.eigenvaluesReal).sort((a, b) => b - a);
      expect(sorted[0]).toBeCloseTo(5, 3);
    });

    it('SM finds smallest magnitude eigenvalues', async () => {
      const diagReal = [1, 2, 3, 4, 5];
      const diagImag = [0, 0, 0, 0, 0];
      const matvec = createComplexDiagonal(diagReal, diagImag);

      const result = await zeigs(matvec, 5, 2, { which: 'SM' });

      expect(result.success).toBe(true);
      // Should get eigenvalues 1 and 2
      const sorted = Array.from(result.eigenvaluesReal).sort((a, b) => a - b);
      expect(sorted[0]).toBeCloseTo(1, 3);
    });
  });
});
