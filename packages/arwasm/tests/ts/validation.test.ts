/**
 * Validation and Verification Tests
 *
 * Tests for verifyEigs, verifyEign, verifySvds,
 * checkSymmetry, checkPositiveDefinite, checkNormalization.
 */

import { describe, it, expect, beforeAll } from 'vitest';
import {
  loadARPACKModule,
  verifyEigs,
  verifyEign,
  verifySvds,
  checkSymmetry,
  checkPositiveDefinite,
  checkNormalization,
  eigs,
  eign,
  svds,
  diagMatvec,
  denseMatvec,
  denseMatvecT,
} from '../../dist/arwasm.mjs';

// Create 1D Laplacian (symmetric)
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

// Create a non-symmetric matrix operator
function createNonSymmetric(n: number): (x: Float64Array) => Float64Array {
  // Upper triangular matrix
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

describe('Validation and Verification', () => {
  beforeAll(async () => {
    await loadARPACKModule();
  });

  describe('verifyEigs - Symmetric eigenvalue verification', () => {
    it('verifies correct eigenpairs', async () => {
      const n = 10;
      const matvec = createLaplacian1D(n);
      const result = await eigs(matvec, n, 3, {
        which: 'SM',
        return_eigenvectors: true,
      });

      expect(result.success).toBe(true);
      expect(result.eigenvectors).toBeDefined();

      const verification = verifyEigs(
        matvec,
        result.eigenvalues,
        result.eigenvectors!
      );

      expect(verification.isValid).toBe(true);
      expect(verification.maxResidual).toBeLessThan(1e-8);
    });

    it('detects incorrect eigenpairs', async () => {
      const n = 10;
      const matvec = createLaplacian1D(n);
      const result = await eigs(matvec, n, 2, {
        which: 'SM',
        return_eigenvectors: true,
      });

      // Perturb eigenvalues to make them wrong
      const wrongEigenvalues = new Float64Array(result.eigenvalues);
      wrongEigenvalues[0] += 1.0; // Significant error

      const verification = verifyEigs(
        matvec,
        wrongEigenvalues,
        result.eigenvectors!
      );

      expect(verification.isValid).toBe(false);
      expect(verification.maxResidual).toBeGreaterThan(0.1);
    });

    it('checks orthogonality of eigenvectors', async () => {
      const diag = [1, 2, 3, 4, 5];
      const matvec = diagMatvec(new Float64Array(diag));
      const result = await eigs(matvec, 5, 3, {
        which: 'LM',
        return_eigenvectors: true,
      });

      expect(result.success).toBe(true);

      const verification = verifyEigs(
        matvec,
        result.eigenvalues,
        result.eigenvectors!
      );

      // Orthogonality should be close to zero (eigenvectors are orthogonal)
      expect(verification.orthogonality).toBeLessThan(1e-8);
    });

    it('computes relative residuals', async () => {
      const n = 8;
      const matvec = createLaplacian1D(n);
      const result = await eigs(matvec, n, 2, {
        which: 'SM',
        return_eigenvectors: true,
      });

      const verification = verifyEigs(
        matvec,
        result.eigenvalues,
        result.eigenvectors!
      );

      expect(verification.residuals.length).toBe(result.nconv);
      expect(verification.relativeResiduals.length).toBe(result.nconv);

      // All residuals should be small
      for (let i = 0; i < result.nconv; i++) {
        expect(verification.residuals[i]).toBeLessThan(1e-8);
      }
    });
  });

  describe('verifyEign - Non-symmetric eigenvalue verification', () => {
    it('verifies correct eigenpairs', async () => {
      const diag = [1, 2, 3, 4, 5];
      const matvec = diagMatvec(new Float64Array(diag));
      const result = await eign(matvec, 5, 2, {
        which: 'LM',
        return_eigenvectors: true,
      });

      expect(result.success).toBe(true);
      expect(result.eigenvectors).toBeDefined();

      const verification = verifyEign(
        matvec,
        result.eigenvaluesReal,
        result.eigenvaluesImag,
        result.eigenvectors!
      );

      expect(verification.isValid).toBe(true);
      expect(verification.maxResidual).toBeLessThan(1e-8);
    });

    it('detects incorrect eigenpairs', async () => {
      const diag = [1, 2, 3, 4, 5];
      const matvec = diagMatvec(new Float64Array(diag));
      const result = await eign(matvec, 5, 2, {
        which: 'LM',
        return_eigenvectors: true,
      });

      expect(result.eigenvectors).toBeDefined();

      // Perturb eigenvalues
      const wrongReal = new Float64Array(result.eigenvaluesReal);
      wrongReal[0] += 2.0;

      const verification = verifyEign(
        matvec,
        wrongReal,
        result.eigenvaluesImag,
        result.eigenvectors!
      );

      expect(verification.isValid).toBe(false);
    });
  });

  describe('verifySvds - SVD verification', () => {
    it('verifies correct SVD', async () => {
      // 3x3 diagonal matrix
      const matrix = new Float64Array([3, 0, 0, 0, 2, 0, 0, 0, 1]);
      const matvec = denseMatvec(matrix, 3, 3, true);
      const matvecT = denseMatvecT(matrix, 3, 3, true);

      const result = await svds(matvec, matvecT, 3, 3, 2, {
        which: 'LM',
        return_singular_vectors: true,
      });

      expect(result.success).toBe(true);
      expect(result.U).toBeDefined();
      expect(result.Vt).toBeDefined();

      const verification = verifySvds(
        matvec,
        matvecT,
        result.U!,
        result.s,
        result.Vt!
      );

      expect(verification.isValid).toBe(true);
      expect(verification.maxResidual).toBeLessThan(1e-6);
    });

    it('checks U orthonormality', async () => {
      const matrix = new Float64Array([1, 2, 3, 4, 5, 6, 7, 8, 9]);
      const matvec = denseMatvec(matrix, 3, 3, true);
      const matvecT = denseMatvecT(matrix, 3, 3, true);

      const result = await svds(matvec, matvecT, 3, 3, 2, {
        which: 'LM',
        return_singular_vectors: true,
      });

      expect(result.success).toBe(true);

      const verification = verifySvds(
        matvec,
        matvecT,
        result.U!,
        result.s,
        result.Vt!
      );

      expect(verification.orthogonalityU).toBeLessThan(1e-8);
    });

    it('checks Vt orthonormality', async () => {
      const matrix = new Float64Array([1, 2, 3, 4, 5, 6, 7, 8, 9]);
      const matvec = denseMatvec(matrix, 3, 3, true);
      const matvecT = denseMatvecT(matrix, 3, 3, true);

      const result = await svds(matvec, matvecT, 3, 3, 2, {
        which: 'LM',
        return_singular_vectors: true,
      });

      const verification = verifySvds(
        matvec,
        matvecT,
        result.U!,
        result.s,
        result.Vt!
      );

      expect(verification.orthogonalityV).toBeLessThan(1e-8);
    });
  });

  describe('checkSymmetry', () => {
    it('identifies symmetric matrix', () => {
      const matvec = createLaplacian1D(10);
      const result = checkSymmetry(matvec, 10, 5);

      expect(result.isSymmetric).toBe(true);
      expect(result.maxAsymmetry).toBeLessThan(1e-10);
    });

    it('identifies non-symmetric matrix', () => {
      const matvec = createNonSymmetric(10);
      const result = checkSymmetry(matvec, 10, 5);

      expect(result.isSymmetric).toBe(false);
      expect(result.maxAsymmetry).toBeGreaterThan(0.1);
    });

    it('uses specified number of probes', () => {
      const matvec = createLaplacian1D(10);

      // With more probes, should still be symmetric
      const result1 = checkSymmetry(matvec, 10, 3);
      const result2 = checkSymmetry(matvec, 10, 10);

      expect(result1.isSymmetric).toBe(true);
      expect(result2.isSymmetric).toBe(true);
    });

    it('handles diagonal matrix', () => {
      const matvec = diagMatvec(new Float64Array([1, 2, 3, 4, 5]));
      const result = checkSymmetry(matvec, 5, 5);

      expect(result.isSymmetric).toBe(true);
    });
  });

  describe('checkPositiveDefinite', () => {
    it('identifies positive definite matrix', async () => {
      // Laplacian is positive semi-definite, shift to make it positive definite
      const n = 10;
      const matvec = (x: Float64Array): Float64Array => {
        const y = new Float64Array(n);
        for (let i = 0; i < n; i++) {
          y[i] = 3 * x[i]; // Larger diagonal
          if (i > 0) y[i] -= x[i - 1];
          if (i < n - 1) y[i] -= x[i + 1];
        }
        return y;
      };

      const result = await checkPositiveDefinite(matvec, n);

      expect(result.isPositiveDefinite).toBe(true);
      expect(result.smallestEigenvalue).toBeGreaterThan(0);
    });

    it('identifies non-positive-definite matrix', async () => {
      // Matrix with negative eigenvalue
      const diag = [-1, 2, 3, 4, 5];
      const matvec = diagMatvec(new Float64Array(diag));

      const result = await checkPositiveDefinite(matvec, 5);

      expect(result.isPositiveDefinite).toBe(false);
      expect(result.smallestEigenvalue).toBeLessThan(0);
    });

    it('handles positive semi-definite (zero eigenvalue)', async () => {
      // Matrix with zero eigenvalue
      const diag = [0, 1, 2, 3, 4];
      const matvec = diagMatvec(new Float64Array(diag));

      const result = await checkPositiveDefinite(matvec, 5);

      // Zero eigenvalue means not strictly positive definite
      expect(result.isPositiveDefinite).toBe(false);
      expect(Math.abs(result.smallestEigenvalue)).toBeLessThan(1e-5);
    });
  });

  describe('checkNormalization', () => {
    it('identifies normalized vectors', async () => {
      const n = 8;
      const matvec = createLaplacian1D(n);
      const result = await eigs(matvec, n, 3, {
        which: 'SM',
        return_eigenvectors: true,
      });

      expect(result.eigenvectors).toBeDefined();

      const normCheck = checkNormalization(result.eigenvectors!);

      expect(normCheck.isNormalized).toBe(true);
      expect(normCheck.maxDeviation).toBeLessThan(1e-8);
    });

    it('identifies non-normalized vectors', () => {
      // Create unnormalized vectors
      const vectors = [
        new Float64Array([2, 0, 0]), // norm = 2
        new Float64Array([0, 3, 0]), // norm = 3
      ];

      const normCheck = checkNormalization(vectors);

      expect(normCheck.isNormalized).toBe(false);
      expect(normCheck.maxDeviation).toBeGreaterThan(0.5);
    });

    it('reports individual norms', () => {
      const vectors = [
        new Float64Array([1, 0, 0]),
        new Float64Array([0, 2, 0]),
      ];

      const normCheck = checkNormalization(vectors);

      expect(normCheck.norms.length).toBe(2);
      expect(normCheck.norms[0]).toBeCloseTo(1, 10);
      expect(normCheck.norms[1]).toBeCloseTo(2, 10);
    });
  });

  describe('Integration: compute and verify', () => {
    it('full workflow: compute eigenvalues, then verify', async () => {
      const n = 15;
      const matvec = createLaplacian1D(n);

      // Check symmetry first
      const symCheck = checkSymmetry(matvec, n, 5);
      expect(symCheck.isSymmetric).toBe(true);

      // Check positive definiteness
      const pdCheck = await checkPositiveDefinite(matvec, n);
      expect(pdCheck.isPositiveDefinite).toBe(true);

      // Compute eigenvalues
      const result = await eigs(matvec, n, 4, {
        which: 'SM',
        return_eigenvectors: true,
      });
      expect(result.success).toBe(true);

      // Verify normalization
      const normCheck = checkNormalization(result.eigenvectors!);
      expect(normCheck.isNormalized).toBe(true);

      // Verify eigenpairs
      const verification = verifyEigs(
        matvec,
        result.eigenvalues,
        result.eigenvectors!
      );
      expect(verification.isValid).toBe(true);
    });
  });
});
