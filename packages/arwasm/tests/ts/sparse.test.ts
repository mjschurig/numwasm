/**
 * Sparse Matrix Helper Tests
 *
 * Tests for csrMatvec, cscMatvec, cooMatvec, denseMatvec,
 * diagMatvec, tridiagMatvec, bandedMatvec.
 */

import { describe, it, expect, beforeAll } from 'vitest';
import {
  loadARPACKModule,
  csrMatvec,
  csrMatvecT,
  cscMatvec,
  cscMatvecT,
  cooMatvec,
  cooMatvecT,
  denseMatvec,
  denseMatvecT,
  diagMatvec,
  diagMatvecInv,
  diagMatvecSqrt,
  diagMatvecInvSqrt,
  tridiagMatvec,
  symTridiagMatvec,
  toeplitzTridiagMatvec,
  bandedMatvec,
  symBandedMatvec,
  eigs,
} from '../../dist/arwasm.mjs';

describe('Sparse Matrix Helpers', () => {
  beforeAll(async () => {
    await loadARPACKModule();
  });

  describe('csrMatvec - CSR format', () => {
    it('computes y = A*x for CSR matrix', () => {
      // Matrix: [[1, 2, 0], [0, 3, 4], [5, 0, 6]]
      // CSR format:
      // data = [1, 2, 3, 4, 5, 6]
      // indices = [0, 1, 1, 2, 0, 2]
      // indptr = [0, 2, 4, 6]
      const indptr = new Int32Array([0, 2, 4, 6]);
      const indices = new Int32Array([0, 1, 1, 2, 0, 2]);
      const data = new Float64Array([1, 2, 3, 4, 5, 6]);

      const matvec = csrMatvec(indptr, indices, data, [3, 3]);
      const x = new Float64Array([1, 2, 3]);
      const y = matvec(x);

      // [1, 2, 0] * [1, 2, 3]^T = 1 + 4 = 5
      // [0, 3, 4] * [1, 2, 3]^T = 6 + 12 = 18
      // [5, 0, 6] * [1, 2, 3]^T = 5 + 18 = 23
      expect(y[0]).toBeCloseTo(5, 10);
      expect(y[1]).toBeCloseTo(18, 10);
      expect(y[2]).toBeCloseTo(23, 10);
    });

    it('computes y = A^T*x for CSR matrix', () => {
      const indptr = new Int32Array([0, 2, 4, 6]);
      const indices = new Int32Array([0, 1, 1, 2, 0, 2]);
      const data = new Float64Array([1, 2, 3, 4, 5, 6]);

      const matvecT = csrMatvecT(indptr, indices, data, [3, 3]);
      const x = new Float64Array([1, 2, 3]);
      const y = matvecT(x);

      // A^T = [[1, 0, 5], [2, 3, 0], [0, 4, 6]]
      // [1, 0, 5] * [1, 2, 3]^T = 1 + 15 = 16
      // [2, 3, 0] * [1, 2, 3]^T = 2 + 6 = 8
      // [0, 4, 6] * [1, 2, 3]^T = 8 + 18 = 26
      expect(y[0]).toBeCloseTo(16, 10);
      expect(y[1]).toBeCloseTo(8, 10);
      expect(y[2]).toBeCloseTo(26, 10);
    });

    it('handles identity matrix', () => {
      // Identity in CSR
      const n = 5;
      const indptr = new Int32Array([0, 1, 2, 3, 4, 5]);
      const indices = new Int32Array([0, 1, 2, 3, 4]);
      const data = new Float64Array([1, 1, 1, 1, 1]);

      const matvec = csrMatvec(indptr, indices, data, [n, n]);
      const x = new Float64Array([1, 2, 3, 4, 5]);
      const y = matvec(x);

      for (let i = 0; i < n; i++) {
        expect(y[i]).toBeCloseTo(x[i], 10);
      }
    });
  });

  describe('cscMatvec - CSC format', () => {
    it('computes y = A*x for CSC matrix', () => {
      // Matrix: [[1, 2, 0], [0, 3, 4], [5, 0, 6]]
      // CSC format (column-major):
      // indptr = [0, 2, 4, 6] (columns)
      // indices = [0, 2, 0, 1, 1, 2] (row indices)
      // data = [1, 5, 2, 3, 4, 6]
      const indptr = new Int32Array([0, 2, 4, 6]);
      const indices = new Int32Array([0, 2, 0, 1, 1, 2]);
      const data = new Float64Array([1, 5, 2, 3, 4, 6]);

      const matvec = cscMatvec(indptr, indices, data, [3, 3]);
      const x = new Float64Array([1, 2, 3]);
      const y = matvec(x);

      // Same result as CSR
      expect(y[0]).toBeCloseTo(5, 10);
      expect(y[1]).toBeCloseTo(18, 10);
      expect(y[2]).toBeCloseTo(23, 10);
    });

    it('computes y = A^T*x for CSC matrix', () => {
      const indptr = new Int32Array([0, 2, 4, 6]);
      const indices = new Int32Array([0, 2, 0, 1, 1, 2]);
      const data = new Float64Array([1, 5, 2, 3, 4, 6]);

      const matvecT = cscMatvecT(indptr, indices, data, [3, 3]);
      const x = new Float64Array([1, 2, 3]);
      const y = matvecT(x);

      expect(y[0]).toBeCloseTo(16, 10);
      expect(y[1]).toBeCloseTo(8, 10);
      expect(y[2]).toBeCloseTo(26, 10);
    });
  });

  describe('cooMatvec - COO format', () => {
    it('computes y = A*x for COO matrix', () => {
      // Matrix: [[1, 2, 0], [0, 3, 4], [5, 0, 6]]
      const rows = new Int32Array([0, 0, 1, 1, 2, 2]);
      const cols = new Int32Array([0, 1, 1, 2, 0, 2]);
      const values = new Float64Array([1, 2, 3, 4, 5, 6]);

      const matvec = cooMatvec(rows, cols, values, [3, 3]);
      const x = new Float64Array([1, 2, 3]);
      const y = matvec(x);

      expect(y[0]).toBeCloseTo(5, 10);
      expect(y[1]).toBeCloseTo(18, 10);
      expect(y[2]).toBeCloseTo(23, 10);
    });

    it('computes y = A^T*x for COO matrix', () => {
      const rows = new Int32Array([0, 0, 1, 1, 2, 2]);
      const cols = new Int32Array([0, 1, 1, 2, 0, 2]);
      const values = new Float64Array([1, 2, 3, 4, 5, 6]);

      const matvecT = cooMatvecT(rows, cols, values, [3, 3]);
      const x = new Float64Array([1, 2, 3]);
      const y = matvecT(x);

      expect(y[0]).toBeCloseTo(16, 10);
      expect(y[1]).toBeCloseTo(8, 10);
      expect(y[2]).toBeCloseTo(26, 10);
    });

    it('handles duplicate entries (sums them)', () => {
      // Two entries at (0,0): 1 and 2 should sum to 3
      const rows = new Int32Array([0, 0, 1]);
      const cols = new Int32Array([0, 0, 1]);
      const values = new Float64Array([1, 2, 4]);

      const matvec = cooMatvec(rows, cols, values, [2, 2]);
      const x = new Float64Array([1, 1]);
      const y = matvec(x);

      // [[3, 0], [0, 4]] * [1, 1]^T = [3, 4]
      expect(y[0]).toBeCloseTo(3, 10);
      expect(y[1]).toBeCloseTo(4, 10);
    });
  });

  describe('denseMatvec - Dense format', () => {
    it('computes y = A*x for row-major dense matrix', () => {
      // Matrix: [[1, 2], [3, 4]]
      const matrix = new Float64Array([1, 2, 3, 4]);
      const matvec = denseMatvec(matrix, 2, 2, true);
      const x = new Float64Array([1, 2]);
      const y = matvec(x);

      // [1, 2] * [1, 2]^T = 5
      // [3, 4] * [1, 2]^T = 11
      expect(y[0]).toBeCloseTo(5, 10);
      expect(y[1]).toBeCloseTo(11, 10);
    });

    it('computes y = A*x for column-major dense matrix', () => {
      // Matrix: [[1, 2], [3, 4]] in column-major = [1, 3, 2, 4]
      const matrix = new Float64Array([1, 3, 2, 4]);
      const matvec = denseMatvec(matrix, 2, 2, false);
      const x = new Float64Array([1, 2]);
      const y = matvec(x);

      expect(y[0]).toBeCloseTo(5, 10);
      expect(y[1]).toBeCloseTo(11, 10);
    });

    it('computes y = A^T*x for dense matrix', () => {
      const matrix = new Float64Array([1, 2, 3, 4]);
      const matvecT = denseMatvecT(matrix, 2, 2, true);
      const x = new Float64Array([1, 2]);
      const y = matvecT(x);

      // A^T = [[1, 3], [2, 4]]
      // [1, 3] * [1, 2]^T = 7
      // [2, 4] * [1, 2]^T = 10
      expect(y[0]).toBeCloseTo(7, 10);
      expect(y[1]).toBeCloseTo(10, 10);
    });

    it('handles rectangular matrices', () => {
      // 2x3 matrix
      const matrix = new Float64Array([1, 2, 3, 4, 5, 6]);
      const matvec = denseMatvec(matrix, 2, 3, true);
      const x = new Float64Array([1, 2, 3]);
      const y = matvec(x);

      // [1, 2, 3] * [1, 2, 3]^T = 14
      // [4, 5, 6] * [1, 2, 3]^T = 32
      expect(y[0]).toBeCloseTo(14, 10);
      expect(y[1]).toBeCloseTo(32, 10);
    });
  });

  describe('diagMatvec - Diagonal format', () => {
    it('computes y = D*x for diagonal matrix', () => {
      const diag = new Float64Array([1, 2, 3, 4]);
      const matvec = diagMatvec(diag);
      const x = new Float64Array([2, 3, 4, 5]);
      const y = matvec(x);

      expect(y[0]).toBeCloseTo(2, 10);
      expect(y[1]).toBeCloseTo(6, 10);
      expect(y[2]).toBeCloseTo(12, 10);
      expect(y[3]).toBeCloseTo(20, 10);
    });

    it('computes y = D^{-1}*x for diagonal matrix', () => {
      const diag = new Float64Array([1, 2, 4, 8]);
      const matvec = diagMatvecInv(diag);
      const x = new Float64Array([1, 2, 4, 8]);
      const y = matvec(x);

      expect(y[0]).toBeCloseTo(1, 10);
      expect(y[1]).toBeCloseTo(1, 10);
      expect(y[2]).toBeCloseTo(1, 10);
      expect(y[3]).toBeCloseTo(1, 10);
    });

    it('computes y = D^{1/2}*x for diagonal matrix', () => {
      const diag = new Float64Array([1, 4, 9, 16]);
      const matvec = diagMatvecSqrt(diag);
      const x = new Float64Array([1, 1, 1, 1]);
      const y = matvec(x);

      expect(y[0]).toBeCloseTo(1, 10);
      expect(y[1]).toBeCloseTo(2, 10);
      expect(y[2]).toBeCloseTo(3, 10);
      expect(y[3]).toBeCloseTo(4, 10);
    });

    it('computes y = D^{-1/2}*x for diagonal matrix', () => {
      const diag = new Float64Array([1, 4, 9, 16]);
      const matvec = diagMatvecInvSqrt(diag);
      const x = new Float64Array([1, 2, 3, 4]);
      const y = matvec(x);

      expect(y[0]).toBeCloseTo(1, 10);
      expect(y[1]).toBeCloseTo(1, 10);
      expect(y[2]).toBeCloseTo(1, 10);
      expect(y[3]).toBeCloseTo(1, 10);
    });
  });

  describe('tridiagMatvec - Tridiagonal format', () => {
    it('computes y = T*x for tridiagonal matrix', () => {
      // Matrix: [[2, -1, 0], [-1, 2, -1], [0, -1, 2]]
      const lower = new Float64Array([-1, -1]);
      const diag = new Float64Array([2, 2, 2]);
      const upper = new Float64Array([-1, -1]);

      const matvec = tridiagMatvec(lower, diag, upper);
      const x = new Float64Array([1, 2, 3]);
      const y = matvec(x);

      // Row 0: 2*1 + (-1)*2 = 0
      // Row 1: (-1)*1 + 2*2 + (-1)*3 = 0
      // Row 2: (-1)*2 + 2*3 = 4
      expect(y[0]).toBeCloseTo(0, 10);
      expect(y[1]).toBeCloseTo(0, 10);
      expect(y[2]).toBeCloseTo(4, 10);
    });

    it('symTridiagMatvec uses symmetric structure', () => {
      const diag = new Float64Array([2, 2, 2]);
      const offDiag = new Float64Array([-1, -1]);

      const matvec = symTridiagMatvec(diag, offDiag);
      const x = new Float64Array([1, 2, 3]);
      const y = matvec(x);

      expect(y[0]).toBeCloseTo(0, 10);
      expect(y[1]).toBeCloseTo(0, 10);
      expect(y[2]).toBeCloseTo(4, 10);
    });

    it('toeplitzTridiagMatvec uses constant diagonals', () => {
      // Constant tridiagonal: lower=-1, diag=2, upper=-1, n=3
      const matvec = toeplitzTridiagMatvec(-1, 2, -1, 3);
      const x = new Float64Array([1, 2, 3]);
      const y = matvec(x);

      // Row 0: 2*1 + (-1)*2 = 0
      // Row 1: (-1)*1 + 2*2 + (-1)*3 = 0
      // Row 2: (-1)*2 + 2*3 = 4
      expect(y[0]).toBeCloseTo(0, 10);
      expect(y[1]).toBeCloseTo(0, 10);
      expect(y[2]).toBeCloseTo(4, 10);
    });
  });

  describe('bandedMatvec - General banded format', () => {
    it('computes y = B*x for banded matrix', () => {
      // Pentadiagonal: offsets [-2, -1, 0, 1, 2]
      // 5x5 matrix with bands
      const n = 5;
      const diag = new Float64Array([4, 4, 4, 4, 4]);
      const sub1 = new Float64Array([-1, -1, -1, -1]); // offset -1
      const sup1 = new Float64Array([-1, -1, -1, -1]); // offset +1
      const sub2 = new Float64Array([0.5, 0.5, 0.5]); // offset -2
      const sup2 = new Float64Array([0.5, 0.5, 0.5]); // offset +2

      const bands = [sub2, sub1, diag, sup1, sup2];
      const offsets = [-2, -1, 0, 1, 2];

      const matvec = bandedMatvec(bands, offsets, n);
      const x = new Float64Array([1, 1, 1, 1, 1]);
      const y = matvec(x);

      // Each interior row: 0.5 - 1 + 4 - 1 + 0.5 = 3
      // But boundary rows have fewer terms
      expect(y[2]).toBeCloseTo(3, 10); // Middle row
    });

    it('symBandedMatvec uses symmetric structure', () => {
      const n = 4;
      const diag = new Float64Array([2, 2, 2, 2]);
      const band1 = new Float64Array([-1, -1, -1]); // offset 1 (and -1), length n-1

      // symBandedMatvec takes (bands, n) where bands[0]=diag, bands[1]=off1, etc.
      const matvec = symBandedMatvec([diag, band1], n);
      const x = new Float64Array([1, 2, 3, 4]);
      const y = matvec(x);

      // This is a symmetric tridiagonal
      // Row 0: 2*1 - 1*2 = 0
      // Row 1: -1*1 + 2*2 - 1*3 = 0
      // Row 2: -1*2 + 2*3 - 1*4 = 0
      // Row 3: -1*3 + 2*4 = 5
      expect(y[0]).toBeCloseTo(0, 10);
      expect(y[1]).toBeCloseTo(0, 10);
      expect(y[2]).toBeCloseTo(0, 10);
      expect(y[3]).toBeCloseTo(5, 10);
    });
  });

  describe('Integration with eigensolvers', () => {
    it('csrMatvec works with eigs', async () => {
      // Symmetric tridiagonal in CSR
      const n = 10;
      const nnz = 3 * n - 2;
      const indptr = new Int32Array(n + 1);
      const indices = new Int32Array(nnz);
      const data = new Float64Array(nnz);

      let idx = 0;
      for (let i = 0; i < n; i++) {
        indptr[i] = idx;
        if (i > 0) {
          indices[idx] = i - 1;
          data[idx] = -1;
          idx++;
        }
        indices[idx] = i;
        data[idx] = 2;
        idx++;
        if (i < n - 1) {
          indices[idx] = i + 1;
          data[idx] = -1;
          idx++;
        }
      }
      indptr[n] = idx;

      const matvec = csrMatvec(indptr, indices, data, [n, n]);
      const result = await eigs(matvec, n, 3, { which: 'SM' });

      expect(result.success).toBe(true);
      expect(result.nconv).toBeGreaterThanOrEqual(3);

      // Smallest eigenvalue of 1D Laplacian
      const expectedSmallest = 2 - 2 * Math.cos(Math.PI / (n + 1));
      const minEig = Math.min(...Array.from(result.eigenvalues));
      expect(minEig).toBeCloseTo(expectedSmallest, 3);
    });

    it('diagMatvec works with eigs', async () => {
      const diag = new Float64Array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10]);
      const matvec = diagMatvec(diag);

      const result = await eigs(matvec, 10, 3, { which: 'LM' });

      expect(result.success).toBe(true);
      const sorted = Array.from(result.eigenvalues).sort((a, b) => b - a);
      expect(sorted[0]).toBeCloseTo(10, 5);
    });

    it('symTridiagMatvec works with eigs', async () => {
      const n = 15;
      const diag = new Float64Array(n).fill(2);
      const offDiag = new Float64Array(n - 1).fill(-1);

      const matvec = symTridiagMatvec(diag, offDiag);
      const result = await eigs(matvec, n, 3, { which: 'SM' });

      expect(result.success).toBe(true);
      const expectedSmallest = 2 - 2 * Math.cos(Math.PI / (n + 1));
      const minEig = Math.min(...Array.from(result.eigenvalues));
      expect(minEig).toBeCloseTo(expectedSmallest, 3);
    });
  });
});
