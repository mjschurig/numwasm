/**
 * Operator Combinator Tests
 *
 * Tests for addMatvec, mulMatvec, shiftMatvec, scaleMatvec,
 * transposeMatvec, symmetrizeMatvec, etc.
 */

import { describe, it, expect, beforeAll } from 'vitest';
import {
  loadARPACKModule,
  addMatvec,
  mulMatvec,
  shiftMatvec,
  scaleMatvec,
  transposeMatvec,
  symmetrizeMatvec,
  identityMatvec,
  negateMatvec,
  powerMatvec,
  blockDiagMatvec,
  diagMatvec,
  denseMatvec,
  denseMatvecT,
  eigs,
} from '../../dist/arwasm.mjs';

// Helper to check array equality
function arrayClose(
  a: Float64Array,
  b: number[],
  tol = 1e-10
): boolean {
  if (a.length !== b.length) return false;
  for (let i = 0; i < a.length; i++) {
    if (Math.abs(a[i] - b[i]) > tol) return false;
  }
  return true;
}

describe('Operator Combinators', () => {
  beforeAll(async () => {
    await loadARPACKModule();
  });

  describe('addMatvec - Linear combination', () => {
    it('computes y = (A + B) * x with default weights', () => {
      const A = diagMatvec(new Float64Array([1, 2, 3]));
      const B = diagMatvec(new Float64Array([4, 5, 6]));

      const combined = addMatvec(A, B);
      const x = new Float64Array([1, 1, 1]);
      const y = combined(x);

      // (1+4, 2+5, 3+6) * (1, 1, 1) = (5, 7, 9)
      expect(arrayClose(y, [5, 7, 9])).toBe(true);
    });

    it('computes y = (αA + βB) * x with custom weights', () => {
      const A = diagMatvec(new Float64Array([1, 2, 3]));
      const B = diagMatvec(new Float64Array([4, 5, 6]));

      const combined = addMatvec(A, B, 2, 3);
      const x = new Float64Array([1, 1, 1]);
      const y = combined(x);

      // (2*1 + 3*4, 2*2 + 3*5, 2*3 + 3*6) * (1, 1, 1) = (14, 19, 24)
      expect(arrayClose(y, [14, 19, 24])).toBe(true);
    });

    it('computes A - B with α=1, β=-1', () => {
      const A = diagMatvec(new Float64Array([5, 6, 7]));
      const B = diagMatvec(new Float64Array([1, 2, 3]));

      const combined = addMatvec(A, B, 1, -1);
      const x = new Float64Array([1, 1, 1]);
      const y = combined(x);

      expect(arrayClose(y, [4, 4, 4])).toBe(true);
    });
  });

  describe('mulMatvec - Operator composition', () => {
    it('computes y = A * B * x', () => {
      const A = diagMatvec(new Float64Array([2, 2, 2]));
      const B = diagMatvec(new Float64Array([1, 2, 3]));

      const combined = mulMatvec(A, B);
      const x = new Float64Array([1, 1, 1]);
      const y = combined(x);

      // B*x = (1, 2, 3), A*(B*x) = (2, 4, 6)
      expect(arrayClose(y, [2, 4, 6])).toBe(true);
    });

    it('order matters for non-commuting operators', () => {
      // A = [[1, 2], [0, 1]], B = [[1, 0], [1, 1]]
      const A = denseMatvec(new Float64Array([1, 2, 0, 1]), 2, 2, true);
      const B = denseMatvec(new Float64Array([1, 0, 1, 1]), 2, 2, true);

      const AB = mulMatvec(A, B);
      const BA = mulMatvec(B, A);
      const x = new Float64Array([1, 1]);

      const yAB = AB(x);
      const yBA = BA(x);

      // They should be different for non-commuting matrices
      // A*B = [[3, 2], [1, 1]], B*A = [[1, 2], [1, 3]]
      // AB*[1,1] = [5, 2], BA*[1,1] = [3, 4]
      expect(arrayClose(yAB, [5, 2])).toBe(true);
      expect(arrayClose(yBA, [3, 4])).toBe(true);
    });
  });

  describe('shiftMatvec - (A - σI)', () => {
    it('computes y = (A - σI) * x', () => {
      const A = diagMatvec(new Float64Array([5, 6, 7]));
      const sigma = 3;

      const shifted = shiftMatvec(A, sigma);
      const x = new Float64Array([1, 1, 1]);
      const y = shifted(x);

      // (5-3, 6-3, 7-3) * (1, 1, 1) = (2, 3, 4)
      expect(arrayClose(y, [2, 3, 4])).toBe(true);
    });

    it('shifts eigenvalues by sigma', async () => {
      const diag = [1, 2, 3, 4, 5];
      const A = diagMatvec(new Float64Array(diag));
      const sigma = 2;

      const shifted = shiftMatvec(A, sigma);
      const result = await eigs(shifted, 5, 3, { which: 'LM' });

      expect(result.success).toBe(true);
      // Original eigenvalues: 1, 2, 3, 4, 5
      // Shifted: -1, 0, 1, 2, 3
      // Largest magnitude: -1, 3, 2
      const sorted = Array.from(result.eigenvalues).sort((a, b) => b - a);
      expect(sorted[0]).toBeCloseTo(3, 3);
    });
  });

  describe('scaleMatvec - αA', () => {
    it('computes y = α * A * x', () => {
      const A = diagMatvec(new Float64Array([1, 2, 3]));
      const alpha = 5;

      const scaled = scaleMatvec(A, alpha);
      const x = new Float64Array([1, 1, 1]);
      const y = scaled(x);

      expect(arrayClose(y, [5, 10, 15])).toBe(true);
    });

    it('scales eigenvalues by alpha', async () => {
      const diag = [1, 2, 3, 4, 5];
      const A = diagMatvec(new Float64Array(diag));
      const alpha = 2;

      const scaled = scaleMatvec(A, alpha);
      const result = await eigs(scaled, 5, 2, { which: 'LM' });

      expect(result.success).toBe(true);
      const sorted = Array.from(result.eigenvalues).sort((a, b) => b - a);
      expect(sorted[0]).toBeCloseTo(10, 3); // 2 * 5 = 10
    });

    it('handles negative scale', () => {
      const A = diagMatvec(new Float64Array([1, 2, 3]));
      const scaled = scaleMatvec(A, -1);
      const x = new Float64Array([1, 1, 1]);
      const y = scaled(x);

      expect(arrayClose(y, [-1, -2, -3])).toBe(true);
    });
  });

  describe('transposeMatvec - A^T', () => {
    it('computes y = A^T * x', () => {
      // A = [[1, 2], [3, 4]]
      const matrix = new Float64Array([1, 2, 3, 4]);
      const A = denseMatvec(matrix, 2, 2, true);
      const AT = denseMatvecT(matrix, 2, 2, true);

      const transposed = transposeMatvec(A, 2, 2);
      const x = new Float64Array([1, 1]);

      const y1 = transposed(x);
      const y2 = AT(x);

      // Both should give same result
      expect(arrayClose(y1, Array.from(y2))).toBe(true);
    });

    it('handles rectangular matrices', () => {
      // A = [[1, 2, 3], [4, 5, 6]] (2x3)
      // A^T is 3x2
      const matrix = new Float64Array([1, 2, 3, 4, 5, 6]);
      const A = denseMatvec(matrix, 2, 3, true);

      const transposed = transposeMatvec(A, 2, 3);
      const x = new Float64Array([1, 1]); // 2-vector for A^T

      const y = transposed(x);
      // A^T = [[1, 4], [2, 5], [3, 6]]
      // A^T * [1, 1] = [5, 7, 9]
      expect(y.length).toBe(3);
      expect(arrayClose(y, [5, 7, 9])).toBe(true);
    });
  });

  describe('symmetrizeMatvec - (A + A^T) / 2', () => {
    it('symmetrizes a matrix', () => {
      // A = [[1, 2], [3, 4]]
      // A^T = [[1, 3], [2, 4]]
      // (A + A^T)/2 = [[1, 2.5], [2.5, 4]]
      const matrix = new Float64Array([1, 2, 3, 4]);
      const A = denseMatvec(matrix, 2, 2, true);
      const AT = denseMatvecT(matrix, 2, 2, true);

      const sym = symmetrizeMatvec(A, AT);
      const x = new Float64Array([1, 0]);
      const y = sym(x);

      // [[1, 2.5], [2.5, 4]] * [1, 0] = [1, 2.5]
      expect(y[0]).toBeCloseTo(1, 10);
      expect(y[1]).toBeCloseTo(2.5, 10);
    });

    it('already symmetric matrix stays same', () => {
      // Symmetric: [[2, 1], [1, 2]]
      const matrix = new Float64Array([2, 1, 1, 2]);
      const A = denseMatvec(matrix, 2, 2, true);
      const AT = denseMatvecT(matrix, 2, 2, true);

      const sym = symmetrizeMatvec(A, AT);
      const x = new Float64Array([1, 1]);

      const yA = A(x);
      const ySym = sym(x);

      expect(arrayClose(ySym, Array.from(yA))).toBe(true);
    });
  });

  describe('identityMatvec - I', () => {
    it('returns input unchanged', () => {
      const I = identityMatvec(5);
      const x = new Float64Array([1, 2, 3, 4, 5]);
      const y = I(x);

      expect(arrayClose(y, [1, 2, 3, 4, 5])).toBe(true);
    });

    it('has eigenvalue 1', async () => {
      const I = identityMatvec(10);
      const result = await eigs(I, 10, 3);

      expect(result.success).toBe(true);
      for (let i = 0; i < result.nconv; i++) {
        expect(result.eigenvalues[i]).toBeCloseTo(1, 5);
      }
    });
  });

  describe('negateMatvec - -A', () => {
    it('negates the operator', () => {
      const A = diagMatvec(new Float64Array([1, 2, 3]));
      const negA = negateMatvec(A);
      const x = new Float64Array([1, 1, 1]);
      const y = negA(x);

      expect(arrayClose(y, [-1, -2, -3])).toBe(true);
    });

    it('equivalent to scale(-1)', () => {
      const A = diagMatvec(new Float64Array([1, 2, 3]));
      const negA = negateMatvec(A);
      const scaledA = scaleMatvec(A, -1);
      const x = new Float64Array([1, 2, 3]);

      const y1 = negA(x);
      const y2 = scaledA(x);

      expect(arrayClose(y1, Array.from(y2))).toBe(true);
    });
  });

  describe('powerMatvec - A^k', () => {
    it('computes A^2', () => {
      const A = diagMatvec(new Float64Array([2, 3, 4]));
      const A2 = powerMatvec(A, 2);
      const x = new Float64Array([1, 1, 1]);
      const y = A2(x);

      // A^2 = diag(4, 9, 16)
      expect(arrayClose(y, [4, 9, 16])).toBe(true);
    });

    it('computes A^3', () => {
      const A = diagMatvec(new Float64Array([2, 3]));
      const A3 = powerMatvec(A, 3);
      const x = new Float64Array([1, 1]);
      const y = A3(x);

      // A^3 = diag(8, 27)
      expect(arrayClose(y, [8, 27])).toBe(true);
    });

    it('A^1 equals A', () => {
      const A = diagMatvec(new Float64Array([1, 2, 3]));
      const A1 = powerMatvec(A, 1);
      const x = new Float64Array([1, 1, 1]);

      const yA = A(x);
      const yA1 = A1(x);

      expect(arrayClose(yA1, Array.from(yA))).toBe(true);
    });

    it('A^0 is identity', () => {
      const A = diagMatvec(new Float64Array([5, 10, 15]));
      const A0 = powerMatvec(A, 0);
      const x = new Float64Array([1, 2, 3]);
      const y = A0(x);

      expect(arrayClose(y, [1, 2, 3])).toBe(true);
    });
  });

  describe('blockDiagMatvec - block diagonal', () => {
    it('combines two diagonal operators', () => {
      const A = diagMatvec(new Float64Array([1, 2]));
      const B = diagMatvec(new Float64Array([3, 4, 5]));

      const block = blockDiagMatvec([A, B], [2, 3]);
      const x = new Float64Array([1, 1, 1, 1, 1]);
      const y = block(x);

      // Block diagonal: [1, 2, 3, 4, 5] applied element-wise
      expect(arrayClose(y, [1, 2, 3, 4, 5])).toBe(true);
    });

    it('preserves block structure', () => {
      // [[A, 0], [0, B]] where A and B are 2x2
      const A = denseMatvec(new Float64Array([1, 2, 3, 4]), 2, 2, true);
      const B = denseMatvec(new Float64Array([5, 6, 7, 8]), 2, 2, true);

      const block = blockDiagMatvec([A, B], [2, 2]);
      const x = new Float64Array([1, 0, 0, 1]);
      const y = block(x);

      // A*[1,0] = [1, 3], B*[0,1] = [6, 8]
      expect(arrayClose(y, [1, 3, 6, 8])).toBe(true);
    });

    it('single block equals original', () => {
      const A = diagMatvec(new Float64Array([1, 2, 3]));
      const block = blockDiagMatvec([A], [3]);
      const x = new Float64Array([1, 1, 1]);

      const yA = A(x);
      const yBlock = block(x);

      expect(arrayClose(yBlock, Array.from(yA))).toBe(true);
    });
  });

  describe('Integration with eigensolvers', () => {
    it('addMatvec result works with eigs', async () => {
      const A = diagMatvec(new Float64Array([1, 2, 3, 4, 5]));
      const B = diagMatvec(new Float64Array([1, 1, 1, 1, 1]));
      const combined = addMatvec(A, B);

      const result = await eigs(combined, 5, 2, { which: 'LM' });

      expect(result.success).toBe(true);
      // Eigenvalues of A+B: 2, 3, 4, 5, 6
      const sorted = Array.from(result.eigenvalues).sort((a, b) => b - a);
      expect(sorted[0]).toBeCloseTo(6, 3);
    });

    it('scaleMatvec and shiftMatvec compose correctly', async () => {
      const A = diagMatvec(new Float64Array([1, 2, 3, 4, 5]));
      // 2A - 3I has eigenvalues: -1, 1, 3, 5, 7
      const scaled = scaleMatvec(A, 2);
      const shifted = shiftMatvec(scaled, 3);

      const result = await eigs(shifted, 5, 2, { which: 'LM' });

      expect(result.success).toBe(true);
      const sorted = Array.from(result.eigenvalues).sort((a, b) => b - a);
      expect(sorted[0]).toBeCloseTo(7, 3);
    });
  });
});
