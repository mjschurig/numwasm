/**
 * Linear Algebra Tests (Phase 13)
 *
 * Tests for numpy.linalg compatible linear algebra operations.
 */

import { describe, it, expect, beforeAll } from 'vitest';
import {
  NDArray,
  loadWasmModule,
  linalg,
  LinAlgError,
  matmul,
  dot,
  vdot,
  inner,
  outer,
  cholesky,
  qr,
  svd,
  svdvals,
  eig,
  eigvals,
  norm,
  det,
  slogdet,
  matrix_rank,
  trace,
  cond,
  solve,
  lstsq,
  inv,
  pinv,
  matrix_power,
  DType,
} from '../../dist/numjs.mjs';

// Helper to check if arrays are close
function allClose(actual: number[], expected: number[], rtol = 1e-5, atol = 1e-8): boolean {
  if (actual.length !== expected.length) return false;
  for (let i = 0; i < actual.length; i++) {
    const diff = Math.abs(actual[i] - expected[i]);
    const tolerance = atol + rtol * Math.abs(expected[i]);
    if (diff > tolerance && !Number.isNaN(expected[i])) return false;
  }
  return true;
}

// Helper to create identity matrix
async function eye(n: number): Promise<NDArray> {
  const data = new Float64Array(n * n);
  for (let i = 0; i < n; i++) {
    data[i * n + i] = 1;
  }
  return NDArray.fromTypedArray(data, [n, n], DType.Float64);
}

describe('Linear Algebra Module', () => {
  beforeAll(async () => {
    await loadWasmModule();
  });

  describe('linalg namespace', () => {
    it('exports all functions', () => {
      expect(linalg).toBeDefined();
      expect(linalg.matmul).toBe(matmul);
      expect(linalg.dot).toBe(dot);
      expect(linalg.solve).toBe(solve);
      expect(linalg.inv).toBe(inv);
      expect(linalg.det).toBe(det);
      expect(linalg.norm).toBe(norm);
      expect(linalg.qr).toBe(qr);
      expect(linalg.svd).toBe(svd);
      expect(linalg.eig).toBe(eig);
      expect(linalg.cholesky).toBe(cholesky);
    });
  });

  describe('LinAlgError', () => {
    it('is a proper Error subclass', () => {
      const err = new LinAlgError('test error');
      expect(err).toBeInstanceOf(Error);
      expect(err).toBeInstanceOf(LinAlgError);
      expect(err.name).toBe('LinAlgError');
      expect(err.message).toBe('test error');
    });
  });

  describe('Matrix multiplication', () => {
    describe('matmul', () => {
      it('multiplies 2x2 matrices', async () => {
        const a = await NDArray.fromTypedArray(
          new Float64Array([1, 2, 3, 4]),
          [2, 2],
          DType.Float64
        );
        const b = await NDArray.fromTypedArray(
          new Float64Array([5, 6, 7, 8]),
          [2, 2],
          DType.Float64
        );
        const c = await matmul(a, b);

        expect(c.shape).toEqual([2, 2]);
        // [[1,2],[3,4]] @ [[5,6],[7,8]] = [[19,22],[43,50]]
        const result = Array.from(await c.toTypedArray());
        expect(allClose(result, [19, 22, 43, 50])).toBe(true);
      });

      it('multiplies non-square matrices', async () => {
        // 2x3 @ 3x2 = 2x2
        const a = await NDArray.fromTypedArray(
          new Float64Array([1, 2, 3, 4, 5, 6]),
          [2, 3],
          DType.Float64
        );
        const b = await NDArray.fromTypedArray(
          new Float64Array([7, 8, 9, 10, 11, 12]),
          [3, 2],
          DType.Float64
        );
        const c = await matmul(a, b);

        expect(c.shape).toEqual([2, 2]);
        // Row 0: 1*7 + 2*9 + 3*11 = 58, 1*8 + 2*10 + 3*12 = 64
        // Row 1: 4*7 + 5*9 + 6*11 = 139, 4*8 + 5*10 + 6*12 = 154
        const result = Array.from(await c.toTypedArray());
        expect(allClose(result, [58, 64, 139, 154])).toBe(true);
      });

      it('throws on incompatible shapes', async () => {
        const a = await NDArray.fromTypedArray(
          new Float64Array([1, 2, 3, 4]),
          [2, 2],
          DType.Float64
        );
        const b = await NDArray.fromTypedArray(
          new Float64Array([1, 2, 3]),
          [3, 1],
          DType.Float64
        );
        await expect(matmul(a, b)).rejects.toThrow();
      });
    });

    describe('dot', () => {
      it('computes dot product of 1D arrays', async () => {
        const a = await NDArray.fromTypedArray(
          new Float64Array([1, 2, 3]),
          [3],
          DType.Float64
        );
        const b = await NDArray.fromTypedArray(
          new Float64Array([4, 5, 6]),
          [3],
          DType.Float64
        );
        const result = await dot(a, b);

        expect(result.shape).toEqual([]);
        // 1*4 + 2*5 + 3*6 = 32
        const data = await result.toTypedArray();
        expect(data[0]).toBeCloseTo(32, 10);
      });

      it('computes matrix dot product (same as matmul)', async () => {
        const a = await NDArray.fromTypedArray(
          new Float64Array([1, 2, 3, 4]),
          [2, 2],
          DType.Float64
        );
        const b = await NDArray.fromTypedArray(
          new Float64Array([5, 6, 7, 8]),
          [2, 2],
          DType.Float64
        );
        const result = await dot(a, b);

        expect(result.shape).toEqual([2, 2]);
        // Same as matmul for 2D arrays
        const data = Array.from(await result.toTypedArray());
        expect(allClose(data, [19, 22, 43, 50])).toBe(true);
      });
    });

    describe('vdot', () => {
      it('computes conjugated dot product for real arrays', async () => {
        const a = await NDArray.fromTypedArray(
          new Float64Array([1, 2, 3]),
          [3],
          DType.Float64
        );
        const b = await NDArray.fromTypedArray(
          new Float64Array([4, 5, 6]),
          [3],
          DType.Float64
        );
        const result = await vdot(a, b);

        // vdot returns a number, not an NDArray
        // 1*4 + 2*5 + 3*6 = 32
        expect(result).toBeCloseTo(32, 10);
      });
    });

    describe('inner', () => {
      it('computes inner product of 1D arrays', async () => {
        const a = await NDArray.fromTypedArray(
          new Float64Array([1, 2, 3]),
          [3],
          DType.Float64
        );
        const b = await NDArray.fromTypedArray(
          new Float64Array([4, 5, 6]),
          [3],
          DType.Float64
        );
        const result = await inner(a, b);

        const data = await result.toTypedArray();
        expect(data[0]).toBeCloseTo(32, 10);
      });
    });

    describe('outer', () => {
      it('computes outer product', async () => {
        const a = await NDArray.fromTypedArray(
          new Float64Array([1, 2, 3]),
          [3],
          DType.Float64
        );
        const b = await NDArray.fromTypedArray(
          new Float64Array([4, 5]),
          [2],
          DType.Float64
        );
        const result = await outer(a, b);

        expect(result.shape).toEqual([3, 2]);
        // [[1*4, 1*5], [2*4, 2*5], [3*4, 3*5]] = [[4,5],[8,10],[12,15]]
        const data = Array.from(await result.toTypedArray());
        expect(allClose(data, [4, 5, 8, 10, 12, 15])).toBe(true);
      });
    });
  });

  describe('Decompositions', () => {
    describe('cholesky', () => {
      it('decomposes positive definite matrix', async () => {
        // A = [[4, 2], [2, 5]] is positive definite
        const a = await NDArray.fromTypedArray(
          new Float64Array([4, 2, 2, 5]),
          [2, 2],
          DType.Float64
        );
        const L = await cholesky(a);

        expect(L.shape).toEqual([2, 2]);
        // L @ L^T should equal A
        const LT = await L.transpose();
        const reconstructed = await matmul(L, LT);
        const origData = Array.from(await a.toTypedArray());
        const reconData = Array.from(await reconstructed.toTypedArray());
        expect(allClose(reconData, origData, 1e-10)).toBe(true);
      });

      it('returns upper triangular with upper=true', async () => {
        const a = await NDArray.fromTypedArray(
          new Float64Array([4, 2, 2, 5]),
          [2, 2],
          DType.Float64
        );
        const U = await cholesky(a, true);

        // U^T @ U should equal A
        const UT = await U.transpose();
        const reconstructed = await matmul(UT, U);
        const origData = Array.from(await a.toTypedArray());
        const reconData = Array.from(await reconstructed.toTypedArray());
        expect(allClose(reconData, origData, 1e-10)).toBe(true);
      });
    });

    describe('qr', () => {
      it('decomposes matrix into Q and R', async () => {
        const a = await NDArray.fromTypedArray(
          new Float64Array([1, 2, 3, 4, 5, 6]),
          [3, 2],
          DType.Float64
        );
        const { Q, R } = await qr(a);

        // Q @ R should equal A
        const reconstructed = await matmul(Q, R);
        const origData = Array.from(await a.toTypedArray());
        const reconData = Array.from(await reconstructed.toTypedArray());
        expect(allClose(reconData, origData, 1e-10)).toBe(true);

        // Q should be orthogonal: Q^T @ Q = I
        const QT = await Q.transpose();
        const QtQ = await matmul(QT, Q);
        const QtQData = Array.from(await QtQ.toTypedArray());
        const I = [1, 0, 0, 1];
        expect(allClose(QtQData, I, 1e-10)).toBe(true);
      });
    });

    describe('svd', () => {
      it('decomposes matrix into U, S, Vh', async () => {
        const a = await NDArray.fromTypedArray(
          new Float64Array([1, 2, 3, 4, 5, 6]),
          [2, 3],
          DType.Float64
        );
        const { U, S, Vh } = await svd(a);

        expect(U.shape[0]).toBe(2);
        expect(S.shape[0]).toBe(2); // min(m, n) singular values
        expect(Vh.shape[1]).toBe(3);

        // Singular values should be non-negative and sorted
        const sData = Array.from(await S.toTypedArray());
        expect(sData[0]).toBeGreaterThan(0);
        expect(sData[0]).toBeGreaterThanOrEqual(sData[1]);
      });
    });

    describe('svdvals', () => {
      it('returns only singular values', async () => {
        const a = await NDArray.fromTypedArray(
          new Float64Array([1, 2, 3, 4, 5, 6]),
          [2, 3],
          DType.Float64
        );
        const s = await svdvals(a);

        expect(s.shape).toEqual([2]);
        const sData = Array.from(await s.toTypedArray());
        expect(sData[0]).toBeGreaterThan(0);
        expect(sData[0]).toBeGreaterThanOrEqual(sData[1]);
      });
    });

    describe('eig', () => {
      it('computes eigenvalues and eigenvectors of symmetric matrix', async () => {
        // Symmetric matrix with known eigenvalues
        const a = await NDArray.fromTypedArray(
          new Float64Array([2, 1, 1, 2]),
          [2, 2],
          DType.Float64
        );
        const { eigenvalues, eigenvectors } = await eig(a);

        expect(eigenvalues.shape).toEqual([2]);
        // Eigenvalues of [[2,1],[1,2]] are 3 and 1
        const evalData = Array.from(await eigenvalues.toTypedArray()).sort((a, b) => b - a);
        expect(evalData[0]).toBeCloseTo(3, 5);
        expect(evalData[1]).toBeCloseTo(1, 5);
      });
    });

    describe('eigvals', () => {
      it('returns only eigenvalues', async () => {
        const a = await NDArray.fromTypedArray(
          new Float64Array([2, 1, 1, 2]),
          [2, 2],
          DType.Float64
        );
        const vals = await eigvals(a);

        expect(vals.shape).toEqual([2]);
        const evalData = Array.from(await vals.toTypedArray()).sort((a, b) => b - a);
        expect(evalData[0]).toBeCloseTo(3, 5);
        expect(evalData[1]).toBeCloseTo(1, 5);
      });
    });
  });

  describe('Norms and related', () => {
    describe('norm', () => {
      it('computes Frobenius norm by default', async () => {
        const a = await NDArray.fromTypedArray(
          new Float64Array([1, 2, 3, 4]),
          [2, 2],
          DType.Float64
        );
        const n = await norm(a);

        // Frobenius norm: sqrt(1 + 4 + 9 + 16) = sqrt(30)
        expect(typeof n).toBe('number');
        expect(n).toBeCloseTo(Math.sqrt(30), 10);
      });

      it('computes L2 norm for vectors', async () => {
        const a = await NDArray.fromTypedArray(
          new Float64Array([3, 4]),
          [2],
          DType.Float64
        );
        const n = await norm(a);

        // L2 norm: sqrt(9 + 16) = 5
        expect(n).toBeCloseTo(5, 10);
      });

      it('computes L1 norm', async () => {
        const a = await NDArray.fromTypedArray(
          new Float64Array([-1, 2, -3]),
          [3],
          DType.Float64
        );
        const n = await norm(a, 1);

        // L1 norm: |(-1)| + |2| + |(-3)| = 6
        expect(n).toBeCloseTo(6, 10);
      });

      it('computes infinity norm', async () => {
        const a = await NDArray.fromTypedArray(
          new Float64Array([-1, 5, -3]),
          [3],
          DType.Float64
        );
        const n = await norm(a, Infinity);

        // Infinity norm: max(|elements|) = 5
        expect(n).toBeCloseTo(5, 10);
      });
    });

    describe('det', () => {
      it('computes determinant of 2x2 matrix', async () => {
        const a = await NDArray.fromTypedArray(
          new Float64Array([1, 2, 3, 4]),
          [2, 2],
          DType.Float64
        );
        const d = await det(a);

        // det([[1,2],[3,4]]) = 1*4 - 2*3 = -2
        expect(d).toBeCloseTo(-2, 10);
      });

      it('computes determinant of 3x3 matrix', async () => {
        const a = await NDArray.fromTypedArray(
          new Float64Array([1, 2, 3, 4, 5, 6, 7, 8, 9]),
          [3, 3],
          DType.Float64
        );
        const d = await det(a);

        // This matrix is singular (rows are linearly dependent)
        expect(d).toBeCloseTo(0, 5);
      });

      it('computes determinant of identity matrix', async () => {
        const I = await eye(3);
        const d = await det(I);
        expect(d).toBeCloseTo(1, 10);
      });
    });

    describe('slogdet', () => {
      it('returns sign and log determinant', async () => {
        const a = await NDArray.fromTypedArray(
          new Float64Array([1, 2, 3, 4]),
          [2, 2],
          DType.Float64
        );
        const { sign, logabsdet } = await slogdet(a);

        // det = -2, so sign = -1, logabsdet = log(2)
        expect(sign).toBe(-1);
        expect(logabsdet).toBeCloseTo(Math.log(2), 10);
      });
    });

    describe('trace', () => {
      it('computes trace of matrix', async () => {
        const a = await NDArray.fromTypedArray(
          new Float64Array([1, 2, 3, 4, 5, 6, 7, 8, 9]),
          [3, 3],
          DType.Float64
        );
        const t = await trace(a);

        // trace = 1 + 5 + 9 = 15
        expect(t).toBeCloseTo(15, 10);
      });
    });

    describe('matrix_rank', () => {
      it('computes rank of full rank matrix', async () => {
        const a = await NDArray.fromTypedArray(
          new Float64Array([1, 2, 3, 4]),
          [2, 2],
          DType.Float64
        );
        const r = await matrix_rank(a);
        expect(r).toBe(2);
      });

      it('computes rank of rank-deficient matrix', async () => {
        // Rows are linearly dependent
        const a = await NDArray.fromTypedArray(
          new Float64Array([1, 2, 2, 4]),
          [2, 2],
          DType.Float64
        );
        const r = await matrix_rank(a);
        expect(r).toBe(1);
      });
    });

    describe('cond', () => {
      it('computes condition number', async () => {
        const a = await NDArray.fromTypedArray(
          new Float64Array([1, 0, 0, 1]),
          [2, 2],
          DType.Float64
        );
        const c = await cond(a);

        // Condition number of identity is 1
        expect(c).toBeCloseTo(1, 5);
      });

      it('returns large value for ill-conditioned matrix', async () => {
        const a = await NDArray.fromTypedArray(
          new Float64Array([1, 0, 0, 1e-10]),
          [2, 2],
          DType.Float64
        );
        const c = await cond(a);

        // Condition number should be very large
        expect(c).toBeGreaterThan(1e9);
      });
    });
  });

  describe('Linear equation solving', () => {
    describe('solve', () => {
      it('solves linear system Ax = b', async () => {
        // A = [[3, 1], [1, 2]], b = [9, 8]
        // Solution: x = [2, 3]
        const A = await NDArray.fromTypedArray(
          new Float64Array([3, 1, 1, 2]),
          [2, 2],
          DType.Float64
        );
        const b = await NDArray.fromTypedArray(
          new Float64Array([9, 8]),
          [2],
          DType.Float64
        );
        const x = await solve(A, b);

        const xData = Array.from(await x.toTypedArray());
        expect(allClose(xData, [2, 3], 1e-10)).toBe(true);
      });

      it('solves system with multiple right-hand sides', async () => {
        const A = await NDArray.fromTypedArray(
          new Float64Array([1, 2, 3, 4]),
          [2, 2],
          DType.Float64
        );
        const B = await NDArray.fromTypedArray(
          new Float64Array([5, 6, 7, 8]),
          [2, 2],
          DType.Float64
        );
        const X = await solve(A, B);

        expect(X.shape).toEqual([2, 2]);
        // Verify A @ X = B
        const check = await matmul(A, X);
        const bData = Array.from(await B.toTypedArray());
        const checkData = Array.from(await check.toTypedArray());
        expect(allClose(checkData, bData, 1e-10)).toBe(true);
      });
    });

    describe('inv', () => {
      it('computes matrix inverse', async () => {
        const a = await NDArray.fromTypedArray(
          new Float64Array([4, 7, 2, 6]),
          [2, 2],
          DType.Float64
        );
        const aInv = await inv(a);

        // A @ A^-1 should be identity
        const product = await matmul(a, aInv);
        const prodData = Array.from(await product.toTypedArray());
        expect(allClose(prodData, [1, 0, 0, 1], 1e-10)).toBe(true);
      });

      it('inverse of identity is identity', async () => {
        const I = await eye(3);
        const Iinv = await inv(I);

        const invData = Array.from(await Iinv.toTypedArray());
        const eyeData = [1, 0, 0, 0, 1, 0, 0, 0, 1];
        expect(allClose(invData, eyeData, 1e-10)).toBe(true);
      });
    });

    describe('pinv', () => {
      it('computes pseudo-inverse', async () => {
        // Test with a square non-singular matrix first
        const a = await NDArray.fromTypedArray(
          new Float64Array([4, 7, 2, 6]),
          [2, 2],
          DType.Float64
        );
        const aPinv = await pinv(a);

        expect(aPinv.shape).toEqual([2, 2]);
        // For non-singular square matrices, pinv = inv
        // A @ A+ should be identity
        const product = await matmul(a, aPinv);
        const prodData = Array.from(await product.toTypedArray());
        expect(allClose(prodData, [1, 0, 0, 1], 1e-5)).toBe(true);
      });

      it('computes pseudo-inverse of rectangular matrix', async () => {
        const a = await NDArray.fromTypedArray(
          new Float64Array([1, 2, 3, 4, 5, 6]),
          [2, 3],
          DType.Float64
        );
        const aPinv = await pinv(a);

        expect(aPinv.shape).toEqual([3, 2]);
        // Check that A @ A+ @ A ≈ A (Moore-Penrose condition)
        const aaPinv = await matmul(a, aPinv);
        const result = await matmul(aaPinv, a);
        const origData = Array.from(await a.toTypedArray());
        const resultData = Array.from(await result.toTypedArray());
        // Use more relaxed tolerance for numerical stability
        expect(allClose(resultData, origData, 1e-3, 1e-6)).toBe(true);
      });
    });

    describe('lstsq', () => {
      it('solves overdetermined system', async () => {
        // Overdetermined system (more equations than unknowns)
        const A = await NDArray.fromTypedArray(
          new Float64Array([1, 1, 1, 2, 1, 3]),
          [3, 2],
          DType.Float64
        );
        const b = await NDArray.fromTypedArray(
          new Float64Array([1, 2, 2]),
          [3],
          DType.Float64
        );
        const { x, residuals, rank, s } = await lstsq(A, b);

        expect(x.shape).toEqual([2]);
        expect(rank).toBe(2);
        expect(s.shape[0]).toBe(2);
      });
    });
  });

  describe('Matrix powers', () => {
    describe('matrix_power', () => {
      it('computes positive integer power', async () => {
        const a = await NDArray.fromTypedArray(
          new Float64Array([1, 2, 3, 4]),
          [2, 2],
          DType.Float64
        );
        const a2 = await matrix_power(a, 2);

        // A^2 = A @ A
        const expected = await matmul(a, a);
        const a2Data = Array.from(await a2.toTypedArray());
        const expectedData = Array.from(await expected.toTypedArray());
        expect(allClose(a2Data, expectedData, 1e-10)).toBe(true);
      });

      it('power of 0 returns identity', async () => {
        const a = await NDArray.fromTypedArray(
          new Float64Array([1, 2, 3, 4]),
          [2, 2],
          DType.Float64
        );
        const a0 = await matrix_power(a, 0);

        const a0Data = Array.from(await a0.toTypedArray());
        expect(allClose(a0Data, [1, 0, 0, 1], 1e-10)).toBe(true);
      });

      it('power of 1 returns copy', async () => {
        const a = await NDArray.fromTypedArray(
          new Float64Array([1, 2, 3, 4]),
          [2, 2],
          DType.Float64
        );
        const a1 = await matrix_power(a, 1);

        const origData = Array.from(await a.toTypedArray());
        const a1Data = Array.from(await a1.toTypedArray());
        expect(allClose(a1Data, origData, 1e-10)).toBe(true);
      });

      it('negative power uses inverse', async () => {
        const a = await NDArray.fromTypedArray(
          new Float64Array([4, 7, 2, 6]),
          [2, 2],
          DType.Float64
        );
        const aNeg1 = await matrix_power(a, -1);

        // A^-1 should be the inverse
        const expected = await inv(a);
        const aNeg1Data = Array.from(await aNeg1.toTypedArray());
        const expectedData = Array.from(await expected.toTypedArray());
        expect(allClose(aNeg1Data, expectedData, 1e-10)).toBe(true);
      });
    });
  });
});
