/**
 * Tests for Phase 25: Advanced Linear Algebra
 */
import { describe, it, expect, beforeAll } from 'vitest';
import {
  NDArray,
  loadWasmModule,
  tensordot,
  multi_dot,
  kron,
  cross,
  tensorsolve,
  tensorinv,
  matrix_norm,
  vector_norm,
  einsum,
  einsum_path,
  matmul,
  dot,
  DType,
  linalg,
} from '../../dist/numjs.mjs';

/**
 * Helper for approximate equality with tolerance.
 */
function allClose(
  actual: number[],
  expected: number[],
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

describe('Advanced Linear Algebra (Phase 25)', () => {
  beforeAll(async () => {
    await loadWasmModule();
  });

  describe('linalg namespace', () => {
    it('exports all Phase 25 functions', () => {
      expect(linalg.tensordot).toBe(tensordot);
      expect(linalg.multi_dot).toBe(multi_dot);
      expect(linalg.kron).toBe(kron);
      expect(linalg.cross).toBe(cross);
      expect(linalg.tensorsolve).toBe(tensorsolve);
      expect(linalg.tensorinv).toBe(tensorinv);
      expect(linalg.matrix_norm).toBe(matrix_norm);
      expect(linalg.vector_norm).toBe(vector_norm);
    });
  });

  describe('tensordot', () => {
    it('with axes=1 equals matmul for 2D arrays', async () => {
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

      const c1 = await tensordot(a, b, 1);
      const c2 = await matmul(a, b);

      expect(c1.shape).toEqual([2, 2]);
      const r1 = Array.from(await c1.toTypedArray());
      const r2 = Array.from(await c2.toTypedArray());
      expect(allClose(r1, r2)).toBe(true);
    });

    it('with axes=0 gives outer product shape', async () => {
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

      const c = await tensordot(a, b, 0);
      expect(c.shape).toEqual([3, 2]);

      const result = Array.from(await c.toTypedArray());
      // outer product: [[1*4, 1*5], [2*4, 2*5], [3*4, 3*5]]
      expect(allClose(result, [4, 5, 8, 10, 12, 15])).toBe(true);
    });

    it('with explicit axes works correctly', async () => {
      const a = await NDArray.fromTypedArray(
        new Float64Array([1, 2, 3, 4, 5, 6]),
        [2, 3],
        DType.Float64
      );
      const b = await NDArray.fromTypedArray(
        new Float64Array([1, 2, 3, 4, 5, 6]),
        [3, 2],
        DType.Float64
      );

      // Contract a's axis 1 with b's axis 0
      const c = await tensordot(a, b, [[1], [0]]);
      expect(c.shape).toEqual([2, 2]);
    });

    it('throws on incompatible shapes', async () => {
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

      await expect(tensordot(a, b, 1)).rejects.toThrow();
    });
  });

  describe('multi_dot', () => {
    it('equals sequential dot product', async () => {
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
      const c = await NDArray.fromTypedArray(
        new Float64Array([9, 10, 11, 12]),
        [2, 2],
        DType.Float64
      );

      const result1 = await multi_dot([a, b, c]);
      const ab = await dot(a, b);
      const result2 = await dot(ab, c);

      const r1 = Array.from(await result1.toTypedArray());
      const r2 = Array.from(await result2.toTypedArray());
      expect(allClose(r1, r2)).toBe(true);
    });

    it('handles 1D first and last arrays', async () => {
      const v1 = await NDArray.fromTypedArray(
        new Float64Array([1, 2]),
        [2],
        DType.Float64
      );
      const m = await NDArray.fromTypedArray(
        new Float64Array([1, 2, 3, 4]),
        [2, 2],
        DType.Float64
      );
      const v2 = await NDArray.fromTypedArray(
        new Float64Array([1, 1]),
        [2],
        DType.Float64
      );

      const result = await multi_dot([v1, m, v2]);
      // v1 @ m @ v2 = [1,2] @ [[1,2],[3,4]] @ [1,1] = [7, 10] @ [1,1] = 17
      expect(result.shape).toEqual([]);
      expect(result.item()).toBeCloseTo(17);
    });
  });

  describe('kron', () => {
    it('produces correct shape', async () => {
      const a = await NDArray.fromTypedArray(
        new Float64Array([1, 2, 3, 4]),
        [2, 2],
        DType.Float64
      );
      const b = await NDArray.fromTypedArray(
        new Float64Array([1, 0, 0, 1]),
        [2, 2],
        DType.Float64
      );

      const c = await kron(a, b);
      expect(c.shape).toEqual([4, 4]);
    });

    it('computes correct values for 1D arrays', async () => {
      const a = await NDArray.fromTypedArray(
        new Float64Array([1, 10, 100]),
        [3],
        DType.Float64
      );
      const b = await NDArray.fromTypedArray(
        new Float64Array([5, 6, 7]),
        [3],
        DType.Float64
      );

      const c = await kron(a, b);
      expect(c.shape).toEqual([9]);

      const result = Array.from(await c.toTypedArray());
      // [1*5, 1*6, 1*7, 10*5, 10*6, 10*7, 100*5, 100*6, 100*7]
      expect(allClose(result, [5, 6, 7, 50, 60, 70, 500, 600, 700])).toBe(true);
    });
  });

  describe('cross', () => {
    it('computes 3D cross product correctly', async () => {
      const a = await NDArray.fromTypedArray(
        new Float64Array([1, 0, 0]),
        [3],
        DType.Float64
      );
      const b = await NDArray.fromTypedArray(
        new Float64Array([0, 1, 0]),
        [3],
        DType.Float64
      );

      const c = await cross(a, b);
      expect(c.shape).toEqual([3]);

      const result = Array.from(await c.toTypedArray());
      // [1, 0, 0] × [0, 1, 0] = [0, 0, 1]
      expect(allClose(result, [0, 0, 1])).toBe(true);
    });

    it('handles batched cross products', async () => {
      const a = await NDArray.fromTypedArray(
        new Float64Array([1, 0, 0, 0, 1, 0]),
        [2, 3],
        DType.Float64
      );
      const b = await NDArray.fromTypedArray(
        new Float64Array([0, 1, 0, 0, 0, 1]),
        [2, 3],
        DType.Float64
      );

      const c = await cross(a, b);
      expect(c.shape).toEqual([2, 3]);

      const result = Array.from(await c.toTypedArray());
      // First: [1,0,0] × [0,1,0] = [0,0,1]
      // Second: [0,1,0] × [0,0,1] = [1,0,0]
      expect(allClose(result, [0, 0, 1, 1, 0, 0])).toBe(true);
    });
  });

  describe('tensorsolve', () => {
    it('solves simple tensor equation', async () => {
      // Create identity-like 4D tensor
      const eye4 = await NDArray.eye(4);
      const a = eye4.reshape([2, 2, 2, 2]);
      const b = await NDArray.fromTypedArray(
        new Float64Array([1, 2, 3, 4]),
        [2, 2],
        DType.Float64
      );

      const x = await tensorsolve(a, b);
      expect(x.shape).toEqual([2, 2]);

      // For identity, x should equal b
      const xData = Array.from(await x.toTypedArray());
      const bData = Array.from(await b.toTypedArray());
      expect(allClose(xData, bData)).toBe(true);
    });
  });

  describe('tensorinv', () => {
    it('computes tensor inverse', async () => {
      const eye4 = await NDArray.eye(4);
      const a = eye4.reshape([2, 2, 2, 2]);

      const ainv = await tensorinv(a, 2);
      expect(ainv.shape).toEqual([2, 2, 2, 2]);

      // For identity tensor, inverse is also identity
      const result = Array.from(await ainv.toTypedArray());
      const expected = Array.from(await a.toTypedArray());
      expect(allClose(result, expected)).toBe(true);
    });

    it('throws on non-positive ind', async () => {
      const a = await NDArray.fromTypedArray(
        new Float64Array([1, 2, 3, 4]),
        [2, 2],
        DType.Float64
      );
      await expect(tensorinv(a, 0)).rejects.toThrow();
    });
  });

  describe('matrix_norm', () => {
    it('computes Frobenius norm by default', async () => {
      const a = await NDArray.fromTypedArray(
        new Float64Array([1, 2, 3, 4]),
        [2, 2],
        DType.Float64
      );

      const n = await matrix_norm(a);
      // Frobenius norm = sqrt(1 + 4 + 9 + 16) = sqrt(30)
      expect(n).toBeCloseTo(Math.sqrt(30), 10);
    });

    it('handles keepdims', async () => {
      const a = await NDArray.fromTypedArray(
        new Float64Array([1, 2, 3, 4]),
        [2, 2],
        DType.Float64
      );

      const n = (await matrix_norm(a, 'fro', true)) as NDArray;
      expect(n.shape).toEqual([1, 1]);
    });
  });

  describe('vector_norm', () => {
    it('computes L2 norm correctly', async () => {
      const v = await NDArray.fromTypedArray(
        new Float64Array([3, 4]),
        [2],
        DType.Float64
      );

      const n = await vector_norm(v);
      expect(n).toBeCloseTo(5, 10);
    });

    it('computes L1 norm correctly', async () => {
      const v = await NDArray.fromTypedArray(
        new Float64Array([3, -4]),
        [2],
        DType.Float64
      );

      const n = await vector_norm(v, 1);
      expect(n).toBeCloseTo(7, 10);
    });

    it('computes Linf norm correctly', async () => {
      const v = await NDArray.fromTypedArray(
        new Float64Array([3, -4, 2]),
        [3],
        DType.Float64
      );

      const n = await vector_norm(v, Infinity);
      expect(n).toBeCloseTo(4, 10);
    });

    it('computes L0 norm (count non-zero)', async () => {
      const v = await NDArray.fromTypedArray(
        new Float64Array([3, 0, 2, 0, 1]),
        [5],
        DType.Float64
      );

      const n = await vector_norm(v, 0);
      expect(n).toBe(3);
    });

    it('handles axis parameter', async () => {
      const m = await NDArray.fromTypedArray(
        new Float64Array([3, 4, 5, 12]),
        [2, 2],
        DType.Float64
      );

      const n = (await vector_norm(m, 2, 1)) as NDArray;
      expect(n.shape).toEqual([2]);

      const result = Array.from(await n.toTypedArray());
      // Row 0: sqrt(9+16) = 5
      // Row 1: sqrt(25+144) = 13
      expect(allClose(result, [5, 13])).toBe(true);
    });
  });

  describe('einsum', () => {
    it('computes matrix multiplication with ij,jk->ik', async () => {
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

      const c = await einsum('ij,jk->ik', a, b);
      const expected = await matmul(a, b);

      const r1 = Array.from(await c.toTypedArray());
      const r2 = Array.from(await expected.toTypedArray());
      expect(allClose(r1, r2)).toBe(true);
    });

    it('computes trace with ii->', async () => {
      const a = await NDArray.fromTypedArray(
        new Float64Array([1, 2, 3, 4]),
        [2, 2],
        DType.Float64
      );

      const tr = await einsum('ii->', a);
      expect(tr.shape).toEqual([]);
      expect(tr.item()).toBeCloseTo(5, 10); // 1 + 4 = 5
    });

    it('computes transpose with ij->ji', async () => {
      const a = await NDArray.fromTypedArray(
        new Float64Array([1, 2, 3, 4]),
        [2, 2],
        DType.Float64
      );

      const t = await einsum('ij->ji', a);
      expect(t.shape).toEqual([2, 2]);

      const result = Array.from(await t.toTypedArray());
      expect(allClose(result, [1, 3, 2, 4])).toBe(true);
    });

    it('computes outer product with i,j->ij', async () => {
      const a = await NDArray.fromTypedArray(
        new Float64Array([1, 2]),
        [2],
        DType.Float64
      );
      const b = await NDArray.fromTypedArray(
        new Float64Array([3, 4, 5]),
        [3],
        DType.Float64
      );

      const c = await einsum('i,j->ij', a, b);
      expect(c.shape).toEqual([2, 3]);

      const result = Array.from(await c.toTypedArray());
      expect(allClose(result, [3, 4, 5, 6, 8, 10])).toBe(true);
    });

    it('computes sum with ij->i', async () => {
      const a = await NDArray.fromTypedArray(
        new Float64Array([1, 2, 3, 4, 5, 6]),
        [2, 3],
        DType.Float64
      );

      const s = await einsum('ij->i', a);
      expect(s.shape).toEqual([2]);

      const result = Array.from(await s.toTypedArray());
      expect(allClose(result, [6, 15])).toBe(true); // 1+2+3, 4+5+6
    });
  });

  describe('einsum_path', () => {
    it('returns valid path structure', async () => {
      const a = await NDArray.fromTypedArray(
        new Float64Array(Array(6).fill(1)),
        [2, 3],
        DType.Float64
      );
      const b = await NDArray.fromTypedArray(
        new Float64Array(Array(12).fill(1)),
        [3, 4],
        DType.Float64
      );
      const c = await NDArray.fromTypedArray(
        new Float64Array(Array(20).fill(1)),
        [4, 5],
        DType.Float64
      );

      const [path, info] = einsum_path('ij,jk,kl->il', a, b, c);

      expect(Array.isArray(path)).toBe(true);
      expect(path.length).toBe(2); // Two contractions for 3 operands
      expect(typeof info).toBe('string');
      expect(info).toContain('Contraction path');
    });
  });
});
