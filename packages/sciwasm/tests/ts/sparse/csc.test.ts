/**
 * Tests for CSCMatrix
 *
 * Ported from scipy/sparse/tests/test_base.py
 */
import { describe, it, expect, beforeAll } from 'vitest';
import { sparse, loadWasmModule } from '../../../src/ts/index.js';

const { CSCMatrix, csc_matrix } = sparse;

beforeAll(async () => {
  await loadWasmModule();
});

describe('CSCMatrix', () => {
  it('constructs from (data, indices, indptr) + shape', () => {
    // Column-compressed: indptr has ncol+1 entries, indices are row indices
    const m = new CSCMatrix(
      {
        data: new Float64Array([1, 4, 5, 2, 3, 6]),
        indices: new Int32Array([0, 2, 2, 0, 1, 2]),
        indptr: new Int32Array([0, 2, 3, 6]),
      },
      { shape: [3, 3] }
    );
    expect(m.shape).toEqual([3, 3]);
    expect(m.nnz).toBe(6);
    expect(m.format).toBe('csc');
  });

  it('constructs from dense 2D array', () => {
    const dense = [
      [1, 0, 2],
      [0, 0, 3],
      [4, 5, 6],
    ];
    const m = new CSCMatrix(dense, { shape: [3, 3] });
    expect(m.shape).toEqual([3, 3]);
    expect(m.nnz).toBe(6);
    expect(m.format).toBe('csc');
  });

  it('toArray() roundtrips from dense', () => {
    const dense = [
      [1, 0, 2],
      [0, 0, 3],
      [4, 5, 6],
    ];
    const m = new CSCMatrix(dense, { shape: [3, 3] });
    expect(m.toArray()).toEqual(dense);
  });

  it('tocsr() and back roundtrips', () => {
    const dense = [
      [1, 0, 2],
      [0, 0, 3],
      [4, 5, 6],
    ];
    const m = new CSCMatrix(dense, { shape: [3, 3] });
    const csr = m.tocsr();
    expect(csr.format).toBe('csr');
    expect(csr.shape).toEqual([3, 3]);
    const back = csr.tocsc();
    expect(back.format).toBe('csc');
    expect(back.toArray()).toEqual(dense);
  });

  it('diagonal(0) extracts main diagonal', () => {
    const m = new CSCMatrix(
      [
        [1, 0, 0],
        [0, 2, 0],
        [0, 0, 3],
      ],
      { shape: [3, 3] }
    );
    const d = m.diagonal(0);
    expect(Array.from(d)).toEqual([1, 2, 3]);
  });

  it('dot() performs matrix-vector multiply', () => {
    const m = new CSCMatrix(
      [
        [1, 0, 0],
        [0, 2, 0],
        [0, 0, 3],
      ],
      { shape: [3, 3] }
    );
    const x = new Float64Array([1, 2, 3]);
    const y = m.dot(x);
    expect(Array.from(y)).toEqual([1, 4, 9]);
  });

  it('add() sums two CSC matrices', () => {
    const a = new CSCMatrix(
      [
        [1, 0],
        [0, 2],
      ],
      { shape: [2, 2] }
    );
    const b = new CSCMatrix(
      [
        [0, 3],
        [4, 0],
      ],
      { shape: [2, 2] }
    );
    const c = a.add(b);
    expect(c.toArray()).toEqual([
      [1, 3],
      [4, 2],
    ]);
  });

  it('matmul() sparse * sparse', () => {
    const a = new CSCMatrix(
      [
        [1, 2],
        [3, 4],
      ],
      { shape: [2, 2] }
    );
    const b = new CSCMatrix(
      [
        [5, 6],
        [7, 8],
      ],
      { shape: [2, 2] }
    );
    const c = a.matmul(b);
    expect(c.toArray()).toEqual([
      [19, 22],
      [43, 50],
    ]);
  });

  it('copy() creates independent copy', () => {
    const m = new CSCMatrix(
      [
        [1, 2],
        [3, 4],
      ],
      { shape: [2, 2] }
    );
    const c = m.copy();
    expect(c.toArray()).toEqual(m.toArray());
    expect(c.format).toBe('csc');
    expect(c).not.toBe(m);
  });

  it('transpose() returns CSR', () => {
    const m = new CSCMatrix(
      [
        [1, 2],
        [3, 4],
      ],
      { shape: [2, 2] }
    );
    const t = m.transpose();
    expect(t.format).toBe('csr');
    expect(t.shape).toEqual([2, 2]);
    expect(t.toArray()).toEqual([
      [1, 3],
      [2, 4],
    ]);
  });

  it('csc_matrix() factory function works', () => {
    const m = csc_matrix(
      [
        [1, 0],
        [0, 1],
      ],
      { shape: [2, 2] }
    );
    expect(m.format).toBe('csc');
    expect(m.toArray()).toEqual([
      [1, 0],
      [0, 1],
    ]);
  });

  it('handles empty matrix', () => {
    const m = new CSCMatrix(
      {
        data: new Float64Array(0),
        indices: new Int32Array(0),
        indptr: new Int32Array([0, 0, 0, 0]),
      },
      { shape: [2, 3] }
    );
    expect(m.nnz).toBe(0);
    expect(m.toArray()).toEqual([
      [0, 0, 0],
      [0, 0, 0],
    ]);
  });
});
