/**
 * Tests for CSRMatrix
 *
 * Ported from scipy/sparse/tests/test_base.py
 */
import { describe, it, expect, beforeAll } from 'vitest';
import { sparse, loadWasmModule } from '../../../src/ts/index.js';

const { CSRMatrix, csr_matrix } = sparse;

beforeAll(async () => {
  await loadWasmModule();
});

describe('CSRMatrix', () => {
  it('constructs from (data, indices, indptr) + shape', () => {
    const m = new CSRMatrix(
      {
        data: new Float64Array([1, 2, 3, 4, 5, 6]),
        indices: new Int32Array([0, 2, 2, 0, 1, 2]),
        indptr: new Int32Array([0, 2, 3, 6]),
      },
      { shape: [3, 3] }
    );
    expect(m.shape).toEqual([3, 3]);
    expect(m.nnz).toBe(6);
    expect(m.format).toBe('csr');
  });

  it('constructs from dense 2D array', () => {
    const dense = [
      [1, 0, 2],
      [0, 0, 3],
      [4, 5, 6],
    ];
    const m = new CSRMatrix(dense, { shape: [3, 3] });
    expect(m.shape).toEqual([3, 3]);
    expect(m.nnz).toBe(6);
    expect(m.format).toBe('csr');
  });

  it('toArray() roundtrips from dense', () => {
    const dense = [
      [1, 0, 2],
      [0, 0, 3],
      [4, 5, 6],
    ];
    const m = new CSRMatrix(dense, { shape: [3, 3] });
    expect(m.toArray()).toEqual(dense);
  });

  it('toArray() roundtrips from (data, indices, indptr)', () => {
    // [ [1, 0, 2],
    //   [0, 0, 3],
    //   [4, 5, 6] ]
    const m = new CSRMatrix(
      {
        data: new Float64Array([1, 2, 3, 4, 5, 6]),
        indices: new Int32Array([0, 2, 2, 0, 1, 2]),
        indptr: new Int32Array([0, 2, 3, 6]),
      },
      { shape: [3, 3] }
    );
    expect(m.toArray()).toEqual([
      [1, 0, 2],
      [0, 0, 3],
      [4, 5, 6],
    ]);
  });

  it('tocsc() and back roundtrips', () => {
    const dense = [
      [1, 0, 2],
      [0, 0, 3],
      [4, 5, 6],
    ];
    const m = new CSRMatrix(dense, { shape: [3, 3] });
    const csc = m.tocsc();
    expect(csc.format).toBe('csc');
    expect(csc.shape).toEqual([3, 3]);
    const back = csc.tocsr();
    expect(back.format).toBe('csr');
    expect(back.toArray()).toEqual(dense);
  });

  it('diagonal(0) extracts main diagonal', () => {
    const m = new CSRMatrix(
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

  it('diagonal(k) extracts off-diagonals', () => {
    const m = new CSRMatrix(
      [
        [1, 2, 3],
        [4, 5, 6],
        [7, 8, 9],
      ],
      { shape: [3, 3] }
    );
    expect(Array.from(m.diagonal(1))).toEqual([2, 6]);
    expect(Array.from(m.diagonal(-1))).toEqual([4, 8]);
  });

  it('dot() performs matrix-vector multiply', () => {
    const m = new CSRMatrix(
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

  it('dot() with non-diagonal matrix', () => {
    const m = new CSRMatrix(
      [
        [1, 2],
        [3, 4],
      ],
      { shape: [2, 2] }
    );
    const x = new Float64Array([1, 1]);
    const y = m.dot(x);
    expect(Array.from(y)).toEqual([3, 7]);
  });

  it('add() sums two sparse matrices', () => {
    const a = new CSRMatrix(
      [
        [1, 0],
        [0, 2],
      ],
      { shape: [2, 2] }
    );
    const b = new CSRMatrix(
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

  it('subtract() differences two sparse matrices', () => {
    const a = new CSRMatrix(
      [
        [5, 3],
        [1, 4],
      ],
      { shape: [2, 2] }
    );
    const b = new CSRMatrix(
      [
        [1, 1],
        [1, 1],
      ],
      { shape: [2, 2] }
    );
    const c = a.subtract(b);
    expect(c.toArray()).toEqual([
      [4, 2],
      [0, 3],
    ]);
  });

  it('multiply() element-wise multiplies', () => {
    const a = new CSRMatrix(
      [
        [1, 2],
        [3, 4],
      ],
      { shape: [2, 2] }
    );
    const b = new CSRMatrix(
      [
        [2, 0],
        [0, 3],
      ],
      { shape: [2, 2] }
    );
    const c = a.multiply(b);
    expect(c.toArray()).toEqual([
      [2, 0],
      [0, 12],
    ]);
  });

  it('matmul() sparse * sparse', () => {
    const a = new CSRMatrix(
      [
        [1, 2],
        [3, 4],
      ],
      { shape: [2, 2] }
    );
    const b = new CSRMatrix(
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

  it('matmul() non-square', () => {
    const a = new CSRMatrix(
      [
        [1, 0, 2],
        [0, 3, 0],
      ],
      { shape: [2, 3] }
    );
    const b = new CSRMatrix(
      [
        [1, 0],
        [0, 1],
        [2, 0],
      ],
      { shape: [3, 2] }
    );
    const c = a.matmul(b);
    expect(c.shape).toEqual([2, 2]);
    expect(c.toArray()).toEqual([
      [5, 0],
      [0, 3],
    ]);
  });

  it('copy() creates independent copy', () => {
    const m = new CSRMatrix(
      [
        [1, 2],
        [3, 4],
      ],
      { shape: [2, 2] }
    );
    const c = m.copy();
    expect(c.toArray()).toEqual(m.toArray());
    expect(c).not.toBe(m);
  });

  it('transpose() returns CSC', () => {
    const m = new CSRMatrix(
      [
        [1, 2],
        [3, 4],
      ],
      { shape: [2, 2] }
    );
    const t = m.transpose();
    expect(t.format).toBe('csc');
    expect(t.shape).toEqual([2, 2]);
    expect(t.toArray()).toEqual([
      [1, 3],
      [2, 4],
    ]);
  });

  it('csr_matrix() factory function works', () => {
    const m = csr_matrix(
      [
        [1, 0],
        [0, 1],
      ],
      { shape: [2, 2] }
    );
    expect(m.format).toBe('csr');
    expect(m.toArray()).toEqual([
      [1, 0],
      [0, 1],
    ]);
  });

  it('sortIndices() sorts column indices', () => {
    // Construct with unsorted indices
    const m = new CSRMatrix(
      {
        data: new Float64Array([2, 1, 3]),
        indices: new Int32Array([1, 0, 2]),
        indptr: new Int32Array([0, 2, 3]),
      },
      { shape: [2, 3] }
    );
    m.sortIndices();
    // After sorting, indices should be [0, 1, 2] and data reordered
    expect(Array.from(m.indices)).toEqual([0, 1, 2]);
    expect(Array.from(m.data)).toEqual([1, 2, 3]);
  });

  it('sumDuplicates() merges duplicate entries', () => {
    const m = new CSRMatrix(
      {
        data: new Float64Array([1, 2, 3]),
        indices: new Int32Array([0, 0, 1]),
        indptr: new Int32Array([0, 2, 3]),
      },
      { shape: [2, 2] }
    );
    m.sumDuplicates();
    expect(m.toArray()).toEqual([
      [3, 0],
      [0, 3],
    ]);
  });

  it('eliminateZeros() removes zero entries', () => {
    const m = new CSRMatrix(
      {
        data: new Float64Array([1, 0, 3]),
        indices: new Int32Array([0, 1, 2]),
        indptr: new Int32Array([0, 2, 3]),
      },
      { shape: [2, 3] }
    );
    m.eliminateZeros();
    expect(m.nnz).toBe(2);
    expect(m.toArray()).toEqual([
      [1, 0, 0],
      [0, 0, 3],
    ]);
  });

  it('handles empty matrix', () => {
    const m = new CSRMatrix(
      {
        data: new Float64Array(0),
        indices: new Int32Array(0),
        indptr: new Int32Array([0, 0, 0]),
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
