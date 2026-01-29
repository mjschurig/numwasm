/**
 * Tests for DIAMatrix (Diagonal format)
 */
import { describe, it, expect, beforeAll } from 'vitest';
import { sparse, loadWasmModule } from '../../../src/ts/index.js';

const { DIAMatrix, dia_matrix } = sparse;

beforeAll(async () => {
  await loadWasmModule();
});

describe('DIAMatrix', () => {
  it('constructs from offsets and data arrays', () => {
    // Tridiagonal matrix:
    // [ 2, -1,  0,  0]
    // [-1,  2, -1,  0]
    // [ 0, -1,  2, -1]
    // [ 0,  0, -1,  2]
    const m = new DIAMatrix(
      {
        offsets: [-1, 0, 1],
        data: [
          [-1, -1, -1, 0],  // subdiagonal (offset -1)
          [2, 2, 2, 2],     // main diagonal (offset 0)
          [0, -1, -1, -1],  // superdiagonal (offset 1)
        ],
      },
      { shape: [4, 4] }
    );
    expect(m.shape).toEqual([4, 4]);
    expect(m.format).toBe('dia');
    expect(m.ndiags).toBe(3);
  });

  it('constructs from dense 2D array', () => {
    const dense = [
      [1, 2, 0],
      [3, 4, 5],
      [0, 6, 7],
    ];
    const m = new DIAMatrix(dense);
    expect(m.shape).toEqual([3, 3]);
    expect(m.format).toBe('dia');
  });

  it('toArray() roundtrips from dense', () => {
    const dense = [
      [1, 2, 0],
      [3, 4, 5],
      [0, 6, 7],
    ];
    const m = new DIAMatrix(dense);
    expect(m.toArray()).toEqual(dense);
  });

  it('toArray() tridiagonal matrix', () => {
    const m = new DIAMatrix(
      {
        offsets: [-1, 0, 1],
        data: [
          [-1, -1, -1, 0],
          [2, 2, 2, 2],
          [0, -1, -1, -1],
        ],
      },
      { shape: [4, 4] }
    );
    expect(m.toArray()).toEqual([
      [2, -1, 0, 0],
      [-1, 2, -1, 0],
      [0, -1, 2, -1],
      [0, 0, -1, 2],
    ]);
  });

  it('dot() matrix-vector multiply for diagonal', () => {
    // Identity matrix via DIA
    const m = new DIAMatrix(
      {
        offsets: [0],
        data: [[1, 1, 1]],
      },
      { shape: [3, 3] }
    );
    const x = new Float64Array([1, 2, 3]);
    const y = m.dot(x);
    expect(Array.from(y)).toEqual([1, 2, 3]);
  });

  it('dot() matrix-vector multiply for tridiagonal', () => {
    // Tridiagonal:
    // [ 2, -1,  0]
    // [-1,  2, -1]
    // [ 0, -1,  2]
    const m = new DIAMatrix(
      {
        offsets: [-1, 0, 1],
        data: [
          [-1, -1, 0],
          [2, 2, 2],
          [0, -1, -1],
        ],
      },
      { shape: [3, 3] }
    );
    // x = [1, 2, 3]
    // y[0] = 2*1 + (-1)*2 = 0
    // y[1] = (-1)*1 + 2*2 + (-1)*3 = 0
    // y[2] = (-1)*2 + 2*3 = 4
    const x = new Float64Array([1, 2, 3]);
    const y = m.dot(x);
    expect(Array.from(y)).toEqual([0, 0, 4]);
  });

  it('tocsr() conversion', () => {
    const dense = [
      [1, 2, 0],
      [3, 4, 5],
      [0, 6, 7],
    ];
    const m = new DIAMatrix(dense);
    const csr = m.tocsr();
    expect(csr.format).toBe('csr');
    expect(csr.toArray()).toEqual(dense);
  });

  it('tocoo() conversion', () => {
    const dense = [
      [1, 2, 0],
      [3, 4, 5],
      [0, 6, 7],
    ];
    const m = new DIAMatrix(dense);
    const coo = m.tocoo();
    expect(coo.format).toBe('coo');
    expect(coo.toArray()).toEqual(dense);
  });

  it('tocsc() conversion (via CSR)', () => {
    const dense = [
      [1, 2, 0],
      [3, 4, 5],
      [0, 6, 7],
    ];
    const m = new DIAMatrix(dense);
    const csc = m.tocsc();
    expect(csc.format).toBe('csc');
    expect(csc.toArray()).toEqual(dense);
  });

  it('diagonal() extracts stored diagonal', () => {
    const m = new DIAMatrix(
      {
        offsets: [-1, 0, 1],
        data: [
          [1, 2, 0],
          [3, 4, 5],
          [0, 6, 7],
        ],
      },
      { shape: [3, 3] }
    );
    const diag = m.diagonal(0);
    expect(Array.from(diag)).toEqual([3, 4, 5]);
  });

  it('diagonal() returns zeros for unstored diagonal', () => {
    const m = new DIAMatrix(
      {
        offsets: [0],  // only main diagonal
        data: [[1, 2, 3]],
      },
      { shape: [3, 3] }
    );
    const diag2 = m.diagonal(2);
    expect(Array.from(diag2)).toEqual([0]);
  });

  it('copy() creates independent copy', () => {
    const m1 = new DIAMatrix(
      {
        offsets: [0],
        data: [[1, 2, 3]],
      },
      { shape: [3, 3] }
    );
    const m2 = m1.copy();
    expect(m2.toArray()).toEqual(m1.toArray());
  });

  it('transpose() swaps shape and negates offsets', () => {
    const dense = [
      [1, 2, 0],
      [3, 4, 5],
      [0, 6, 7],
    ];
    const m = new DIAMatrix(dense);
    const mt = m.transpose();
    expect(mt.shape).toEqual([3, 3]);
    expect(mt.toArray()).toEqual([
      [1, 3, 0],
      [2, 4, 6],
      [0, 5, 7],
    ]);
  });

  it('transpose() rectangular matrix', () => {
    const dense = [
      [1, 2, 0, 0],
      [3, 4, 5, 0],
    ];
    const m = new DIAMatrix(dense);
    const mt = m.transpose();
    expect(mt.shape).toEqual([4, 2]);
    expect(mt.toArray()).toEqual([
      [1, 3],
      [2, 4],
      [0, 5],
      [0, 0],
    ]);
  });

  it('handles empty (zero) matrix', () => {
    const m = new DIAMatrix(
      {
        offsets: [],
        data: [],
      },
      { shape: [3, 3] }
    );
    expect(m.nnz).toBe(0);
    expect(m.toArray()).toEqual([
      [0, 0, 0],
      [0, 0, 0],
      [0, 0, 0],
    ]);
  });

  it('handles single diagonal (identity)', () => {
    const m = new DIAMatrix(
      {
        offsets: [0],
        data: [[1, 1, 1, 1]],
      },
      { shape: [4, 4] }
    );
    expect(m.toArray()).toEqual([
      [1, 0, 0, 0],
      [0, 1, 0, 0],
      [0, 0, 1, 0],
      [0, 0, 0, 1],
    ]);
  });

  it('handles rectangular matrices (more rows)', () => {
    const m = new DIAMatrix(
      {
        offsets: [0],
        data: [[1, 2, 3, 4]],
      },
      { shape: [4, 3] }
    );
    expect(m.toArray()).toEqual([
      [1, 0, 0],
      [0, 2, 0],
      [0, 0, 3],
      [0, 0, 0],
    ]);
  });

  it('handles rectangular matrices (more cols)', () => {
    const m = new DIAMatrix(
      {
        offsets: [0],
        data: [[1, 2, 3, 4]],
      },
      { shape: [3, 4] }
    );
    expect(m.toArray()).toEqual([
      [1, 0, 0, 0],
      [0, 2, 0, 0],
      [0, 0, 3, 0],
    ]);
  });

  it('add() delegates to CSR', () => {
    const m1 = new DIAMatrix([[1, 0], [0, 2]]);
    const m2 = new DIAMatrix([[3, 0], [0, 4]]);
    const result = m1.add(m2);
    expect(result.toArray()).toEqual([
      [4, 0],
      [0, 6],
    ]);
  });

  it('subtract() delegates to CSR', () => {
    const m1 = new DIAMatrix([[5, 0], [0, 8]]);
    const m2 = new DIAMatrix([[2, 0], [0, 3]]);
    const result = m1.subtract(m2);
    expect(result.toArray()).toEqual([
      [3, 0],
      [0, 5],
    ]);
  });

  it('multiply() element-wise delegates to CSR', () => {
    const m1 = new DIAMatrix([[2, 0], [0, 3]]);
    const m2 = new DIAMatrix([[4, 0], [0, 5]]);
    const result = m1.multiply(m2);
    expect(result.toArray()).toEqual([
      [8, 0],
      [0, 15],
    ]);
  });

  it('matmul() matrix multiply delegates to CSR', () => {
    const m1 = new DIAMatrix([[1, 2], [3, 4]]);
    const m2 = new DIAMatrix([[5, 6], [7, 8]]);
    const result = m1.matmul(m2);
    expect(result.toArray()).toEqual([
      [19, 22],
      [43, 50],
    ]);
  });

  it('factory function dia_matrix() works', () => {
    const m = dia_matrix([[1, 0], [0, 2]]);
    expect(m.format).toBe('dia');
    expect(m.shape).toEqual([2, 2]);
  });

  it('offsets getter returns copy', () => {
    const m = new DIAMatrix(
      {
        offsets: [-1, 0, 1],
        data: [[1, 2, 0], [3, 4, 5], [0, 6, 7]],
      },
      { shape: [3, 3] }
    );
    const offsets = m.offsets;
    expect(Array.from(offsets)).toEqual([-1, 0, 1]);
  });

  it('nnz counts actual non-zeros', () => {
    const m = new DIAMatrix(
      {
        offsets: [-1, 0, 1],
        data: [
          [-1, -1, 0],  // 2 non-zeros in valid range
          [2, 2, 2],    // 3 non-zeros
          [0, -1, -1],  // 2 non-zeros in valid range
        ],
      },
      { shape: [3, 3] }
    );
    expect(m.nnz).toBe(7);
  });
});
