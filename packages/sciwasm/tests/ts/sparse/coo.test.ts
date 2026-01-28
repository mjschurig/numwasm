/**
 * Tests for COOMatrix (Coordinate format)
 *
 * Ported from scipy/sparse/tests/test_coo.py
 */
import { describe, it, expect, beforeAll } from 'vitest';
import { sparse, loadWasmModule } from '../../../src/ts/index.js';

const { COOMatrix, coo_matrix } = sparse;

beforeAll(async () => {
  await loadWasmModule();
});

describe('COOMatrix', () => {
  it('constructs from (row, col, data) + shape', () => {
    const m = new COOMatrix(
      {
        row: new Int32Array([0, 0, 1, 2, 2, 2]),
        col: new Int32Array([0, 2, 2, 0, 1, 2]),
        data: new Float64Array([1, 2, 3, 4, 5, 6]),
      },
      { shape: [3, 3] }
    );
    expect(m.shape).toEqual([3, 3]);
    expect(m.nnz).toBe(6);
    expect(m.format).toBe('coo');
  });

  it('constructs from dense 2D array', () => {
    const dense = [
      [1, 0, 2],
      [0, 0, 3],
      [4, 5, 6],
    ];
    const m = new COOMatrix(dense, { shape: [3, 3] });
    expect(m.shape).toEqual([3, 3]);
    expect(m.nnz).toBe(6);
    expect(m.format).toBe('coo');
  });

  it('toArray() roundtrips from dense', () => {
    const dense = [
      [1, 0, 2],
      [0, 0, 3],
      [4, 5, 6],
    ];
    const m = new COOMatrix(dense, { shape: [3, 3] });
    expect(m.toArray()).toEqual(dense);
  });

  it('toArray() roundtrips from (row, col, data)', () => {
    // [ [1, 0, 2],
    //   [0, 0, 3],
    //   [4, 5, 6] ]
    const m = new COOMatrix(
      {
        row: new Int32Array([0, 0, 1, 2, 2, 2]),
        col: new Int32Array([0, 2, 2, 0, 1, 2]),
        data: new Float64Array([1, 2, 3, 4, 5, 6]),
      },
      { shape: [3, 3] }
    );
    expect(m.toArray()).toEqual([
      [1, 0, 2],
      [0, 0, 3],
      [4, 5, 6],
    ]);
  });

  it('tocsr() conversion', () => {
    const m = new COOMatrix(
      {
        row: new Int32Array([0, 0, 1, 2, 2, 2]),
        col: new Int32Array([0, 2, 2, 0, 1, 2]),
        data: new Float64Array([1, 2, 3, 4, 5, 6]),
      },
      { shape: [3, 3] }
    );
    const csr = m.tocsr();
    expect(csr.format).toBe('csr');
    expect(csr.shape).toEqual([3, 3]);
    expect(csr.nnz).toBe(6);
    expect(csr.toArray()).toEqual([
      [1, 0, 2],
      [0, 0, 3],
      [4, 5, 6],
    ]);
  });

  it('tocsc() conversion (via CSR)', () => {
    const m = new COOMatrix(
      {
        row: new Int32Array([0, 0, 1, 2, 2, 2]),
        col: new Int32Array([0, 2, 2, 0, 1, 2]),
        data: new Float64Array([1, 2, 3, 4, 5, 6]),
      },
      { shape: [3, 3] }
    );
    const csc = m.tocsc();
    expect(csc.format).toBe('csc');
    expect(csc.shape).toEqual([3, 3]);
    expect(csc.toArray()).toEqual([
      [1, 0, 2],
      [0, 0, 3],
      [4, 5, 6],
    ]);
  });

  it('roundtrip COO -> CSR -> COO preserves values', () => {
    const dense = [
      [1, 0, 2],
      [0, 0, 3],
      [4, 5, 6],
    ];
    const m1 = new COOMatrix(dense, { shape: [3, 3] });
    const csr = m1.tocsr();
    const m2 = new COOMatrix(csr.toArray(), { shape: [3, 3] });
    expect(m2.toArray()).toEqual(dense);
  });

  it('dot() matrix-vector multiply', () => {
    // Matrix:
    // [ [1, 0, 2],
    //   [0, 0, 3],
    //   [4, 5, 6] ]
    // Vector: [1, 2, 3]
    // Result: [1*1 + 0*2 + 2*3, 0*1 + 0*2 + 3*3, 4*1 + 5*2 + 6*3] = [7, 9, 32]
    const m = new COOMatrix(
      {
        row: new Int32Array([0, 0, 1, 2, 2, 2]),
        col: new Int32Array([0, 2, 2, 0, 1, 2]),
        data: new Float64Array([1, 2, 3, 4, 5, 6]),
      },
      { shape: [3, 3] }
    );
    const x = new Float64Array([1, 2, 3]);
    const y = m.dot(x);
    expect(Array.from(y)).toEqual([7, 9, 32]);
  });

  it('transpose() swaps row and col', () => {
    const m = new COOMatrix(
      {
        row: new Int32Array([0, 0, 1, 2, 2, 2]),
        col: new Int32Array([0, 2, 2, 0, 1, 2]),
        data: new Float64Array([1, 2, 3, 4, 5, 6]),
      },
      { shape: [3, 3] }
    );
    const mt = m.transpose();
    expect(mt.shape).toEqual([3, 3]);
    expect(mt.format).toBe('coo');

    // Original matrix:
    // [ [1, 0, 2],
    //   [0, 0, 3],
    //   [4, 5, 6] ]
    // Transpose:
    // [ [1, 0, 4],
    //   [0, 0, 5],
    //   [2, 3, 6] ]
    expect(mt.toArray()).toEqual([
      [1, 0, 4],
      [0, 0, 5],
      [2, 3, 6],
    ]);
  });

  it('copy() creates independent copy', () => {
    const m1 = new COOMatrix(
      {
        row: new Int32Array([0, 1, 2]),
        col: new Int32Array([0, 1, 2]),
        data: new Float64Array([1, 2, 3]),
      },
      { shape: [3, 3] }
    );
    const m2 = m1.copy();

    // Modify m2's data
    m2.data[0] = 999;

    // m1 should be unchanged
    expect(m1.data[0]).toBe(1);
    expect(m2.data[0]).toBe(999);
  });

  it('handles empty matrix', () => {
    const m = new COOMatrix(
      {
        row: new Int32Array([]),
        col: new Int32Array([]),
        data: new Float64Array([]),
      },
      { shape: [5, 5] }
    );
    expect(m.nnz).toBe(0);
    expect(m.shape).toEqual([5, 5]);
    expect(m.toArray()).toEqual([
      [0, 0, 0, 0, 0],
      [0, 0, 0, 0, 0],
      [0, 0, 0, 0, 0],
      [0, 0, 0, 0, 0],
      [0, 0, 0, 0, 0],
    ]);
  });

  it('handles single element', () => {
    const m = new COOMatrix(
      {
        row: new Int32Array([1]),
        col: new Int32Array([2]),
        data: new Float64Array([42]),
      },
      { shape: [3, 3] }
    );
    expect(m.nnz).toBe(1);
    expect(m.toArray()).toEqual([
      [0, 0, 0],
      [0, 0, 42],
      [0, 0, 0],
    ]);
  });

  it('handles rectangular matrices', () => {
    const m = new COOMatrix(
      {
        row: new Int32Array([0, 1, 1]),
        col: new Int32Array([0, 1, 3]),
        data: new Float64Array([1, 2, 3]),
      },
      { shape: [2, 4] }
    );
    expect(m.shape).toEqual([2, 4]);
    expect(m.toArray()).toEqual([
      [1, 0, 0, 0],
      [0, 2, 0, 3],
    ]);
  });

  it('factory function coo_matrix() works', () => {
    const m = coo_matrix(
      {
        row: new Int32Array([0, 1, 2]),
        col: new Int32Array([0, 1, 2]),
        data: new Float64Array([1, 2, 3]),
      },
      { shape: [3, 3] }
    );
    expect(m.format).toBe('coo');
    expect(m.shape).toEqual([3, 3]);
  });

  it('exposes row, col, data arrays', () => {
    const m = new COOMatrix(
      {
        row: new Int32Array([0, 1, 2]),
        col: new Int32Array([0, 1, 2]),
        data: new Float64Array([1, 2, 3]),
      },
      { shape: [3, 3] }
    );
    expect(m.row).toEqual(new Int32Array([0, 1, 2]));
    expect(m.col).toEqual(new Int32Array([0, 1, 2]));
    expect(m.data).toEqual(new Float64Array([1, 2, 3]));
  });

  it('validates array lengths match', () => {
    expect(() => {
      new COOMatrix(
        {
          row: new Int32Array([0, 1]),       // length 2
          col: new Int32Array([0, 1, 2]),    // length 3 - mismatch!
          data: new Float64Array([1, 2]),
        },
        { shape: [3, 3] }
      );
    }).toThrow('row, col, and data arrays must have the same length');
  });

  it('requires shape when constructing from arrays', () => {
    expect(() => {
      new COOMatrix({
        row: new Int32Array([0]),
        col: new Int32Array([0]),
        data: new Float64Array([1]),
      });
    }).toThrow('shape is required');
  });

  it('add() delegates to CSR', () => {
    const m1 = new COOMatrix(
      [[1, 0], [0, 2]],
      { shape: [2, 2] }
    );
    const m2 = new COOMatrix(
      [[3, 0], [0, 4]],
      { shape: [2, 2] }
    );
    const result = m1.add(m2);
    expect(result.toArray()).toEqual([
      [4, 0],
      [0, 6],
    ]);
  });

  it('subtract() delegates to CSR', () => {
    const m1 = new COOMatrix(
      [[5, 0], [0, 8]],
      { shape: [2, 2] }
    );
    const m2 = new COOMatrix(
      [[2, 0], [0, 3]],
      { shape: [2, 2] }
    );
    const result = m1.subtract(m2);
    expect(result.toArray()).toEqual([
      [3, 0],
      [0, 5],
    ]);
  });

  it('multiply() element-wise delegates to CSR', () => {
    const m1 = new COOMatrix(
      [[2, 0], [0, 3]],
      { shape: [2, 2] }
    );
    const m2 = new COOMatrix(
      [[4, 0], [0, 5]],
      { shape: [2, 2] }
    );
    const result = m1.multiply(m2);
    expect(result.toArray()).toEqual([
      [8, 0],
      [0, 15],
    ]);
  });

  it('matmul() matrix multiply delegates to CSR', () => {
    const m1 = new COOMatrix(
      [[1, 2], [3, 4]],
      { shape: [2, 2] }
    );
    const m2 = new COOMatrix(
      [[5, 6], [7, 8]],
      { shape: [2, 2] }
    );
    // [1, 2]   [5, 6]   [1*5 + 2*7, 1*6 + 2*8]   [19, 22]
    // [3, 4] × [7, 8] = [3*5 + 4*7, 3*6 + 4*8] = [43, 50]
    const result = m1.matmul(m2);
    expect(result.toArray()).toEqual([
      [19, 22],
      [43, 50],
    ]);
  });

  it('diagonal() delegates to CSR', () => {
    const m = new COOMatrix(
      [[1, 0, 0], [0, 2, 0], [0, 0, 3]],
      { shape: [3, 3] }
    );
    const diag = m.diagonal(0);
    expect(Array.from(diag)).toEqual([1, 2, 3]);
  });
});
