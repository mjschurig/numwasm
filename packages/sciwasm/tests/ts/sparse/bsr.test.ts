/**
 * Tests for BSRMatrix (Block Sparse Row format)
 */
import { describe, it, expect, beforeAll } from 'vitest';
import { sparse, loadWasmModule } from '../../../src/ts/index.js';

const { BSRMatrix, bsr_matrix } = sparse;

beforeAll(async () => {
  await loadWasmModule();
});

describe('BSRMatrix', () => {
  it('constructs from BSR arrays with blocksize', () => {
    // 4x4 matrix with 2x2 blocks:
    // [[1, 2, 0, 0],
    //  [3, 4, 0, 0],
    //  [0, 0, 5, 6],
    //  [0, 0, 7, 8]]
    const m = new BSRMatrix(
      {
        data: new Float64Array([1, 2, 3, 4, 5, 6, 7, 8]),  // 2 blocks, each 2x2
        indices: new Int32Array([0, 1]),  // block columns
        indptr: new Int32Array([0, 1, 2]),  // block row pointers
        blocksize: [2, 2],
      },
      { shape: [4, 4] }
    );
    expect(m.shape).toEqual([4, 4]);
    expect(m.blocksize).toEqual([2, 2]);
    expect(m.format).toBe('bsr');
    expect(m.nnzb).toBe(2);
  });

  it('constructs from dense 2D array with blocksize', () => {
    const dense = [
      [1, 2, 0, 0],
      [3, 4, 0, 0],
      [0, 0, 5, 6],
      [0, 0, 7, 8],
    ];
    const m = new BSRMatrix(dense, { blocksize: [2, 2] });
    expect(m.shape).toEqual([4, 4]);
    expect(m.blocksize).toEqual([2, 2]);
    expect(m.nnzb).toBe(2);
  });

  it('toArray() roundtrips from dense', () => {
    const dense = [
      [1, 2, 0, 0],
      [3, 4, 0, 0],
      [0, 0, 5, 6],
      [0, 0, 7, 8],
    ];
    const m = new BSRMatrix(dense, { blocksize: [2, 2] });
    expect(m.toArray()).toEqual(dense);
  });

  it('dot() matrix-vector multiply', () => {
    // 4x4 diagonal blocks:
    // [[1, 2, 0, 0],
    //  [3, 4, 0, 0],
    //  [0, 0, 5, 6],
    //  [0, 0, 7, 8]]
    // x = [1, 1, 1, 1]
    // y[0] = 1*1 + 2*1 = 3
    // y[1] = 3*1 + 4*1 = 7
    // y[2] = 5*1 + 6*1 = 11
    // y[3] = 7*1 + 8*1 = 15
    const dense = [
      [1, 2, 0, 0],
      [3, 4, 0, 0],
      [0, 0, 5, 6],
      [0, 0, 7, 8],
    ];
    const m = new BSRMatrix(dense, { blocksize: [2, 2] });
    const x = new Float64Array([1, 1, 1, 1]);
    const y = m.dot(x);
    expect(Array.from(y)).toEqual([3, 7, 11, 15]);
  });

  it('dot() with general vector', () => {
    const dense = [
      [1, 2, 0, 0],
      [3, 4, 0, 0],
      [0, 0, 5, 6],
      [0, 0, 7, 8],
    ];
    const m = new BSRMatrix(dense, { blocksize: [2, 2] });
    // x = [1, 2, 3, 4]
    // y[0] = 1*1 + 2*2 = 5
    // y[1] = 3*1 + 4*2 = 11
    // y[2] = 5*3 + 6*4 = 39
    // y[3] = 7*3 + 8*4 = 53
    const x = new Float64Array([1, 2, 3, 4]);
    const y = m.dot(x);
    expect(Array.from(y)).toEqual([5, 11, 39, 53]);
  });

  it('tocsr() conversion', () => {
    const dense = [
      [1, 2, 0, 0],
      [3, 4, 0, 0],
      [0, 0, 5, 6],
      [0, 0, 7, 8],
    ];
    const m = new BSRMatrix(dense, { blocksize: [2, 2] });
    const csr = m.tocsr();
    expect(csr.format).toBe('csr');
    expect(csr.toArray()).toEqual(dense);
  });

  it('tocoo() conversion', () => {
    const dense = [
      [1, 2, 0, 0],
      [3, 4, 0, 0],
      [0, 0, 5, 6],
      [0, 0, 7, 8],
    ];
    const m = new BSRMatrix(dense, { blocksize: [2, 2] });
    const coo = m.tocoo();
    expect(coo.format).toBe('coo');
    expect(coo.toArray()).toEqual(dense);
  });

  it('tocsc() conversion (via CSR)', () => {
    const dense = [
      [1, 2, 0, 0],
      [3, 4, 0, 0],
      [0, 0, 5, 6],
      [0, 0, 7, 8],
    ];
    const m = new BSRMatrix(dense, { blocksize: [2, 2] });
    const csc = m.tocsc();
    expect(csc.format).toBe('csc');
    expect(csc.toArray()).toEqual(dense);
  });

  it('diagonal() extracts main diagonal', () => {
    const dense = [
      [1, 2, 0, 0],
      [3, 4, 0, 0],
      [0, 0, 5, 6],
      [0, 0, 7, 8],
    ];
    const m = new BSRMatrix(dense, { blocksize: [2, 2] });
    const diag = m.diagonal(0);
    expect(Array.from(diag)).toEqual([1, 4, 5, 8]);
  });

  it('diagonal() extracts superdiagonal', () => {
    const dense = [
      [1, 2, 0, 0],
      [3, 4, 0, 0],
      [0, 0, 5, 6],
      [0, 0, 7, 8],
    ];
    const m = new BSRMatrix(dense, { blocksize: [2, 2] });
    const diag = m.diagonal(1);
    expect(Array.from(diag)).toEqual([2, 0, 6]);
  });

  it('transpose() swaps dimensions and transposes blocks', () => {
    const dense = [
      [1, 2, 0, 0],
      [3, 4, 0, 0],
      [0, 0, 5, 6],
      [0, 0, 7, 8],
    ];
    const m = new BSRMatrix(dense, { blocksize: [2, 2] });
    const mt = m.transpose();
    expect(mt.shape).toEqual([4, 4]);
    expect(mt.blocksize).toEqual([2, 2]);
    expect(mt.toArray()).toEqual([
      [1, 3, 0, 0],
      [2, 4, 0, 0],
      [0, 0, 5, 7],
      [0, 0, 6, 8],
    ]);
  });

  it('copy() creates independent copy', () => {
    const dense = [
      [1, 2, 0, 0],
      [3, 4, 0, 0],
      [0, 0, 5, 6],
      [0, 0, 7, 8],
    ];
    const m1 = new BSRMatrix(dense, { blocksize: [2, 2] });
    const m2 = m1.copy();
    expect(m2.toArray()).toEqual(m1.toArray());
    expect(m2.blocksize).toEqual(m1.blocksize);
  });

  it('handles single block', () => {
    const dense = [
      [1, 2],
      [3, 4],
    ];
    const m = new BSRMatrix(dense, { blocksize: [2, 2] });
    expect(m.nnzb).toBe(1);
    expect(m.toArray()).toEqual(dense);
  });

  it('handles empty blocks', () => {
    const dense = [
      [0, 0, 0, 0],
      [0, 0, 0, 0],
      [0, 0, 1, 2],
      [0, 0, 3, 4],
    ];
    const m = new BSRMatrix(dense, { blocksize: [2, 2] });
    expect(m.nnzb).toBe(1);  // only one non-zero block
    expect(m.toArray()).toEqual(dense);
  });

  it('handles rectangular matrix', () => {
    const dense = [
      [1, 2, 3, 4, 0, 0],
      [5, 6, 7, 8, 0, 0],
    ];
    const m = new BSRMatrix(dense, { blocksize: [2, 2] });
    expect(m.shape).toEqual([2, 6]);
    expect(m.nbrows).toBe(1);
    expect(m.nbcols).toBe(3);
    expect(m.toArray()).toEqual(dense);
  });

  it('add() delegates to CSR', () => {
    const m1 = new BSRMatrix([[1, 0], [0, 2]], { blocksize: [1, 1] });
    const m2 = new BSRMatrix([[3, 0], [0, 4]], { blocksize: [1, 1] });
    const result = m1.add(m2);
    expect(result.toArray()).toEqual([
      [4, 0],
      [0, 6],
    ]);
  });

  it('subtract() delegates to CSR', () => {
    const m1 = new BSRMatrix([[5, 0], [0, 8]], { blocksize: [1, 1] });
    const m2 = new BSRMatrix([[2, 0], [0, 3]], { blocksize: [1, 1] });
    const result = m1.subtract(m2);
    expect(result.toArray()).toEqual([
      [3, 0],
      [0, 5],
    ]);
  });

  it('multiply() element-wise delegates to CSR', () => {
    const m1 = new BSRMatrix([[2, 0], [0, 3]], { blocksize: [1, 1] });
    const m2 = new BSRMatrix([[4, 0], [0, 5]], { blocksize: [1, 1] });
    const result = m1.multiply(m2);
    expect(result.toArray()).toEqual([
      [8, 0],
      [0, 15],
    ]);
  });

  it('matmul() matrix multiply delegates to CSR', () => {
    const m1 = new BSRMatrix([[1, 2], [3, 4]], { blocksize: [1, 1] });
    const m2 = new BSRMatrix([[5, 6], [7, 8]], { blocksize: [1, 1] });
    const result = m1.matmul(m2);
    expect(result.toArray()).toEqual([
      [19, 22],
      [43, 50],
    ]);
  });

  it('factory function bsr_matrix() works', () => {
    const m = bsr_matrix([[1, 0], [0, 2]], { blocksize: [1, 1] });
    expect(m.format).toBe('bsr');
    expect(m.shape).toEqual([2, 2]);
  });

  it('throws on shape not divisible by blocksize', () => {
    expect(() => {
      new BSRMatrix([[1, 2, 3], [4, 5, 6]], { blocksize: [2, 2] });
    }).toThrow('not divisible by blocksize');
  });

  it('nnz counts actual non-zeros in blocks', () => {
    // Blocks may contain zeros
    const m = new BSRMatrix(
      {
        data: new Float64Array([1, 0, 0, 2]),  // block with 2 zeros
        indices: new Int32Array([0]),
        indptr: new Int32Array([0, 1]),
        blocksize: [2, 2],
      },
      { shape: [2, 2] }
    );
    expect(m.nnzb).toBe(1);  // 1 block
    expect(m.nnz).toBe(2);   // 2 actual non-zeros
  });

  it('nbrows and nbcols are correct', () => {
    const m = new BSRMatrix(
      {
        data: new Float64Array([1, 2, 3, 4]),
        indices: new Int32Array([0]),
        indptr: new Int32Array([0, 1, 1, 1]),  // 3 block rows, only first has data
        blocksize: [2, 2],
      },
      { shape: [6, 4] }
    );
    expect(m.nbrows).toBe(3);
    expect(m.nbcols).toBe(2);
  });
});
