/**
 * Tests for LILMatrix (List of Lists format)
 */
import { describe, it, expect, beforeAll } from 'vitest';
import { sparse, loadWasmModule } from '../../../src/ts/index.js';

const { LILMatrix, lil_matrix } = sparse;

beforeAll(async () => {
  await loadWasmModule();
});

describe('LILMatrix', () => {
  it('constructs from shape tuple (empty matrix)', () => {
    const m = new LILMatrix([3, 3]);
    expect(m.shape).toEqual([3, 3]);
    expect(m.nnz).toBe(0);
    expect(m.format).toBe('lil');
  });

  it('constructs from dense 2D array', () => {
    const dense = [
      [1, 0, 2],
      [0, 0, 3],
      [4, 5, 6],
    ];
    const m = new LILMatrix(dense);
    expect(m.shape).toEqual([3, 3]);
    expect(m.nnz).toBe(6);
    expect(m.format).toBe('lil');
  });

  it('toArray() roundtrips from dense', () => {
    const dense = [
      [1, 0, 2],
      [0, 0, 3],
      [4, 5, 6],
    ];
    const m = new LILMatrix(dense);
    expect(m.toArray()).toEqual(dense);
  });

  it('get() retrieves correct values', () => {
    const dense = [
      [1, 0, 2],
      [0, 0, 3],
      [4, 5, 6],
    ];
    const m = new LILMatrix(dense);
    expect(m.get(0, 0)).toBe(1);
    expect(m.get(0, 1)).toBe(0);
    expect(m.get(0, 2)).toBe(2);
    expect(m.get(2, 2)).toBe(6);
  });

  it('get() supports negative indices', () => {
    const dense = [
      [1, 0, 2],
      [0, 0, 3],
      [4, 5, 6],
    ];
    const m = new LILMatrix(dense);
    expect(m.get(-1, -1)).toBe(6);  // last row, last col
    expect(m.get(-3, 0)).toBe(1);   // first row
  });

  it('set() updates values', () => {
    const m = new LILMatrix([3, 3]);
    m.set(0, 0, 1);
    m.set(1, 1, 2);
    m.set(2, 2, 3);
    expect(m.get(0, 0)).toBe(1);
    expect(m.get(1, 1)).toBe(2);
    expect(m.get(2, 2)).toBe(3);
    expect(m.nnz).toBe(3);
  });

  it('set() maintains sorted column order', () => {
    const m = new LILMatrix([1, 5]);
    // Insert out of order (use non-zero values)
    m.set(0, 3, 30);
    m.set(0, 1, 10);
    m.set(0, 4, 40);
    m.set(0, 0, 5);
    m.set(0, 2, 20);

    const row = m.getrow(0);
    expect(row.indices).toEqual([0, 1, 2, 3, 4]);
    expect(row.data).toEqual([5, 10, 20, 30, 40]);
  });

  it('set() with zero removes entry', () => {
    const m = new LILMatrix([[1, 0], [0, 2]]);
    expect(m.nnz).toBe(2);
    m.set(0, 0, 0);
    expect(m.nnz).toBe(1);
    expect(m.get(0, 0)).toBe(0);
  });

  it('getrow() returns row data', () => {
    const m = new LILMatrix([
      [1, 0, 2],
      [0, 0, 3],
      [4, 5, 6],
    ]);
    const row0 = m.getrow(0);
    expect(row0.indices).toEqual([0, 2]);
    expect(row0.data).toEqual([1, 2]);

    const row2 = m.getrow(2);
    expect(row2.indices).toEqual([0, 1, 2]);
    expect(row2.data).toEqual([4, 5, 6]);
  });

  it('getrow() supports negative index', () => {
    const m = new LILMatrix([
      [1, 0, 2],
      [0, 0, 3],
      [4, 5, 6],
    ]);
    const lastRow = m.getrow(-1);
    expect(lastRow.indices).toEqual([0, 1, 2]);
    expect(lastRow.data).toEqual([4, 5, 6]);
  });

  it('tocsr() conversion (efficient)', () => {
    const m = new LILMatrix([
      [1, 0, 2],
      [0, 0, 3],
      [4, 5, 6],
    ]);
    const csr = m.tocsr();
    expect(csr.format).toBe('csr');
    expect(csr.shape).toEqual([3, 3]);
    expect(csr.toArray()).toEqual([
      [1, 0, 2],
      [0, 0, 3],
      [4, 5, 6],
    ]);
  });

  it('tocoo() conversion', () => {
    const m = new LILMatrix([
      [1, 0, 2],
      [0, 0, 3],
      [4, 5, 6],
    ]);
    const coo = m.tocoo();
    expect(coo.format).toBe('coo');
    expect(coo.toArray()).toEqual([
      [1, 0, 2],
      [0, 0, 3],
      [4, 5, 6],
    ]);
  });

  it('tocsc() conversion (via CSR)', () => {
    const m = new LILMatrix([
      [1, 0, 2],
      [0, 0, 3],
      [4, 5, 6],
    ]);
    const csc = m.tocsc();
    expect(csc.format).toBe('csc');
    expect(csc.toArray()).toEqual([
      [1, 0, 2],
      [0, 0, 3],
      [4, 5, 6],
    ]);
  });

  it('diagonal() extracts main diagonal', () => {
    const m = new LILMatrix([
      [1, 0, 2],
      [0, 3, 0],
      [4, 5, 6],
    ]);
    const diag = m.diagonal(0);
    expect(Array.from(diag)).toEqual([1, 3, 6]);
  });

  it('diagonal() extracts superdiagonal', () => {
    const m = new LILMatrix([
      [1, 2, 0],
      [0, 3, 4],
      [0, 0, 5],
    ]);
    const diag = m.diagonal(1);
    expect(Array.from(diag)).toEqual([2, 4]);
  });

  it('diagonal() extracts subdiagonal', () => {
    const m = new LILMatrix([
      [1, 0, 0],
      [2, 3, 0],
      [0, 4, 5],
    ]);
    const diag = m.diagonal(-1);
    expect(Array.from(diag)).toEqual([2, 4]);
  });

  it('copy() creates independent copy', () => {
    const m1 = new LILMatrix([[1, 0], [0, 2]]);
    const m2 = m1.copy();
    m2.set(0, 0, 999);
    expect(m1.get(0, 0)).toBe(1);
    expect(m2.get(0, 0)).toBe(999);
  });

  it('transpose() swaps rows and cols', () => {
    const m = new LILMatrix([
      [1, 2, 3],
      [4, 5, 6],
    ]);
    const mt = m.transpose();
    expect(mt.shape).toEqual([3, 2]);
    expect(mt.toArray()).toEqual([
      [1, 4],
      [2, 5],
      [3, 6],
    ]);
  });

  it('handles empty matrix', () => {
    const m = new LILMatrix([5, 5]);
    expect(m.nnz).toBe(0);
    expect(m.toArray()).toEqual([
      [0, 0, 0, 0, 0],
      [0, 0, 0, 0, 0],
      [0, 0, 0, 0, 0],
      [0, 0, 0, 0, 0],
      [0, 0, 0, 0, 0],
    ]);
  });

  it('handles rectangular matrices', () => {
    const m = new LILMatrix([[1, 2, 3, 4], [5, 6, 7, 8]]);
    expect(m.shape).toEqual([2, 4]);
    expect(m.toArray()).toEqual([
      [1, 2, 3, 4],
      [5, 6, 7, 8],
    ]);
  });

  it('add() delegates to CSR', () => {
    const m1 = new LILMatrix([[1, 0], [0, 2]]);
    const m2 = new LILMatrix([[3, 0], [0, 4]]);
    const result = m1.add(m2);
    expect(result.toArray()).toEqual([
      [4, 0],
      [0, 6],
    ]);
  });

  it('subtract() delegates to CSR', () => {
    const m1 = new LILMatrix([[5, 0], [0, 8]]);
    const m2 = new LILMatrix([[2, 0], [0, 3]]);
    const result = m1.subtract(m2);
    expect(result.toArray()).toEqual([
      [3, 0],
      [0, 5],
    ]);
  });

  it('multiply() element-wise delegates to CSR', () => {
    const m1 = new LILMatrix([[2, 0], [0, 3]]);
    const m2 = new LILMatrix([[4, 0], [0, 5]]);
    const result = m1.multiply(m2);
    expect(result.toArray()).toEqual([
      [8, 0],
      [0, 15],
    ]);
  });

  it('matmul() matrix multiply delegates to CSR', () => {
    const m1 = new LILMatrix([[1, 2], [3, 4]]);
    const m2 = new LILMatrix([[5, 6], [7, 8]]);
    const result = m1.matmul(m2);
    expect(result.toArray()).toEqual([
      [19, 22],
      [43, 50],
    ]);
  });

  it('dot() matrix-vector multiply delegates to CSR', () => {
    const m = new LILMatrix([
      [1, 0, 2],
      [0, 0, 3],
      [4, 5, 6],
    ]);
    const x = new Float64Array([1, 2, 3]);
    const y = m.dot(x);
    expect(Array.from(y)).toEqual([7, 9, 32]);
  });

  it('factory function lil_matrix() works', () => {
    const m = lil_matrix([[1, 0], [0, 2]]);
    expect(m.format).toBe('lil');
    expect(m.shape).toEqual([2, 2]);
  });

  it('factory function lil_matrix() with shape', () => {
    const m = lil_matrix([4, 4]);
    expect(m.format).toBe('lil');
    expect(m.shape).toEqual([4, 4]);
    expect(m.nnz).toBe(0);
  });

  it('throws on out of bounds access', () => {
    const m = new LILMatrix([2, 2]);
    expect(() => m.get(5, 0)).toThrow();
    expect(() => m.set(0, 5, 1)).toThrow();
  });

  it('throws on out of bounds getrow', () => {
    const m = new LILMatrix([2, 2]);
    expect(() => m.getrow(5)).toThrow();
  });
});
