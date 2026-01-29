/**
 * Tests for DOKMatrix (Dictionary of Keys format)
 */
import { describe, it, expect, beforeAll } from 'vitest';
import { sparse, loadWasmModule } from '../../../src/ts/index.js';

const { DOKMatrix, dok_matrix } = sparse;

beforeAll(async () => {
  await loadWasmModule();
});

describe('DOKMatrix', () => {
  it('constructs from shape tuple (empty matrix)', () => {
    const m = new DOKMatrix([3, 3]);
    expect(m.shape).toEqual([3, 3]);
    expect(m.nnz).toBe(0);
    expect(m.format).toBe('dok');
  });

  it('constructs from dense 2D array', () => {
    const dense = [
      [1, 0, 2],
      [0, 0, 3],
      [4, 5, 6],
    ];
    const m = new DOKMatrix(dense);
    expect(m.shape).toEqual([3, 3]);
    expect(m.nnz).toBe(6);
    expect(m.format).toBe('dok');
  });

  it('toArray() roundtrips from dense', () => {
    const dense = [
      [1, 0, 2],
      [0, 0, 3],
      [4, 5, 6],
    ];
    const m = new DOKMatrix(dense);
    expect(m.toArray()).toEqual(dense);
  });

  it('get() retrieves correct values', () => {
    const dense = [
      [1, 0, 2],
      [0, 0, 3],
      [4, 5, 6],
    ];
    const m = new DOKMatrix(dense);
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
    const m = new DOKMatrix(dense);
    expect(m.get(-1, -1)).toBe(6);  // last row, last col
    expect(m.get(-3, 0)).toBe(1);   // first row
  });

  it('set() updates values', () => {
    const m = new DOKMatrix([3, 3]);
    m.set(0, 0, 1);
    m.set(1, 1, 2);
    m.set(2, 2, 3);
    expect(m.get(0, 0)).toBe(1);
    expect(m.get(1, 1)).toBe(2);
    expect(m.get(2, 2)).toBe(3);
    expect(m.nnz).toBe(3);
  });

  it('set() with zero removes entry', () => {
    const m = new DOKMatrix([[1, 0], [0, 2]]);
    expect(m.nnz).toBe(2);
    m.set(0, 0, 0);
    expect(m.nnz).toBe(1);
    expect(m.get(0, 0)).toBe(0);
  });

  it('has() checks if entry exists', () => {
    const m = new DOKMatrix([[1, 0], [0, 2]]);
    expect(m.has(0, 0)).toBe(true);
    expect(m.has(0, 1)).toBe(false);
    expect(m.has(1, 1)).toBe(true);
  });

  it('delete() removes entry', () => {
    const m = new DOKMatrix([[1, 0], [0, 2]]);
    expect(m.delete(0, 0)).toBe(true);
    expect(m.nnz).toBe(1);
    expect(m.delete(0, 0)).toBe(false);  // already deleted
  });

  it('keys() iterates over indices', () => {
    const m = new DOKMatrix([[1, 0], [0, 2]]);
    const keys = Array.from(m.keys());
    expect(keys.length).toBe(2);
    expect(keys).toContainEqual([0, 0]);
    expect(keys).toContainEqual([1, 1]);
  });

  it('values() iterates over values', () => {
    const m = new DOKMatrix([[1, 0], [0, 2]]);
    const values = Array.from(m.values());
    expect(values.sort()).toEqual([1, 2]);
  });

  it('entries() iterates over key-value pairs', () => {
    const m = new DOKMatrix([[1, 0], [0, 2]]);
    const entries = Array.from(m.entries());
    expect(entries.length).toBe(2);
  });

  it('tocoo() conversion', () => {
    const m = new DOKMatrix([
      [1, 0, 2],
      [0, 0, 3],
      [4, 5, 6],
    ]);
    const coo = m.tocoo();
    expect(coo.format).toBe('coo');
    expect(coo.shape).toEqual([3, 3]);
    expect(coo.toArray()).toEqual([
      [1, 0, 2],
      [0, 0, 3],
      [4, 5, 6],
    ]);
  });

  it('tocsr() conversion', () => {
    const m = new DOKMatrix([
      [1, 0, 2],
      [0, 0, 3],
      [4, 5, 6],
    ]);
    const csr = m.tocsr();
    expect(csr.format).toBe('csr');
    expect(csr.toArray()).toEqual([
      [1, 0, 2],
      [0, 0, 3],
      [4, 5, 6],
    ]);
  });

  it('tocsc() conversion', () => {
    const m = new DOKMatrix([
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
    const m = new DOKMatrix([
      [1, 0, 2],
      [0, 3, 0],
      [4, 5, 6],
    ]);
    const diag = m.diagonal(0);
    expect(Array.from(diag)).toEqual([1, 3, 6]);
  });

  it('diagonal() extracts superdiagonal', () => {
    const m = new DOKMatrix([
      [1, 2, 0],
      [0, 3, 4],
      [0, 0, 5],
    ]);
    const diag = m.diagonal(1);
    expect(Array.from(diag)).toEqual([2, 4]);
  });

  it('diagonal() extracts subdiagonal', () => {
    const m = new DOKMatrix([
      [1, 0, 0],
      [2, 3, 0],
      [0, 4, 5],
    ]);
    const diag = m.diagonal(-1);
    expect(Array.from(diag)).toEqual([2, 4]);
  });

  it('copy() creates independent copy', () => {
    const m1 = new DOKMatrix([[1, 0], [0, 2]]);
    const m2 = m1.copy();
    m2.set(0, 0, 999);
    expect(m1.get(0, 0)).toBe(1);
    expect(m2.get(0, 0)).toBe(999);
  });

  it('transpose() swaps rows and cols', () => {
    const m = new DOKMatrix([
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
    const m = new DOKMatrix([5, 5]);
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
    const m = new DOKMatrix([[1, 2, 3, 4], [5, 6, 7, 8]]);
    expect(m.shape).toEqual([2, 4]);
    expect(m.toArray()).toEqual([
      [1, 2, 3, 4],
      [5, 6, 7, 8],
    ]);
  });

  it('add() delegates to CSR', () => {
    const m1 = new DOKMatrix([[1, 0], [0, 2]]);
    const m2 = new DOKMatrix([[3, 0], [0, 4]]);
    const result = m1.add(m2);
    expect(result.toArray()).toEqual([
      [4, 0],
      [0, 6],
    ]);
  });

  it('subtract() delegates to CSR', () => {
    const m1 = new DOKMatrix([[5, 0], [0, 8]]);
    const m2 = new DOKMatrix([[2, 0], [0, 3]]);
    const result = m1.subtract(m2);
    expect(result.toArray()).toEqual([
      [3, 0],
      [0, 5],
    ]);
  });

  it('multiply() element-wise delegates to CSR', () => {
    const m1 = new DOKMatrix([[2, 0], [0, 3]]);
    const m2 = new DOKMatrix([[4, 0], [0, 5]]);
    const result = m1.multiply(m2);
    expect(result.toArray()).toEqual([
      [8, 0],
      [0, 15],
    ]);
  });

  it('matmul() matrix multiply delegates to CSR', () => {
    const m1 = new DOKMatrix([[1, 2], [3, 4]]);
    const m2 = new DOKMatrix([[5, 6], [7, 8]]);
    const result = m1.matmul(m2);
    expect(result.toArray()).toEqual([
      [19, 22],
      [43, 50],
    ]);
  });

  it('dot() matrix-vector multiply delegates to CSR', () => {
    const m = new DOKMatrix([
      [1, 0, 2],
      [0, 0, 3],
      [4, 5, 6],
    ]);
    const x = new Float64Array([1, 2, 3]);
    const y = m.dot(x);
    expect(Array.from(y)).toEqual([7, 9, 32]);
  });

  it('factory function dok_matrix() works', () => {
    const m = dok_matrix([[1, 0], [0, 2]]);
    expect(m.format).toBe('dok');
    expect(m.shape).toEqual([2, 2]);
  });

  it('throws on out of bounds access', () => {
    const m = new DOKMatrix([2, 2]);
    expect(() => m.get(5, 0)).toThrow();
    expect(() => m.set(0, 5, 1)).toThrow();
  });
});
