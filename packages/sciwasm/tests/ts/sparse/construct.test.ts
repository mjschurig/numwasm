/**
 * Tests for sparse construction functions: eye(), diags()
 *
 * Ported from scipy/sparse/tests/test_construct.py
 */
import { describe, it, expect, beforeAll } from 'vitest';
import { sparse, loadWasmModule } from '../../../src/ts/index.js';

const { eye, diags } = sparse;

beforeAll(async () => {
  await loadWasmModule();
});

describe('eye()', () => {
  it('eye(1) gives [[1]]', () => {
    const m = eye(1);
    expect(m.shape).toEqual([1, 1]);
    expect(m.toArray()).toEqual([[1]]);
  });

  it('eye(2) gives 2x2 identity', () => {
    const m = eye(2);
    expect(m.toArray()).toEqual([
      [1, 0],
      [0, 1],
    ]);
  });

  it('eye(3) gives 3x3 identity', () => {
    const m = eye(3);
    expect(m.shape).toEqual([3, 3]);
    expect(m.nnz).toBe(3);
    expect(m.toArray()).toEqual([
      [1, 0, 0],
      [0, 1, 0],
      [0, 0, 1],
    ]);
  });

  it('eye(2, 3) gives rectangular identity', () => {
    const m = eye(2, 3);
    expect(m.shape).toEqual([2, 3]);
    expect(m.toArray()).toEqual([
      [1, 0, 0],
      [0, 1, 0],
    ]);
  });

  it('eye(3, 2) gives tall rectangular identity', () => {
    const m = eye(3, 2);
    expect(m.shape).toEqual([3, 2]);
    expect(m.toArray()).toEqual([
      [1, 0],
      [0, 1],
      [0, 0],
    ]);
  });

  it('eye(3, 3, 1) gives superdiagonal', () => {
    const m = eye(3, 3, 1);
    expect(m.toArray()).toEqual([
      [0, 1, 0],
      [0, 0, 1],
      [0, 0, 0],
    ]);
  });

  it('eye(3, 3, -1) gives subdiagonal', () => {
    const m = eye(3, 3, -1);
    expect(m.toArray()).toEqual([
      [0, 0, 0],
      [1, 0, 0],
      [0, 1, 0],
    ]);
  });

  it('eye(3, 3, 2) gives second superdiagonal', () => {
    const m = eye(3, 3, 2);
    expect(m.toArray()).toEqual([
      [0, 0, 1],
      [0, 0, 0],
      [0, 0, 0],
    ]);
  });

  it('eye(3, 3, -2) gives second subdiagonal', () => {
    const m = eye(3, 3, -2);
    expect(m.toArray()).toEqual([
      [0, 0, 0],
      [0, 0, 0],
      [1, 0, 0],
    ]);
  });

  it('eye(4, 3, 0) rectangular main diagonal', () => {
    const m = eye(4, 3);
    expect(m.shape).toEqual([4, 3]);
    expect(m.toArray()).toEqual([
      [1, 0, 0],
      [0, 1, 0],
      [0, 0, 1],
      [0, 0, 0],
    ]);
  });

  it('eye with format=csc returns CSC', () => {
    const m = eye(3, undefined, 0, 'csc');
    expect(m.format).toBe('csc');
    expect(m.toArray()).toEqual([
      [1, 0, 0],
      [0, 1, 0],
      [0, 0, 1],
    ]);
  });

  it('eye(3, 3, 3) gives empty (offset beyond matrix)', () => {
    const m = eye(3, 3, 3);
    expect(m.nnz).toBe(0);
    expect(m.toArray()).toEqual([
      [0, 0, 0],
      [0, 0, 0],
      [0, 0, 0],
    ]);
  });
});

describe('diags()', () => {
  it('diags with single scalar diagonal', () => {
    const m = diags([1], [0], [3, 3]);
    expect(m.toArray()).toEqual([
      [1, 0, 0],
      [0, 1, 0],
      [0, 0, 1],
    ]);
  });

  it('diags with array diagonal', () => {
    const m = diags([[1, 2, 3]], [0], [3, 3]);
    expect(m.toArray()).toEqual([
      [1, 0, 0],
      [0, 2, 0],
      [0, 0, 3],
    ]);
  });

  it('diags with super-diagonal', () => {
    const m = diags([[1, 2]], [1], [3, 3]);
    expect(m.toArray()).toEqual([
      [0, 1, 0],
      [0, 0, 2],
      [0, 0, 0],
    ]);
  });

  it('diags with sub-diagonal', () => {
    const m = diags([[4, 5]], [-1], [3, 3]);
    expect(m.toArray()).toEqual([
      [0, 0, 0],
      [4, 0, 0],
      [0, 5, 0],
    ]);
  });

  it('diags tridiagonal with scalars', () => {
    const m = diags([1, -2, 1], [-1, 0, 1], [4, 4]);
    expect(m.toArray()).toEqual([
      [-2, 1, 0, 0],
      [1, -2, 1, 0],
      [0, 1, -2, 1],
      [0, 0, 1, -2],
    ]);
  });

  it('diags tridiagonal with arrays', () => {
    const m = diags(
      [
        [1, 1, 1],
        [-2, -2, -2, -2],
        [1, 1, 1],
      ],
      [-1, 0, 1],
      [4, 4]
    );
    expect(m.toArray()).toEqual([
      [-2, 1, 0, 0],
      [1, -2, 1, 0],
      [0, 1, -2, 1],
      [0, 0, 1, -2],
    ]);
  });

  it('diags rectangular', () => {
    const m = diags([[1, 2, 3]], [0], [3, 4]);
    expect(m.shape).toEqual([3, 4]);
    expect(m.toArray()).toEqual([
      [1, 0, 0, 0],
      [0, 2, 0, 0],
      [0, 0, 3, 0],
    ]);
  });

  it('diags with offset=2', () => {
    const m = diags([[6]], [2], [3, 3]);
    expect(m.toArray()).toEqual([
      [0, 0, 6],
      [0, 0, 0],
      [0, 0, 0],
    ]);
  });

  it('diags single offset (number, not array)', () => {
    const m = diags([[1, 2, 3]], 0, [3, 3]);
    expect(m.toArray()).toEqual([
      [1, 0, 0],
      [0, 2, 0],
      [0, 0, 3],
    ]);
  });

  it('diags with format=csc returns CSC', () => {
    const m = diags([1], [0], [3, 3], 'csc');
    expect(m.format).toBe('csc');
    expect(m.toArray()).toEqual([
      [1, 0, 0],
      [0, 1, 0],
      [0, 0, 1],
    ]);
  });

  it('diags infers shape from first diagonal', () => {
    // shape not provided: size = len(diag[0]) + abs(offset[0])
    const m = diags([[1, 2, 3]]);
    expect(m.shape).toEqual([3, 3]);
    expect(m.toArray()).toEqual([
      [1, 0, 0],
      [0, 2, 0],
      [0, 0, 3],
    ]);
  });

  it('diags throws on mismatched lengths', () => {
    expect(() => diags([[1, 2]], [0, 1], [3, 3])).toThrow(
      /does not match/
    );
  });

  it('diags zero values are not stored', () => {
    const m = diags([[0, 0, 0]], [0], [3, 3]);
    expect(m.nnz).toBe(0);
    expect(m.toArray()).toEqual([
      [0, 0, 0],
      [0, 0, 0],
      [0, 0, 0],
    ]);
  });
});
