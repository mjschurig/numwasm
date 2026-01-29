/**
 * Tests for sparse construction functions: eye(), diags(), tril(), triu(),
 * hstack(), vstack(), block_diag(), kron(), kronsum(), random(), issparse()
 *
 * Ported from scipy/sparse/tests/test_construct.py
 */
import { describe, it, expect, beforeAll } from 'vitest';
import { sparse, loadWasmModule } from '../../../src/ts/index.js';

const {
  eye,
  diags,
  issparse,
  tril,
  triu,
  hstack,
  vstack,
  block_diag,
  kron,
  kronsum,
  random,
  csr_matrix,
  coo_matrix,
} = sparse;

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

describe('issparse()', () => {
  it('returns true for CSR matrix', () => {
    const m = csr_matrix([[1, 0], [0, 2]]);
    expect(issparse(m)).toBe(true);
  });

  it('returns true for COO matrix', () => {
    const m = coo_matrix([[1, 0], [0, 2]]);
    expect(issparse(m)).toBe(true);
  });

  it('returns false for dense array', () => {
    expect(issparse([[1, 0], [0, 2]])).toBe(false);
  });

  it('returns false for number', () => {
    expect(issparse(42)).toBe(false);
  });

  it('returns false for null/undefined', () => {
    expect(issparse(null)).toBe(false);
    expect(issparse(undefined)).toBe(false);
  });
});

describe('tril()', () => {
  it('tril extracts lower triangle (k=0)', () => {
    const A = csr_matrix([
      [1, 2, 3],
      [4, 5, 6],
      [7, 8, 9],
    ]);
    const L = tril(A);
    expect(L.toArray()).toEqual([
      [1, 0, 0],
      [4, 5, 0],
      [7, 8, 9],
    ]);
  });

  it('tril with k=1 includes first superdiagonal', () => {
    const A = csr_matrix([
      [1, 2, 3],
      [4, 5, 6],
      [7, 8, 9],
    ]);
    const L = tril(A, 1);
    expect(L.toArray()).toEqual([
      [1, 2, 0],
      [4, 5, 6],
      [7, 8, 9],
    ]);
  });

  it('tril with k=-1 excludes main diagonal', () => {
    const A = csr_matrix([
      [1, 2, 3],
      [4, 5, 6],
      [7, 8, 9],
    ]);
    const L = tril(A, -1);
    expect(L.toArray()).toEqual([
      [0, 0, 0],
      [4, 0, 0],
      [7, 8, 0],
    ]);
  });

  it('tril works on dense input', () => {
    const L = tril([[1, 2], [3, 4]]);
    expect(L.toArray()).toEqual([
      [1, 0],
      [3, 4],
    ]);
  });

  it('tril respects format parameter', () => {
    const L = tril([[1, 2], [3, 4]], 0, 'csr');
    expect(L.format).toBe('csr');
  });

  it('tril on empty matrix', () => {
    const A = csr_matrix([
      [0, 0],
      [0, 0],
    ]);
    const L = tril(A);
    expect(L.nnz).toBe(0);
  });
});

describe('triu()', () => {
  it('triu extracts upper triangle (k=0)', () => {
    const A = csr_matrix([
      [1, 2, 3],
      [4, 5, 6],
      [7, 8, 9],
    ]);
    const U = triu(A);
    expect(U.toArray()).toEqual([
      [1, 2, 3],
      [0, 5, 6],
      [0, 0, 9],
    ]);
  });

  it('triu with k=1 excludes main diagonal', () => {
    const A = csr_matrix([
      [1, 2, 3],
      [4, 5, 6],
      [7, 8, 9],
    ]);
    const U = triu(A, 1);
    expect(U.toArray()).toEqual([
      [0, 2, 3],
      [0, 0, 6],
      [0, 0, 0],
    ]);
  });

  it('triu with k=-1 includes first subdiagonal', () => {
    const A = csr_matrix([
      [1, 2, 3],
      [4, 5, 6],
      [7, 8, 9],
    ]);
    const U = triu(A, -1);
    expect(U.toArray()).toEqual([
      [1, 2, 3],
      [4, 5, 6],
      [0, 8, 9],
    ]);
  });

  it('triu works on dense input', () => {
    const U = triu([[1, 2], [3, 4]]);
    expect(U.toArray()).toEqual([
      [1, 2],
      [0, 4],
    ]);
  });
});

describe('hstack()', () => {
  it('hstack two matrices', () => {
    const A = csr_matrix([[1, 2], [3, 4]]);
    const B = csr_matrix([[5, 6], [7, 8]]);
    const C = hstack([A, B]);
    expect(C.shape).toEqual([2, 4]);
    expect(C.toArray()).toEqual([
      [1, 2, 5, 6],
      [3, 4, 7, 8],
    ]);
  });

  it('hstack three matrices', () => {
    const A = csr_matrix([[1], [2]]);
    const B = csr_matrix([[3], [4]]);
    const C = csr_matrix([[5], [6]]);
    const D = hstack([A, B, C]);
    expect(D.shape).toEqual([2, 3]);
    expect(D.toArray()).toEqual([
      [1, 3, 5],
      [2, 4, 6],
    ]);
  });

  it('hstack single matrix returns copy', () => {
    const A = csr_matrix([[1, 2], [3, 4]]);
    const B = hstack([A]);
    expect(B.toArray()).toEqual(A.toArray());
  });

  it('hstack throws on empty input', () => {
    expect(() => hstack([])).toThrow(/at least one/);
  });

  it('hstack throws on mismatched rows', () => {
    const A = csr_matrix([[1, 2]]);
    const B = csr_matrix([[3], [4]]);
    expect(() => hstack([A, B])).toThrow(/row/i);
  });

  it('hstack works with dense arrays', () => {
    const C = hstack([[[1], [2]], [[3], [4]]]);
    expect(C.toArray()).toEqual([
      [1, 3],
      [2, 4],
    ]);
  });
});

describe('vstack()', () => {
  it('vstack two matrices', () => {
    const A = csr_matrix([[1, 2], [3, 4]]);
    const B = csr_matrix([[5, 6], [7, 8]]);
    const C = vstack([A, B]);
    expect(C.shape).toEqual([4, 2]);
    expect(C.toArray()).toEqual([
      [1, 2],
      [3, 4],
      [5, 6],
      [7, 8],
    ]);
  });

  it('vstack three matrices', () => {
    const A = csr_matrix([[1, 2]]);
    const B = csr_matrix([[3, 4]]);
    const C = csr_matrix([[5, 6]]);
    const D = vstack([A, B, C]);
    expect(D.shape).toEqual([3, 2]);
    expect(D.toArray()).toEqual([
      [1, 2],
      [3, 4],
      [5, 6],
    ]);
  });

  it('vstack single matrix returns copy', () => {
    const A = csr_matrix([[1, 2], [3, 4]]);
    const B = vstack([A]);
    expect(B.toArray()).toEqual(A.toArray());
  });

  it('vstack throws on empty input', () => {
    expect(() => vstack([])).toThrow(/at least one/);
  });

  it('vstack throws on mismatched columns', () => {
    const A = csr_matrix([[1, 2]]);
    const B = csr_matrix([[3, 4, 5]]);
    expect(() => vstack([A, B])).toThrow(/column/i);
  });
});

describe('block_diag()', () => {
  it('block_diag two matrices', () => {
    const A = csr_matrix([[1, 2], [3, 4]]);
    const B = csr_matrix([[5, 6], [7, 8]]);
    const C = block_diag([A, B]);
    expect(C.shape).toEqual([4, 4]);
    expect(C.toArray()).toEqual([
      [1, 2, 0, 0],
      [3, 4, 0, 0],
      [0, 0, 5, 6],
      [0, 0, 7, 8],
    ]);
  });

  it('block_diag three matrices of different sizes', () => {
    const A = csr_matrix([[1]]);
    const B = csr_matrix([[2, 3], [4, 5]]);
    const C = csr_matrix([[6]]);
    const D = block_diag([A, B, C]);
    expect(D.shape).toEqual([4, 4]);
    expect(D.toArray()).toEqual([
      [1, 0, 0, 0],
      [0, 2, 3, 0],
      [0, 4, 5, 0],
      [0, 0, 0, 6],
    ]);
  });

  it('block_diag single matrix', () => {
    const A = csr_matrix([[1, 2], [3, 4]]);
    const B = block_diag([A]);
    expect(B.toArray()).toEqual(A.toArray());
  });

  it('block_diag empty list', () => {
    const B = block_diag([]);
    expect(B.shape).toEqual([0, 0]);
    expect(B.nnz).toBe(0);
  });

  it('block_diag respects format', () => {
    const A = csr_matrix([[1]]);
    const B = block_diag([A], 'csc');
    expect(B.format).toBe('csc');
  });
});

describe('kron()', () => {
  it('kron identity x identity', () => {
    const I2 = eye(2);
    const K = kron(I2, I2);
    expect(K.shape).toEqual([4, 4]);
    expect(K.toArray()).toEqual([
      [1, 0, 0, 0],
      [0, 1, 0, 0],
      [0, 0, 1, 0],
      [0, 0, 0, 1],
    ]);
  });

  it('kron with scalar matrices', () => {
    const A = csr_matrix([[2]]);
    const B = csr_matrix([[3]]);
    const K = kron(A, B);
    expect(K.toArray()).toEqual([[6]]);
  });

  it('kron 2x2 x 2x2', () => {
    const A = csr_matrix([[1, 2], [3, 4]]);
    const B = csr_matrix([[0, 5], [6, 7]]);
    const K = kron(A, B);
    expect(K.shape).toEqual([4, 4]);
    expect(K.toArray()).toEqual([
      [0, 5, 0, 10],
      [6, 7, 12, 14],
      [0, 15, 0, 20],
      [18, 21, 24, 28],
    ]);
  });

  it('kron with empty matrix', () => {
    const A = csr_matrix([[1, 2]]);
    const B = csr_matrix({ data: [], indices: [], indptr: [0] }, { shape: [1, 2] });
    const K = kron(A, B);
    expect(K.shape).toEqual([1, 4]);
    expect(K.nnz).toBe(0);
  });

  it('kron respects format', () => {
    const A = eye(2);
    const B = eye(2);
    const K = kron(A, B, 'coo');
    expect(K.format).toBe('coo');
  });
});

describe('kronsum()', () => {
  it('kronsum of 2x2 matrices', () => {
    // kronsum(A, B) = kron(A, I_b) + kron(I_a, B)
    const A = csr_matrix([[1, 0], [0, 2]]);
    const B = csr_matrix([[3, 0], [0, 4]]);
    const K = kronsum(A, B);
    expect(K.shape).toEqual([4, 4]);
    // kron(A, I2) + kron(I2, B)
    // = [[1,0,0,0],[0,1,0,0],[0,0,2,0],[0,0,0,2]] + [[3,0,0,0],[0,4,0,0],[0,0,3,0],[0,0,0,4]]
    // = [[4,0,0,0],[0,5,0,0],[0,0,5,0],[0,0,0,6]]
    expect(K.toArray()).toEqual([
      [4, 0, 0, 0],
      [0, 5, 0, 0],
      [0, 0, 5, 0],
      [0, 0, 0, 6],
    ]);
  });

  it('kronsum throws on non-square A', () => {
    const A = csr_matrix([[1, 2, 3]]);
    const B = csr_matrix([[1]]);
    expect(() => kronsum(A, B)).toThrow(/square/);
  });

  it('kronsum throws on non-square B', () => {
    const A = csr_matrix([[1]]);
    const B = csr_matrix([[1, 2], [3, 4], [5, 6]]);
    expect(() => kronsum(A, B)).toThrow(/square/);
  });
});

describe('random()', () => {
  it('random generates correct shape', () => {
    const R = random(10, 20, 0.1);
    expect(R.shape).toEqual([10, 20]);
  });

  it('random with density=0 gives empty matrix', () => {
    const R = random(10, 10, 0);
    expect(R.nnz).toBe(0);
  });

  it('random with density=1 gives full matrix', () => {
    const R = random(3, 3, 1);
    expect(R.nnz).toBe(9);
  });

  it('random respects approximate density', () => {
    const R = random(100, 100, 0.1);
    // Should have approximately 1000 non-zeros (10% of 10000)
    expect(R.nnz).toBeGreaterThan(800);
    expect(R.nnz).toBeLessThan(1200);
  });

  it('random with custom rng is deterministic', () => {
    let seed = 12345;
    const rng = () => {
      seed = (seed * 1103515245 + 12345) % 2147483648;
      return seed / 2147483648;
    };
    const R1 = random(5, 5, 0.5, 'csr', rng);

    seed = 12345;
    const R2 = random(5, 5, 0.5, 'csr', rng);

    expect(R1.toArray()).toEqual(R2.toArray());
  });

  it('random throws on invalid density', () => {
    expect(() => random(10, 10, -0.1)).toThrow(/density/);
    expect(() => random(10, 10, 1.5)).toThrow(/density/);
  });

  it('random respects format', () => {
    const R = random(5, 5, 0.5, 'coo');
    expect(R.format).toBe('coo');
  });

  it('random values are in [0, 1)', () => {
    const R = random(10, 10, 0.5);
    const arr = R.toArray();
    for (const row of arr) {
      for (const val of row) {
        if (val !== 0) {
          expect(val).toBeGreaterThanOrEqual(0);
          expect(val).toBeLessThan(1);
        }
      }
    }
  });
});
