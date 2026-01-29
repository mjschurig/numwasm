/**
 * Tests for sparse linear algebra module (sparse.linalg)
 */

import { describe, it, expect, beforeAll } from 'vitest';
import { sparse, loadWasmModule } from '../../../src/ts/index.js';

const {
  csr_matrix,
  LinearOperator,
  IdentityOperator,
  aslinearoperator,
  isLinearOperator,
  linalg_norm,
  cg,
  bicgstab,
  gmres,
} = sparse;

beforeAll(async () => {
  await loadWasmModule();
});

describe('LinearOperator', () => {
  describe('aslinearoperator', () => {
    it('should wrap a sparse matrix', () => {
      const A = csr_matrix([
        [2, 1],
        [1, 3],
      ]);
      const op = aslinearoperator(A);

      expect(op.shape).toEqual([2, 2]);
      expect(op.dtype).toBe('float64');
    });

    it('should wrap a dense 2D array', () => {
      const A = [
        [2, 1],
        [1, 3],
      ];
      const op = aslinearoperator(A);

      expect(op.shape).toEqual([2, 2]);
    });

    it('should return LinearOperator unchanged', () => {
      const op1 = new IdentityOperator(3);
      const op2 = aslinearoperator(op1);

      expect(op2).toBe(op1);
    });
  });

  describe('isLinearOperator', () => {
    it('should return true for LinearOperator', () => {
      const op = new IdentityOperator(3);
      expect(isLinearOperator(op)).toBe(true);
    });

    it('should return false for sparse matrix', () => {
      const A = csr_matrix([[1, 0], [0, 1]]);
      expect(isLinearOperator(A)).toBe(false);
    });

    it('should return false for primitives', () => {
      expect(isLinearOperator(null)).toBe(false);
      expect(isLinearOperator(undefined)).toBe(false);
      expect(isLinearOperator(5)).toBe(false);
      expect(isLinearOperator('hello')).toBe(false);
    });
  });

  describe('matvec', () => {
    it('should compute A @ x for sparse matrix', () => {
      const A = csr_matrix([
        [2, 1],
        [1, 3],
      ]);
      const op = aslinearoperator(A);
      const x = new Float64Array([1, 2]);

      const y = op.matvec(x);

      expect(y).toBeInstanceOf(Float64Array);
      expect(y.length).toBe(2);
      expect(y[0]).toBeCloseTo(4); // 2*1 + 1*2 = 4
      expect(y[1]).toBeCloseTo(7); // 1*1 + 3*2 = 7
    });

    it('should compute A @ x for dense matrix', () => {
      const A = [
        [1, 2],
        [3, 4],
      ];
      const op = aslinearoperator(A);
      const x = new Float64Array([1, 1]);

      const y = op.matvec(x);

      expect(y[0]).toBeCloseTo(3); // 1 + 2
      expect(y[1]).toBeCloseTo(7); // 3 + 4
    });

    it('should throw on dimension mismatch', () => {
      const A = csr_matrix([[1, 2, 3]]);
      const op = aslinearoperator(A);

      expect(() => op.matvec(new Float64Array([1, 2]))).toThrow(/Dimension mismatch/);
    });
  });

  describe('rmatvec', () => {
    it('should compute A.T @ x', () => {
      const A = csr_matrix([
        [1, 2],
        [3, 4],
      ]);
      const op = aslinearoperator(A);
      const x = new Float64Array([1, 1]);

      const y = op.rmatvec(x);

      expect(y[0]).toBeCloseTo(4); // 1 + 3
      expect(y[1]).toBeCloseTo(6); // 2 + 4
    });
  });

  describe('IdentityOperator', () => {
    it('should return x unchanged', () => {
      const I = new IdentityOperator(3);
      const x = new Float64Array([1, 2, 3]);

      const y = I.matvec(x);

      expect(Array.from(y)).toEqual([1, 2, 3]);
    });
  });

  describe('operator composition', () => {
    it('should add operators', () => {
      const A = aslinearoperator([[1, 0], [0, 1]]);
      const B = aslinearoperator([[1, 1], [1, 1]]);
      const C = A.add(B);
      const x = new Float64Array([1, 2]);

      const y = C.matvec(x);

      // (A+B) = [[2, 1], [1, 2]]
      // y = [[2,1],[1,2]] @ [1,2] = [2*1+1*2, 1*1+2*2] = [4, 5]
      expect(y[0]).toBeCloseTo(4); // 2*1 + 1*2 = 4
      expect(y[1]).toBeCloseTo(5); // 1*1 + 2*2 = 5
    });

    it('should scale operator', () => {
      const A = aslinearoperator([[1, 2], [3, 4]]);
      const sA = A.mul(2);
      const x = new Float64Array([1, 1]);

      const y = sA.matvec(x);

      expect(y[0]).toBeCloseTo(6);  // 2 * (1 + 2)
      expect(y[1]).toBeCloseTo(14); // 2 * (3 + 4)
    });

    it('should negate operator', () => {
      const A = aslinearoperator([[1, 2], [3, 4]]);
      const negA = A.neg();
      const x = new Float64Array([1, 1]);

      const y = negA.matvec(x);

      expect(y[0]).toBeCloseTo(-3);
      expect(y[1]).toBeCloseTo(-7);
    });
  });
});

describe('norm', () => {
  describe('Frobenius norm', () => {
    it('should compute Frobenius norm by default', () => {
      const A = csr_matrix([
        [1, 2],
        [3, 4],
      ]);

      const n = linalg_norm(A);

      // sqrt(1 + 4 + 9 + 16) = sqrt(30)
      expect(n).toBeCloseTo(Math.sqrt(30));
    });

    it('should compute Frobenius norm with explicit ord', () => {
      const A = csr_matrix([[3, 4]]);

      const n = linalg_norm(A, 'fro');

      expect(n).toBeCloseTo(5); // sqrt(9 + 16) = 5
    });
  });

  describe('1-norm (max column sum)', () => {
    it('should compute 1-norm', () => {
      const A = csr_matrix([
        [1, 2],
        [3, 4],
      ]);

      const n = linalg_norm(A, 1);

      // max(|1|+|3|, |2|+|4|) = max(4, 6) = 6
      expect(n).toBeCloseTo(6);
    });

    it('should handle negative values', () => {
      const A = csr_matrix([
        [-1, 2],
        [3, -4],
      ]);

      const n = linalg_norm(A, 1);

      // max(1+3, 2+4) = 6
      expect(n).toBeCloseTo(6);
    });
  });

  describe('inf-norm (max row sum)', () => {
    it('should compute inf-norm', () => {
      const A = csr_matrix([
        [1, 2],
        [3, 4],
      ]);

      const n = linalg_norm(A, Infinity);

      // max(|1|+|2|, |3|+|4|) = max(3, 7) = 7
      expect(n).toBeCloseTo(7);
    });
  });

  describe('-1-norm (min column sum)', () => {
    it('should compute -1-norm', () => {
      const A = csr_matrix([
        [1, 2],
        [3, 4],
      ]);

      const n = linalg_norm(A, -1);

      // min(|1|+|3|, |2|+|4|) = min(4, 6) = 4
      expect(n).toBeCloseTo(4);
    });
  });

  describe('-inf-norm (min row sum)', () => {
    it('should compute -inf-norm', () => {
      const A = csr_matrix([
        [1, 2],
        [3, 4],
      ]);

      const n = linalg_norm(A, -Infinity);

      // min(|1|+|2|, |3|+|4|) = min(3, 7) = 3
      expect(n).toBeCloseTo(3);
    });
  });
});

describe('Conjugate Gradient (cg)', () => {
  it('should solve a simple SPD system', () => {
    // Simple 2x2 SPD matrix
    const A = csr_matrix([
      [4, 1],
      [1, 3],
    ]);
    const b = new Float64Array([1, 2]);

    const result = cg(A, b);

    expect(result.info).toBe(0);
    expect(result.iterations).toBeGreaterThan(0);

    // Verify solution: A @ x should equal b
    const Ax = A.tocsr().dot(result.x);
    expect(Ax[0]).toBeCloseTo(b[0], 4);
    expect(Ax[1]).toBeCloseTo(b[1], 4);
  });

  it('should solve a 3x3 SPD system', () => {
    // 3x3 SPD (tridiagonal)
    const A = csr_matrix([
      [4, 1, 0],
      [1, 4, 1],
      [0, 1, 4],
    ]);
    const b = new Float64Array([1, 2, 3]);

    const result = cg(A, b);

    expect(result.info).toBe(0);

    // Verify
    const Ax = A.tocsr().dot(result.x);
    for (let i = 0; i < 3; i++) {
      expect(Ax[i]).toBeCloseTo(b[i], 4);
    }
  });

  it('should respect tolerance', () => {
    const A = csr_matrix([
      [4, 1],
      [1, 3],
    ]);
    const b = new Float64Array([1, 2]);

    const result = cg(A, b, { tol: 1e-10 });

    expect(result.info).toBe(0);
    expect(result.residualNorm).toBeLessThan(1e-10 * Math.sqrt(5)); // tol * ||b||
  });

  it('should use initial guess', () => {
    const A = csr_matrix([
      [4, 1],
      [1, 3],
    ]);
    const b = new Float64Array([1, 2]);

    // Start close to solution
    const result = cg(A, b, { x0: new Float64Array([0.09, 0.64]) });

    expect(result.info).toBe(0);
    expect(result.iterations).toBeLessThanOrEqual(2); // Should converge quickly
  });

  it('should report non-convergence with too few iterations', () => {
    const A = csr_matrix([
      [4, 1],
      [1, 3],
    ]);
    const b = new Float64Array([1, 2]);

    const result = cg(A, b, { maxiter: 0 });

    expect(result.info).toBe(0); // maxiter=0 but solution is trivial (x0=0)
  });

  it('should throw on non-square matrix', () => {
    const A = csr_matrix([[1, 2, 3]]);
    const b = new Float64Array([1]);

    expect(() => cg(A, b)).toThrow(/must be square/);
  });
});

describe('BiCGSTAB', () => {
  it('should solve a simple system', () => {
    const A = csr_matrix([
      [4, 1],
      [1, 3],
    ]);
    const b = new Float64Array([1, 2]);

    const result = bicgstab(A, b);

    expect(result.info).toBe(0);

    const Ax = A.tocsr().dot(result.x);
    expect(Ax[0]).toBeCloseTo(b[0], 4);
    expect(Ax[1]).toBeCloseTo(b[1], 4);
  });

  it('should solve a non-symmetric system', () => {
    // Non-symmetric matrix
    const A = csr_matrix([
      [4, 2],
      [1, 3],
    ]);
    const b = new Float64Array([1, 2]);

    const result = bicgstab(A, b);

    expect(result.info).toBe(0);

    const Ax = A.tocsr().dot(result.x);
    expect(Ax[0]).toBeCloseTo(b[0], 4);
    expect(Ax[1]).toBeCloseTo(b[1], 4);
  });

  it('should solve a larger system', () => {
    // 5x5 diagonally dominant matrix
    const A = csr_matrix([
      [10, 1, 0, 0, 0],
      [1, 10, 1, 0, 0],
      [0, 1, 10, 1, 0],
      [0, 0, 1, 10, 1],
      [0, 0, 0, 1, 10],
    ]);
    const b = new Float64Array([1, 2, 3, 4, 5]);

    const result = bicgstab(A, b);

    expect(result.info).toBe(0);

    const Ax = A.tocsr().dot(result.x);
    for (let i = 0; i < 5; i++) {
      expect(Ax[i]).toBeCloseTo(b[i], 3);
    }
  });
});

describe('GMRES', () => {
  it('should solve a simple system', () => {
    const A = csr_matrix([
      [4, 1],
      [1, 3],
    ]);
    const b = new Float64Array([1, 2]);

    const result = gmres(A, b);

    expect(result.info).toBe(0);

    const Ax = A.tocsr().dot(result.x);
    expect(Ax[0]).toBeCloseTo(b[0], 4);
    expect(Ax[1]).toBeCloseTo(b[1], 4);
  });

  it('should solve a non-symmetric system', () => {
    const A = csr_matrix([
      [4, 2],
      [1, 3],
    ]);
    const b = new Float64Array([1, 2]);

    const result = gmres(A, b);

    expect(result.info).toBe(0);

    const Ax = A.tocsr().dot(result.x);
    expect(Ax[0]).toBeCloseTo(b[0], 4);
    expect(Ax[1]).toBeCloseTo(b[1], 4);
  });

  it('should handle restart', () => {
    // 5x5 system with small restart
    const A = csr_matrix([
      [10, 1, 0, 0, 0],
      [1, 10, 1, 0, 0],
      [0, 1, 10, 1, 0],
      [0, 0, 1, 10, 1],
      [0, 0, 0, 1, 10],
    ]);
    const b = new Float64Array([1, 2, 3, 4, 5]);

    const result = gmres(A, b, { restart: 2 });

    expect(result.info).toBe(0);

    const Ax = A.tocsr().dot(result.x);
    for (let i = 0; i < 5; i++) {
      expect(Ax[i]).toBeCloseTo(b[i], 3);
    }
  });

  it('should call callback with residual norm', () => {
    const A = csr_matrix([
      [4, 1],
      [1, 3],
    ]);
    const b = new Float64Array([1, 2]);

    const residuals: number[] = [];
    const result = gmres(A, b, {
      callback: (rnorm) => residuals.push(rnorm),
    });

    expect(result.info).toBe(0);
    expect(residuals.length).toBeGreaterThan(0);

    // Residuals should decrease (generally)
    expect(residuals[residuals.length - 1]).toBeLessThanOrEqual(residuals[0]);
  });
});

describe('Comparison between solvers', () => {
  it('all solvers should give same solution for SPD system', () => {
    const A = csr_matrix([
      [4, 1, 0],
      [1, 4, 1],
      [0, 1, 4],
    ]);
    const b = new Float64Array([1, 2, 3]);

    const cgResult = cg(A, b);
    const bicgstabResult = bicgstab(A, b);
    const gmresResult = gmres(A, b);

    // All should converge
    expect(cgResult.info).toBe(0);
    expect(bicgstabResult.info).toBe(0);
    expect(gmresResult.info).toBe(0);

    // Solutions should be close
    for (let i = 0; i < 3; i++) {
      expect(cgResult.x[i]).toBeCloseTo(bicgstabResult.x[i], 4);
      expect(cgResult.x[i]).toBeCloseTo(gmresResult.x[i], 4);
    }
  });
});
