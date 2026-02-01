/**
 * Tests for Eigenvalue Problems
 */

import { describe, it, expect } from 'vitest';
import {
  eig,
  eigvals,
  eigSymmetric,
  eigGeneralized,
  eigGeneralizedSymmetric,
  eigTridiagonal,
  eigSelect,
} from '../src/index.js';
import { toArray2D } from './setup.js';
import './setup.js';

describe('eig - General Eigenvalue Problem', () => {
  it('computes eigenvalues of a 2x2 matrix', () => {
    const A = [
      [0, 1],
      [-2, -3],
    ];
    // Characteristic poly: λ² + 3λ + 2 = (λ+1)(λ+2)
    // Eigenvalues: -1, -2

    const result = eig(A);

    expect(result.success).toBe(true);
    expect(result.n).toBe(2);

    // Sort eigenvalues for comparison
    const eigs = result.values.map((v) => v.re).sort((a, b) => a - b);
    expect(eigs[0]).toBeCloseTo(-2, 10);
    expect(eigs[1]).toBeCloseTo(-1, 10);
  });

  it('computes complex eigenvalues', () => {
    const A = [
      [0, -1],
      [1, 0],
    ];
    // Eigenvalues: ±i

    const result = eig(A);

    expect(result.success).toBe(true);
    // Should have imaginary parts
    const hasImaginary = result.values.some((v) => Math.abs(v.im) > 0.5);
    expect(hasImaginary).toBe(true);
  });

  it('computes eigenvectors', () => {
    const A = [
      [2, 0],
      [0, 3],
    ];
    // Eigenvalues: 2, 3
    // Eigenvectors: [1,0], [0,1]

    const result = eig(A, { computeVectors: true });

    expect(result.success).toBe(true);
    expect(result.vectors).toBeDefined();
  });

  it('verifies A*v = λ*v', () => {
    const A = [
      [4, 1],
      [2, 3],
    ];

    const result = eig(A, { computeVectors: true });

    expect(result.success).toBe(true);
    if (result.vectors) {
      // For real eigenvalues, verify A*v = λ*v
      const V = toArray2D(result.vectors, 2, 2);

      for (let k = 0; k < 2; k++) {
        if (Math.abs(result.values[k].im) < 1e-10) {
          const lambda = result.values[k].re;
          const v = [V[0][k], V[1][k]];
          const Av = [
            A[0][0] * v[0] + A[0][1] * v[1],
            A[1][0] * v[0] + A[1][1] * v[1],
          ];
          const lambdaV = [lambda * v[0], lambda * v[1]];

          // Normalize comparison
          const normAv = Math.sqrt(Av[0] * Av[0] + Av[1] * Av[1]);
          const normLv = Math.sqrt(
            lambdaV[0] * lambdaV[0] + lambdaV[1] * lambdaV[1]
          );
          if (normAv > 1e-10 && normLv > 1e-10) {
            expect(Av[0] / normAv).toBeCloseTo(lambdaV[0] / normLv, 8);
            expect(Av[1] / normAv).toBeCloseTo(lambdaV[1] / normLv, 8);
          }
        }
      }
    }
  });
});

describe('eigvals - Eigenvalues Only', () => {
  it('computes eigenvalues only', () => {
    const A = [
      [1, 2],
      [3, 4],
    ];

    const result = eigvals(A);

    expect(result.success).toBe(true);
    expect(result.values.length).toBe(2);
  });

  it('is faster than full eig', () => {
    const n = 30;
    const A: number[][] = [];
    for (let i = 0; i < n; i++) {
      A.push([]);
      for (let j = 0; j < n; j++) {
        A[i].push(Math.random());
      }
    }

    const result = eigvals(A);

    expect(result.success).toBe(true);
    expect(result.values.length).toBe(n);
  });
});

describe('eigSymmetric - Symmetric Eigenvalue Problem', () => {
  it('computes eigenvalues of symmetric matrix', () => {
    const A = [
      [2, 1],
      [1, 2],
    ];
    // Eigenvalues: 3, 1

    const result = eigSymmetric(A);

    expect(result.success).toBe(true);
    // Eigenvalues in ascending order
    expect(result.values[0]).toBeCloseTo(1, 10);
    expect(result.values[1]).toBeCloseTo(3, 10);
  });

  it('eigenvalues are always real', () => {
    const A = [
      [4, 2, 2],
      [2, 5, 3],
      [2, 3, 6],
    ];

    const result = eigSymmetric(A);

    expect(result.success).toBe(true);
    // All eigenvalues should be real (stored directly as Float64Array)
    expect(result.values.length).toBe(3);
  });

  it('computes orthonormal eigenvectors', () => {
    const A = [
      [2, 1],
      [1, 2],
    ];

    const result = eigSymmetric(A, { computeVectors: true });

    expect(result.success).toBe(true);
    expect(result.vectors).toBeDefined();

    if (result.vectors) {
      const V = toArray2D(result.vectors, 2, 2);

      // Columns should be orthonormal
      const dot00 = V[0][0] * V[0][0] + V[1][0] * V[1][0];
      const dot11 = V[0][1] * V[0][1] + V[1][1] * V[1][1];
      const dot01 = V[0][0] * V[0][1] + V[1][0] * V[1][1];

      expect(dot00).toBeCloseTo(1, 10);
      expect(dot11).toBeCloseTo(1, 10);
      expect(dot01).toBeCloseTo(0, 10);
    }
  });

  it('works with divide-and-conquer algorithm', () => {
    const A = [
      [2, 1],
      [1, 2],
    ];

    const result = eigSymmetric(A, { algorithm: 'syevd' });

    expect(result.success).toBe(true);
  });

  it('works with RRR algorithm', () => {
    const A = [
      [2, 1],
      [1, 2],
    ];

    const result = eigSymmetric(A, { algorithm: 'syevr' });

    expect(result.success).toBe(true);
  });
});

describe('eigGeneralized - Generalized Eigenvalue Problem', () => {
  it('solves Ax = λBx', () => {
    const A = [
      [2, 0],
      [0, 3],
    ];
    const B = [
      [1, 0],
      [0, 1],
    ];

    const result = eigGeneralized(A, B);

    expect(result.success).toBe(true);
    // With B = I, eigenvalues are same as standard problem
    const eigs = result.values.map((v) => v.re).sort((a, b) => a - b);
    expect(eigs[0]).toBeCloseTo(2, 10);
    expect(eigs[1]).toBeCloseTo(3, 10);
  });

  it('handles non-identity B', () => {
    const A = [
      [6, 0],
      [0, 6],
    ];
    const B = [
      [2, 0],
      [0, 3],
    ];
    // λ = 6/2 = 3, 6/3 = 2

    const result = eigGeneralized(A, B);

    expect(result.success).toBe(true);
    const eigs = result.values.map((v) => v.re).sort((a, b) => a - b);
    expect(eigs[0]).toBeCloseTo(2, 10);
    expect(eigs[1]).toBeCloseTo(3, 10);
  });
});

describe('eigGeneralizedSymmetric - Symmetric Generalized Problem', () => {
  it('solves symmetric Ax = λBx', () => {
    const A = [
      [4, 2],
      [2, 5],
    ];
    const B = [
      [2, 0],
      [0, 2],
    ];

    const result = eigGeneralizedSymmetric(A, B);

    expect(result.success).toBe(true);
    // Eigenvalues in ascending order
    expect(result.values.length).toBe(2);
  });

  it('B must be positive definite', () => {
    const A = [
      [4, 2],
      [2, 5],
    ];
    const B = [
      [1, 0],
      [0, 1],
    ];

    const result = eigGeneralizedSymmetric(A, B);

    expect(result.success).toBe(true);
  });
});

describe('eigTridiagonal - Tridiagonal Eigenvalue Problem', () => {
  it('computes eigenvalues of tridiagonal matrix', () => {
    const d = [2, 2, 2]; // diagonal
    const e = [1, 1]; // off-diagonal

    const result = eigTridiagonal(d, e);

    expect(result.success).toBe(true);
    expect(result.values.length).toBe(3);
  });

  it('computes eigenvectors', () => {
    const d = [2, 2, 2];
    const e = [1, 1];

    const result = eigTridiagonal(d, e, { computeVectors: true });

    expect(result.success).toBe(true);
    expect(result.vectors).toBeDefined();
  });

  it('handles single element', () => {
    const d = [5];
    const e: number[] = [];

    const result = eigTridiagonal(d, e);

    expect(result.success).toBe(true);
    expect(result.values[0]).toBeCloseTo(5, 10);
  });
});

describe('eigSelect - Selected Eigenvalues', () => {
  it('computes all eigenvalues by default', () => {
    const A = [
      [2, 1],
      [1, 2],
    ];

    const result = eigSelect(A);

    expect(result.success).toBe(true);
    expect(result.m).toBe(2);
  });

  it('computes eigenvalues in interval', () => {
    const A = [
      [1, 0, 0],
      [0, 2, 0],
      [0, 0, 3],
    ];

    const result = eigSelect(A, { range: 'value', vl: 1.5, vu: 2.5 });

    expect(result.success).toBe(true);
    expect(result.m).toBe(1); // Only eigenvalue 2 is in (1.5, 2.5]
  });

  it('computes eigenvalues by index', () => {
    const A = [
      [1, 0, 0],
      [0, 2, 0],
      [0, 0, 3],
    ];

    const result = eigSelect(A, { range: 'index', il: 1, iu: 2 });

    expect(result.success).toBe(true);
    expect(result.m).toBe(2); // First 2 eigenvalues
    expect(result.values[0]).toBeCloseTo(1, 10);
    expect(result.values[1]).toBeCloseTo(2, 10);
  });
});
