/**
 * Symbolic matrix operations.
 * @module matrices
 */

import type { Expr } from '../core/index.js';
import { NotImplementedError } from '../errors.js';

/**
 * Symbolic matrix.
 * Mirrors sympy.Matrix.
 */
export class Matrix {
  readonly rows: number;
  readonly cols: number;

  constructor(_data: Expr[][] | number[][]) {
    this.rows = 0;
    this.cols = 0;
    throw new NotImplementedError('symwasm.matrices.Matrix');
  }

  /** Get element at (i, j). */
  get(_i: number, _j: number): Expr {
    throw new NotImplementedError('symwasm.matrices.Matrix.get');
  }

  /** Compute the determinant. */
  det(): Expr {
    throw new NotImplementedError('symwasm.matrices.Matrix.det');
  }

  /** Compute the inverse. */
  inv(): Matrix {
    throw new NotImplementedError('symwasm.matrices.Matrix.inv');
  }

  /** Compute the transpose. */
  transpose(): Matrix {
    throw new NotImplementedError('symwasm.matrices.Matrix.transpose');
  }

  /** Compute eigenvalues. */
  eigenvals(): Map<Expr, number> {
    throw new NotImplementedError('symwasm.matrices.Matrix.eigenvals');
  }

  /** Compute eigenvectors. */
  eigenvects(): Array<{ eigenval: Expr; multiplicity: number; vectors: Matrix[] }> {
    throw new NotImplementedError('symwasm.matrices.Matrix.eigenvects');
  }

  /** Row-reduce to echelon form. */
  rref(): { matrix: Matrix; pivots: number[] } {
    throw new NotImplementedError('symwasm.matrices.Matrix.rref');
  }
}

/** Create an identity matrix. */
export function eye(_n: number): Matrix {
  throw new NotImplementedError('symwasm.matrices.eye');
}

/** Create a zero matrix. */
export function zeros(_rows: number, _cols: number): Matrix {
  throw new NotImplementedError('symwasm.matrices.zeros');
}

/** Create a matrix of ones. */
export function ones(_rows: number, _cols: number): Matrix {
  throw new NotImplementedError('symwasm.matrices.ones');
}

/** Create a diagonal matrix. */
export function diag(..._values: (Expr | number)[]): Matrix {
  throw new NotImplementedError('symwasm.matrices.diag');
}
