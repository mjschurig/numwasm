/**
 * LinearOperator - Abstract interface for matrix-free linear algebra
 *
 * Based on scipy.sparse.linalg.LinearOperator
 */

import type { SparseMatrix } from '../base.js';
import type { LinearOperator as ILinearOperator } from './types.js';

/**
 * Check if a value is a LinearOperator
 */
export function isLinearOperator(x: unknown): x is ILinearOperator {
  return (
    typeof x === 'object' &&
    x !== null &&
    'shape' in x &&
    'matvec' in x &&
    typeof (x as ILinearOperator).matvec === 'function'
  );
}

/**
 * Convert input to Float64Array
 */
function asFloat64Array(x: Float64Array | { data: ArrayLike<number> }): Float64Array {
  if (x instanceof Float64Array) {
    return x;
  }
  // Object with data property
  const data = x.data;
  if (data instanceof Float64Array) {
    return data;
  }
  return new Float64Array(data);
}

/**
 * Vector operations for internal use (pure TypeScript, synchronous)
 */
function vecAdd(a: Float64Array, b: Float64Array): Float64Array {
  const result = new Float64Array(a.length);
  for (let i = 0; i < a.length; i++) {
    result[i] = a[i] + b[i];
  }
  return result;
}

function vecSub(a: Float64Array, b: Float64Array): Float64Array {
  const result = new Float64Array(a.length);
  for (let i = 0; i < a.length; i++) {
    result[i] = a[i] - b[i];
  }
  return result;
}

function vecScale(scalar: number, a: Float64Array): Float64Array {
  const result = new Float64Array(a.length);
  for (let i = 0; i < a.length; i++) {
    result[i] = scalar * a[i];
  }
  return result;
}

function vecNeg(a: Float64Array): Float64Array {
  const result = new Float64Array(a.length);
  for (let i = 0; i < a.length; i++) {
    result[i] = -a[i];
  }
  return result;
}

/**
 * Base LinearOperator class for matrix-free operations
 */
export class LinearOperator implements ILinearOperator {
  readonly shape: [number, number];
  readonly dtype: string = 'float64';

  protected _matvec: (x: Float64Array) => Float64Array;
  protected _rmatvec?: (x: Float64Array) => Float64Array;

  constructor(
    shape: [number, number],
    matvec: (x: Float64Array) => Float64Array,
    rmatvec?: (x: Float64Array) => Float64Array
  ) {
    this.shape = shape;
    this._matvec = matvec;
    this._rmatvec = rmatvec;
  }

  /**
   * Matrix-vector product: y = A @ x
   */
  matvec(x: Float64Array | { data: ArrayLike<number> }): Float64Array {
    const xarr = asFloat64Array(x);
    if (xarr.length !== this.shape[1]) {
      throw new Error(
        `Dimension mismatch: operator is ${this.shape[0]}x${this.shape[1]}, vector has length ${xarr.length}`
      );
    }
    return this._matvec(xarr);
  }

  /**
   * Transpose matrix-vector product: y = A.T @ x
   */
  rmatvec(x: Float64Array | { data: ArrayLike<number> }): Float64Array {
    const xarr = asFloat64Array(x);
    if (xarr.length !== this.shape[0]) {
      throw new Error(
        `Dimension mismatch: operator.T is ${this.shape[1]}x${this.shape[0]}, vector has length ${xarr.length}`
      );
    }
    if (this._rmatvec) {
      return this._rmatvec(xarr);
    }
    throw new Error('rmatvec is not defined for this operator');
  }

  /**
   * Matrix-matrix product: Y = A @ X
   */
  matmat(X: Float64Array, nCols: number): Float64Array {
    const [m, n] = this.shape;
    const nRows = X.length / nCols;
    if (nRows !== n) {
      throw new Error(`Dimension mismatch: operator is ${m}x${n}, matrix is ${nRows}x${nCols}`);
    }
    const result = new Float64Array(m * nCols);

    for (let j = 0; j < nCols; j++) {
      const col = new Float64Array(n);
      for (let i = 0; i < n; i++) {
        col[i] = X[i * nCols + j];
      }
      const y = this._matvec(col);
      for (let i = 0; i < m; i++) {
        result[i * nCols + j] = y[i];
      }
    }
    return result;
  }

  /**
   * Return the adjoint (conjugate transpose) operator
   * For real matrices, this is just the transpose
   */
  adjoint(): LinearOperator {
    return new _AdjointLinearOperator(this);
  }

  /**
   * Alias for adjoint
   */
  get H(): LinearOperator {
    return this.adjoint();
  }

  /**
   * Transpose operator
   */
  get T(): LinearOperator {
    return this.adjoint(); // Same as adjoint for real matrices
  }

  /**
   * Add two operators
   */
  add(other: ILinearOperator): LinearOperator {
    return new _SumLinearOperator(this, other);
  }

  /**
   * Subtract operators
   */
  sub(other: ILinearOperator): LinearOperator {
    return new _DiffLinearOperator(this, other);
  }

  /**
   * Multiply by scalar
   */
  mul(scalar: number): LinearOperator {
    return new _ScaledLinearOperator(this, scalar);
  }

  /**
   * Negate operator
   */
  neg(): LinearOperator {
    return new _NegatedLinearOperator(this);
  }

  /**
   * Compose operators: (A @ B) @ x = A @ (B @ x)
   */
  compose(other: ILinearOperator): LinearOperator {
    return new _ProductLinearOperator(this, other);
  }
}

/**
 * LinearOperator wrapping a sparse matrix
 */
class _SparseLinearOperator extends LinearOperator {
  constructor(matrix: SparseMatrix) {
    super(
      matrix.shape,
      (x: Float64Array) => {
        // Use the sparse matrix's dot method
        const csr = matrix.tocsr() as import('../compressed.js').CompressedSparseMatrix;
        return csr.dot(x);
      },
      (x: Float64Array) => {
        // Transpose matvec: A.T @ x
        // transpose() returns a CSC with swapped shape that represents A.T
        const At = matrix.transpose() as import('../compressed.js').CompressedSparseMatrix;
        return At.dot(x);
      }
    );
  }
}

/**
 * LinearOperator wrapping a dense 2D array
 */
class _DenseLinearOperator extends LinearOperator {
  constructor(data: number[][]) {
    const m = data.length;
    const n = m > 0 ? data[0].length : 0;

    super(
      [m, n],
      (x: Float64Array) => {
        const result = new Float64Array(m);
        for (let i = 0; i < m; i++) {
          let sum = 0;
          for (let j = 0; j < n; j++) {
            sum += data[i][j] * x[j];
          }
          result[i] = sum;
        }
        return result;
      },
      (x: Float64Array) => {
        const result = new Float64Array(n);
        for (let j = 0; j < n; j++) {
          let sum = 0;
          for (let i = 0; i < m; i++) {
            sum += data[i][j] * x[i];
          }
          result[j] = sum;
        }
        return result;
      }
    );
  }
}

/**
 * Adjoint (transpose) of a LinearOperator
 */
class _AdjointLinearOperator extends LinearOperator {
  constructor(base: ILinearOperator) {
    super(
      [base.shape[1], base.shape[0]],
      (x: Float64Array) => asFloat64Array(base.rmatvec(x)),
      (x: Float64Array) => asFloat64Array(base.matvec(x))
    );
  }
}

/**
 * Sum of two LinearOperators
 */
class _SumLinearOperator extends LinearOperator {
  constructor(A: ILinearOperator, B: ILinearOperator) {
    if (A.shape[0] !== B.shape[0] || A.shape[1] !== B.shape[1]) {
      throw new Error(`Shape mismatch: ${A.shape} vs ${B.shape}`);
    }
    super(
      A.shape,
      (x: Float64Array) => {
        const ya = asFloat64Array(A.matvec(x));
        const yb = asFloat64Array(B.matvec(x));
        return vecAdd(ya, yb);
      },
      (x: Float64Array) => {
        const ya = asFloat64Array(A.rmatvec(x));
        const yb = asFloat64Array(B.rmatvec(x));
        return vecAdd(ya, yb);
      }
    );
  }
}

/**
 * Difference of two LinearOperators
 */
class _DiffLinearOperator extends LinearOperator {
  constructor(A: ILinearOperator, B: ILinearOperator) {
    if (A.shape[0] !== B.shape[0] || A.shape[1] !== B.shape[1]) {
      throw new Error(`Shape mismatch: ${A.shape} vs ${B.shape}`);
    }
    super(
      A.shape,
      (x: Float64Array) => {
        const ya = asFloat64Array(A.matvec(x));
        const yb = asFloat64Array(B.matvec(x));
        return vecSub(ya, yb);
      },
      (x: Float64Array) => {
        const ya = asFloat64Array(A.rmatvec(x));
        const yb = asFloat64Array(B.rmatvec(x));
        return vecSub(ya, yb);
      }
    );
  }
}

/**
 * Scaled LinearOperator
 */
class _ScaledLinearOperator extends LinearOperator {
  constructor(base: ILinearOperator, scalar: number) {
    super(
      base.shape,
      (x: Float64Array) => {
        const y = asFloat64Array(base.matvec(x));
        return vecScale(scalar, y);
      },
      (x: Float64Array) => {
        const y = asFloat64Array(base.rmatvec(x));
        return vecScale(scalar, y);
      }
    );
  }
}

/**
 * Negated LinearOperator
 */
class _NegatedLinearOperator extends LinearOperator {
  constructor(base: ILinearOperator) {
    super(
      base.shape,
      (x: Float64Array) => vecNeg(asFloat64Array(base.matvec(x))),
      (x: Float64Array) => vecNeg(asFloat64Array(base.rmatvec(x)))
    );
  }
}

/**
 * Product of two LinearOperators: (A @ B)
 */
class _ProductLinearOperator extends LinearOperator {
  constructor(A: ILinearOperator, B: ILinearOperator) {
    if (A.shape[1] !== B.shape[0]) {
      throw new Error(`Shape mismatch for composition: ${A.shape} @ ${B.shape}`);
    }
    super(
      [A.shape[0], B.shape[1]],
      (x: Float64Array) => {
        const y = asFloat64Array(B.matvec(x));
        return asFloat64Array(A.matvec(y));
      },
      (x: Float64Array) => {
        const y = asFloat64Array(A.rmatvec(x));
        return asFloat64Array(B.rmatvec(y));
      }
    );
  }
}

/**
 * Identity LinearOperator
 */
export class IdentityOperator extends LinearOperator {
  constructor(n: number) {
    super(
      [n, n],
      (x: Float64Array) => new Float64Array(x),
      (x: Float64Array) => new Float64Array(x)
    );
  }
}

/**
 * Convert various types to a LinearOperator
 *
 * @param A - Sparse matrix, dense array, or LinearOperator
 * @returns LinearOperator wrapping the input
 */
export function aslinearoperator(A: ILinearOperator | SparseMatrix | number[][]): LinearOperator {
  // Already a LinearOperator
  if (A instanceof LinearOperator) {
    return A;
  }

  // Check if it implements the interface
  if (isLinearOperator(A)) {
    // Wrap it
    return new LinearOperator(
      A.shape,
      (x: Float64Array) => asFloat64Array(A.matvec(x)),
      A.rmatvec ? (x: Float64Array) => asFloat64Array(A.rmatvec(x)) : undefined
    );
  }

  // Check for SparseMatrix (duck typing - has shape and tocsr)
  if (
    typeof A === 'object' &&
    A !== null &&
    'shape' in A &&
    'tocsr' in A &&
    typeof (A as SparseMatrix).tocsr === 'function'
  ) {
    return new _SparseLinearOperator(A as SparseMatrix);
  }

  // Dense 2D array
  if (Array.isArray(A) && A.length > 0 && Array.isArray(A[0])) {
    return new _DenseLinearOperator(A as number[][]);
  }

  throw new Error('Cannot convert input to LinearOperator');
}
