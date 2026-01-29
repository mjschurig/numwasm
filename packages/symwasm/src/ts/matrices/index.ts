/**
 * Symbolic matrix operations.
 * @module matrices
 */

import { Expr, Integer, exprFromWasm } from '../core/index.js';
import { getWasmModule } from '../wasm-loader.js';
import {
  DenseMatrixObject,
  createDenseMatrix,
  SymEngineVec,
  createBasic,
  checkException,
} from '../wasm-memory.js';
import { NotImplementedError } from '../errors.js';

/**
 * Helper to convert a value to an Expr
 * @internal
 */
function toExpr(value: Expr | number): Expr {
  if (typeof value === 'number') {
    return new Integer(value);
  }
  return value;
}

/**
 * Symbolic dense matrix.
 * Mirrors sympy.Matrix.
 *
 * Matrix is backed by SymEngine's CDenseMatrix and contains symbolic expressions.
 * Unlike Expr, Matrix does not extend the expression hierarchy since matrices
 * are containers for expressions, not expressions themselves.
 */
export class Matrix {
  /** @internal WASM pointer wrapper */
  private _obj: DenseMatrixObject;

  /**
   * Create a matrix from a 2D array of expressions or numbers.
   * @param data Nested array where each inner array is a row
   *
   * @example
   * const m = new Matrix([[1, 2], [3, 4]]);
   * const symbolic = new Matrix([[x, y], [z, w]]);
   */
  constructor(data: (Expr | number)[][]) {
    if (!data || data.length === 0) {
      throw new Error('Matrix data cannot be empty');
    }

    const rows = data.length;
    const cols = data[0].length;

    if (cols === 0) {
      throw new Error('Matrix rows cannot be empty');
    }

    // Validate rectangular
    for (const row of data) {
      if (row.length !== cols) {
        throw new Error('All rows must have the same number of columns');
      }
    }

    const wasm = getWasmModule();

    // Build a CVecBasic with all elements in row-major order
    const vec = new SymEngineVec();
    try {
      for (let i = 0; i < rows; i++) {
        for (let j = 0; j < cols; j++) {
          const expr = toExpr(data[i][j]);
          wasm._vecbasic_push_back(vec.getPtr(), expr.getWasmPtr());
        }
      }

      // Create matrix from vector
      const ptr = wasm._dense_matrix_new_vec(rows, cols, vec.getPtr());
      this._obj = new DenseMatrixObject(ptr);
    } finally {
      vec.free();
    }
  }

  /**
   * @internal Create a Matrix from an existing DenseMatrixObject
   */
  static _fromDenseMatrixObject(obj: DenseMatrixObject): Matrix {
    const m = Object.create(Matrix.prototype) as Matrix;
    m._obj = obj;
    return m;
  }

  /**
   * Create a matrix from a flat array.
   * @param flat Array of elements in row-major order
   * @param rows Number of rows
   * @param cols Number of columns
   *
   * @example
   * const m = Matrix.fromFlat([1, 2, 3, 4, 5, 6], 2, 3);
   * // Creates [[1, 2, 3], [4, 5, 6]]
   */
  static fromFlat(flat: (Expr | number)[], rows: number, cols: number): Matrix {
    if (flat.length !== rows * cols) {
      throw new Error(
        `Flat array length (${flat.length}) does not match rows*cols (${rows * cols})`
      );
    }

    const wasm = getWasmModule();
    const vec = new SymEngineVec();

    try {
      for (const elem of flat) {
        const expr = toExpr(elem);
        wasm._vecbasic_push_back(vec.getPtr(), expr.getWasmPtr());
      }

      const ptr = wasm._dense_matrix_new_vec(rows, cols, vec.getPtr());
      return Matrix._fromDenseMatrixObject(new DenseMatrixObject(ptr));
    } finally {
      vec.free();
    }
  }

  /**
   * @internal Get the underlying WASM pointer
   */
  getWasmPtr(): number {
    return this._obj.getPtr();
  }

  /**
   * Number of rows in the matrix.
   */
  get rows(): number {
    return this._obj.rows();
  }

  /**
   * Number of columns in the matrix.
   */
  get cols(): number {
    return this._obj.cols();
  }

  /**
   * Shape of the matrix as [rows, cols] tuple.
   */
  get shape(): [number, number] {
    return [this.rows, this.cols];
  }

  /**
   * Get element at position (i, j).
   * @param i Row index (0-based)
   * @param j Column index (0-based)
   * @returns The expression at the specified position
   *
   * @example
   * const m = new Matrix([[1, 2], [3, 4]]);
   * m.get(0, 1);  // Returns Integer(2)
   */
  get(i: number, j: number): Expr {
    if (i < 0 || i >= this.rows || j < 0 || j >= this.cols) {
      throw new RangeError(
        `Index (${i}, ${j}) out of bounds for matrix of shape (${this.rows}, ${this.cols})`
      );
    }

    const wasm = getWasmModule();
    const result = createBasic();

    try {
      const code = wasm._dense_matrix_get_basic(
        result.getPtr(),
        this._obj.getPtr(),
        i,
        j
      );
      checkException(code);
      return exprFromWasm(result);
    } catch (e) {
      result.free();
      throw e;
    }
  }

  /**
   * Set element at position (i, j).
   * @param i Row index (0-based)
   * @param j Column index (0-based)
   * @param value The expression or number to set
   *
   * @example
   * const m = zeros(2, 2);
   * m.set(0, 0, 5);
   */
  set(i: number, j: number, value: Expr | number): void {
    if (i < 0 || i >= this.rows || j < 0 || j >= this.cols) {
      throw new RangeError(
        `Index (${i}, ${j}) out of bounds for matrix of shape (${this.rows}, ${this.cols})`
      );
    }

    const wasm = getWasmModule();
    const expr = toExpr(value);
    const code = wasm._dense_matrix_set_basic(
      this._obj.getPtr(),
      i,
      j,
      expr.getWasmPtr()
    );
    checkException(code);
  }

  /**
   * String representation of the matrix.
   */
  toString(): string {
    return this._obj.toString();
  }

  /**
   * Check structural equality with another matrix.
   * Two matrices are equal if they have the same dimensions and all elements are equal.
   */
  equals(other: Matrix): boolean {
    if (!(other instanceof Matrix)) {
      return false;
    }
    const wasm = getWasmModule();
    return wasm._dense_matrix_eq(this._obj.getPtr(), other._obj.getPtr()) !== 0;
  }

  /**
   * Free the underlying WASM memory.
   * After calling this, the matrix becomes invalid.
   */
  free(): void {
    this._obj.free();
  }

  // ============================================================================
  // Stubs for Phase 3.1 Part 2 (Basic Operations)
  // ============================================================================

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
}

// ============================================================================
// Factory Functions
// ============================================================================

/**
 * Create an identity matrix.
 * @param n Number of rows (and columns if m not specified)
 * @param m Number of columns (optional, defaults to n)
 * @param k Diagonal offset (default 0, positive for superdiagonals, negative for subdiagonals)
 * @returns Identity-like matrix
 *
 * @example
 * eye(3);           // 3x3 identity
 * eye(3, 4);        // 3x4 matrix with 1s on main diagonal
 * eye(3, 3, 1);     // 3x3 matrix with 1s on first superdiagonal
 */
export function eye(n: number, m?: number, k: number = 0): Matrix {
  const rows = n;
  const cols = m ?? n;

  const wasm = getWasmModule();
  const matObj = createDenseMatrix();

  try {
    const code = wasm._dense_matrix_eye(matObj.getPtr(), rows, cols, k);
    checkException(code);
    return Matrix._fromDenseMatrixObject(matObj);
  } catch (e) {
    matObj.free();
    throw e;
  }
}

/**
 * Create a matrix of zeros.
 * @param rows Number of rows
 * @param cols Number of columns
 * @returns Matrix filled with zeros
 *
 * @example
 * zeros(2, 3);  // 2x3 zero matrix
 */
export function zeros(rows: number, cols: number): Matrix {
  const wasm = getWasmModule();
  const matObj = createDenseMatrix();

  try {
    const code = wasm._dense_matrix_zeros(matObj.getPtr(), rows, cols);
    checkException(code);
    return Matrix._fromDenseMatrixObject(matObj);
  } catch (e) {
    matObj.free();
    throw e;
  }
}

/**
 * Create a matrix of ones.
 * @param rows Number of rows
 * @param cols Number of columns
 * @returns Matrix filled with ones
 *
 * @example
 * ones(2, 3);  // 2x3 matrix of ones
 */
export function ones(rows: number, cols: number): Matrix {
  const wasm = getWasmModule();
  const matObj = createDenseMatrix();

  try {
    const code = wasm._dense_matrix_ones(matObj.getPtr(), rows, cols);
    checkException(code);
    return Matrix._fromDenseMatrixObject(matObj);
  } catch (e) {
    matObj.free();
    throw e;
  }
}

/**
 * Create a diagonal matrix from values.
 * @param values Elements for the diagonal
 * @param k Diagonal offset (default 0, positive for superdiagonals, negative for subdiagonals)
 * @returns Diagonal matrix
 *
 * @example
 * diag([1, 2, 3]);      // 3x3 diagonal with 1,2,3
 * diag([1, 2, 3], 1);   // 4x4 with 1,2,3 on superdiagonal
 */
export function diag(values: (Expr | number)[], k: number = 0): Matrix {
  const wasm = getWasmModule();
  const vec = new SymEngineVec();
  const matObj = createDenseMatrix();

  try {
    // Build vector of diagonal elements
    for (const val of values) {
      const expr = toExpr(val);
      wasm._vecbasic_push_back(vec.getPtr(), expr.getWasmPtr());
    }

    const code = wasm._dense_matrix_diag(matObj.getPtr(), vec.getPtr(), k);
    checkException(code);
    return Matrix._fromDenseMatrixObject(matObj);
  } catch (e) {
    matObj.free();
    throw e;
  } finally {
    vec.free();
  }
}
