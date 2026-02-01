/**
 * Symbolic matrix operations.
 * @module matrices
 */

import { Expr, Integer, exprFromWasm } from '../core/index.js';
import { getWasmModule } from '../wasm-loader.js';
import {
  DenseMatrixObject,
  createDenseMatrix,
  SparseMatrixObject,
  createSparseMatrixWithSize,
  SymEngineVec,
  createBasic,
  checkException,
} from '../wasm-memory.js';

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
  // Submatrix Operations
  // ============================================================================

  /**
   * Extract a submatrix from this matrix.
   * @param r1 Starting row index (inclusive)
   * @param c1 Starting column index (inclusive)
   * @param r2 Ending row index (inclusive)
   * @param c2 Ending column index (inclusive)
   * @param rStep Row step size (default 1)
   * @param cStep Column step size (default 1)
   * @returns A new Matrix containing the extracted submatrix
   *
   * @example
   * const m = new Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]]);
   * m.submatrix(0, 0, 1, 1);  // Returns [[1, 2], [4, 5]]
   */
  submatrix(
    r1: number,
    c1: number,
    r2: number,
    c2: number,
    rStep: number = 1,
    cStep: number = 1
  ): Matrix {
    if (r1 < 0 || r1 >= this.rows || r2 < 0 || r2 >= this.rows) {
      throw new RangeError(
        `Row indices (${r1}, ${r2}) out of bounds for matrix with ${this.rows} rows`
      );
    }
    if (c1 < 0 || c1 >= this.cols || c2 < 0 || c2 >= this.cols) {
      throw new RangeError(
        `Column indices (${c1}, ${c2}) out of bounds for matrix with ${this.cols} columns`
      );
    }
    if (rStep <= 0 || cStep <= 0) {
      throw new Error('Step sizes must be positive');
    }

    const wasm = getWasmModule();
    const resultObj = createDenseMatrix();

    try {
      const code = wasm._dense_matrix_submatrix(
        resultObj.getPtr(),
        this._obj.getPtr(),
        r1,
        c1,
        r2,
        c2,
        rStep,
        cStep
      );
      checkException(code);
      return Matrix._fromDenseMatrixObject(resultObj);
    } catch (e) {
      resultObj.free();
      throw e;
    }
  }

  /**
   * Horizontally stack another matrix to the right of this one.
   * Both matrices must have the same number of rows.
   * @param other The matrix to append horizontally
   * @returns A new Matrix with columns from both matrices
   *
   * @example
   * const a = new Matrix([[1, 2], [3, 4]]);
   * const b = new Matrix([[5], [6]]);
   * a.rowJoin(b);  // Returns [[1, 2, 5], [3, 4, 6]]
   */
  rowJoin(other: Matrix): Matrix {
    if (this.rows !== other.rows) {
      throw new Error(
        `Cannot row_join matrices with different row counts: ${this.rows} vs ${other.rows}`
      );
    }

    const wasm = getWasmModule();
    const resultObj = createDenseMatrix();

    try {
      // Copy this matrix to result first (indices are inclusive)
      const code1 = wasm._dense_matrix_submatrix(
        resultObj.getPtr(),
        this._obj.getPtr(),
        0,
        0,
        this.rows - 1,
        this.cols - 1,
        1,
        1
      );
      checkException(code1);

      // Join with other
      const code2 = wasm._dense_matrix_row_join(resultObj.getPtr(), other._obj.getPtr());
      checkException(code2);

      return Matrix._fromDenseMatrixObject(resultObj);
    } catch (e) {
      resultObj.free();
      throw e;
    }
  }

  /**
   * Vertically stack another matrix below this one.
   * Both matrices must have the same number of columns.
   * @param other The matrix to append vertically
   * @returns A new Matrix with rows from both matrices
   *
   * @example
   * const a = new Matrix([[1, 2], [3, 4]]);
   * const b = new Matrix([[5, 6]]);
   * a.colJoin(b);  // Returns [[1, 2], [3, 4], [5, 6]]
   */
  colJoin(other: Matrix): Matrix {
    if (this.cols !== other.cols) {
      throw new Error(
        `Cannot col_join matrices with different column counts: ${this.cols} vs ${other.cols}`
      );
    }

    const wasm = getWasmModule();
    const resultObj = createDenseMatrix();

    try {
      // Copy this matrix to result first (indices are inclusive)
      const code1 = wasm._dense_matrix_submatrix(
        resultObj.getPtr(),
        this._obj.getPtr(),
        0,
        0,
        this.rows - 1,
        this.cols - 1,
        1,
        1
      );
      checkException(code1);

      // Join with other
      const code2 = wasm._dense_matrix_col_join(resultObj.getPtr(), other._obj.getPtr());
      checkException(code2);

      return Matrix._fromDenseMatrixObject(resultObj);
    } catch (e) {
      resultObj.free();
      throw e;
    }
  }

  /**
   * Delete a row from this matrix (in-place).
   * @param k The row index to delete (0-based)
   *
   * @example
   * const m = new Matrix([[1, 2], [3, 4], [5, 6]]);
   * m.rowDel(1);  // m becomes [[1, 2], [5, 6]]
   */
  rowDel(k: number): void {
    if (k < 0 || k >= this.rows) {
      throw new RangeError(
        `Row index ${k} out of bounds for matrix with ${this.rows} rows`
      );
    }

    const wasm = getWasmModule();
    const code = wasm._dense_matrix_row_del(this._obj.getPtr(), k);
    checkException(code);
  }

  /**
   * Delete a column from this matrix (in-place).
   * @param k The column index to delete (0-based)
   *
   * @example
   * const m = new Matrix([[1, 2, 3], [4, 5, 6]]);
   * m.colDel(1);  // m becomes [[1, 3], [4, 6]]
   */
  colDel(k: number): void {
    if (k < 0 || k >= this.cols) {
      throw new RangeError(
        `Column index ${k} out of bounds for matrix with ${this.cols} columns`
      );
    }

    const wasm = getWasmModule();
    const code = wasm._dense_matrix_col_del(this._obj.getPtr(), k);
    checkException(code);
  }

  // ============================================================================
  // Basic Operations (Phase 3.1 Part 2)
  // ============================================================================

  /**
   * Compute the determinant.
   * Only valid for square matrices.
   * @returns The determinant as a symbolic expression
   *
   * @example
   * const m = new Matrix([[1, 2], [3, 4]]);
   * m.det();  // Returns Integer(-2)
   */
  det(): Expr {
    if (this.rows !== this.cols) {
      throw new Error(
        `Determinant requires square matrix, got (${this.rows}, ${this.cols})`
      );
    }

    const wasm = getWasmModule();
    const result = createBasic();

    try {
      const code = wasm._dense_matrix_det(result.getPtr(), this._obj.getPtr());
      checkException(code);
      return exprFromWasm(result);
    } catch (e) {
      result.free();
      throw e;
    }
  }

  /**
   * Compute the inverse.
   * Only valid for square matrices with non-zero determinant.
   * @returns The inverse matrix
   *
   * @example
   * const m = new Matrix([[1, 2], [3, 4]]);
   * m.inv();  // Returns inverse matrix
   */
  inv(): Matrix {
    if (this.rows !== this.cols) {
      throw new Error(
        `Inverse requires square matrix, got (${this.rows}, ${this.cols})`
      );
    }

    const wasm = getWasmModule();
    const result = createDenseMatrix();

    try {
      const code = wasm._dense_matrix_inv(result.getPtr(), this._obj.getPtr());
      checkException(code);
      return Matrix._fromDenseMatrixObject(result);
    } catch (e) {
      result.free();
      throw e;
    }
  }

  /**
   * Compute the transpose.
   * @returns The transposed matrix
   *
   * @example
   * const m = new Matrix([[1, 2, 3], [4, 5, 6]]);
   * m.transpose();  // Returns 3x2 matrix [[1, 4], [2, 5], [3, 6]]
   */
  transpose(): Matrix {
    const wasm = getWasmModule();
    const result = createDenseMatrix();

    try {
      const code = wasm._dense_matrix_transpose(
        result.getPtr(),
        this._obj.getPtr()
      );
      checkException(code);
      return Matrix._fromDenseMatrixObject(result);
    } catch (e) {
      result.free();
      throw e;
    }
  }

  /**
   * Add another matrix.
   * Both matrices must have the same dimensions.
   * @param other The matrix to add
   * @returns The sum of the two matrices
   *
   * @example
   * const a = new Matrix([[1, 2], [3, 4]]);
   * const b = new Matrix([[5, 6], [7, 8]]);
   * a.add(b);  // Returns [[6, 8], [10, 12]]
   */
  add(other: Matrix): Matrix {
    if (this.rows !== other.rows || this.cols !== other.cols) {
      throw new Error(
        `Matrix dimensions must match for addition: (${this.rows}, ${this.cols}) vs (${other.rows}, ${other.cols})`
      );
    }

    const wasm = getWasmModule();
    const result = createDenseMatrix();

    try {
      const code = wasm._dense_matrix_add_matrix(
        result.getPtr(),
        this._obj.getPtr(),
        other._obj.getPtr()
      );
      checkException(code);
      return Matrix._fromDenseMatrixObject(result);
    } catch (e) {
      result.free();
      throw e;
    }
  }

  /**
   * Multiply by another matrix.
   * The number of columns in this matrix must equal the number of rows in other.
   * @param other The matrix to multiply by
   * @returns The product matrix
   *
   * @example
   * const a = new Matrix([[1, 2], [3, 4]]);
   * const b = new Matrix([[5, 6], [7, 8]]);
   * a.mul(b);  // Returns [[19, 22], [43, 50]]
   */
  mul(other: Matrix): Matrix {
    if (this.cols !== other.rows) {
      throw new Error(
        `Matrix dimensions incompatible for multiplication: (${this.rows}, ${this.cols}) x (${other.rows}, ${other.cols})`
      );
    }

    const wasm = getWasmModule();
    const result = createDenseMatrix();

    try {
      const code = wasm._dense_matrix_mul_matrix(
        result.getPtr(),
        this._obj.getPtr(),
        other._obj.getPtr()
      );
      checkException(code);
      return Matrix._fromDenseMatrixObject(result);
    } catch (e) {
      result.free();
      throw e;
    }
  }

  /**
   * Add a scalar to every element.
   * @param k The scalar to add (expression or number)
   * @returns New matrix with scalar added to each element
   *
   * @example
   * const m = new Matrix([[1, 2], [3, 4]]);
   * m.addScalar(10);  // Returns [[11, 12], [13, 14]]
   */
  addScalar(k: Expr | number): Matrix {
    const wasm = getWasmModule();
    const result = createDenseMatrix();
    const scalar = toExpr(k);

    try {
      const code = wasm._dense_matrix_add_scalar(
        result.getPtr(),
        this._obj.getPtr(),
        scalar.getWasmPtr()
      );
      checkException(code);
      return Matrix._fromDenseMatrixObject(result);
    } catch (e) {
      result.free();
      throw e;
    }
  }

  /**
   * Multiply every element by a scalar.
   * @param k The scalar to multiply by (expression or number)
   * @returns New matrix with each element multiplied by scalar
   *
   * @example
   * const m = new Matrix([[1, 2], [3, 4]]);
   * m.mulScalar(2);  // Returns [[2, 4], [6, 8]]
   */
  mulScalar(k: Expr | number): Matrix {
    const wasm = getWasmModule();
    const result = createDenseMatrix();
    const scalar = toExpr(k);

    try {
      const code = wasm._dense_matrix_mul_scalar(
        result.getPtr(),
        this._obj.getPtr(),
        scalar.getWasmPtr()
      );
      checkException(code);
      return Matrix._fromDenseMatrixObject(result);
    } catch (e) {
      result.free();
      throw e;
    }
  }

  // ============================================================================
  // Factorizations
  // ============================================================================

  /**
   * LU factorization of the matrix.
   * Decomposes matrix A into lower triangular (L) and upper triangular (U) matrices.
   * A = L*U
   *
   * @returns Tuple [L, U] where L is lower triangular, U is upper triangular
   * @throws Error if decomposition fails
   *
   * @example
   * const A = new Matrix([[1, 2], [3, 4]]);
   * const [L, U] = A.lu();
   */
  lu(): [Matrix, Matrix] {
    const wasm = getWasmModule();
    const lObj = createDenseMatrix();
    const uObj = createDenseMatrix();

    try {
      const code = wasm._dense_matrix_LU(
        lObj.getPtr(),
        uObj.getPtr(),
        this._obj.getPtr()
      );
      checkException(code);
      return [
        Matrix._fromDenseMatrixObject(lObj),
        Matrix._fromDenseMatrixObject(uObj),
      ];
    } catch (e) {
      lObj.free();
      uObj.free();
      throw e;
    }
  }

  /**
   * LDL factorization of the matrix.
   * Decomposes a symmetric matrix A into L*D*L^T where L is lower triangular
   * and D is diagonal.
   *
   * @returns Tuple [L, D] where L is lower triangular, D is diagonal
   * @throws Error if decomposition fails (matrix must be symmetric)
   *
   * @example
   * const A = new Matrix([[4, 2], [2, 5]]);
   * const [L, D] = A.ldl();
   */
  ldl(): [Matrix, Matrix] {
    const wasm = getWasmModule();
    const lObj = createDenseMatrix();
    const dObj = createDenseMatrix();

    try {
      const code = wasm._dense_matrix_LDL(
        lObj.getPtr(),
        dObj.getPtr(),
        this._obj.getPtr()
      );
      checkException(code);
      return [
        Matrix._fromDenseMatrixObject(lObj),
        Matrix._fromDenseMatrixObject(dObj),
      ];
    } catch (e) {
      lObj.free();
      dObj.free();
      throw e;
    }
  }

  /**
   * Fraction-free LU factorization of the matrix.
   * Returns a single matrix containing the LU decomposition without fractions.
   * Useful for symbolic computation to avoid expression swell.
   *
   * @returns Matrix containing the fraction-free LU decomposition
   * @throws Error if decomposition fails
   *
   * @example
   * const A = new Matrix([[1, 2], [3, 4]]);
   * const LU = A.fflu();
   */
  fflu(): Matrix {
    const wasm = getWasmModule();
    const luObj = createDenseMatrix();

    try {
      const code = wasm._dense_matrix_FFLU(luObj.getPtr(), this._obj.getPtr());
      checkException(code);
      return Matrix._fromDenseMatrixObject(luObj);
    } catch (e) {
      luObj.free();
      throw e;
    }
  }

  /**
   * Fraction-free LDU factorization of the matrix.
   * Decomposes matrix A into L*D*U where L is lower triangular,
   * D is diagonal, and U is upper triangular, all without fractions.
   *
   * @returns Tuple [L, D, U] where L is lower triangular, D is diagonal, U is upper triangular
   * @throws Error if decomposition fails
   *
   * @example
   * const A = new Matrix([[1, 2], [3, 4]]);
   * const [L, D, U] = A.ffldu();
   */
  ffldu(): [Matrix, Matrix, Matrix] {
    const wasm = getWasmModule();
    const lObj = createDenseMatrix();
    const dObj = createDenseMatrix();
    const uObj = createDenseMatrix();

    try {
      const code = wasm._dense_matrix_FFLDU(
        lObj.getPtr(),
        dObj.getPtr(),
        uObj.getPtr(),
        this._obj.getPtr()
      );
      checkException(code);
      return [
        Matrix._fromDenseMatrixObject(lObj),
        Matrix._fromDenseMatrixObject(dObj),
        Matrix._fromDenseMatrixObject(uObj),
      ];
    } catch (e) {
      lObj.free();
      dObj.free();
      uObj.free();
      throw e;
    }
  }

  /**
   * Solve the linear system Ax = b using LU decomposition.
   * This matrix (A) is the coefficient matrix, b is the right-hand side.
   *
   * @param b Right-hand side matrix or vector
   * @returns Solution matrix x such that A*x = b
   * @throws Error if the system cannot be solved
   *
   * @example
   * const A = new Matrix([[2, 1], [1, 3]]);
   * const b = new Matrix([[1], [2]]);
   * const x = A.luSolve(b);  // Solves Ax = b
   */
  luSolve(b: Matrix): Matrix {
    const wasm = getWasmModule();
    const xObj = createDenseMatrix();

    try {
      const code = wasm._dense_matrix_LU_solve(
        xObj.getPtr(),
        this._obj.getPtr(),
        b._obj.getPtr()
      );
      checkException(code);
      return Matrix._fromDenseMatrixObject(xObj);
    } catch (e) {
      xObj.free();
      throw e;
    }
  }

  // ============================================================================
  // Calculus
  // ============================================================================

  /**
   * Compute the elementwise derivative of the matrix with respect to a symbol.
   * Each element of the matrix is differentiated with respect to x.
   *
   * @param x The symbol to differentiate with respect to
   * @returns New matrix with each element differentiated
   *
   * @example
   * const x = new Symbol('x');
   * const m = new Matrix([[x, pow(x, 2)], [sin(x), exp(x)]]);
   * m.diff(x);  // Returns [[1, 2*x], [cos(x), exp(x)]]
   */
  diff(x: Expr): Matrix {
    const wasm = getWasmModule();
    const result = createDenseMatrix();

    try {
      const code = wasm._dense_matrix_diff(
        result.getPtr(),
        this._obj.getPtr(),
        x.getWasmPtr()
      );
      checkException(code);
      return Matrix._fromDenseMatrixObject(result);
    } catch (e) {
      result.free();
      throw e;
    }
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

/**
 * Compute the Jacobian matrix of a vector of functions with respect to variables.
 *
 * The Jacobian is a matrix of partial derivatives where element (i,j) is
 * the partial derivative of the i-th function with respect to the j-th variable.
 *
 * @param funcs Column matrix of functions (expressions)
 * @param vars Column matrix of variables (symbols) to differentiate with respect to
 * @returns The Jacobian matrix with shape (len(funcs), len(vars))
 *
 * @remarks
 * If funcs = [f1, f2, ..., fm] and vars = [x1, x2, ..., xn], then the Jacobian is:
 * ```
 * | ∂f1/∂x1  ∂f1/∂x2  ...  ∂f1/∂xn |
 * | ∂f2/∂x1  ∂f2/∂x2  ...  ∂f2/∂xn |
 * |   ...      ...    ...    ...   |
 * | ∂fm/∂x1  ∂fm/∂x2  ...  ∂fm/∂xn |
 * ```
 *
 * @example
 * ```typescript
 * const x = new Symbol('x');
 * const y = new Symbol('y');
 *
 * // Functions: [x^2 + y, x*y]
 * const funcs = new Matrix([[pow(x, 2).add(y)], [mul(x, y)]]);
 * const vars = new Matrix([[x], [y]]);
 *
 * jacobian(funcs, vars);
 * // Returns: [[2*x, 1], [y, x]]
 * ```
 *
 * @see {@link Matrix.diff} - Differentiate matrix elements
 */
export function jacobian(funcs: Matrix, vars: Matrix): Matrix {
  const wasm = getWasmModule();
  const result = createDenseMatrix();

  try {
    const code = wasm._dense_matrix_jacobian(
      result.getPtr(),
      funcs.getWasmPtr(),
      vars.getWasmPtr()
    );
    checkException(code);
    return Matrix._fromDenseMatrixObject(result);
  } catch (e) {
    result.free();
    throw e;
  }
}

// ============================================================================
// Sparse Matrix
// ============================================================================

/**
 * Symbolic sparse matrix in CSR (Compressed Sparse Row) format.
 *
 * Sparse matrices are efficient for storing matrices where most elements are zero.
 * They use less memory and can perform operations faster on sparse data.
 *
 * @remarks
 * SymEngine uses CSR format internally, which stores:
 * - Non-zero values in a contiguous array
 * - Column indices for each non-zero value
 * - Row pointers indicating where each row starts
 *
 * Use SparseMatrix when:
 * - Your matrix has many zero elements (>50% sparse)
 * - You need memory-efficient storage for large matrices
 * - You're working with structured sparsity patterns
 *
 * @example
 * ```typescript
 * // Create a 3x3 sparse matrix
 * const s = new SparseMatrix(3, 3);
 *
 * // Set non-zero elements
 * s.set(0, 0, 1);
 * s.set(1, 1, 2);
 * s.set(2, 2, 3);
 *
 * // Get elements
 * s.get(0, 0);  // 1
 * s.get(0, 1);  // 0 (default for unset elements)
 * ```
 */
export class SparseMatrix {
  /** @internal WASM pointer wrapper */
  private _obj: SparseMatrixObject;
  private _rows: number;
  private _cols: number;

  /**
   * Create a sparse matrix with given dimensions.
   * @param rows Number of rows
   * @param cols Number of columns
   *
   * @example
   * const s = new SparseMatrix(100, 100);  // 100x100 sparse matrix
   */
  constructor(rows: number, cols: number) {
    if (rows <= 0 || cols <= 0) {
      throw new Error('Matrix dimensions must be positive');
    }
    this._rows = rows;
    this._cols = cols;
    this._obj = createSparseMatrixWithSize(rows, cols);
  }

  /**
   * @internal Create a SparseMatrix from an existing SparseMatrixObject
   */
  static _fromSparseMatrixObject(
    obj: SparseMatrixObject,
    rows: number,
    cols: number
  ): SparseMatrix {
    const m = Object.create(SparseMatrix.prototype) as SparseMatrix;
    m._obj = obj;
    m._rows = rows;
    m._cols = cols;
    return m;
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
    return this._rows;
  }

  /**
   * Number of columns in the matrix.
   */
  get cols(): number {
    return this._cols;
  }

  /**
   * Shape of the matrix as [rows, cols] tuple.
   */
  get shape(): [number, number] {
    return [this._rows, this._cols];
  }

  /**
   * Get element at position (i, j).
   * Returns 0 for elements that haven't been set.
   * @param i Row index (0-based)
   * @param j Column index (0-based)
   * @returns The expression at the specified position
   *
   * @example
   * const s = new SparseMatrix(3, 3);
   * s.set(0, 0, 5);
   * s.get(0, 0);  // Returns Integer(5)
   * s.get(0, 1);  // Returns Integer(0)
   */
  get(i: number, j: number): Expr {
    if (i < 0 || i >= this._rows || j < 0 || j >= this._cols) {
      throw new RangeError(
        `Index (${i}, ${j}) out of bounds for matrix of shape (${this._rows}, ${this._cols})`
      );
    }

    const wasm = getWasmModule();
    const result = createBasic();

    try {
      const code = wasm._sparse_matrix_get_basic(
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
   * const s = new SparseMatrix(3, 3);
   * s.set(0, 0, 5);
   * s.set(1, 1, x);  // Can set symbolic values
   */
  set(i: number, j: number, value: Expr | number): void {
    if (i < 0 || i >= this._rows || j < 0 || j >= this._cols) {
      throw new RangeError(
        `Index (${i}, ${j}) out of bounds for matrix of shape (${this._rows}, ${this._cols})`
      );
    }

    const wasm = getWasmModule();
    const expr = toExpr(value);
    const code = wasm._sparse_matrix_set_basic(
      this._obj.getPtr(),
      i,
      j,
      expr.getWasmPtr()
    );
    checkException(code);
  }

  /**
   * String representation of the sparse matrix.
   */
  toString(): string {
    return this._obj.toString();
  }

  /**
   * Check structural equality with another sparse matrix.
   * Two sparse matrices are equal if they have the same dimensions and all elements are equal.
   */
  equals(other: SparseMatrix): boolean {
    if (!(other instanceof SparseMatrix)) {
      return false;
    }
    if (this._rows !== other._rows || this._cols !== other._cols) {
      return false;
    }
    const wasm = getWasmModule();
    return wasm._sparse_matrix_eq(this._obj.getPtr(), other._obj.getPtr()) !== 0;
  }

  /**
   * Free the underlying WASM memory.
   * After calling this, the matrix becomes invalid.
   */
  free(): void {
    this._obj.free();
  }
}
