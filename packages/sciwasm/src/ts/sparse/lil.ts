import { SparseMatrix } from './base.js';
import type { CompressedSparseMatrix } from './compressed.js';
import type { SparseFormat, LILRow, LILConstructorArrays } from './types.js';
import { createCSR, createCOO, registerLILFactory } from './_factory.js';

/**
 * Binary search to find insertion point (bisect_left)
 * Returns index where value should be inserted to maintain sorted order
 */
function bisectLeft(arr: number[], value: number): number {
  let lo = 0;
  let hi = arr.length;

  while (lo < hi) {
    const mid = (lo + hi) >>> 1;
    if (arr[mid] < value) {
      lo = mid + 1;
    } else {
      hi = mid;
    }
  }

  return lo;
}

/**
 * LIL (List of Lists) sparse matrix format
 *
 * Stores matrix as an array of rows, where each row contains:
 * - indices: sorted array of column indices with non-zero values
 * - data: corresponding non-zero values
 *
 * This format is optimal for:
 * - Incremental row-wise construction
 * - Flexible slicing and element access
 * - Efficient conversion to CSR (rows already sorted)
 *
 * Performance characteristics:
 * - Element lookup: O(log n) via binary search
 * - Element insertion: O(n) due to array shifting
 * - Row access: O(1)
 */
export class LILMatrix extends SparseMatrix {
  protected _rows: LILRow[];
  protected _shape: [number, number];

  constructor(
    arg: LILConstructorArrays | number[][] | [number, number],
    options?: { shape?: [number, number] }
  ) {
    super();

    if (typeof arg === 'object' && 'rows' in arg) {
      // Construct from LILConstructorArrays (used by factory)
      if (!options?.shape) {
        throw new Error('shape is required when constructing LILMatrix from arrays');
      }
      const arrays = arg as LILConstructorArrays;
      this._shape = options.shape;
      this._rows = arrays.rows.map(row => ({
        indices: [...row.indices],
        data: [...row.data],
      }));
    } else if (Array.isArray(arg) && arg.length === 2 && typeof arg[0] === 'number') {
      // Construct empty matrix from shape tuple [nrow, ncol]
      this._shape = arg as [number, number];
      this._rows = [];
      for (let i = 0; i < this._shape[0]; i++) {
        this._rows.push({ indices: [], data: [] });
      }
    } else if (Array.isArray(arg)) {
      // Construct from dense 2D array
      const dense = arg as number[][];
      const nrow = dense.length;
      const ncol = nrow > 0 ? dense[0].length : 0;
      this._shape = options?.shape ?? [nrow, ncol];
      this._rows = [];

      for (let i = 0; i < nrow; i++) {
        const indices: number[] = [];
        const data: number[] = [];

        for (let j = 0; j < ncol; j++) {
          if (dense[i][j] !== 0) {
            indices.push(j);
            data.push(dense[i][j]);
          }
        }

        this._rows.push({ indices, data });
      }

      // Pad with empty rows if shape is larger
      for (let i = nrow; i < this._shape[0]; i++) {
        this._rows.push({ indices: [], data: [] });
      }
    } else {
      throw new Error('Invalid argument for LILMatrix constructor');
    }
  }

  get format(): SparseFormat {
    return 'lil';
  }

  get shape(): [number, number] {
    return [this._shape[0], this._shape[1]];
  }

  get nnz(): number {
    let count = 0;
    for (const row of this._rows) {
      count += row.data.length;
    }
    return count;
  }

  /**
   * Get the row data (indices and values) for row i
   */
  getrow(i: number): LILRow {
    // Handle negative index
    if (i < 0) i += this._shape[0];

    if (i < 0 || i >= this._shape[0]) {
      throw new Error(`Row index ${i} out of bounds for shape ${this._shape}`);
    }

    // Return a copy
    return {
      indices: [...this._rows[i].indices],
      data: [...this._rows[i].data],
    };
  }

  /**
   * Get element at (i, j). Returns 0 if not stored.
   */
  get(i: number, j: number): number {
    // Handle negative indices
    if (i < 0) i += this._shape[0];
    if (j < 0) j += this._shape[1];

    // Bounds check
    if (i < 0 || i >= this._shape[0] || j < 0 || j >= this._shape[1]) {
      throw new Error(`Index (${i}, ${j}) out of bounds for shape ${this._shape}`);
    }

    const row = this._rows[i];
    const pos = bisectLeft(row.indices, j);

    if (pos < row.indices.length && row.indices[pos] === j) {
      return row.data[pos];
    }

    return 0;
  }

  /**
   * Set element at (i, j). If value is 0, removes the entry.
   */
  set(i: number, j: number, val: number): void {
    // Handle negative indices
    if (i < 0) i += this._shape[0];
    if (j < 0) j += this._shape[1];

    // Bounds check
    if (i < 0 || i >= this._shape[0] || j < 0 || j >= this._shape[1]) {
      throw new Error(`Index (${i}, ${j}) out of bounds for shape ${this._shape}`);
    }

    const row = this._rows[i];
    const pos = bisectLeft(row.indices, j);

    if (val === 0) {
      // Remove entry if it exists
      if (pos < row.indices.length && row.indices[pos] === j) {
        row.indices.splice(pos, 1);
        row.data.splice(pos, 1);
      }
    } else {
      if (pos < row.indices.length && row.indices[pos] === j) {
        // Update existing entry
        row.data[pos] = val;
      } else {
        // Insert new entry
        row.indices.splice(pos, 0, j);
        row.data.splice(pos, 0, val);
      }
    }
  }

  /**
   * Convert to CSR format (efficient - rows already sorted)
   */
  tocsr(): SparseMatrix {
    const [nrow] = this._shape;
    const nnz = this.nnz;

    const indptr = new Int32Array(nrow + 1);
    const indices = new Int32Array(nnz);
    const data = new Float64Array(nnz);

    let ptr = 0;
    for (let i = 0; i < nrow; i++) {
      indptr[i] = ptr;
      const row = this._rows[i];
      for (let k = 0; k < row.indices.length; k++) {
        indices[ptr] = row.indices[k];
        data[ptr] = row.data[k];
        ptr++;
      }
    }
    indptr[nrow] = ptr;

    return createCSR({ data, indices, indptr }, this._shape);
  }

  /**
   * Convert to COO format
   */
  tocoo(): SparseMatrix {
    const nnz = this.nnz;
    const row = new Int32Array(nnz);
    const col = new Int32Array(nnz);
    const data = new Float64Array(nnz);

    let ptr = 0;
    for (let i = 0; i < this._rows.length; i++) {
      const rowData = this._rows[i];
      for (let k = 0; k < rowData.indices.length; k++) {
        row[ptr] = i;
        col[ptr] = rowData.indices[k];
        data[ptr] = rowData.data[k];
        ptr++;
      }
    }

    return createCOO({ row, col, data }, this._shape);
  }

  /**
   * Convert to CSC format (via CSR)
   */
  tocsc(): SparseMatrix {
    return this.tocsr().tocsc();
  }

  /**
   * Convert to dense 2D array
   */
  toArray(): number[][] {
    const [nrow, ncol] = this._shape;
    const result: number[][] = [];

    for (let i = 0; i < nrow; i++) {
      const rowArr: number[] = new Array(ncol).fill(0);
      const row = this._rows[i];

      for (let k = 0; k < row.indices.length; k++) {
        rowArr[row.indices[k]] = row.data[k];
      }

      result.push(rowArr);
    }

    return result;
  }

  /**
   * Extract k-th diagonal
   */
  diagonal(k: number = 0): Float64Array {
    const [nrow, ncol] = this._shape;

    // Calculate diagonal length
    let diagLen: number;
    if (k >= 0) {
      diagLen = Math.min(nrow, ncol - k);
    } else {
      diagLen = Math.min(nrow + k, ncol);
    }

    if (diagLen <= 0) {
      return new Float64Array(0);
    }

    const result = new Float64Array(diagLen);

    // Starting position
    const rowStart = k >= 0 ? 0 : -k;
    const colStart = k >= 0 ? k : 0;

    for (let d = 0; d < diagLen; d++) {
      result[d] = this.get(rowStart + d, colStart + d);
    }

    return result;
  }

  /**
   * Create a deep copy
   */
  copy(): LILMatrix {
    const rows: LILRow[] = this._rows.map(row => ({
      indices: [...row.indices],
      data: [...row.data],
    }));

    return new LILMatrix({ rows }, { shape: [this._shape[0], this._shape[1]] });
  }

  /**
   * Transpose the matrix
   * Note: Transposing LIL is not efficient, returns COO-based result
   */
  transpose(): SparseMatrix {
    // Build transposed as LIL
    const [nrow, ncol] = this._shape;
    const newRows: LILRow[] = [];
    for (let j = 0; j < ncol; j++) {
      newRows.push({ indices: [], data: [] });
    }

    for (let i = 0; i < nrow; i++) {
      const row = this._rows[i];
      for (let k = 0; k < row.indices.length; k++) {
        const j = row.indices[k];
        const val = row.data[k];
        // Insert into transposed row j at position i
        const pos = bisectLeft(newRows[j].indices, i);
        newRows[j].indices.splice(pos, 0, i);
        newRows[j].data.splice(pos, 0, val);
      }
    }

    return new LILMatrix({ rows: newRows }, { shape: [ncol, nrow] });
  }

  /**
   * Add another sparse matrix (delegates to CSR)
   */
  add(other: SparseMatrix): SparseMatrix {
    return (this.tocsr() as CompressedSparseMatrix).add(other.tocsr() as CompressedSparseMatrix);
  }

  /**
   * Subtract another sparse matrix (delegates to CSR)
   */
  subtract(other: SparseMatrix): SparseMatrix {
    return (this.tocsr() as CompressedSparseMatrix).subtract(other.tocsr() as CompressedSparseMatrix);
  }

  /**
   * Element-wise multiply (delegates to CSR)
   */
  multiply(other: SparseMatrix): SparseMatrix {
    return (this.tocsr() as CompressedSparseMatrix).multiply(other.tocsr() as CompressedSparseMatrix);
  }

  /**
   * Matrix-matrix multiplication (delegates to CSR)
   */
  matmul(other: SparseMatrix): SparseMatrix {
    return (this.tocsr() as CompressedSparseMatrix).matmul(other.tocsr() as CompressedSparseMatrix);
  }

  /**
   * Matrix-vector multiplication (delegates to CSR)
   */
  dot(x: Float64Array): Float64Array {
    return (this.tocsr() as CompressedSparseMatrix).dot(x);
  }
}

/**
 * Factory function to create a LIL matrix
 */
export function lil_matrix(
  arg: number[][] | [number, number],
  options?: { shape?: [number, number] }
): LILMatrix {
  return new LILMatrix(arg, options);
}

// Register factory so other modules can create LIL without circular imports
registerLILFactory((arrays, opts) => new LILMatrix(arrays, opts));
