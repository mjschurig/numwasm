import { SparseMatrix } from './base.js';
import type { CompressedSparseMatrix } from './compressed.js';
import type { SparseFormat } from './types.js';
import { createCOO, registerDOKFactory } from './_factory.js';

/**
 * DOK (Dictionary Of Keys) sparse matrix format
 *
 * Stores matrix entries in a Map with string keys of the form `${row},${col}`.
 * This format is optimal for:
 * - Incremental matrix construction with random access patterns
 * - Efficient element lookup and modification (O(1) average case)
 * - Small to medium matrices where random access is frequent
 *
 * Not optimal for arithmetic operations - convert to CSR/CSC first.
 */
export class DOKMatrix extends SparseMatrix {
  protected _dict: Map<string, number>;
  protected _shape: [number, number];

  constructor(
    arg: number[][] | [number, number] | Map<string, number>,
    options?: { shape?: [number, number] }
  ) {
    super();

    if (arg instanceof Map) {
      // Construct from existing Map (used by factory)
      if (!options?.shape) {
        throw new Error('shape is required when constructing DOKMatrix from Map');
      }
      this._dict = new Map(arg);
      this._shape = options.shape;
    } else if (Array.isArray(arg) && arg.length === 2 && typeof arg[0] === 'number') {
      // Construct empty matrix from shape tuple [nrow, ncol]
      this._shape = arg as [number, number];
      this._dict = new Map();
    } else if (Array.isArray(arg)) {
      // Construct from dense 2D array
      const dense = arg as number[][];
      const nrow = dense.length;
      const ncol = nrow > 0 ? dense[0].length : 0;
      this._shape = options?.shape ?? [nrow, ncol];
      this._dict = new Map();

      for (let i = 0; i < nrow; i++) {
        for (let j = 0; j < ncol; j++) {
          if (dense[i][j] !== 0) {
            this._dict.set(this._key(i, j), dense[i][j]);
          }
        }
      }
    } else {
      throw new Error('Invalid argument for DOKMatrix constructor');
    }
  }

  /**
   * Generate string key for (row, col) pair
   */
  private _key(i: number, j: number): string {
    return `${i},${j}`;
  }

  /**
   * Parse string key back to (row, col) pair
   */
  private _parseKey(key: string): [number, number] {
    const [i, j] = key.split(',').map(Number);
    return [i, j];
  }

  get format(): SparseFormat {
    return 'dok';
  }

  get shape(): [number, number] {
    return [this._shape[0], this._shape[1]];
  }

  get nnz(): number {
    return this._dict.size;
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

    return this._dict.get(this._key(i, j)) ?? 0;
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

    const key = this._key(i, j);
    if (val === 0) {
      this._dict.delete(key);
    } else {
      this._dict.set(key, val);
    }
  }

  /**
   * Check if element at (i, j) is stored (non-zero)
   */
  has(i: number, j: number): boolean {
    // Handle negative indices
    if (i < 0) i += this._shape[0];
    if (j < 0) j += this._shape[1];

    return this._dict.has(this._key(i, j));
  }

  /**
   * Delete element at (i, j)
   */
  delete(i: number, j: number): boolean {
    // Handle negative indices
    if (i < 0) i += this._shape[0];
    if (j < 0) j += this._shape[1];

    return this._dict.delete(this._key(i, j));
  }

  /**
   * Iterate over stored (row, col) pairs
   */
  *keys(): IterableIterator<[number, number]> {
    for (const key of this._dict.keys()) {
      yield this._parseKey(key);
    }
  }

  /**
   * Iterate over stored values
   */
  values(): IterableIterator<number> {
    return this._dict.values();
  }

  /**
   * Iterate over stored entries as [[row, col], value]
   */
  *entries(): IterableIterator<[[number, number], number]> {
    for (const [key, value] of this._dict.entries()) {
      yield [this._parseKey(key), value];
    }
  }

  /**
   * Convert to COO format (natural conversion for DOK)
   */
  tocoo(): SparseMatrix {
    const nnz = this.nnz;
    const row = new Int32Array(nnz);
    const col = new Int32Array(nnz);
    const data = new Float64Array(nnz);

    let idx = 0;
    for (const [key, value] of this._dict.entries()) {
      const [i, j] = this._parseKey(key);
      row[idx] = i;
      col[idx] = j;
      data[idx] = value;
      idx++;
    }

    return createCOO({ row, col, data }, this._shape);
  }

  /**
   * Convert to CSR format (via COO)
   */
  tocsr(): SparseMatrix {
    return this.tocoo().tocsr();
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
      const row: number[] = new Array(ncol).fill(0);
      result.push(row);
    }

    for (const [key, value] of this._dict.entries()) {
      const [i, j] = this._parseKey(key);
      result[i][j] = value;
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
  copy(): DOKMatrix {
    return new DOKMatrix(new Map(this._dict), { shape: [this._shape[0], this._shape[1]] });
  }

  /**
   * Transpose the matrix (swap row and col in keys)
   */
  transpose(): DOKMatrix {
    const newDict = new Map<string, number>();

    for (const [key, value] of this._dict.entries()) {
      const [i, j] = this._parseKey(key);
      newDict.set(`${j},${i}`, value);
    }

    return new DOKMatrix(newDict, { shape: [this._shape[1], this._shape[0]] });
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
 * Factory function to create a DOK matrix
 */
export function dok_matrix(
  arg: number[][] | [number, number],
  options?: { shape?: [number, number] }
): DOKMatrix {
  return new DOKMatrix(arg, options);
}

// Register factory so other modules can create DOK without circular imports
registerDOKFactory((dict, opts) => new DOKMatrix(dict, opts));
