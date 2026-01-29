import { SparseMatrix } from './base.js';
import type { SparseFormat, COOConstructorArrays } from './types.js';
import { getWasmModule } from '../wasm-loader.js';
import type { SciWasmModule } from '../wasm-types.js';
import { createCSR, registerCOOFactory } from './_factory.js';
import type { CSRMatrix } from './csr.js';

/**
 * Helper: copy a typed array into WASM heap and return the pointer.
 */
function toWasm(wasm: SciWasmModule, arr: Int32Array | Float64Array): number {
  const ptr = wasm._malloc(arr.byteLength);
  if (arr instanceof Float64Array) {
    wasm.HEAPF64.set(arr, ptr >> 3);
  } else {
    wasm.HEAP32.set(arr, ptr >> 2);
  }
  return ptr;
}

function fromWasmF64(wasm: SciWasmModule, ptr: number, len: number): Float64Array {
  const result = new Float64Array(len);
  result.set(wasm.HEAPF64.subarray(ptr >> 3, (ptr >> 3) + len));
  return result;
}

function fromWasmI32(wasm: SciWasmModule, ptr: number, len: number): Int32Array {
  const result = new Int32Array(len);
  result.set(wasm.HEAP32.subarray(ptr >> 2, (ptr >> 2) + len));
  return result;
}

/**
 * COO (COOrdinate) sparse matrix format
 *
 * Stores matrix as three arrays: row indices, column indices, and data values.
 * This is a simple format for constructing sparse matrices, but less efficient
 * for arithmetic operations compared to CSR/CSC.
 *
 * Allows duplicate (i,j) entries which are summed together.
 */
export class COOMatrix extends SparseMatrix {
  protected _row: Int32Array;
  protected _col: Int32Array;
  protected _data: Float64Array;
  protected _shape: [number, number];

  constructor(
    arg: COOConstructorArrays | number[][],
    options?: { shape?: [number, number] }
  ) {
    super();

    if (Array.isArray(arg)) {
      // Construct from dense 2D array
      const dense = arg as number[][];
      const nrow = dense.length;
      const ncol = nrow > 0 ? dense[0].length : 0;
      this._shape = options?.shape ?? [nrow, ncol];

      // Collect COO triplets
      const rows: number[] = [];
      const cols: number[] = [];
      const vals: number[] = [];
      for (let i = 0; i < nrow; i++) {
        for (let j = 0; j < ncol; j++) {
          if (dense[i][j] !== 0) {
            rows.push(i);
            cols.push(j);
            vals.push(dense[i][j]);
          }
        }
      }

      this._row = new Int32Array(rows);
      this._col = new Int32Array(cols);
      this._data = new Float64Array(vals);
    } else {
      // Construct from (row, col, data) arrays
      const arrays = arg as COOConstructorArrays;
      if (!options?.shape) {
        throw new Error('shape is required when constructing COOMatrix from arrays');
      }
      this._row = new Int32Array(arrays.row);
      this._col = new Int32Array(arrays.col);
      this._data = new Float64Array(arrays.data);
      this._shape = options.shape;
    }

    // Validate
    if (this._row.length !== this._col.length || this._row.length !== this._data.length) {
      throw new Error('row, col, and data arrays must have the same length');
    }
  }

  get format(): SparseFormat {
    return 'coo';
  }

  get shape(): [number, number] {
    return [this._shape[0], this._shape[1]];
  }

  get nnz(): number {
    return this._data.length;
  }

  /**
   * Get row indices array
   */
  get row(): Int32Array {
    return this._row;
  }

  /**
   * Get column indices array
   */
  get col(): Int32Array {
    return this._col;
  }

  /**
   * Get data values array
   */
  get data(): Float64Array {
    return this._data;
  }

  /**
   * Convert to CSR (Compressed Sparse Row) format using WASM
   */
  tocsr(): SparseMatrix {
    const wasm = getWasmModule();
    const [nrow, ncol] = this.shape;
    const nnz = this.nnz;

    // Allocate input arrays in WASM
    const rowPtr = toWasm(wasm, this._row);
    const colPtr = toWasm(wasm, this._col);
    const dataPtr = toWasm(wasm, this._data);

    // Allocate output arrays
    const indptrPtr = wasm._malloc((nrow + 1) * 4);
    const indicesPtr = wasm._malloc(nnz * 4);
    const csrDataPtr = wasm._malloc(nnz * 8);

    try {
      // Call WASM function to convert COO -> CSR
      wasm._sp_coo_tocsr_f64(
        nrow, ncol, nnz,
        rowPtr, colPtr, dataPtr,
        indptrPtr, indicesPtr, csrDataPtr
      );

      // Extract results
      const indptr = fromWasmI32(wasm, indptrPtr, nrow + 1);
      const indices = fromWasmI32(wasm, indicesPtr, nnz);
      const csrData = fromWasmF64(wasm, csrDataPtr, nnz);

      return createCSR({ data: csrData, indices, indptr }, this.shape);
    } finally {
      // Free all allocated memory
      wasm._free(rowPtr);
      wasm._free(colPtr);
      wasm._free(dataPtr);
      wasm._free(indptrPtr);
      wasm._free(indicesPtr);
      wasm._free(csrDataPtr);
    }
  }

  /**
   * Convert to CSC format (via CSR)
   */
  tocsc(): SparseMatrix {
    return this.tocsr().tocsc();
  }

  /**
   * Convert to dense 2D array using WASM
   */
  toArray(): number[][] {
    const wasm = getWasmModule();
    const [nrow, ncol] = this.shape;
    const nnz = this.nnz;

    // Allocate input arrays in WASM
    const rowPtr = toWasm(wasm, this._row);
    const colPtr = toWasm(wasm, this._col);
    const dataPtr = toWasm(wasm, this._data);

    // Allocate output dense array and initialize to zero
    const densePtr = wasm._malloc(nrow * ncol * 8);
    // Zero the output array - coo_todense accumulates values
    wasm.HEAPF64.fill(0, densePtr >> 3, (densePtr >> 3) + nrow * ncol);

    try {
      // Call WASM function to convert COO -> dense
      wasm._sp_coo_todense_f64(
        nrow, ncol, nnz,
        rowPtr, colPtr, dataPtr,
        densePtr, 0
      );

      // Extract dense array and convert to 2D
      const dense = fromWasmF64(wasm, densePtr, nrow * ncol);
      const result: number[][] = [];
      for (let i = 0; i < nrow; i++) {
        result.push(Array.from(dense.subarray(i * ncol, (i + 1) * ncol)));
      }
      return result;
    } finally {
      wasm._free(rowPtr);
      wasm._free(colPtr);
      wasm._free(dataPtr);
      wasm._free(densePtr);
    }
  }

  /**
   * Matrix-vector multiplication: y = A * x using WASM
   */
  dot(x: Float64Array): Float64Array {
    const wasm = getWasmModule();
    const [nrow, ncol] = this.shape;
    const nnz = this.nnz;

    if (x.length !== ncol) {
      throw new Error(`dimension mismatch: matrix has ${ncol} columns, vector has ${x.length} elements`);
    }

    // Allocate input arrays in WASM
    const rowPtr = toWasm(wasm, this._row);
    const colPtr = toWasm(wasm, this._col);
    const dataPtr = toWasm(wasm, this._data);
    const xPtr = toWasm(wasm, x);

    // Allocate output vector (initialized to zero)
    const yPtr = wasm._malloc(nrow * 8);
    for (let i = 0; i < nrow; i++) {
      wasm.HEAPF64[(yPtr >> 3) + i] = 0;
    }

    try {
      // Call WASM function for COO matrix-vector multiply
      wasm._sp_coo_matvec_f64(
        nnz,
        rowPtr, colPtr, dataPtr,
        xPtr, yPtr
      );

      return fromWasmF64(wasm, yPtr, nrow);
    } finally {
      wasm._free(rowPtr);
      wasm._free(colPtr);
      wasm._free(dataPtr);
      wasm._free(xPtr);
      wasm._free(yPtr);
    }
  }

  /**
   * Transpose the matrix (swap row and col arrays)
   */
  transpose(): COOMatrix {
    return new COOMatrix(
      {
        row: new Int32Array(this._col),  // Swap row <-> col
        col: new Int32Array(this._row),
        data: new Float64Array(this._data),
      },
      { shape: [this._shape[1], this._shape[0]] }  // Swap shape
    );
  }

  /**
   * Create a copy of this matrix
   */
  copy(): COOMatrix {
    return new COOMatrix(
      {
        row: new Int32Array(this._row),
        col: new Int32Array(this._col),
        data: new Float64Array(this._data),
      },
      { shape: [this._shape[0], this._shape[1]] }
    );
  }

  /**
   * Add another sparse matrix (delegates to CSR)
   */
  add(other: SparseMatrix): SparseMatrix {
    return (this.tocsr() as CSRMatrix).add(other.tocsr() as CSRMatrix);
  }

  /**
   * Subtract another sparse matrix (delegates to CSR)
   */
  subtract(other: SparseMatrix): SparseMatrix {
    return (this.tocsr() as CSRMatrix).subtract(other.tocsr() as CSRMatrix);
  }

  /**
   * Element-wise multiply (delegates to CSR)
   */
  multiply(other: SparseMatrix): SparseMatrix {
    return (this.tocsr() as CSRMatrix).multiply(other.tocsr() as CSRMatrix);
  }

  /**
   * Matrix-matrix multiplication (delegates to CSR)
   */
  matmul(other: SparseMatrix): SparseMatrix {
    return (this.tocsr() as CSRMatrix).matmul(other.tocsr() as CSRMatrix);
  }

  /**
   * Extract diagonal (delegates to CSR)
   */
  diagonal(k: number = 0): Float64Array {
    return this.tocsr().diagonal(k);
  }
}

/**
 * Factory function to create a COO matrix
 */
export function coo_matrix(
  arg: COOConstructorArrays | number[][],
  options?: { shape?: [number, number] }
): COOMatrix {
  return new COOMatrix(arg, options);
}

// Register factory so other modules can create COO without circular imports
registerCOOFactory((arrays, opts) => new COOMatrix(arrays, opts));
