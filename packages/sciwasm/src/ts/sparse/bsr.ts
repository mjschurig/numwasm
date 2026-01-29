import { SparseMatrix } from './base.js';
import type { CompressedSparseMatrix } from './compressed.js';
import type { SparseFormat, BSRConstructorArrays } from './types.js';
import { getWasmModule } from '../wasm-loader.js';
import { createCSR, createCOO, registerBSRFactory } from './_factory.js';

function toWasm(wasm: any, arr: Int32Array | Float64Array): number {
  const ptr = wasm._malloc(arr.byteLength);
  if (arr instanceof Float64Array) {
    wasm.HEAPF64.set(arr, ptr >> 3);
  } else {
    wasm.HEAP32.set(arr, ptr >> 2);
  }
  return ptr;
}

function fromWasmF64(wasm: any, ptr: number, len: number): Float64Array {
  const result = new Float64Array(len);
  result.set(wasm.HEAPF64.subarray(ptr >> 3, (ptr >> 3) + len));
  return result;
}

function fromWasmI32(wasm: any, ptr: number, len: number): Int32Array {
  const result = new Int32Array(len);
  result.set(wasm.HEAP32.subarray(ptr >> 2, (ptr >> 2) + len));
  return result;
}

/**
 * BSR (Block Sparse Row) sparse matrix format
 *
 * Like CSR but stores dense blocks instead of individual elements.
 * Each non-zero "entry" is actually an R×C dense block.
 *
 * Internal storage:
 * - data: Float64Array of size [nnzb * R * C] (flattened blocks in row-major order)
 * - indices: Block column indices
 * - indptr: Block row pointers
 * - blocksize: [R, C] block dimensions
 *
 * This format is optimal for:
 * - Block-structured matrices (FEM, discretized PDEs)
 * - Matrices with dense subblocks on the diagonal
 * - Efficient matrix-vector multiply for block structure
 *
 * Constraints:
 * - Matrix shape must be divisible by blocksize
 * - No random element access (use tocsr() first)
 */
export class BSRMatrix extends SparseMatrix {
  protected _data: Float64Array; // Flattened: [nnzb * R * C]
  protected _indices: Int32Array; // Block column indices
  protected _indptr: Int32Array; // Block row pointers
  protected _blocksize: [number, number];
  protected _shape: [number, number];

  constructor(
    arg: BSRConstructorArrays | number[][],
    options?: { shape?: [number, number]; blocksize?: [number, number] }
  ) {
    super();

    if (typeof arg === 'object' && 'blocksize' in arg) {
      // Construct from BSRConstructorArrays
      if (!options?.shape) {
        throw new Error('shape is required when constructing BSRMatrix from arrays');
      }
      this._shape = options.shape;
      this._blocksize = arg.blocksize;

      // Validate shape divisibility
      if (this._shape[0] % this._blocksize[0] !== 0 ||
          this._shape[1] % this._blocksize[1] !== 0) {
        throw new Error(
          `Shape ${this._shape} is not divisible by blocksize ${this._blocksize}`
        );
      }

      this._data = new Float64Array(arg.data);
      this._indices = new Int32Array(arg.indices);
      this._indptr = new Int32Array(arg.indptr);
    } else if (Array.isArray(arg)) {
      // Construct from dense 2D array
      const dense = arg as number[][];
      const nrow = dense.length;
      const ncol = nrow > 0 ? dense[0].length : 0;
      this._shape = options?.shape ?? [nrow, ncol];
      this._blocksize = options?.blocksize ?? [1, 1];

      const [R, C] = this._blocksize;

      // Validate shape divisibility
      if (this._shape[0] % R !== 0 || this._shape[1] % C !== 0) {
        throw new Error(
          `Shape ${this._shape} is not divisible by blocksize ${this._blocksize}`
        );
      }

      const n_brow = this._shape[0] / R;
      const n_bcol = this._shape[1] / C;

      // Build BSR structure by scanning blocks
      const indptr: number[] = [0];
      const indices: number[] = [];
      const blocks: number[][] = [];

      for (let bi = 0; bi < n_brow; bi++) {
        for (let bj = 0; bj < n_bcol; bj++) {
          // Check if block (bi, bj) has any non-zeros
          let hasNonzero = false;
          const block: number[] = [];

          for (let r = 0; r < R; r++) {
            for (let c = 0; c < C; c++) {
              const i = bi * R + r;
              const j = bj * C + c;
              const val = i < nrow && j < ncol ? dense[i][j] : 0;
              block.push(val);
              if (val !== 0) hasNonzero = true;
            }
          }

          if (hasNonzero) {
            indices.push(bj);
            blocks.push(block);
          }
        }
        indptr.push(indices.length);
      }

      this._indptr = new Int32Array(indptr);
      this._indices = new Int32Array(indices);
      this._data = new Float64Array(blocks.flat());
    } else {
      throw new Error('Invalid argument for BSRMatrix constructor');
    }
  }

  get format(): SparseFormat {
    return 'bsr';
  }

  get shape(): [number, number] {
    return [this._shape[0], this._shape[1]];
  }

  get blocksize(): [number, number] {
    return [this._blocksize[0], this._blocksize[1]];
  }

  /**
   * Number of stored blocks
   */
  get nnzb(): number {
    return this._indices.length;
  }

  /**
   * Number of stored elements (nnzb * R * C)
   * Note: This counts all elements in blocks, including zeros within blocks
   */
  get nnz(): number {
    // Count actual non-zeros in all blocks
    let count = 0;
    for (let i = 0; i < this._data.length; i++) {
      if (this._data[i] !== 0) count++;
    }
    return count;
  }

  /**
   * Number of block rows
   */
  get nbrows(): number {
    return this._shape[0] / this._blocksize[0];
  }

  /**
   * Number of block columns
   */
  get nbcols(): number {
    return this._shape[1] / this._blocksize[1];
  }

  /**
   * Matrix-vector multiplication using WASM
   */
  dot(x: Float64Array): Float64Array {
    const [nrow, ncol] = this._shape;
    const [R, C] = this._blocksize;
    const n_brow = nrow / R;
    const n_bcol = ncol / C;

    if (x.length !== ncol) {
      throw new Error(`Vector length ${x.length} doesn't match matrix columns ${ncol}`);
    }

    const wasm = getWasmModule();

    const indptrPtr = toWasm(wasm, this._indptr);
    const indicesPtr = toWasm(wasm, this._indices);
    const dataPtr = toWasm(wasm, this._data);
    const xPtr = toWasm(wasm, x);
    const yPtr = wasm._malloc(nrow * 8);

    // Initialize y to zero
    for (let i = 0; i < nrow; i++) {
      wasm.HEAPF64[(yPtr >> 3) + i] = 0;
    }

    wasm._sp_bsr_matvec_f64(n_brow, n_bcol, R, C,
      indptrPtr, indicesPtr, dataPtr, xPtr, yPtr);

    const result = fromWasmF64(wasm, yPtr, nrow);

    wasm._free(indptrPtr);
    wasm._free(indicesPtr);
    wasm._free(dataPtr);
    wasm._free(xPtr);
    wasm._free(yPtr);

    return result;
  }

  /**
   * Convert to CSR format using WASM
   */
  tocsr(): SparseMatrix {
    const [nrow, ncol] = this._shape;
    const [R, C] = this._blocksize;
    const n_brow = nrow / R;
    const n_bcol = ncol / C;

    // CSR nnz = nnzb * R * C
    const csrNnz = this.nnzb * R * C;

    const wasm = getWasmModule();

    const indptrPtr = toWasm(wasm, this._indptr);
    const indicesPtr = toWasm(wasm, this._indices);
    const dataPtr = toWasm(wasm, this._data);

    const bpPtr = wasm._malloc((nrow + 1) * 4);
    const bjPtr = wasm._malloc(csrNnz * 4);
    const bxPtr = wasm._malloc(csrNnz * 8);

    wasm._sp_bsr_tocsr_f64(n_brow, n_bcol, R, C,
      indptrPtr, indicesPtr, dataPtr,
      bpPtr, bjPtr, bxPtr);

    const Bp = fromWasmI32(wasm, bpPtr, nrow + 1);
    const Bj = fromWasmI32(wasm, bjPtr, csrNnz);
    const Bx = fromWasmF64(wasm, bxPtr, csrNnz);

    wasm._free(indptrPtr);
    wasm._free(indicesPtr);
    wasm._free(dataPtr);
    wasm._free(bpPtr);
    wasm._free(bjPtr);
    wasm._free(bxPtr);

    return createCSR({ data: Bx, indices: Bj, indptr: Bp }, this._shape);
  }

  /**
   * Convert to COO format
   */
  tocoo(): SparseMatrix {
    const [R, C] = this._blocksize;
    const RC = R * C;

    // Build COO arrays
    const rows: number[] = [];
    const cols: number[] = [];
    const data: number[] = [];

    const n_brow = this._indptr.length - 1;

    for (let bi = 0; bi < n_brow; bi++) {
      for (let jj = this._indptr[bi]; jj < this._indptr[bi + 1]; jj++) {
        const bj = this._indices[jj];
        const blockStart = jj * RC;

        // Expand block to individual elements
        for (let r = 0; r < R; r++) {
          for (let c = 0; c < C; c++) {
            const val = this._data[blockStart + r * C + c];
            if (val !== 0) {
              rows.push(bi * R + r);
              cols.push(bj * C + c);
              data.push(val);
            }
          }
        }
      }
    }

    return createCOO(
      { row: new Int32Array(rows), col: new Int32Array(cols), data: new Float64Array(data) },
      this._shape
    );
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
    const [R, C] = this._blocksize;
    const RC = R * C;

    const result: number[][] = [];
    for (let i = 0; i < nrow; i++) {
      result.push(new Array(ncol).fill(0));
    }

    const n_brow = this._indptr.length - 1;

    for (let bi = 0; bi < n_brow; bi++) {
      for (let jj = this._indptr[bi]; jj < this._indptr[bi + 1]; jj++) {
        const bj = this._indices[jj];
        const blockStart = jj * RC;

        // Copy block to result
        for (let r = 0; r < R; r++) {
          for (let c = 0; c < C; c++) {
            result[bi * R + r][bj * C + c] = this._data[blockStart + r * C + c];
          }
        }
      }
    }

    return result;
  }

  /**
   * Extract k-th diagonal using WASM
   */
  diagonal(k: number = 0): Float64Array {
    const [nrow, ncol] = this._shape;
    const [R, C] = this._blocksize;
    const n_brow = nrow / R;
    const n_bcol = ncol / C;

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

    const wasm = getWasmModule();

    const indptrPtr = toWasm(wasm, this._indptr);
    const indicesPtr = toWasm(wasm, this._indices);
    const dataPtr = toWasm(wasm, this._data);
    const yPtr = wasm._malloc(diagLen * 8);

    // Initialize to zero
    for (let i = 0; i < diagLen; i++) {
      wasm.HEAPF64[(yPtr >> 3) + i] = 0;
    }

    wasm._sp_bsr_diagonal_f64(k, n_brow, n_bcol, R, C,
      indptrPtr, indicesPtr, dataPtr, yPtr);

    const result = fromWasmF64(wasm, yPtr, diagLen);

    wasm._free(indptrPtr);
    wasm._free(indicesPtr);
    wasm._free(dataPtr);
    wasm._free(yPtr);

    return result;
  }

  /**
   * Transpose the matrix using WASM
   */
  transpose(): BSRMatrix {
    const [nrow, ncol] = this._shape;
    const [R, C] = this._blocksize;
    const n_brow = nrow / R;
    const n_bcol = ncol / C;

    const wasm = getWasmModule();

    const indptrPtr = toWasm(wasm, this._indptr);
    const indicesPtr = toWasm(wasm, this._indices);
    const dataPtr = toWasm(wasm, this._data);

    // Output arrays - transposed has same nnzb but swapped dimensions
    const bpPtr = wasm._malloc((n_bcol + 1) * 4);
    const biPtr = wasm._malloc(this.nnzb * 4);
    const bxPtr = wasm._malloc(this._data.length * 8);

    wasm._sp_bsr_transpose_f64(n_brow, n_bcol, R, C,
      indptrPtr, indicesPtr, dataPtr,
      bpPtr, biPtr, bxPtr);

    const Bp = fromWasmI32(wasm, bpPtr, n_bcol + 1);
    const Bi = fromWasmI32(wasm, biPtr, this.nnzb);
    const Bx = fromWasmF64(wasm, bxPtr, this._data.length);

    wasm._free(indptrPtr);
    wasm._free(indicesPtr);
    wasm._free(dataPtr);
    wasm._free(bpPtr);
    wasm._free(biPtr);
    wasm._free(bxPtr);

    return new BSRMatrix(
      { data: Bx, indices: Bi, indptr: Bp, blocksize: [C, R] },
      { shape: [ncol, nrow] }
    );
  }

  /**
   * Create a deep copy
   */
  copy(): BSRMatrix {
    return new BSRMatrix(
      {
        data: new Float64Array(this._data),
        indices: new Int32Array(this._indices),
        indptr: new Int32Array(this._indptr),
        blocksize: [this._blocksize[0], this._blocksize[1]],
      },
      { shape: [this._shape[0], this._shape[1]] }
    );
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
}

/**
 * Factory function to create a BSR matrix
 */
export function bsr_matrix(
  arg: BSRConstructorArrays | number[][],
  options?: { shape?: [number, number]; blocksize?: [number, number] }
): BSRMatrix {
  return new BSRMatrix(arg, options);
}

// Register factory so other modules can create BSR without circular imports
registerBSRFactory((arrays, opts) => new BSRMatrix(arrays, opts));
