import { SparseMatrix } from './base.js';
import type { CompressedSparseMatrix } from './compressed.js';
import type { SparseFormat, DIAConstructorArrays } from './types.js';
import { getWasmModule } from '../wasm-loader.js';
import { createCSR, createCOO, registerDIAFactory } from './_factory.js';

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
 * DIA (DIAgonal) sparse matrix format
 *
 * Stores matrix diagonals in a dense 2D array where each row corresponds
 * to a diagonal. The offsets array specifies which diagonal each row represents
 * (0 = main diagonal, positive = superdiagonals, negative = subdiagonals).
 *
 * This format is optimal for:
 * - Banded matrices (tridiagonal, pentadiagonal, etc.)
 * - Efficient matrix-vector multiply for banded structure
 * - Matrices where non-zeros are concentrated along diagonals
 *
 * Not optimal for:
 * - Random access to individual elements
 * - Sparse matrices with scattered non-zeros
 */
export class DIAMatrix extends SparseMatrix {
  protected _data: Float64Array; // Flattened: [n_diags * L]
  protected _offsets: Int32Array; // Diagonal offsets
  protected _shape: [number, number];
  protected _L: number; // Length of each diagonal in data array

  constructor(
    arg: DIAConstructorArrays | number[][],
    options?: { shape?: [number, number] }
  ) {
    super();

    if (typeof arg === 'object' && 'offsets' in arg) {
      // Construct from DIAConstructorArrays
      if (!options?.shape) {
        throw new Error('shape is required when constructing DIAMatrix from arrays');
      }
      this._shape = options.shape;

      // Convert offsets
      if (arg.offsets instanceof Int32Array) {
        this._offsets = new Int32Array(arg.offsets);
      } else {
        this._offsets = new Int32Array(arg.offsets);
      }

      const n_diags = this._offsets.length;
      // L is max(nrow, ncol) to accommodate any diagonal position
      this._L = Math.max(this._shape[0], this._shape[1]);

      // Convert data - can be array of arrays or array of Float64Arrays
      if (Array.isArray(arg.data) && arg.data.length > 0) {
        if (arg.data[0] instanceof Float64Array) {
          // Array of Float64Array
          this._data = new Float64Array(n_diags * this._L);
          for (let d = 0; d < n_diags; d++) {
            const diag = arg.data[d] as Float64Array;
            this._data.set(diag.subarray(0, Math.min(diag.length, this._L)), d * this._L);
          }
        } else {
          // Array of number[]
          this._data = new Float64Array(n_diags * this._L);
          for (let d = 0; d < n_diags; d++) {
            const diag = arg.data[d] as number[];
            for (let j = 0; j < Math.min(diag.length, this._L); j++) {
              this._data[d * this._L + j] = diag[j];
            }
          }
        }
      } else {
        this._data = new Float64Array(n_diags * this._L);
      }
    } else if (Array.isArray(arg)) {
      // Construct from dense 2D array
      const dense = arg as number[][];
      const nrow = dense.length;
      const ncol = nrow > 0 ? dense[0].length : 0;
      this._shape = options?.shape ?? [nrow, ncol];
      this._L = Math.max(this._shape[0], this._shape[1]);

      // Find non-zero diagonals
      // DIA format stores data[d, j] = A[j - offset[d], j]
      const diagValues = new Map<number, number[]>();

      for (let i = 0; i < nrow; i++) {
        for (let j = 0; j < ncol; j++) {
          if (dense[i][j] !== 0) {
            const k = j - i; // diagonal offset: A[i,j] -> k = j - i
            if (!diagValues.has(k)) {
              diagValues.set(k, new Array(this._L).fill(0));
            }
            // In DIA format, data[d, j] = A[j - k, j], so position is j
            diagValues.get(k)![j] = dense[i][j];
          }
        }
      }

      // Sort offsets and build arrays
      const sortedOffsets = Array.from(diagValues.keys()).sort((a, b) => a - b);
      const n_diags = sortedOffsets.length;

      this._offsets = new Int32Array(sortedOffsets);
      this._data = new Float64Array(n_diags * this._L);

      for (let d = 0; d < n_diags; d++) {
        const diag = diagValues.get(sortedOffsets[d])!;
        for (let j = 0; j < this._L; j++) {
          this._data[d * this._L + j] = diag[j];
        }
      }
    } else {
      throw new Error('Invalid argument for DIAMatrix constructor');
    }
  }

  get format(): SparseFormat {
    return 'dia';
  }

  get shape(): [number, number] {
    return [this._shape[0], this._shape[1]];
  }

  get nnz(): number {
    // Count actual non-zeros in the valid range of each diagonal
    let count = 0;
    const [nrow, ncol] = this._shape;

    for (let d = 0; d < this._offsets.length; d++) {
      const k = this._offsets[d];
      // Valid range: j >= max(0, k) and j < min(ncol, nrow + k)
      const jStart = Math.max(0, k);
      const jEnd = Math.min(ncol, nrow + k);

      for (let j = jStart; j < jEnd; j++) {
        if (this._data[d * this._L + j] !== 0) {
          count++;
        }
      }
    }

    return count;
  }

  /**
   * Get the number of diagonals stored
   */
  get ndiags(): number {
    return this._offsets.length;
  }

  /**
   * Get the diagonal offsets
   */
  get offsets(): Int32Array {
    return new Int32Array(this._offsets);
  }

  /**
   * Matrix-vector multiplication using WASM
   */
  dot(x: Float64Array): Float64Array {
    const [nrow, ncol] = this._shape;

    if (x.length !== ncol) {
      throw new Error(`Vector length ${x.length} doesn't match matrix columns ${ncol}`);
    }

    const wasm = getWasmModule();

    const offsetsPtr = toWasm(wasm, this._offsets);
    const dataPtr = toWasm(wasm, this._data);
    const xPtr = toWasm(wasm, x);
    const yPtr = wasm._malloc(nrow * 8);

    // Initialize y to zero
    for (let i = 0; i < nrow; i++) {
      wasm.HEAPF64[(yPtr >> 3) + i] = 0;
    }

    wasm._sp_dia_matvec_f64(nrow, ncol, this._offsets.length, this._L,
      offsetsPtr, dataPtr, xPtr, yPtr);

    const result = fromWasmF64(wasm, yPtr, nrow);

    wasm._free(offsetsPtr);
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
    const n_diags = this._offsets.length;

    // Estimate nnz (upper bound)
    const nnzEstimate = this.nnz;

    const wasm = getWasmModule();

    const offsetsPtr = toWasm(wasm, this._offsets);
    const dataPtr = toWasm(wasm, this._data);

    // Create diagonal order array (for sorted output)
    const order = new Int32Array(n_diags);
    for (let i = 0; i < n_diags; i++) order[i] = i;
    const orderPtr = toWasm(wasm, order);

    const indptrPtr = wasm._malloc((nrow + 1) * 4);
    const indicesPtr = wasm._malloc(nnzEstimate * 4);
    const csrDataPtr = wasm._malloc(nnzEstimate * 8);

    const actualNnz = wasm._sp_dia_tocsr_f64(nrow, ncol, n_diags, this._L,
      offsetsPtr, dataPtr, orderPtr,
      csrDataPtr, indicesPtr, indptrPtr);

    const indptr = fromWasmI32(wasm, indptrPtr, nrow + 1);
    const indices = fromWasmI32(wasm, indicesPtr, actualNnz);
    const csrData = fromWasmF64(wasm, csrDataPtr, actualNnz);

    wasm._free(offsetsPtr);
    wasm._free(dataPtr);
    wasm._free(orderPtr);
    wasm._free(indptrPtr);
    wasm._free(indicesPtr);
    wasm._free(csrDataPtr);

    return createCSR({ data: csrData, indices, indptr }, this._shape);
  }

  /**
   * Convert to COO format
   */
  tocoo(): SparseMatrix {
    const [nrow, ncol] = this._shape;
    const nnzCount = this.nnz;

    const rowArr = new Int32Array(nnzCount);
    const colArr = new Int32Array(nnzCount);
    const dataArr = new Float64Array(nnzCount);

    let ptr = 0;
    for (let d = 0; d < this._offsets.length; d++) {
      const k = this._offsets[d];
      // Valid range: j >= max(0, k) and j < min(ncol, nrow + k)
      const jStart = Math.max(0, k);
      const jEnd = Math.min(ncol, nrow + k);

      for (let j = jStart; j < jEnd; j++) {
        const i = j - k;  // row index
        const val = this._data[d * this._L + j];
        if (val !== 0) {
          rowArr[ptr] = i;
          colArr[ptr] = j;
          dataArr[ptr] = val;
          ptr++;
        }
      }
    }

    return createCOO({ row: rowArr, col: colArr, data: dataArr }, this._shape);
  }

  /**
   * Convert to CSC format (via CSR)
   */
  tocsc(): SparseMatrix {
    return this.tocsr().tocsc();
  }

  /**
   * Convert to dense 2D array
   *
   * DIA format stores data[d, j] = A[j - offset[d], j] for valid positions.
   * For superdiagonals (k > 0), first k positions are padding.
   * For subdiagonals (k < 0), last |k| positions are padding.
   */
  toArray(): number[][] {
    const [nrow, ncol] = this._shape;
    const result: number[][] = [];

    for (let i = 0; i < nrow; i++) {
      result.push(new Array(ncol).fill(0));
    }

    for (let d = 0; d < this._offsets.length; d++) {
      const k = this._offsets[d];
      // For diagonal k: A[i, j] where j = i + k, so i = j - k
      // Data is stored by column j: data[d, j] = A[j-k, j]
      // Valid range: j >= max(0, k) and j < min(ncol, nrow + k)
      const jStart = Math.max(0, k);
      const jEnd = Math.min(ncol, nrow + k);

      for (let j = jStart; j < jEnd; j++) {
        const i = j - k;  // row index
        const val = this._data[d * this._L + j];
        result[i][j] = val;
      }
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

    // Check if we have this diagonal stored
    const diagIdx = this._offsets.indexOf(k);
    if (diagIdx >= 0) {
      // For diagonal k: element i of the result is A[i, i+k] (if k>=0)
      // or A[i-k, i] (if k<0), which is stored at data[diagIdx, i+k] or data[diagIdx, i]
      const jStart = Math.max(0, k);

      for (let i = 0; i < diagLen; i++) {
        const j = jStart + i;
        result[i] = this._data[diagIdx * this._L + j];
      }
    }
    // If diagonal not stored, result is all zeros (already initialized)

    return result;
  }

  /**
   * Create a deep copy
   */
  copy(): DIAMatrix {
    const data: Float64Array[] = [];
    for (let d = 0; d < this._offsets.length; d++) {
      data.push(new Float64Array(this._data.subarray(d * this._L, (d + 1) * this._L)));
    }

    return new DIAMatrix(
      { data, offsets: new Int32Array(this._offsets) },
      { shape: [this._shape[0], this._shape[1]] }
    );
  }

  /**
   * Transpose the matrix
   *
   * For original A with offset k: A[i, j] where j = i + k, stored at data[d, j]
   * Transposed B = A^T: B[j, i] = A[i, j], which has offset -k
   * For B with offset -k: B[p, q] where q = p - k, stored at data[d, q]
   * Since B[j, i] = A[i, j] and i = j + (-k) = j - k, we have q = i
   * So data_B[d, i] = A[i, j] = data_A[d, j] where j = i + k
   */
  transpose(): DIAMatrix {
    const [nrow, ncol] = this._shape;
    const newOffsets = new Int32Array(this._offsets.length);
    const newL = Math.max(ncol, nrow);
    const newData: Float64Array[] = [];

    for (let d = 0; d < this._offsets.length; d++) {
      const k = this._offsets[d];
      const newK = -k;
      newOffsets[d] = newK;
      const newDiag = new Float64Array(newL);

      // Original: valid columns are j in [max(0,k), min(ncol, nrow+k))
      // For each valid j in original, the value A[j-k, j] goes to B[j, j-k]
      // In B (with offset newK=-k), the value at (j, j-k) is stored at data_B[d, j-k]
      // Since i = j - k, we have data_B[d, i] = data_A[d, j] where j = i + k
      const jStart = Math.max(0, k);
      const jEnd = Math.min(ncol, nrow + k);

      for (let j = jStart; j < jEnd; j++) {
        const i = j - k;  // i is the row in A, becomes column in B
        // In B, the position (j, i) has value A[i, j]
        // For B with offset -k, this is stored at data_B[d, i]
        // (since for offset -k: B[p, q] at q = p + (-k), so p = q + k, and element is at data[d, q])
        // Here p = j, q = i, so data_B[d, i] = A[i, j]
        newDiag[i] = this._data[d * this._L + j];
      }

      newData.push(newDiag);
    }

    return new DIAMatrix(
      { data: newData, offsets: newOffsets },
      { shape: [ncol, nrow] }
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
 * Factory function to create a DIA matrix
 */
export function dia_matrix(
  arg: DIAConstructorArrays | number[][],
  options?: { shape?: [number, number] }
): DIAMatrix {
  return new DIAMatrix(arg, options);
}

// Register factory so other modules can create DIA without circular imports
registerDIAFactory((arrays, opts) => new DIAMatrix(arrays, opts));
