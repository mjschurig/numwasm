import { SparseMatrix } from './base.js';
import type { SparseConstructorArrays } from './types.js';
import { getWasmModule } from '../wasm-loader.js';
import type { SciWasmModule } from '../wasm-types.js';
import { createCSR } from './_factory.js';

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

export abstract class CompressedSparseMatrix extends SparseMatrix {
  protected _data: Float64Array;
  protected _indices: Int32Array;
  protected _indptr: Int32Array;
  protected _shape: [number, number];

  constructor(
    arg: SparseConstructorArrays | number[][],
    options?: { shape?: [number, number] }
  ) {
    super();

    if (Array.isArray(arg)) {
      // Construct from dense 2D array
      const dense = arg as number[][];
      const nrow = dense.length;
      const ncol = nrow > 0 ? dense[0].length : 0;
      this._shape = options?.shape ?? [nrow, ncol];

      // Collect COO triples
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

      // For CSR: rows are major, cols are minor
      // For CSC: cols are major, rows are minor
      const [majorIdx, minorIdx] = this._swapArrays(rows, cols);
      const [M] = this._swapPair(nrow, ncol);

      // Build compressed format directly in JS (no WASM needed for small construction)
      const nnz = vals.length;
      this._indptr = new Int32Array(M + 1);
      this._indices = new Int32Array(nnz);
      this._data = new Float64Array(nnz);

      // Count per major axis
      for (let n = 0; n < nnz; n++) {
        this._indptr[majorIdx[n] + 1]++;
      }
      // Cumsum
      for (let i = 1; i <= M; i++) {
        this._indptr[i] += this._indptr[i - 1];
      }
      // Fill
      const offset = new Int32Array(M);
      for (let n = 0; n < nnz; n++) {
        const major = majorIdx[n];
        const dest = this._indptr[major] + offset[major];
        this._indices[dest] = minorIdx[n];
        this._data[dest] = vals[n];
        offset[major]++;
      }
    } else {
      // Construct from (data, indices, indptr) arrays
      const arrays = arg as SparseConstructorArrays;
      this._data = options?.shape ? new Float64Array(arrays.data) : arrays.data;
      this._indices = options?.shape ? new Int32Array(arrays.indices) : arrays.indices;
      this._indptr = options?.shape ? new Int32Array(arrays.indptr) : arrays.indptr;
      if (!options?.shape) {
        throw new Error('shape is required when constructing from arrays');
      }
      this._shape = options.shape;
    }
  }

  // Swap pair for CSR vs CSC semantics
  protected abstract _swapPair(a: number, b: number): [number, number];
  protected abstract _swapArrays(a: number[], b: number[]): [number[], number[]];

  get data(): Float64Array { return this._data; }
  get indices(): Int32Array { return this._indices; }
  get indptr(): Int32Array { return this._indptr; }
  get shape(): [number, number] { return this._shape; }
  get nnz(): number { return this._data.length; }

  toArray(): number[][] {
    const wasm = getWasmModule();
    const [nrow, ncol] = this.shape;
    // Always convert to CSR for todense
    const csr = this.tocsr() as CompressedSparseMatrix;

    const resultSize = nrow * ncol;
    const resultPtr = wasm._malloc(resultSize * 8);
    // Zero the output
    wasm.HEAPF64.fill(0, resultPtr >> 3, (resultPtr >> 3) + resultSize);

    const apPtr = toWasm(wasm, csr._indptr);
    const ajPtr = toWasm(wasm, csr._indices);
    const axPtr = toWasm(wasm, csr._data);

    wasm._sp_csr_todense_f64(nrow, ncol, apPtr, ajPtr, axPtr, resultPtr);

    const flat = fromWasmF64(wasm, resultPtr, resultSize);
    wasm._free(apPtr);
    wasm._free(ajPtr);
    wasm._free(axPtr);
    wasm._free(resultPtr);

    // Convert to 2D array
    const result: number[][] = [];
    for (let i = 0; i < nrow; i++) {
      result.push(Array.from(flat.subarray(i * ncol, (i + 1) * ncol)));
    }
    return result;
  }

  diagonal(k: number = 0): Float64Array {
    const wasm = getWasmModule();
    const [nrow, ncol] = this.shape;
    const csr = this.tocsr() as CompressedSparseMatrix;

    const firstRow = k >= 0 ? 0 : -k;
    const firstCol = k >= 0 ? k : 0;
    const diagLen = Math.min(nrow - firstRow, ncol - firstCol);
    if (diagLen <= 0) return new Float64Array(0);

    const apPtr = toWasm(wasm, csr._indptr);
    const ajPtr = toWasm(wasm, csr._indices);
    const axPtr = toWasm(wasm, csr._data);
    const yPtr = wasm._malloc(diagLen * 8);
    wasm.HEAPF64.fill(0, yPtr >> 3, (yPtr >> 3) + diagLen);

    wasm._sp_csr_diagonal_f64(k, nrow, ncol, apPtr, ajPtr, axPtr, yPtr);

    const result = fromWasmF64(wasm, yPtr, diagLen);
    wasm._free(apPtr);
    wasm._free(ajPtr);
    wasm._free(axPtr);
    wasm._free(yPtr);
    return result;
  }

  /**
   * Matrix-vector multiply: y = A * x
   */
  dot(x: Float64Array): Float64Array {
    const wasm = getWasmModule();
    const [nrow, ncol] = this.shape;

    if (x.length !== ncol) {
      throw new Error(`Dimension mismatch: matrix is ${nrow}x${ncol}, vector has length ${x.length}`);
    }

    const apPtr = toWasm(wasm, this._indptr);
    const ajPtr = toWasm(wasm, this._indices);
    const axPtr = toWasm(wasm, this._data);
    const xxPtr = toWasm(wasm, x);
    const yxPtr = wasm._malloc(nrow * 8);
    wasm.HEAPF64.fill(0, yxPtr >> 3, (yxPtr >> 3) + nrow);

    if (this.format === 'csr') {
      wasm._sp_csr_matvec_f64(nrow, ncol, apPtr, ajPtr, axPtr, xxPtr, yxPtr);
    } else {
      wasm._sp_csc_matvec_f64(nrow, ncol, apPtr, ajPtr, axPtr, xxPtr, yxPtr);
    }

    const result = fromWasmF64(wasm, yxPtr, nrow);
    wasm._free(apPtr);
    wasm._free(ajPtr);
    wasm._free(axPtr);
    wasm._free(xxPtr);
    wasm._free(yxPtr);
    return result;
  }

  /**
   * Add two sparse matrices
   */
  add(other: CompressedSparseMatrix): CompressedSparseMatrix {
    return this._binop(other, '_sp_csr_plus_csr_f64');
  }

  /**
   * Subtract two sparse matrices
   */
  subtract(other: CompressedSparseMatrix): CompressedSparseMatrix {
    return this._binop(other, '_sp_csr_minus_csr_f64');
  }

  /**
   * Element-wise multiply
   */
  multiply(other: CompressedSparseMatrix): CompressedSparseMatrix {
    return this._binop(other, '_sp_csr_elmul_csr_f64');
  }

  protected _binop(
    other: CompressedSparseMatrix,
    wasmFn: '_sp_csr_plus_csr_f64' | '_sp_csr_minus_csr_f64' | '_sp_csr_elmul_csr_f64' | '_sp_csr_eldiv_csr_f64'
  ): CompressedSparseMatrix {
    const wasm = getWasmModule();
    const [nrow, ncol] = this.shape;

    // Convert both to CSR
    const a = this.tocsr() as CompressedSparseMatrix;
    const b = other.tocsr() as CompressedSparseMatrix;

    const maxNnz = a.nnz + b.nnz;

    const apPtr = toWasm(wasm, a._indptr);
    const ajPtr = toWasm(wasm, a._indices);
    const axPtr = toWasm(wasm, a._data);
    const bpPtr = toWasm(wasm, b._indptr);
    const bjPtr = toWasm(wasm, b._indices);
    const bxPtr = toWasm(wasm, b._data);
    const cpPtr = wasm._malloc((nrow + 1) * 4);
    const cjPtr = wasm._malloc(maxNnz * 4);
    const cxPtr = wasm._malloc(maxNnz * 8);

    wasm[wasmFn](nrow, ncol, apPtr, ajPtr, axPtr, bpPtr, bjPtr, bxPtr, cpPtr, cjPtr, cxPtr);

    const Cp = fromWasmI32(wasm, cpPtr, nrow + 1);
    const actualNnz = Cp[nrow];
    const Cj = fromWasmI32(wasm, cjPtr, actualNnz);
    const Cx = fromWasmF64(wasm, cxPtr, actualNnz);

    wasm._free(apPtr); wasm._free(ajPtr); wasm._free(axPtr);
    wasm._free(bpPtr); wasm._free(bjPtr); wasm._free(bxPtr);
    wasm._free(cpPtr); wasm._free(cjPtr); wasm._free(cxPtr);

    return createCSR(
      { data: Cx, indices: Cj, indptr: Cp },
      [nrow, ncol]
    );
  }

  /**
   * Matrix-matrix multiply (sparse * sparse -> sparse)
   */
  matmul(other: CompressedSparseMatrix): CompressedSparseMatrix {
    const wasm = getWasmModule();
    const [nrowA, ncolA] = this.shape;
    const [nrowB, ncolB] = other.shape;

    if (ncolA !== nrowB) {
      throw new Error(`Dimension mismatch: ${nrowA}x${ncolA} @ ${nrowB}x${ncolB}`);
    }

    const a = this.tocsr() as CompressedSparseMatrix;
    const b = other.tocsr() as CompressedSparseMatrix;

    const apPtr = toWasm(wasm, a._indptr);
    const ajPtr = toWasm(wasm, a._indices);
    const axPtr = toWasm(wasm, a._data);
    const bpPtr = toWasm(wasm, b._indptr);
    const bjPtr = toWasm(wasm, b._indices);
    const bxPtr = toWasm(wasm, b._data);

    // Phase 1: compute nnz(C)
    const maxNnz = wasm._sp_csr_matmat_maxnnz(nrowA, ncolB, apPtr, ajPtr, bpPtr, bjPtr);

    // Phase 2: compute C
    const cpPtr = wasm._malloc((nrowA + 1) * 4);
    const cjPtr = wasm._malloc(maxNnz * 4);
    const cxPtr = wasm._malloc(maxNnz * 8);

    wasm._sp_csr_matmat_f64(nrowA, ncolB, apPtr, ajPtr, axPtr, bpPtr, bjPtr, bxPtr, cpPtr, cjPtr, cxPtr);

    const Cp = fromWasmI32(wasm, cpPtr, nrowA + 1);
    const actualNnz = Cp[nrowA];
    const Cj = fromWasmI32(wasm, cjPtr, actualNnz);
    const Cx = fromWasmF64(wasm, cxPtr, actualNnz);

    wasm._free(apPtr); wasm._free(ajPtr); wasm._free(axPtr);
    wasm._free(bpPtr); wasm._free(bjPtr); wasm._free(bxPtr);
    wasm._free(cpPtr); wasm._free(cjPtr); wasm._free(cxPtr);

    return createCSR(
      { data: Cx, indices: Cj, indptr: Cp },
      [nrowA, ncolB]
    );
  }

  /**
   * Sort indices in place
   */
  sortIndices(): void {
    const wasm = getWasmModule();
    const [M] = this._swapPair(this._shape[0], this._shape[1]);

    const apPtr = toWasm(wasm, this._indptr);
    const ajPtr = toWasm(wasm, this._indices);
    const axPtr = toWasm(wasm, this._data);

    wasm._sp_csr_sort_indices_f64(M, apPtr, ajPtr, axPtr);

    // Read back modified indices and data
    this._indices = fromWasmI32(wasm, ajPtr, this._indices.length);
    this._data = fromWasmF64(wasm, axPtr, this._data.length);

    wasm._free(apPtr);
    wasm._free(ajPtr);
    wasm._free(axPtr);
  }

  /**
   * Sum duplicate entries
   */
  sumDuplicates(): void {
    const wasm = getWasmModule();
    const [M, N] = this._swapPair(this._shape[0], this._shape[1]);

    const apPtr = toWasm(wasm, this._indptr);
    const ajPtr = toWasm(wasm, this._indices);
    const axPtr = toWasm(wasm, this._data);

    wasm._sp_csr_sum_duplicates_f64(M, N, apPtr, ajPtr, axPtr);

    // Read back — indptr tells us the new nnz
    const newIndptr = fromWasmI32(wasm, apPtr, M + 1);
    const newNnz = newIndptr[M];
    const newIndices = fromWasmI32(wasm, ajPtr, newNnz);
    const newData = fromWasmF64(wasm, axPtr, newNnz);

    this._indptr = newIndptr;
    this._indices = newIndices;
    this._data = newData;

    wasm._free(apPtr);
    wasm._free(ajPtr);
    wasm._free(axPtr);
  }

  /**
   * Remove zero entries
   */
  eliminateZeros(): void {
    const wasm = getWasmModule();
    const [M, N] = this._swapPair(this._shape[0], this._shape[1]);

    const apPtr = toWasm(wasm, this._indptr);
    const ajPtr = toWasm(wasm, this._indices);
    const axPtr = toWasm(wasm, this._data);

    wasm._sp_csr_eliminate_zeros_f64(M, N, apPtr, ajPtr, axPtr);

    const newIndptr = fromWasmI32(wasm, apPtr, M + 1);
    const newNnz = newIndptr[M];
    const newIndices = fromWasmI32(wasm, ajPtr, newNnz);
    const newData = fromWasmF64(wasm, axPtr, newNnz);

    this._indptr = newIndptr;
    this._indices = newIndices;
    this._data = newData;

    wasm._free(apPtr);
    wasm._free(ajPtr);
    wasm._free(axPtr);
  }
}
