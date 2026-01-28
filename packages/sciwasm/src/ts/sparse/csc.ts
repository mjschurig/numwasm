import { CompressedSparseMatrix } from './compressed.js';
import type { SparseFormat, SparseConstructorArrays } from './types.js';
import { getWasmModule } from '../wasm-loader.js';
import { registerCSCFactory, createCSR } from './_factory.js';

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

export class CSCMatrix extends CompressedSparseMatrix {
  get format(): SparseFormat { return 'csc'; }

  protected _swapPair(a: number, b: number): [number, number] {
    return [b, a]; // swap for CSC — columns are major axis
  }

  protected _swapArrays(a: number[], b: number[]): [number[], number[]] {
    return [b, a]; // swap for CSC
  }

  tocsc(): CSCMatrix {
    return this;
  }

  tocsr(): CompressedSparseMatrix {
    // CSC->CSR: csr_tocsc with swapped dimensions
    const wasm = getWasmModule();
    const [nrow, ncol] = this.shape;
    const nnz = this.nnz;

    const apPtr = toWasm(wasm, this._indptr);
    const ajPtr = toWasm(wasm, this._indices);
    const axPtr = toWasm(wasm, this._data);
    const bpPtr = wasm._malloc((nrow + 1) * 4);
    const biPtr = wasm._malloc(nnz * 4);
    const bxPtr = wasm._malloc(nnz * 8);

    // Treat CSC(nrow,ncol) as CSR(ncol,nrow) and convert
    wasm._sp_csr_tocsc_f64(ncol, nrow, apPtr, ajPtr, axPtr, bpPtr, biPtr, bxPtr);

    const Bp = fromWasmI32(wasm, bpPtr, nrow + 1);
    const Bi = fromWasmI32(wasm, biPtr, nnz);
    const Bx = fromWasmF64(wasm, bxPtr, nnz);

    wasm._free(apPtr); wasm._free(ajPtr); wasm._free(axPtr);
    wasm._free(bpPtr); wasm._free(biPtr); wasm._free(bxPtr);

    return createCSR(
      { data: Bx, indices: Bi, indptr: Bp },
      [nrow, ncol]
    );
  }

  copy(): CSCMatrix {
    return new CSCMatrix(
      {
        data: new Float64Array(this._data),
        indices: new Int32Array(this._indices),
        indptr: new Int32Array(this._indptr),
      },
      { shape: [this._shape[0], this._shape[1]] }
    );
  }

  transpose(): CompressedSparseMatrix {
    // CSC transpose = CSR with same arrays, swapped shape
    return createCSR(
      { data: this._data, indices: this._indices, indptr: this._indptr },
      [this._shape[1], this._shape[0]]
    );
  }
}

export function csc_matrix(
  arg: SparseConstructorArrays | number[][],
  options?: { shape?: [number, number] }
): CSCMatrix {
  return new CSCMatrix(arg, options);
}

// Register factory so other modules can create CSC without circular imports
registerCSCFactory((arrays, opts) => new CSCMatrix(arrays, opts));
