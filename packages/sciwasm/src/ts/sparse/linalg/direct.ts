/**
 * Direct sparse linear algebra solvers using SuperLU
 *
 * Provides direct (non-iterative) methods for solving sparse linear systems:
 * - spsolve: Solve Ax = b
 * - splu: Sparse LU decomposition
 * - spilu: Incomplete LU decomposition (for preconditioning)
 * - inv: Sparse matrix inverse
 */

import type { SparseMatrix } from '../base.js';
import type { CompressedSparseMatrix } from '../compressed.js';
import type { LinearOperator } from './types.js';
import { getSuperLUModule, loadSuperLUModule } from '../../superlu-loader.js';
import type { SuperLUModule } from '../../superlu-types.js';
import { SuperLUColPerm, SuperLUStype, SuperLUDtype, SuperLUMtype, SuperLUTrans, SuperLUMiluT } from '../../superlu-types.js';
import { SingularMatrixError, DimensionMismatchError } from '../../errors.js';
import { createCSC } from '../_factory.js';

// ============================================================
// WASM Memory Helpers
// ============================================================

function toWasmF64(wasm: SuperLUModule, arr: Float64Array): number {
  const ptr = wasm._malloc(arr.byteLength);
  wasm.HEAPF64.set(arr, ptr >> 3);
  return ptr;
}

function toWasmI32(wasm: SuperLUModule, arr: Int32Array): number {
  const ptr = wasm._malloc(arr.byteLength);
  wasm.HEAP32.set(arr, ptr >> 2);
  return ptr;
}

function fromWasmF64(wasm: SuperLUModule, ptr: number, len: number): Float64Array {
  const result = new Float64Array(len);
  result.set(wasm.HEAPF64.subarray(ptr >> 3, (ptr >> 3) + len));
  return result;
}

function fromWasmI32(wasm: SuperLUModule, ptr: number, len: number): Int32Array {
  const result = new Int32Array(len);
  result.set(wasm.HEAP32.subarray(ptr >> 2, (ptr >> 2) + len));
  return result;
}

// ============================================================
// SuperLU Structure Sizes (estimated for WASM 32-bit)
// ============================================================

const SUPERMATRIX_SIZE = 40;        // SuperMatrix struct
const OPTIONS_SIZE = 80;            // superlu_options_t
const STAT_SIZE = 64;               // SuperLUStat_t
const GLU_SIZE = 128;               // GlobalLU_t

// ============================================================
// Options Types
// ============================================================

export type ColPermSpec = 'NATURAL' | 'MMD_ATA' | 'MMD_AT_PLUS_A' | 'COLAMD';
export type TransSpec = 'N' | 'T' | 'C';

export interface SpsolveOptions {
  /** Transpose mode: 'N' for Ax=b, 'T' for A^T x=b, 'C' for A^H x=b */
  trans?: TransSpec;
  /** Column permutation strategy (default: 'COLAMD') */
  permc_spec?: ColPermSpec;
}

export interface SpLUOptions {
  /** Column permutation strategy (default: 'COLAMD') */
  permc_spec?: ColPermSpec;
  /** Threshold for partial pivoting (default: 1.0) */
  diag_pivot_thresh?: number;
}

export interface SpILUOptions {
  /** Drop tolerance (default: 1e-4) */
  drop_tol?: number;
  /** Fill factor (default: 10) */
  fill_factor?: number;
  /** MILU variant: 'silu', 'smilu_1', 'smilu_2', 'smilu_3' */
  milu?: 'silu' | 'smilu_1' | 'smilu_2' | 'smilu_3';
  /** Column permutation strategy (default: 'COLAMD') */
  permc_spec?: ColPermSpec;
}

export interface InvOptions {
  /** Output format: 'sparse' or 'dense' (default: 'dense') */
  output?: 'sparse' | 'dense';
  /** Column permutation strategy (default: 'COLAMD') */
  permc_spec?: ColPermSpec;
}

// ============================================================
// Result Types
// ============================================================

/**
 * Result of sparse LU decomposition: A = P_r * L * U * P_c
 */
export interface SpLUResult {
  /** Lower triangular factor L */
  readonly L: SparseMatrix;
  /** Upper triangular factor U */
  readonly U: SparseMatrix;
  /** Row permutation vector */
  readonly perm_r: Int32Array;
  /** Column permutation vector */
  readonly perm_c: Int32Array;
  /** Shape of the original matrix */
  readonly shape: [number, number];

  /**
   * Solve A @ x = b (or A.T @ x = b) using the LU factorization
   */
  solve(b: Float64Array, trans?: TransSpec): Float64Array;

  /**
   * Free WASM resources. Call when done with the factorization.
   */
  dispose(): void;
}

/**
 * Result of incomplete LU decomposition.
 * Can be used as a preconditioner for iterative solvers.
 */
export interface SpILUResult extends LinearOperator {
  /** Approximate lower triangular factor */
  readonly L: SparseMatrix;
  /** Approximate upper triangular factor */
  readonly U: SparseMatrix;
  /** Row permutation vector */
  readonly perm_r: Int32Array;
  /** Column permutation vector */
  readonly perm_c: Int32Array;

  /**
   * Solve the approximate system (ILU) @ x = b
   */
  solve(b: Float64Array): Float64Array;

  /**
   * Free WASM resources. Call when done with the factorization.
   */
  dispose(): void;
}

// ============================================================
// Helper Functions
// ============================================================

function getColPermCode(spec?: ColPermSpec): number {
  switch (spec) {
    case 'NATURAL': return SuperLUColPerm.NATURAL;
    case 'MMD_ATA': return SuperLUColPerm.MMD_ATA;
    case 'MMD_AT_PLUS_A': return SuperLUColPerm.MMD_AT_PLUS_A;
    case 'COLAMD':
    default: return SuperLUColPerm.COLAMD;
  }
}

function getTransCode(trans?: TransSpec): number {
  switch (trans) {
    case 'T': return SuperLUTrans.TRANS;
    case 'C': return SuperLUTrans.CONJ;
    case 'N':
    default: return SuperLUTrans.NOTRANS;
  }
}

function getMiluCode(milu?: string): number {
  switch (milu) {
    case 'smilu_1': return SuperLUMiluT.SMILU_1;
    case 'smilu_2': return SuperLUMiluT.SMILU_2;
    case 'smilu_3': return SuperLUMiluT.SMILU_3;
    case 'silu':
    default: return SuperLUMiluT.SILU;
  }
}

/**
 * Convert a sparse matrix to CSC format and extract arrays
 */
function toCSCArrays(A: SparseMatrix): {
  data: Float64Array;
  indices: Int32Array;
  indptr: Int32Array;
  shape: [number, number];
} {
  const csc = A.tocsc() as CompressedSparseMatrix;
  return {
    data: csc.data,
    indices: csc.indices,
    indptr: csc.indptr,
    shape: csc.shape,
  };
}

// ============================================================
// spsolve - Solve Ax = b
// ============================================================

/**
 * Solve the sparse linear system A @ x = b using SuperLU.
 *
 * @param A - Sparse square matrix (n x n)
 * @param b - Right-hand side vector (length n) or matrix (n x nrhs)
 * @param options - Solver options
 * @returns Solution vector x with same shape as b
 * @throws {SingularMatrixError} if A is singular
 * @throws {DimensionMismatchError} if dimensions don't match
 *
 * @example
 * ```typescript
 * const A = csr_matrix([[4, 1], [1, 3]]);
 * const b = new Float64Array([1, 2]);
 * const x = await spsolve(A, b);
 * // A @ x ≈ b
 * ```
 */
export async function spsolve(
  A: SparseMatrix,
  b: Float64Array | Float64Array[],
  options?: SpsolveOptions
): Promise<Float64Array | Float64Array[]> {
  await loadSuperLUModule();
  const slu = getSuperLUModule();

  const { data, indices, indptr, shape } = toCSCArrays(A);
  const [m, n] = shape;

  if (m !== n) {
    throw new DimensionMismatchError(`Matrix must be square, got ${m}x${n}`);
  }

  // Handle single vector or multiple RHS
  const isMultipleRHS = Array.isArray(b);
  const bVectors = isMultipleRHS ? b : [b];
  const nrhs = bVectors.length;

  // Validate dimensions
  for (const bVec of bVectors) {
    if (bVec.length !== m) {
      throw new DimensionMismatchError(
        `RHS vector length ${bVec.length} doesn't match matrix size ${m}`
      );
    }
  }

  // Flatten b into column-major format for SuperLU
  const bFlat = new Float64Array(m * nrhs);
  for (let j = 0; j < nrhs; j++) {
    bFlat.set(bVectors[j], j * m);
  }

  const nnz = data.length;

  // Allocate WASM memory
  const nzvalPtr = toWasmF64(slu, data);
  const rowindPtr = toWasmI32(slu, indices);
  const colptrPtr = toWasmI32(slu, indptr);
  const bDataPtr = toWasmF64(slu, bFlat);

  const AMatPtr = slu._malloc(SUPERMATRIX_SIZE);
  const BMatPtr = slu._malloc(SUPERMATRIX_SIZE);
  const LMatPtr = slu._malloc(SUPERMATRIX_SIZE);
  const UMatPtr = slu._malloc(SUPERMATRIX_SIZE);
  const optionsPtr = slu._malloc(OPTIONS_SIZE);
  const statPtr = slu._malloc(STAT_SIZE);
  const permCPtr = slu._malloc(n * 4);
  const permRPtr = slu._malloc(m * 4);
  const infoPtr = slu._malloc(4);

  try {
    // Create SuperMatrix for A (compressed column format)
    slu._dCreate_CompCol_Matrix(
      AMatPtr, m, n, nnz,
      nzvalPtr, rowindPtr, colptrPtr,
      SuperLUStype.SLU_NC, SuperLUDtype.SLU_D, SuperLUMtype.SLU_GE
    );

    // Create SuperMatrix for B (dense, solution overwrites B)
    slu._dCreate_Dense_Matrix(
      BMatPtr, m, nrhs, bDataPtr, m,
      SuperLUStype.SLU_DN, SuperLUDtype.SLU_D, SuperLUMtype.SLU_GE
    );

    // Initialize options
    slu._set_default_options(optionsPtr);

    // Set column permutation strategy
    const colPermCode = getColPermCode(options?.permc_spec);
    // Options struct: first int is Fact, second is Equil, third is ColPerm
    slu.setValue(optionsPtr + 8, colPermCode, 'i32');

    // Set transpose if needed
    const transCode = getTransCode(options?.trans);
    // Trans is at offset 20 in options struct
    slu.setValue(optionsPtr + 20, transCode, 'i32');

    // Initialize statistics
    slu._StatInit(statPtr);

    // Solve: A @ x = b
    slu._dgssv(
      optionsPtr,
      AMatPtr,
      permCPtr,
      permRPtr,
      LMatPtr,
      UMatPtr,
      BMatPtr,
      statPtr,
      infoPtr
    );

    // Check for errors
    const info = slu.getValue(infoPtr, 'i32');
    if (info > 0) {
      throw new SingularMatrixError(info);
    } else if (info < 0) {
      throw new Error(`SuperLU argument error at position ${-info}`);
    }

    // Extract solution from B (solution is stored in B's data)
    const solution = fromWasmF64(slu, bDataPtr, m * nrhs);

    // Clean up SuperLU structures
    slu._Destroy_SuperMatrix_Store(AMatPtr);
    slu._Destroy_SuperMatrix_Store(BMatPtr);
    slu._Destroy_SuperNode_Matrix(LMatPtr);
    slu._Destroy_CompCol_Matrix(UMatPtr);
    slu._StatFree(statPtr);

    // Return in appropriate format
    if (isMultipleRHS) {
      const result: Float64Array[] = [];
      for (let j = 0; j < nrhs; j++) {
        result.push(solution.slice(j * m, (j + 1) * m));
      }
      return result;
    } else {
      return solution.slice(0, m);
    }
  } finally {
    // Free WASM memory
    slu._free(nzvalPtr);
    slu._free(rowindPtr);
    slu._free(colptrPtr);
    slu._free(bDataPtr);
    slu._free(AMatPtr);
    slu._free(BMatPtr);
    slu._free(LMatPtr);
    slu._free(UMatPtr);
    slu._free(optionsPtr);
    slu._free(statPtr);
    slu._free(permCPtr);
    slu._free(permRPtr);
    slu._free(infoPtr);
  }
}

// ============================================================
// splu - Sparse LU Decomposition
// ============================================================

/**
 * Internal class implementing SpLUResult
 */
class SpLUResultImpl implements SpLUResult {
  private _disposed = false;
  private _slu: SuperLUModule;

  // WASM pointers (kept for solve() and dispose())
  private _LMatPtr: number;
  private _UMatPtr: number;
  private _permCPtr: number;
  private _permRPtr: number;
  private _statPtr: number;

  // Cached data
  private _L: SparseMatrix | null = null;
  private _U: SparseMatrix | null = null;
  private _perm_r: Int32Array;
  private _perm_c: Int32Array;
  readonly shape: [number, number];

  constructor(
    slu: SuperLUModule,
    LMatPtr: number,
    UMatPtr: number,
    permCPtr: number,
    permRPtr: number,
    statPtr: number,
    shape: [number, number],
    perm_r: Int32Array,
    perm_c: Int32Array
  ) {
    this._slu = slu;
    this._LMatPtr = LMatPtr;
    this._UMatPtr = UMatPtr;
    this._permCPtr = permCPtr;
    this._permRPtr = permRPtr;
    this._statPtr = statPtr;
    this.shape = shape;
    this._perm_r = perm_r;
    this._perm_c = perm_c;
  }

  get L(): SparseMatrix {
    if (this._disposed) throw new Error('SpLUResult has been disposed');
    if (!this._L) {
      // TODO: Extract L from SuperLU's supernode format
      // For now, return a placeholder
      throw new Error('L extraction not yet implemented');
    }
    return this._L;
  }

  get U(): SparseMatrix {
    if (this._disposed) throw new Error('SpLUResult has been disposed');
    if (!this._U) {
      // TODO: Extract U from SuperLU's compressed column format
      throw new Error('U extraction not yet implemented');
    }
    return this._U;
  }

  get perm_r(): Int32Array {
    if (this._disposed) throw new Error('SpLUResult has been disposed');
    return this._perm_r;
  }

  get perm_c(): Int32Array {
    if (this._disposed) throw new Error('SpLUResult has been disposed');
    return this._perm_c;
  }

  solve(b: Float64Array, trans?: TransSpec): Float64Array {
    if (this._disposed) throw new Error('SpLUResult has been disposed');

    const slu = this._slu;
    const m = this.shape[0];

    if (b.length !== m) {
      throw new DimensionMismatchError(
        `RHS vector length ${b.length} doesn't match matrix size ${m}`
      );
    }

    // Create dense matrix for B
    const bCopy = new Float64Array(b);
    const bDataPtr = toWasmF64(slu, bCopy);
    const BMatPtr = slu._malloc(SUPERMATRIX_SIZE);
    const infoPtr = slu._malloc(4);

    try {
      slu._dCreate_Dense_Matrix(
        BMatPtr, m, 1, bDataPtr, m,
        SuperLUStype.SLU_DN, SuperLUDtype.SLU_D, SuperLUMtype.SLU_GE
      );

      // Solve using triangular solve
      slu._dgstrs(
        getTransCode(trans),
        this._LMatPtr,
        this._UMatPtr,
        this._permCPtr,
        this._permRPtr,
        BMatPtr,
        this._statPtr,
        infoPtr
      );

      const info = slu.getValue(infoPtr, 'i32');
      if (info !== 0) {
        throw new Error(`SuperLU triangular solve error: ${info}`);
      }

      // Extract solution
      const result = fromWasmF64(slu, bDataPtr, m);
      slu._Destroy_SuperMatrix_Store(BMatPtr);
      return result;
    } finally {
      slu._free(bDataPtr);
      slu._free(BMatPtr);
      slu._free(infoPtr);
    }
  }

  dispose(): void {
    if (this._disposed) return;
    this._disposed = true;

    const slu = this._slu;
    slu._Destroy_SuperNode_Matrix(this._LMatPtr);
    slu._Destroy_CompCol_Matrix(this._UMatPtr);
    slu._StatFree(this._statPtr);
    slu._free(this._LMatPtr);
    slu._free(this._UMatPtr);
    slu._free(this._permCPtr);
    slu._free(this._permRPtr);
    slu._free(this._statPtr);
  }
}

/**
 * Compute sparse LU decomposition of matrix A.
 *
 * Returns a factorization object that can be reused for solving
 * multiple systems with the same matrix A.
 *
 * @param A - Sparse square matrix (n x n)
 * @param options - Factorization options
 * @returns SpLUResult with L, U factors and solve method
 * @throws {SingularMatrixError} if A is singular
 *
 * @example
 * ```typescript
 * const A = csr_matrix([[4, 1], [1, 3]]);
 * const lu = await splu(A);
 *
 * const x1 = lu.solve(new Float64Array([1, 2]));
 * const x2 = lu.solve(new Float64Array([3, 4]));
 *
 * lu.dispose(); // Free resources when done
 * ```
 */
export async function splu(A: SparseMatrix, options?: SpLUOptions): Promise<SpLUResult> {
  await loadSuperLUModule();
  const slu = getSuperLUModule();

  const { data, indices, indptr, shape } = toCSCArrays(A);
  const [m, n] = shape;

  if (m !== n) {
    throw new DimensionMismatchError(`Matrix must be square, got ${m}x${n}`);
  }

  const nnz = data.length;

  // Allocate WASM memory
  const nzvalPtr = toWasmF64(slu, data);
  const rowindPtr = toWasmI32(slu, indices);
  const colptrPtr = toWasmI32(slu, indptr);

  const AMatPtr = slu._malloc(SUPERMATRIX_SIZE);
  const ACPtr = slu._malloc(SUPERMATRIX_SIZE);  // Permuted matrix
  const LMatPtr = slu._malloc(SUPERMATRIX_SIZE);
  const UMatPtr = slu._malloc(SUPERMATRIX_SIZE);
  const optionsPtr = slu._malloc(OPTIONS_SIZE);
  const statPtr = slu._malloc(STAT_SIZE);
  const GluPtr = slu._malloc(GLU_SIZE);
  const permCPtr = slu._malloc(n * 4);
  const permRPtr = slu._malloc(m * 4);
  const etreePtr = slu._malloc(n * 4);
  const infoPtr = slu._malloc(4);

  try {
    // Create SuperMatrix for A
    slu._dCreate_CompCol_Matrix(
      AMatPtr, m, n, nnz,
      nzvalPtr, rowindPtr, colptrPtr,
      SuperLUStype.SLU_NC, SuperLUDtype.SLU_D, SuperLUMtype.SLU_GE
    );

    // Initialize options
    slu._set_default_options(optionsPtr);

    // Set column permutation strategy
    const colPermCode = getColPermCode(options?.permc_spec);
    slu.setValue(optionsPtr + 8, colPermCode, 'i32');

    // Initialize statistics
    slu._StatInit(statPtr);

    // Get column permutation
    slu._get_perm_c(colPermCode, AMatPtr, permCPtr);

    // Symbolic factorization (preorder)
    slu._sp_preorder(optionsPtr, AMatPtr, permCPtr, etreePtr, ACPtr);

    // Get panel size and relaxation parameter
    const panelSize = slu._sp_ienv(1);
    const relax = slu._sp_ienv(2);

    // Numeric factorization
    slu._dgstrf(
      optionsPtr,
      ACPtr,
      relax,
      panelSize,
      etreePtr,
      0,  // work
      0,  // lwork (0 = let SuperLU allocate)
      permCPtr,
      permRPtr,
      LMatPtr,
      UMatPtr,
      GluPtr,
      statPtr,
      infoPtr
    );

    // Check for errors
    const info = slu.getValue(infoPtr, 'i32');
    if (info > 0) {
      throw new SingularMatrixError(info);
    } else if (info < 0) {
      throw new Error(`SuperLU argument error at position ${-info}`);
    }

    // Read permutation vectors
    const perm_r = fromWasmI32(slu, permRPtr, m);
    const perm_c = fromWasmI32(slu, permCPtr, n);

    // Clean up intermediate structures (but keep L, U, perms for solving)
    slu._Destroy_SuperMatrix_Store(ACPtr);
    slu._free(AMatPtr);
    slu._free(ACPtr);
    slu._free(nzvalPtr);
    slu._free(rowindPtr);
    slu._free(colptrPtr);
    slu._free(optionsPtr);
    slu._free(GluPtr);
    slu._free(etreePtr);
    slu._free(infoPtr);

    return new SpLUResultImpl(
      slu,
      LMatPtr,
      UMatPtr,
      permCPtr,
      permRPtr,
      statPtr,
      [m, n],
      perm_r,
      perm_c
    );
  } catch (e) {
    // Clean up on error
    slu._free(nzvalPtr);
    slu._free(rowindPtr);
    slu._free(colptrPtr);
    slu._free(AMatPtr);
    slu._free(ACPtr);
    slu._free(LMatPtr);
    slu._free(UMatPtr);
    slu._free(optionsPtr);
    slu._free(statPtr);
    slu._free(GluPtr);
    slu._free(permCPtr);
    slu._free(permRPtr);
    slu._free(etreePtr);
    slu._free(infoPtr);
    throw e;
  }
}

// ============================================================
// spilu - Incomplete LU Decomposition
// ============================================================

/**
 * Internal class implementing SpILUResult
 */
class SpILUResultImpl implements SpILUResult {
  private _disposed = false;
  private _slu: SuperLUModule;

  // WASM pointers
  private _LMatPtr: number;
  private _UMatPtr: number;
  private _permCPtr: number;
  private _permRPtr: number;
  private _statPtr: number;

  // Cached data
  private _L: SparseMatrix | null = null;
  private _U: SparseMatrix | null = null;
  private _perm_r: Int32Array;
  private _perm_c: Int32Array;
  readonly shape: [number, number];
  readonly dtype = 'float64';

  constructor(
    slu: SuperLUModule,
    LMatPtr: number,
    UMatPtr: number,
    permCPtr: number,
    permRPtr: number,
    statPtr: number,
    shape: [number, number],
    perm_r: Int32Array,
    perm_c: Int32Array
  ) {
    this._slu = slu;
    this._LMatPtr = LMatPtr;
    this._UMatPtr = UMatPtr;
    this._permCPtr = permCPtr;
    this._permRPtr = permRPtr;
    this._statPtr = statPtr;
    this.shape = shape;
    this._perm_r = perm_r;
    this._perm_c = perm_c;
  }

  get L(): SparseMatrix {
    if (this._disposed) throw new Error('SpILUResult has been disposed');
    if (!this._L) {
      throw new Error('L extraction not yet implemented');
    }
    return this._L;
  }

  get U(): SparseMatrix {
    if (this._disposed) throw new Error('SpILUResult has been disposed');
    if (!this._U) {
      throw new Error('U extraction not yet implemented');
    }
    return this._U;
  }

  get perm_r(): Int32Array {
    if (this._disposed) throw new Error('SpILUResult has been disposed');
    return this._perm_r;
  }

  get perm_c(): Int32Array {
    if (this._disposed) throw new Error('SpILUResult has been disposed');
    return this._perm_c;
  }

  /**
   * LinearOperator interface: apply preconditioner (solve ILU system)
   */
  matvec(x: Float64Array): Float64Array {
    return this.solve(x);
  }

  /**
   * LinearOperator interface: transpose solve
   */
  rmatvec(x: Float64Array): Float64Array {
    return this.solve(x);  // Approximate for ILU
  }

  solve(b: Float64Array): Float64Array {
    if (this._disposed) throw new Error('SpILUResult has been disposed');

    const slu = this._slu;
    const m = this.shape[0];

    if (b.length !== m) {
      throw new DimensionMismatchError(
        `RHS vector length ${b.length} doesn't match matrix size ${m}`
      );
    }

    // Create dense matrix for B
    const bCopy = new Float64Array(b);
    const bDataPtr = toWasmF64(slu, bCopy);
    const BMatPtr = slu._malloc(SUPERMATRIX_SIZE);
    const infoPtr = slu._malloc(4);

    try {
      slu._dCreate_Dense_Matrix(
        BMatPtr, m, 1, bDataPtr, m,
        SuperLUStype.SLU_DN, SuperLUDtype.SLU_D, SuperLUMtype.SLU_GE
      );

      // Solve using triangular solve
      slu._dgstrs(
        SuperLUTrans.NOTRANS,
        this._LMatPtr,
        this._UMatPtr,
        this._permCPtr,
        this._permRPtr,
        BMatPtr,
        this._statPtr,
        infoPtr
      );

      const info = slu.getValue(infoPtr, 'i32');
      if (info !== 0) {
        throw new Error(`SuperLU ILU solve error: ${info}`);
      }

      const result = fromWasmF64(slu, bDataPtr, m);
      slu._Destroy_SuperMatrix_Store(BMatPtr);
      return result;
    } finally {
      slu._free(bDataPtr);
      slu._free(BMatPtr);
      slu._free(infoPtr);
    }
  }

  dispose(): void {
    if (this._disposed) return;
    this._disposed = true;

    const slu = this._slu;
    slu._Destroy_SuperNode_Matrix(this._LMatPtr);
    slu._Destroy_CompCol_Matrix(this._UMatPtr);
    slu._StatFree(this._statPtr);
    slu._free(this._LMatPtr);
    slu._free(this._UMatPtr);
    slu._free(this._permCPtr);
    slu._free(this._permRPtr);
    slu._free(this._statPtr);
  }
}

/**
 * Compute incomplete LU factorization for use as a preconditioner.
 *
 * The result implements LinearOperator and can be passed as the M
 * parameter to iterative solvers (cg, gmres, bicgstab).
 *
 * @param A - Sparse square matrix (n x n)
 * @param options - ILU options including drop tolerance
 * @returns SpILUResult that can be used with iterative solvers
 *
 * @example
 * ```typescript
 * const A = csr_matrix([[4, 1], [1, 3]]);
 * const M = await spilu(A, { drop_tol: 1e-4 });
 *
 * // Use as preconditioner
 * const result = cg(A, b, { M });
 *
 * M.dispose(); // Free resources when done
 * ```
 */
export async function spilu(A: SparseMatrix, options?: SpILUOptions): Promise<SpILUResult> {
  await loadSuperLUModule();
  const slu = getSuperLUModule();

  const { data, indices, indptr, shape } = toCSCArrays(A);
  const [m, n] = shape;

  if (m !== n) {
    throw new DimensionMismatchError(`Matrix must be square, got ${m}x${n}`);
  }

  const nnz = data.length;
  const drop_tol = options?.drop_tol ?? 1e-4;
  const fill_factor = options?.fill_factor ?? 10;

  // Allocate WASM memory
  const nzvalPtr = toWasmF64(slu, data);
  const rowindPtr = toWasmI32(slu, indices);
  const colptrPtr = toWasmI32(slu, indptr);

  const AMatPtr = slu._malloc(SUPERMATRIX_SIZE);
  const ACPtr = slu._malloc(SUPERMATRIX_SIZE);
  const LMatPtr = slu._malloc(SUPERMATRIX_SIZE);
  const UMatPtr = slu._malloc(SUPERMATRIX_SIZE);
  const optionsPtr = slu._malloc(OPTIONS_SIZE);
  const statPtr = slu._malloc(STAT_SIZE);
  const GluPtr = slu._malloc(GLU_SIZE);
  const permCPtr = slu._malloc(n * 4);
  const permRPtr = slu._malloc(m * 4);
  const etreePtr = slu._malloc(n * 4);
  const infoPtr = slu._malloc(4);

  try {
    // Create SuperMatrix for A
    slu._dCreate_CompCol_Matrix(
      AMatPtr, m, n, nnz,
      nzvalPtr, rowindPtr, colptrPtr,
      SuperLUStype.SLU_NC, SuperLUDtype.SLU_D, SuperLUMtype.SLU_GE
    );

    // Initialize options with ILU settings
    slu._set_default_options(optionsPtr);

    // Set column permutation
    const colPermCode = getColPermCode(options?.permc_spec);
    slu.setValue(optionsPtr + 8, colPermCode, 'i32');

    // Set ILU-specific options in the options struct
    // ILU_DropTol is at offset 52 (double)
    slu.setValue(optionsPtr + 52, drop_tol, 'double');
    // ILU_FillFactor is at offset 60 (double)
    slu.setValue(optionsPtr + 60, fill_factor, 'double');
    // ILU_MILU is at offset 48 (int)
    slu.setValue(optionsPtr + 48, getMiluCode(options?.milu), 'i32');

    // Initialize statistics
    slu._StatInit(statPtr);

    // Get column permutation
    slu._get_perm_c(colPermCode, AMatPtr, permCPtr);

    // Symbolic factorization
    slu._sp_preorder(optionsPtr, AMatPtr, permCPtr, etreePtr, ACPtr);

    // Get panel size and relaxation parameter
    const panelSize = slu._sp_ienv(1);
    const relax = slu._sp_ienv(2);

    // ILU factorization
    slu._dgsitrf(
      optionsPtr,
      ACPtr,
      relax,
      panelSize,
      etreePtr,
      0,  // work
      0,  // lwork
      permCPtr,
      permRPtr,
      LMatPtr,
      UMatPtr,
      GluPtr,
      statPtr,
      infoPtr
    );

    const info = slu.getValue(infoPtr, 'i32');
    if (info > 0) {
      throw new SingularMatrixError(info);
    } else if (info < 0) {
      throw new Error(`SuperLU ILU argument error at position ${-info}`);
    }

    // Read permutation vectors
    const perm_r = fromWasmI32(slu, permRPtr, m);
    const perm_c = fromWasmI32(slu, permCPtr, n);

    // Clean up intermediate structures
    slu._Destroy_SuperMatrix_Store(ACPtr);
    slu._free(AMatPtr);
    slu._free(ACPtr);
    slu._free(nzvalPtr);
    slu._free(rowindPtr);
    slu._free(colptrPtr);
    slu._free(optionsPtr);
    slu._free(GluPtr);
    slu._free(etreePtr);
    slu._free(infoPtr);

    return new SpILUResultImpl(
      slu,
      LMatPtr,
      UMatPtr,
      permCPtr,
      permRPtr,
      statPtr,
      [m, n],
      perm_r,
      perm_c
    );
  } catch (e) {
    // Clean up on error
    slu._free(nzvalPtr);
    slu._free(rowindPtr);
    slu._free(colptrPtr);
    slu._free(AMatPtr);
    slu._free(ACPtr);
    slu._free(LMatPtr);
    slu._free(UMatPtr);
    slu._free(optionsPtr);
    slu._free(statPtr);
    slu._free(GluPtr);
    slu._free(permCPtr);
    slu._free(permRPtr);
    slu._free(etreePtr);
    slu._free(infoPtr);
    throw e;
  }
}

// ============================================================
// inv - Sparse Matrix Inverse
// ============================================================

/**
 * Compute the inverse of a sparse matrix.
 *
 * Warning: The inverse of a sparse matrix is typically dense!
 * Consider using splu + solve for solving systems instead.
 *
 * @param A - Sparse square matrix (n x n)
 * @param options - Options including output format
 * @returns Inverse matrix (dense by default)
 * @throws {SingularMatrixError} if A is singular
 *
 * @example
 * ```typescript
 * const A = csr_matrix([[4, 1], [1, 3]]);
 * const Ainv = await inv(A);
 * // A @ Ainv ≈ I
 * ```
 */
export async function inv(
  A: SparseMatrix,
  options?: InvOptions
): Promise<SparseMatrix | number[][]> {
  const [m, n] = A.shape;

  if (m !== n) {
    throw new DimensionMismatchError(`Matrix must be square, got ${m}x${n}`);
  }

  // Compute LU factorization
  const lu = await splu(A, { permc_spec: options?.permc_spec });

  try {
    // Solve for each column of identity matrix
    const invColumns: Float64Array[] = [];
    const eyeCol = new Float64Array(n);

    for (let j = 0; j < n; j++) {
      // Create j-th column of identity
      eyeCol.fill(0);
      eyeCol[j] = 1;

      // Solve A @ x_j = e_j
      const xj = lu.solve(eyeCol);
      invColumns.push(xj);
    }

    // Return as dense 2D array by default
    if (options?.output !== 'sparse') {
      const result: number[][] = [];
      for (let i = 0; i < n; i++) {
        const row: number[] = [];
        for (let j = 0; j < n; j++) {
          row.push(invColumns[j][i]);
        }
        result.push(row);
      }
      return result;
    }

    // Build sparse matrix from columns (COO format, then convert to CSC)
    const rows: number[] = [];
    const cols: number[] = [];
    const vals: number[] = [];
    const tol = 1e-15;  // Drop near-zero values

    for (let j = 0; j < n; j++) {
      for (let i = 0; i < n; i++) {
        if (Math.abs(invColumns[j][i]) > tol) {
          rows.push(i);
          cols.push(j);
          vals.push(invColumns[j][i]);
        }
      }
    }

    // Build CSC directly
    const nnz = vals.length;
    const indptr = new Int32Array(n + 1);
    const indices = new Int32Array(nnz);
    const data = new Float64Array(nnz);

    // Count entries per column
    for (let k = 0; k < nnz; k++) {
      indptr[cols[k] + 1]++;
    }
    // Cumsum
    for (let j = 1; j <= n; j++) {
      indptr[j] += indptr[j - 1];
    }
    // Fill
    const offset = new Int32Array(n);
    for (let k = 0; k < nnz; k++) {
      const j = cols[k];
      const dest = indptr[j] + offset[j];
      indices[dest] = rows[k];
      data[dest] = vals[k];
      offset[j]++;
    }

    return createCSC({ data, indices, indptr }, [n, n]);
  } finally {
    lu.dispose();
  }
}
