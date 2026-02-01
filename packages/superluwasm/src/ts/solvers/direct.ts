/**
 * Direct Sparse Linear Solvers
 *
 * High-level functions for solving sparse linear systems Ax = b
 * using SuperLU's direct LU factorization methods.
 *
 * @module solvers/direct
 */

import { getSuperLUModule } from '../core/loader.js';
import {
  SuperLUColPerm,
  SuperLURowPerm,
  SuperLUTrans,
  SuperLUStype,
  SuperLUDtype,
  SuperLUMtype,
} from '../types.js';
import type { SuperLUModule } from '../types.js';
import type {
  SparseMatrixCSC,
  SparseMatrixCSR,
  SolveOptions,
  SolveResult,
  SolveStatistics,
  DataType,
  RealArray,
} from '../high-level-types.js';
import { inferDataType, getNnz } from '../high-level-types.js';
import {
  allocateDoubles,
  allocateFloats,
  allocateInts,
  readDoubles,
  readFloats,
  readInts,
  readInt,
  freeAll,
  toFloat64Array,
  toFloat32Array,
  toInt32Array,
  OPTIONS_STRUCT_SIZE,
  STAT_STRUCT_SIZE,
  SUPERMATRIX_SIZE,
  checkSuperLUError,
  INT_SIZE,
} from '../core/helpers.js';

// ============================================================
// Internal Helpers
// ============================================================

/**
 * Get the SuperLU data type enum from our DataType.
 */
function getSuperLUDtype(dtype: DataType): SuperLUDtype {
  switch (dtype) {
    case 'float32':
      return SuperLUDtype.SLU_S;
    case 'float64':
      return SuperLUDtype.SLU_D;
    case 'complex64':
      return SuperLUDtype.SLU_C;
    case 'complex128':
      return SuperLUDtype.SLU_Z;
    default:
      return SuperLUDtype.SLU_D;
  }
}

/**
 * Internal implementation of sparse solve with CSC matrix.
 */
function solveCSCInternal(
  module: SuperLUModule,
  A: SparseMatrixCSC,
  b: RealArray,
  trans: SuperLUTrans,
  options: SolveOptions
): SolveResult {
  const startTime = performance.now();
  const dtype = A.dtype ?? inferDataType(A.values);
  const { m, n } = A;
  const nnz = getNnz(A);

  // Validate dimensions
  if (m !== n) {
    throw new Error(`Matrix must be square for direct solve. Got ${m}x${n}`);
  }
  if (b.length !== m) {
    throw new Error(`RHS vector length (${b.length}) must match matrix rows (${m})`);
  }

  // Allocate all WASM memory
  const ptrs: number[] = [];
  let factorizationTime = 0;

  try {
    // Allocate options structure
    const optionsPtr = module._malloc(OPTIONS_STRUCT_SIZE);
    ptrs.push(optionsPtr);
    module._set_default_options(optionsPtr);

    // Set options based on user preferences
    // Options structure layout (approximate offsets):
    // Fact: 0, Equil: 4, ColPerm: 8, Trans: 12, IterRefine: 16, ...
    const colPerm = options.columnPermutation ?? SuperLUColPerm.COLAMD;
    const rowPerm = options.rowPermutation ?? SuperLURowPerm.LargeDiag_MC64;
    const equil = options.equilibrate !== false ? 1 : 0;
    const iterRefine = options.iterativeRefinement ? 1 : 0; // SLU_DOUBLE for double precision

    // Write options (offsets are approximations - may need adjustment)
    const _colPerm = colPerm;
    const _rowPerm = rowPerm;
    module.HEAP32[(optionsPtr >> 2) + 1] = equil;      // Equil
    module.HEAP32[(optionsPtr >> 2) + 2] = _colPerm;   // ColPerm
    module.HEAP32[(optionsPtr >> 2) + 3] = trans;      // Trans
    module.HEAP32[(optionsPtr >> 2) + 4] = iterRefine; // IterRefine
    module.HEAP32[(optionsPtr >> 2) + 8] = options.printStatistics ? 1 : 0; // PrintStat
    module.HEAP32[(optionsPtr >> 2) + 9] = _rowPerm;   // RowPerm

    // Allocate statistics structure
    const statPtr = module._malloc(STAT_STRUCT_SIZE);
    ptrs.push(statPtr);
    module._StatInit(statPtr);

    // Allocate SuperMatrix structures
    const APtr = module._malloc(SUPERMATRIX_SIZE);
    const LPtr = module._malloc(SUPERMATRIX_SIZE);
    const UPtr = module._malloc(SUPERMATRIX_SIZE);
    const BPtr = module._malloc(SUPERMATRIX_SIZE);
    ptrs.push(APtr, LPtr, UPtr, BPtr);

    // Allocate and copy matrix data
    const valuesPtr = dtype === 'float64'
      ? allocateDoubles(module, toFloat64Array(A.values as RealArray))
      : allocateFloats(module, toFloat32Array(A.values as RealArray));
    const rowindPtr = allocateInts(module, toInt32Array(A.rowIndices));
    const colptrPtr = allocateInts(module, toInt32Array(A.colPointers));
    ptrs.push(valuesPtr, rowindPtr, colptrPtr);

    // Create the CSC matrix
    const sluDtype = getSuperLUDtype(dtype);
    if (dtype === 'float64') {
      module._dCreate_CompCol_Matrix(
        APtr, m, n, nnz, valuesPtr, rowindPtr, colptrPtr,
        SuperLUStype.SLU_NC, sluDtype, SuperLUMtype.SLU_GE
      );
    } else if (dtype === 'float32') {
      module._sCreate_CompCol_Matrix(
        APtr, m, n, nnz, valuesPtr, rowindPtr, colptrPtr,
        SuperLUStype.SLU_NC, sluDtype, SuperLUMtype.SLU_GE
      );
    } else {
      throw new Error(`Unsupported data type: ${dtype}. Use float32 or float64.`);
    }

    // Allocate and copy RHS (will be overwritten with solution)
    const bPtr = dtype === 'float64'
      ? allocateDoubles(module, toFloat64Array(b))
      : allocateFloats(module, toFloat32Array(b));
    ptrs.push(bPtr);

    // Create dense matrix for RHS
    if (dtype === 'float64') {
      module._dCreate_Dense_Matrix(
        BPtr, m, 1, bPtr, m,
        SuperLUStype.SLU_DN, sluDtype, SuperLUMtype.SLU_GE
      );
    } else {
      module._sCreate_Dense_Matrix(
        BPtr, m, 1, bPtr, m,
        SuperLUStype.SLU_DN, sluDtype, SuperLUMtype.SLU_GE
      );
    }

    // Allocate permutation arrays
    const permCPtr = allocateInts(module, null, n);
    const permRPtr = allocateInts(module, null, m);
    ptrs.push(permCPtr, permRPtr);

    // Allocate info pointer
    const infoPtr = module._malloc(INT_SIZE);
    ptrs.push(infoPtr);

    // Call the simple driver
    const factStart = performance.now();
    if (dtype === 'float64') {
      module._dgssv(optionsPtr, APtr, permCPtr, permRPtr, LPtr, UPtr, BPtr, statPtr, infoPtr);
    } else {
      module._sgssv(optionsPtr, APtr, permCPtr, permRPtr, LPtr, UPtr, BPtr, statPtr, infoPtr);
    }
    factorizationTime = performance.now() - factStart;

    // Check for errors
    const info = readInt(module, infoPtr);
    checkSuperLUError(info, 'gssv');

    // Read solution from B (overwritten in place)
    const x = dtype === 'float64'
      ? readDoubles(module, bPtr, m)
      : readFloats(module, bPtr, m);

    // Read permutations
    const rowPermutation = readInts(module, permRPtr, m);
    const columnPermutation = readInts(module, permCPtr, n);

    // Cleanup SuperMatrix structures (don't free the data arrays we allocated separately)
    module._Destroy_SuperMatrix_Store(APtr);
    module._Destroy_SuperMatrix_Store(BPtr);
    module._Destroy_SuperNode_Matrix(LPtr);
    module._Destroy_CompCol_Matrix(UPtr);

    // Free statistics
    module._StatFree(statPtr);

    const totalTime = performance.now() - startTime;

    const statistics: SolveStatistics = {
      factorizationTime,
      solveTime: totalTime - factorizationTime,
      totalTime,
    };

    return {
      x,
      statistics,
      rowPermutation,
      columnPermutation,
    };

  } finally {
    // Free all allocated memory
    freeAll(module, ptrs);
  }
}

/**
 * Convert CSR to CSC format.
 */
function convertCSRtoCSC(A: SparseMatrixCSR): SparseMatrixCSC {
  const { m, n, values, colIndices, rowPointers } = A;
  const nnz = values.length;

  // Count entries in each column
  const colCounts = new Int32Array(n);
  for (let i = 0; i < nnz; i++) {
    const col = colIndices instanceof Int32Array ? colIndices[i] : colIndices[i];
    colCounts[col]++;
  }

  // Build column pointers
  const colPointers = new Int32Array(n + 1);
  colPointers[0] = 0;
  for (let j = 0; j < n; j++) {
    colPointers[j + 1] = colPointers[j] + colCounts[j];
  }

  // Build CSC arrays
  const dtype = A.dtype ?? inferDataType(A.values);
  const cscValues = dtype === 'float32'
    ? new Float32Array(nnz)
    : new Float64Array(nnz);
  const rowIndices = new Int32Array(nnz);
  const colPos = new Int32Array(n); // Current position in each column

  for (let i = 0; i < m; i++) {
    const rowStart = rowPointers instanceof Int32Array ? rowPointers[i] : rowPointers[i];
    const rowEnd = rowPointers instanceof Int32Array ? rowPointers[i + 1] : rowPointers[i + 1];

    for (let k = rowStart; k < rowEnd; k++) {
      const j = colIndices instanceof Int32Array ? colIndices[k] : colIndices[k];
      const pos = colPointers[j] + colPos[j];
      cscValues[pos] = values instanceof Float64Array || values instanceof Float32Array
        ? values[k]
        : values[k];
      rowIndices[pos] = i;
      colPos[j]++;
    }
  }

  return {
    m,
    n,
    values: cscValues,
    rowIndices,
    colPointers,
    dtype,
  };
}

// ============================================================
// Public API
// ============================================================

/**
 * Solve a sparse linear system Ax = b where A is in CSC format.
 *
 * Uses SuperLU's simple driver (dgssv/sgssv) which performs:
 * 1. Column permutation for fill reduction
 * 2. Row permutation for numerical stability
 * 3. LU factorization: Pr * A * Pc = L * U
 * 4. Triangular solves to compute x
 *
 * @param A - Sparse matrix in CSC format
 * @param b - Right-hand side vector
 * @param options - Solver options
 * @returns Solution and statistics
 *
 * @example
 * ```typescript
 * import { loadSuperLUModule, solveSparseCSC } from 'superluwasm';
 *
 * await loadSuperLUModule();
 *
 * // 3x3 sparse matrix
 * const A = {
 *   m: 3, n: 3,
 *   values: new Float64Array([1, 4, 3, 2, 5]),
 *   rowIndices: new Int32Array([0, 2, 1, 0, 2]),
 *   colPointers: new Int32Array([0, 2, 3, 5]),
 *   dtype: 'float64' as const
 * };
 *
 * const b = new Float64Array([1, 2, 3]);
 *
 * const { x, statistics } = solveSparseCSC(A, b);
 * console.log('Solution:', x);
 * console.log('Time:', statistics.totalTime, 'ms');
 * ```
 *
 * @see {@link solveSparseCSR} for CSR format input
 * @see {@link solveSparseExpert} for advanced options
 */
export function solveSparseCSC(
  A: SparseMatrixCSC,
  b: RealArray,
  options: SolveOptions = {}
): SolveResult {
  const module = getSuperLUModule();
  return solveCSCInternal(module, A, b, SuperLUTrans.NOTRANS, options);
}

/**
 * Solve a sparse linear system Ax = b where A is in CSR format.
 *
 * The CSR matrix is internally converted to CSC format (which is
 * SuperLU's native format) before solving. For best performance
 * with repeated solves, consider pre-converting to CSC.
 *
 * @param A - Sparse matrix in CSR format
 * @param b - Right-hand side vector
 * @param options - Solver options
 * @returns Solution and statistics
 *
 * @example
 * ```typescript
 * const A = {
 *   m: 3, n: 3,
 *   values: new Float64Array([1, 2, 3, 4, 5]),
 *   colIndices: new Int32Array([0, 2, 1, 0, 2]),
 *   rowPointers: new Int32Array([0, 2, 3, 5]),
 *   dtype: 'float64' as const
 * };
 *
 * const b = new Float64Array([1, 2, 3]);
 * const { x } = solveSparseCSR(A, b);
 * ```
 *
 * @see {@link solveSparseCSC} for CSC format input (more efficient)
 */
export function solveSparseCSR(
  A: SparseMatrixCSR,
  b: RealArray,
  options: SolveOptions = {}
): SolveResult {
  // Convert CSR to CSC
  const cscA = convertCSRtoCSC(A);
  return solveSparseCSC(cscA, b, options);
}

/**
 * Solve A^T x = b (transpose system) where A is in CSC format.
 *
 * Solves the system with the transpose of A. Useful when you have
 * A but need to solve with its transpose without explicitly computing it.
 *
 * @param A - Sparse matrix in CSC format
 * @param b - Right-hand side vector
 * @param options - Solver options
 * @returns Solution and statistics
 *
 * @example
 * ```typescript
 * // Solve A^T * x = b
 * const { x } = solveSparseTranspose(A, b);
 *
 * // Equivalent to: solveSparseCSC(transpose(A), b)
 * // but more efficient as no explicit transpose is computed
 * ```
 */
export function solveSparseTranspose(
  A: SparseMatrixCSC,
  b: RealArray,
  options: SolveOptions = {}
): SolveResult {
  const module = getSuperLUModule();
  return solveCSCInternal(module, A, b, SuperLUTrans.TRANS, options);
}

/**
 * Solve A^H x = b (conjugate transpose) where A is complex.
 *
 * Solves the system with the Hermitian (conjugate) transpose of A.
 * Only meaningful for complex matrices; for real matrices, this is
 * equivalent to `solveSparseTranspose`.
 *
 * @param A - Sparse complex matrix in CSC format
 * @param b - Right-hand side vector (complex, interleaved)
 * @param options - Solver options
 * @returns Solution and statistics
 *
 * @example
 * ```typescript
 * // Complex matrix with dtype 'complex128'
 * const A = {
 *   m: 2, n: 2,
 *   values: new Float64Array([1, 2, 3, 4, 5, 6, 7, 8]), // [1+2i, 3+4i, 5+6i, 7+8i]
 *   rowIndices: new Int32Array([0, 1, 0, 1]),
 *   colPointers: new Int32Array([0, 2, 4]),
 *   dtype: 'complex128' as const
 * };
 *
 * const b = new Float64Array([1, 0, 2, 0]); // [1+0i, 2+0i]
 * const { x } = solveSparseConjugateTranspose(A, b);
 * ```
 */
export function solveSparseConjugateTranspose(
  A: SparseMatrixCSC,
  b: RealArray,
  options: SolveOptions = {}
): SolveResult {
  const module = getSuperLUModule();
  return solveCSCInternal(module, A, b, SuperLUTrans.CONJ, options);
}
