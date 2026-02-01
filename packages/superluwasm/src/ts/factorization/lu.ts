/**
 * Sparse LU Factorization
 *
 * Compute the LU factorization of a sparse matrix: Pr * A * Pc = L * U
 * where L is lower triangular, U is upper triangular, and Pr, Pc are
 * permutation matrices.
 *
 * @module factorization/lu
 */

import { getSuperLUModule } from '../core/loader.js';
import {
  SuperLUColPerm,
  SuperLURowPerm,
  SuperLUStype,
  SuperLUDtype,
  SuperLUMtype,
} from '../types.js';
import type {
  SparseMatrixCSC,
  LUFactorizationOptions,
  LUFactorization,
  LUFactorHandle,
  SolveStatistics,
  RealArray,
} from '../high-level-types.js';
import { inferDataType, getNnz } from '../high-level-types.js';
import {
  allocateDoubles,
  allocateInts,
  readDoubles,
  readInts,
  readInt,
  freeAll,
  toFloat64Array,
  toInt32Array,
  OPTIONS_STRUCT_SIZE,
  STAT_STRUCT_SIZE,
  SUPERMATRIX_SIZE,
  GLOBALLU_SIZE,
  checkSuperLUError,
  INT_SIZE,
} from '../core/helpers.js';

// ============================================================
// Internal Helpers
// ============================================================

/**
 * Create an LU factor handle.
 */
function createFactorHandle(ptr: number): LUFactorHandle {
  return {
    _ptr: ptr,
    _disposed: false,
  };
}

// ============================================================
// Public API
// ============================================================

/**
 * Compute the LU factorization of a sparse matrix.
 *
 * Computes the factorization Pr * A * Pc = L * U where:
 * - L is unit lower triangular (stored in SuperLU's supernode format)
 * - U is upper triangular (stored in compressed column format)
 * - Pr is the row permutation matrix
 * - Pc is the column permutation matrix
 *
 * The factorization can be reused to solve multiple systems with the
 * same matrix but different right-hand sides.
 *
 * @param A - Sparse matrix in CSC format
 * @param options - Factorization options
 * @returns LU factorization that can be reused for solving
 *
 * @example Basic factorization
 * ```typescript
 * import { loadSuperLUModule, sparseLU } from 'superluwasm';
 *
 * await loadSuperLUModule();
 *
 * const A = {
 *   m: 3, n: 3,
 *   values: new Float64Array([4, 1, 1, 3, 2, 5]),
 *   rowIndices: new Int32Array([0, 1, 1, 2, 0, 2]),
 *   colPointers: new Int32Array([0, 2, 4, 6]),
 *   dtype: 'float64' as const
 * };
 *
 * const lu = sparseLU(A);
 *
 * console.log('Row permutation:', lu.rowPermutation);
 * console.log('Column permutation:', lu.columnPermutation);
 * console.log('Factorization time:', lu.statistics.totalTime, 'ms');
 *
 * // Use for solving (with solveSparseExpert)
 * // ...
 *
 * // Clean up when done
 * lu.dispose();
 * ```
 *
 * @example With options
 * ```typescript
 * const lu = sparseLU(A, {
 *   columnPermutation: SuperLUColPerm.METIS_AT_PLUS_A, // Better for large matrices
 *   rowPermutation: SuperLURowPerm.LargeDiag_MC64,
 *   equilibrate: true,
 *   printStatistics: true
 * });
 * ```
 *
 * @see {@link sparseILU} for incomplete factorization (preconditioning)
 * @see {@link solveSparseExpert} to use factorization for solving
 */
export function sparseLU(
  A: SparseMatrixCSC,
  options: LUFactorizationOptions = {}
): LUFactorization {
  const module = getSuperLUModule();
  const startTime = performance.now();
  const dtype = A.dtype ?? inferDataType(A.values);
  const { m, n } = A;
  const nnz = getNnz(A);

  // Only support double precision for now
  if (dtype !== 'float64') {
    throw new Error(`LU factorization currently only supports float64. Got: ${dtype}`);
  }

  // Validate dimensions
  if (m !== n) {
    throw new Error(`Matrix must be square for LU factorization. Got ${m}x${n}`);
  }

  const ptrs: number[] = [];
  // Track pointers that should NOT be freed (they're part of the factorization)
  const factorPtrs: number[] = [];

  try {
    // Allocate options structure
    const optionsPtr = module._malloc(OPTIONS_STRUCT_SIZE);
    ptrs.push(optionsPtr);
    module._set_default_options(optionsPtr);

    // Set options
    const base = optionsPtr >> 2;
    module.HEAP32[base + 1] = options.equilibrate !== false ? 1 : 0; // Equil
    module.HEAP32[base + 2] = options.columnPermutation ?? SuperLUColPerm.COLAMD; // ColPerm
    module.HEAP32[base + 9] = options.rowPermutation ?? SuperLURowPerm.LargeDiag_MC64; // RowPerm
    module.HEAP32[base + 24] = options.printStatistics ? 1 : 0; // PrintStat

    // Initialize statistics
    const statPtr = module._malloc(STAT_STRUCT_SIZE);
    ptrs.push(statPtr);
    module._StatInit(statPtr);

    // Allocate SuperMatrix structures
    const APtr = module._malloc(SUPERMATRIX_SIZE);
    const ACPtr = module._malloc(SUPERMATRIX_SIZE);
    const LPtr = module._malloc(SUPERMATRIX_SIZE);
    const UPtr = module._malloc(SUPERMATRIX_SIZE);
    ptrs.push(APtr, ACPtr);
    factorPtrs.push(LPtr, UPtr);

    // Allocate GlobalLU structure
    const GluPtr = module._malloc(GLOBALLU_SIZE);
    factorPtrs.push(GluPtr);

    // Allocate and copy matrix data
    const valuesPtr = allocateDoubles(module, toFloat64Array(A.values as RealArray));
    const rowindPtr = allocateInts(module, toInt32Array(A.rowIndices));
    const colptrPtr = allocateInts(module, toInt32Array(A.colPointers));
    ptrs.push(valuesPtr, rowindPtr, colptrPtr);

    // Create CSC matrix
    module._dCreate_CompCol_Matrix(
      APtr, m, n, nnz, valuesPtr, rowindPtr, colptrPtr,
      SuperLUStype.SLU_NC, SuperLUDtype.SLU_D, SuperLUMtype.SLU_GE
    );

    // Allocate permutation arrays
    const permCPtr = allocateInts(module, null, n);
    const permRPtr = allocateInts(module, null, m);
    factorPtrs.push(permCPtr, permRPtr);

    // Allocate elimination tree
    const etreePtr = allocateInts(module, null, n);
    factorPtrs.push(etreePtr);

    // Allocate scaling factors
    const RPtr = allocateDoubles(module, null, m);
    const CPtr = allocateDoubles(module, null, n);
    factorPtrs.push(RPtr, CPtr);

    // Get column permutation
    const colPerm = options.columnPermutation ?? SuperLUColPerm.COLAMD;
    module._get_perm_c(colPerm, APtr, permCPtr);

    // Symbolic factorization (preordering)
    module._sp_preorder(optionsPtr, APtr, permCPtr, etreePtr, ACPtr);

    // Get tuning parameters
    const panelSize = options.panelSize ?? module._sp_ienv(1);
    const relax = options.relaxation ?? module._sp_ienv(2);

    // Allocate info pointer
    const infoPtr = module._malloc(INT_SIZE);
    ptrs.push(infoPtr);

    // Numerical factorization
    const factStart = performance.now();
    module._dgstrf(
      optionsPtr,
      ACPtr,
      relax,
      panelSize,
      etreePtr,
      0,      // work
      0,      // lwork (0 = automatic allocation)
      permCPtr,
      permRPtr,
      LPtr,
      UPtr,
      GluPtr,
      statPtr,
      infoPtr
    );
    const factorizationTime = performance.now() - factStart;

    // Check for errors
    const info = readInt(module, infoPtr);
    if (info < 0) {
      checkSuperLUError(info, 'dgstrf');
    }

    // Read permutations
    const rowPermutation = readInts(module, permRPtr, m);
    const columnPermutation = readInts(module, permCPtr, n);
    const eliminationTree = readInts(module, etreePtr, n);

    // Read scaling factors
    const rowScale = readDoubles(module, RPtr, m);
    const columnScale = readDoubles(module, CPtr, n);

    // Cleanup temporary structures (but not factor structures)
    module._Destroy_SuperMatrix_Store(APtr);
    module._Destroy_CompCol_Matrix(ACPtr); // Permuted matrix
    module._StatFree(statPtr);

    const totalTime = performance.now() - startTime;

    const statistics: SolveStatistics = {
      factorizationTime,
      totalTime,
    };

    // Create dispose function that cleans up factor memory
    let disposed = false;
    const dispose = () => {
      if (disposed) return;
      disposed = true;

      // Free L and U matrices
      module._Destroy_SuperNode_Matrix(LPtr);
      module._Destroy_CompCol_Matrix(UPtr);

      // Free other factor-related memory
      freeAll(module, [permCPtr, permRPtr, etreePtr, RPtr, CPtr, GluPtr]);
    };

    // Create factor handles
    const L = createFactorHandle(LPtr);
    const U = createFactorHandle(UPtr);

    const result: LUFactorization = {
      L,
      U,
      rowPermutation,
      columnPermutation,
      eliminationTree,
      rowScale,
      columnScale,
      equilibration: 'N', // TODO: Read actual equilibration status
      m,
      n,
      dtype,
      statistics,
      dispose,
    };

    // Warn if matrix was singular
    if (info > 0) {
      console.warn(`Matrix may be singular: U(${info},${info}) is exactly zero`);
    }

    return result;

  } catch (error) {
    // On error, free factor pointers too
    freeAll(module, factorPtrs);
    throw error;
  } finally {
    // Always free temporary pointers
    freeAll(module, ptrs);
  }
}
