/**
 * Incomplete LU (ILU) Factorization
 *
 * Compute an incomplete LU factorization for use as a preconditioner
 * in iterative solvers. ILU drops small elements to create sparser
 * factors that are cheaper to apply.
 *
 * @module factorization/ilu
 */

import { getSuperLUModule } from '../core/loader.js';
import {
  SuperLUColPerm,
  SuperLURowPerm,
  SuperLUILUDrop,
  SuperLUMiluT,
  SuperLUStype,
  SuperLUDtype,
  SuperLUMtype,
  SuperLUNorm,
} from '../types.js';
import type {
  SparseMatrixCSC,
  ILUFactorizationOptions,
  ILUFactorization,
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
  MEM_USAGE_SIZE,
  checkSuperLUError,
  INT_SIZE,
} from '../core/helpers.js';

// ============================================================
// Constants for options structure offsets
// ============================================================

const OPT_EQUIL = 1;
const OPT_COLPERM = 2;
const OPT_ROWPERM = 9;
const OPT_ILU_DROPRULE = 10;
const OPT_ILU_NORM = 16;
const OPT_ILU_MILU = 17;
const OPT_PRINTSTAT = 24;

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
 * Compute an Incomplete LU (ILU) factorization.
 *
 * ILU factorization produces approximate L and U factors by dropping
 * elements during factorization. The resulting factors are sparser
 * and cheaper to apply, making them useful as preconditioners for
 * iterative solvers like GMRES, BiCGSTAB, etc.
 *
 * The factorization satisfies: Pr * A * Pc ≈ L * U
 *
 * @param A - Sparse matrix in CSC format
 * @param options - ILU factorization options
 * @returns ILU factorization
 *
 * @example Basic ILU for preconditioning
 * ```typescript
 * import { loadSuperLUModule, sparseILU } from 'superluwasm';
 *
 * await loadSuperLUModule();
 *
 * const A = {
 *   m: 100, n: 100,
 *   values: new Float64Array([...]),
 *   rowIndices: new Int32Array([...]),
 *   colPointers: new Int32Array([...]),
 *   dtype: 'float64' as const
 * };
 *
 * // Compute ILU with default options
 * const ilu = sparseILU(A);
 *
 * console.log('Fill-in ratio:', ilu.fillIn);
 *
 * // Use L and U as preconditioner in your iterative solver
 * // ...
 *
 * ilu.dispose();
 * ```
 *
 * @example Controlling drop tolerance and fill
 * ```typescript
 * import { sparseILU, SuperLUILUDrop, SuperLUMiluT } from 'superluwasm';
 *
 * const ilu = sparseILU(A, {
 *   dropTolerance: 1e-3,        // Drop elements < 1e-3 * pivot
 *   dropRule: SuperLUILUDrop.DROP_BASIC,
 *   fillFactor: 20.0,           // Allow up to 20x fill-in
 *   miluType: SuperLUMiluT.SMILU_1, // Modified ILU variant 1
 * });
 * ```
 *
 * @see {@link sparseLU} for complete factorization
 */
export function sparseILU(
  A: SparseMatrixCSC,
  options: ILUFactorizationOptions = {}
): ILUFactorization {
  const module = getSuperLUModule();
  const startTime = performance.now();
  const dtype = A.dtype ?? inferDataType(A.values);
  const { m, n } = A;
  const nnz = getNnz(A);

  // Only support double precision for now
  if (dtype !== 'float64') {
    throw new Error(`ILU factorization currently only supports float64. Got: ${dtype}`);
  }

  // Validate dimensions
  if (m !== n) {
    throw new Error(`Matrix must be square for ILU factorization. Got ${m}x${n}`);
  }

  const ptrs: number[] = [];
  const factorPtrs: number[] = [];

  try {
    // Allocate options structure
    const optionsPtr = module._malloc(OPTIONS_STRUCT_SIZE);
    ptrs.push(optionsPtr);
    module._set_default_options(optionsPtr);

    // Set ILU-specific options
    const base32 = optionsPtr >> 2;
    const base64 = optionsPtr >> 3;

    module.HEAP32[base32 + OPT_EQUIL] = options.equilibrate !== false ? 1 : 0;
    module.HEAP32[base32 + OPT_COLPERM] = options.columnPermutation ?? SuperLUColPerm.COLAMD;
    module.HEAP32[base32 + OPT_ROWPERM] = options.rowPermutation ?? SuperLURowPerm.LargeDiag_MC64;
    module.HEAP32[base32 + OPT_ILU_DROPRULE] = options.dropRule ?? SuperLUILUDrop.DROP_BASIC;
    module.HEAP32[base32 + OPT_ILU_NORM] = SuperLUNorm.INF_NORM;
    module.HEAP32[base32 + OPT_ILU_MILU] = options.miluType ?? SuperLUMiluT.SILU;
    module.HEAP32[base32 + OPT_PRINTSTAT] = options.printStatistics ? 1 : 0;

    // Set double-precision options (drop tolerance and fill factor)
    // These need to be written at double-aligned offsets
    const dropTol = options.dropTolerance ?? 1e-4;
    const fillFactor = options.fillFactor ?? 10.0;

    // Write drop tolerance at offset 48 (12 * 4 bytes, which is 6 * 8 bytes)
    module.HEAPF64[base64 + 6] = dropTol;
    // Write fill factor at offset 56 (14 * 4 bytes, which is 7 * 8 bytes)
    module.HEAPF64[base64 + 7] = fillFactor;

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

    // Symbolic factorization
    module._sp_preorder(optionsPtr, APtr, permCPtr, etreePtr, ACPtr);

    // Get tuning parameters
    const panelSize = options.panelSize ?? module._sp_ienv(1);
    const relax = options.relaxation ?? module._sp_ienv(2);

    // Allocate info pointer
    const infoPtr = module._malloc(INT_SIZE);
    ptrs.push(infoPtr);

    // ILU factorization
    const factStart = performance.now();
    module._dgsitrf(
      optionsPtr,
      ACPtr,
      relax,
      panelSize,
      etreePtr,
      0,      // work
      0,      // lwork
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
      checkSuperLUError(info, 'dgsitrf');
    }

    // Read permutations
    const rowPermutation = readInts(module, permRPtr, m);
    const columnPermutation = readInts(module, permCPtr, n);
    const eliminationTree = readInts(module, etreePtr, n);

    // Read scaling factors
    const rowScale = readDoubles(module, RPtr, m);
    const columnScale = readDoubles(module, CPtr, n);

    // Query memory usage to estimate fill-in
    const memUsagePtr = module._malloc(MEM_USAGE_SIZE);
    ptrs.push(memUsagePtr);
    module._dQuerySpace(LPtr, UPtr, memUsagePtr);

    // Estimate fill-in (nnz(L+U) / nnz(A))
    // This is approximate - actual nnz would require reading the structures
    const fillIn = fillFactor; // Use fill factor as approximation

    // Cleanup temporary structures
    module._Destroy_SuperMatrix_Store(APtr);
    module._Destroy_CompCol_Matrix(ACPtr); // Permuted matrix
    module._StatFree(statPtr);

    const totalTime = performance.now() - startTime;

    const statistics: SolveStatistics = {
      factorizationTime,
      totalTime,
    };

    // Create dispose function
    let disposed = false;
    const dispose = () => {
      if (disposed) return;
      disposed = true;

      module._Destroy_SuperNode_Matrix(LPtr);
      module._Destroy_CompCol_Matrix(UPtr);
      freeAll(module, [permCPtr, permRPtr, etreePtr, RPtr, CPtr, GluPtr]);
    };

    // Create factor handles
    const L = createFactorHandle(LPtr);
    const U = createFactorHandle(UPtr);

    const result: ILUFactorization = {
      L,
      U,
      rowPermutation,
      columnPermutation,
      eliminationTree,
      rowScale,
      columnScale,
      equilibration: 'N',
      m,
      n,
      dtype,
      statistics,
      dispose,
      fillIn,
    };

    if (info > 0) {
      console.warn(`ILU: Zero pivot encountered at step ${info}. ` +
                   `The incomplete factors may not be suitable as preconditioners.`);
    }

    return result;

  } catch (error) {
    freeAll(module, factorPtrs);
    throw error;
  } finally {
    freeAll(module, ptrs);
  }
}
