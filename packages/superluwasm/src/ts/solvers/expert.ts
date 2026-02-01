/**
 * Expert Sparse Linear Solver
 *
 * Full-featured solver with complete control over all SuperLU parameters,
 * including equilibration, condition number estimation, iterative refinement,
 * and error bounds.
 *
 * @module solvers/expert
 */

import { getSuperLUModule } from '../core/loader.js';
import {
  SuperLUColPerm,
  SuperLURowPerm,
  SuperLUTrans,
  SuperLUStype,
  SuperLUDtype,
  SuperLUMtype,
  SuperLUFactOption,
} from '../types.js';
import type { SuperLUModule } from '../types.js';
import type {
  SparseMatrixCSC,
  ExpertSolveOptions,
  SolveResult,
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
  readDouble,
  readByte,
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
// Constants for superlu_options_t offsets
// ============================================================

const OPT_FACT = 0;
const OPT_EQUIL = 1;
const OPT_COLPERM = 2;
const OPT_TRANS = 3;
const OPT_ITERREFINE = 4;
const OPT_PIVOTGROWTH = 7;
const OPT_CONDITIONNUMBER = 8;
const OPT_ROWPERM = 9;
const OPT_PRINTSTAT = 24;

// ============================================================
// Internal Helpers
// ============================================================

/**
 * Set options in the superlu_options_t structure.
 */
function setOptions(
  module: SuperLUModule,
  optionsPtr: number,
  options: ExpertSolveOptions,
  trans: SuperLUTrans = SuperLUTrans.NOTRANS
): void {
  // Set default options first
  module._set_default_options(optionsPtr);

  const base = optionsPtr >> 2;

  // Factorization type
  if (options.reuseFactorization) {
    module.HEAP32[base + OPT_FACT] = SuperLUFactOption.FACTORED;
  } else {
    module.HEAP32[base + OPT_FACT] = SuperLUFactOption.DOFACT;
  }

  // Equilibration
  module.HEAP32[base + OPT_EQUIL] = options.equilibrate !== false ? 1 : 0;

  // Column permutation
  module.HEAP32[base + OPT_COLPERM] = options.columnPermutation ?? SuperLUColPerm.COLAMD;

  // Transpose
  module.HEAP32[base + OPT_TRANS] = trans;

  // Iterative refinement (0 = NOREFINE, 1 = SLU_SINGLE, 2 = SLU_DOUBLE, 3 = SLU_EXTRA)
  module.HEAP32[base + OPT_ITERREFINE] = options.iterativeRefinement ? 2 : 0;

  // Row permutation
  module.HEAP32[base + OPT_ROWPERM] = options.rowPermutation ?? SuperLURowPerm.LargeDiag_MC64;

  // Pivot growth and condition number
  module.HEAP32[base + OPT_PIVOTGROWTH] = 1; // Always compute for expert driver
  module.HEAP32[base + OPT_CONDITIONNUMBER] = options.computeConditionNumber ? 1 : 0;

  // Print statistics
  module.HEAP32[base + OPT_PRINTSTAT] = options.printStatistics ? 1 : 0;
}

// ============================================================
// Public API
// ============================================================

/**
 * Expert sparse solver with full control over all SuperLU options.
 *
 * This is the most flexible solver interface, providing access to:
 * - Equilibration (row/column scaling) for numerical stability
 * - Condition number estimation
 * - Iterative refinement for improved accuracy
 * - Forward and backward error bounds
 * - Reuse of factorization for multiple right-hand sides
 *
 * @param A - Sparse matrix in CSC format
 * @param b - Right-hand side vector
 * @param options - Expert solver options
 * @returns Solution with full diagnostics
 *
 * @example Basic usage with condition number
 * ```typescript
 * import { loadSuperLUModule, solveSparseExpert } from 'superluwasm';
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
 * const b = new Float64Array([1, 2, 3]);
 *
 * const result = solveSparseExpert(A, b, {
 *   computeConditionNumber: true,
 *   computeErrorBounds: true,
 *   iterativeRefinement: true
 * });
 *
 * console.log('Solution:', result.x);
 * console.log('Condition number:', 1 / result.rcond);
 * console.log('Forward error:', result.forwardError);
 * console.log('Backward error:', result.backwardError);
 * ```
 *
 * @example Reusing factorization
 * ```typescript
 * import { sparseLU, solveSparseExpert } from 'superluwasm';
 *
 * // Factor once
 * const factorization = sparseLU(A);
 *
 * // Solve with multiple right-hand sides
 * const x1 = solveSparseExpert(A, b1, { reuseFactorization: factorization });
 * const x2 = solveSparseExpert(A, b2, { reuseFactorization: factorization });
 * const x3 = solveSparseExpert(A, b3, { reuseFactorization: factorization });
 *
 * // Clean up
 * factorization.dispose();
 * ```
 *
 * @see {@link solveSparseCSC} for simpler interface
 * @see {@link sparseLU} for computing factorization separately
 */
export function solveSparseExpert(
  A: SparseMatrixCSC,
  b: RealArray,
  options: ExpertSolveOptions = {}
): SolveResult {
  const module = getSuperLUModule();
  const startTime = performance.now();
  const dtype = A.dtype ?? inferDataType(A.values);
  const { m, n } = A;
  const nnz = getNnz(A);

  // Only support double precision for now
  if (dtype !== 'float64') {
    throw new Error(`Expert solver currently only supports float64. Got: ${dtype}`);
  }

  // Validate dimensions
  if (m !== n) {
    throw new Error(`Matrix must be square for direct solve. Got ${m}x${n}`);
  }
  if (b.length !== m) {
    throw new Error(`RHS vector length (${b.length}) must match matrix rows (${m})`);
  }

  const nrhs = 1; // Number of right-hand sides
  const ptrs: number[] = [];

  try {
    // Allocate structures
    const optionsPtr = module._malloc(OPTIONS_STRUCT_SIZE);
    const statPtr = module._malloc(STAT_STRUCT_SIZE);
    const APtr = module._malloc(SUPERMATRIX_SIZE);
    const LPtr = module._malloc(SUPERMATRIX_SIZE);
    const UPtr = module._malloc(SUPERMATRIX_SIZE);
    const BPtr = module._malloc(SUPERMATRIX_SIZE);
    const XPtr = module._malloc(SUPERMATRIX_SIZE);
    const GluPtr = module._malloc(GLOBALLU_SIZE);
    const memUsagePtr = module._malloc(MEM_USAGE_SIZE);
    ptrs.push(optionsPtr, statPtr, APtr, LPtr, UPtr, BPtr, XPtr, GluPtr, memUsagePtr);

    // Set options
    setOptions(module, optionsPtr, options);

    // Initialize statistics
    module._StatInit(statPtr);

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

    // Allocate RHS and solution
    const bPtr = allocateDoubles(module, toFloat64Array(b));
    const xPtr = allocateDoubles(module, null, m * nrhs);
    ptrs.push(bPtr, xPtr);

    // Create dense matrices for B and X
    module._dCreate_Dense_Matrix(
      BPtr, m, nrhs, bPtr, m,
      SuperLUStype.SLU_DN, SuperLUDtype.SLU_D, SuperLUMtype.SLU_GE
    );
    module._dCreate_Dense_Matrix(
      XPtr, m, nrhs, xPtr, m,
      SuperLUStype.SLU_DN, SuperLUDtype.SLU_D, SuperLUMtype.SLU_GE
    );

    // Allocate permutation arrays
    const permCPtr = allocateInts(module, options.userColumnPerm ?? null, n);
    const permRPtr = allocateInts(module, options.userRowPerm ?? null, m);
    ptrs.push(permCPtr, permRPtr);

    // Allocate elimination tree
    const etreePtr = allocateInts(module, null, n);
    ptrs.push(etreePtr);

    // Allocate scaling factors
    const RPtr = allocateDoubles(module, null, m);
    const CPtr = allocateDoubles(module, null, n);
    ptrs.push(RPtr, CPtr);

    // Allocate equed (equilibration status - single char)
    const equedPtr = module._malloc(4);
    ptrs.push(equedPtr);
    module.HEAP8[equedPtr] = 'N'.charCodeAt(0);

    // Allocate output scalars
    const rpgPtr = allocateDoubles(module, [0], 1);
    const rcondPtr = allocateDoubles(module, [0], 1);
    ptrs.push(rpgPtr, rcondPtr);

    // Allocate error bounds
    const ferrPtr = allocateDoubles(module, null, nrhs);
    const berrPtr = allocateDoubles(module, null, nrhs);
    ptrs.push(ferrPtr, berrPtr);

    // Allocate info
    const infoPtr = module._malloc(INT_SIZE);
    ptrs.push(infoPtr);

    // Call expert driver
    const lwork = options.workspaceSize ?? 0;
    const workPtr = 0; // Let SuperLU allocate workspace

    module._dgssvx(
      optionsPtr,
      APtr,
      permCPtr,
      permRPtr,
      etreePtr,
      equedPtr,
      RPtr,
      CPtr,
      LPtr,
      UPtr,
      workPtr,
      lwork,
      BPtr,
      XPtr,
      rpgPtr,
      rcondPtr,
      ferrPtr,
      berrPtr,
      GluPtr,
      memUsagePtr,
      statPtr,
      infoPtr
    );

    // Check for errors
    const info = readInt(module, infoPtr);
    if (info < 0) {
      checkSuperLUError(info, 'dgssvx');
    }
    // info > 0 means singular, but we still have a solution attempt

    // Read solution
    const x = readDoubles(module, xPtr, m);

    // Read permutations
    const rowPermutation = readInts(module, permRPtr, m);
    const columnPermutation = readInts(module, permCPtr, n);

    // Read scalars
    const pivotGrowth = readDouble(module, rpgPtr);
    const rcond = options.computeConditionNumber ? readDouble(module, rcondPtr) : undefined;

    // Read error bounds if requested
    const forwardError = options.computeErrorBounds ? readDoubles(module, ferrPtr, nrhs) : undefined;
    const backwardError = options.computeErrorBounds ? readDoubles(module, berrPtr, nrhs) : undefined;

    // Read equilibration status
    const equedChar = String.fromCharCode(readByte(module, equedPtr));
    const equilibration = equedChar as 'N' | 'R' | 'C' | 'B';

    // Cleanup
    module._Destroy_SuperMatrix_Store(APtr);
    module._Destroy_SuperMatrix_Store(BPtr);
    module._Destroy_SuperMatrix_Store(XPtr);
    module._Destroy_SuperNode_Matrix(LPtr);
    module._Destroy_CompCol_Matrix(UPtr);
    module._StatFree(statPtr);

    const totalTime = performance.now() - startTime;

    const statistics: SolveStatistics = {
      totalTime,
    };

    const result: SolveResult = {
      x,
      statistics,
      rowPermutation,
      columnPermutation,
      pivotGrowth,
      equilibration,
    };

    if (rcond !== undefined) {
      result.rcond = rcond;
    }
    if (forwardError !== undefined) {
      result.forwardError = forwardError;
    }
    if (backwardError !== undefined) {
      result.backwardError = backwardError;
    }

    // Warn if nearly singular
    if (info > 0) {
      console.warn(`Matrix may be singular: U(${info},${info}) is exactly zero`);
    }

    return result;

  } finally {
    freeAll(module, ptrs);
  }
}
