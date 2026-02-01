# SuperLU WASM TypeScript API TODO

Complete high-level TypeScript API for sparse direct solver functionality.

---

## 1. Module Loading & Configuration

| Function | Status | Description | File |
|----------|--------|-------------|------|
| `loadSuperLUModule()` | ✅ | Async WASM module initialization | `ts/core/loader.ts` |
| `getSuperLUModule()` | ✅ | Get loaded module synchronously | `ts/core/loader.ts` |
| `isSuperLULoaded()` | ✅ | Check module status | `ts/core/loader.ts` |
| `isSuperLULoading()` | ✅ | Check if module is loading | `ts/core/loader.ts` |
| `resetSuperLUModule()` | ✅ | Unload and reset | `ts/core/loader.ts` |
| `configureSuperLU(config)` | ✅ | Configure WASM asset URLs | `ts/core/loader.ts` |
| `getSuperLUConfig()` | ✅ | Get current configuration | `ts/core/loader.ts` |

---

## 2. Core Sparse Linear Solvers

| Function | Status | Description | File |
|----------|--------|-------------|------|
| `solveSparseCSC(A, b, options?)` | ✅ | Solve Ax=b with CSC matrix | `ts/solvers/direct.ts` |
| `solveSparseCSR(A, b, options?)` | ✅ | Solve Ax=b with CSR matrix | `ts/solvers/direct.ts` |
| `solveSparseTranspose(A, b, options?)` | ✅ | Solve A^T x = b | `ts/solvers/direct.ts` |
| `solveSparseConjugateTranspose(A, b, options?)` | ✅ | Solve A^H x = b (complex) | `ts/solvers/direct.ts` |
| `solveSparseExpert(A, b, options?)` | ✅ | Full control over all SuperLU options | `ts/solvers/expert.ts` |

---

## 3. Sparse LU Factorization

| Function | Status | Description | File |
|----------|--------|-------------|------|
| `sparseLU(A, options?)` | ✅ | Compute L, U, P, Q factorization | `ts/factorization/lu.ts` |
| `sparseILU(A, options?)` | ✅ | Incomplete LU for preconditioning | `ts/factorization/ilu.ts` |

### ILU Options
- `dropTolerance` - Drop strategy threshold
- `dropRule` - DROP_BASIC, DROP_PROWS, DROP_COLUMN, DROP_AREA, DROP_DYNAMIC
- `miluType` - SILU, SMILU_1, SMILU_2, SMILU_3
- `fillFactor` - Maximum allowed fill-in ratio

---

## 4. Triangular Solvers (Using LU Factors)

| Function | Status | Description | File |
|----------|--------|-------------|------|
| `solveTriangularLU(L, U, permC, permR, b, options?)` | 🔲 | Solve with precomputed LU | `ts/solvers/triangular.ts` |
| `solveTriangularLUTranspose(L, U, permC, permR, b)` | 🔲 | Transpose solve | `ts/solvers/triangular.ts` |
| `solveTriangularLUConjugateTranspose(L, U, permC, permR, b)` | 🔲 | Hermitian transpose | `ts/solvers/triangular.ts` |
| `solveTriangularCSC(A, b, lower?, unitDiagonal?, transpose?)` | 🔲 | Direct triangular CSC | `ts/solvers/triangular.ts` |
| `solveTriangularCSR(A, b, lower?, unitDiagonal?, transpose?)` | 🔲 | Direct triangular CSR | `ts/solvers/triangular.ts` |

---

## 5. Multiple Right-Hand Sides

| Function | Status | Description | File |
|----------|--------|-------------|------|
| `solveMultipleRHS(A, B, options?)` | 🔲 | Solve AX=B (B is matrix) | `ts/solvers/batch.ts` |
| `solveReusedFactorization(L, U, permC, permR, B)` | 🔲 | Reuse LU for multiple RHS | `ts/solvers/batch.ts` |

---

## 6. Condition Number & Stability Analysis

| Function | Status | Description | File |
|----------|--------|-------------|------|
| `estimateConditionNumber(L, U, normType?, originalMatrixNorm?)` | 🔲 | Reciprocal condition number | `ts/analysis/condition.ts` |
| `computeMatrixNorm(A, normType)` | 🔲 | ONE_NORM, INF_NORM, TWO_NORM, FROB_NORM | `ts/analysis/condition.ts` |
| `iterativelyRefine(A, L, U, permC, permR, b, x, options?)` | 🔲 | Improve solution accuracy | `ts/analysis/refinement.ts` |

---

## 7. Equilibration (Scaling)

| Function | Status | Description | File |
|----------|--------|-------------|------|
| `computeEquilibrationScaling(A)` | 🔲 | Get row/column scaling factors | `ts/equilibration/scaling.ts` |
| `applyEquilibrationScaling(A, rowScaling, columnScaling)` | 🔲 | Scale matrix | `ts/equilibration/scaling.ts` |
| `undoEquilibrationScaling(x, columnScaling)` | 🔲 | Unscale solution | `ts/equilibration/scaling.ts` |

---

## 8. Permutation Utilities

| Function | Status | Description | File |
|----------|--------|-------------|------|
| `computeColumnPermutation(A, strategy)` | 🔲 | NATURAL, MMD_ATA, MMD_AT_PLUS_A, COLAMD, METIS | `ts/permutation/ordering.ts` |
| `computeRowPermutation(A, strategy)` | 🔲 | NOROWPERM, LargeDiag_MC64, LargeDiag_AWPM | `ts/permutation/ordering.ts` |
| `permuteSparseMatrix(A, rowPerm, colPerm, format?)` | 🔲 | Apply P*A*Q | `ts/permutation/ordering.ts` |

---

## 9. Elimination Tree & Symbolic Analysis

| Function | Status | Description | File |
|----------|--------|-------------|------|
| `getEliminationTree(A, columnPermutation?)` | 🔲 | Factorization structure | `ts/analysis/symbolic.ts` |
| `getEliminationTreeStatistics(etree, A, columnPermutation?)` | 🔲 | Predict fill-in | `ts/analysis/symbolic.ts` |
| `predictFillIn(A, perm)` | 🔲 | Estimate fill before factorization | `ts/analysis/symbolic.ts` |
| `getFlopsEstimate(A)` | 🔲 | Estimate floating point operations | `ts/analysis/symbolic.ts` |

---

## 10. Sparse Matrix Construction

| Function | Status | Description | File |
|----------|--------|-------------|------|
| `createSparseCSC(m, n, rowIndices, colPointers, values, dtype?)` | 🔲 | Create CSC matrix | `ts/matrix/construction.ts` |
| `createSparseCSR(m, n, colIndices, rowPointers, values, dtype?)` | 🔲 | Create CSR matrix | `ts/matrix/construction.ts` |
| `createSparseCOO(m, n, rows, cols, values, dtype?)` | 🔲 | Create COO matrix | `ts/matrix/construction.ts` |
| `createDenseMatrix(data, m, n, rowMajor?)` | 🔲 | Create dense matrix | `ts/matrix/construction.ts` |

---

## 11. Format Conversion

| Function | Status | Description | File |
|----------|--------|-------------|------|
| `convertCSRtoCSC(A)` | 🔲 | CSR → CSC | `ts/matrix/conversion.ts` |
| `convertCSCtoCSR(A)` | 🔲 | CSC → CSR | `ts/matrix/conversion.ts` |
| `toCOO(matrix)` | 🔲 | Any format → COO | `ts/matrix/conversion.ts` |
| `toCSC(matrix)` | 🔲 | Any format → CSC | `ts/matrix/conversion.ts` |
| `toCSR(matrix)` | 🔲 | Any format → CSR | `ts/matrix/conversion.ts` |
| `getTransposeCSC(A)` | 🔲 | Efficient transpose | `ts/matrix/conversion.ts` |
| `sparseToDense(sparseMatrix)` | 🔲 | Sparse → dense array | `ts/matrix/conversion.ts` |
| `denseToSparse(denseArray, format)` | 🔲 | Dense → sparse | `ts/matrix/conversion.ts` |

---

## 12. Matrix Utilities

| Function | Status | Description | File |
|----------|--------|-------------|------|
| `copySparseMatrix(A, format?)` | 🔲 | Deep copy | `ts/matrix/utilities.ts` |
| `getMatrixNonzeros(A)` | 🔲 | Count/list nonzero elements | `ts/matrix/utilities.ts` |
| `getMatrixStatistics(A)` | 🔲 | Sparsity %, fill pattern analysis | `ts/matrix/utilities.ts` |
| `getSparsityPattern(A)` | 🔲 | Visualization data | `ts/matrix/utilities.ts` |
| `getSupernodeStructure(A, perm)` | 🔲 | Supernode organization | `ts/matrix/utilities.ts` |

---

## 13. Matrix Properties & Validation

| Function | Status | Description | File |
|----------|--------|-------------|------|
| `validateSparseMatrix(A)` | 🔲 | Check format consistency | `ts/matrix/validation.ts` |
| `isSymmetricSparse(A, tolerance?)` | 🔲 | Check symmetry | `ts/matrix/validation.ts` |
| `isHermitianSparse(A, tolerance?)` | 🔲 | Check Hermitian property | `ts/matrix/validation.ts` |
| `isPositiveDefiniteSparse(A)` | 🔲 | Check positive definiteness | `ts/matrix/validation.ts` |
| `detectMatrixProperties(A)` | 🔲 | Full property analysis | `ts/matrix/validation.ts` |
| `checkFactorization(A, L, U, P, Q, tolerance?)` | 🔲 | Verify ‖PAQ - LU‖ | `ts/matrix/validation.ts` |

---

## 14. Specialized Solvers

| Function | Status | Description | File |
|----------|--------|-------------|------|
| `solveSymmetricSparse(A, b, options?)` | 🔲 | Symmetric solver (uses only triangle) | `ts/solvers/symmetric.ts` |
| `solveHermitianSparse(A, b, options?)` | 🔲 | Complex Hermitian | `ts/solvers/symmetric.ts` |
| `solveLeastSquaresSparse(A, b, options?)` | 🔲 | Overdetermined systems | `ts/solvers/least-squares.ts` |
| `solveWithIterativeRefinement(A, L, U, permC, permR, b, maxIter?, tol?)` | 🔲 | Refined solve | `ts/solvers/refined.ts` |
| `solveShiftedSystem(A, shift, b, options?)` | 🔲 | (A - σI)x = b for eigenvalue work | `ts/solvers/shifted.ts` |
| `solveWithDeflation(A, b, deflationVectors, options?)` | 🔲 | Deflated solve | `ts/solvers/deflation.ts` |

---

## 15. Memory Management & Statistics

| Function | Status | Description | File |
|----------|--------|-------------|------|
| `getMemoryUsage(L, U)` | 🔲 | L_bytes, U_bytes, total | `ts/core/memory.ts` |
| `estimateMemoryRequired(A, options?)` | 🔲 | Predict before factorization | `ts/core/memory.ts` |
| `queryFactorStatistics(L, U)` | 🔲 | nnz_L, nnz_U, flops | `ts/core/memory.ts` |
| `freeFactorization(L, U)` | 🔲 | Release LU factor memory | `ts/core/memory.ts` |
| `allocateWorkspace(size)` | 🔲 | User-supplied workspace | `ts/core/memory.ts` |
| `deallocateWorkspace(ptr)` | 🔲 | Release workspace | `ts/core/memory.ts` |
| `setWorkspaceSize(bytes)` | 🔲 | Global workspace limit | `ts/core/memory.ts` |

---

## 16. Error Handling & Diagnostics

| Function | Status | Description | File |
|----------|--------|-------------|------|
| `getSuperLUErrorMessage(errorCode)` | 🔲 | Human-readable errors | `ts/core/errors.ts` |
| `getLastSolveStatistics()` | 🔲 | factorTime, solveTime, refinementTime | `ts/core/errors.ts` |
| `enableDetailedLogging(level?)` | 🔲 | 0=off, 1=basic, 2=detailed, 3=verbose | `ts/core/errors.ts` |

---

## 17. Test Matrix Generation

| Function | Status | Description | File |
|----------|--------|-------------|------|
| `createTestMatrix(type, n)` | 🔲 | Standard test matrices | `ts/testing/matrices.ts` |

### Test Matrix Types
- `'tridiagonal'` - Tridiagonal matrix
- `'laplacian_2d'` - 2D Laplacian (5-point stencil)
- `'random'` - Random dense matrix
- `'sparse_random'` - Random sparse matrix
- `'poisson'` - Poisson problem matrix

---

## 18. Types & Interfaces

| Type | Status | Description | File |
|------|--------|-------------|------|
| `SparseMatrix` | 🔲 | Sparse matrix type (CSC/CSR/COO) | `ts/types.ts` |
| `LUFactorization` | 🔲 | L, U, P, Q result | `ts/types.ts` |
| `SolveResult` | 🔲 | Solution + statistics | `ts/types.ts` |
| `SolveOptions` | 🔲 | Solver configuration | `ts/types.ts` |
| `ILUOptions` | 🔲 | ILU configuration | `ts/types.ts` |
| `FactorStatistics` | 🔲 | Factorization timing/memory | `ts/types.ts` |
| `SolveStatistics` | 🔲 | Solve timing/iterations | `ts/types.ts` |
| `PermutationStrategy` | 🔲 | Column/row ordering enum | `ts/types.ts` |
| `DropStrategy` | 🔲 | ILU drop rule enum | `ts/types.ts` |
| `MILUType` | 🔲 | Modified ILU type enum | `ts/types.ts` |
| `NormType` | 🔲 | Matrix norm type enum | `ts/types.ts` |

---

## Precision Support

All solver functions support:
- `float32` - Single precision real
- `float64` - Double precision real
- `complex64` - Single precision complex
- `complex128` - Double precision complex

---

## Folder Structure

```
src/ts/
├── index.ts                    # Main exports
├── types.ts                    # Type definitions
├── core/
│   ├── loader.ts               # Module loading
│   ├── memory.ts               # Memory management
│   └── errors.ts               # Error handling
├── solvers/
│   ├── direct.ts               # Direct solvers (CSC/CSR)
│   ├── expert.ts               # Expert solver with full options
│   ├── triangular.ts           # Triangular solvers
│   ├── batch.ts                # Multiple RHS
│   ├── symmetric.ts            # Symmetric/Hermitian solvers
│   ├── least-squares.ts        # Least squares
│   ├── refined.ts              # Iterative refinement
│   ├── shifted.ts              # Shift-invert
│   └── deflation.ts            # Deflated solvers
├── factorization/
│   ├── lu.ts                   # Sparse LU
│   └── ilu.ts                  # Incomplete LU
├── matrix/
│   ├── construction.ts         # Matrix creation
│   ├── conversion.ts           # Format conversion
│   ├── utilities.ts            # Matrix utilities
│   └── validation.ts           # Property checking
├── analysis/
│   ├── condition.ts            # Condition number
│   ├── refinement.ts           # Iterative refinement
│   └── symbolic.ts             # Elimination tree
├── equilibration/
│   └── scaling.ts              # Row/column scaling
├── permutation/
│   └── ordering.ts             # Permutation strategies
└── testing/
    └── matrices.ts             # Test matrix generation
```

---

## Summary

| Category | Count | Status |
|----------|-------|--------|
| Module Loading | 7 | ✅ |
| Core Solvers | 5 | ✅ |
| Factorization | 2 | ✅ |
| Triangular Solvers | 5 | 🔲 |
| Multiple RHS | 2 | 🔲 |
| Condition/Stability | 3 | 🔲 |
| Equilibration | 3 | 🔲 |
| Permutations | 3 | 🔲 |
| Elimination Tree | 4 | 🔲 |
| Matrix Construction | 4 | 🔲 |
| Format Conversion | 8 | 🔲 |
| Matrix Utilities | 5 | 🔲 |
| Properties/Validation | 6 | 🔲 |
| Specialized Solvers | 6 | 🔲 |
| Memory/Statistics | 7 | 🔲 |
| Error Handling | 3 | 🔲 |
| Test Utilities | 1 | 🔲 |
| Types | 11 | 🔲 |
| **Total** | **83** | 🔲 |
