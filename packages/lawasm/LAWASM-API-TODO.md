# LAWasm High-Level TypeScript API - TODO

A comprehensive list of high-level TypeScript API functions for lawasm (LAPACK/BLAS WebAssembly bindings).

---

## 1. Linear System Solvers

- [x] `solve(A, b)` - Solve Ax = b for general matrix A
- [x] `solveTriangular(A, b, opts)` - Solve triangular system (upper/lower)
- [x] `solveBanded(A, b, kl, ku)` - Solve banded system
- [x] `solveSymmetric(A, b)` - Solve symmetric positive definite system
- [x] `solveHermitian(A, b)` - Solve complex Hermitian system
- [x] `solveTridiagonal(dl, d, du, b)` - Solve tridiagonal system

---

## 2. Matrix Factorizations

- [x] `lu(A)` - LU factorization with partial pivoting → {L, U, P}
- [x] `cholesky(A, upper?)` - Cholesky factorization → L or U
- [x] `qr(A, mode?)` - QR factorization → {Q, R} or {R} (reduced/full)
- [x] `qrPivoted(A)` - QR with column pivoting → {Q, R, P}
- [x] `lq(A)` - LQ factorization → {L, Q}
- [x] `ldl(A)` - LDL^T factorization for symmetric matrices
- [x] `schur(A)` - Schur decomposition → {T, Z}
- [x] `hessenberg(A)` - Hessenberg reduction → {H, Q}

---

## 3. Singular Value Decomposition

- [x] `svd(A, full?)` - Full SVD → {U, S, Vt}
- [x] `svdvals(A)` - Singular values only (faster)
- [x] `svdCompact(A)` - Compact/thin SVD
- [x] `svdRank(A, tol?)` - SVD-based rank computation

---

## 4. Eigenvalue Problems

- [x] `eig(A)` - Eigenvalues and eigenvectors of general matrix
- [x] `eigvals(A)` - Eigenvalues only (faster)
- [x] `eigSymmetric(A, opts?)` - Symmetric eigenvalue problem
- [x] `eigHermitian(A, opts?)` - Hermitian eigenvalue problem
- [x] `eigGeneralized(A, B)` - Generalized eigenvalue Ax = λBx
- [x] `eigGeneralizedSymmetric(A, B)` - Symmetric generalized eigenvalue
- [x] `eigBanded(A, kl, ku)` - Eigenvalues of banded matrix
- [x] `eigTridiagonal(d, e)` - Eigenvalues of tridiagonal matrix
- [x] `eigSelect(A, select)` - Selected eigenvalues/eigenvectors

---

## 5. Least Squares & Minimum Norm

- [x] `lstsq(A, b)` - Least squares solution via QR
- [x] `lstsqSVD(A, b, rcond?)` - Least squares via SVD (rank-deficient safe)
- [x] `lstsqGelsy(A, b, rcond?)` - Least squares via QR with pivoting
- [x] `constrainedLstSq(A, b, C, d)` - Equality-constrained least squares
- [x] `generalizedLstSq(A, B, d)` - Generalized least squares

---

## 6. Matrix Inverses & Pseudoinverse

- [x] `inv(A)` - General matrix inverse
- [x] `invTriangular(A, upper?)` - Triangular matrix inverse
- [x] `invSymmetric(A)` - Symmetric positive definite inverse
- [x] `pinv(A, rcond?)` - Moore-Penrose pseudoinverse

---

## 7. Matrix Norms & Condition Numbers

- [x] `norm(A, ord?)` - Matrix norm (1, 2, inf, 'fro')
- [x] `cond(A, ord?)` - Condition number
- [x] `condEst(A)` - Fast condition number estimate
- [x] `rcond(A)` - Reciprocal condition number

---

## 8. Determinants & Rank

- [x] `det(A)` - Determinant
- [x] `logdet(A)` - Log-determinant (for large values)
- [x] `slogdet(A)` - Sign and log-determinant
- [x] `rank(A, tol?)` - Matrix rank

---

## 9. BLAS Level 3 (Matrix-Matrix)

- [x] `matmul(A, B, opts?)` - General matrix multiplication C = αAB + βC
- [x] `matmulTriangular(A, B, side, upper?)` - Triangular matrix multiply
- [x] `solveMatrixTriangular(A, B, side, upper?)` - Solve AX=B or XA=B for triangular A
- [x] `syrk(A, C?, alpha?, beta?)` - Symmetric rank-k update C = αAA^T + βC
- [x] `syr2k(A, B, C?)` - Symmetric rank-2k update
- [x] `herk(A, C?)` - Hermitian rank-k update
- [x] `her2k(A, B, C?)` - Hermitian rank-2k update

---

## 10. BLAS Level 2 (Matrix-Vector)

- [x] `matvec(A, x, opts?)` - Matrix-vector multiply y = αAx + βy
- [x] `matvecTriangular(A, x, upper?)` - Triangular matrix-vector multiply
- [x] `solveVectorTriangular(A, b, upper?)` - Solve Ax=b for triangular A
- [x] `ger(x, y, A?)` - Rank-1 update A = A + αxy^T
- [x] `syr(x, A?)` - Symmetric rank-1 update
- [x] `her(x, A?)` - Hermitian rank-1 update

---

## 11. BLAS Level 1 (Vector)

- [x] `dot(x, y)` - Dot product
- [x] `dotc(x, y)` - Conjugate dot product (complex)
- [x] `axpy(a, x, y)` - y = αx + y
- [x] `scal(a, x)` - x = αx
- [x] `copy(x, y)` - y = x
- [x] `swap(x, y)` - Swap x and y
- [x] `nrm2(x)` - Euclidean norm
- [x] `asum(x)` - Sum of absolute values
- [x] `iamax(x)` - Index of max absolute value

---

## 12. Matrix Utilities

- [x] `transpose(A)` - Matrix transpose
- [x] `conjugate(A)` - Complex conjugate
- [x] `hermitian(A)` - Conjugate transpose
- [x] `triu(A, k?)` - Upper triangular part
- [x] `tril(A, k?)` - Lower triangular part
- [x] `diag(A)` - Extract/create diagonal
- [x] `trace(A)` - Matrix trace
- [x] `balance(A)` - Balance matrix for eigenvalue computation

---

## 13. Matrix Properties & Tests

- [x] `isSymmetric(A, tol?)` - Test for symmetry
- [x] `isHermitian(A, tol?)` - Test for Hermitian
- [x] `isPositiveDefinite(A)` - Test positive definiteness
- [x] `isOrthogonal(A, tol?)` - Test orthogonality
- [x] `isUnitary(A, tol?)` - Test unitarity
- [x] `isSingular(A, tol?)` - Test singularity

---

## 14. Matrix Functions

- [x] `expm(A)` - Matrix exponential
- [x] `logm(A)` - Matrix logarithm
- [x] `sqrtm(A)` - Matrix square root
- [x] `powm(A, n)` - Matrix power
- [x] `funm(A, f)` - General matrix function

---

## 15. Specialized Decompositions

- [x] `polarDecomposition(A)` - A = UP (unitary × positive semidefinite)
- [x] `rrqr(A, tol?)` - Rank-revealing QR
- [x] `csd(X, p, q)` - Cosine-Sine decomposition

---

## 16. Complex-Specific Functions

- [x] `eigHermitian(A)` - Hermitian eigenvalue problem (already in eigenvalues module)
- [x] `choleskyHermitian(A)` - Cholesky for Hermitian matrices
- [x] `solveHermitian(A, b)` - Solve Hermitian system (already in linear-solvers module)

---

## Type Definitions Needed

- [x] `LUResult` - { L: Matrix; U: Matrix; P: number[]; } (in factorizations/types.ts)
- [x] `QRResult` - { Q: Matrix; R: Matrix; } (in factorizations/types.ts)
- [x] `SVDResult` - { U: Matrix; S: number[]; Vt: Matrix; rank?: number; } (in svd/types.ts)
- [x] `EigResult` - { values: Complex[]; vectors?: Matrix; } (in eigenvalues/types.ts)
- [x] `LstSqResult` - { x: Matrix; residuals?: number[]; rank: number; s?: number[]; } (in least-squares/types.ts)
- [x] `CholeskyResult` - { L: Matrix; } (in factorizations/types.ts)
- [x] `SchurResult` - { T: Matrix; Z: Matrix; } (in factorizations/types.ts)
- [x] `HessenbergResult` - { H: Matrix; Q: Matrix; } (in factorizations/types.ts)
- [x] `LDLResult` - { L: Matrix; D: number[]; } (in factorizations/types.ts)
- [x] `PolarResult` - { U: Matrix; P: Matrix; } (in matrix-functions/types.ts as PolarDecompositionResult)

---

## Infrastructure

- [x] Helper functions module (workspace sizing, error messages) - in helpers.ts
- [x] Memory management utilities (auto alloc/free) - allocateDoubles, allocateInts, freeAll in helpers.ts
- [x] Array conversion utilities (TypedArray ↔ WASM memory) - readDoubles, readInts in helpers.ts
- [x] Column-major ↔ row-major conversion helpers - toColumnMajor, fromColumnMajor in helpers.ts
- [x] Complex number utilities - Complex type and handling throughout
- [x] Input validation and error handling - getLapackErrorMessage, validation in each function

---

## Summary

| Category | Count |
|----------|-------|
| Linear System Solvers | 6 |
| Matrix Factorizations | 8 |
| SVD | 4 |
| Eigenvalues | 9 |
| Least Squares | 5 |
| Inverses | 4 |
| Norms & Condition | 4 |
| Determinants & Rank | 4 |
| BLAS Level 3 | 7 |
| BLAS Level 2 | 6 |
| BLAS Level 1 | 9 |
| Matrix Utilities | 8 |
| Matrix Properties | 6 |
| Matrix Functions | 5 |
| Specialized | 3 |
| Complex-Specific | 3 |
| Type Definitions | 10 |
| Infrastructure | 6 |
| **Total** | **~107 items** |
