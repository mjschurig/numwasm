# LAWasm High-Level TypeScript API - TODO

A comprehensive list of high-level TypeScript API functions for lawasm (LAPACK/BLAS WebAssembly bindings).

---

## 1. Linear System Solvers

- [ ] `solve(A, b)` - Solve Ax = b for general matrix A
- [ ] `solveTriangular(A, b, opts)` - Solve triangular system (upper/lower)
- [ ] `solveBanded(A, b, kl, ku)` - Solve banded system
- [ ] `solveSymmetric(A, b)` - Solve symmetric positive definite system
- [ ] `solveHermitian(A, b)` - Solve complex Hermitian system
- [ ] `solveTridiagonal(dl, d, du, b)` - Solve tridiagonal system

---

## 2. Matrix Factorizations

- [ ] `lu(A)` - LU factorization with partial pivoting → {L, U, P}
- [ ] `cholesky(A, upper?)` - Cholesky factorization → L or U
- [ ] `qr(A, mode?)` - QR factorization → {Q, R} or {R} (reduced/full)
- [ ] `qrPivoted(A)` - QR with column pivoting → {Q, R, P}
- [ ] `lq(A)` - LQ factorization → {L, Q}
- [ ] `ldl(A)` - LDL^T factorization for symmetric matrices
- [ ] `schur(A)` - Schur decomposition → {T, Z}
- [ ] `hessenberg(A)` - Hessenberg reduction → {H, Q}

---

## 3. Singular Value Decomposition

- [ ] `svd(A, full?)` - Full SVD → {U, S, Vt}
- [ ] `svdvals(A)` - Singular values only (faster)
- [ ] `svdCompact(A)` - Compact/thin SVD
- [ ] `svdRank(A, tol?)` - SVD-based rank computation

---

## 4. Eigenvalue Problems

- [ ] `eig(A)` - Eigenvalues and eigenvectors of general matrix
- [ ] `eigvals(A)` - Eigenvalues only (faster)
- [ ] `eigSymmetric(A, opts?)` - Symmetric eigenvalue problem
- [ ] `eigHermitian(A, opts?)` - Hermitian eigenvalue problem
- [ ] `eigGeneralized(A, B)` - Generalized eigenvalue Ax = λBx
- [ ] `eigGeneralizedSymmetric(A, B)` - Symmetric generalized eigenvalue
- [ ] `eigBanded(A, kl, ku)` - Eigenvalues of banded matrix
- [ ] `eigTridiagonal(d, e)` - Eigenvalues of tridiagonal matrix
- [ ] `eigSelect(A, select)` - Selected eigenvalues/eigenvectors

---

## 5. Least Squares & Minimum Norm

- [ ] `lstsq(A, b)` - Least squares solution via QR
- [ ] `lstsqSVD(A, b, rcond?)` - Least squares via SVD (rank-deficient safe)
- [ ] `lstsqGelsy(A, b, rcond?)` - Least squares via QR with pivoting
- [ ] `constrainedLstSq(A, b, C, d)` - Equality-constrained least squares
- [ ] `generalizedLstSq(A, B, d)` - Generalized least squares

---

## 6. Matrix Inverses & Pseudoinverse

- [ ] `inv(A)` - General matrix inverse
- [ ] `invTriangular(A, upper?)` - Triangular matrix inverse
- [ ] `invSymmetric(A)` - Symmetric positive definite inverse
- [ ] `pinv(A, rcond?)` - Moore-Penrose pseudoinverse

---

## 7. Matrix Norms & Condition Numbers

- [ ] `norm(A, ord?)` - Matrix norm (1, 2, inf, 'fro')
- [ ] `cond(A, ord?)` - Condition number
- [ ] `condEst(A)` - Fast condition number estimate
- [ ] `rcond(A)` - Reciprocal condition number

---

## 8. Determinants & Rank

- [ ] `det(A)` - Determinant
- [ ] `logdet(A)` - Log-determinant (for large values)
- [ ] `slogdet(A)` - Sign and log-determinant
- [ ] `rank(A, tol?)` - Matrix rank

---

## 9. BLAS Level 3 (Matrix-Matrix)

- [ ] `matmul(A, B, opts?)` - General matrix multiplication C = αAB + βC
- [ ] `matmulTriangular(A, B, side, upper?)` - Triangular matrix multiply
- [ ] `solveMatrixTriangular(A, B, side, upper?)` - Solve AX=B or XA=B for triangular A
- [ ] `syrk(A, C?, alpha?, beta?)` - Symmetric rank-k update C = αAA^T + βC
- [ ] `syr2k(A, B, C?)` - Symmetric rank-2k update
- [ ] `herk(A, C?)` - Hermitian rank-k update
- [ ] `her2k(A, B, C?)` - Hermitian rank-2k update

---

## 10. BLAS Level 2 (Matrix-Vector)

- [ ] `matvec(A, x, opts?)` - Matrix-vector multiply y = αAx + βy
- [ ] `matvecTriangular(A, x, upper?)` - Triangular matrix-vector multiply
- [ ] `solveVectorTriangular(A, b, upper?)` - Solve Ax=b for triangular A
- [ ] `ger(x, y, A?)` - Rank-1 update A = A + αxy^T
- [ ] `syr(x, A?)` - Symmetric rank-1 update
- [ ] `her(x, A?)` - Hermitian rank-1 update

---

## 11. BLAS Level 1 (Vector)

- [ ] `dot(x, y)` - Dot product
- [ ] `dotc(x, y)` - Conjugate dot product (complex)
- [ ] `axpy(a, x, y)` - y = αx + y
- [ ] `scal(a, x)` - x = αx
- [ ] `copy(x, y)` - y = x
- [ ] `swap(x, y)` - Swap x and y
- [ ] `nrm2(x)` - Euclidean norm
- [ ] `asum(x)` - Sum of absolute values
- [ ] `iamax(x)` - Index of max absolute value

---

## 12. Matrix Utilities

- [ ] `transpose(A)` - Matrix transpose
- [ ] `conjugate(A)` - Complex conjugate
- [ ] `hermitian(A)` - Conjugate transpose
- [ ] `triu(A, k?)` - Upper triangular part
- [ ] `tril(A, k?)` - Lower triangular part
- [ ] `diag(A)` - Extract/create diagonal
- [ ] `trace(A)` - Matrix trace
- [ ] `balance(A)` - Balance matrix for eigenvalue computation

---

## 13. Matrix Properties & Tests

- [ ] `isSymmetric(A, tol?)` - Test for symmetry
- [ ] `isHermitian(A, tol?)` - Test for Hermitian
- [ ] `isPositiveDefinite(A)` - Test positive definiteness
- [ ] `isOrthogonal(A, tol?)` - Test orthogonality
- [ ] `isUnitary(A, tol?)` - Test unitarity
- [ ] `isSingular(A, tol?)` - Test singularity

---

## 14. Matrix Functions

- [ ] `expm(A)` - Matrix exponential
- [ ] `logm(A)` - Matrix logarithm
- [ ] `sqrtm(A)` - Matrix square root
- [ ] `powm(A, n)` - Matrix power
- [ ] `funm(A, f)` - General matrix function

---

## 15. Specialized Decompositions

- [ ] `polarDecomposition(A)` - A = UP (unitary × positive semidefinite)
- [ ] `rrqr(A, tol?)` - Rank-revealing QR
- [ ] `csd(X, p, q)` - Cosine-Sine decomposition

---

## 16. Complex-Specific Functions

- [ ] `eigHermitian(A)` - Hermitian eigenvalue problem
- [ ] `choleskyHermitian(A)` - Cholesky for Hermitian matrices
- [ ] `solveHermitian(A, b)` - Solve Hermitian system

---

## Type Definitions Needed

- [ ] `LUResult` - { L: Matrix; U: Matrix; P: number[]; }
- [ ] `QRResult` - { Q: Matrix; R: Matrix; }
- [ ] `SVDResult` - { U: Matrix; S: number[]; Vt: Matrix; rank?: number; }
- [ ] `EigResult` - { values: Complex[]; vectors?: Matrix; }
- [ ] `LstSqResult` - { x: Matrix; residuals?: number[]; rank: number; s?: number[]; }
- [ ] `CholeskyResult` - { L: Matrix; } (or U depending on upper flag)
- [ ] `SchurResult` - { T: Matrix; Z: Matrix; }
- [ ] `HessenbergResult` - { H: Matrix; Q: Matrix; }
- [ ] `LDLResult` - { L: Matrix; D: number[]; }
- [ ] `PolarResult` - { U: Matrix; P: Matrix; }

---

## Infrastructure

- [ ] Helper functions module (workspace sizing, error messages)
- [ ] Memory management utilities (auto alloc/free)
- [ ] Array conversion utilities (TypedArray ↔ WASM memory)
- [ ] Column-major ↔ row-major conversion helpers
- [ ] Complex number utilities
- [ ] Input validation and error handling

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
