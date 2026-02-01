# LINWASM Full API Reference

A comprehensive list of high-level TypeScript API functions for the linwasm package (LINPACK WebAssembly).

LINPACK is the classic Fortran library for numerical linear algebra. This document describes all high-level TypeScript wrapper functions that should be implemented.

---

## Summary

| Category | Count |
|----------|-------|
| Linear System Solvers | 14 |
| Matrix Factorizations | 10 |
| Matrix Inverses | 8 |
| Determinants | 6 |
| Condition Number Estimation | 8 |
| SVD | 4 |
| QR Decomposition | 4 |
| Least Squares | 4 |
| Cholesky Updates | 4 |
| Triangular Operations | 4 |
| Matrix Inertia | 2 |
| BLAS Level 1 (Vector) | 16 |
| Matrix Utilities | 12 |
| Matrix Properties | 8 |
| **Total** | **~104** |

---

## 1. Linear System Solvers

### General Dense Matrices
- [ ] `solve(A, b, options?)` - Solve Ax = b for general matrix (LU factorization)
- [ ] `solveTranspose(A, b, options?)` - Solve A'x = b (transpose system)
- [ ] `solveComplex(A, b, options?)` - Solve complex system (single precision)
- [ ] `solveComplexDouble(A, b, options?)` - Solve complex system (double precision)

### Banded Matrices
- [ ] `solveBanded(A, b, ml, mu, options?)` - Solve banded system Ax = b
- [ ] `solveBandedTranspose(A, b, ml, mu, options?)` - Solve banded A'x = b

### Symmetric/Hermitian Positive Definite
- [ ] `solvePositiveDefinite(A, b, options?)` - Solve SPD system (Cholesky)
- [ ] `solvePositiveDefinitePacked(A, b, options?)` - Solve SPD system (packed storage)
- [ ] `solvePositiveDefiniteBanded(A, b, m, options?)` - Solve banded SPD system
- [ ] `solveHermitian(A, b, options?)` - Solve Hermitian positive definite (complex)

### Symmetric Indefinite
- [ ] `solveSymmetricIndefinite(A, b, options?)` - Solve symmetric indefinite (Bunch-Kaufman)
- [ ] `solveSymmetricIndefinitePacked(A, b, options?)` - Solve symmetric indefinite (packed)

### Triangular Systems
- [ ] `solveTriangular(A, b, options?)` - Solve triangular system Tx = b or T'x = b

### Tridiagonal Systems
- [ ] `solveTridiagonal(d, e, f, b, options?)` - Solve general tridiagonal system
- [ ] `solveTridiagonalPositiveDefinite(d, e, b, options?)` - Solve SPD tridiagonal

---

## 2. Matrix Factorizations

### LU Factorization
- [ ] `lu(A, options?)` - LU factorization with partial pivoting → {L, U, P, info}
- [ ] `luBanded(A, ml, mu, options?)` - LU factorization for banded matrix

### Cholesky Factorization
- [ ] `cholesky(A, options?)` - Cholesky factorization A = R'R → {R, info}
- [ ] `choleskyPacked(A, options?)` - Cholesky for packed storage
- [ ] `choleskyBanded(A, m, options?)` - Cholesky for banded SPD matrix
- [ ] `choleskyPivoted(A, options?)` - Cholesky with column pivoting (DCHDC)

### Symmetric Indefinite Factorization
- [ ] `ldl(A, options?)` - Bunch-Kaufman LDL' factorization → {L, D, P, info}
- [ ] `ldlPacked(A, options?)` - Bunch-Kaufman for packed storage

### QR Factorization
- [ ] `qr(A, options?)` - QR factorization with optional pivoting → {Q, R, P?, info}

### SVD
- [ ] `svd(A, options?)` - Singular value decomposition → {U, S, V, info}

---

## 3. Matrix Inverses

### General Matrix Inverse
- [ ] `inv(A, options?)` - General matrix inverse (via LU)
- [ ] `invComplex(A, options?)` - Complex matrix inverse

### Positive Definite Inverse
- [ ] `invPositiveDefinite(A, options?)` - SPD matrix inverse (via Cholesky)
- [ ] `invPositiveDefinitePacked(A, options?)` - SPD inverse (packed storage)

### Symmetric Indefinite Inverse
- [ ] `invSymmetricIndefinite(A, options?)` - Symmetric indefinite inverse (Bunch-Kaufman)
- [ ] `invSymmetricIndefinitePacked(A, options?)` - Symmetric indefinite inverse (packed)

### Triangular Inverse
- [ ] `invTriangular(A, options?)` - Triangular matrix inverse

### Pseudoinverse
- [ ] `pinv(A, rcond?, options?)` - Moore-Penrose pseudoinverse (via SVD)

---

## 4. Determinants

- [ ] `det(A, options?)` - Determinant of general matrix (via LU)
- [ ] `detBanded(A, ml, mu, options?)` - Determinant of banded matrix
- [ ] `detPositiveDefinite(A, options?)` - Determinant of SPD matrix (via Cholesky)
- [ ] `detPositiveDefinitePacked(A, options?)` - Determinant of SPD (packed)
- [ ] `detPositiveDefiniteBanded(A, m, options?)` - Determinant of banded SPD
- [ ] `detSymmetricIndefinite(A, options?)` - Determinant of symmetric indefinite
- [ ] `slogdet(A, options?)` - Sign and log-determinant (numerically stable)

---

## 5. Condition Number Estimation

- [ ] `rcond(A, options?)` - Reciprocal condition number of general matrix
- [ ] `rcondBanded(A, ml, mu, options?)` - Reciprocal condition of banded matrix
- [ ] `rcondPositiveDefinite(A, options?)` - Reciprocal condition of SPD matrix
- [ ] `rcondPositiveDefinitePacked(A, options?)` - Reciprocal condition of SPD (packed)
- [ ] `rcondPositiveDefiniteBanded(A, m, options?)` - Reciprocal condition of banded SPD
- [ ] `rcondSymmetricIndefinite(A, options?)` - Reciprocal condition of symmetric indefinite
- [ ] `rcondSymmetricIndefinitePacked(A, options?)` - Reciprocal condition (packed)
- [ ] `rcondTriangular(A, options?)` - Reciprocal condition of triangular matrix
- [ ] `cond(A, options?)` - Condition number (1/rcond)

---

## 6. Singular Value Decomposition

- [ ] `svd(A, options?)` - Full SVD: A = U * S * V' → {U, S, V, info}
- [ ] `svdvals(A, options?)` - Singular values only → {S, info}
- [ ] `svdCompact(A, options?)` - Compact SVD (economy size)
- [ ] `svdRank(A, tol?, options?)` - Numerical rank via SVD

---

## 7. QR Decomposition

- [ ] `qr(A, options?)` - QR factorization → {Q, R, info}
- [ ] `qrPivoted(A, options?)` - QR with column pivoting → {Q, R, P, info}
- [ ] `qrApply(qrFactors, y, options?)` - Apply Q or Q' to vector
- [ ] `qrSolve(qrFactors, b, options?)` - Solve least squares via QR

---

## 8. Least Squares

- [ ] `lstsq(A, b, options?)` - Least squares solution min||Ax - b||₂ (via QR)
- [ ] `lstsqSVD(A, b, rcond?, options?)` - Least squares via SVD (rank-deficient safe)
- [ ] `lstsqResidual(A, b, options?)` - Least squares with residual
- [ ] `lstsqFitted(A, b, options?)` - Least squares with fitted values

---

## 9. Cholesky Update Routines

- [ ] `choleskyUpdate(R, x, options?)` - Rank-1 update: R'R + xx' (DCHUD)
- [ ] `choleskyDowndate(R, x, options?)` - Rank-1 downdate: R'R - xx' (DCHDD)
- [ ] `choleskyExchange(R, k, l, options?)` - Reorder columns in Cholesky factor (DCHEX)
- [ ] `choleskyPivoted(A, options?)` - Cholesky with pivoting (DCHDC)

---

## 10. Triangular Matrix Operations

- [ ] `solveTriangular(T, b, options?)` - Solve Tx = b or T'x = b
- [ ] `invTriangular(T, options?)` - Inverse of triangular matrix
- [ ] `rcondTriangular(T, options?)` - Condition number of triangular matrix
- [ ] `detTriangular(T, options?)` - Determinant of triangular matrix

---

## 11. Matrix Inertia

- [ ] `inertia(A, options?)` - Inertia of symmetric matrix: counts of (neg, zero, pos) eigenvalues
- [ ] `inertiaPacked(A, options?)` - Inertia of symmetric matrix (packed storage)

---

## 12. BLAS Level 1 (Vector Operations)

### Double Precision
- [ ] `dot(x, y)` - Dot product x'y
- [ ] `nrm2(x)` - Euclidean norm ||x||₂
- [ ] `asum(x)` - Sum of absolute values ||x||₁
- [ ] `iamax(x)` - Index of max absolute value
- [ ] `axpy(alpha, x, y)` - y := αx + y
- [ ] `scal(alpha, x)` - x := αx
- [ ] `copy(x, y)` - y := x
- [ ] `swap(x, y)` - Swap x ↔ y
- [ ] `rotg(a, b)` - Generate Givens rotation
- [ ] `rot(x, y, c, s)` - Apply Givens rotation

### Complex
- [ ] `dotc(x, y)` - Conjugated dot product conj(x)'y
- [ ] `dotu(x, y)` - Unconjugated dot product x'y
- [ ] `scnrm2(x)` - Complex Euclidean norm (single)
- [ ] `dznrm2(x)` - Complex Euclidean norm (double)
- [ ] `icamax(x)` - Index of max |x_i| (complex single)
- [ ] `izamax(x)` - Index of max |x_i| (complex double)

---

## 13. Matrix Utilities

### Creation
- [ ] `zeros(m, n, dtype?)` - Zero matrix
- [ ] `ones(m, n, dtype?)` - Matrix of ones
- [ ] `eye(n, dtype?)` - Identity matrix
- [ ] `diag(v, k?)` - Create diagonal matrix from vector

### Extraction
- [ ] `diagonal(A, k?)` - Extract diagonal
- [ ] `triu(A, k?)` - Extract upper triangular
- [ ] `tril(A, k?)` - Extract lower triangular

### Transformations
- [ ] `transpose(A)` - Matrix transpose
- [ ] `hermitian(A)` - Conjugate transpose (complex)
- [ ] `conjugate(A)` - Complex conjugate

### Conversion
- [ ] `toPacked(A)` - Convert to packed storage (upper triangle)
- [ ] `fromPacked(ap, n)` - Convert from packed storage

---

## 14. Matrix Properties & Tests

- [ ] `isSymmetric(A, tol?)` - Test for symmetry
- [ ] `isHermitian(A, tol?)` - Test for Hermitian (complex)
- [ ] `isPositiveDefinite(A)` - Test for positive definiteness (via Cholesky)
- [ ] `isTriangular(A, upper?, tol?)` - Test for triangular
- [ ] `isSingular(A, tol?)` - Test for singularity
- [ ] `bandwidth(A, tol?)` - Compute bandwidth {ml, mu}
- [ ] `rank(A, tol?)` - Numerical rank
- [ ] `trace(A)` - Matrix trace

---

## 15. Norms

- [ ] `norm(A, ord?)` - Matrix norm (1, inf, fro)
- [ ] `normFrobenius(A)` - Frobenius norm
- [ ] `norm1(A)` - 1-norm (max column sum)
- [ ] `normInf(A)` - Infinity norm (max row sum)
- [ ] `vectorNorm(x, ord?)` - Vector norm (1, 2, inf)

---

## Options Interfaces

```typescript
interface SolveOptions {
  overwriteA?: boolean;      // Allow modifying input matrix
  overwriteB?: boolean;      // Allow modifying input RHS
  transpose?: boolean;       // Solve A'x = b instead of Ax = b
  dtype?: 'float32' | 'float64' | 'complex64' | 'complex128';
}

interface FactorizationOptions {
  overwriteA?: boolean;      // Allow modifying input matrix
  pivoting?: boolean;        // Use column pivoting (for QR)
  computeQ?: boolean;        // Compute Q matrix explicitly
  dtype?: 'float32' | 'float64' | 'complex64' | 'complex128';
}

interface SVDOptions {
  overwriteA?: boolean;
  computeU?: boolean | 'full' | 'reduced';
  computeV?: boolean | 'full' | 'reduced';
  dtype?: 'float32' | 'float64' | 'complex64' | 'complex128';
}

interface TriangularOptions {
  upper?: boolean;           // Upper (true) or lower (false) triangular
  unitDiagonal?: boolean;    // Assume diagonal is all 1s
  transpose?: boolean;       // Use transpose
}
```

---

## Result Interfaces

```typescript
interface SolveResult {
  x: Float64Array;           // Solution vector
  info: number;              // LINPACK status code
  success: boolean;          // info === 0
  message: string;           // Human-readable status
}

interface LUResult {
  L: Float64Array;           // Lower triangular factor
  U: Float64Array;           // Upper triangular factor
  P: Int32Array;             // Pivot permutation
  info: number;
  success: boolean;
}

interface CholeskyResult {
  R: Float64Array;           // Upper triangular factor (A = R'R)
  info: number;
  success: boolean;
  message: string;
}

interface SVDResult {
  U?: Float64Array;          // Left singular vectors (m×k)
  S: Float64Array;           // Singular values (descending)
  V?: Float64Array;          // Right singular vectors (n×k)
  info: number;
  success: boolean;
}

interface QRResult {
  Q: Float64Array;           // Orthogonal matrix
  R: Float64Array;           // Upper triangular
  P?: Int32Array;            // Column permutation (if pivoting)
  qraux: Float64Array;       // Auxiliary data for Q reconstruction
  info: number;
  success: boolean;
}

interface ConditionResult {
  rcond: number;             // Reciprocal condition number
  cond: number;              // Condition number (1/rcond)
  isWellConditioned: boolean; // rcond > machine epsilon
}

interface InertiaResult {
  negative: number;          // Count of negative eigenvalues
  zero: number;              // Count of zero eigenvalues
  positive: number;          // Count of positive eigenvalues
}
```

---

## Precision Support

LINPACK provides four precision variants:

| Prefix | Type | TypedArray | Bytes | Precision |
|--------|------|------------|-------|-----------|
| `s` | single real | Float32Array | 4 | ~7 digits |
| `d` | double real | Float64Array | 8 | ~16 digits |
| `c` | single complex | Float32Array (pairs) | 8 | ~7 digits |
| `z` | double complex | Float64Array (pairs) | 16 | ~16 digits |

High-level API functions should support a `dtype` option to select precision.

---

## Matrix Storage Formats

### Column-Major (Default)
Standard Fortran/LINPACK storage. Element A[i,j] at index `i + j*lda`.

### Packed Storage
For symmetric matrices. Upper triangle stored by columns:
```
A[i,j] (i ≤ j) at index k = i + j*(j+1)/2
```

### Band Storage
For banded matrices with `ml` subdiagonals and `mu` superdiagonals:
```
A[i,j] stored at ABD[i-j+ml+mu, j]
Leading dimension = 2*ml + mu + 1
```
