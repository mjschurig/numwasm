# ARWasm High-Level TypeScript API - Complete Feature List

A comprehensive list of all high-level TypeScript API functions that arwasm could provide for sparse eigenvalue computation.

---

## Currently Implemented ✅

- [x] `eigs()` - Symmetric eigenvalue solver (Implicitly Restarted Lanczos)
- [x] `eign()` - Non-symmetric eigenvalue solver (Implicitly Restarted Arnoldi)
- [x] `eigsh()` - Unified dispatcher that selects eigs/eign based on problem type
- [x] `isEigsResult()` / `isEignResult()` - Type guards for result types
- [x] Module loader functions (`loadARPACKModule`, `getARPACKModule`, etc.)
- [x] Helper utilities (`dsaupdWorklSize`, `defaultNcv`, error message functions)

---

## Core Eigenvalue Solvers

### Complex Eigenvalue Solvers (High Priority) ✅

- [x] `zeigs()` - Complex non-symmetric eigenvalue problem (wraps znaupd/zneupd)
  ```typescript
  async function zeigs(
    matvec: ComplexMatVecFunction,
    n: number,
    nev: number,
    options?: ZeigsOptions
  ): Promise<ZeigsResult>
  ```

- [x] `zeigsh()` - Complex Hermitian eigenvalue problem
  ```typescript
  async function zeigsh(
    matvec: ComplexMatVecFunction,
    n: number,
    nev: number,
    options?: ZeigshOptions
  ): Promise<ZeigshResult>
  ```

### Single Precision Variants (Medium Priority)

- [ ] `seigs()` - Single precision symmetric (wraps ssaupd/sseupd)
  ```typescript
  async function seigs(
    matvec: Float32MatVecFunction,
    n: number,
    nev: number,
    options?: SeigsOptions
  ): Promise<SeigsResult>
  ```

- [ ] `seign()` - Single precision non-symmetric (wraps snaupd/sneupd)
  ```typescript
  async function seign(
    matvec: Float32MatVecFunction,
    n: number,
    nev: number,
    options?: SeignOptions
  ): Promise<SeignResult>
  ```

---

## Singular Value Decomposition (High Priority) ✅

- [x] `svds()` - Partial SVD via ARPACK (like scipy.sparse.linalg.svds)
  ```typescript
  async function svds(
    matvec: MatVecFunction,           // A @ x
    matvecT: MatVecFunction,          // A.T @ x
    m: number,                        // rows
    n: number,                        // cols
    k: number,                        // number of singular values
    options?: SvdsOptions
  ): Promise<SvdsResult>

  interface SvdsResult {
    U: Float64Array[];      // Left singular vectors (m x k)
    s: Float64Array;        // Singular values (k)
    Vt: Float64Array[];     // Right singular vectors (k x n)
    nconv: number;
    niter: number;
  }

  interface SvdsOptions {
    which?: 'LM' | 'SM';              // Largest or smallest singular values
    tol?: number;
    maxiter?: number;
    v0?: Float64Array;
    return_singular_vectors?: boolean | 'u' | 'vh' | 'both';
  }
  ```

---

## Specialized Eigenvalue Problems

### Generalized Eigenvalue Problem (High Priority) ✅

- [x] `geigs()` - Explicit generalized eigenvalue: A @ x = λ * B @ x
  ```typescript
  async function geigs(
    Amatvec: MatVecFunction,
    Bmatvec: MatVecFunction,
    n: number,
    nev: number,
    options?: GeigsOptions
  ): Promise<EigsResult>
  ```

### Shift-Invert Convenience Functions (High Priority) ✅

- [x] `eigsNear()` - Find eigenvalues near sigma (symmetric)
  ```typescript
  async function eigsNear(
    matvec: MatVecFunction,
    solveShifted: OperatorFunction,   // Solves (A - σI)x = b
    n: number,
    nev: number,
    sigma: number,
    options?: EigsOptions
  ): Promise<EigsResult>
  ```

- [x] `eignNear()` - Find eigenvalues near sigma (non-symmetric)
  ```typescript
  async function eignNear(
    matvec: MatVecFunction,
    solveShifted: OperatorFunction,
    n: number,
    nev: number,
    sigma: number | Complex,
    options?: EignOptions
  ): Promise<EignResult>
  ```

### Advanced ARPACK Modes (Low Priority) ✅

- [x] `bucklingEigs()` - Buckling mode: K @ x = λ * Kg @ x (Mode 4)
  ```typescript
  async function bucklingEigs(
    Kmatvec: MatVecFunction,           // Stiffness matrix
    Kgmatvec: MatVecFunction,          // Geometric stiffness matrix
    solveK: OperatorFunction,
    n: number,
    nev: number,
    options?: BucklingOptions
  ): Promise<EigsResult>
  ```

- [x] `cayleyEigs()` - Cayley transform mode (Mode 5)
  ```typescript
  async function cayleyEigs(
    Amatvec: MatVecFunction,
    Bmatvec: MatVecFunction,
    solveCayley: OperatorFunction,
    n: number,
    nev: number,
    sigma: number,
    options?: CayleyOptions
  ): Promise<EigsResult>
  ```

- [x] `eigsInInterval()` - Find eigenvalues in an interval
- [x] `criticalBucklingLoad()` - Find critical buckling load

---

## Application-Specific Functions

### Graph/Network Analysis (Medium Priority) ✅

- [x] `laplacianEigs()` - Graph Laplacian eigenvalues for spectral clustering
  ```typescript
  async function laplacianEigs(
    adjacency: SparseMatrix | MatVecFunction,
    n: number,
    nev: number,
    options?: LaplacianEigsOptions
  ): Promise<LaplacianEigsResult>

  interface LaplacianEigsOptions {
    normalized?: boolean;             // Use normalized Laplacian (default: true)
    which?: 'SM' | 'LM';             // Usually want smallest for clustering
    dropFirst?: boolean;              // Drop trivial eigenvalue 0
  }

  interface LaplacianEigsResult extends EigsResult {
    fiedlerValue?: number;            // Second smallest eigenvalue
    fiedlerVector?: Float64Array;     // Corresponding eigenvector
  }
  ```

- [x] `pagerank()` - PageRank via power iteration / ARPACK
  ```typescript
  async function pagerank(
    outlinks: MatVecFunction,         // Sparse adjacency (column-stochastic)
    n: number,
    options?: PagerankOptions
  ): Promise<PagerankResult>

  interface PagerankOptions {
    damping?: number;                 // Damping factor (default: 0.85)
    tol?: number;
    maxiter?: number;
  }

  interface PagerankResult {
    scores: Float64Array;
    niter: number;
    converged: boolean;
  }
  ```

### Dimensionality Reduction (Medium Priority) ✅

- [x] `spectralEmbedding()` - Laplacian eigenmaps for manifold learning
  ```typescript
  async function spectralEmbedding(
    affinity: MatVecFunction | SparseMatrix,
    n: number,
    nComponents: number,
    options?: SpectralEmbeddingOptions
  ): Promise<SpectralEmbeddingResult>

  interface SpectralEmbeddingOptions {
    normalized?: boolean;
    dropFirst?: boolean;              // Drop first trivial component
  }

  interface SpectralEmbeddingResult {
    embedding: Float64Array[];        // nComponents x n
    eigenvalues: Float64Array;
  }
  ```

- [x] `truncatedPCA()` - Truncated PCA via eigendecomposition of covariance
  ```typescript
  async function truncatedPCA(
    covariance: MatVecFunction | SparseMatrix,
    n: number,
    nComponents: number,
    options?: TruncatedPCAOptions
  ): Promise<TruncatedPCAResult>

  interface TruncatedPCAResult {
    components: Float64Array[];       // Principal components
    explainedVariance: Float64Array;  // Eigenvalues
    explainedVarianceRatio: Float64Array;
  }
  ```

---

## Linear Algebra Utilities

### Matrix Properties (Medium Priority) ✅

- [x] `spectralRadius()` - Largest magnitude eigenvalue
  ```typescript
  async function spectralRadius(
    matvec: MatVecFunction,
    n: number,
    options?: SpectralRadiusOptions
  ): Promise<SpectralRadiusResult>
  ```

- [x] `spectralNorm()` - Largest singular value (2-norm)
  ```typescript
  async function spectralNorm(
    matvec: MatVecFunction,
    matvecT: MatVecFunction,
    m: number,
    n: number,
    options?: SpectralNormOptions
  ): Promise<SpectralNormResult>
  ```

- [x] `condest()` - Condition number estimate
  ```typescript
  async function condest(
    matvec: MatVecFunction,
    matvecT: MatVecFunction,
    n: number,
    options?: CondestOptions
  ): Promise<CondestResult>
  ```

- [x] `nuclearNormApprox()` - Approximate nuclear norm (sum of singular values)
  ```typescript
  async function nuclearNormApprox(
    matvec: MatVecFunction,
    matvecT: MatVecFunction,
    m: number,
    n: number,
    k: number,                        // Number of singular values to sum
    options?: NuclearNormOptions
  ): Promise<NuclearNormResult>
  ```

### Matrix Functions (Low Priority) ✅

- [x] `expmv()` - Matrix exponential times vector: exp(t*A) @ v
  ```typescript
  function expmv(
    matvec: MatVecFunction,
    v: Float64Array,
    t: number,
    options?: ExpmvOptions
  ): ExpmvResult
  ```

- [x] `sqrtmv()` - Matrix square root times vector (for SPD matrices)
  ```typescript
  async function sqrtmv(
    matvec: MatVecFunction,
    v: Float64Array,
    options?: SqrtmvOptions
  ): Promise<SqrtmvResult>
  ```

- [x] `invsqrtmv()` - Inverse matrix square root times vector
- [x] `matpowv()` - Matrix power times vector: A^p @ v
- [x] `expmvMultiple()` - Matrix exponential for multiple time points

---

## Iterative Refinement & Continuation

### Deflation & Hot Start (Low Priority) ✅

- [x] `eigsContinue()` - Continue from previous computation
  ```typescript
  async function eigsContinue(
    matvec: MatVecFunction,
    n: number,
    additionalNev: number,
    previousResult: EigsResult,
    options?: EigsDeflatedOptions
  ): Promise<EigsResult>
  ```

- [x] `eigsDeflated()` - Find eigenvalues orthogonal to known eigenvectors
  ```typescript
  async function eigsDeflated(
    matvec: MatVecFunction,
    n: number,
    nev: number,
    knownEigenvectors: Float64Array[],
    knownEigenvalues?: Float64Array,
    options?: EigsDeflatedOptions
  ): Promise<EigsResult>
  ```

- [x] `eignDeflated()` - Non-symmetric deflated eigenvalue solver

### Block Methods (Low Priority)

- [ ] `eigsBlock()` - LOBPCG-style block eigenvalue solver
  ```typescript
  async function eigsBlock(
    matvec: BlockMatVecFunction,
    n: number,
    nev: number,
    blockSize?: number,
    options?: EigsBlockOptions
  ): Promise<EigsResult>

  type BlockMatVecFunction = (X: Float64Array[]) => Float64Array[];
  ```

---

## Sparse Matrix Helpers (High Priority) ✅

### Create MatVec from Sparse Formats

- [x] `csrMatvec()` - CSR (Compressed Sparse Row) format
  ```typescript
  function csrMatvec(
    indptr: Int32Array,
    indices: Int32Array,
    data: Float64Array,
    shape: [number, number]
  ): MatVecFunction
  ```

- [x] `cscMatvec()` - CSC (Compressed Sparse Column) format
  ```typescript
  function cscMatvec(
    indptr: Int32Array,
    indices: Int32Array,
    data: Float64Array,
    shape: [number, number]
  ): MatVecFunction
  ```

- [x] `cooMatvec()` - COO (Coordinate) format
  ```typescript
  function cooMatvec(
    rows: Int32Array,
    cols: Int32Array,
    values: Float64Array,
    shape: [number, number]
  ): MatVecFunction
  ```

- [x] `denseMatvec()` - Dense array (for small problems)
  ```typescript
  function denseMatvec(
    matrix: Float64Array,
    m: number,
    n: number,
    rowMajor?: boolean
  ): MatVecFunction
  ```

- [x] `diagMatvec()` - Diagonal matrix
  ```typescript
  function diagMatvec(diagonal: Float64Array): MatVecFunction
  ```

- [x] `tridiagMatvec()` - Tridiagonal matrix
  ```typescript
  function tridiagMatvec(
    lower: Float64Array,
    diag: Float64Array,
    upper: Float64Array
  ): MatVecFunction
  ```

- [x] `bandedMatvec()` - General banded matrix
  ```typescript
  function bandedMatvec(
    bands: Float64Array[],
    offsets: number[],
    n: number
  ): MatVecFunction
  ```

---

## Operator Combinators (Medium Priority) ✅

- [x] `addMatvec()` - Linear combination: (α*A + β*B) @ x
  ```typescript
  function addMatvec(
    A: MatVecFunction,
    B: MatVecFunction,
    alpha?: number,
    beta?: number
  ): MatVecFunction
  ```

- [x] `mulMatvec()` - Product: A @ B @ x
  ```typescript
  function mulMatvec(
    A: MatVecFunction,
    B: MatVecFunction
  ): MatVecFunction
  ```

- [x] `shiftMatvec()` - Shift: (A - σI) @ x
  ```typescript
  function shiftMatvec(
    matvec: MatVecFunction,
    sigma: number
  ): MatVecFunction
  ```

- [x] `scaleMatvec()` - Scale: α * A @ x
  ```typescript
  function scaleMatvec(
    matvec: MatVecFunction,
    alpha: number
  ): MatVecFunction
  ```

- [x] `transposeMatvec()` - Transpose via explicit computation
  ```typescript
  function transposeMatvec(
    matvec: MatVecFunction,
    m: number,
    n: number
  ): MatVecFunction
  ```

- [x] `symmetrizeMatvec()` - (A + A.T) / 2
  ```typescript
  function symmetrizeMatvec(
    matvec: MatVecFunction,
    matvecT: MatVecFunction
  ): MatVecFunction
  ```

Also implemented:
- [x] `identityMatvec()` - Identity operator
- [x] `negateMatvec()` - Negation operator
- [x] `powerMatvec()` - Power operator A^k
- [x] `blockDiagMatvec()` - Block diagonal operator

---

## Validation & Diagnostics (Medium Priority) ✅

### Result Verification

- [x] `verifyEigs()` - Verify eigenvalue/eigenvector pairs
  ```typescript
  function verifyEigs(
    matvec: MatVecFunction,
    eigenvalues: Float64Array,
    eigenvectors: Float64Array[],
    options?: VerifyOptions
  ): VerifyEigsResult

  interface VerifyEigsResult {
    residuals: Float64Array;          // ||Av - λv|| for each pair
    maxResidual: number;
    relativeResiduals: Float64Array;  // ||Av - λv|| / |λ|
    orthogonality: number;            // max |v_i . v_j| for i ≠ j
    isValid: boolean;
  }
  ```

- [x] `verifyEign()` - Verify non-symmetric eigenvalue/eigenvector pairs

- [x] `verifySvds()` - Verify singular value decomposition
  ```typescript
  function verifySvds(
    matvec: MatVecFunction,
    matvecT: MatVecFunction,
    U: Float64Array[],
    s: Float64Array,
    Vt: Float64Array[],
    options?: VerifyOptions
  ): VerifySvdsResult
  ```

### Matrix Properties

- [x] `checkSymmetry()` - Check if matrix is symmetric via random probing
  ```typescript
  function checkSymmetry(
    matvec: MatVecFunction,
    n: number,
    nProbes?: number
  ): SymmetryCheckResult
  ```

- [x] `checkPositiveDefinite()` - Check if matrix is positive definite
  ```typescript
  async function checkPositiveDefinite(
    matvec: MatVecFunction,
    n: number,
    tol?: number
  ): Promise<PositiveDefiniteCheckResult>
  ```

- [x] `checkNormalization()` - Check eigenvector normalization

### Convergence Monitoring (Not implemented - requires ARPACK internals)

- [ ] Add `onIteration` callback to all solver options
  ```typescript
  interface ConvergenceCallback {
    (info: {
      iteration: number;
      ritzValues: Float64Array;
      residualNorms: Float64Array;
      converged: number;
    }): boolean | void;               // Return false to stop early
  }

  interface EigsOptionsWithCallback extends EigsOptions {
    onIteration?: ConvergenceCallback;
    returnHistory?: boolean;
  }
  ```

- [ ] Add convergence history to results
  ```typescript
  interface EigsResultWithHistory extends EigsResult {
    history?: {
      iterations: number[];
      ritzValues: Float64Array[];
      residualNorms: Float64Array[];
    };
  }
  ```

---

## Type Definitions

### Core Types (to add)

```typescript
// Complex number support
interface Complex {
  re: number;
  im: number;
}

type ComplexArray = {
  re: Float64Array;
  im: Float64Array;
};

// Function types
type MatVecFunction = (x: Float64Array) => Float64Array;
type Float32MatVecFunction = (x: Float32Array) => Float32Array;
type ComplexMatVecFunction = (x: ComplexArray) => ComplexArray;
type BlockMatVecFunction = (X: Float64Array[]) => Float64Array[];
type OperatorFunction = (b: Float64Array) => Float64Array;

// Sparse matrix type (for convenience wrappers)
interface SparseMatrix {
  format: 'csr' | 'csc' | 'coo';
  data: Float64Array;
  indices: Int32Array;
  indptr?: Int32Array;                // For CSR/CSC
  rows?: Int32Array;                  // For COO
  cols?: Int32Array;                  // For COO
  shape: [number, number];
}

// Complex result types
interface ZeigsResult {
  eigenvalues: ComplexArray;
  eigenvectors?: ComplexArray[];
  nconv: number;
  niter: number;
  info: number;
}
```

---

## Priority Summary

| Priority | Category | Functions | Status |
|----------|----------|-----------|--------|
| **High** | Core Eigensolvers | `zeigs`, `zeigsh` | ✅ Done |
| **High** | SVD | `svds` | ✅ Done |
| **High** | Generalized/Shift-Invert | `geigs`, `eigsNear`, `eignNear` | ✅ Done |
| **High** | Sparse Helpers | `csrMatvec`, `cscMatvec`, `cooMatvec`, `denseMatvec`, `diagMatvec`, `tridiagMatvec`, `bandedMatvec` | ✅ Done |
| **Medium** | Single Precision | `seigs`, `seign` | ⏳ Pending |
| **Medium** | Graph/ML | `laplacianEigs`, `spectralEmbedding`, `pagerank` | ✅ Done |
| **Medium** | Matrix Properties | `spectralRadius`, `spectralNorm`, `condest`, `nuclearNormApprox` | ✅ Done |
| **Medium** | Operator Combinators | `addMatvec`, `shiftMatvec`, `scaleMatvec`, `mulMatvec`, `symmetrizeMatvec`, `transposeMatvec` | ✅ Done |
| **Medium** | Validation | `verifyEigs`, `verifySvds`, `checkSymmetry`, `checkPositiveDefinite` | ✅ Done |
| **Low** | Advanced Modes | `bucklingEigs`, `cayleyEigs`, `eigsInInterval`, `criticalBucklingLoad` | ✅ Done |
| **Low** | Matrix Functions | `expmv`, `expmvMultiple`, `sqrtmv`, `invsqrtmv`, `matpowv` | ✅ Done |
| **Low** | Continuation | `eigsContinue`, `eigsDeflated`, `eignDeflated` | ✅ Done |
| **Low** | Block Methods | `eigsBlock` | ⏳ Pending |

---

## Implementation Notes

1. **WASM already exports**: The current WASM build exports `znaupd/zneupd` and `cnaupd/cneupd` for complex eigenvalues, so implementing `zeigs`/`zeigsh` is straightforward.

2. **SVD implementation**: Can be implemented by computing eigenvalues of A^T A or A A^T using the existing `eigs` function.

3. **Sparse matrix helpers**: These are pure TypeScript and don't require any WASM changes.

4. **Memory considerations**: Complex and block operations will need careful WASM heap management for interleaved real/imaginary arrays.

5. **Integration with numwasm**: Consider having sparse helpers accept numwasm NDArray types directly.

---

## References

- [SciPy sparse.linalg Reference](https://docs.scipy.org/doc/scipy/reference/sparse.linalg.html)
- [SciPy ARPACK Tutorial](https://docs.scipy.org/doc/scipy/tutorial/arpack.html)
- [Spectra C++ Library](https://spectralib.org/)
- [Julia Arpack.jl API](https://julialinearalgebra.github.io/Arpack.jl/stable/api/)
- [ARPACK-NG GitHub](https://github.com/opencollab/arpack-ng)
- [scikit-learn SpectralClustering](https://scikit-learn.org/stable/modules/generated/sklearn.cluster.SpectralClustering.html)
