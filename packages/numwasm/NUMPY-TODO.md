# NumWasm NumPy Implementation TODO

This document tracks the implementation status of NumPy functions in NumWasm, based on the official NumPy 2.x API.

**Goal:** Export all NumPy package/function names from this package.

## Legend

| Symbol | Meaning |
|--------|---------|
| ✅ | TypeScript API implemented |
| 🔧 | Has WASM/C backing |
| ⬜ | Not implemented |

---

## Summary by Module

| Module | Implemented | Total | Coverage |
|--------|-------------|-------|----------|
| **numpy (top-level)** | ~280 | ~350 | ~80% |
| **numpy.linalg** | 28 | 29 | 97% |
| **numpy.fft** | 18 | 18 | 100% |
| **numpy.random** | 45 | 56 | 80% |
| **numpy.polynomial** | 6 classes + functions | 6 classes | 100% (classes) |
| **numpy.ma** | Core class | Full module | Partial |
| **numpy.strings** | 40 | 40 | 100% |
| **numpy.rec** | Core functions | Full module | Partial |
| **numpy.testing** | 20 | 20 | 100% |

---

## Implementation Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                    TypeScript API Layer                          │
│  (packages/numwasm/src/ts/)                                     │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                    WASM Binary (numjs.wasm)                      │
│  Compiled from: packages/numwasm/src/wasm/                      │
│  Size: ~290KB (standalone, self-contained)                      │
└─────────────────────────────────────────────────────────────────┘
```

### WASM C Implementation (Self-Contained)

| Category | Lines | Status |
|----------|-------|--------|
| Core NDArray (`ndarray.c/h`) | ~1300 | ✅ Complete |
| DTypes (`dtype.c/h`) | ~350 | ✅ Complete |
| Broadcasting (`broadcast.c/h`) | ~200 | ✅ Complete |
| Indexing (`indexing.c/h`) | ~650 | ✅ Complete |
| Ufuncs (`ufunc*.c/h`) | ~1800 | ✅ Complete |
| Logic (`logic.c/h`) | ~500 | ✅ Complete |
| Sorting (`sorting.c/h`) | ~400 | ✅ Complete |
| Searching (`searching.c/h`) | ~250 | ✅ Complete |
| Statistics (`statistics.c/h`) | ~800 | ✅ Complete |
| Set Operations (`setops.c/h`) | ~1000 | ✅ Complete |
| BLAS (`blas.c/h`) | ~1400 | ✅ Standalone |
| LAPACK (`lapack.c/h`) | ~1800 | ✅ Core functions |
| Linalg (`linalg.c/h`) | ~1100 | ✅ Complete |
| FFT (`fft.c/h`) | ~650 | ✅ Standalone |
| Random (`random/*.c`) | ~2500 | ✅ Complete |
| **Total** | **~14,700** | **~95% core** |

---

## 1. numpy (Top-Level Namespace)

Based on NumPy's `__init__.py` exports from `_core` and `lib` modules.

### 1.1 Array Creation

| Function | TS | WASM | Notes |
|----------|:--:|:----:|-------|
| `array` | ✅ | 🔧 | |
| `zeros` | ✅ | 🔧 | Uses ndarray_create |
| `ones` | ✅ | 🔧 | Uses ndarray_full |
| `empty` | ✅ | 🔧 | |
| `full` | ✅ | 🔧 | |
| `arange` | ✅ | 🔧 | |
| `linspace` | ✅ | 🔧 | |
| `logspace` | ✅ | 🔧 | Uses linspace |
| `geomspace` | ✅ | 🔧 | Uses linspace |
| `eye` | ✅ | 🔧 | |
| `identity` | ✅ | 🔧 | Uses eye |
| `diag` | ✅ | 🔧 | |
| `diagflat` | ✅ | 🔧 | Uses diag |
| `tri` | ✅ | 🔧 | |
| `tril` | ✅ | 🔧 | |
| `triu` | ✅ | 🔧 | |
| `vander` | ✅ | 🔧 | |
| `zeros_like` | ✅ | 🔧 | Uses zeros |
| `ones_like` | ✅ | 🔧 | Uses ones |
| `empty_like` | ✅ | 🔧 | Uses empty |
| `full_like` | ✅ | 🔧 | Uses full |
| `fromfunction` | ✅ | ⬜ | Requires callback |
| `fromiter` | ✅ | ⬜ | Requires JS iterator |
| `fromstring` | ✅ | ⬜ | Parses text data |
| `frombuffer` | ✅ | 🔧 | |
| `fromfile` | ✅ | 🔧 | |
| `fromregex` | ✅ | ⬜ | Requires JS regex |
| `copy` | ✅ | 🔧 | |
| `asarray` | ✅ | 🔧 | |
| `asanyarray` | ✅ | 🔧 | Alias for asarray |
| `ascontiguousarray` | ✅ | 🔧 | |
| `asfortranarray` | ✅ | 🔧 | |
| `asarray_chkfinite` | ✅ | 🔧 | Uses asarray + isfinite |
| `require` | ✅ | ⬜ | Ensures array meets requirements |
| `astype` | ✅ | ⬜ | Top-level function wrapping method |

### 1.2 Array Manipulation

| Function | TS | WASM | Notes |
|----------|:--:|:----:|-------|
| `reshape` | ✅ | 🔧 | |
| `ravel` | ✅ | 🔧 | |
| `flatten` | ✅ | 🔧 | |
| `squeeze` | ✅ | 🔧 | |
| `expand_dims` | ✅ | 🔧 | |
| `transpose` | ✅ | 🔧 | |
| `permute_dims` | ✅ | 🔧 | NumPy 2.0 alias for transpose |
| `swapaxes` | ✅ | 🔧 | |
| `moveaxis` | ✅ | 🔧 | Uses transpose |
| `rollaxis` | ✅ | 🔧 | Uses transpose |
| `atleast_1d` | ✅ | 🔧 | Uses reshape/expandDims |
| `atleast_2d` | ✅ | 🔧 | Uses reshape/expandDims |
| `atleast_3d` | ✅ | 🔧 | Uses reshape/expandDims |
| `concatenate` | ✅ | 🔧 | |
| `concat` | ✅ | 🔧 | NumPy 2.0 alias |
| `stack` | ✅ | 🔧 | Uses expandDims + concatenate |
| `unstack` | ✅ | 🔧 | Uses slice |
| `vstack` | ✅ | 🔧 | Uses atleast_2d + concatenate |
| `hstack` | ✅ | 🔧 | Uses concatenate |
| `dstack` | ✅ | 🔧 | Uses atleast_3d + concatenate |
| `column_stack` | ✅ | 🔧 | Uses atleast_2d + concatenate |
| `block` | ✅ | 🔧 | Uses concatenate |
| `split` | ✅ | 🔧 | Uses slice + copy |
| `array_split` | ✅ | 🔧 | Uses slice + copy |
| `vsplit` | ✅ | 🔧 | Uses split |
| `hsplit` | ✅ | 🔧 | Uses split |
| `dsplit` | ✅ | 🔧 | Uses split |
| `tile` | ✅ | 🔧 | Uses concatenate |
| `repeat` | ✅ | 🔧 | Uses concatenate |
| `pad` | ✅ | 🔧 | Uses copy/fill |
| `flip` | ✅ | 🔧 | Uses slice with negative stride |
| `fliplr` | ✅ | 🔧 | Uses flip |
| `flipud` | ✅ | 🔧 | Uses flip |
| `roll` | ✅ | 🔧 | Uses slice + concatenate |
| `rot90` | ✅ | 🔧 | Uses flip + swapaxes |
| `resize` | ✅ | 🔧 | Uses concatenate + slice |
| `trim_zeros` | ✅ | 🔧 | Uses slice (WASM) |
| `insert` | ✅ | 🔧 | Uses slice + concatenate |
| `delete` | ✅ | 🔧 | Exported as `deleteArr`, uses slice + concatenate |
| `append` | ✅ | 🔧 | Uses concatenate |
| `copyto` | ✅ | 🔧 | Uses memcpy for contiguous arrays |

### 1.3 Mathematical Functions (Ufuncs)

#### Arithmetic

| Function | TS | WASM | Notes |
|----------|:--:|:----:|-------|
| `add` | ✅ | 🔧 | |
| `subtract` | ✅ | 🔧 | |
| `multiply` | ✅ | 🔧 | |
| `divide` | ✅ | 🔧 | |
| `true_divide` | ✅ | 🔧 | |
| `floor_divide` | ✅ | 🔧 | |
| `negative` | ✅ | 🔧 | |
| `positive` | ✅ | 🔧 | |
| `power` | ✅ | 🔧 | |
| `pow` | ✅ | 🔧 | Alias for power |
| `float_power` | ✅ | 🔧 | Same as power (uses float64) |
| `remainder` | ✅ | 🔧 | |
| `mod` | ✅ | 🔧 | |
| `fmod` | ✅ | 🔧 | |
| `divmod` | ✅ | 🔧 | |
| `absolute` | ✅ | 🔧 | |
| `abs` | ✅ | 🔧 | Alias |
| `fabs` | ✅ | 🔧 | Alias for absolute |
| `sign` | ✅ | 🔧 | |
| `reciprocal` | ✅ | 🔧 | |
| `sqrt` | ✅ | 🔧 | |
| `square` | ✅ | 🔧 | |
| `cbrt` | ✅ | 🔧 | |

#### Trigonometric

| Function | TS | WASM | Notes |
|----------|:--:|:----:|-------|
| `sin` | ✅ | 🔧 | |
| `cos` | ✅ | 🔧 | |
| `tan` | ✅ | 🔧 | |
| `arcsin` | ✅ | 🔧 | |
| `arccos` | ✅ | 🔧 | |
| `arctan` | ✅ | 🔧 | |
| `arctan2` | ✅ | 🔧 | |
| `hypot` | ✅ | 🔧 | |
| `sinh` | ✅ | 🔧 | |
| `cosh` | ✅ | 🔧 | |
| `tanh` | ✅ | 🔧 | |
| `arcsinh` | ✅ | 🔧 | |
| `arccosh` | ✅ | 🔧 | |
| `arctanh` | ✅ | 🔧 | |
| `degrees` | ✅ | 🔧 | |
| `radians` | ✅ | 🔧 | |
| `deg2rad` | ✅ | 🔧 | |
| `rad2deg` | ✅ | 🔧 | |
| `asin` | ✅ | 🔧 | NumPy 2.0 alias for arcsin |
| `acos` | ✅ | 🔧 | NumPy 2.0 alias for arccos |
| `atan` | ✅ | 🔧 | NumPy 2.0 alias for arctan |
| `atan2` | ✅ | 🔧 | NumPy 2.0 alias for arctan2 |
| `asinh` | ✅ | 🔧 | NumPy 2.0 alias for arcsinh |
| `acosh` | ✅ | 🔧 | NumPy 2.0 alias for arccosh |
| `atanh` | ✅ | 🔧 | NumPy 2.0 alias for arctanh |

#### Exponential & Logarithmic

| Function | TS | WASM | Notes |
|----------|:--:|:----:|-------|
| `exp` | ✅ | 🔧 | |
| `exp2` | ✅ | 🔧 | |
| `expm1` | ✅ | 🔧 | |
| `log` | ✅ | 🔧 | |
| `log2` | ✅ | 🔧 | |
| `log10` | ✅ | 🔧 | |
| `log1p` | ✅ | 🔧 | |
| `logaddexp` | ✅ | 🔧 | |
| `logaddexp2` | ✅ | 🔧 | |

#### Rounding

| Function | TS | WASM | Notes |
|----------|:--:|:----:|-------|
| `floor` | ✅ | 🔧 | |
| `ceil` | ✅ | 🔧 | |
| `trunc` | ✅ | 🔧 | |
| `rint` | ✅ | 🔧 | |
| `round` | ✅ | 🔧 | |
| `around` | ✅ | 🔧 | Alias for round |
| `fix` | ✅ | 🔧 | Alias for trunc |

#### Comparison

| Function | TS | WASM | Notes |
|----------|:--:|:----:|-------|
| `equal` | ✅ | 🔧 | |
| `not_equal` | ✅ | 🔧 | |
| `less` | ✅ | 🔧 | |
| `less_equal` | ✅ | 🔧 | |
| `greater` | ✅ | 🔧 | |
| `greater_equal` | ✅ | 🔧 | |
| `maximum` | ✅ | 🔧 | |
| `minimum` | ✅ | 🔧 | |
| `fmax` | ✅ | 🔧 | |
| `fmin` | ✅ | 🔧 | |

#### Logical

| Function | TS | WASM | Notes |
|----------|:--:|:----:|-------|
| `logical_and` | ✅ | 🔧 | |
| `logical_or` | ✅ | 🔧 | |
| `logical_xor` | ✅ | 🔧 | |
| `logical_not` | ✅ | 🔧 | |

#### Bitwise

| Function | TS | WASM | Notes |
|----------|:--:|:----:|-------|
| `bitwise_and` | ✅ | 🔧 | |
| `bitwise_or` | ✅ | 🔧 | |
| `bitwise_xor` | ✅ | 🔧 | |
| `bitwise_not` | ✅ | 🔧 | |
| `bitwise_invert` | ✅ | 🔧 | NumPy 2.0 alias for invert |
| `invert` | ✅ | 🔧 | |
| `left_shift` | ✅ | 🔧 | |
| `right_shift` | ✅ | 🔧 | |
| `bitwise_left_shift` | ✅ | 🔧 | NumPy 2.0 alias for left_shift |
| `bitwise_right_shift` | ✅ | 🔧 | NumPy 2.0 alias for right_shift |
| `bitwise_count` | ✅ | 🔧 | |

#### Special Math

| Function | TS | WASM | Notes |
|----------|:--:|:----:|-------|
| `copysign` | ✅ | 🔧 | |
| `signbit` | ✅ | 🔧 | |
| `heaviside` | ✅ | 🔧 | |
| `sinc` | ✅ | 🔧 | |
| `frexp` | ✅ | 🔧 | |
| `ldexp` | ✅ | 🔧 | |
| `nextafter` | ✅ | 🔧 | |
| `spacing` | ✅ | 🔧 | |
| `modf` | ✅ | 🔧 | |
| `gcd` | ✅ | 🔧 | |
| `lcm` | ✅ | 🔧 | |

#### Complex Numbers

| Function | TS | WASM | Notes |
|----------|:--:|:----:|-------|
| `conj` | ✅ | 🔧 | |
| `conjugate` | ✅ | 🔧 | |
| `real` | ✅ | ⬜ | |
| `imag` | ✅ | ⬜ | |
| `angle` | ✅ | ⬜ | |

### 1.4 Statistics & Reductions

| Function | TS | WASM | Notes |
|----------|:--:|:----:|-------|
| `sum` | ✅ | 🔧 | |
| `prod` | ✅ | 🔧 | |
| `mean` | ✅ | 🔧 | |
| `std` | ✅ | 🔧 | |
| `var` | ✅ | 🔧 | Exported as `var_` |
| `min` | ✅ | 🔧 | |
| `max` | ✅ | 🔧 | |
| `amin` | ✅ | 🔧 | Alias for min |
| `amax` | ✅ | 🔧 | Alias for max |
| `argmin` | ✅ | 🔧 | |
| `argmax` | ✅ | 🔧 | |
| `ptp` | ✅ | 🔧 | |
| `median` | ✅ | 🔧 | |
| `average` | ✅ | 🔧 | Uses multiply+sum WASM |
| `percentile` | ✅ | 🔧 | |
| `quantile` | ✅ | 🔧 | |
| `corrcoef` | ✅ | 🔧 | |
| `cov` | ✅ | 🔧 | |
| `histogram` | ✅ | ⬜ | |
| `histogram2d` | ✅ | ⬜ | |
| `histogramdd` | ✅ | ⬜ | |
| `histogram_bin_edges` | ✅ | ⬜ | |
| `bincount` | ✅ | ⬜ | |
| `digitize` | ✅ | ⬜ | |
| `count_nonzero` | ✅ | 🔧 | |

#### Cumulative Operations

| Function | TS | WASM | Notes |
|----------|:--:|:----:|-------|
| `cumsum` | ✅ | 🔧 | |
| `cumprod` | ✅ | 🔧 | |
| `cumulative_sum` | ✅ | 🔧 | NumPy 2.0 alias for cumsum |
| `cumulative_prod` | ✅ | 🔧 | NumPy 2.0 alias for cumprod |

#### NaN-handling Functions

| Function | TS | WASM | Notes |
|----------|:--:|:----:|-------|
| `nansum` | ✅ | 🔧 | Uses sum WASM |
| `nanprod` | ✅ | ⬜ | Pure TS loop |
| `nanmin` | ✅ | 🔧 | Uses min WASM |
| `nanmax` | ✅ | 🔧 | Uses max WASM |
| `nanargmin` | ✅ | 🔧 | Uses argmin WASM |
| `nanargmax` | ✅ | 🔧 | Uses argmax WASM |
| `nanmean` | ✅ | 🔧 | Uses sum WASM |
| `nanstd` | ✅ | 🔧 | Uses sum/ufuncs WASM |
| `nanvar` | ✅ | 🔧 | Uses sum/ufuncs WASM |
| `nanmedian` | ✅ | 🔧 | Uses median/extract WASM |
| `nanpercentile` | ✅ | 🔧 | Uses sort/extract WASM |
| `nanquantile` | ✅ | 🔧 | Uses sort/extract WASM |
| `nancumsum` | ✅ | 🔧 | Direct WASM |
| `nancumprod` | ✅ | 🔧 | Direct WASM |
| `nan_to_num` | ✅ | 🔧 | Uses where/isnan WASM |

### 1.5 Sorting & Searching

| Function | TS | WASM | Notes |
|----------|:--:|:----:|-------|
| `sort` | ✅ | 🔧 | |
| `argsort` | ✅ | 🔧 | |
| `lexsort` | ✅ | ⬜ | Pure TS stable sort |
| `partition` | ✅ | 🔧 | |
| `argpartition` | ✅ | 🔧 | |
| `sort_complex` | ✅ | ⬜ | |
| `searchsorted` | ✅ | 🔧 | |
| `extract` | ✅ | 🔧 | |

### 1.6 Logic & Predicates

| Function | TS | WASM | Notes |
|----------|:--:|:----:|-------|
| `all` | ✅ | 🔧 | |
| `any` | ✅ | 🔧 | |
| `isfinite` | ✅ | 🔧 | |
| `isinf` | ✅ | 🔧 | |
| `isnan` | ✅ | 🔧 | |
| `isnat` | ⬜ | ⬜ | datetime64 specific |
| `isneginf` | ✅ | 🔧 | |
| `isposinf` | ✅ | 🔧 | |
| `iscomplex` | ✅ | 🔧 | |
| `isreal` | ✅ | 🔧 | |
| `iscomplexobj` | ✅ | ⬜ | Pure TS dtype check |
| `isrealobj` | ✅ | ⬜ | Pure TS dtype check |
| `isfortran` | ✅ | ⬜ | Pure TS flags check |
| `isscalar` | ✅ | ⬜ | Pure TS typeof check |
| `isclose` | ✅ | 🔧 | |
| `allclose` | ✅ | 🔧 | |
| `array_equal` | ✅ | 🔧 | |
| `array_equiv` | ✅ | 🔧 | |
| `real_if_close` | ✅ | 🔧 | Uses absolute/isclose WASM |

### 1.7 Set Operations

| Function | TS | WASM | Notes |
|----------|:--:|:----:|-------|
| `unique` | ✅ | 🔧 | |
| `unique_all` | ✅ | 🔧 | Uses unique WASM |
| `unique_counts` | ✅ | 🔧 | Uses unique WASM |
| `unique_inverse` | ✅ | 🔧 | Uses unique WASM |
| `unique_values` | ✅ | 🔧 | |
| `union1d` | ✅ | 🔧 | |
| `intersect1d` | ✅ | 🔧 | |
| `setdiff1d` | ✅ | 🔧 | |
| `setxor1d` | ✅ | 🔧 | |
| `isin` | ✅ | 🔧 | |
| `in1d` | ✅ | 🔧 | Deprecated alias |
| `ediff1d` | ✅ | ⬜ | Pure TS loop |

### 1.8 Indexing

| Function | TS | WASM | Notes |
|----------|:--:|:----:|-------|
| `where` | ✅ | 🔧 | |
| `nonzero` | ✅ | 🔧 | |
| `flatnonzero` | ✅ | 🔧 | |
| `argwhere` | ✅ | 🔧 | Uses nonzero WASM |
| `take` | ✅ | 🔧 | |
| `take_along_axis` | ✅ | ⬜ | Pure TS loop |
| `put` | ✅ | 🔧 | |
| `put_along_axis` | ✅ | ⬜ | Pure TS loop |
| `putmask` | ✅ | ⬜ | Pure TS loop |
| `place` | ✅ | ⬜ | Uses putmask |
| `choose` | ✅ | 🔧 | |
| `compress` | ✅ | 🔧 | |
| `select` | ✅ | ⬜ | Pure TS loop |
| `diagonal` | ✅ | 🔧 | |
| `indices` | ✅ | ⬜ | Pure TS grid gen |
| `ix_` | ✅ | ⬜ | Pure TS reshape |
| `diag_indices` | ✅ | ⬜ | Uses arange |
| `diag_indices_from` | ✅ | ⬜ | Uses diag_indices |
| `tril_indices` | ✅ | ⬜ | Pure TS loop |
| `triu_indices` | ✅ | ⬜ | Pure TS loop |
| `tril_indices_from` | ✅ | ⬜ | Uses tril_indices |
| `triu_indices_from` | ✅ | ⬜ | Uses triu_indices |
| `mask_indices` | ✅ | ⬜ | Pure TS implementation |
| `fill_diagonal` | ✅ | ⬜ | In-place diagonal fill |
| `unravel_index` | ✅ | ⬜ | Pure TS arithmetic |
| `ravel_multi_index` | ✅ | ⬜ | Pure TS arithmetic |
| `meshgrid` | ✅ | ⬜ | Pure TS loop |

### 1.9 Special Functions (lib._function_base_impl)

| Function | TS | WASM | Notes |
|----------|:--:|:----:|-------|
| `clip` | ✅ | 🔧 | Uses maximum/minimum WASM |
| `diff` | ✅ | 🔧 | |
| `gradient` | ✅ | 🔧 | |
| `convolve` | ✅ | 🔧 | |
| `correlate` | ✅ | 🔧 | |
| `interp` | ✅ | ⬜ | Pure TS (binary search + linear interp) |
| `trapezoid` | ✅ | ⬜ | Pure TS |
| `unwrap` | ✅ | ⬜ | Pure TS |
| `i0` | ✅ | ⬜ | Bessel function |

### 1.10 Window Functions

| Function | TS | WASM | Notes |
|----------|:--:|:----:|-------|
| `bartlett` | ✅ | ⬜ | |
| `blackman` | ✅ | ⬜ | |
| `hamming` | ✅ | ⬜ | |
| `hanning` | ✅ | ⬜ | |
| `kaiser` | ✅ | ⬜ | |

### 1.11 Functional Programming

| Function | TS | WASM | Notes |
|----------|:--:|:----:|-------|
| `apply_along_axis` | ✅ | ⬜ | |
| `apply_over_axes` | ✅ | ⬜ | |
| `vectorize` | ✅ | ⬜ | |
| `frompyfunc` | ✅ | ⬜ | |
| `piecewise` | ✅ | ⬜ | |
| `iterable` | ✅ | ⬜ | Check if object is iterable |

### 1.12 Broadcasting

| Function | TS | WASM | Notes |
|----------|:--:|:----:|-------|
| `broadcast_to` | ✅ | ⬜ | |
| `broadcast_arrays` | ✅ | ⬜ | |
| `broadcast_shapes` | ✅ | ⬜ | |
| `broadcast` | ✅ | 🔧 | Iterator class for broadcasting |

### 1.13 I/O Functions

| Function | TS | WASM | Notes |
|----------|:--:|:----:|-------|
| `save` | ✅ | ⬜ | NPY format |
| `load` | ✅ | ⬜ | NPY format |
| `savez` | ✅ | ⬜ | NPZ format (ZIP archive) |
| `savez_compressed` | ✅ | ⬜ | NPZ format (DEFLATE) |
| `savetxt` | ✅ | ⬜ | |
| `loadtxt` | ✅ | ⬜ | |
| `genfromtxt` | ✅ | ⬜ | |
| `packbits` | ✅ | ⬜ | Bit packing |
| `unpackbits` | ✅ | ⬜ | Bit unpacking |

### 1.14 String & Base Representation

| Function | TS | WASM | Notes |
|----------|:--:|:----:|-------|
| `binary_repr` | ✅ | ⬜ | |
| `base_repr` | ✅ | ⬜ | |
| `array2string` | ✅ | ⬜ | |
| `array_repr` | ✅ | ⬜ | |
| `array_str` | ✅ | ⬜ | |
| `format_float_positional` | ✅ | ⬜ | |
| `format_float_scientific` | ✅ | ⬜ | |
| `set_printoptions` | ✅ | ⬜ | |
| `get_printoptions` | ✅ | ⬜ | |
| `printoptions` | ⬜ | ⬜ | Context manager |

### 1.15 Type Information

| Function | TS | WASM | Notes |
|----------|:--:|:----:|-------|
| `dtype` | ✅ | ⬜ | Class |
| `finfo` | ✅ | ⬜ | |
| `iinfo` | ✅ | ⬜ | |
| `can_cast` | ✅ | ⬜ | Type casting check |
| `result_type` | ✅ | ⬜ | Type promotion |
| `promote_types` | ✅ | ⬜ | |
| `min_scalar_type` | ✅ | ⬜ | Returns minimum dtype for scalar |
| `issubdtype` | ✅ | ⬜ | Dtype hierarchy check |
| `isdtype` | ✅ | ⬜ | NumPy 2.0 dtype check |
| `common_type` | ✅ | ⬜ | Type promotion |
| `mintypecode` | ⬜ | ⬜ | Low priority |
| `typename` | ⬜ | ⬜ | Low priority |

### 1.16 Constants

| Constant | TS | Notes |
|----------|:--:|-------|
| `e` | ✅ | Euler's number |
| `pi` | ✅ | |
| `euler_gamma` | ✅ | Euler-Mascheroni constant |
| `inf` | ✅ | |
| `nan` | ✅ | |
| `newaxis` | ✅ | Alias for `None` in indexing |
| `PINF` | ✅ | |
| `NINF` | ✅ | |
| `PZERO` | ✅ | |
| `NZERO` | ✅ | |

### 1.17 Linear Algebra (Top-Level)

| Function | TS | WASM | Notes |
|----------|:--:|:----:|-------|
| `dot` | ✅ | 🔧 | |
| `vdot` | ✅ | 🔧 | |
| `inner` | ✅ | 🔧 | |
| `outer` | ✅ | 🔧 | Uses matmul WASM |
| `matmul` | ✅ | 🔧 | |
| `tensordot` | ✅ | 🔧 | Uses matmul WASM |
| `einsum` | ✅ | 🔧 | Uses matmul/sum WASM |
| `einsum_path` | ✅ | ⬜ | Pure TS planning |
| `kron` | ✅ | 🔧 | Uses ufunc_multiply WASM |
| `cross` | ✅ | ⬜ | Pure TS loop |
| `trace` | ✅ | 🔧 | Uses diagonal/sum WASM |
| `vecdot` | ✅ | 🔧 | NumPy 2.0, uses multiply+sum WASM |
| `matvec` | ✅ | 🔧 | NumPy 2.0, uses matmul WASM |
| `vecmat` | ✅ | 🔧 | NumPy 2.0, uses matmul WASM |
| `matrix_transpose` | ✅ | ⬜ | NumPy 2.0, uses swapaxes |

### 1.18 Polynomial (Top-Level Legacy)

| Function | TS | WASM | Notes |
|----------|:--:|:----:|-------|
| `poly` | ✅ | ⬜ | Legacy function (descending powers) |
| `poly1d` | ✅ | ⬜ | Class for polynomial operations |
| `polyadd` | ✅ | ⬜ | |
| `polyder` | ✅ | ⬜ | |
| `polydiv` | ✅ | ⬜ | |
| `polyfit` | ✅ | ⬜ | |
| `polyint` | ✅ | ⬜ | |
| `polymul` | ✅ | ⬜ | |
| `polysub` | ✅ | ⬜ | |
| `polyval` | ✅ | ⬜ | |
| `roots` | ✅ | ⬜ | Legacy function (descending powers) |

---

## 2. numpy.linalg

| Function | TS | WASM | Notes |
|----------|:--:|:----:|-------|
| `LinAlgError` | ✅ | ⬜ | Exception class |
| **Matrix Products** | | | |
| `cross` | ✅ | ⬜ | Pure TS loop |
| `multi_dot` | ✅ | 🔧 | Uses dot WASM |
| `matrix_power` | ✅ | 🔧 | Uses inv/matmul WASM |
| `tensordot` | ✅ | 🔧 | Uses dot WASM |
| `matmul` | ✅ | 🔧 | |
| `outer` | ✅ | 🔧 | Uses matmul WASM |
| **Decompositions** | | | |
| `cholesky` | ✅ | 🔧 | |
| `qr` | ✅ | 🔧 | |
| `svd` | ✅ | 🔧 | |
| `svdvals` | ✅ | 🔧 | Uses svd WASM |
| **Eigenvalues** | | | |
| `eig` | ✅ | 🔧 | |
| `eigh` | ✅ | 🔧 | Uses eig WASM |
| `eigvals` | ✅ | 🔧 | Uses eig WASM |
| `eigvalsh` | ✅ | 🔧 | Uses eigh WASM |
| **Norms** | | | |
| `norm` | ✅ | 🔧 | |
| `matrix_norm` | ✅ | 🔧 | Uses norm WASM |
| `vector_norm` | ✅ | ⬜ | Pure TS loop |
| `cond` | ✅ | 🔧 | Uses svdvals/norm/inv WASM |
| `det` | ✅ | 🔧 | |
| `matrix_rank` | ✅ | 🔧 | Uses svdvals WASM |
| `slogdet` | ✅ | 🔧 | Uses det WASM |
| `trace` | ✅ | 🔧 | Uses diagonal/sum WASM |
| **Solving** | | | |
| `solve` | ✅ | 🔧 | |
| `tensorsolve` | ✅ | 🔧 | Uses solve WASM |
| `lstsq` | ✅ | 🔧 | Uses pinv/matmul WASM |
| `inv` | ✅ | 🔧 | |
| `pinv` | ✅ | 🔧 | Uses svd WASM |
| `tensorinv` | ✅ | 🔧 | Uses inv WASM |
| **Other** | | | |
| `diagonal` | ✅ | 🔧 | |
| `matrix_transpose` | ✅ | ⬜ | NumPy 2.0, uses swapaxes |

---

## 3. numpy.fft - Complete ✅

| Function | TS | WASM | Notes |
|----------|:--:|:----:|-------|
| **Standard FFTs** | | | |
| `fft` | ✅ | 🔧 | |
| `ifft` | ✅ | 🔧 | |
| `fft2` | ✅ | 🔧 | |
| `ifft2` | ✅ | 🔧 | |
| `fftn` | ✅ | 🔧 | |
| `ifftn` | ✅ | 🔧 | |
| **Real FFTs** | | | |
| `rfft` | ✅ | 🔧 | |
| `irfft` | ✅ | 🔧 | |
| `rfft2` | ✅ | 🔧 | |
| `irfft2` | ✅ | 🔧 | |
| `rfftn` | ✅ | 🔧 | |
| `irfftn` | ✅ | 🔧 | |
| **Hermitian FFTs** | | | |
| `hfft` | ✅ | 🔧 | |
| `ihfft` | ✅ | 🔧 | |
| **Helpers** | | | |
| `fftfreq` | ✅ | ⬜ | |
| `rfftfreq` | ✅ | ⬜ | |
| `fftshift` | ✅ | ⬜ | |
| `ifftshift` | ✅ | ⬜ | |

---

## 4. numpy.random

### Classes & Utilities

| Item | TS | Notes |
|------|:--:|-------|
| `Generator` | ✅ | Main random class |
| `RandomState` | ⬜ | Legacy, low priority |
| `SeedSequence` | ✅ | |
| `BitGenerator` | ✅ | Base class |
| `default_rng` | ✅ | Factory function |

### BitGenerators

| BitGenerator | TS | Notes |
|--------------|:--:|-------|
| `MT19937` | ✅ | Mersenne Twister |
| `PCG64` | ✅ | Default |
| `PCG64DXSM` | ✅ | Pure TS 128-bit arithmetic |
| `Philox` | ✅ | |
| `SFC64` | ✅ | |

### Utility Functions

| Function | TS | WASM | Notes |
|----------|:--:|:----:|-------|
| `seed` | ✅ | ⬜ | |
| `random` | ✅ | 🔧 | |
| `randn` | ✅ | 🔧 | |
| `randint` | ✅ | 🔧 | |
| `rand` | ✅ | ⬜ | Legacy (deprecated) |
| `ranf` | ✅ | ⬜ | Legacy (deprecated) |
| `random_sample` | ✅ | ⬜ | Legacy (deprecated) |
| `random_integers` | ✅ | ⬜ | Deprecated |
| `sample` | ✅ | ⬜ | Legacy (deprecated) |
| `bytes` | ✅ | ⬜ | Top-level function |
| `choice` | ✅ | ⬜ | |
| `shuffle` | ✅ | ⬜ | |
| `permutation` | ✅ | ⬜ | |
| `get_state` | ⬜ | ⬜ | Legacy |
| `set_state` | ⬜ | ⬜ | Legacy |

### Distributions (Generator methods)

| Distribution | TS | WASM | Notes |
|--------------|:--:|:----:|-------|
| `beta` | ✅ | 🔧 | Generator method |
| `binomial` | ✅ | 🔧 | Generator method |
| `chisquare` | ✅ | 🔧 | Generator method |
| `dirichlet` | ✅ | ⬜ | Generator method |
| `exponential` | ✅ | 🔧 | Generator method |
| `f` | ✅ | 🔧 | Generator method |
| `gamma` | ✅ | 🔧 | Generator method |
| `geometric` | ✅ | 🔧 | Generator method |
| `gumbel` | ✅ | 🔧 | Generator method |
| `hypergeometric` | ✅ | 🔧 | Generator method |
| `laplace` | ✅ | 🔧 | Generator method |
| `logistic` | ✅ | 🔧 | Generator method |
| `lognormal` | ✅ | 🔧 | Generator method |
| `logseries` | ✅ | 🔧 | Generator method |
| `multinomial` | ✅ | ⬜ | Generator method |
| `multivariate_normal` | ✅ | ⬜ | Generator method |
| `negative_binomial` | ✅ | 🔧 | Generator method |
| `noncentral_chisquare` | ✅ | 🔧 | Generator method |
| `noncentral_f` | ✅ | 🔧 | Generator method |
| `normal` | ✅ | 🔧 | Generator method |
| `pareto` | ✅ | 🔧 | Generator method |
| `poisson` | ✅ | 🔧 | Generator method |
| `power` | ✅ | 🔧 | Generator method |
| `rayleigh` | ✅ | 🔧 | Generator method |
| `standard_cauchy` | ✅ | 🔧 | Generator method |
| `standard_exponential` | ✅ | 🔧 | Generator method |
| `standard_gamma` | ✅ | 🔧 | Generator method |
| `standard_normal` | ✅ | 🔧 | Generator method |
| `standard_t` | ✅ | 🔧 | Generator method |
| `triangular` | ✅ | 🔧 | Generator method |
| `uniform` | ✅ | 🔧 | Generator method |
| `vonmises` | ✅ | 🔧 | Generator method |
| `wald` | ✅ | 🔧 | Generator method |
| `weibull` | ✅ | 🔧 | Generator method |
| `zipf` | ✅ | 🔧 | Generator method |

---

## 5. numpy.polynomial

### Classes - All Implemented ✅

| Class | TS | Notes |
|-------|:--:|-------|
| `Polynomial` | ✅ | Power series |
| `Chebyshev` | ✅ | Chebyshev series |
| `Legendre` | ✅ | Legendre series |
| `Hermite` | ✅ | Hermite (physicist's) |
| `HermiteE` | ✅ | Hermite (probabilist's) |
| `Laguerre` | ✅ | Laguerre series |

### Utility Functions

| Function | TS | Notes |
|----------|:--:|-------|
| `set_default_printstyle` | ⬜ | Low priority |

Each polynomial class has associated functions (e.g., `polyval`, `chebval`, `legval`, etc.) which are implemented.

---

## 6. numpy.ma (Masked Arrays)

| Item | TS | Notes |
|------|:--:|-------|
| `MaskedArray` | ✅ | Core class |
| `ma` namespace | ✅ | Module object |
| Core operations | ✅ | Via delegation |

Most mathematical operations are implemented via delegation to the main numpy functions. Full implementation of all `ma` functions is lower priority.

---

## 7. numpy.strings

| Function | TS | Notes |
|----------|:--:|-------|
| **Comparison** | | |
| `equal` | ✅ | |
| `not_equal` | ✅ | |
| `less` | ✅ | |
| `less_equal` | ✅ | |
| `greater` | ✅ | |
| `greater_equal` | ✅ | |
| **Property Testing** | | |
| `isalpha` | ✅ | |
| `isdigit` | ✅ | |
| `isalnum` | ✅ | |
| `isspace` | ✅ | |
| `islower` | ✅ | |
| `isupper` | ✅ | |
| `istitle` | ✅ | |
| `isdecimal` | ✅ | |
| `isnumeric` | ✅ | |
| `str_len` | ✅ | |
| **Search** | | |
| `find` | ✅ | |
| `rfind` | ✅ | |
| `index` | ✅ | |
| `rindex` | ✅ | |
| `count` | ✅ | |
| `startswith` | ✅ | |
| `endswith` | ✅ | |
| **Manipulation** | | |
| `add` | ✅ | |
| `multiply` | ✅ | |
| `upper` | ✅ | |
| `lower` | ✅ | |
| `swapcase` | ✅ | |
| `capitalize` | ✅ | |
| `title` | ✅ | |
| `strip` | ✅ | |
| `lstrip` | ✅ | |
| `rstrip` | ✅ | |
| `replace` | ✅ | |
| `center` | ✅ | |
| `ljust` | ✅ | |
| `rjust` | ✅ | |
| `zfill` | ✅ | |
| `expandtabs` | ✅ | |
| `partition` | ✅ | |
| `rpartition` | ✅ | |
| `encode` | ✅ | |
| `decode` | ✅ | |
| `mod` | ✅ | Printf-style formatting |
| `translate` | ✅ | Character translation |
| `slice` | ✅ | NumPy 2.0 string slicing |

---

## 8. numpy.rec (Record Arrays)

| Item | TS | Notes |
|------|:--:|-------|
| `recarray` | ✅ | Core class |
| `record` | ✅ | Record type |
| `fromarrays` | ✅ | |
| `fromrecords` | ✅ | |
| `fromstring` | ✅ | |
| `fromfile` | ✅ | |
| `array` | ✅ | |
| `format_parser` | ✅ | |
| `find_duplicate` | ✅ | |

---

## 9. numpy.testing

| Item | TS | Notes |
|------|:--:|-------|
| `SkipTest` | ✅ | |
| `KnownFailureException` | ✅ | |
| `AssertionError` | ✅ | |
| `assert_` | ✅ | Basic assertion |
| `assert_equal` | ✅ | Scalar/array equality |
| `assert_almost_equal` | ✅ | Decimal places tolerance |
| `assert_approx_equal` | ✅ | Significant figures tolerance |
| `assert_array_equal` | ✅ | Array element-wise equality |
| `assert_array_almost_equal` | ✅ | Array decimal tolerance |
| `assert_allclose` | ✅ | rtol/atol tolerance |
| `assert_array_less` | ✅ | Element-wise less than |
| `assert_array_max_ulp` | ✅ | ULP difference check |
| `assert_array_compare` | ✅ | Custom comparison |
| `assert_raises` | ✅ | Exception assertion |
| `assert_raises_regex` | ✅ | Exception + message match |
| `assert_warns` | ✅ | Warning assertion |
| `assert_string_equal` | ✅ | String equality with diff |
| `assert_no_warnings` | ✅ | No-warning assertion |
| `measure` | ✅ | Timing utility |
| `print_assert_equal` | ✅ | Debug printing |

---

## Not Planned for NumWasm

These items are Python-specific, deprecated, or low priority:

### Python/Environment Specific
- `datetime64`, `timedelta64` - Use JavaScript Date
- `busday_count`, `busday_offset`, `busdaycalendar`, `is_busday` - Business day functions
- `errstate`, `geterr`, `seterr`, `geterrcall`, `seterrcall` - Error state management
- `getbufsize`, `setbufsize` - Buffer management
- `show_config`, `show_runtime`, `get_include` - Environment introspection
- `may_share_memory`, `shares_memory` - Memory introspection
- `from_dlpack` - DLPack interop
- `nested_iters`, `nditer` class - Advanced iteration
- `flatiter` class - Iterator class
- `ufunc` class - Cannot create custom ufuncs in TypeScript
- `memmap` - Memory-mapped files (use browser File API)

### Deprecated/Legacy
- `matrix` class - Deprecated in NumPy
- All `*_` scalar types (`int_`, `float_`, `bool_`, etc.) - Use TypeScript types
- `random.RandomState` - Legacy API
- `random.rand`, `random.ranf`, `random.sample` - Legacy aliases

### Low Priority
- `packbits`, `unpackbits` - Bit packing
- Most `numpy.testing` functions - Testing utilities
- `ctypeslib` module - ctypes interop
- `f2py` module - Fortran-to-Python (not applicable)

---

## Priority TODO Items

### Priority 1: High-Impact Missing Functions

| Function | Category | Notes |
|----------|----------|-------|
| ~~`prod`~~ | Statistics | ✅ DONE |
| ~~`clip`~~ | Math | ✅ DONE |
| ~~`diff`~~ | Math | ✅ DONE |
| ~~`gradient`~~ | Math | ✅ DONE |
| ~~`convolve`~~ | Math | ✅ DONE |
| ~~`correlate`~~ | Math | ✅ DONE |
| ~~`cov`~~ | Statistics | ✅ DONE |
| ~~`corrcoef`~~ | Statistics | ✅ DONE |
| ~~`average`~~ | Statistics | ✅ DONE |
| ~~`ptp`~~ | Statistics | ✅ DONE |
| ~~`amin`~~ | Statistics | ✅ DONE (alias) |
| ~~`amax`~~ | Statistics | ✅ DONE (alias) |
| ~~`percentile`~~ | Statistics | ✅ DONE |
| ~~`quantile`~~ | Statistics | ✅ DONE |

### Priority 2: Export Functions with WASM Backing

These have C implementations but need TypeScript exports:
- ~~`reshape`, `ravel`, `squeeze`, `expand_dims`, `transpose`, `swapaxes`~~ ✅ Already exported
- ~~`copy`, `ascontiguousarray`, `asfortranarray`~~ ✅ Already exported
- ~~`percentile`, `quantile`~~ ✅ DONE
- ~~`tri`, `tril`, `triu`~~ ✅ Already exported

### Priority 3: Complete numpy.random

Most distributions have WASM backing but need TypeScript API:
- All 30+ distributions listed in section 4
- ~~`choice`, `shuffle`, `permutation`~~ ✅ DONE
- ~~`dirichlet`, `multinomial`, `multivariate_normal`~~ ✅ DONE
- ~~Legacy functions: `rand`, `ranf`, `random_sample`, `sample`, `random_integers`~~ ✅ DONE
- ~~`PCG64DXSM` BitGenerator~~ ✅ DONE

### Priority 4: NumPy 2.0 Compatibility

Add aliases and new functions from NumPy 2.0:
- ~~`cumulative_sum`, `cumulative_prod`~~ ✅ DONE
- ~~`matrix_transpose`, `vecdot`, `matvec`, `vecmat`~~ ✅ DONE
- ~~`concat` (alias for `concatenate`)~~ ✅ Already exported
- ~~Various trigonometric aliases (`asin`, `acos`, etc.)~~ ✅ Already exported

---

## Action Items Summary

### Immediate (High Impact)
1. ~~Implement `prod`, `clip`, `average`, `ptp`~~ ✅ DONE
2. ~~Implement `diff`, `gradient`, `convolve`, `correlate`~~ ✅ DONE
3. ~~Implement `cov`, `corrcoef`~~ ✅ DONE
4. ~~Export `reshape`, `transpose`, etc. as top-level functions~~ ✅ Already exported

### Short-term
5. Complete numpy.random TypeScript API for all distributions
6. ~~Add `choice`, `shuffle`, `permutation`~~ ✅ DONE
7. ~~Add `diagflat`, `tril`, `triu`, `vander`~~ ✅ Already exported

### Medium-term
8. ~~Add NumPy 2.0 aliases and new functions~~ ✅ DONE (cumulative_sum/prod, matrix_transpose, vecdot, matvec, vecmat)
9. Complete numpy.ma module
10. ~~Add missing numpy.testing assertions~~ ✅ DONE (all assertions implemented)
11. ~~Add numpy.strings mod/translate/slice~~ ✅ DONE
12. ~~Add savez/savez_compressed I/O~~ ✅ DONE
13. ~~Add packbits/unpackbits~~ ✅ DONE
14. ~~Add fromstring~~ ✅ DONE
