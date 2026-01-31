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
| **numpy.random** | 15 | 56 | 27% |
| **numpy.polynomial** | 6 classes + functions | 6 classes | 100% (classes) |
| **numpy.ma** | Core class | Full module | Partial |
| **numpy.strings** | 38 | 40 | 95% |
| **numpy.rec** | Core functions | Full module | Partial |
| **numpy.testing** | 3 | 30+ | 10% |

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
| `fromstring` | ⬜ | ⬜ | Low priority |
| `frombuffer` | ✅ | 🔧 | |
| `fromfile` | ✅ | 🔧 | |
| `fromregex` | ✅ | ⬜ | Requires JS regex |
| `copy` | ✅ | 🔧 | |
| `asarray` | ✅ | 🔧 | |
| `asanyarray` | ✅ | 🔧 | Alias for asarray |
| `ascontiguousarray` | ✅ | 🔧 | |
| `asfortranarray` | ✅ | 🔧 | |
| `asarray_chkfinite` | ✅ | 🔧 | Uses asarray + isfinite |
| `require` | ⬜ | ⬜ | Low priority |
| `astype` | ⬜ | ⬜ | Method exists on NDArray |

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
| `bitwise_invert` | ⬜ | ⬜ | NumPy 2.0 alias |
| `invert` | ✅ | 🔧 | |
| `left_shift` | ✅ | 🔧 | |
| `right_shift` | ✅ | 🔧 | |
| `bitwise_left_shift` | ⬜ | ⬜ | NumPy 2.0 alias |
| `bitwise_right_shift` | ⬜ | ⬜ | NumPy 2.0 alias |
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
| `prod` | ⬜ | ⬜ | **HIGH PRIORITY** |
| `mean` | ✅ | 🔧 | |
| `std` | ✅ | 🔧 | |
| `var` | ✅ | 🔧 | Exported as `var_` |
| `min` | ✅ | 🔧 | |
| `max` | ✅ | 🔧 | |
| `amin` | ⬜ | ⬜ | Alias for min |
| `amax` | ⬜ | ⬜ | Alias for max |
| `argmin` | ✅ | 🔧 | |
| `argmax` | ✅ | 🔧 | |
| `ptp` | ⬜ | ⬜ | **HIGH PRIORITY** |
| `median` | ✅ | 🔧 | |
| `average` | ⬜ | ⬜ | **HIGH PRIORITY** |
| `percentile` | ⬜ | 🔧 | Has WASM |
| `quantile` | ⬜ | 🔧 | Has WASM |
| `corrcoef` | ⬜ | ⬜ | **HIGH PRIORITY** |
| `cov` | ⬜ | ⬜ | **HIGH PRIORITY** |
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
| `cumulative_sum` | ⬜ | ⬜ | NumPy 2.0 |
| `cumulative_prod` | ⬜ | ⬜ | NumPy 2.0 |

#### NaN-handling Functions

| Function | TS | WASM | Notes |
|----------|:--:|:----:|-------|
| `nansum` | ✅ | 🔧 | |
| `nanprod` | ✅ | ⬜ | |
| `nanmin` | ✅ | ⬜ | |
| `nanmax` | ✅ | ⬜ | |
| `nanargmin` | ✅ | ⬜ | |
| `nanargmax` | ✅ | ⬜ | |
| `nanmean` | ✅ | 🔧 | |
| `nanstd` | ✅ | 🔧 | |
| `nanvar` | ✅ | 🔧 | |
| `nanmedian` | ✅ | ⬜ | |
| `nanpercentile` | ✅ | ⬜ | |
| `nanquantile` | ✅ | ⬜ | |
| `nancumsum` | ✅ | 🔧 | |
| `nancumprod` | ✅ | 🔧 | |
| `nan_to_num` | ✅ | ⬜ | |

### 1.5 Sorting & Searching

| Function | TS | WASM | Notes |
|----------|:--:|:----:|-------|
| `sort` | ✅ | 🔧 | |
| `argsort` | ✅ | 🔧 | |
| `lexsort` | ⬜ | ⬜ | TODO |
| `partition` | ✅ | 🔧 | |
| `argpartition` | ✅ | 🔧 | |
| `sort_complex` | ✅ | ⬜ | |
| `searchsorted` | ✅ | 🔧 | |
| `extract` | ✅ | ⬜ | |

### 1.6 Logic & Predicates

| Function | TS | WASM | Notes |
|----------|:--:|:----:|-------|
| `all` | ✅ | 🔧 | |
| `any` | ✅ | 🔧 | |
| `isfinite` | ✅ | 🔧 | |
| `isinf` | ✅ | 🔧 | |
| `isnan` | ✅ | 🔧 | |
| `isnat` | ⬜ | ⬜ | datetime64 specific |
| `isneginf` | ✅ | ⬜ | |
| `isposinf` | ✅ | ⬜ | |
| `iscomplex` | ✅ | ⬜ | |
| `isreal` | ✅ | ⬜ | |
| `iscomplexobj` | ✅ | ⬜ | |
| `isrealobj` | ✅ | ⬜ | |
| `isfortran` | ✅ | ⬜ | |
| `isscalar` | ✅ | ⬜ | |
| `isclose` | ✅ | ⬜ | |
| `allclose` | ✅ | ⬜ | |
| `array_equal` | ✅ | ⬜ | |
| `array_equiv` | ✅ | ⬜ | |
| `real_if_close` | ✅ | ⬜ | |

### 1.7 Set Operations

| Function | TS | WASM | Notes |
|----------|:--:|:----:|-------|
| `unique` | ✅ | 🔧 | |
| `unique_all` | ✅ | ⬜ | NumPy 2.0 |
| `unique_counts` | ✅ | ⬜ | NumPy 2.0 |
| `unique_inverse` | ✅ | ⬜ | NumPy 2.0 |
| `unique_values` | ✅ | ⬜ | NumPy 2.0 |
| `union1d` | ✅ | ⬜ | |
| `intersect1d` | ✅ | ⬜ | |
| `setdiff1d` | ✅ | ⬜ | |
| `setxor1d` | ✅ | ⬜ | |
| `isin` | ✅ | ⬜ | |
| `in1d` | ✅ | ⬜ | Deprecated alias |
| `ediff1d` | ✅ | ⬜ | |

### 1.8 Indexing

| Function | TS | WASM | Notes |
|----------|:--:|:----:|-------|
| `where` | ✅ | 🔧 | |
| `nonzero` | ✅ | 🔧 | |
| `flatnonzero` | ✅ | ⬜ | |
| `argwhere` | ✅ | ⬜ | |
| `take` | ✅ | 🔧 | |
| `take_along_axis` | ✅ | ⬜ | |
| `put` | ✅ | ⬜ | |
| `put_along_axis` | ✅ | ⬜ | |
| `putmask` | ✅ | ⬜ | |
| `place` | ✅ | ⬜ | |
| `choose` | ✅ | ⬜ | |
| `compress` | ✅ | ⬜ | |
| `select` | ✅ | ⬜ | |
| `diagonal` | ✅ | 🔧 | |
| `indices` | ✅ | ⬜ | |
| `ix_` | ✅ | ⬜ | |
| `diag_indices` | ✅ | ⬜ | |
| `diag_indices_from` | ⬜ | ⬜ | TODO |
| `tril_indices` | ✅ | ⬜ | |
| `triu_indices` | ✅ | ⬜ | |
| `tril_indices_from` | ⬜ | ⬜ | TODO |
| `triu_indices_from` | ⬜ | ⬜ | TODO |
| `mask_indices` | ⬜ | ⬜ | TODO |
| `fill_diagonal` | ⬜ | ⬜ | TODO |
| `unravel_index` | ✅ | ⬜ | |
| `ravel_multi_index` | ✅ | ⬜ | |
| `meshgrid` | ✅ | ⬜ | |

### 1.9 Special Functions (lib._function_base_impl)

| Function | TS | WASM | Notes |
|----------|:--:|:----:|-------|
| `clip` | ⬜ | ⬜ | **HIGH PRIORITY** |
| `diff` | ⬜ | ⬜ | **HIGH PRIORITY** |
| `gradient` | ⬜ | ⬜ | **HIGH PRIORITY** |
| `convolve` | ⬜ | ⬜ | **HIGH PRIORITY** |
| `correlate` | ⬜ | ⬜ | **HIGH PRIORITY** |
| `interp` | ⬜ | ⬜ | TODO |
| `trapezoid` | ⬜ | ⬜ | TODO |
| `unwrap` | ⬜ | ⬜ | TODO |
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
| `iterable` | ⬜ | ⬜ | Low priority |

### 1.12 Broadcasting

| Function | TS | WASM | Notes |
|----------|:--:|:----:|-------|
| `broadcast_to` | ✅ | ⬜ | |
| `broadcast_arrays` | ✅ | ⬜ | |
| `broadcast_shapes` | ✅ | ⬜ | |
| `broadcast` | ⬜ | ⬜ | Class, low priority |

### 1.13 I/O Functions

| Function | TS | WASM | Notes |
|----------|:--:|:----:|-------|
| `save` | ✅ | ⬜ | NPY format |
| `load` | ✅ | ⬜ | NPY format |
| `savez` | ⬜ | ⬜ | TODO |
| `savez_compressed` | ⬜ | ⬜ | TODO |
| `savetxt` | ✅ | ⬜ | |
| `loadtxt` | ✅ | ⬜ | |
| `genfromtxt` | ✅ | ⬜ | |
| `packbits` | ⬜ | ⬜ | Low priority |
| `unpackbits` | ⬜ | ⬜ | Low priority |

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
| `can_cast` | ⬜ | ⬜ | TODO |
| `result_type` | ⬜ | ⬜ | TODO |
| `promote_types` | ✅ | ⬜ | |
| `min_scalar_type` | ⬜ | ⬜ | Low priority |
| `issubdtype` | ⬜ | ⬜ | TODO |
| `isdtype` | ⬜ | ⬜ | NumPy 2.0 |
| `common_type` | ⬜ | ⬜ | TODO |
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
| `outer` | ✅ | ⬜ | |
| `matmul` | ✅ | 🔧 | |
| `tensordot` | ✅ | ⬜ | |
| `einsum` | ✅ | ⬜ | |
| `einsum_path` | ✅ | ⬜ | |
| `kron` | ✅ | ⬜ | |
| `cross` | ✅ | ⬜ | |
| `trace` | ✅ | ⬜ | |
| `vecdot` | ⬜ | ⬜ | NumPy 2.0 |
| `matvec` | ⬜ | ⬜ | NumPy 2.0 |
| `vecmat` | ⬜ | ⬜ | NumPy 2.0 |
| `matrix_transpose` | ⬜ | ⬜ | NumPy 2.0 |

### 1.18 Polynomial (Top-Level Legacy)

| Function | TS | WASM | Notes |
|----------|:--:|:----:|-------|
| `poly` | ⬜ | ⬜ | TODO |
| `poly1d` | ⬜ | ⬜ | TODO |
| `polyadd` | ✅ | ⬜ | |
| `polyder` | ✅ | ⬜ | |
| `polydiv` | ✅ | ⬜ | |
| `polyfit` | ✅ | ⬜ | |
| `polyint` | ✅ | ⬜ | |
| `polymul` | ✅ | ⬜ | |
| `polysub` | ✅ | ⬜ | |
| `polyval` | ✅ | ⬜ | |
| `roots` | ⬜ | ⬜ | TODO |

---

## 2. numpy.linalg

| Function | TS | WASM | Notes |
|----------|:--:|:----:|-------|
| `LinAlgError` | ✅ | ⬜ | Exception class |
| **Matrix Products** | | | |
| `cross` | ✅ | ⬜ | |
| `multi_dot` | ✅ | ⬜ | |
| `matrix_power` | ✅ | ⬜ | |
| `tensordot` | ✅ | ⬜ | |
| `matmul` | ✅ | 🔧 | |
| `outer` | ✅ | ⬜ | |
| **Decompositions** | | | |
| `cholesky` | ✅ | 🔧 | |
| `qr` | ✅ | 🔧 | |
| `svd` | ✅ | 🔧 | |
| `svdvals` | ✅ | ⬜ | |
| **Eigenvalues** | | | |
| `eig` | ✅ | 🔧 | |
| `eigh` | ✅ | ⬜ | |
| `eigvals` | ✅ | ⬜ | |
| `eigvalsh` | ✅ | ⬜ | |
| **Norms** | | | |
| `norm` | ✅ | 🔧 | |
| `matrix_norm` | ✅ | ⬜ | NumPy 2.0 |
| `vector_norm` | ✅ | ⬜ | NumPy 2.0 |
| `cond` | ✅ | ⬜ | |
| `det` | ✅ | 🔧 | |
| `matrix_rank` | ✅ | ⬜ | |
| `slogdet` | ✅ | ⬜ | |
| `trace` | ✅ | ⬜ | |
| **Solving** | | | |
| `solve` | ✅ | 🔧 | |
| `tensorsolve` | ✅ | ⬜ | |
| `lstsq` | ✅ | ⬜ | |
| `inv` | ✅ | 🔧 | |
| `pinv` | ✅ | ⬜ | |
| `tensorinv` | ✅ | ⬜ | |
| **Other** | | | |
| `diagonal` | ✅ | 🔧 | |
| `matrix_transpose` | ⬜ | ⬜ | NumPy 2.0 |

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
| `PCG64DXSM` | ⬜ | TODO |
| `Philox` | ✅ | |
| `SFC64` | ✅ | |

### Utility Functions

| Function | TS | WASM | Notes |
|----------|:--:|:----:|-------|
| `seed` | ✅ | ⬜ | |
| `random` | ✅ | 🔧 | |
| `randn` | ✅ | 🔧 | |
| `randint` | ✅ | 🔧 | |
| `rand` | ⬜ | ⬜ | Legacy |
| `ranf` | ⬜ | ⬜ | Legacy |
| `random_sample` | ⬜ | ⬜ | Legacy |
| `random_integers` | ⬜ | ⬜ | Deprecated |
| `sample` | ⬜ | ⬜ | Legacy |
| `bytes` | ⬜ | ⬜ | TODO |
| `choice` | ⬜ | ⬜ | **HIGH PRIORITY** |
| `shuffle` | ⬜ | ⬜ | **HIGH PRIORITY** |
| `permutation` | ⬜ | ⬜ | **HIGH PRIORITY** |
| `get_state` | ⬜ | ⬜ | Legacy |
| `set_state` | ⬜ | ⬜ | Legacy |

### Distributions (Generator methods)

| Distribution | TS | WASM | Notes |
|--------------|:--:|:----:|-------|
| `beta` | ⬜ | 🔧 | Has WASM |
| `binomial` | ⬜ | 🔧 | Has WASM |
| `chisquare` | ⬜ | 🔧 | Has WASM |
| `dirichlet` | ⬜ | ⬜ | TODO |
| `exponential` | ⬜ | 🔧 | Has WASM |
| `f` | ⬜ | 🔧 | Has WASM |
| `gamma` | ⬜ | 🔧 | Has WASM |
| `geometric` | ⬜ | 🔧 | Has WASM |
| `gumbel` | ⬜ | 🔧 | Has WASM |
| `hypergeometric` | ⬜ | 🔧 | Has WASM |
| `laplace` | ⬜ | 🔧 | Has WASM |
| `logistic` | ⬜ | 🔧 | Has WASM |
| `lognormal` | ⬜ | 🔧 | Has WASM |
| `logseries` | ⬜ | 🔧 | Has WASM |
| `multinomial` | ⬜ | ⬜ | TODO |
| `multivariate_normal` | ⬜ | ⬜ | TODO |
| `negative_binomial` | ⬜ | 🔧 | Has WASM |
| `noncentral_chisquare` | ⬜ | 🔧 | Has WASM |
| `noncentral_f` | ⬜ | 🔧 | Has WASM |
| `normal` | ⬜ | 🔧 | Has WASM |
| `pareto` | ⬜ | 🔧 | Has WASM |
| `poisson` | ⬜ | 🔧 | Has WASM |
| `power` | ✅ | 🔧 | |
| `rayleigh` | ⬜ | 🔧 | Has WASM |
| `standard_cauchy` | ⬜ | 🔧 | Has WASM |
| `standard_exponential` | ⬜ | 🔧 | Has WASM |
| `standard_gamma` | ⬜ | 🔧 | Has WASM |
| `standard_normal` | ⬜ | 🔧 | Has WASM |
| `standard_t` | ⬜ | 🔧 | Has WASM |
| `triangular` | ⬜ | 🔧 | Has WASM |
| `uniform` | ⬜ | 🔧 | Has WASM |
| `vonmises` | ⬜ | 🔧 | Has WASM |
| `wald` | ⬜ | 🔧 | Has WASM |
| `weibull` | ⬜ | 🔧 | Has WASM |
| `zipf` | ⬜ | 🔧 | Has WASM |

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
| `mod` | ⬜ | TODO |
| `translate` | ⬜ | TODO |
| `slice` | ⬜ | NumPy 2.0 |

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
| `assert_equal` | ⬜ | TODO |
| `assert_almost_equal` | ⬜ | TODO |
| `assert_array_equal` | ⬜ | TODO |
| `assert_array_almost_equal` | ⬜ | TODO |
| `assert_allclose` | ⬜ | TODO |
| `assert_raises` | ⬜ | TODO |
| `assert_warns` | ⬜ | TODO |
| Other assertions | ⬜ | Low priority |

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
| `prod` | Statistics | Product reduction |
| `clip` | Math | Value clipping - very common |
| `diff` | Math | N-th differences |
| `gradient` | Math | Numerical gradient |
| `convolve` | Math | 1D convolution |
| `correlate` | Math | 1D correlation |
| `cov` | Statistics | Covariance matrix |
| `corrcoef` | Statistics | Correlation coefficients |
| `average` | Statistics | Weighted average |
| `ptp` | Statistics | Peak-to-peak (max-min) |

### Priority 2: Export Functions with WASM Backing

These have C implementations but need TypeScript exports:
- `reshape`, `ravel`, `squeeze`, `expand_dims`, `transpose`, `swapaxes`
- `copy`, `ascontiguousarray`, `asfortranarray`
- `percentile`, `quantile`
- `tri`, `tril`, `triu`

### Priority 3: Complete numpy.random

Most distributions have WASM backing but need TypeScript API:
- All 30+ distributions listed in section 4
- `choice`, `shuffle`, `permutation`
- `dirichlet`, `multinomial`, `multivariate_normal`

### Priority 4: NumPy 2.0 Compatibility

Add aliases and new functions from NumPy 2.0:
- `cumulative_sum`, `cumulative_prod`
- `matrix_transpose`, `vecdot`, `matvec`, `vecmat`
- `concat` (alias for `concatenate`)
- Various trigonometric aliases (`asin`, `acos`, etc.)

---

## Action Items Summary

### Immediate (High Impact)
1. Implement `prod`, `clip`, `average`, `ptp`
2. Implement `diff`, `gradient`, `convolve`, `correlate`
3. Implement `cov`, `corrcoef`
4. Export `reshape`, `transpose`, etc. as top-level functions

### Short-term
5. Complete numpy.random TypeScript API for all distributions
6. Add `choice`, `shuffle`, `permutation`
7. Add `diagflat`, `tril`, `triu`, `vander`

### Medium-term
8. Add NumPy 2.0 aliases and new functions
9. Complete numpy.ma module
10. Add missing numpy.testing assertions
