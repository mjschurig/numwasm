# NumWasm NumPy Implementation TODO

This document tracks the implementation status of NumPy functions in NumWasm.

**Goal:** Export all NumPy package/function names from this package.

## Legend

| Symbol | Meaning |
|--------|---------|
| ✅ | Implemented (TypeScript API exists) |
| 🔧 | Has WASM backing (C implementation) |
| ⬜ | Not implemented |

## Summary

| Module | Implemented | Total | Coverage |
|--------|-------------|-------|----------|
| **numpy (main)** | 295 | 484 | 61% |
| **numpy.linalg** | 30 | 32 | 94% |
| **numpy.fft** | 18 | 18 | 100% |
| **numpy.random** | 5 | 50 | 10% |
| **numpy.ma** | 124 | 224 | 55% |
| **numpy.polynomial** | 111 | 200 | 56% |
| **numpy.strings** | 38 | 46 | 83% |
| **numpy.testing** | 2 | 48 | 4% |
| **Total** | ~623 | ~1102 | ~57% |

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
└─────────────────────────────────────────────────────────────────┘
```

### WASM C Implementation Status

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

## Priority TODO Items

### Priority 1: High-Impact Missing Functions

These are commonly-used NumPy functions that should be added:

| Function | Category | Notes |
|----------|----------|-------|
| `prod` | Statistics | Product reduction - common operation |
| `clip` | Math | Value clipping - very common |
| `diff` | Math | N-th differences |
| `gradient` | Math | Numerical gradient |
| `convolve` | Math | 1D convolution |
| `correlate` | Math | 1D correlation |
| `cov` | Statistics | Covariance matrix |
| `corrcoef` | Statistics | Correlation coefficients |
| `average` | Statistics | Weighted average |
| `ptp` | Statistics | Peak-to-peak (max-min) |
| `percentile` | Statistics | Percentile calculation |
| `quantile` | Statistics | Quantile calculation |

### Priority 2: numpy.random Module (5/50 implemented)

The random module has WASM backing for most distributions but lacks TypeScript API:

**Already have WASM, need TypeScript API:**
- `beta`, `binomial`, `chisquare`, `exponential`, `f`
- `gamma`, `geometric`, `gumbel`, `hypergeometric`, `laplace`
- `logistic`, `lognormal`, `logseries`, `negative_binomial`
- `noncentral_chisquare`, `noncentral_f`, `normal`, `pareto`
- `poisson`, `rayleigh`, `standard_cauchy`, `standard_exponential`
- `standard_gamma`, `standard_normal`, `standard_t`, `triangular`
- `uniform`, `vonmises`, `wald`, `weibull`, `zipf`

**Need both WASM and TypeScript:**
- `choice`, `shuffle`, `permutation`
- `dirichlet`, `multinomial`, `multivariate_normal`

### Priority 3: Array Manipulation Functions

These functions are marked as not implemented but are commonly used:

| Function | Notes |
|----------|-------|
| `reshape` | Top-level function (method exists) |
| `ravel` | Top-level function (method exists) |
| `squeeze` | Top-level function (method exists) |
| `expand_dims` | Top-level function (method exists) |
| `transpose` | Top-level function (method exists) |
| `swapaxes` | Top-level function (method exists) |
| `moveaxis` | Move array axis |
| `rollaxis` | Roll axis |
| `copy` | Deep copy |

### Priority 4: Type System Functions

| Function | Notes |
|----------|-------|
| `can_cast` | Check type casting |
| `result_type` | Determine result type |
| `astype` | Top-level type conversion |
| `issubdtype` | Check dtype hierarchy |

---

## numpy (Main Namespace) - Detailed Status

### Array Creation (20/35)

| Function | Status | WASM | Notes |
|----------|--------|------|-------|
| `array` | ✅ | 🔧 | |
| `zeros` | ✅ | ⬜ | |
| `ones` | ✅ | ⬜ | |
| `empty` | ✅ | 🔧 | |
| `full` | ✅ | 🔧 | |
| `arange` | ✅ | ⬜ | |
| `linspace` | ✅ | ⬜ | |
| `logspace` | ✅ | ⬜ | |
| `geomspace` | ✅ | ⬜ | |
| `eye` | ✅ | ⬜ | |
| `identity` | ✅ | ⬜ | |
| `diag` | ✅ | 🔧 | |
| `diagflat` | ⬜ | ⬜ | TODO |
| `tri` | ⬜ | 🔧 | Has WASM |
| `tril` | ⬜ | ⬜ | TODO |
| `triu` | ⬜ | ⬜ | TODO |
| `vander` | ⬜ | ⬜ | TODO |
| `zeros_like` | ✅ | ⬜ | |
| `ones_like` | ✅ | ⬜ | |
| `empty_like` | ✅ | ⬜ | |
| `full_like` | ✅ | ⬜ | |
| `fromfunction` | ⬜ | ⬜ | TODO |
| `fromiter` | ⬜ | ⬜ | TODO |
| `fromstring` | ⬜ | ⬜ | TODO |
| `frombuffer` | ✅ | ⬜ | |
| `fromfile` | ✅ | ⬜ | |
| `fromregex` | ✅ | ⬜ | |
| `copy` | ⬜ | 🔧 | Has WASM, needs TS |
| `asarray` | ✅ | ⬜ | |
| `asanyarray` | ⬜ | ⬜ | TODO |
| `ascontiguousarray` | ⬜ | 🔧 | Has WASM |
| `asfortranarray` | ⬜ | 🔧 | Has WASM |
| `asarray_chkfinite` | ⬜ | ⬜ | TODO |
| `require` | ⬜ | ⬜ | TODO |
| `from_dlpack` | ⬜ | ⬜ | Low priority |

### Array Manipulation (31/39)

| Function | Status | WASM | Notes |
|----------|--------|------|-------|
| `reshape` | ⬜ | 🔧 | Has WASM, needs TS export |
| `ravel` | ⬜ | 🔧 | Has WASM, needs TS export |
| `squeeze` | ⬜ | 🔧 | Has WASM, needs TS export |
| `expand_dims` | ⬜ | 🔧 | Has WASM, needs TS export |
| `transpose` | ⬜ | 🔧 | Has WASM, needs TS export |
| `swapaxes` | ⬜ | 🔧 | Has WASM, needs TS export |
| `moveaxis` | ⬜ | ⬜ | TODO |
| `rollaxis` | ⬜ | ⬜ | TODO |
| `atleast_1d` | ✅ | ⬜ | |
| `atleast_2d` | ✅ | ⬜ | |
| `atleast_3d` | ✅ | ⬜ | |
| `concatenate` | ✅ | 🔧 | |
| `stack` | ✅ | ⬜ | |
| `vstack` | ✅ | ⬜ | |
| `hstack` | ✅ | ⬜ | |
| `dstack` | ✅ | ⬜ | |
| `column_stack` | ✅ | ⬜ | |
| `split` | ✅ | ⬜ | |
| `array_split` | ✅ | ⬜ | |
| `vsplit` | ✅ | ⬜ | |
| `hsplit` | ✅ | ⬜ | |
| `dsplit` | ✅ | ⬜ | |
| `unstack` | ✅ | ⬜ | |
| `tile` | ✅ | ⬜ | |
| `repeat` | ✅ | ⬜ | |
| `pad` | ✅ | ⬜ | |
| `flip` | ✅ | ⬜ | |
| `fliplr` | ✅ | ⬜ | |
| `flipud` | ✅ | ⬜ | |
| `roll` | ✅ | ⬜ | |
| `rot90` | ✅ | ⬜ | |
| `resize` | ✅ | ⬜ | |
| `trim_zeros` | ✅ | ⬜ | |
| `insert` | ✅ | ⬜ | |
| `delete` | ⬜ | ⬜ | Exported as `deleteArr` |
| `append` | ✅ | ⬜ | |
| `copyto` | ✅ | ⬜ | |
| `block` | ✅ | ⬜ | |

### Math - Basic (20/22)

| Function | Status | WASM | Notes |
|----------|--------|------|-------|
| `add` | ✅ | 🔧 | |
| `subtract` | ✅ | 🔧 | |
| `multiply` | ✅ | 🔧 | |
| `divide` | ✅ | 🔧 | |
| `true_divide` | ✅ | 🔧 | |
| `floor_divide` | ✅ | 🔧 | |
| `negative` | ✅ | 🔧 | |
| `positive` | ✅ | 🔧 | |
| `power` | ✅ | 🔧 | |
| `pow` | ⬜ | 🔧 | Has WASM, needs TS alias |
| `remainder` | ✅ | 🔧 | |
| `mod` | ✅ | 🔧 | |
| `fmod` | ✅ | 🔧 | |
| `divmod` | ✅ | 🔧 | |
| `absolute` | ✅ | 🔧 | |
| `abs` | ✅ | 🔧 | |
| `fabs` | ⬜ | ⬜ | TODO |
| `sign` | ✅ | 🔧 | |
| `reciprocal` | ✅ | 🔧 | |
| `sqrt` | ✅ | 🔧 | |
| `square` | ✅ | 🔧 | |
| `cbrt` | ✅ | 🔧 | |

### Math - Special (3/10) - HIGH PRIORITY

| Function | Status | WASM | Notes |
|----------|--------|------|-------|
| `clip` | ⬜ | ⬜ | **HIGH PRIORITY** |
| `convolve` | ⬜ | ⬜ | **HIGH PRIORITY** |
| `correlate` | ⬜ | ⬜ | **HIGH PRIORITY** |
| `diff` | ⬜ | ⬜ | **HIGH PRIORITY** |
| `gradient` | ⬜ | ⬜ | **HIGH PRIORITY** |
| `interp` | ⬜ | ⬜ | TODO |
| `trapezoid` | ⬜ | ⬜ | TODO |
| `sinc` | ✅ | 🔧 | |
| `i0` | ✅ | ⬜ | |
| `heaviside` | ✅ | 🔧 | |

### Statistics (34/43) - NEEDS WORK

| Function | Status | WASM | Notes |
|----------|--------|------|-------|
| `sum` | ✅ | 🔧 | |
| `prod` | ⬜ | ⬜ | **HIGH PRIORITY** |
| `mean` | ✅ | 🔧 | |
| `std` | ✅ | 🔧 | |
| `var` | ✅ | 🔧 | Exported as `var_` |
| `min` | ✅ | 🔧 | |
| `max` | ✅ | 🔧 | |
| `argmin` | ✅ | 🔧 | |
| `argmax` | ✅ | 🔧 | |
| `nanmin` | ✅ | ⬜ | |
| `nanmax` | ✅ | ⬜ | |
| `nansum` | ✅ | 🔧 | |
| `nanprod` | ✅ | ⬜ | |
| `nanmean` | ✅ | 🔧 | |
| `nanstd` | ✅ | 🔧 | |
| `nanvar` | ✅ | 🔧 | |
| `nanargmin` | ✅ | ⬜ | |
| `nanargmax` | ✅ | ⬜ | |
| `nanmedian` | ✅ | ⬜ | |
| `nanpercentile` | ✅ | ⬜ | |
| `nanquantile` | ✅ | ⬜ | |
| `median` | ✅ | 🔧 | |
| `average` | ⬜ | ⬜ | **HIGH PRIORITY** |
| `percentile` | ⬜ | 🔧 | Has WASM |
| `quantile` | ⬜ | 🔧 | Has WASM |
| `ptp` | ⬜ | ⬜ | **HIGH PRIORITY** |
| `corrcoef` | ⬜ | ⬜ | **HIGH PRIORITY** |
| `cov` | ⬜ | ⬜ | **HIGH PRIORITY** |
| `histogram` | ✅ | ⬜ | |
| `histogram2d` | ✅ | ⬜ | |
| `histogramdd` | ✅ | ⬜ | |
| `histogram_bin_edges` | ✅ | ⬜ | |
| `bincount` | ✅ | ⬜ | |
| `digitize` | ✅ | ⬜ | |
| `count_nonzero` | ✅ | 🔧 | |
| `cumsum` | ✅ | 🔧 | |
| `cumprod` | ✅ | 🔧 | |
| `nancumsum` | ✅ | 🔧 | |
| `nancumprod` | ✅ | 🔧 | |
| `cumulative_sum` | ⬜ | ⬜ | NumPy 2.0 |
| `cumulative_prod` | ⬜ | ⬜ | NumPy 2.0 |

---

## numpy.linalg (30/32) - Nearly Complete

| Function | Status | WASM | Notes |
|----------|--------|------|-------|
| `LinAlgError` | ✅ | ⬜ | |
| `cholesky` | ✅ | 🔧 | |
| `cond` | ✅ | ⬜ | |
| `cross` | ✅ | ⬜ | |
| `det` | ✅ | 🔧 | |
| `diagonal` | ✅ | 🔧 | |
| `eig` | ✅ | 🔧 | |
| `eigh` | ✅ | ⬜ | |
| `eigvals` | ✅ | ⬜ | |
| `eigvalsh` | ✅ | ⬜ | |
| `inv` | ✅ | 🔧 | |
| `lstsq` | ✅ | ⬜ | |
| `matmul` | ✅ | 🔧 | |
| `matrix_norm` | ✅ | ⬜ | |
| `matrix_power` | ✅ | ⬜ | |
| `matrix_rank` | ✅ | ⬜ | |
| `matrix_transpose` | ⬜ | ⬜ | TODO |
| `multi_dot` | ✅ | ⬜ | |
| `norm` | ✅ | 🔧 | |
| `outer` | ✅ | ⬜ | |
| `pinv` | ✅ | ⬜ | |
| `qr` | ✅ | 🔧 | |
| `slogdet` | ✅ | ⬜ | |
| `solve` | ✅ | 🔧 | |
| `svd` | ✅ | 🔧 | |
| `svdvals` | ✅ | ⬜ | |
| `tensordot` | ✅ | ⬜ | |
| `tensorinv` | ✅ | ⬜ | |
| `tensorsolve` | ✅ | ⬜ | |
| `trace` | ✅ | ⬜ | |
| `vecdot` | ⬜ | ⬜ | TODO |
| `vector_norm` | ✅ | ⬜ | |

---

## numpy.fft (18/18) - Complete ✅

All FFT functions are implemented:
- `fft`, `ifft`, `fft2`, `ifft2`, `fftn`, `ifftn`
- `rfft`, `irfft`, `rfft2`, `irfft2`, `rfftn`, `irfftn`
- `hfft`, `ihfft`
- `fftfreq`, `rfftfreq`, `fftshift`, `ifftshift`

---

## numpy.random (5/50) - Needs Major Work

### Currently Implemented
- `power`, `randint`, `randn`, `random`, `seed`

### Has WASM, Needs TypeScript API
| Distribution | WASM Status |
|--------------|-------------|
| `beta` | 🔧 |
| `binomial` | 🔧 |
| `chisquare` | 🔧 |
| `exponential` | 🔧 |
| `f` | 🔧 |
| `gamma` | 🔧 |
| `geometric` | 🔧 |
| `gumbel` | 🔧 |
| `hypergeometric` | 🔧 |
| `laplace` | 🔧 |
| `logistic` | 🔧 |
| `lognormal` | 🔧 |
| `logseries` | 🔧 |
| `negative_binomial` | 🔧 |
| `noncentral_chisquare` | 🔧 |
| `noncentral_f` | 🔧 |
| `normal` | 🔧 |
| `pareto` | 🔧 |
| `poisson` | 🔧 |
| `rayleigh` | 🔧 |
| `standard_cauchy` | 🔧 |
| `standard_exponential` | 🔧 |
| `standard_gamma` | 🔧 |
| `standard_normal` | 🔧 |
| `standard_t` | 🔧 |
| `triangular` | 🔧 |
| `uniform` | 🔧 |
| `vonmises` | 🔧 |
| `wald` | 🔧 |
| `weibull` | 🔧 |
| `zipf` | 🔧 |

### Needs Implementation
- `choice`, `shuffle`, `permutation`
- `dirichlet`, `multinomial`, `multivariate_normal`
- `bytes`, `rand`, `random_integers`, `random_sample`, `ranf`, `sample`
- `get_state`, `set_state`

---

## numpy.ma (Masked Arrays) (124/224)

Masked array module is ~55% complete. Most math operations are implemented through delegation to the main numpy functions.

### Key Missing Functions
- `anom`, `anomalies`, `average`
- `clip`, `compressed`, `convolve`, `correlate`, `cov`, `corrcoef`
- `diff`, `filled`, `fix_invalid`
- `flatten_mask`, `flatten_structured_array`
- `getdata`, `getmask`, `getmaskarray`
- `make_mask`, `make_mask_descr`, `make_mask_none`
- `masked_*` family (equal, greater, inside, invalid, less, outside, values, where)
- `notmasked_*` family
- `prod`, `ptp`
- `soften_mask`, `harden_mask`

---

## numpy.polynomial (111/200)

Polynomial module is ~56% complete. Core polynomial classes and operations are implemented.

### Implemented Classes
- `Polynomial`, `Chebyshev`, `Hermite`, `HermiteE`, `Laguerre`, `Legendre`

### Key Missing Functions
- `*domain`, `*line`, `*one`, `*zero`, `*x` constants for each family
- `*gauss` (Gaussian quadrature)
- `*grid2d`, `*grid3d` (grid evaluation)
- `*vander2d`, `*vander3d` (2D/3D Vandermonde)
- `*weight` (weight functions)
- `set_default_printstyle`

---

## numpy.strings (38/46)

String operations are ~83% complete.

### Missing
- `capitalize`, `center`, `count`, `find`, `index`
- `replace`, `title`, `translate`

---

## numpy.testing (2/48)

Testing module is minimally implemented. Most functions are low priority for a runtime library.

### Implemented
- `KnownFailureException`, `SkipTest`

---

## Not Planned for NumWasm

These are Python-specific or low-priority items:

### Python-Specific
- `datetime64`, `timedelta64` - Use JavaScript Date
- `busday_*`, `busdaycalendar` - Business day functions
- All `*_` scalar types (`int_`, `float_`, `bool_`, etc.)
- `matrix` class (deprecated in NumPy)
- `ufunc` class registration

### Environment-Specific
- `errstate`, `geterr`, `seterr` - Error state management
- `getbufsize`, `setbufsize` - Buffer management
- `show_config`, `show_runtime`, `get_include`
- `may_share_memory`, `shares_memory` - Memory introspection

### Low Priority
- `from_dlpack` - DLPack interop
- `nested_iters` - Advanced iteration
- `packbits`, `unpackbits` - Bit packing
- Most `numpy.testing` functions

---

## Action Items Summary

### Immediate (High Impact)
1. Add TypeScript exports for functions with WASM backing:
   - `reshape`, `ravel`, `squeeze`, `expand_dims`, `transpose`, `swapaxes`
   - `copy`, `ascontiguousarray`, `asfortranarray`
   - `percentile`, `quantile`
   - `pow` (alias for `power`)

2. Implement high-priority missing functions:
   - `prod`, `clip`, `average`, `ptp`
   - `diff`, `gradient`
   - `cov`, `corrcoef`
   - `convolve`, `correlate`

### Short-term
3. Complete numpy.random TypeScript API for all distributions with WASM backing

4. Add remaining array manipulation functions:
   - `moveaxis`, `rollaxis`
   - `diagflat`, `tril`, `triu`, `vander`

### Medium-term
5. Complete numpy.ma module
6. Complete numpy.polynomial module
7. Add remaining string functions
