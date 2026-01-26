# NumJS-WASM Implementation Plan

Complete implementation roadmap for a TypeScript/WebAssembly NumPy clone with proper dependency ordering.

---

## ⚠️ Implementation Guideline

**The original NumPy source code is available in `/numpy`.**

Implementation should follow the original NumPy code as closely as possible to ensure:
- API compatibility with NumPy
- Consistent behavior and edge case handling
- Numerical accuracy matching NumPy's algorithms
- Familiar patterns for users coming from Python/NumPy

---

## NumPy Reference Source Tree

```
/numpy/
├── numpy/                              # Main package
│   ├── __init__.py                     # Entry point (~400+ exports)
│   │
│   ├── _core/                          # ⭐ CORE IMPLEMENTATION (most important)
│   │   ├── multiarray.py               # Array construction/manipulation (~1,740 lines)
│   │   ├── numeric.py                  # Numerical operations (~2,711 lines)
│   │   ├── fromnumeric.py              # Array creation functions (~4,233 lines)
│   │   ├── umath.py                    # Universal math functions
│   │   ├── numerictypes.py             # Data type handling
│   │   ├── _dtype.py                   # dtype representation
│   │   ├── einsumfunc.py               # Einstein summation (~1,650 lines)
│   │   ├── arrayprint.py               # Array printing (~1,779 lines)
│   │   ├── shape_base.py               # Shape manipulation
│   │   ├── function_base.py            # Utility functions
│   │   │
│   │   └── src/                        # ⭐ C/Cython implementations
│   │       ├── multiarray/             # Core array (150+ .c files)
│   │       │   ├── arrayobject.c       # ndarray object implementation
│   │       │   ├── methods.c           # Array methods
│   │       │   ├── mapping.c           # Indexing/slicing
│   │       │   ├── item_selection.c    # Item access
│   │       │   ├── shape.c             # Shape operations
│   │       │   ├── conversion_utils.c  # Type conversion
│   │       │   ├── dtype_transfer.c    # dtype handling
│   │       │   ├── descriptor.c        # dtype descriptors
│   │       │   ├── iterator.c          # Iterator implementation
│   │       │   ├── nditer_*.c          # N-dimensional iteration
│   │       │   └── sorting/            # Sort algorithms
│   │       │
│   │       ├── umath/                  # Universal math functions
│   │       │   ├── umathmodule.c       # Module init
│   │       │   ├── ufunc_object.c      # Ufunc object
│   │       │   └── reduction.c         # Reduction operations
│   │       │
│   │       ├── common/                 # Common utilities
│   │       ├── npymath/                # Math library
│   │       ├── npysort/                # Sorting algorithms
│   │       └── _simd/                  # SIMD optimizations
│   │
│   ├── lib/                            # Library functions (~50+ files)
│   │   ├── _arraypad_impl.py           # Array padding
│   │   ├── _arraysetops_impl.py        # Set operations
│   │   ├── _nanfunctions_impl.py       # NaN handling
│   │   ├── _histograms_impl.py         # Histograms
│   │   └── _index_tricks_impl.py       # Indexing tricks
│   │
│   ├── linalg/                         # Linear algebra
│   │   ├── _linalg.py                  # Main interface
│   │   └── lapack_lite/                # LAPACK routines
│   │
│   ├── fft/                            # Fast Fourier Transforms
│   │   ├── _pocketfft.py               # FFT interface
│   │   └── pocketfft/                  # FFT implementation
│   │
│   ├── random/                         # Random number generation
│   │   ├── _generator.pyx              # Generator class (~200KB)
│   │   ├── _pcg64.pyx                  # PCG64 generator
│   │   └── src/                        # C implementations
│   │
│   ├── polynomial/                     # Polynomial operations
│   ├── ma/                             # Masked arrays
│   └── testing/                        # Testing utilities
```

### Key Files to Reference

| Component | NumPy Reference File |
|-----------|---------------------|
| ndarray core | `numpy/_core/src/multiarray/arrayobject.c` |
| dtype system | `numpy/_core/src/multiarray/descriptor.c` |
| Indexing/slicing | `numpy/_core/src/multiarray/mapping.c` |
| Shape operations | `numpy/_core/src/multiarray/shape.c` |
| Ufunc infrastructure | `numpy/_core/src/umath/ufunc_object.c` |
| Reductions | `numpy/_core/src/umath/reduction.c` |
| Type promotion | `numpy/_core/numerictypes.py` |
| Array creation | `numpy/_core/fromnumeric.py` |

---

## Implementation Status

Legend: ✅ Complete | 🔶 Partial | ❌ Not Started

---

## Dependency Tree

```
LEVEL 0: FOUNDATION
├── ✅ Memory Management (existing)
├── ✅ DType System
│   ├── ✅ Core Types (Float32, Float64, Int32, Int64)
│   ├── ✅ Extended Types (Bool, UInt8-64, Int8, Int16, Float16, Complex64, Complex128)
│   └── ✅ Type Promotion Rules
└── ✅ NDArray Core Structure

LEVEL 1: CORE OPERATIONS
├── ✅ Element Access (get/set by index)
├── ✅ Iteration Infrastructure
└── ✅ Shape Manipulation (reshape, transpose, ravel, squeeze, expand_dims, moveaxis, atleast_Nd)

LEVEL 2: VIEWS & BROADCASTING
├── ✅ View System (shared data, different shape/strides)
├── ✅ Slicing (basic, integer array, boolean array)
├── ✅ Broadcasting (shape compatibility, stride adjustment)
└── ✅ Index Generation (indices, ix_, diag_indices, tril/triu_indices, argwhere)

LEVEL 3: UNIVERSAL FUNCTIONS
├── ❌ Ufunc Infrastructure
├── ❌ Unary Ufuncs (neg, abs, sqrt, exp, log, trig, rounding, predicates)
├── ❌ Binary Ufuncs (arithmetic, comparison, logical, bitwise)
└── 🔶 Reductions with Axis (sum only, without axis support)

LEVEL 4: MANIPULATION & STATISTICS
├── ❌ Array Manipulation (concat, stack, split, flip, roll, tile)
├── 🔶 Sorting & Searching (argsort 1D only, argmax/argmin partial)
├── ❌ Statistics (median, percentile, histogram, corrcoef)
└── ❌ Set Operations (unique, intersect1d, union1d)

LEVEL 5: HIGHER-LEVEL MODULES
├── ❌ numpy.linalg (dot, matmul, solve, inv, eig, svd)
├── ❌ numpy.fft (fft, ifft, fft2, fftn, fftfreq)
└── ❌ numpy.random (Generator, PCG64, distributions)
```

---

## Phase 1: Foundation Enhancement

### 1.1 Extended DType System ✅ COMPLETE
```
numpy.dtypes
├── ✅ BoolDType
├── ✅ Int8DType, Int16DType, Int32DType, Int64DType
├── ✅ UInt8DType, UInt16DType, UInt32DType, UInt64DType
├── ✅ Float16DType, Float32DType, Float64DType
├── ✅ Complex64DType, Complex128DType
└── ✅ Type Promotion
    ├── ✅ can_cast(from, to) → canCast()
    ├── ✅ result_type(*arrays_or_dtypes) → promoteTypes()
    └── ✅ common_type(*arrays) → commonType()
```

**Implementation:** `src/ts/types.ts`, `src/ts/dtype.ts`, `src/wasm/dtype.c`

### 1.2 Element Access ✅ COMPLETE
```
NDArray
├── ✅ get(...indices) → scalar
├── ✅ set(...indices, value)
├── ✅ item() → scalar (for 0-d or single element)
├── ✅ itemset(value)
├── ✅ getFlat(index), setFlat(index, value)
├── ✅ getComplex(), setComplex() (for complex types)
├── ✅ toArray() → number[]
├── ✅ toTypedArray()
└── ✅ flat → FlatIterator
```

**Implementation:** `src/ts/NDArray.ts`, `src/wasm/ndarray.c`

### 1.3 Array Creation ✅ COMPLETE
```
Creation Functions
├── ✅ zeros(shape, dtype)
├── ✅ ones(shape, dtype)
├── ✅ empty(shape, dtype)
├── ✅ full(shape, fill_value, dtype)
├── ✅ zerosLike(arr), onesLike(arr), emptyLike(arr), fullLike(arr)
├── ✅ fromArray(data, shape, dtype) - supports nested arrays, shape inference
├── ✅ arange(start, stop, step, dtype)
├── ✅ linspace(start, stop, num, endpoint, dtype)
├── ✅ logspace(start, stop, num, endpoint, base, dtype)
├── ✅ geomspace(start, stop, num, endpoint, dtype)
├── ✅ eye(N, M, k, dtype)
├── ✅ identity(n, dtype)
├── ✅ diag(v, k)
├── ✅ tri(N, M, k, dtype)
├── ✅ tril(arr, k), triu(arr, k)
└── ✅ meshgrid(*xi, indexing)
```

**Implementation:** `src/ts/NDArray.ts`

---

## Phase 2: Iteration & Shape Manipulation

### 2.1 Iterator Infrastructure ✅ COMPLETE
```
Iteration
├── ✅ nditer(arr) → Iterator
├── ✅ ndenumerate(arr) → Iterator<[index, value]>
├── ✅ ndindex(*shape) → Iterator<index>
└── ✅ FlatIterator → Iterator<scalar>
```

**Implementation:** `src/ts/iterators.ts`

### 2.2 Shape Operations ✅ COMPLETE
```
Shape Manipulation
├── ✅ reshape(arr, newshape) - supports -1 dimension
├── ✅ ravel(arr) - view when contiguous
├── ✅ flatten(arr) → always copy
├── ✅ transpose(arr, axes) - custom axis permutation
├── ✅ moveaxis(arr, source, destination)
├── ✅ swapaxes(arr, axis1, axis2)
├── ✅ atleast_1d(*arrs), atleast_2d(*arrs), atleast_3d(*arrs)
├── ✅ expand_dims(arr, axis)
├── ✅ squeeze(arr, axis)
└── ✅ Properties
    ├── ✅ ndim
    ├── ✅ shape
    ├── ✅ size
    ├── ✅ T (transpose)
    └── ✅ flat
```

**Implementation:** `src/ts/NDArray.ts`, `src/ts/indexing.ts`, `src/wasm/ndarray.c`

---

## Phase 3: Views, Slicing & Broadcasting

### 3.1 View System ✅ COMPLETE
```
Views
├── ✅ view(shape, strides) → view with custom layout
├── ✅ viewDtype(dtype) → view with different dtype
├── ✅ copy() → deep copy
├── ✅ ascontiguousarray(arr)
├── ✅ asfortranarray(arr)
└── ✅ flags: c_contiguous, f_contiguous, owndata, writeable, aligned
```

**Implementation:** `src/ts/NDArray.ts`, `src/wasm/ndarray.c`

### 3.2 Indexing & Slicing ✅ COMPLETE
```
Indexing
├── ✅ Basic Slicing
│   ├── ✅ arr.at(i) → element or subarray
│   ├── ✅ arr.slice([start:stop:step])
│   ├── ✅ arr.slice([..., i]) (ellipsis)
│   └── ✅ newaxis support (expand dims)
├── ❌ Advanced Indexing (not yet WASM-accelerated)
│   ├── ❌ arr[int_array] → fancy indexing
│   └── ❌ arr[bool_array] → boolean masking
├── ✅ Index Functions
│   ├── ✅ take(arr, indices, axis)
│   ├── ✅ take_along_axis(arr, indices, axis)
│   ├── ✅ put(arr, indices, values)
│   ├── ✅ put_along_axis(arr, indices, values, axis)
│   ├── ✅ putmask(arr, mask, values)
│   ├── ✅ place(arr, mask, vals)
│   ├── ✅ compress(condition, arr, axis)
│   ├── ✅ extract(condition, arr)
│   ├── ✅ choose(a, choices)
│   └── ✅ select(condlist, choicelist, default)
└── ✅ Index Generation
    ├── ✅ indices(dimensions)
    ├── ✅ ix_(*args)
    ├── ❌ ogrid, mgrid (lower priority, complex slice-based API)
    ├── ✅ diag_indices(n, ndim)
    ├── ✅ tril_indices(n, k, m)
    ├── ✅ triu_indices(n, k, m)
    ├── ❌ mask_indices(n, mask_func, k)
    ├── ✅ nonzero(arr)
    ├── ✅ flatnonzero(arr)
    ├── ✅ argwhere(arr)
    ├── ✅ where(condition, x, y)
    ├── ✅ ravelMultiIndex(multi_index, dims)
    └── ✅ unravelIndex(indices, shape)
```

**Implementation:** `src/ts/slice.ts`, `src/ts/indexing.ts`, `src/wasm/indexing.c`

### 3.3 Broadcasting ✅ COMPLETE
```
Broadcasting
├── ✅ broadcastTo(arr, shape)
├── ✅ broadcastArrays(*args)
├── ✅ broadcastShapes(*shapes)
├── ✅ broadcastShapesMulti(shapes[])
├── ✅ shapesAreBroadcastable()
└── ✅ computeBroadcastStrides() (internal)
```

**Implementation:** `src/ts/broadcast.ts`, `src/wasm/broadcast.c`

---

## Phase 4: Universal Functions (Ufuncs) ❌ NOT STARTED

### 4.1 Ufunc Infrastructure ❌
```
Ufunc System
├── ❌ Ufunc class
│   ├── ❌ __call__(inputs, out, where, casting, dtype)
│   ├── ❌ reduce(arr, axis, dtype, out, keepdims, initial)
│   ├── ❌ accumulate(arr, axis, dtype, out)
│   ├── ❌ reduceat(arr, indices, axis, dtype, out)
│   ├── ❌ outer(a, b, out)
│   └── ❌ at(arr, indices, b)
├── ❌ Inner loop dispatch (by dtype)
├── ❌ Output allocation
└── ❌ Broadcasting integration
```

### 4.2 Math Ufuncs - Unary ❌
```
Unary Math
├── Arithmetic
│   ├── ❌ negative(x), positive(x)
│   ├── ❌ absolute(x), fabs(x)
│   ├── ❌ sign(x)
│   ├── ❌ sqrt(x), square(x), cbrt(x)
│   └── ❌ reciprocal(x)
├── Exponents & Logarithms
│   ├── ❌ exp(x), exp2(x), expm1(x)
│   ├── ❌ log(x), log2(x), log10(x), log1p(x)
│   └── ❌ logaddexp(x1, x2), logaddexp2(x1, x2)
├── Trigonometric
│   ├── ❌ sin(x), cos(x), tan(x)
│   ├── ❌ arcsin(x), arccos(x), arctan(x)
│   ├── ❌ degrees(x), radians(x)
│   └── ❌ deg2rad(x), rad2deg(x)
├── Hyperbolic
│   ├── ❌ sinh(x), cosh(x), tanh(x)
│   └── ❌ arcsinh(x), arccosh(x), arctanh(x)
├── Rounding
│   ├── ❌ round(x, decimals), around(x, decimals)
│   ├── ❌ rint(x), fix(x)
│   └── ❌ floor(x), ceil(x), trunc(x)
├── Floating Point
│   ├── ❌ signbit(x), copysign(x1, x2)
│   ├── ❌ frexp(x), ldexp(x1, x2)
│   ├── ❌ nextafter(x1, x2), spacing(x)
│   └── ❌ nan_to_num(x, nan, posinf, neginf)
├── Special
│   ├── ❌ i0(x) → Bessel
│   ├── ❌ sinc(x)
│   └── ❌ heaviside(x1, x2)
└── Predicates
    ├── ❌ isnan(x), isinf(x), isfinite(x)
    ├── ❌ isneginf(x), isposinf(x)
    └── ❌ isnat(x)
```

### 4.3 Math Ufuncs - Binary ❌
```
Binary Math
├── Arithmetic
│   ├── ❌ add(x1, x2), subtract(x1, x2)
│   ├── ❌ multiply(x1, x2), divide(x1, x2)
│   ├── ❌ true_divide(x1, x2), floor_divide(x1, x2)
│   ├── ❌ power(x1, x2), float_power(x1, x2)
│   ├── ❌ mod(x1, x2), remainder(x1, x2), fmod(x1, x2)
│   └── ❌ divmod(x1, x2), modf(x)
├── Comparison
│   ├── ❌ greater(x1, x2), greater_equal(x1, x2)
│   ├── ❌ less(x1, x2), less_equal(x1, x2)
│   ├── ❌ equal(x1, x2), not_equal(x1, x2)
│   └── ❌ maximum(x1, x2), minimum(x1, x2)
│       fmax(x1, x2), fmin(x1, x2)
├── Logical
│   ├── ❌ logical_and(x1, x2)
│   ├── ❌ logical_or(x1, x2)
│   ├── ❌ logical_xor(x1, x2)
│   └── ❌ logical_not(x)
├── Bitwise
│   ├── ❌ bitwise_and(x1, x2), bitwise_or(x1, x2)
│   ├── ❌ bitwise_xor(x1, x2), invert(x)
│   ├── ❌ left_shift(x1, x2), right_shift(x1, x2)
│   └── ❌ bitwise_count(x)
├── Trigonometric
│   ├── ❌ arctan2(x1, x2)
│   └── ❌ hypot(x1, x2)
└── Rational
    ├── ❌ gcd(x1, x2)
    └── ❌ lcm(x1, x2)
```

### 4.4 Reductions (with axis support) 🔶 PARTIAL
```
Reductions
├── 🔶 sum(arr) - implemented WITHOUT axis support, uses pairwise summation
├── ❌ prod(arr, axis, dtype, out, keepdims, initial)
├── ❌ nansum(arr, axis, ...), nanprod(arr, axis, ...)
├── ❌ cumsum(arr, axis, dtype, out), cumprod(arr, axis, dtype, out)
├── ❌ nancumsum(arr, ...), nancumprod(arr, ...)
├── ❌ diff(arr, n, axis, prepend, append)
├── ❌ ediff1d(arr, to_end, to_begin)
├── ❌ gradient(f, *varargs, axis, edge_order)
├── ❌ cross(a, b, axisa, axisb, axisc, axis)
├── ❌ trapezoid(y, x, dx, axis)
└── Aggregations
    ├── ❌ min(arr, axis, ...), max(arr, axis, ...)
    ├── ❌ amin(arr, ...), amax(arr, ...)
    ├── ❌ nanmin(arr, ...), nanmax(arr, ...)
    ├── ❌ ptp(arr, axis) → peak to peak
    ├── 🔶 argmin(arr), argmax(arr) - TypeScript only, limited axis support
    ├── ❌ nanargmin(arr, ...), nanargmax(arr, ...)
    ├── ❌ all(arr, axis), any(arr, axis)
    └── ✅ countNonzero(arr, axis)
```

**Existing Implementation:** `src/wasm/pairwise_sum.c` (accurate summation algorithm)

### 4.5 Complex Numbers ❌
```
Complex
├── ❌ real(x), imag(x)
├── ❌ conj(x), conjugate(x)
└── ❌ angle(z, deg)
```

---

## Phase 5: Array Manipulation ❌ NOT STARTED

### 5.1 Joining Arrays ❌
```
Joining
├── ❌ concatenate(arrays, axis, out, dtype, casting)
├── ❌ stack(arrays, axis, out, dtype, casting)
├── ❌ vstack(tup), row_stack(tup)
├── ❌ hstack(tup)
├── ❌ dstack(tup)
├── ❌ column_stack(tup)
├── ❌ block(arrays)
└── ❌ append(arr, values, axis)
```

### 5.2 Splitting Arrays ❌
```
Splitting
├── ❌ split(arr, indices_or_sections, axis)
├── ❌ array_split(arr, indices_or_sections, axis)
├── ❌ vsplit(arr, indices_or_sections)
├── ❌ hsplit(arr, indices_or_sections)
├── ❌ dsplit(arr, indices_or_sections)
└── ❌ unstack(x, axis)
```

### 5.3 Tiling & Repeating ❌
```
Tiling
├── ❌ tile(arr, reps)
├── ❌ repeat(arr, repeats, axis)
└── ❌ pad(arr, pad_width, mode, ...)
```

### 5.4 Rearranging ❌
```
Rearranging
├── ❌ flip(arr, axis)
├── ❌ fliplr(arr), flipud(arr)
├── ❌ roll(arr, shift, axis)
├── ❌ rot90(arr, k, axes)
├── ❌ resize(arr, new_shape)
├── ❌ trim_zeros(filt, trim)
├── ❌ insert(arr, obj, values, axis)
└── ❌ delete(arr, obj, axis)
```

### 5.5 Copying 🔶 PARTIAL
```
Copying
├── ✅ copy(a)
├── ❌ copyto(dst, src, casting, where)
└── ❌ asarray(a, dtype, order, ...)
```

---

## Phase 6: Sorting, Searching & Statistics 🔶 PARTIAL

### 6.1 Sorting 🔶 PARTIAL
```
Sorting
├── ❌ sort(arr, axis, kind, order)
├── 🔶 argsort(arr) - TypeScript only, 1D arrays only
├── ❌ lexsort(keys, axis)
├── ❌ sort_complex(arr)
├── ❌ partition(arr, kth, axis, kind, order)
├── ❌ argpartition(arr, kth, axis, kind, order)
└── ❌ msort(arr) → sort along first axis
```

### 6.2 Searching 🔶 PARTIAL
```
Searching
├── 🔶 argmax(arr) - TypeScript only, limited axis support
├── 🔶 argmin(arr) - TypeScript only, limited axis support
├── ❌ nanargmax(arr, axis, out, keepdims)
├── ❌ nanargmin(arr, axis, out, keepdims)
├── ✅ nonzero(arr)
├── ✅ flatnonzero(arr)
├── ❌ argwhere(arr)
├── ✅ where(condition, x, y)
├── ❌ searchsorted(a, v, side, sorter)
└── ✅ extract(condition, arr)
```

### 6.3 Statistics ❌ NOT STARTED
```
Statistics
├── Averages & Variances
│   ├── ❌ mean(arr, axis, dtype, out, keepdims)
│   ├── ❌ average(arr, axis, weights, returned, keepdims)
│   ├── ❌ std(arr, axis, dtype, out, ddof, keepdims)
│   ├── ❌ var(arr, axis, dtype, out, ddof, keepdims)
│   ├── ❌ nanmean(...), nanstd(...), nanvar(...)
│   └── ❌ median(arr, axis, out, overwrite_input, keepdims)
│       nanmedian(...)
├── Order Statistics
│   ├── ❌ amin(arr, ...), amax(arr, ...)
│   ├── ❌ ptp(arr, axis, out, keepdims)
│   ├── ❌ percentile(arr, q, axis, out, ...)
│   ├── ❌ quantile(arr, q, axis, out, ...)
│   └── ❌ nanpercentile(...), nanquantile(...)
├── Correlating
│   ├── ❌ corrcoef(x, y, rowvar, bias, ddof, dtype)
│   ├── ❌ correlate(a, v, mode)
│   └── ❌ cov(m, y, rowvar, bias, ddof, fweights, aweights, dtype)
└── Histograms
    ├── ❌ histogram(a, bins, range, density, weights)
    ├── ❌ histogram2d(x, y, bins, range, density, weights)
    ├── ❌ histogramdd(sample, bins, range, density, weights)
    ├── ❌ histogram_bin_edges(a, bins, range, weights)
    ├── ❌ bincount(x, weights, minlength)
    └── ❌ digitize(x, bins, right)
```

---

## Phase 7: Logic & Comparison ❌ NOT STARTED

### 7.1 Truth Testing ❌
```
Truth Testing
├── ❌ all(a, axis, out, keepdims, where)
├── ❌ any(a, axis, out, keepdims, where)
├── ❌ allclose(a, b, rtol, atol, equal_nan)
├── ❌ isclose(a, b, rtol, atol, equal_nan)
├── ❌ array_equal(a1, a2, equal_nan)
└── ❌ array_equiv(a1, a2)
```

### 7.2 Array Contents ❌
```
Array Contents
├── ❌ isfinite(x), isinf(x), isnan(x)
├── ❌ isnat(x)
├── ❌ isneginf(x), isposinf(x)
├── ❌ iscomplex(x), iscomplexobj(x)
├── ❌ isreal(x), isrealobj(x)
├── ❌ isfortran(a)
└── ❌ isscalar(element)
```

---

## Phase 8: Set Operations ❌ NOT STARTED

```
Set Operations
├── ❌ unique(ar, return_index, return_inverse, return_counts, axis, equal_nan)
├── ❌ unique_all(x), unique_counts(x)
├── ❌ unique_inverse(x), unique_values(x)
├── ❌ in1d(ar1, ar2, assume_unique, invert, kind)
├── ❌ isin(element, test_elements, assume_unique, invert, kind)
├── ❌ intersect1d(ar1, ar2, assume_unique, return_indices)
├── ❌ setdiff1d(ar1, ar2, assume_unique)
├── ❌ setxor1d(ar1, ar2, assume_unique)
└── ❌ union1d(ar1, ar2)
```

---

## Phase 9: I/O Operations ❌ NOT STARTED

```
Input/Output
├── Binary Files
│   ├── ❌ save(file, arr, allow_pickle, fix_imports)
│   ├── ❌ load(file, mmap_mode, allow_pickle, fix_imports, encoding)
│   ├── ❌ savez(file, *args, **kwds)
│   └── ❌ savez_compressed(file, *args, **kwds)
├── Text Files
│   ├── ❌ loadtxt(fname, dtype, comments, delimiter, ...)
│   ├── ❌ savetxt(fname, X, fmt, delimiter, newline, ...)
│   ├── ❌ genfromtxt(fname, dtype, comments, delimiter, ...)
│   └── ❌ fromregex(file, regexp, dtype, encoding)
├── Raw Binary
│   ├── ❌ fromfile(file, dtype, count, sep, offset, like)
│   └── ❌ tofile(fid, sep, format)
├── String Formatting
│   ├── ❌ array2string(a, max_line_width, precision, ...)
│   ├── ❌ array_repr(arr, max_line_width, precision, ...)
│   ├── ❌ array_str(a, max_line_width, precision, ...)
│   ├── ❌ format_float_positional(x, precision, ...)
│   └── ❌ format_float_scientific(x, precision, ...)
├── Memory Mapping
│   └── ❌ memmap(filename, dtype, mode, offset, shape, order)
├── Print Options
│   ├── ❌ set_printoptions(precision, threshold, edgeitems, ...)
│   ├── ❌ get_printoptions()
│   └── ❌ printoptions(*args, **kwargs)
└── Base Conversion
    ├── ❌ binary_repr(num, width)
    └── ❌ base_repr(number, base, padding)
```

---

## Phase 10: Functional Programming ❌ NOT STARTED

```
Functional
├── ❌ apply_along_axis(func1d, axis, arr, *args, **kwargs)
├── ❌ apply_over_axes(func, a, axes)
├── ❌ vectorize(pyfunc, otypes, doc, excluded, cache, signature)
├── ❌ frompyfunc(func, nin, nout, identity)
└── ❌ piecewise(x, condlist, funclist, *args, **kw)
```

---

## Phase 11: Window Functions ❌ NOT STARTED

```
Window Functions
├── ❌ bartlett(M)
├── ❌ blackman(M)
├── ❌ hamming(M)
├── ❌ hanning(M)
└── ❌ kaiser(M, beta)
```

---

## Phase 12: Constants ❌ NOT STARTED

```
Constants
├── ❌ e → 2.71828...
├── ❌ euler_gamma → 0.57721...
├── ❌ inf → positive infinity
├── ❌ nan → Not a Number
├── ❌ newaxis → None (for indexing)
└── ❌ pi → 3.14159...
```

---

## Phase 13: numpy.linalg ❌ NOT STARTED

```
numpy.linalg
├── Matrix Products
│   ├── ❌ dot(a, b, out)
│   ├── ❌ vdot(a, b)
│   ├── ❌ inner(a, b)
│   ├── ❌ outer(a, b, out)
│   ├── ❌ matmul(x1, x2, out) / @ operator
│   ├── ❌ tensordot(a, b, axes)
│   ├── ❌ einsum(subscripts, *operands, out, ...)
│   ├── ❌ einsum_path(subscripts, *operands, optimize)
│   ├── ❌ kron(a, b)
│   ├── ❌ cross(a, b, axisa, axisb, axisc, axis)
│   └── ❌ multi_dot(arrays, out)
├── Decompositions
│   ├── ❌ cholesky(a)
│   ├── ❌ qr(a, mode)
│   ├── ❌ svd(a, full_matrices, compute_uv, hermitian)
│   └── ❌ svdvals(x)
├── Eigenvalues
│   ├── ❌ eig(a)
│   ├── ❌ eigh(a, UPLO)
│   ├── ❌ eigvals(a)
│   └── ❌ eigvalsh(a, UPLO)
├── Norms & Numbers
│   ├── ❌ norm(x, ord, axis, keepdims)
│   ├── ❌ matrix_norm(x, ord, keepdims)
│   ├── ❌ vector_norm(x, ord, axis, keepdims)
│   ├── ❌ cond(x, p)
│   ├── ❌ det(a)
│   ├── ❌ slogdet(a)
│   ├── ❌ matrix_rank(A, tol, hermitian, rtol)
│   └── ❌ trace(a, offset, axis1, axis2, dtype, out)
├── Solving & Inverting
│   ├── ❌ solve(a, b)
│   ├── ❌ tensorsolve(a, b, axes)
│   ├── ❌ lstsq(a, b, rcond)
│   ├── ❌ inv(a)
│   ├── ❌ pinv(a, rcond, hermitian, rtol)
│   └── ❌ tensorinv(a, ind)
├── Matrix Operations
│   ├── ❌ matrix_power(a, n)
│   ├── ✅ diagonal(a, offset, axis1, axis2)
│   └── ❌ matrix_transpose(x)
└── Exception
    └── ❌ LinAlgError
```

---

## Phase 14: numpy.fft ❌ NOT STARTED

```
numpy.fft
├── Standard FFTs
│   ├── ❌ fft(a, n, axis, norm, out)
│   ├── ❌ ifft(a, n, axis, norm, out)
│   ├── ❌ fft2(a, s, axes, norm, out)
│   ├── ❌ ifft2(a, s, axes, norm, out)
│   ├── ❌ fftn(a, s, axes, norm, out)
│   └── ❌ ifftn(a, s, axes, norm, out)
├── Real FFTs
│   ├── ❌ rfft(a, n, axis, norm, out)
│   ├── ❌ irfft(a, n, axis, norm, out)
│   ├── ❌ rfft2(a, s, axes, norm, out)
│   ├── ❌ irfft2(a, s, axes, norm, out)
│   ├── ❌ rfftn(a, s, axes, norm, out)
│   └── ❌ irfftn(a, s, axes, norm, out)
├── Hermitian FFTs
│   ├── ❌ hfft(a, n, axis, norm, out)
│   └── ❌ ihfft(a, n, axis, norm, out)
└── Helper Functions
    ├── ❌ fftfreq(n, d, device)
    ├── ❌ rfftfreq(n, d, device)
    ├── ❌ fftshift(x, axes)
    └── ❌ ifftshift(x, axes)
```

---

## Phase 15: numpy.random ❌ NOT STARTED

```
numpy.random
├── Generator Class
│   ├── ❌ default_rng(seed) → Generator
│   └── Generator Methods
│       ├── ❌ random(size, dtype, out)
│       ├── ❌ integers(low, high, size, dtype, endpoint)
│       ├── ❌ uniform(low, high, size)
│       ├── ❌ normal(loc, scale, size)
│       ├── ❌ standard_normal(size, dtype, out)
│       ├── ❌ exponential(scale, size)
│       ├── ❌ poisson(lam, size)
│       ├── ❌ binomial(n, p, size)
│       ├── ❌ beta(a, b, size)
│       ├── ❌ gamma(shape, scale, size)
│       ├── ❌ chisquare(df, size)
│       ├── ❌ choice(a, size, replace, p, axis, shuffle)
│       ├── ❌ shuffle(x, axis)
│       ├── ❌ permutation(x, axis)
│       └── ... (many more distributions)
├── BitGenerator Infrastructure
│   ├── ❌ PCG64 (default)
│   ├── ❌ MT19937
│   ├── ❌ Philox
│   └── ❌ SFC64
├── SeedSequence
│   └── ❌ SeedSequence(entropy, spawn_key, pool_size)
└── Legacy (RandomState)
    └── ❌ Backward compatibility functions
```

---

## Phase 16: Additional Modules ❌ NOT STARTED

### numpy.strings (2.0+) ❌
### numpy.polynomial ❌
### numpy.ma (Masked Arrays) ❌
### numpy.rec (Record Arrays) ❌
### numpy.testing ❌

---

## Phase 17: Error Handling & Configuration ❌ NOT STARTED

### Error Handling ❌
### Exceptions ❌

---

## Current Implementation Summary

### TypeScript Files (`src/ts/`)
| File | Lines | Description |
|------|-------|-------------|
| `NDArray.ts` | ~1,527 | Core NDArray class with all methods |
| `types.ts` | ~400 | Type definitions and DType system |
| `dtype.ts` | ~200 | DType utilities and conversion |
| `broadcast.ts` | ~150 | Broadcasting functions |
| `indexing.ts` | ~350 | Index operations |
| `slice.ts` | ~200 | Slicing utilities |
| `iterators.ts` | ~150 | Iterator implementations |
| `wasm-loader.ts` | ~100 | WASM module management |
| `index.ts` | ~50 | Main exports |

### C/WASM Files (`src/wasm/`)
| File | Lines | Description |
|------|-------|-------------|
| `ndarray.c` | ~1,500 | Core array operations |
| `dtype.c` | ~318 | Type system |
| `broadcast.c` | ~212 | Broadcasting |
| `indexing.c` | ~626 | Index operations |
| `pairwise_sum.c` | ~151 | Accurate summation algorithm |

**Total: ~5,300 lines of implementation code**

---

## Implementation Priority Summary

```
CRITICAL PATH (Enables Everything Else):
1. ✅ Extended DTypes + Type Promotion
2. ✅ Element Access (get/set)
3. ✅ Iterators
4. ✅ Views + reshape/transpose
5. ✅ Slicing (basic)
6. ✅ Broadcasting
7. ❌ Ufunc Infrastructure          ← NEXT PRIORITY
8. ❌ Core Ufuncs (add, subtract, multiply, divide, comparison)
9. ❌ Reductions with axis

HIGH VALUE:
10. ❌ Array manipulation (concat, stack, split)
11. ❌ Sorting & Searching
12. ❌ Statistics (mean, std, var, median)
13. ❌ numpy.linalg (dot, matmul, solve, inv)
14. ❌ numpy.random (Generator, basic distributions)

MEDIUM VALUE:
15. ❌ numpy.fft
16. ❌ Advanced indexing (fancy, boolean)
17. ❌ Set operations
18. ❌ Window functions
19. ❌ I/O operations

LOWER PRIORITY:
20. ❌ numpy.ma (masked arrays)
21. ❌ numpy.polynomial
22. ❌ numpy.strings
23. ❌ numpy.rec
24. ❌ numpy.testing
```

---

## Verification

After each phase:
1. Run existing tests: `npm test`
2. Add new tests for implemented features
3. Run benchmarks: `npm run benchmark`
4. Validate against NumPy test vectors (extend `tests/python/generate_test_cases.py`)
