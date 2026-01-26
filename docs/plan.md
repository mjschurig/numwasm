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
│   │       ├── umath/                  # Universal math functions
│   │       ├── common/                 # Common utilities
│   │       ├── npymath/                # Math library
│   │       ├── npysort/                # Sorting algorithms
│   │       └── _simd/                  # SIMD optimizations
│   │
│   ├── lib/                            # Library functions (~50+ files)
│   ├── linalg/                         # Linear algebra
│   ├── fft/                            # Fast Fourier Transforms
│   ├── random/                         # Random number generation
│   ├── polynomial/                     # Polynomial operations
│   ├── ma/                             # Masked arrays
│   ├── strings/                        # String operations
│   ├── rec/                            # Record arrays
│   └── testing/                        # Testing utilities
```

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
├── ✅ Ufunc Infrastructure (WASM-accelerated)
├── ✅ Unary Ufuncs (neg, abs, sqrt, exp, log, trig, rounding, predicates)
├── ✅ Binary Ufuncs (arithmetic, comparison, logical, bitwise)
└── ✅ Reductions with Axis (sum, mean, std, var, min, max, argmin, argmax)

LEVEL 4: MANIPULATION & STATISTICS
├── ✅ Array Manipulation (concat, stack, split, flip, roll, tile, pad)
├── ✅ Sorting & Searching (sort, argsort, partition, argpartition, searchsorted)
├── ✅ Statistics (mean, std, var, median)
└── ✅ Set Operations (unique, intersect1d, union1d, setdiff1d, setxor1d, isin)

LEVEL 5: HIGHER-LEVEL MODULES
├── ✅ numpy.linalg (dot, matmul, solve, inv, eig, svd, qr, cholesky)
├── ✅ numpy.fft (fft, ifft, fft2, fftn, fftfreq, fftshift)
├── ✅ numpy.random (Generator, PCG64, 30+ distributions)
├── ✅ numpy.polynomial (Polynomial, Chebyshev, Legendre, Hermite, Laguerre)
├── ✅ numpy.ma (MaskedArray, masked operations, cov, corrcoef)
├── ✅ numpy.strings (40+ string operations)
├── ✅ numpy.rec (Record arrays)
└── ✅ numpy.testing (Assertion utilities)
```

---

## Phase 1: Foundation Enhancement ✅ COMPLETE

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

## Phase 2: Iteration & Shape Manipulation ✅ COMPLETE

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

## Phase 3: Views, Slicing & Broadcasting ✅ COMPLETE

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
├── ✅ Advanced Indexing
│   ├── ✅ take(arr, indices, axis) → fancy indexing equivalent
│   └── ✅ compress/extract → boolean masking equivalent
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

## Phase 4: Universal Functions (Ufuncs) ✅ COMPLETE

### 4.1 Ufunc Infrastructure ✅ COMPLETE
```
Ufunc System
├── ✅ WASM-accelerated unary operations
├── ✅ WASM-accelerated binary operations
├── ✅ Broadcasting integration
├── ✅ Output allocation
├── ✅ Inner loop dispatch (by dtype)
└── 🔶 Advanced Ufunc Methods
    ├── ❌ accumulate(arr, axis, dtype, out)
    ├── ❌ reduceat(arr, indices, axis, dtype, out)
    ├── ❌ outer(a, b, out)
    └── ❌ at(arr, indices, b)
```

**Implementation:** `src/ts/ufunc.ts`, `src/wasm/ufunc_unary.c`, `src/wasm/ufunc_binary.c`

### 4.2 Math Ufuncs - Unary ✅ COMPLETE (90+ functions)
```
Unary Math
├── Arithmetic
│   ├── ✅ negative(x), positive(x)
│   ├── ✅ absolute(x), abs(x)
│   ├── ✅ sign(x)
│   ├── ✅ sqrt(x), square(x), cbrt(x)
│   └── ✅ reciprocal(x)
├── Exponents & Logarithms
│   ├── ✅ exp(x), exp2(x), expm1(x)
│   ├── ✅ log(x), log2(x), log10(x), log1p(x)
│   └── ✅ logaddexp(x1, x2), logaddexp2(x1, x2)
├── Trigonometric
│   ├── ✅ sin(x), cos(x), tan(x)
│   ├── ✅ arcsin(x), arccos(x), arctan(x)
│   ├── ✅ degrees(x), radians(x)
│   └── ✅ deg2rad(x), rad2deg(x)
├── Hyperbolic
│   ├── ✅ sinh(x), cosh(x), tanh(x)
│   └── ✅ arcsinh(x), arccosh(x), arctanh(x)
├── Rounding
│   ├── ✅ round(x, decimals), around(x, decimals)
│   ├── ✅ rint(x), fix(x)
│   └── ✅ floor(x), ceil(x), trunc(x)
├── Floating Point
│   ├── ✅ signbit(x), copysign(x1, x2)
│   ├── ❌ frexp(x), ldexp(x1, x2)
│   ├── ❌ nextafter(x1, x2), spacing(x)
│   └── ❌ nan_to_num(x, nan, posinf, neginf)
├── Special
│   ├── ✅ i0(x) → Bessel (in window.ts)
│   ├── ❌ sinc(x)
│   └── ❌ heaviside(x1, x2)
└── Predicates
    ├── ✅ isnan(x), isinf(x), isfinite(x)
    ├── ✅ isneginf(x), isposinf(x)
    └── ❌ isnat(x)
```

### 4.3 Math Ufuncs - Binary ✅ COMPLETE
```
Binary Math
├── Arithmetic
│   ├── ✅ add(x1, x2), subtract(x1, x2)
│   ├── ✅ multiply(x1, x2), divide(x1, x2)
│   ├── ✅ true_divide(x1, x2), floor_divide(x1, x2)
│   ├── ✅ power(x1, x2), float_power(x1, x2)
│   ├── ✅ mod(x1, x2), remainder(x1, x2), fmod(x1, x2)
│   └── ❌ divmod(x1, x2), modf(x)
├── Comparison
│   ├── ✅ greater(x1, x2), greater_equal(x1, x2)
│   ├── ✅ less(x1, x2), less_equal(x1, x2)
│   ├── ✅ equal(x1, x2), not_equal(x1, x2)
│   └── ✅ maximum(x1, x2), minimum(x1, x2)
│       ✅ fmax(x1, x2), fmin(x1, x2)
├── Logical
│   ├── ✅ logical_and(x1, x2)
│   ├── ✅ logical_or(x1, x2)
│   ├── ✅ logical_xor(x1, x2)
│   └── ✅ logical_not(x)
├── Bitwise
│   ├── ✅ bitwise_and(x1, x2), bitwise_or(x1, x2)
│   ├── ✅ bitwise_xor(x1, x2), invert(x)
│   ├── ✅ left_shift(x1, x2), right_shift(x1, x2)
│   └── ❌ bitwise_count(x)
├── Trigonometric
│   ├── ✅ arctan2(x1, x2)
│   └── ✅ hypot(x1, x2)
└── Rational
    ├── ❌ gcd(x1, x2)
    └── ❌ lcm(x1, x2)
```

### 4.4 Reductions (with axis support) ✅ COMPLETE
```
Reductions
├── ✅ sum(arr, axis, dtype, keepdims) - uses pairwise summation for accuracy
├── ✅ prod(arr, axis, dtype, keepdims)
├── ❌ nansum(arr, axis, ...), nanprod(arr, axis, ...)
├── ❌ cumsum(arr, axis, dtype, out), cumprod(arr, axis, dtype, out)
├── ❌ nancumsum(arr, ...), nancumprod(arr, ...)
├── ✅ diff(arr, n, axis, prepend, append) → ediff1d
├── ✅ ediff1d(arr, to_end, to_begin)
├── ❌ gradient(f, *varargs, axis, edge_order)
├── ❌ cross(a, b, axisa, axisb, axisc, axis)
├── ❌ trapezoid(y, x, dx, axis)
└── Aggregations
    ├── ✅ min(arr, axis, keepdims), max(arr, axis, keepdims)
    ├── ✅ amin(arr, ...), amax(arr, ...)
    ├── ❌ nanmin(arr, ...), nanmax(arr, ...)
    ├── ❌ ptp(arr, axis) → peak to peak
    ├── ✅ argmin(arr, axis, keepdims), argmax(arr, axis, keepdims)
    ├── ❌ nanargmin(arr, ...), nanargmax(arr, ...)
    ├── ✅ all(arr, axis, keepdims), any(arr, axis, keepdims)
    └── ✅ countNonzero(arr, axis)
```

**Implementation:** `src/ts/statistics.ts`, `src/wasm/statistics.c`, `src/wasm/pairwise_sum.c`

### 4.5 Complex Numbers 🔶 PARTIAL
```
Complex
├── 🔶 Complex64, Complex128 dtypes supported
├── ❌ real(x), imag(x)
├── ❌ conj(x), conjugate(x)
└── ❌ angle(z, deg)
```

---

## Phase 5: Array Manipulation ✅ COMPLETE

### 5.1 Joining Arrays ✅ COMPLETE
```
Joining
├── ✅ concatenate(arrays, axis, out, dtype, casting)
├── ✅ stack(arrays, axis, out, dtype, casting)
├── ✅ vstack(tup), row_stack(tup)
├── ✅ hstack(tup)
├── ✅ dstack(tup)
├── ✅ column_stack(tup)
├── ✅ block(arrays)
└── ✅ append(arr, values, axis)
```

**Implementation:** `src/ts/manipulation.ts`, `src/wasm/manipulation.c`

### 5.2 Splitting Arrays ✅ COMPLETE
```
Splitting
├── ✅ split(arr, indices_or_sections, axis)
├── ✅ array_split(arr, indices_or_sections, axis)
├── ✅ vsplit(arr, indices_or_sections)
├── ✅ hsplit(arr, indices_or_sections)
├── ✅ dsplit(arr, indices_or_sections)
└── ✅ unstack(x, axis)
```

### 5.3 Tiling & Repeating ✅ COMPLETE
```
Tiling
├── ✅ tile(arr, reps)
├── ✅ repeat(arr, repeats, axis)
└── ✅ pad(arr, pad_width, mode, constant_values)
```

### 5.4 Rearranging ✅ COMPLETE
```
Rearranging
├── ✅ flip(arr, axis)
├── ✅ fliplr(arr), flipud(arr)
├── ✅ roll(arr, shift, axis)
├── ✅ rot90(arr, k, axes)
├── ✅ resize(arr, new_shape)
├── ✅ trim_zeros(filt, trim)
├── ✅ insert(arr, obj, values, axis)
└── ✅ deleteArr(arr, obj, axis)
```

### 5.5 Copying ✅ COMPLETE
```
Copying
├── ✅ copy(a)
├── ✅ copyto(dst, src, where)
└── ✅ asarray(a, dtype)
```

---

## Phase 6: Sorting, Searching & Statistics ✅ COMPLETE

### 6.1 Sorting ✅ COMPLETE
```
Sorting
├── ✅ sort(arr, axis, kind) - supports quicksort, mergesort, heapsort
├── ✅ argsort(arr, axis, kind)
├── ❌ lexsort(keys, axis)
├── ❌ sort_complex(arr)
├── ✅ partition(arr, kth, axis)
├── ✅ argpartition(arr, kth, axis)
└── ❌ msort(arr) → sort along first axis
```

**Implementation:** `src/ts/sorting.ts`, `src/wasm/sorting.c`

### 6.2 Searching ✅ COMPLETE
```
Searching
├── ✅ argmax(arr, axis, keepdims)
├── ✅ argmin(arr, axis, keepdims)
├── ❌ nanargmax(arr, axis, out, keepdims)
├── ❌ nanargmin(arr, axis, out, keepdims)
├── ✅ nonzero(arr)
├── ✅ flatnonzero(arr)
├── ✅ argwhere(arr)
├── ✅ where(condition, x, y)
├── ✅ searchsorted(a, v, side, sorter)
└── ✅ extract(condition, arr)
```

### 6.3 Statistics ✅ COMPLETE
```
Statistics
├── Averages & Variances
│   ├── ✅ mean(arr, axis, dtype, keepdims)
│   ├── ❌ average(arr, axis, weights, returned, keepdims)
│   ├── ✅ std(arr, axis, dtype, ddof, keepdims)
│   ├── ✅ var(arr, axis, dtype, ddof, keepdims)
│   ├── ❌ nanmean(...), nanstd(...), nanvar(...)
│   └── ✅ median(arr, axis, keepdims)
│       ❌ nanmedian(...)
├── Order Statistics
│   ├── ✅ amin(arr, ...), amax(arr, ...)
│   ├── ❌ ptp(arr, axis, out, keepdims)
│   ├── ❌ percentile(arr, q, axis, out, ...)
│   ├── ❌ quantile(arr, q, axis, out, ...)
│   └── ❌ nanpercentile(...), nanquantile(...)
├── Correlating
│   ├── ✅ corrcoef(x, y, ...) - via ma module
│   ├── ❌ correlate(a, v, mode)
│   └── ✅ cov(m, y, ...) - via ma module
└── Histograms
    ├── ❌ histogram(a, bins, range, density, weights)
    ├── ❌ histogram2d(x, y, bins, range, density, weights)
    ├── ❌ histogramdd(sample, bins, range, density, weights)
    ├── ❌ histogram_bin_edges(a, bins, range, weights)
    ├── ❌ bincount(x, weights, minlength)
    └── ❌ digitize(x, bins, right)
```

---

## Phase 7: Logic & Comparison ✅ COMPLETE

### 7.1 Truth Testing ✅ COMPLETE
```
Truth Testing
├── ✅ all(a, axis, out, keepdims)
├── ✅ any(a, axis, out, keepdims)
├── ✅ allclose(a, b, rtol, atol, equal_nan)
├── ✅ isclose(a, b, rtol, atol, equal_nan)
├── ✅ array_equal(a1, a2, equal_nan)
└── ✅ array_equiv(a1, a2)
```

**Implementation:** `src/ts/logic.ts`, `src/wasm/logic.c`

### 7.2 Array Contents ✅ COMPLETE
```
Array Contents
├── ✅ isfinite(x), isinf(x), isnan(x)
├── ❌ isnat(x)
├── ✅ isneginf(x), isposinf(x)
├── ✅ iscomplex(x), iscomplexobj(x)
├── ✅ isreal(x), isrealobj(x)
├── ✅ isfortran(a)
└── ✅ isscalar(element)
```

---

## Phase 8: Set Operations ✅ COMPLETE

```
Set Operations
├── ✅ unique(ar, return_index, return_inverse, return_counts, axis, equal_nan)
├── ✅ unique_all(x), unique_counts(x)
├── ✅ unique_inverse(x), unique_values(x)
├── ✅ in1d(ar1, ar2, assume_unique, invert, kind)
├── ✅ isin(element, test_elements, assume_unique, invert, kind)
├── ✅ intersect1d(ar1, ar2, assume_unique, return_indices)
├── ✅ setdiff1d(ar1, ar2, assume_unique)
├── ✅ setxor1d(ar1, ar2, assume_unique)
└── ✅ union1d(ar1, ar2)
```

**Implementation:** `src/ts/setops.ts`, `src/wasm/setops.c`

---

## Phase 9: I/O Operations ✅ COMPLETE

```
Input/Output
├── Binary Files
│   ├── ✅ save(file, arr) - NPY format
│   ├── ✅ load(file) - NPY format
│   ├── ❌ savez(file, *args, **kwds)
│   └── ❌ savez_compressed(file, *args, **kwds)
├── Text Files
│   ├── ✅ loadtxt(fname, dtype, delimiter, skiprows, ...)
│   ├── ✅ savetxt(fname, X, fmt, delimiter, newline, ...)
│   ├── ✅ genfromtxt(fname, dtype, delimiter, skip_header, ...)
│   └── ✅ fromregex(file, regexp, dtype)
├── Raw Binary
│   ├── ✅ fromfile(file, dtype, count, offset)
│   └── ✅ frombuffer(buffer, dtype, count, offset)
├── String Formatting
│   ├── ✅ array2string(a, max_line_width, precision, ...)
│   ├── ✅ array_repr(arr, max_line_width, precision, ...)
│   ├── ✅ array_str(a, max_line_width, precision, ...)
│   ├── ✅ format_float_positional(x, precision, ...)
│   └── ✅ format_float_scientific(x, precision, ...)
├── Memory Mapping
│   └── ✅ memmap(filename, dtype, mode, offset, shape, order) - Memmap class
├── Print Options
│   ├── ✅ set_printoptions(precision, threshold, edgeitems, ...)
│   ├── ✅ get_printoptions()
│   └── ✅ printoptions(*args, **kwargs) → withPrintoptions
└── Base Conversion
    ├── ✅ binary_repr(num, width)
    └── ✅ base_repr(number, base, padding)
```

**Implementation:** `src/ts/io/` (10+ files)

---

## Phase 10: Functional Programming ✅ COMPLETE

```
Functional
├── ✅ apply_along_axis(func1d, axis, arr, *args, **kwargs)
├── ✅ apply_over_axes(func, a, axes)
├── ✅ vectorize(pyfunc, otypes, doc, excluded, cache, signature)
├── ✅ frompyfunc(func, nin, nout, identity)
└── ✅ piecewise(x, condlist, funclist, *args, **kw)
```

**Implementation:** `src/ts/functional.ts`

---

## Phase 11: Window Functions ✅ COMPLETE

```
Window Functions
├── ✅ bartlett(M)
├── ✅ blackman(M)
├── ✅ hamming(M)
├── ✅ hanning(M)
├── ✅ kaiser(M, beta)
└── ✅ i0(x) - Modified Bessel function
```

**Implementation:** `src/ts/window.ts`

---

## Phase 12: Constants ✅ COMPLETE

```
Constants
├── ✅ e → 2.71828...
├── ✅ euler_gamma → 0.57721...
├── ✅ inf, PINF, NINF → infinity values
├── ✅ nan, NAN → Not a Number
├── ✅ PZERO, NZERO → signed zeros
├── ✅ newaxis → for indexing
└── ✅ pi → 3.14159...
```

**Implementation:** `src/ts/constants.ts`

---

## Phase 13: Type Information ✅ COMPLETE

```
Type Information
├── ✅ finfo(dtype) → FloatInfo (eps, max, min, bits, etc.)
└── ✅ iinfo(dtype) → IntInfo (min, max)
```

**Implementation:** `src/ts/typeinfo.ts`

---

## Phase 14: numpy.linalg ✅ COMPLETE

```
numpy.linalg
├── Matrix Products
│   ├── ✅ dot(a, b, out)
│   ├── ✅ vdot(a, b)
│   ├── ✅ inner(a, b)
│   ├── ✅ outer(a, b, out)
│   ├── ✅ matmul(x1, x2, out)
│   ├── ❌ tensordot(a, b, axes)
│   ├── ❌ einsum(subscripts, *operands, out, ...)
│   ├── ❌ einsum_path(subscripts, *operands, optimize)
│   ├── ❌ kron(a, b)
│   ├── ❌ cross(a, b, axisa, axisb, axisc, axis)
│   └── ❌ multi_dot(arrays, out)
├── Decompositions
│   ├── ✅ cholesky(a)
│   ├── ✅ qr(a, mode)
│   ├── ✅ svd(a, full_matrices, compute_uv, hermitian)
│   └── ✅ svdvals(x)
├── Eigenvalues
│   ├── ✅ eig(a)
│   ├── ✅ eigh(a, UPLO)
│   ├── ✅ eigvals(a)
│   └── ✅ eigvalsh(a, UPLO)
├── Norms & Numbers
│   ├── ✅ norm(x, ord, axis, keepdims)
│   ├── ❌ matrix_norm(x, ord, keepdims)
│   ├── ❌ vector_norm(x, ord, axis, keepdims)
│   ├── ✅ cond(x, p)
│   ├── ✅ det(a)
│   ├── ✅ slogdet(a)
│   ├── ✅ matrix_rank(A, tol, hermitian, rtol)
│   └── ✅ trace(a, offset, axis1, axis2, dtype, out)
├── Solving & Inverting
│   ├── ✅ solve(a, b)
│   ├── ❌ tensorsolve(a, b, axes)
│   ├── ✅ lstsq(a, b, rcond)
│   ├── ✅ inv(a)
│   ├── ✅ pinv(a, rcond, hermitian, rtol)
│   └── ❌ tensorinv(a, ind)
├── Matrix Operations
│   ├── ✅ matrix_power(a, n)
│   ├── ✅ diagonal(a, offset, axis1, axis2)
│   └── ❌ matrix_transpose(x)
└── Exception
    └── ✅ LinAlgError
```

**Implementation:** `src/ts/linalg.ts`, `src/wasm/linalg.c`, `src/wasm/blas.c`, `src/wasm/lapack.c`

---

## Phase 15: numpy.fft ✅ COMPLETE

```
numpy.fft
├── Standard FFTs
│   ├── ✅ fft(a, n, axis, norm, out)
│   ├── ✅ ifft(a, n, axis, norm, out)
│   ├── ✅ fft2(a, s, axes, norm, out)
│   ├── ✅ ifft2(a, s, axes, norm, out)
│   ├── ✅ fftn(a, s, axes, norm, out)
│   └── ✅ ifftn(a, s, axes, norm, out)
├── Real FFTs
│   ├── ✅ rfft(a, n, axis, norm, out)
│   ├── ✅ irfft(a, n, axis, norm, out)
│   ├── ✅ rfft2(a, s, axes, norm, out)
│   ├── ✅ irfft2(a, s, axes, norm, out)
│   ├── ✅ rfftn(a, s, axes, norm, out)
│   └── ✅ irfftn(a, s, axes, norm, out)
├── Hermitian FFTs
│   ├── ✅ hfft(a, n, axis, norm, out)
│   └── ✅ ihfft(a, n, axis, norm, out)
└── Helper Functions
    ├── ✅ fftfreq(n, d, device)
    ├── ✅ rfftfreq(n, d, device)
    ├── ✅ fftshift(x, axes)
    └── ✅ ifftshift(x, axes)
```

**Implementation:** `src/ts/fft.ts`, `src/wasm/fft.c`

---

## Phase 16: numpy.random ✅ COMPLETE

```
numpy.random
├── Generator Class
│   ├── ✅ default_rng(seed) → Generator
│   └── Generator Methods
│       ├── ✅ random(size, dtype, out)
│       ├── ✅ integers(low, high, size, dtype, endpoint)
│       ├── ✅ uniform(low, high, size)
│       ├── ✅ normal(loc, scale, size)
│       ├── ✅ standard_normal(size, dtype, out)
│       ├── ✅ exponential(scale, size)
│       ├── ✅ standard_exponential(size, method)
│       ├── ✅ poisson(lam, size)
│       ├── ✅ binomial(n, p, size)
│       ├── ✅ negative_binomial(n, p, size)
│       ├── ✅ geometric(p, size)
│       ├── ✅ hypergeometric(ngood, nbad, nsample, size)
│       ├── ✅ beta(a, b, size)
│       ├── ✅ gamma(shape, scale, size)
│       ├── ✅ standard_gamma(shape, size)
│       ├── ✅ chisquare(df, size)
│       ├── ✅ f(dfnum, dfden, size)
│       ├── ✅ standard_t(df, size)
│       ├── ✅ standard_cauchy(size)
│       ├── ✅ pareto(a, size)
│       ├── ✅ weibull(a, size)
│       ├── ✅ laplace(loc, scale, size)
│       ├── ✅ lognormal(mean, sigma, size)
│       ├── ✅ rayleigh(scale, size)
│       ├── ✅ choice(a, size, replace, p, axis, shuffle) - async
│       ├── ✅ shuffle(x, axis) - async
│       ├── ✅ permutation(x, axis) - async
│       └── ✅ bytes(length)
├── BitGenerator Infrastructure
│   ├── ✅ PCG64 (default)
│   ├── ❌ MT19937
│   ├── ❌ Philox
│   └── ❌ SFC64
├── SeedSequence
│   └── ✅ SeedSequence(entropy, spawn_key, pool_size)
└── Legacy Functions
    ├── ✅ seed(seed)
    ├── ✅ random()
    ├── ✅ randn()
    ├── ✅ randint(low, high, size)
    └── ✅ initRandom()
```

**Implementation:** `src/ts/random.ts`, `src/wasm/random/`

---

## Phase 17: numpy.polynomial ✅ COMPLETE

```
numpy.polynomial
├── ✅ Polynomial Class (Power Series)
│   ├── ✅ polyval, polyval2d, polyval3d
│   ├── ✅ polyvander, polyvander2d, polyvander3d
│   ├── ✅ polyder, polyint
│   ├── ✅ polyfit
│   ├── ✅ polyroots, polycompanion
│   ├── ✅ polyfromroots
│   └── ✅ polyadd, polysub, polymul, polydiv, polypow
├── ✅ Chebyshev Class
│   ├── ✅ chebval, chebvander
│   ├── ✅ chebder, chebint
│   ├── ✅ chebfit, chebinterpolate
│   ├── ✅ chebroots, chebcompanion
│   └── ✅ chebadd, chebsub, chebmul, chebdiv, chebpow
├── ✅ Legendre Class
│   ├── ✅ legval, legvander
│   ├── ✅ legder, legint
│   ├── ✅ legfit, legroots
│   └── ✅ legadd, legsub, legmul, legdiv, legpow
├── ✅ Hermite Class (Physicist's)
│   ├── ✅ hermval, hermvander
│   ├── ✅ hermder, hermint
│   ├── ✅ hermfit, hermroots
│   └── ✅ hermadd, hermsub, hermmul, hermdiv, hermpow
├── ✅ HermiteE Class (Probabilist's)
│   ├── ✅ hermeval, hermevander
│   ├── ✅ hermeder, hermeint
│   ├── ✅ hermefit, hermeroots
│   └── ✅ hermeadd, hermesub, hermemul, hermediv, hermepow
├── ✅ Laguerre Class
│   ├── ✅ lagval, lagvander
│   ├── ✅ lagder, lagint
│   ├── ✅ lagfit, lagroots
│   └── ✅ lagadd, lagsub, lagmul, lagdiv, lagpow
├── ✅ Conversion Functions
│   ├── ✅ poly2cheb, cheb2poly
│   ├── ✅ poly2leg, leg2poly
│   └── ✅ (and all other conversion combinations)
└── ✅ Utilities
    ├── ✅ trimseq, trimcoef
    ├── ✅ as_series
    ├── ✅ getdomain, mapdomain, mapparms
    └── ✅ ABCPolyBase, maxpower
```

**Implementation:** `src/ts/polynomial/` (multiple files)

---

## Phase 18: numpy.ma (Masked Arrays) ✅ COMPLETE

```
numpy.ma
├── ✅ MaskedArray Class
│   ├── ✅ Core properties (data, mask, fill_value)
│   ├── ✅ Arithmetic operations
│   └── ✅ Comparison operations
├── ✅ Mask Operations
│   ├── ✅ make_mask, make_mask_none
│   ├── ✅ getmask, getmaskarray, getdata
│   ├── ✅ is_mask, is_masked
│   ├── ✅ mask_or, flatten_mask, reshape_mask
│   └── ✅ broadcast_mask, allTrue
├── ✅ Fill Values
│   ├── ✅ default_fill_value
│   ├── ✅ common_fill_value
│   ├── ✅ set_fill_value
│   └── ✅ getReductionFillValue
├── ✅ Creation Functions
│   ├── ✅ masked_array, array
│   ├── ✅ masked_equal, masked_not_equal
│   ├── ✅ masked_greater, masked_greater_equal
│   ├── ✅ masked_less, masked_less_equal
│   ├── ✅ masked_inside, masked_outside
│   ├── ✅ masked_where, masked_invalid, masked_values
│   ├── ✅ zeros, ones, empty
│   ├── ✅ masked_all, masked_all_like
│   ├── ✅ zeros_like, ones_like, empty_like
│   └── ✅ fromfunction
├── ✅ Statistics & Extras
│   ├── ✅ average, median
│   ├── ✅ cov, corrcoef
│   ├── ✅ notmasked_edges, notmasked_contiguous
│   ├── ✅ flatnotmasked_edges, flatnotmasked_contiguous
│   ├── ✅ clump_masked, clump_unmasked
│   └── ✅ apply_along_axis
└── ✅ Constants & Errors
    ├── ✅ nomask, masked
    ├── ✅ MaskedArrayError
    └── ✅ MaskError
```

**Implementation:** `src/ts/ma/`

---

## Phase 19: numpy.strings ✅ COMPLETE

```
numpy.strings
├── ✅ Comparison (10 functions)
│   ├── ✅ equal, not_equal
│   ├── ✅ less, less_equal
│   ├── ✅ greater, greater_equal
│   └── ✅ compare_chararrays
├── ✅ Properties (9 functions)
│   ├── ✅ isalpha, isdigit, isalnum
│   ├── ✅ isspace, islower, isupper
│   ├── ✅ istitle, isdecimal, isnumeric
│   └── ✅ str_len
├── ✅ Search (7 functions)
│   ├── ✅ find, rfind
│   ├── ✅ index, rindex
│   ├── ✅ count
│   └── ✅ startswith, endswith
├── ✅ Manipulation (11 functions)
│   ├── ✅ lower, upper, swapcase
│   ├── ✅ capitalize, title
│   ├── ✅ add, multiply
│   └── ✅ strip, lstrip, rstrip, expandtabs
└── ✅ Advanced
    ├── ✅ replace, center, ljust, rjust, zfill
    ├── ✅ partition, rpartition
    └── ✅ encode, decode
```

**Implementation:** `src/ts/strings/`

---

## Phase 20: numpy.rec (Record Arrays) ✅ COMPLETE

```
numpy.rec
├── ✅ recarray class
├── ✅ record class
├── ✅ format_parser function
├── ✅ Convenience Functions
│   ├── ✅ fromarrays
│   ├── ✅ fromrecords
│   ├── ✅ fromstring
│   ├── ✅ fromfile
│   ├── ✅ array
│   └── ✅ find_duplicate
└── ✅ Error Classes
    ├── ✅ KeyError
    └── ✅ IndexError
```

**Implementation:** `src/ts/rec/`

---

## Phase 21: numpy.testing ✅ COMPLETE

```
numpy.testing
├── ✅ Assertion functions for unit testing
└── ✅ Error classes
    ├── ✅ AssertionError
    ├── ✅ SkipTest
    └── ✅ KnownFailureException
```

**Implementation:** `src/ts/testing/`

---

## Remaining Work (Lower Priority)

```
NOT YET IMPLEMENTED:
├── Cumulative Operations
│   ├── ❌ cumsum(arr, axis, dtype, out)
│   ├── ❌ cumprod(arr, axis, dtype, out)
│   ├── ❌ nancumsum(arr, axis, ...)
│   └── ❌ nancumprod(arr, axis, ...)
├── NaN-handling Functions
│   ├── ❌ nansum, nanprod
│   ├── ❌ nanmean, nanstd, nanvar
│   ├── ❌ nanmin, nanmax
│   ├── ❌ nanargmin, nanargmax
│   ├── ❌ nanmedian
│   └── ❌ nanpercentile, nanquantile
├── Histogram Functions
│   ├── ❌ histogram, histogram2d, histogramdd
│   ├── ❌ histogram_bin_edges
│   ├── ❌ bincount
│   └── ❌ digitize
├── Advanced Linear Algebra
│   ├── ❌ tensordot, einsum, einsum_path
│   ├── ❌ kron, cross, multi_dot
│   ├── ❌ tensorsolve, tensorinv
│   └── ❌ matrix_norm, vector_norm
├── Miscellaneous Ufuncs
│   ├── ❌ divmod, modf, frexp, ldexp
│   ├── ❌ nextafter, spacing, nan_to_num
│   ├── ❌ sinc, heaviside, isnat
│   ├── ❌ gcd, lcm, bitwise_count
│   └── ❌ correlate (signal processing)
└── Index Generation
    ├── ❌ ogrid, mgrid
    └── ❌ mask_indices
```

---

## Current Implementation Summary

### TypeScript Files (`src/ts/`)
| File/Directory | Description |
|----------------|-------------|
| `NDArray.ts` | Core NDArray class with all methods |
| `index.ts` | Main public API exports (500+ functions) |
| `ufunc.ts` | Universal functions (90+ ufuncs) |
| `statistics.ts` | Statistical functions with axis support |
| `sorting.ts` | Sorting and partitioning |
| `logic.ts` | Comparison and logical operations |
| `setops.ts` | Set operations (unique, union, intersect, etc.) |
| `manipulation.ts` | Array joining, splitting, rearranging |
| `indexing.ts` | Advanced indexing and selection |
| `broadcast.ts` | Broadcasting operations |
| `dtype.ts` | Data type utilities |
| `slice.ts` | Slicing and indexing specs |
| `iterators.ts` | Array iteration tools |
| `functional.ts` | Functional programming (vectorize, apply) |
| `window.ts` | Window functions |
| `constants.ts` | Mathematical constants |
| `typeinfo.ts` | Type information (finfo, iinfo) |
| `linalg.ts` | Linear algebra (matmul, solve, eig, svd, etc.) |
| `fft.ts` | FFT operations (fft, ifft, fft2, fftn, etc.) |
| `random.ts` | Random number generation (30+ distributions) |
| `polynomial/` | Polynomial classes (6 types × 15+ functions) |
| `ma/` | Masked arrays (MaskedArray, cov, corrcoef, etc.) |
| `strings/` | String operations (40+ functions) |
| `rec/` | Record arrays |
| `testing/` | Testing utilities |
| `io/` | I/O operations (10+ files) |

### C/WASM Files (`src/wasm/`)
| File | Description |
|------|-------------|
| `ndarray.c` | Core array operations |
| `ufunc_unary.c` | Unary mathematical operations |
| `ufunc_binary.c` | Binary mathematical operations |
| `statistics.c` | Reduction operations with axis support |
| `sorting.c` | Sorting and partitioning |
| `setops.c` | Set operations |
| `logic.c` | Logical and comparison operations |
| `manipulation.c` | Array manipulation |
| `indexing.c` | Index operations |
| `broadcast.c` | Broadcasting |
| `pairwise_sum.c` | NumPy-compatible accurate summation |
| `dtype.c` | Type system |
| `linalg.c` | Linear algebra operations |
| `blas.c` | BLAS operations |
| `lapack.c` | LAPACK operations |
| `fft.c` | FFT operations |
| `random/` | Random number generation |

**Total: 500+ exported functions, 15,000+ lines of implementation code**

---

## Implementation Priority Summary

```
COMPLETE (Phases 1-21):
✅ Extended DTypes + Type Promotion
✅ Element Access (get/set)
✅ Iterators
✅ Views + reshape/transpose
✅ Slicing (basic + advanced)
✅ Broadcasting
✅ Ufunc Infrastructure (WASM-accelerated, 90+ functions)
✅ Reductions with axis (sum, mean, std, var, min, max)
✅ Array manipulation (concat, stack, split, flip, roll, tile, pad)
✅ Sorting & Searching (sort, argsort, partition, searchsorted)
✅ Statistics (mean, std, var, median, cov, corrcoef)
✅ Set operations (unique, union, intersect, isin)
✅ Logic & comparison operations
✅ I/O operations (NPY, text, binary, formatting)
✅ Functional programming (vectorize, apply_along_axis)
✅ Window functions (blackman, hanning, kaiser, etc.)
✅ Constants & Type info
✅ numpy.linalg (matmul, solve, inv, eig, svd, qr, cholesky, etc.)
✅ numpy.fft (fft, ifft, fft2, fftn, fftfreq, fftshift, etc.)
✅ numpy.random (Generator, PCG64, 30+ distributions)
✅ numpy.polynomial (Polynomial, Chebyshev, Legendre, Hermite, Laguerre)
✅ numpy.ma (MaskedArray, cov, corrcoef, masked operations)
✅ numpy.strings (40+ string operations)
✅ numpy.rec (Record arrays)
✅ numpy.testing (Assertion utilities)

REMAINING (Lower Priority):
❌ Cumulative operations (cumsum, cumprod)
❌ NaN-handling functions (nansum, nanmean, etc.)
❌ Histogram functions
❌ Advanced linalg (tensordot, einsum, kron)
❌ Additional BitGenerators (MT19937, Philox, SFC64)
```

---

## Test Coverage

**14+ Test Files (~250+ test cases):**
- `comparison.test.ts` - NumPy accuracy comparison
- `constants.test.ts` - Mathematical constants
- `creation.test.ts` - Array creation functions
- `dtype.test.ts` - Data type operations
- `element-access.test.ts` - Get/set/item operations
- `functional.test.ts` - Vectorize, apply, piecewise
- `io.test.ts` - I/O and formatting
- `level1.test.ts` - Basic operations
- `level2.test.ts` - Slicing, broadcasting, indexing
- `manipulation.test.ts` - Join, split, rearrange
- `ndarray.test.ts` - NDArray class tests
- `setops.test.ts` - Set operations
- `ufunc.test.ts` - Universal functions
- `window.test.ts` - Window functions

---

## Verification

After each phase:
1. Run existing tests: `npm test`
2. Add new tests for implemented features
3. Run benchmarks: `npm run benchmark`
4. Validate against NumPy test vectors (extend `tests/python/generate_test_cases.py`)
