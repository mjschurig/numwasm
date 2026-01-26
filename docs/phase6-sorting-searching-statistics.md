# Phase 6: Sorting, Searching & Statistics Implementation Plan

Phase 6 extends NumJS-WASM with three major components: **Sorting Functions** (sort, argsort, partition), **Searching Functions** (argmax, argmin, searchsorted), and **Statistics Functions** (mean, std, var, median). This layer enables data analysis, array ordering, and statistical computations critical for numerical computing.

---

## Current State (Levels 1-2 Complete)

```
src/wasm/
├── ndarray.h          # NDArray struct, views, slicing, shape ops
├── ndarray.c          # Core array operations (~800 lines)
├── dtype.h/c          # DType utilities
├── pairwise_sum.h/c   # Accurate pairwise summation (O(lg n) error)
├── broadcast.h/c      # Shape broadcasting
└── indexing.h/c       # take, put, nonzero, where, compress, extract (~626 lines)

src/ts/
├── types.ts           # DType enum, WasmModule interface, flags
├── NDArray.ts         # Factory methods, element access, shape ops
├── iterators.ts       # FlatIterator, nditer, ndenumerate, ndindex
├── dtype.ts           # Type promotion utilities
├── slice.ts           # Slice class, ellipsis, newaxis
├── broadcast.ts       # broadcastTo, broadcastArrays, broadcastShapes
├── indexing.ts        # Partial argsort/argmax/argmin (1D only, TypeScript)
├── wasm-loader.ts     # WASM module loading
└── index.ts           # Public exports
```

**Existing Partial Implementations (in src/ts/indexing.ts):**
- `argsort()` - Lines 386-400: TypeScript only, 1D arrays only
- `argmax()` - Lines 410-445: TypeScript only, axis=0 for 2D only
- `argmin()` - Lines 455-490: Same limitations as argmax

**Existing Infrastructure:**
- `pairwise_sum.c` (~151 lines) - Accurate summation algorithm with O(lg n) rounding error
- Broadcasting system fully implemented
- View system with slicing fully implemented

---

## Phase 6 Implementation Tree

```
PHASE 6: SORTING, SEARCHING & STATISTICS
│
├── 6.1 Sorting Functions (C + TypeScript)
│   ├── 6.1.1 C: Introsort implementation (quicksort + heapsort fallback)
│   ├── 6.1.2 C: ndarray_sort() → in-place sort along axis
│   ├── 6.1.3 C: ndarray_sort_copy() → sorted copy
│   ├── 6.1.4 C: ndarray_argsort() → indices that would sort
│   ├── 6.1.5 C: ndarray_partition() → partial sort around kth
│   ├── 6.1.6 C: ndarray_argpartition() → indices for partition
│   ├── 6.1.7 TS: sort() wrapper with axis support
│   ├── 6.1.8 TS: argsort() wrapper (replace existing)
│   └── 6.1.9 TS: partition(), argpartition() wrappers
│
│   Dependencies: Level 1 views, Level 2 slicing
│
├── 6.2 Searching Functions (C + TypeScript)
│   ├── 6.2.1 C: ndarray_argmax() → index of maximum along axis
│   ├── 6.2.2 C: ndarray_argmin() → index of minimum along axis
│   ├── 6.2.3 C: ndarray_searchsorted() → binary search
│   ├── 6.2.4 C: ndarray_nanargmax() → argmax ignoring NaN
│   ├── 6.2.5 C: ndarray_nanargmin() → argmin ignoring NaN
│   ├── 6.2.6 TS: argmax() wrapper with keepdims (replace existing)
│   ├── 6.2.7 TS: argmin() wrapper with keepdims (replace existing)
│   └── 6.2.8 TS: searchsorted() wrapper
│
│   Dependencies: 6.1.* (sorted arrays for searchsorted)
│
├── 6.3 Basic Statistics (C + TypeScript)
│   ├── 6.3.1 C: ndarray_sum_axis() → sum along axis (uses pairwise_sum)
│   ├── 6.3.2 C: ndarray_mean_axis() → mean along axis
│   ├── 6.3.3 C: ndarray_var_axis() → variance with ddof
│   ├── 6.3.4 C: ndarray_std_axis() → standard deviation
│   ├── 6.3.5 C: ndarray_min_axis() → minimum along axis
│   ├── 6.3.6 C: ndarray_max_axis() → maximum along axis
│   ├── 6.3.7 TS: sum(), mean(), var(), std() wrappers
│   └── 6.3.8 TS: min(), max() wrappers
│
│   Dependencies: pairwise_sum.c, Level 2 broadcasting
│
├── 6.4 Advanced Statistics (C + TypeScript)
│   ├── 6.4.1 C: ndarray_median() → median (uses partition)
│   ├── 6.4.2 C: ndarray_percentile() → nth percentile
│   ├── 6.4.3 C: ndarray_quantile() → quantile
│   └── 6.4.4 TS: median(), percentile(), quantile() wrappers
│
│   Dependencies: 6.1.* (partition), 6.3.*
│
├── 6.5 NaN-aware Statistics (C + TypeScript)
│   ├── 6.5.1 C: ndarray_nansum() → sum ignoring NaN
│   ├── 6.5.2 C: ndarray_nanmean() → mean ignoring NaN
│   ├── 6.5.3 C: ndarray_nanvar() → variance ignoring NaN
│   ├── 6.5.4 C: ndarray_nanstd() → std ignoring NaN
│   ├── 6.5.5 C: ndarray_nanmin() → min ignoring NaN
│   ├── 6.5.6 C: ndarray_nanmax() → max ignoring NaN
│   └── 6.5.7 TS: nansum(), nanmean(), nanvar(), nanstd(), nanmin(), nanmax()
│
│   Dependencies: 6.3.*
│
└── 6.6 Reduction Utilities (C + TypeScript)
    ├── 6.6.1 C: Helper for axis reduction with keepdims
    ├── 6.6.2 C: Helper for multi-axis reduction
    ├── 6.6.3 TS: prod() → product of elements
    └── 6.6.4 TS: ptp() → peak-to-peak (max - min)

    Dependencies: 6.3.*
```

---

## Detailed Implementation Specifications

### 6.1 Sorting Functions (C)

**File:** `src/wasm/sorting.h` (new file)

```c
#ifndef NUMJS_SORTING_H
#define NUMJS_SORTING_H

#include "ndarray.h"

/* Sort Kind Constants */
#define SORT_QUICKSORT  0  /* Introsort (quicksort + heapsort fallback) */
#define SORT_MERGESORT  1  /* Stable merge sort */
#define SORT_HEAPSORT   2  /* Heapsort */
#define SORT_STABLE     3  /* Stable sort (alias for mergesort) */

/**
 * Sort array in-place along specified axis.
 *
 * @param arr   Array to sort (modified in place)
 * @param axis  Axis along which to sort (-1 for last axis)
 * @param kind  Sorting algorithm (SORT_* constant)
 * @return      0 on success, -1 on error
 */
int ndarray_sort(NDArray* arr, int32_t axis, int32_t kind);

/**
 * Return a sorted copy of array.
 *
 * @param arr   Source array
 * @param axis  Axis along which to sort (-1 for last, INT32_MIN to flatten)
 * @param kind  Sorting algorithm
 * @return      New sorted array or NULL on error
 */
NDArray* ndarray_sort_copy(NDArray* arr, int32_t axis, int32_t kind);

/**
 * Return indices that would sort array along axis.
 *
 * @param arr   Source array
 * @param axis  Axis along which to sort (-1 for last, INT32_MIN to flatten)
 * @param kind  Sorting algorithm
 * @return      Array of indices (Int32) or NULL on error
 */
NDArray* ndarray_argsort(NDArray* arr, int32_t axis, int32_t kind);

/**
 * Return a partitioned copy of array.
 * Elements before kth are smaller, elements after are larger.
 *
 * @param arr   Source array
 * @param kth   Index to partition around
 * @param axis  Axis along which to partition
 * @return      Partitioned array or NULL on error
 */
NDArray* ndarray_partition(NDArray* arr, int32_t kth, int32_t axis);

/**
 * Return indices that would partition array.
 */
NDArray* ndarray_argpartition(NDArray* arr, int32_t kth, int32_t axis);

#endif /* NUMJS_SORTING_H */
```

**File:** `src/wasm/sorting.c` (new file) - Core sorting implementation

```c
#include "sorting.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#define EXPORT EMSCRIPTEN_KEEPALIVE
#else
#define EXPORT
#endif

/* Swap two elements at byte offsets */
static void swap_elements(void* data, size_t off1, size_t off2, size_t elem_size) {
    uint8_t temp[16];
    uint8_t* ptr = (uint8_t*)data;
    memcpy(temp, ptr + off1, elem_size);
    memcpy(ptr + off1, ptr + off2, elem_size);
    memcpy(ptr + off2, temp, elem_size);
}

/* Compare elements by dtype. NaN sorts to end. */
static int compare_elements(const void* a, const void* b, DType dtype) {
    switch (dtype) {
        case DTYPE_FLOAT64: {
            double va = *(const double*)a, vb = *(const double*)b;
            if (isnan(va) && isnan(vb)) return 0;
            if (isnan(va)) return 1;
            if (isnan(vb)) return -1;
            return (va > vb) - (va < vb);
        }
        case DTYPE_FLOAT32: {
            float va = *(const float*)a, vb = *(const float*)b;
            if (isnan(va) && isnan(vb)) return 0;
            if (isnan(va)) return 1;
            if (isnan(vb)) return -1;
            return (va > vb) - (va < vb);
        }
        case DTYPE_INT32: {
            int32_t va = *(const int32_t*)a, vb = *(const int32_t*)b;
            return (va > vb) - (va < vb);
        }
        case DTYPE_INT64: {
            int64_t va = *(const int64_t*)a, vb = *(const int64_t*)b;
            return (va > vb) - (va < vb);
        }
        default: return 0;
    }
}

/* Heapsort for guaranteed O(n log n) */
static void heapsort_1d(void* data, size_t n, size_t elem_size, DType dtype) {
    if (n < 2) return;
    uint8_t* arr = (uint8_t*)data;

    /* Build max heap */
    for (size_t i = n / 2; i > 0; i--) {
        size_t root = i - 1, child = 2 * root + 1;
        while (child < n) {
            if (child + 1 < n &&
                compare_elements(arr + child * elem_size,
                               arr + (child + 1) * elem_size, dtype) < 0)
                child++;
            if (compare_elements(arr + root * elem_size,
                               arr + child * elem_size, dtype) >= 0)
                break;
            swap_elements(data, root * elem_size, child * elem_size, elem_size);
            root = child;
            child = 2 * root + 1;
        }
    }

    /* Extract from heap */
    for (size_t i = n - 1; i > 0; i--) {
        swap_elements(data, 0, i * elem_size, elem_size);
        size_t root = 0, child = 1;
        while (child < i) {
            if (child + 1 < i &&
                compare_elements(arr + child * elem_size,
                               arr + (child + 1) * elem_size, dtype) < 0)
                child++;
            if (compare_elements(arr + root * elem_size,
                               arr + child * elem_size, dtype) >= 0)
                break;
            swap_elements(data, root * elem_size, child * elem_size, elem_size);
            root = child;
            child = 2 * root + 1;
        }
    }
}

/* Quicksort partition with median-of-three pivot */
static size_t qs_partition(void* data, size_t lo, size_t hi,
                           size_t elem_size, DType dtype) {
    uint8_t* arr = (uint8_t*)data;
    size_t mid = lo + (hi - lo) / 2;

    /* Median-of-three */
    if (compare_elements(arr + mid * elem_size, arr + lo * elem_size, dtype) < 0)
        swap_elements(data, lo * elem_size, mid * elem_size, elem_size);
    if (compare_elements(arr + hi * elem_size, arr + lo * elem_size, dtype) < 0)
        swap_elements(data, lo * elem_size, hi * elem_size, elem_size);
    if (compare_elements(arr + mid * elem_size, arr + hi * elem_size, dtype) < 0)
        swap_elements(data, mid * elem_size, hi * elem_size, elem_size);

    size_t i = lo;
    for (size_t j = lo; j < hi; j++) {
        if (compare_elements(arr + j * elem_size, arr + hi * elem_size, dtype) < 0) {
            swap_elements(data, i * elem_size, j * elem_size, elem_size);
            i++;
        }
    }
    swap_elements(data, i * elem_size, hi * elem_size, elem_size);
    return i;
}

/* Introsort: quicksort with heapsort fallback */
static void introsort_1d(void* data, size_t n, size_t elem_size,
                         DType dtype, int depth_limit) {
    while (n > 16) {
        if (depth_limit == 0) {
            heapsort_1d(data, n, elem_size, dtype);
            return;
        }
        depth_limit--;
        size_t p = qs_partition(data, 0, n - 1, elem_size, dtype);
        if (p < n - 1 - p) {
            introsort_1d(data, p, elem_size, dtype, depth_limit);
            data = (uint8_t*)data + (p + 1) * elem_size;
            n = n - p - 1;
        } else {
            introsort_1d((uint8_t*)data + (p + 1) * elem_size,
                        n - p - 1, elem_size, dtype, depth_limit);
            n = p;
        }
    }

    /* Insertion sort for small arrays */
    uint8_t* arr = (uint8_t*)data;
    uint8_t temp[16];
    for (size_t i = 1; i < n; i++) {
        memcpy(temp, arr + i * elem_size, elem_size);
        size_t j = i;
        while (j > 0 && compare_elements(arr + (j-1) * elem_size, temp, dtype) > 0) {
            memcpy(arr + j * elem_size, arr + (j-1) * elem_size, elem_size);
            j--;
        }
        memcpy(arr + j * elem_size, temp, elem_size);
    }
}

static void sort_1d(void* data, size_t n, size_t elem_size, DType dtype, int32_t kind) {
    if (n < 2) return;
    if (kind == SORT_HEAPSORT) {
        heapsort_1d(data, n, elem_size, dtype);
    } else {
        int depth = 0;
        for (size_t m = n; m > 0; m >>= 1) depth++;
        introsort_1d(data, n, elem_size, dtype, depth * 2);
    }
}

EXPORT int ndarray_sort(NDArray* arr, int32_t axis, int32_t kind) {
    if (!arr || !arr->data || !(arr->flags & NDARRAY_WRITEABLE)) return -1;

    if (axis < 0) axis += arr->ndim;
    if (axis < 0 || axis >= arr->ndim) return -1;

    size_t elem_size = dtype_size(arr->dtype);
    int32_t axis_size = arr->shape[axis];
    if (axis_size <= 1) return 0;

    size_t outer_size = 1, inner_size = 1;
    for (int i = 0; i < axis; i++) outer_size *= arr->shape[i];
    for (int i = axis + 1; i < arr->ndim; i++) inner_size *= arr->shape[i];

    if (inner_size == 1 && arr->strides[axis] == (int32_t)elem_size) {
        /* Contiguous along axis */
        for (size_t outer = 0; outer < outer_size; outer++) {
            size_t offset = outer * axis_size * elem_size;
            sort_1d((uint8_t*)arr->data + offset, axis_size, elem_size, arr->dtype, kind);
        }
    } else {
        /* Non-contiguous: use temp buffer */
        void* temp = malloc(axis_size * elem_size);
        if (!temp) return -1;

        for (size_t outer = 0; outer < outer_size; outer++) {
            for (size_t inner = 0; inner < inner_size; inner++) {
                for (int32_t i = 0; i < axis_size; i++) {
                    size_t idx = outer * axis_size * inner_size + i * inner_size + inner;
                    memcpy((uint8_t*)temp + i * elem_size,
                           (uint8_t*)arr->data + idx * elem_size, elem_size);
                }
                sort_1d(temp, axis_size, elem_size, arr->dtype, kind);
                for (int32_t i = 0; i < axis_size; i++) {
                    size_t idx = outer * axis_size * inner_size + i * inner_size + inner;
                    memcpy((uint8_t*)arr->data + idx * elem_size,
                           (uint8_t*)temp + i * elem_size, elem_size);
                }
            }
        }
        free(temp);
    }
    return 0;
}

EXPORT NDArray* ndarray_sort_copy(NDArray* arr, int32_t axis, int32_t kind) {
    if (!arr) return NULL;

    NDArray* result;
    if (axis == -2147483648) {  /* INT32_MIN: flatten */
        result = ndarray_flatten(arr);
        if (!result) return NULL;
        axis = 0;
    } else {
        result = ndarray_copy(arr);
        if (!result) return NULL;
    }

    if (ndarray_sort(result, axis, kind) != 0) {
        ndarray_free(result);
        return NULL;
    }
    return result;
}

EXPORT NDArray* ndarray_argsort(NDArray* arr, int32_t axis, int32_t kind) {
    if (!arr) return NULL;

    NDArray* work = arr;
    bool free_work = false;

    if (axis == -2147483648) {
        work = ndarray_flatten(arr);
        if (!work) return NULL;
        free_work = true;
        axis = 0;
    }

    if (axis < 0) axis += work->ndim;
    if (axis < 0 || axis >= work->ndim) {
        if (free_work) ndarray_free(work);
        return NULL;
    }

    NDArray* result = ndarray_empty(work->ndim, work->shape, DTYPE_INT32);
    if (!result) {
        if (free_work) ndarray_free(work);
        return NULL;
    }

    int32_t axis_size = work->shape[axis];
    size_t outer_size = 1, inner_size = 1;
    for (int i = 0; i < axis; i++) outer_size *= work->shape[i];
    for (int i = axis + 1; i < work->ndim; i++) inner_size *= work->shape[i];

    typedef struct { double val; int32_t idx; } IdxVal;
    IdxVal* indexed = malloc(axis_size * sizeof(IdxVal));
    if (!indexed) {
        ndarray_free(result);
        if (free_work) ndarray_free(work);
        return NULL;
    }

    for (size_t outer = 0; outer < outer_size; outer++) {
        for (size_t inner = 0; inner < inner_size; inner++) {
            for (int32_t i = 0; i < axis_size; i++) {
                size_t idx = outer * axis_size * inner_size + i * inner_size + inner;
                indexed[i].val = ndarray_get_flat(work, idx);
                indexed[i].idx = i;
            }
            /* Insertion sort for stability */
            for (int32_t i = 1; i < axis_size; i++) {
                IdxVal key = indexed[i];
                int32_t j = i - 1;
                while (j >= 0 && indexed[j].val > key.val) {
                    indexed[j + 1] = indexed[j];
                    j--;
                }
                indexed[j + 1] = key;
            }
            int32_t* res_data = (int32_t*)result->data;
            for (int32_t i = 0; i < axis_size; i++) {
                size_t idx = outer * axis_size * inner_size + i * inner_size + inner;
                res_data[idx] = indexed[i].idx;
            }
        }
    }

    free(indexed);
    if (free_work) ndarray_free(work);
    return result;
}

EXPORT NDArray* ndarray_partition(NDArray* arr, int32_t kth, int32_t axis) {
    NDArray* result = ndarray_copy(arr);
    if (!result) return NULL;

    if (axis < 0) axis += result->ndim;
    if (axis < 0 || axis >= result->ndim) {
        ndarray_free(result);
        return NULL;
    }

    int32_t axis_size = result->shape[axis];
    if (kth < 0) kth += axis_size;
    if (kth < 0 || kth >= axis_size) {
        ndarray_free(result);
        return NULL;
    }

    /* For simplicity, sort and return (true partition would be O(n)) */
    ndarray_sort(result, axis, SORT_QUICKSORT);
    return result;
}

EXPORT NDArray* ndarray_argpartition(NDArray* arr, int32_t kth, int32_t axis) {
    /* Use argsort for now */
    return ndarray_argsort(arr, axis, SORT_QUICKSORT);
}
```

---

### 6.2 Searching Functions (C)

**File:** `src/wasm/searching.h` (new file)

```c
#ifndef NUMJS_SEARCHING_H
#define NUMJS_SEARCHING_H

#include "ndarray.h"

#define SEARCH_LEFT   0
#define SEARCH_RIGHT  1

/**
 * Return indices of maximum values along axis.
 */
NDArray* ndarray_argmax(NDArray* arr, int32_t axis, bool keepdims);

/**
 * Return indices of minimum values along axis.
 */
NDArray* ndarray_argmin(NDArray* arr, int32_t axis, bool keepdims);

/**
 * Find indices where elements should be inserted to maintain order.
 */
NDArray* ndarray_searchsorted(NDArray* sorted_arr, NDArray* values,
                               int32_t side, NDArray* sorter);

#endif /* NUMJS_SEARCHING_H */
```

**File:** `src/wasm/searching.c` (new file)

```c
#include "searching.h"
#include <stdlib.h>
#include <math.h>

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#define EXPORT EMSCRIPTEN_KEEPALIVE
#else
#define EXPORT
#endif

static void get_reduced_shape(const int32_t* shape, int32_t ndim, int32_t axis,
                              bool keepdims, int32_t* out_shape, int32_t* out_ndim) {
    if (keepdims) {
        *out_ndim = ndim;
        for (int i = 0; i < ndim; i++)
            out_shape[i] = (i == axis) ? 1 : shape[i];
    } else {
        *out_ndim = (ndim > 1) ? ndim - 1 : 1;
        int j = 0;
        for (int i = 0; i < ndim; i++)
            if (i != axis) out_shape[j++] = shape[i];
        if (*out_ndim == 1 && j == 0) out_shape[0] = 1;
    }
}

EXPORT NDArray* ndarray_argmax(NDArray* arr, int32_t axis, bool keepdims) {
    if (!arr || !arr->data) return NULL;

    /* Flatten case */
    if (axis == -2147483648) {
        size_t max_idx = 0;
        double max_val = ndarray_get_flat(arr, 0);
        for (size_t i = 1; i < arr->size; i++) {
            double val = ndarray_get_flat(arr, i);
            if (!isnan(val) && (isnan(max_val) || val > max_val)) {
                max_val = val;
                max_idx = i;
            }
        }
        int32_t shape[1] = {1};
        NDArray* result = ndarray_create(keepdims ? arr->ndim : 1,
                                         keepdims ? NULL : shape, DTYPE_INT32);
        if (keepdims) {
            for (int i = 0; i < arr->ndim; i++) result->shape[i] = 1;
        }
        ((int32_t*)result->data)[0] = (int32_t)max_idx;
        return result;
    }

    if (axis < 0) axis += arr->ndim;
    if (axis < 0 || axis >= arr->ndim) return NULL;

    int32_t out_shape[32], out_ndim;
    get_reduced_shape(arr->shape, arr->ndim, axis, keepdims, out_shape, &out_ndim);

    NDArray* result = ndarray_create(out_ndim, out_shape, DTYPE_INT32);
    if (!result) return NULL;

    int32_t axis_size = arr->shape[axis];
    size_t outer_size = 1, inner_size = 1;
    for (int i = 0; i < axis; i++) outer_size *= arr->shape[i];
    for (int i = axis + 1; i < arr->ndim; i++) inner_size *= arr->shape[i];

    int32_t* res_data = (int32_t*)result->data;
    size_t res_idx = 0;

    for (size_t outer = 0; outer < outer_size; outer++) {
        for (size_t inner = 0; inner < inner_size; inner++) {
            int32_t max_idx = 0;
            double max_val = ndarray_get_flat(arr,
                outer * axis_size * inner_size + inner);

            for (int32_t i = 1; i < axis_size; i++) {
                size_t idx = outer * axis_size * inner_size + i * inner_size + inner;
                double val = ndarray_get_flat(arr, idx);
                if (!isnan(val) && (isnan(max_val) || val > max_val)) {
                    max_val = val;
                    max_idx = i;
                }
            }
            res_data[res_idx++] = max_idx;
        }
    }
    return result;
}

EXPORT NDArray* ndarray_argmin(NDArray* arr, int32_t axis, bool keepdims) {
    if (!arr || !arr->data) return NULL;

    if (axis == -2147483648) {
        size_t min_idx = 0;
        double min_val = ndarray_get_flat(arr, 0);
        for (size_t i = 1; i < arr->size; i++) {
            double val = ndarray_get_flat(arr, i);
            if (!isnan(val) && (isnan(min_val) || val < min_val)) {
                min_val = val;
                min_idx = i;
            }
        }
        int32_t shape[1] = {1};
        NDArray* result = ndarray_create(keepdims ? arr->ndim : 1,
                                         keepdims ? NULL : shape, DTYPE_INT32);
        if (keepdims) {
            for (int i = 0; i < arr->ndim; i++) result->shape[i] = 1;
        }
        ((int32_t*)result->data)[0] = (int32_t)min_idx;
        return result;
    }

    if (axis < 0) axis += arr->ndim;
    if (axis < 0 || axis >= arr->ndim) return NULL;

    int32_t out_shape[32], out_ndim;
    get_reduced_shape(arr->shape, arr->ndim, axis, keepdims, out_shape, &out_ndim);

    NDArray* result = ndarray_create(out_ndim, out_shape, DTYPE_INT32);
    if (!result) return NULL;

    int32_t axis_size = arr->shape[axis];
    size_t outer_size = 1, inner_size = 1;
    for (int i = 0; i < axis; i++) outer_size *= arr->shape[i];
    for (int i = axis + 1; i < arr->ndim; i++) inner_size *= arr->shape[i];

    int32_t* res_data = (int32_t*)result->data;
    size_t res_idx = 0;

    for (size_t outer = 0; outer < outer_size; outer++) {
        for (size_t inner = 0; inner < inner_size; inner++) {
            int32_t min_idx = 0;
            double min_val = ndarray_get_flat(arr,
                outer * axis_size * inner_size + inner);

            for (int32_t i = 1; i < axis_size; i++) {
                size_t idx = outer * axis_size * inner_size + i * inner_size + inner;
                double val = ndarray_get_flat(arr, idx);
                if (!isnan(val) && (isnan(min_val) || val < min_val)) {
                    min_val = val;
                    min_idx = i;
                }
            }
            res_data[res_idx++] = min_idx;
        }
    }
    return result;
}

EXPORT NDArray* ndarray_searchsorted(NDArray* sorted_arr, NDArray* values,
                                      int32_t side, NDArray* sorter) {
    if (!sorted_arr || !values || sorted_arr->ndim != 1) return NULL;

    NDArray* result = ndarray_empty(values->ndim, values->shape, DTYPE_INT32);
    if (!result) return NULL;

    int32_t* res_data = (int32_t*)result->data;

    for (size_t i = 0; i < values->size; i++) {
        double val = ndarray_get_flat(values, i);
        size_t lo = 0, hi = sorted_arr->size;

        while (lo < hi) {
            size_t mid = lo + (hi - lo) / 2;
            size_t idx = sorter ? (size_t)ndarray_get_flat(sorter, mid) : mid;
            double arr_val = ndarray_get_flat(sorted_arr, idx);

            if (side == SEARCH_RIGHT ? arr_val <= val : arr_val < val)
                lo = mid + 1;
            else
                hi = mid;
        }
        res_data[i] = (int32_t)lo;
    }
    return result;
}
```

---

### 6.3 Statistics Functions (C)

**File:** `src/wasm/statistics.h` (new file)

```c
#ifndef NUMJS_STATISTICS_H
#define NUMJS_STATISTICS_H

#include "ndarray.h"

/* Basic reductions with axis support */
NDArray* ndarray_sum_axis(NDArray* arr, int32_t axis, bool keepdims, DType dtype);
NDArray* ndarray_mean_axis(NDArray* arr, int32_t axis, bool keepdims, DType dtype);
NDArray* ndarray_var_axis(NDArray* arr, int32_t axis, bool keepdims, DType dtype, int32_t ddof);
NDArray* ndarray_std_axis(NDArray* arr, int32_t axis, bool keepdims, DType dtype, int32_t ddof);
NDArray* ndarray_min_axis(NDArray* arr, int32_t axis, bool keepdims);
NDArray* ndarray_max_axis(NDArray* arr, int32_t axis, bool keepdims);

/* Advanced statistics */
NDArray* ndarray_median(NDArray* arr, int32_t axis, bool keepdims);
NDArray* ndarray_percentile(NDArray* arr, double q, int32_t axis, bool keepdims);

/* NaN-aware versions */
NDArray* ndarray_nansum(NDArray* arr, int32_t axis, bool keepdims, DType dtype);
NDArray* ndarray_nanmean(NDArray* arr, int32_t axis, bool keepdims, DType dtype);
NDArray* ndarray_nanvar(NDArray* arr, int32_t axis, bool keepdims, DType dtype, int32_t ddof);
NDArray* ndarray_nanstd(NDArray* arr, int32_t axis, bool keepdims, DType dtype, int32_t ddof);

#endif /* NUMJS_STATISTICS_H */
```

**File:** `src/wasm/statistics.c` (new file)

```c
#include "statistics.h"
#include "pairwise_sum.h"
#include <stdlib.h>
#include <math.h>

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#define EXPORT EMSCRIPTEN_KEEPALIVE
#else
#define EXPORT
#endif

static void get_reduced_shape(const int32_t* shape, int32_t ndim, int32_t axis,
                              bool keepdims, int32_t* out_shape, int32_t* out_ndim) {
    if (axis == -2147483648) {
        if (keepdims) {
            *out_ndim = ndim;
            for (int i = 0; i < ndim; i++) out_shape[i] = 1;
        } else {
            *out_ndim = 0;
        }
        return;
    }
    if (keepdims) {
        *out_ndim = ndim;
        for (int i = 0; i < ndim; i++)
            out_shape[i] = (i == axis) ? 1 : shape[i];
    } else {
        *out_ndim = (ndim > 1) ? ndim - 1 : 0;
        int j = 0;
        for (int i = 0; i < ndim; i++)
            if (i != axis) out_shape[j++] = shape[i];
    }
}

static DType get_reduction_dtype(DType input, DType requested) {
    if (requested != (DType)-1) return requested;
    switch (input) {
        case DTYPE_INT8: case DTYPE_INT16: case DTYPE_INT32: case DTYPE_INT64:
        case DTYPE_UINT8: case DTYPE_UINT16: case DTYPE_UINT32: case DTYPE_UINT64:
        case DTYPE_BOOL:
            return DTYPE_FLOAT64;
        default:
            return input;
    }
}

EXPORT NDArray* ndarray_sum_axis(NDArray* arr, int32_t axis, bool keepdims, DType dtype) {
    if (!arr || !arr->data) return NULL;
    DType out_dtype = get_reduction_dtype(arr->dtype, dtype);

    if (axis == -2147483648) {
        double total = ndarray_sum(arr);  /* Uses pairwise sum */
        if (keepdims) {
            int32_t* shape = malloc(arr->ndim * sizeof(int32_t));
            for (int i = 0; i < arr->ndim; i++) shape[i] = 1;
            NDArray* r = ndarray_full(arr->ndim, shape, out_dtype, total);
            free(shape);
            return r;
        }
        return ndarray_scalar(total, out_dtype);
    }

    if (axis < 0) axis += arr->ndim;
    if (axis < 0 || axis >= arr->ndim) return NULL;

    int32_t out_shape[32], out_ndim;
    get_reduced_shape(arr->shape, arr->ndim, axis, keepdims, out_shape, &out_ndim);

    NDArray* result = out_ndim == 0 ? ndarray_scalar(0, out_dtype) :
                      ndarray_create(out_ndim, out_shape, out_dtype);
    if (!result) return NULL;

    int32_t axis_size = arr->shape[axis];
    size_t outer = 1, inner = 1;
    for (int i = 0; i < axis; i++) outer *= arr->shape[i];
    for (int i = axis + 1; i < arr->ndim; i++) inner *= arr->shape[i];

    double* temp = malloc(axis_size * sizeof(double));
    size_t res_idx = 0;

    for (size_t o = 0; o < outer; o++) {
        for (size_t n = 0; n < inner; n++) {
            for (int32_t i = 0; i < axis_size; i++) {
                size_t idx = o * axis_size * inner + i * inner + n;
                temp[i] = ndarray_get_flat(arr, idx);
            }
            double sum = pairwise_sum_f64(temp, axis_size, 1);
            ndarray_set_flat(result, res_idx++, sum);
        }
    }
    free(temp);
    return result;
}

EXPORT NDArray* ndarray_mean_axis(NDArray* arr, int32_t axis, bool keepdims, DType dtype) {
    if (!arr) return NULL;

    if (axis == -2147483648) {
        double mean_val = ndarray_sum(arr) / (double)arr->size;
        DType out_dtype = get_reduction_dtype(arr->dtype, dtype);
        if (keepdims) {
            int32_t* shape = malloc(arr->ndim * sizeof(int32_t));
            for (int i = 0; i < arr->ndim; i++) shape[i] = 1;
            NDArray* r = ndarray_full(arr->ndim, shape, out_dtype, mean_val);
            free(shape);
            return r;
        }
        return ndarray_scalar(mean_val, out_dtype);
    }

    NDArray* sum_res = ndarray_sum_axis(arr, axis, keepdims, dtype);
    if (!sum_res) return NULL;

    int32_t norm_axis = axis < 0 ? axis + arr->ndim : axis;
    int32_t count = arr->shape[norm_axis];

    for (size_t i = 0; i < sum_res->size; i++) {
        double val = ndarray_get_flat(sum_res, i);
        ndarray_set_flat(sum_res, i, val / count);
    }
    return sum_res;
}

EXPORT NDArray* ndarray_var_axis(NDArray* arr, int32_t axis, bool keepdims,
                                  DType dtype, int32_t ddof) {
    if (!arr) return NULL;
    DType out_dtype = get_reduction_dtype(arr->dtype, dtype);

    if (axis == -2147483648) {
        double mean_val = ndarray_sum(arr) / (double)arr->size;
        double sum_sq = 0;
        for (size_t i = 0; i < arr->size; i++) {
            double diff = ndarray_get_flat(arr, i) - mean_val;
            sum_sq += diff * diff;
        }
        double var_val = sum_sq / (arr->size - ddof);
        if (keepdims) {
            int32_t* shape = malloc(arr->ndim * sizeof(int32_t));
            for (int i = 0; i < arr->ndim; i++) shape[i] = 1;
            NDArray* r = ndarray_full(arr->ndim, shape, out_dtype, var_val);
            free(shape);
            return r;
        }
        return ndarray_scalar(var_val, out_dtype);
    }

    if (axis < 0) axis += arr->ndim;
    if (axis < 0 || axis >= arr->ndim) return NULL;

    NDArray* mean_res = ndarray_mean_axis(arr, axis, true, dtype);
    if (!mean_res) return NULL;

    int32_t out_shape[32], out_ndim;
    get_reduced_shape(arr->shape, arr->ndim, axis, keepdims, out_shape, &out_ndim);

    NDArray* result = out_ndim == 0 ? ndarray_scalar(0, out_dtype) :
                      ndarray_create(out_ndim, out_shape, out_dtype);
    if (!result) { ndarray_free(mean_res); return NULL; }

    int32_t axis_size = arr->shape[axis];
    size_t outer = 1, inner = 1;
    for (int i = 0; i < axis; i++) outer *= arr->shape[i];
    for (int i = axis + 1; i < arr->ndim; i++) inner *= arr->shape[i];

    size_t res_idx = 0, mean_idx = 0;
    for (size_t o = 0; o < outer; o++) {
        for (size_t n = 0; n < inner; n++) {
            double mean_val = ndarray_get_flat(mean_res, mean_idx++);
            double sum_sq = 0;
            for (int32_t i = 0; i < axis_size; i++) {
                size_t idx = o * axis_size * inner + i * inner + n;
                double diff = ndarray_get_flat(arr, idx) - mean_val;
                sum_sq += diff * diff;
            }
            ndarray_set_flat(result, res_idx++, sum_sq / (axis_size - ddof));
        }
    }
    ndarray_free(mean_res);
    return result;
}

EXPORT NDArray* ndarray_std_axis(NDArray* arr, int32_t axis, bool keepdims,
                                  DType dtype, int32_t ddof) {
    NDArray* var_res = ndarray_var_axis(arr, axis, keepdims, dtype, ddof);
    if (!var_res) return NULL;
    for (size_t i = 0; i < var_res->size; i++) {
        ndarray_set_flat(var_res, i, sqrt(ndarray_get_flat(var_res, i)));
    }
    return var_res;
}

EXPORT NDArray* ndarray_min_axis(NDArray* arr, int32_t axis, bool keepdims) {
    if (!arr) return NULL;

    if (axis == -2147483648) {
        double min_val = ndarray_get_flat(arr, 0);
        for (size_t i = 1; i < arr->size; i++) {
            double v = ndarray_get_flat(arr, i);
            if (v < min_val) min_val = v;
        }
        if (keepdims) {
            int32_t* shape = malloc(arr->ndim * sizeof(int32_t));
            for (int i = 0; i < arr->ndim; i++) shape[i] = 1;
            NDArray* r = ndarray_full(arr->ndim, shape, arr->dtype, min_val);
            free(shape);
            return r;
        }
        return ndarray_scalar(min_val, arr->dtype);
    }

    if (axis < 0) axis += arr->ndim;
    if (axis < 0 || axis >= arr->ndim) return NULL;

    int32_t out_shape[32], out_ndim;
    get_reduced_shape(arr->shape, arr->ndim, axis, keepdims, out_shape, &out_ndim);

    NDArray* result = out_ndim == 0 ? ndarray_scalar(0, arr->dtype) :
                      ndarray_create(out_ndim, out_shape, arr->dtype);
    if (!result) return NULL;

    int32_t axis_size = arr->shape[axis];
    size_t outer = 1, inner = 1;
    for (int i = 0; i < axis; i++) outer *= arr->shape[i];
    for (int i = axis + 1; i < arr->ndim; i++) inner *= arr->shape[i];

    size_t res_idx = 0;
    for (size_t o = 0; o < outer; o++) {
        for (size_t n = 0; n < inner; n++) {
            double min_val = ndarray_get_flat(arr, o * axis_size * inner + n);
            for (int32_t i = 1; i < axis_size; i++) {
                double v = ndarray_get_flat(arr, o * axis_size * inner + i * inner + n);
                if (v < min_val) min_val = v;
            }
            ndarray_set_flat(result, res_idx++, min_val);
        }
    }
    return result;
}

EXPORT NDArray* ndarray_max_axis(NDArray* arr, int32_t axis, bool keepdims) {
    if (!arr) return NULL;

    if (axis == -2147483648) {
        double max_val = ndarray_get_flat(arr, 0);
        for (size_t i = 1; i < arr->size; i++) {
            double v = ndarray_get_flat(arr, i);
            if (v > max_val) max_val = v;
        }
        if (keepdims) {
            int32_t* shape = malloc(arr->ndim * sizeof(int32_t));
            for (int i = 0; i < arr->ndim; i++) shape[i] = 1;
            NDArray* r = ndarray_full(arr->ndim, shape, arr->dtype, max_val);
            free(shape);
            return r;
        }
        return ndarray_scalar(max_val, arr->dtype);
    }

    if (axis < 0) axis += arr->ndim;
    if (axis < 0 || axis >= arr->ndim) return NULL;

    int32_t out_shape[32], out_ndim;
    get_reduced_shape(arr->shape, arr->ndim, axis, keepdims, out_shape, &out_ndim);

    NDArray* result = out_ndim == 0 ? ndarray_scalar(0, arr->dtype) :
                      ndarray_create(out_ndim, out_shape, arr->dtype);
    if (!result) return NULL;

    int32_t axis_size = arr->shape[axis];
    size_t outer = 1, inner = 1;
    for (int i = 0; i < axis; i++) outer *= arr->shape[i];
    for (int i = axis + 1; i < arr->ndim; i++) inner *= arr->shape[i];

    size_t res_idx = 0;
    for (size_t o = 0; o < outer; o++) {
        for (size_t n = 0; n < inner; n++) {
            double max_val = ndarray_get_flat(arr, o * axis_size * inner + n);
            for (int32_t i = 1; i < axis_size; i++) {
                double v = ndarray_get_flat(arr, o * axis_size * inner + i * inner + n);
                if (v > max_val) max_val = v;
            }
            ndarray_set_flat(result, res_idx++, max_val);
        }
    }
    return result;
}

EXPORT NDArray* ndarray_median(NDArray* arr, int32_t axis, bool keepdims) {
    if (!arr) return NULL;

    if (axis == -2147483648) {
        NDArray* flat = ndarray_flatten(arr);
        ndarray_sort(flat, 0, 0);
        size_t n = flat->size, mid = n / 2;
        double median_val = (n % 2 == 0) ?
            (ndarray_get_flat(flat, mid-1) + ndarray_get_flat(flat, mid)) / 2.0 :
            ndarray_get_flat(flat, mid);
        ndarray_free(flat);

        if (keepdims) {
            int32_t* shape = malloc(arr->ndim * sizeof(int32_t));
            for (int i = 0; i < arr->ndim; i++) shape[i] = 1;
            NDArray* r = ndarray_full(arr->ndim, shape, DTYPE_FLOAT64, median_val);
            free(shape);
            return r;
        }
        return ndarray_scalar(median_val, DTYPE_FLOAT64);
    }

    /* Axis-specific median implementation */
    if (axis < 0) axis += arr->ndim;
    if (axis < 0 || axis >= arr->ndim) return NULL;

    int32_t out_shape[32], out_ndim;
    get_reduced_shape(arr->shape, arr->ndim, axis, keepdims, out_shape, &out_ndim);

    NDArray* result = out_ndim == 0 ? ndarray_scalar(0, DTYPE_FLOAT64) :
                      ndarray_create(out_ndim, out_shape, DTYPE_FLOAT64);
    if (!result) return NULL;

    int32_t axis_size = arr->shape[axis];
    size_t mid = axis_size / 2;
    size_t outer = 1, inner = 1;
    for (int i = 0; i < axis; i++) outer *= arr->shape[i];
    for (int i = axis + 1; i < arr->ndim; i++) inner *= arr->shape[i];

    double* temp = malloc(axis_size * sizeof(double));
    size_t res_idx = 0;

    for (size_t o = 0; o < outer; o++) {
        for (size_t n = 0; n < inner; n++) {
            for (int32_t i = 0; i < axis_size; i++) {
                temp[i] = ndarray_get_flat(arr, o * axis_size * inner + i * inner + n);
            }
            /* Sort temp */
            for (int32_t i = 1; i < axis_size; i++) {
                double key = temp[i];
                int32_t j = i - 1;
                while (j >= 0 && temp[j] > key) { temp[j+1] = temp[j]; j--; }
                temp[j+1] = key;
            }
            double median_val = (axis_size % 2 == 0) ?
                (temp[mid-1] + temp[mid]) / 2.0 : temp[mid];
            ndarray_set_flat(result, res_idx++, median_val);
        }
    }
    free(temp);
    return result;
}

EXPORT NDArray* ndarray_percentile(NDArray* arr, double q, int32_t axis, bool keepdims) {
    return ndarray_quantile(arr, q / 100.0, axis, keepdims);
}

EXPORT NDArray* ndarray_quantile(NDArray* arr, double q, int32_t axis, bool keepdims) {
    if (!arr || q < 0 || q > 1) return NULL;

    if (axis == -2147483648) {
        NDArray* flat = ndarray_flatten(arr);
        ndarray_sort(flat, 0, 0);
        size_t n = flat->size;
        double idx = q * (n - 1);
        size_t lo = (size_t)floor(idx), hi = (size_t)ceil(idx);
        double frac = idx - lo;
        double qval = (lo == hi) ? ndarray_get_flat(flat, lo) :
            ndarray_get_flat(flat, lo) * (1 - frac) + ndarray_get_flat(flat, hi) * frac;
        ndarray_free(flat);

        if (keepdims) {
            int32_t* shape = malloc(arr->ndim * sizeof(int32_t));
            for (int i = 0; i < arr->ndim; i++) shape[i] = 1;
            NDArray* r = ndarray_full(arr->ndim, shape, DTYPE_FLOAT64, qval);
            free(shape);
            return r;
        }
        return ndarray_scalar(qval, DTYPE_FLOAT64);
    }
    return NULL;  /* Axis-specific implementation TODO */
}

/* NaN-aware implementations */
EXPORT NDArray* ndarray_nansum(NDArray* arr, int32_t axis, bool keepdims, DType dtype) {
    if (!arr) return NULL;
    DType out_dtype = get_reduction_dtype(arr->dtype, dtype);

    if (axis == -2147483648) {
        double total = 0;
        for (size_t i = 0; i < arr->size; i++) {
            double v = ndarray_get_flat(arr, i);
            if (!isnan(v)) total += v;
        }
        if (keepdims) {
            int32_t* shape = malloc(arr->ndim * sizeof(int32_t));
            for (int i = 0; i < arr->ndim; i++) shape[i] = 1;
            NDArray* r = ndarray_full(arr->ndim, shape, out_dtype, total);
            free(shape);
            return r;
        }
        return ndarray_scalar(total, out_dtype);
    }
    return NULL;  /* Axis-specific TODO */
}

EXPORT NDArray* ndarray_nanmean(NDArray* arr, int32_t axis, bool keepdims, DType dtype) {
    if (!arr) return NULL;
    DType out_dtype = get_reduction_dtype(arr->dtype, dtype);

    if (axis == -2147483648) {
        double total = 0;
        size_t count = 0;
        for (size_t i = 0; i < arr->size; i++) {
            double v = ndarray_get_flat(arr, i);
            if (!isnan(v)) { total += v; count++; }
        }
        double mean_val = count > 0 ? total / count : NAN;
        if (keepdims) {
            int32_t* shape = malloc(arr->ndim * sizeof(int32_t));
            for (int i = 0; i < arr->ndim; i++) shape[i] = 1;
            NDArray* r = ndarray_full(arr->ndim, shape, out_dtype, mean_val);
            free(shape);
            return r;
        }
        return ndarray_scalar(mean_val, out_dtype);
    }
    return NULL;
}

EXPORT NDArray* ndarray_nanvar(NDArray* arr, int32_t axis, bool keepdims,
                                DType dtype, int32_t ddof) {
    return NULL;  /* TODO */
}

EXPORT NDArray* ndarray_nanstd(NDArray* arr, int32_t axis, bool keepdims,
                                DType dtype, int32_t ddof) {
    return NULL;  /* TODO */
}
```

---

## TypeScript Wrappers

**File:** `src/ts/sorting.ts` (new file)

```typescript
import { NDArray } from './NDArray.js';
import { DType } from './types.js';

export type SortKind = 'quicksort' | 'mergesort' | 'heapsort' | 'stable';
const KIND_MAP: Record<SortKind, number> = {
  quicksort: 0, mergesort: 1, heapsort: 2, stable: 3
};

export async function sort(a: NDArray, axis: number | null = -1,
                           kind: SortKind = 'quicksort'): Promise<NDArray> {
  const axisVal = axis === null ? -2147483648 : axis;
  const resultPtr = a._wasmModule._ndarray_sort_copy(a._wasmPtr, axisVal, KIND_MAP[kind]);
  if (resultPtr === 0) throw new Error('sort failed');
  return NDArray._fromPtr(resultPtr, a._wasmModule);
}

export async function argsort(a: NDArray, axis: number | null = -1,
                              kind: SortKind = 'quicksort'): Promise<NDArray> {
  const axisVal = axis === null ? -2147483648 : axis;
  const resultPtr = a._wasmModule._ndarray_argsort(a._wasmPtr, axisVal, KIND_MAP[kind]);
  if (resultPtr === 0) throw new Error('argsort failed');
  return NDArray._fromPtr(resultPtr, a._wasmModule);
}

export async function partition(a: NDArray, kth: number,
                                axis: number = -1): Promise<NDArray> {
  const resultPtr = a._wasmModule._ndarray_partition(a._wasmPtr, kth, axis);
  if (resultPtr === 0) throw new Error('partition failed');
  return NDArray._fromPtr(resultPtr, a._wasmModule);
}

export async function argpartition(a: NDArray, kth: number,
                                   axis: number = -1): Promise<NDArray> {
  const resultPtr = a._wasmModule._ndarray_argpartition(a._wasmPtr, kth, axis);
  if (resultPtr === 0) throw new Error('argpartition failed');
  return NDArray._fromPtr(resultPtr, a._wasmModule);
}
```

**File:** `src/ts/statistics.ts` (new file)

```typescript
import { NDArray } from './NDArray.js';
import { DType } from './types.js';

export async function sum(a: NDArray, axis: number | null = null,
                          keepdims = false, dtype?: DType): Promise<NDArray | number> {
  const axisVal = axis === null ? -2147483648 : axis;
  const resultPtr = a._wasmModule._ndarray_sum_axis(a._wasmPtr, axisVal, keepdims, dtype ?? -1);
  if (resultPtr === 0) throw new Error('sum failed');
  const result = NDArray._fromPtr(resultPtr, a._wasmModule);
  if (axis === null && !keepdims) { const v = result.item(); result.dispose(); return v; }
  return result;
}

export async function mean(a: NDArray, axis: number | null = null,
                           keepdims = false, dtype?: DType): Promise<NDArray | number> {
  const axisVal = axis === null ? -2147483648 : axis;
  const resultPtr = a._wasmModule._ndarray_mean_axis(a._wasmPtr, axisVal, keepdims, dtype ?? -1);
  if (resultPtr === 0) throw new Error('mean failed');
  const result = NDArray._fromPtr(resultPtr, a._wasmModule);
  if (axis === null && !keepdims) { const v = result.item(); result.dispose(); return v; }
  return result;
}

export async function variance(a: NDArray, axis: number | null = null, ddof = 0,
                               keepdims = false, dtype?: DType): Promise<NDArray | number> {
  const axisVal = axis === null ? -2147483648 : axis;
  const resultPtr = a._wasmModule._ndarray_var_axis(a._wasmPtr, axisVal, keepdims, dtype ?? -1, ddof);
  if (resultPtr === 0) throw new Error('variance failed');
  const result = NDArray._fromPtr(resultPtr, a._wasmModule);
  if (axis === null && !keepdims) { const v = result.item(); result.dispose(); return v; }
  return result;
}
export const var_ = variance;

export async function std(a: NDArray, axis: number | null = null, ddof = 0,
                          keepdims = false, dtype?: DType): Promise<NDArray | number> {
  const axisVal = axis === null ? -2147483648 : axis;
  const resultPtr = a._wasmModule._ndarray_std_axis(a._wasmPtr, axisVal, keepdims, dtype ?? -1, ddof);
  if (resultPtr === 0) throw new Error('std failed');
  const result = NDArray._fromPtr(resultPtr, a._wasmModule);
  if (axis === null && !keepdims) { const v = result.item(); result.dispose(); return v; }
  return result;
}

export async function min(a: NDArray, axis: number | null = null,
                          keepdims = false): Promise<NDArray | number> {
  const axisVal = axis === null ? -2147483648 : axis;
  const resultPtr = a._wasmModule._ndarray_min_axis(a._wasmPtr, axisVal, keepdims);
  if (resultPtr === 0) throw new Error('min failed');
  const result = NDArray._fromPtr(resultPtr, a._wasmModule);
  if (axis === null && !keepdims) { const v = result.item(); result.dispose(); return v; }
  return result;
}

export async function max(a: NDArray, axis: number | null = null,
                          keepdims = false): Promise<NDArray | number> {
  const axisVal = axis === null ? -2147483648 : axis;
  const resultPtr = a._wasmModule._ndarray_max_axis(a._wasmPtr, axisVal, keepdims);
  if (resultPtr === 0) throw new Error('max failed');
  const result = NDArray._fromPtr(resultPtr, a._wasmModule);
  if (axis === null && !keepdims) { const v = result.item(); result.dispose(); return v; }
  return result;
}

export async function median(a: NDArray, axis: number | null = null,
                             keepdims = false): Promise<NDArray | number> {
  const axisVal = axis === null ? -2147483648 : axis;
  const resultPtr = a._wasmModule._ndarray_median(a._wasmPtr, axisVal, keepdims);
  if (resultPtr === 0) throw new Error('median failed');
  const result = NDArray._fromPtr(resultPtr, a._wasmModule);
  if (axis === null && !keepdims) { const v = result.item(); result.dispose(); return v; }
  return result;
}

export async function argmax(a: NDArray, axis: number | null = null,
                             keepdims = false): Promise<NDArray | number> {
  const axisVal = axis === null ? -2147483648 : axis;
  const resultPtr = a._wasmModule._ndarray_argmax(a._wasmPtr, axisVal, keepdims);
  if (resultPtr === 0) throw new Error('argmax failed');
  const result = NDArray._fromPtr(resultPtr, a._wasmModule);
  if (axis === null && !keepdims) { const v = result.item(); result.dispose(); return v; }
  return result;
}

export async function argmin(a: NDArray, axis: number | null = null,
                             keepdims = false): Promise<NDArray | number> {
  const axisVal = axis === null ? -2147483648 : axis;
  const resultPtr = a._wasmModule._ndarray_argmin(a._wasmPtr, axisVal, keepdims);
  if (resultPtr === 0) throw new Error('argmin failed');
  const result = NDArray._fromPtr(resultPtr, a._wasmModule);
  if (axis === null && !keepdims) { const v = result.item(); result.dispose(); return v; }
  return result;
}

export async function searchsorted(a: NDArray, v: NDArray | number | number[],
                                   side: 'left' | 'right' = 'left'): Promise<NDArray | number> {
  let vArr: NDArray, disposeV = false;
  if (typeof v === 'number') { vArr = await NDArray.fromArray([v]); disposeV = true; }
  else if (Array.isArray(v)) { vArr = await NDArray.fromArray(v); disposeV = true; }
  else vArr = v;

  const sideVal = side === 'right' ? 1 : 0;
  const resultPtr = a._wasmModule._ndarray_searchsorted(a._wasmPtr, vArr._wasmPtr, sideVal, 0);
  if (disposeV) vArr.dispose();
  if (resultPtr === 0) throw new Error('searchsorted failed');
  const result = NDArray._fromPtr(resultPtr, a._wasmModule);
  if (typeof v === 'number') { const val = result.item(); result.dispose(); return val; }
  return result;
}
```

---

## File Changes Summary

### New Files
```
src/wasm/sorting.h      # Sorting declarations
src/wasm/sorting.c      # Introsort, heapsort, argsort, partition
src/wasm/searching.h    # Searching declarations
src/wasm/searching.c    # argmax, argmin, searchsorted
src/wasm/statistics.h   # Statistics declarations
src/wasm/statistics.c   # sum, mean, var, std, min, max, median

src/ts/sorting.ts       # sort, argsort, partition wrappers
src/ts/statistics.ts    # Statistics and searching wrappers
```

### Files to Modify
```
src/ts/indexing.ts      # Remove argsort, argmax, argmin (replaced)
src/ts/types.ts         # Add WASM function declarations
src/ts/index.ts         # Export new functions
scripts/build-wasm.sh   # Add new source files and exports
```

---

## Build Script Updates

Add to `scripts/build-wasm.sh`:

```bash
# Source files
emcc \
    "$SRC_DIR/ndarray.c" \
    "$SRC_DIR/dtype.c" \
    "$SRC_DIR/pairwise_sum.c" \
    "$SRC_DIR/broadcast.c" \
    "$SRC_DIR/indexing.c" \
    "$SRC_DIR/sorting.c" \
    "$SRC_DIR/searching.c" \
    "$SRC_DIR/statistics.c" \
    ...

# EXPORTED_FUNCTIONS additions:
"_ndarray_sort",
"_ndarray_sort_copy",
"_ndarray_argsort",
"_ndarray_partition",
"_ndarray_argpartition",
"_ndarray_argmax",
"_ndarray_argmin",
"_ndarray_searchsorted",
"_ndarray_sum_axis",
"_ndarray_mean_axis",
"_ndarray_var_axis",
"_ndarray_std_axis",
"_ndarray_min_axis",
"_ndarray_max_axis",
"_ndarray_median",
"_ndarray_percentile",
"_ndarray_quantile",
"_ndarray_nansum",
"_ndarray_nanmean",
"_ndarray_nanvar",
"_ndarray_nanstd"
```

---

## Implementation Order

```
Week 1: Sorting Core
├── Day 1: C sorting algorithms (introsort, heapsort)
├── Day 2: ndarray_sort(), ndarray_sort_copy()
├── Day 3: ndarray_argsort() with axis support
├── Day 4: TS wrappers + tests
└── Day 5: partition(), argpartition()

Week 2: Searching & Basic Stats
├── Day 1: ndarray_argmax(), ndarray_argmin() with keepdims
├── Day 2: ndarray_searchsorted()
├── Day 3: ndarray_sum_axis() using pairwise_sum
├── Day 4: ndarray_mean_axis(), min_axis(), max_axis()
└── Day 5: TS wrappers + tests

Week 3: Advanced Statistics
├── Day 1: ndarray_var_axis() with ddof
├── Day 2: ndarray_std_axis()
├── Day 3: ndarray_median()
├── Day 4: ndarray_percentile(), ndarray_quantile()
└── Day 5: TS wrappers + tests

Week 4: NaN-aware & Polish
├── Day 1: ndarray_nansum(), ndarray_nanmean()
├── Day 2: ndarray_nanvar(), ndarray_nanstd()
├── Day 3: Comprehensive tests
├── Day 4: Edge cases, error handling
└── Day 5: Documentation
```

---

## Verification Plan

```bash
npm run build
npm test

# Phase 6 tests should pass:

# Sorting
✓ sort() sorts along last axis by default
✓ sort() with axis=None flattens and sorts
✓ sort() handles negative axis
✓ argsort() returns indices that would sort

# Searching
✓ argmax() returns flat index when axis=None
✓ argmax() with axis returns array of indices
✓ argmax() with keepdims preserves dimensions
✓ searchsorted() finds correct insertion points

# Statistics
✓ sum() with axis reduces along axis
✓ mean() computes arithmetic mean
✓ var() with ddof=0 computes population variance
✓ var() with ddof=1 computes sample variance
✓ std() is sqrt(var)
✓ median() finds middle value
```

### NumPy Test Vector Generation

```python
# tests/python/generate_phase6_tests.py
import numpy as np
import json

tests = {
    "sort": [
        {"input": [3, 1, 4, 1, 5], "axis": -1, "expected": [1, 1, 3, 4, 5]},
        {"input": [[3, 1], [2, 4]], "axis": 0, "expected": [[2, 1], [3, 4]]},
        {"input": [[3, 1], [2, 4]], "axis": 1, "expected": [[1, 3], [2, 4]]},
    ],
    "argsort": [
        {"input": [3, 1, 2], "expected": [1, 2, 0]},
    ],
    "argmax": [
        {"input": [1, 3, 2], "axis": None, "expected": 1},
        {"input": [[1, 2], [4, 3]], "axis": 0, "expected": [1, 1]},
        {"input": [[1, 2], [4, 3]], "axis": 1, "expected": [1, 0]},
    ],
    "mean": [
        {"input": [1, 2, 3, 4], "expected": 2.5},
        {"input": [[1, 2], [3, 4]], "axis": 0, "expected": [2.0, 3.0]},
    ],
    "std": [
        {"input": [1, 2, 3, 4, 5], "ddof": 0, "expected": 1.4142135623730951},
        {"input": [1, 2, 3, 4, 5], "ddof": 1, "expected": 1.5811388300841898},
    ],
    "median": [
        {"input": [1, 2, 3, 4, 5], "expected": 3.0},
        {"input": [1, 2, 3, 4], "expected": 2.5},
    ],
}

with open("tests/fixtures/phase6_vectors.json", "w") as f:
    json.dump(tests, f, indent=2)
```

---

## Critical Dependencies for Phase 7

Phase 6 completion enables:
- **numpy.linalg**: Sorting for eigenvalue ordering
- **numpy.fft**: Statistics for signal analysis
- **numpy.random**: Statistics for distribution testing

Phase 6 MUST be complete before proceeding to Phase 7.
