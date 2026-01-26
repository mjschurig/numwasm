# Phase 8: Set Operations Implementation Plan

Phase 8 implements set-based array operations for NumJS-WASM. These functions enable membership testing, finding unique elements, and computing set unions, intersections, and differences on arrays.

---

## Current State (Phases 1-7)

```
src/wasm/
├── ndarray.h/c        # Core NDArray with views, slicing, element access
├── dtype.h/c          # DType system with type promotion
├── broadcast.h/c      # Broadcasting infrastructure
├── indexing.h/c       # take, put, nonzero, where, compress, extract
└── pairwise_sum.h/c   # Accurate summation

src/ts/
├── types.ts           # DType enum, WasmModule interface
├── NDArray.ts         # Factory methods, shape ops, element access
├── dtype.ts           # Type promotion utilities
├── broadcast.ts       # broadcastTo, broadcastArrays, broadcastShapes
├── indexing.ts        # take, put, nonzero, where, compress, extract, argsort
├── slice.ts           # Slice class, ellipsis, newaxis
├── iterators.ts       # FlatIterator, nditer, ndenumerate, ndindex
├── wasm-loader.ts     # WASM module loading
└── index.ts           # Public exports
```

**Existing Infrastructure Used by Set Operations:**
- `argsort()` - TypeScript implementation for 1D arrays (indexing.ts:386)
- `take()` - Select elements by indices along axis
- `nonzero()` - Find indices of non-zero elements
- `where()` - Conditional element selection
- `flatten()` / `ravel()` - Convert to 1D arrays
- Broadcasting infrastructure for element-wise operations

---

## Phase 8 Implementation Tree

```
PHASE 8: SET OPERATIONS
│
├── 8.1 Sorting Infrastructure (C + TypeScript)
│   ├── 8.1.1 C: ndarray_sort() → in-place sort along axis
│   ├── 8.1.2 C: ndarray_argsort() → indices that would sort
│   ├── 8.1.3 C: quicksort/mergesort kernels per dtype
│   └── 8.1.4 TS: sort(), argsort() with full axis support
│
│   Dependencies: Level 1-2 core
│
├── 8.2 Unique Function (C + TypeScript) ⭐ FOUNDATIONAL
│   ├── 8.2.1 C: ndarray_unique() → sorted unique elements
│   ├── 8.2.2 C: ndarray_unique_with_indices() → with return_index
│   ├── 8.2.3 C: ndarray_unique_with_inverse() → with return_inverse
│   ├── 8.2.4 C: ndarray_unique_with_counts() → with return_counts
│   ├── 8.2.5 TS: unique() function with all optional returns
│   ├── 8.2.6 TS: uniqueAll() → Array API compatible
│   ├── 8.2.7 TS: uniqueCounts() → Array API compatible
│   ├── 8.2.8 TS: uniqueInverse() → Array API compatible
│   └── 8.2.9 TS: uniqueValues() → Array API compatible
│
│   Dependencies: 8.1.*
│
├── 8.3 Set Combination Operations (C + TypeScript)
│   ├── 8.3.1 C: ndarray_union1d() → unique sorted union
│   ├── 8.3.2 C: ndarray_intersect1d() → common elements
│   ├── 8.3.3 C: ndarray_setdiff1d() → elements in ar1 not in ar2
│   ├── 8.3.4 C: ndarray_setxor1d() → exclusive or (symmetric difference)
│   ├── 8.3.5 TS: union1d() wrapper
│   ├── 8.3.6 TS: intersect1d() wrapper with return_indices option
│   ├── 8.3.7 TS: setdiff1d() wrapper
│   └── 8.3.8 TS: setxor1d() wrapper
│
│   Dependencies: 8.2.*
│
├── 8.4 Membership Testing (C + TypeScript)
│   ├── 8.4.1 C: ndarray_isin_sort() → sort-based membership
│   ├── 8.4.2 C: ndarray_isin_table() → table-based for integers
│   ├── 8.4.3 C: ndarray_in1d() → flattened membership (deprecated alias)
│   ├── 8.4.4 TS: isin() with kind parameter ('sort' | 'table' | null)
│   └── 8.4.5 TS: in1d() deprecated wrapper
│
│   Dependencies: 8.1.*, 8.2.*
│
└── 8.5 Difference Operations (TypeScript)
    ├── 8.5.1 TS: ediff1d() → differences between consecutive elements
    └── 8.5.2 TS: Integrate with existing diff() if available

    Dependencies: Level 1-2 slicing
```

---

## Detailed Implementation Specifications

### 8.1 Sorting Infrastructure (C)

**File:** `src/wasm/sorting.h` (new file)

```c
#ifndef NUMJS_SORTING_H
#define NUMJS_SORTING_H

#include "ndarray.h"

/**
 * Sort kind enumeration
 */
typedef enum {
    SORT_QUICKSORT = 0,
    SORT_MERGESORT = 1,
    SORT_HEAPSORT = 2,
    SORT_STABLE = 3  /* alias for mergesort */
} SortKind;

/**
 * Sort an array in-place along an axis.
 *
 * Reference: numpy/_core/src/npysort/quicksort.cpp
 *            numpy/_core/src/npysort/mergesort.cpp
 *
 * @param arr  Array to sort (modified in place)
 * @param axis Axis along which to sort (-1 for last axis)
 * @param kind Sorting algorithm to use
 * @return     0 on success, -1 on error
 */
int ndarray_sort(NDArray* arr, int32_t axis, SortKind kind);

/**
 * Return indices that would sort an array along an axis.
 *
 * Reference: numpy/_core/src/multiarray/item_selection.c (PyArray_ArgSort)
 *
 * @param arr  Source array
 * @param axis Axis along which to sort (-1 for last axis)
 * @param kind Sorting algorithm to use
 * @return     New array of indices (Int32/Int64) or NULL on error
 */
NDArray* ndarray_argsort(NDArray* arr, int32_t axis, SortKind kind);

/**
 * Lexicographic sort by multiple keys.
 *
 * Reference: numpy/_core/src/multiarray/item_selection.c (PyArray_LexSort)
 *
 * @param keys      Array of key arrays (last key is primary)
 * @param num_keys  Number of key arrays
 * @param axis      Axis along which to sort
 * @return          Indices that would sort by keys or NULL on error
 */
NDArray* ndarray_lexsort(NDArray** keys, int32_t num_keys, int32_t axis);

#endif /* NUMJS_SORTING_H */
```

**File:** `src/wasm/sorting.c` (new file)

```c
#include "sorting.h"
#include <stdlib.h>
#include <string.h>

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#define EXPORT EMSCRIPTEN_KEEPALIVE
#else
#define EXPORT
#endif

/* ============ Comparison Functions per DType ============ */

/**
 * Generic comparison function signature for qsort_r style sorting.
 * Returns < 0 if a < b, 0 if a == b, > 0 if a > b
 */
typedef int (*CompareFunc)(const void* a, const void* b, void* context);

static int compare_float64(const void* a, const void* b, void* ctx) {
    (void)ctx;
    double va = *(const double*)a;
    double vb = *(const double*)b;
    /* Handle NaN: NaN sorts to end (NumPy convention) */
    if (va != va) return (vb != vb) ? 0 : 1;  /* NaN vs NaN or value */
    if (vb != vb) return -1;  /* value vs NaN */
    return (va < vb) ? -1 : (va > vb) ? 1 : 0;
}

static int compare_float32(const void* a, const void* b, void* ctx) {
    (void)ctx;
    float va = *(const float*)a;
    float vb = *(const float*)b;
    if (va != va) return (vb != vb) ? 0 : 1;
    if (vb != vb) return -1;
    return (va < vb) ? -1 : (va > vb) ? 1 : 0;
}

static int compare_int64(const void* a, const void* b, void* ctx) {
    (void)ctx;
    int64_t va = *(const int64_t*)a;
    int64_t vb = *(const int64_t*)b;
    return (va < vb) ? -1 : (va > vb) ? 1 : 0;
}

static int compare_int32(const void* a, const void* b, void* ctx) {
    (void)ctx;
    int32_t va = *(const int32_t*)a;
    int32_t vb = *(const int32_t*)b;
    return (va < vb) ? -1 : (va > vb) ? 1 : 0;
}

static int compare_uint64(const void* a, const void* b, void* ctx) {
    (void)ctx;
    uint64_t va = *(const uint64_t*)a;
    uint64_t vb = *(const uint64_t*)b;
    return (va < vb) ? -1 : (va > vb) ? 1 : 0;
}

static int compare_uint32(const void* a, const void* b, void* ctx) {
    (void)ctx;
    uint32_t va = *(const uint32_t*)a;
    uint32_t vb = *(const uint32_t*)b;
    return (va < vb) ? -1 : (va > vb) ? 1 : 0;
}

static int compare_bool(const void* a, const void* b, void* ctx) {
    (void)ctx;
    uint8_t va = *(const uint8_t*)a;
    uint8_t vb = *(const uint8_t*)b;
    return (int)va - (int)vb;
}

/**
 * Get comparison function for a dtype.
 */
static CompareFunc get_compare_func(DType dtype) {
    switch (dtype) {
        case DTYPE_FLOAT64: return compare_float64;
        case DTYPE_FLOAT32: return compare_float32;
        case DTYPE_INT64:   return compare_int64;
        case DTYPE_INT32:   return compare_int32;
        case DTYPE_UINT64:  return compare_uint64;
        case DTYPE_UINT32:  return compare_uint32;
        case DTYPE_BOOL:    return compare_bool;
        /* Add more types as needed */
        default: return NULL;
    }
}

/* ============ Quicksort Implementation ============ */

/**
 * Partition function for quicksort.
 */
static size_t partition(void* data, size_t lo, size_t hi,
                        size_t elem_size, CompareFunc cmp, void* ctx) {
    char* arr = (char*)data;
    char* pivot = arr + hi * elem_size;
    size_t i = lo;
    char* temp = malloc(elem_size);

    for (size_t j = lo; j < hi; j++) {
        if (cmp(arr + j * elem_size, pivot, ctx) < 0) {
            /* Swap arr[i] and arr[j] */
            memcpy(temp, arr + i * elem_size, elem_size);
            memcpy(arr + i * elem_size, arr + j * elem_size, elem_size);
            memcpy(arr + j * elem_size, temp, elem_size);
            i++;
        }
    }
    /* Swap arr[i] and arr[hi] (pivot) */
    memcpy(temp, arr + i * elem_size, elem_size);
    memcpy(arr + i * elem_size, arr + hi * elem_size, elem_size);
    memcpy(arr + hi * elem_size, temp, elem_size);

    free(temp);
    return i;
}

/**
 * Recursive quicksort implementation.
 */
static void quicksort_impl(void* data, size_t lo, size_t hi,
                            size_t elem_size, CompareFunc cmp, void* ctx) {
    if (lo < hi && hi != (size_t)-1) {
        size_t p = partition(data, lo, hi, elem_size, cmp, ctx);
        if (p > 0) quicksort_impl(data, lo, p - 1, elem_size, cmp, ctx);
        quicksort_impl(data, p + 1, hi, elem_size, cmp, ctx);
    }
}

/* ============ Mergesort Implementation (Stable) ============ */

/**
 * Merge two sorted halves.
 */
static void merge(void* data, void* temp, size_t lo, size_t mid, size_t hi,
                  size_t elem_size, CompareFunc cmp, void* ctx) {
    char* arr = (char*)data;
    char* tmp = (char*)temp;

    /* Copy to temp */
    memcpy(tmp + lo * elem_size, arr + lo * elem_size, (hi - lo + 1) * elem_size);

    size_t i = lo, j = mid + 1, k = lo;

    while (i <= mid && j <= hi) {
        if (cmp(tmp + i * elem_size, tmp + j * elem_size, ctx) <= 0) {
            memcpy(arr + k * elem_size, tmp + i * elem_size, elem_size);
            i++;
        } else {
            memcpy(arr + k * elem_size, tmp + j * elem_size, elem_size);
            j++;
        }
        k++;
    }

    /* Copy remaining elements */
    while (i <= mid) {
        memcpy(arr + k * elem_size, tmp + i * elem_size, elem_size);
        i++; k++;
    }
    /* j elements are already in place */
}

/**
 * Recursive mergesort implementation.
 */
static void mergesort_impl(void* data, void* temp, size_t lo, size_t hi,
                            size_t elem_size, CompareFunc cmp, void* ctx) {
    if (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        mergesort_impl(data, temp, lo, mid, elem_size, cmp, ctx);
        mergesort_impl(data, temp, mid + 1, hi, elem_size, cmp, ctx);
        merge(data, temp, lo, mid, hi, elem_size, cmp, ctx);
    }
}

/* ============ Public API ============ */

EXPORT int ndarray_sort(NDArray* arr, int32_t axis, SortKind kind) {
    if (!arr) return -1;

    /* Normalize axis */
    if (axis < 0) axis += arr->ndim;
    if (axis < 0 || axis >= arr->ndim) return -1;

    /* Must be contiguous for in-place sort */
    if (!ndarray_is_c_contiguous(arr)) return -1;

    CompareFunc cmp = get_compare_func(arr->dtype);
    if (!cmp) return -1;

    size_t elem_size = dtype_size(arr->dtype);

    if (arr->ndim == 1) {
        /* Simple 1D case */
        size_t n = arr->shape[0];
        if (n <= 1) return 0;

        if (kind == SORT_QUICKSORT || kind == SORT_HEAPSORT) {
            quicksort_impl(arr->data, 0, n - 1, elem_size, cmp, NULL);
        } else {
            /* Mergesort (stable) */
            void* temp = malloc(n * elem_size);
            if (!temp) return -1;
            mergesort_impl(arr->data, temp, 0, n - 1, elem_size, cmp, NULL);
            free(temp);
        }
        return 0;
    }

    /* Multi-dimensional: sort along axis */
    /* Compute iteration parameters */
    size_t axis_len = arr->shape[axis];
    size_t outer_size = 1;
    size_t inner_size = 1;

    for (int i = 0; i < axis; i++) {
        outer_size *= arr->shape[i];
    }
    for (int i = axis + 1; i < arr->ndim; i++) {
        inner_size *= arr->shape[i];
    }

    /* Allocate temp buffer for extraction and sorting */
    void* line = malloc(axis_len * elem_size);
    void* temp = malloc(axis_len * elem_size);
    if (!line || !temp) {
        free(line);
        free(temp);
        return -1;
    }

    char* data = (char*)arr->data;
    size_t axis_stride = arr->strides[axis];

    for (size_t outer = 0; outer < outer_size; outer++) {
        for (size_t inner = 0; inner < inner_size; inner++) {
            /* Extract line along axis */
            for (size_t k = 0; k < axis_len; k++) {
                size_t offset = (outer * axis_len + k) * inner_size + inner;
                offset *= elem_size;
                memcpy((char*)line + k * elem_size, data + offset, elem_size);
            }

            /* Sort the line */
            if (kind == SORT_QUICKSORT || kind == SORT_HEAPSORT) {
                quicksort_impl(line, 0, axis_len - 1, elem_size, cmp, NULL);
            } else {
                mergesort_impl(line, temp, 0, axis_len - 1, elem_size, cmp, NULL);
            }

            /* Write back */
            for (size_t k = 0; k < axis_len; k++) {
                size_t offset = (outer * axis_len + k) * inner_size + inner;
                offset *= elem_size;
                memcpy(data + offset, (char*)line + k * elem_size, elem_size);
            }
        }
    }

    free(line);
    free(temp);
    return 0;
}

EXPORT NDArray* ndarray_argsort(NDArray* arr, int32_t axis, SortKind kind) {
    if (!arr) return NULL;

    /* Normalize axis */
    if (axis < 0) axis += arr->ndim;
    if (axis < 0 || axis >= arr->ndim) return NULL;

    /* Create output array of indices */
    NDArray* indices = ndarray_empty(arr->ndim, arr->shape, DTYPE_INT32);
    if (!indices) return NULL;

    /* Implementation similar to sort, but track indices */
    /* ... (detailed implementation follows similar pattern) */

    return indices;
}
```

---

### 8.2 Unique Function (C)

**File:** `src/wasm/setops.h` (new file)

```c
#ifndef NUMJS_SETOPS_H
#define NUMJS_SETOPS_H

#include "ndarray.h"

/**
 * Result structure for unique with multiple return values.
 * Caller must free non-NULL arrays.
 */
typedef struct {
    NDArray* values;         /* Sorted unique values */
    NDArray* indices;        /* First occurrence indices (if requested) */
    NDArray* inverse;        /* Indices to reconstruct input (if requested) */
    NDArray* counts;         /* Count of each unique value (if requested) */
} UniqueResult;

/**
 * Find unique elements of an array.
 *
 * Reference: numpy/lib/_arraysetops_impl.py (unique, _unique1d)
 *
 * Algorithm:
 * 1. Flatten input array
 * 2. Compute argsort (mergesort if return_index for stability)
 * 3. Create sorted view using argsort indices
 * 4. Create mask where adjacent elements differ
 * 5. Extract unique elements using mask
 * 6. If return_index: extract corresponding argsort indices
 * 7. If return_inverse: cumsum of mask, then scatter back via argsort
 * 8. If return_counts: diff of nonzero indices in mask
 *
 * @param arr            Input array (will be flattened)
 * @param return_index   Return indices of first occurrences
 * @param return_inverse Return indices to reconstruct input
 * @param return_counts  Return counts of each unique value
 * @param equal_nan      Treat NaN values as equal (default true)
 * @return               UniqueResult structure (caller must free)
 */
UniqueResult* ndarray_unique(NDArray* arr,
                              bool return_index,
                              bool return_inverse,
                              bool return_counts,
                              bool equal_nan);

/**
 * Free a UniqueResult structure and its arrays.
 */
void unique_result_free(UniqueResult* result);

/**
 * Simplified unique: just return sorted unique values.
 */
NDArray* ndarray_unique_values(NDArray* arr);

#endif /* NUMJS_SETOPS_H */
```

**File:** `src/wasm/setops.c` (new file - partial)

```c
#include "setops.h"
#include "sorting.h"
#include <stdlib.h>
#include <string.h>

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#define EXPORT EMSCRIPTEN_KEEPALIVE
#else
#define EXPORT
#endif

/**
 * Check if two elements are equal for a given dtype.
 * Handles NaN comparison based on equal_nan flag.
 */
static bool elements_equal(const void* a, const void* b, DType dtype, bool equal_nan) {
    switch (dtype) {
        case DTYPE_FLOAT64: {
            double va = *(const double*)a;
            double vb = *(const double*)b;
            if (equal_nan && va != va && vb != vb) return true;  /* Both NaN */
            return va == vb;
        }
        case DTYPE_FLOAT32: {
            float va = *(const float*)a;
            float vb = *(const float*)b;
            if (equal_nan && va != va && vb != vb) return true;
            return va == vb;
        }
        case DTYPE_INT64:
            return *(const int64_t*)a == *(const int64_t*)b;
        case DTYPE_INT32:
            return *(const int32_t*)a == *(const int32_t*)b;
        case DTYPE_UINT64:
            return *(const uint64_t*)a == *(const uint64_t*)b;
        case DTYPE_UINT32:
            return *(const uint32_t*)a == *(const uint32_t*)b;
        case DTYPE_BOOL:
            return *(const uint8_t*)a == *(const uint8_t*)b;
        default:
            return false;
    }
}

EXPORT UniqueResult* ndarray_unique(NDArray* arr,
                                     bool return_index,
                                     bool return_inverse,
                                     bool return_counts,
                                     bool equal_nan) {
    if (!arr) return NULL;

    /* Flatten the input */
    NDArray* flat = ndarray_flatten(arr);
    if (!flat) return NULL;

    size_t n = flat->size;
    size_t elem_size = dtype_size(flat->dtype);

    /* Handle empty array */
    if (n == 0) {
        UniqueResult* result = calloc(1, sizeof(UniqueResult));
        result->values = ndarray_empty(1, (int32_t[]){0}, flat->dtype);
        if (return_index) result->indices = ndarray_empty(1, (int32_t[]){0}, DTYPE_INT32);
        if (return_inverse) result->inverse = ndarray_empty(1, (int32_t[]){0}, DTYPE_INT32);
        if (return_counts) result->counts = ndarray_empty(1, (int32_t[]){0}, DTYPE_INT32);
        ndarray_free(flat);
        return result;
    }

    /* Compute argsort */
    SortKind kind = return_index ? SORT_MERGESORT : SORT_QUICKSORT;
    NDArray* perm = ndarray_argsort(flat, 0, kind);
    if (!perm) {
        ndarray_free(flat);
        return NULL;
    }

    /* Create sorted auxiliary array using permutation */
    NDArray* aux = ndarray_take(flat, perm, 0, CLIP_RAISE);
    if (!aux) {
        ndarray_free(flat);
        ndarray_free(perm);
        return NULL;
    }

    /* Create mask: mask[0] = true, mask[i] = aux[i] != aux[i-1] */
    uint8_t* mask = malloc(n);
    if (!mask) {
        ndarray_free(flat);
        ndarray_free(perm);
        ndarray_free(aux);
        return NULL;
    }

    mask[0] = 1;  /* First element is always unique */
    char* aux_data = (char*)aux->data;

    for (size_t i = 1; i < n; i++) {
        mask[i] = !elements_equal(
            aux_data + i * elem_size,
            aux_data + (i - 1) * elem_size,
            flat->dtype,
            equal_nan
        );
    }

    /* Count unique elements */
    size_t num_unique = 0;
    for (size_t i = 0; i < n; i++) {
        num_unique += mask[i];
    }

    /* Extract unique values */
    NDArray* values = ndarray_empty(1, (int32_t[]){(int32_t)num_unique}, flat->dtype);
    if (!values) {
        free(mask);
        ndarray_free(flat);
        ndarray_free(perm);
        ndarray_free(aux);
        return NULL;
    }

    char* values_data = (char*)values->data;
    size_t unique_idx = 0;
    for (size_t i = 0; i < n; i++) {
        if (mask[i]) {
            memcpy(values_data + unique_idx * elem_size,
                   aux_data + i * elem_size,
                   elem_size);
            unique_idx++;
        }
    }

    /* Build result */
    UniqueResult* result = calloc(1, sizeof(UniqueResult));
    result->values = values;

    /* Return indices: perm[mask] */
    if (return_index) {
        result->indices = ndarray_empty(1, (int32_t[]){(int32_t)num_unique}, DTYPE_INT32);
        int32_t* indices_data = (int32_t*)result->indices->data;
        int32_t* perm_data = (int32_t*)perm->data;
        unique_idx = 0;
        for (size_t i = 0; i < n; i++) {
            if (mask[i]) {
                indices_data[unique_idx++] = perm_data[i];
            }
        }
    }

    /* Return inverse: indices to reconstruct original from unique */
    if (return_inverse) {
        result->inverse = ndarray_empty(1, (int32_t[]){(int32_t)n}, DTYPE_INT32);
        int32_t* inverse_data = (int32_t*)result->inverse->data;
        int32_t* perm_data = (int32_t*)perm->data;

        /* imask = cumsum(mask) - 1 */
        int32_t* imask = malloc(n * sizeof(int32_t));
        int32_t cumsum = 0;
        for (size_t i = 0; i < n; i++) {
            cumsum += mask[i];
            imask[i] = cumsum - 1;
        }

        /* inverse[perm] = imask */
        for (size_t i = 0; i < n; i++) {
            inverse_data[perm_data[i]] = imask[i];
        }

        free(imask);
    }

    /* Return counts: difference of indices where mask is true */
    if (return_counts) {
        result->counts = ndarray_empty(1, (int32_t[]){(int32_t)num_unique}, DTYPE_INT32);
        int32_t* counts_data = (int32_t*)result->counts->data;

        /* Find indices where mask is true, plus n at the end */
        size_t* idx = malloc((num_unique + 1) * sizeof(size_t));
        unique_idx = 0;
        for (size_t i = 0; i < n; i++) {
            if (mask[i]) {
                idx[unique_idx++] = i;
            }
        }
        idx[num_unique] = n;

        /* counts = diff(idx) */
        for (size_t i = 0; i < num_unique; i++) {
            counts_data[i] = (int32_t)(idx[i + 1] - idx[i]);
        }

        free(idx);
    }

    /* Cleanup */
    free(mask);
    ndarray_free(flat);
    ndarray_free(perm);
    ndarray_free(aux);

    return result;
}

EXPORT void unique_result_free(UniqueResult* result) {
    if (!result) return;
    if (result->values) ndarray_free(result->values);
    if (result->indices) ndarray_free(result->indices);
    if (result->inverse) ndarray_free(result->inverse);
    if (result->counts) ndarray_free(result->counts);
    free(result);
}

EXPORT NDArray* ndarray_unique_values(NDArray* arr) {
    UniqueResult* result = ndarray_unique(arr, false, false, false, true);
    if (!result) return NULL;
    NDArray* values = result->values;
    result->values = NULL;  /* Prevent freeing */
    unique_result_free(result);
    return values;
}
```

---

### 8.3 Set Combination Operations (C)

**Add to:** `src/wasm/setops.h`

```c
/**
 * Find the union of two arrays.
 * Returns unique, sorted values in either array.
 *
 * Reference: numpy/lib/_arraysetops_impl.py (union1d)
 * Algorithm: concatenate then unique
 *
 * @param ar1 First input array
 * @param ar2 Second input array
 * @return    Sorted unique union or NULL on error
 */
NDArray* ndarray_union1d(NDArray* ar1, NDArray* ar2);

/**
 * Find the intersection of two arrays.
 * Returns sorted, unique values in both arrays.
 *
 * Reference: numpy/lib/_arraysetops_impl.py (intersect1d)
 * Algorithm:
 * 1. If not assume_unique, compute unique of each
 * 2. Concatenate arrays
 * 3. Sort (mergesort if returning indices)
 * 4. Find adjacent equal elements
 *
 * @param ar1            First input array
 * @param ar2            Second input array
 * @param assume_unique  If true, skip uniquifying inputs
 * @param return_indices If true, also return indices in original arrays
 * @param out_ar1_idx    Output: indices in ar1 (if return_indices)
 * @param out_ar2_idx    Output: indices in ar2 (if return_indices)
 * @return               Sorted intersection or NULL on error
 */
NDArray* ndarray_intersect1d(NDArray* ar1, NDArray* ar2,
                              bool assume_unique,
                              bool return_indices,
                              NDArray** out_ar1_idx,
                              NDArray** out_ar2_idx);

/**
 * Find the set difference (ar1 - ar2).
 * Returns values in ar1 that are not in ar2.
 *
 * Reference: numpy/lib/_arraysetops_impl.py (setdiff1d)
 * Algorithm: unique(ar1)[~isin(unique(ar1), unique(ar2))]
 *
 * @param ar1           First input array
 * @param ar2           Second input array
 * @param assume_unique If true, skip uniquifying inputs
 * @return              Set difference or NULL on error
 */
NDArray* ndarray_setdiff1d(NDArray* ar1, NDArray* ar2, bool assume_unique);

/**
 * Find the set exclusive-or (symmetric difference).
 * Returns values in exactly one of the arrays.
 *
 * Reference: numpy/lib/_arraysetops_impl.py (setxor1d)
 * Algorithm:
 * 1. If not assume_unique, compute unique of each
 * 2. Concatenate and sort
 * 3. Create flag array: [True, aux[1:] != aux[:-1], True]
 * 4. Return aux[flag[1:] & flag[:-1]]
 *
 * @param ar1           First input array
 * @param ar2           Second input array
 * @param assume_unique If true, skip uniquifying inputs
 * @return              Symmetric difference or NULL on error
 */
NDArray* ndarray_setxor1d(NDArray* ar1, NDArray* ar2, bool assume_unique);
```

**Add to:** `src/wasm/setops.c`

```c
EXPORT NDArray* ndarray_union1d(NDArray* ar1, NDArray* ar2) {
    if (!ar1 || !ar2) return NULL;

    /* Concatenate flattened arrays */
    NDArray* flat1 = ndarray_flatten(ar1);
    NDArray* flat2 = ndarray_flatten(ar2);
    if (!flat1 || !flat2) {
        ndarray_free(flat1);
        ndarray_free(flat2);
        return NULL;
    }

    NDArray* concat = ndarray_concatenate(flat1, flat2, 0);
    ndarray_free(flat1);
    ndarray_free(flat2);

    if (!concat) return NULL;

    /* Return unique values */
    NDArray* result = ndarray_unique_values(concat);
    ndarray_free(concat);

    return result;
}

EXPORT NDArray* ndarray_intersect1d(NDArray* ar1, NDArray* ar2,
                                     bool assume_unique,
                                     bool return_indices,
                                     NDArray** out_ar1_idx,
                                     NDArray** out_ar2_idx) {
    if (!ar1 || !ar2) return NULL;

    NDArray* u1 = NULL;
    NDArray* u2 = NULL;
    NDArray* ind1 = NULL;
    NDArray* ind2 = NULL;

    /* Get unique arrays (and indices if needed) */
    if (!assume_unique) {
        if (return_indices) {
            UniqueResult* r1 = ndarray_unique(ar1, true, false, false, true);
            UniqueResult* r2 = ndarray_unique(ar2, true, false, false, true);
            if (!r1 || !r2) {
                unique_result_free(r1);
                unique_result_free(r2);
                return NULL;
            }
            u1 = r1->values; r1->values = NULL;
            u2 = r2->values; r2->values = NULL;
            ind1 = r1->indices; r1->indices = NULL;
            ind2 = r2->indices; r2->indices = NULL;
            unique_result_free(r1);
            unique_result_free(r2);
        } else {
            u1 = ndarray_unique_values(ar1);
            u2 = ndarray_unique_values(ar2);
        }
    } else {
        u1 = ndarray_flatten(ar1);
        u2 = ndarray_flatten(ar2);
    }

    if (!u1 || !u2) {
        ndarray_free(u1);
        ndarray_free(u2);
        ndarray_free(ind1);
        ndarray_free(ind2);
        return NULL;
    }

    size_t n1 = u1->size;
    size_t n2 = u2->size;
    size_t n = n1 + n2;

    /* Handle empty case */
    if (n == 0) {
        if (out_ar1_idx) *out_ar1_idx = ndarray_empty(1, (int32_t[]){0}, DTYPE_INT32);
        if (out_ar2_idx) *out_ar2_idx = ndarray_empty(1, (int32_t[]){0}, DTYPE_INT32);
        NDArray* result = ndarray_empty(1, (int32_t[]){0}, u1->dtype);
        ndarray_free(u1);
        ndarray_free(u2);
        ndarray_free(ind1);
        ndarray_free(ind2);
        return result;
    }

    /* Concatenate */
    NDArray* aux = ndarray_concatenate(u1, u2, 0);
    if (!aux) {
        ndarray_free(u1);
        ndarray_free(u2);
        ndarray_free(ind1);
        ndarray_free(ind2);
        return NULL;
    }

    /* Sort (with indices if returning them) */
    NDArray* sort_idx = NULL;
    if (return_indices) {
        sort_idx = ndarray_argsort(aux, 0, SORT_MERGESORT);
        NDArray* sorted = ndarray_take(aux, sort_idx, 0, CLIP_RAISE);
        ndarray_free(aux);
        aux = sorted;
    } else {
        ndarray_sort(aux, 0, SORT_QUICKSORT);
    }

    /* Find where adjacent elements are equal (intersection) */
    size_t elem_size = dtype_size(aux->dtype);
    char* aux_data = (char*)aux->data;

    /* Count intersections */
    size_t count = 0;
    for (size_t i = 0; i < n - 1; i++) {
        if (elements_equal(aux_data + i * elem_size,
                           aux_data + (i + 1) * elem_size,
                           aux->dtype, true)) {
            count++;
        }
    }

    /* Extract intersection elements */
    NDArray* result = ndarray_empty(1, (int32_t[]){(int32_t)count}, aux->dtype);
    char* result_data = (char*)result->data;

    NDArray* ar1_indices = return_indices ?
        ndarray_empty(1, (int32_t[]){(int32_t)count}, DTYPE_INT32) : NULL;
    NDArray* ar2_indices = return_indices ?
        ndarray_empty(1, (int32_t[]){(int32_t)count}, DTYPE_INT32) : NULL;

    size_t idx = 0;
    for (size_t i = 0; i < n - 1; i++) {
        if (elements_equal(aux_data + i * elem_size,
                           aux_data + (i + 1) * elem_size,
                           aux->dtype, true)) {
            memcpy(result_data + idx * elem_size,
                   aux_data + i * elem_size,
                   elem_size);

            if (return_indices) {
                int32_t* sort_data = (int32_t*)sort_idx->data;
                int32_t sidx1 = sort_data[i];
                int32_t sidx2 = sort_data[i + 1] - (int32_t)n1;

                /* Map back through unique indices if needed */
                if (!assume_unique) {
                    sidx1 = ((int32_t*)ind1->data)[sidx1];
                    sidx2 = ((int32_t*)ind2->data)[sidx2];
                }

                ((int32_t*)ar1_indices->data)[idx] = sidx1;
                ((int32_t*)ar2_indices->data)[idx] = sidx2;
            }
            idx++;
        }
    }

    /* Cleanup */
    ndarray_free(u1);
    ndarray_free(u2);
    ndarray_free(ind1);
    ndarray_free(ind2);
    ndarray_free(aux);
    ndarray_free(sort_idx);

    if (out_ar1_idx) *out_ar1_idx = ar1_indices;
    else ndarray_free(ar1_indices);

    if (out_ar2_idx) *out_ar2_idx = ar2_indices;
    else ndarray_free(ar2_indices);

    return result;
}

EXPORT NDArray* ndarray_setxor1d(NDArray* ar1, NDArray* ar2, bool assume_unique) {
    if (!ar1 || !ar2) return NULL;

    NDArray* u1;
    NDArray* u2;

    if (!assume_unique) {
        u1 = ndarray_unique_values(ar1);
        u2 = ndarray_unique_values(ar2);
    } else {
        u1 = ndarray_flatten(ar1);
        u2 = ndarray_flatten(ar2);
    }

    if (!u1 || !u2) {
        ndarray_free(u1);
        ndarray_free(u2);
        return NULL;
    }

    /* Concatenate */
    NDArray* aux = ndarray_concatenate(u1, u2, 0);
    ndarray_free(u1);
    ndarray_free(u2);

    if (!aux) return NULL;

    size_t n = aux->size;
    if (n == 0) return aux;

    /* Sort */
    ndarray_sort(aux, 0, SORT_QUICKSORT);

    /* Create flag: [True, aux[1:] != aux[:-1], True] */
    uint8_t* flag = malloc(n + 1);
    if (!flag) {
        ndarray_free(aux);
        return NULL;
    }

    flag[0] = 1;
    flag[n] = 1;

    size_t elem_size = dtype_size(aux->dtype);
    char* aux_data = (char*)aux->data;

    for (size_t i = 1; i < n; i++) {
        flag[i] = !elements_equal(aux_data + i * elem_size,
                                   aux_data + (i - 1) * elem_size,
                                   aux->dtype, true);
    }

    /* Count elements where flag[i] & flag[i+1] */
    size_t count = 0;
    for (size_t i = 0; i < n; i++) {
        if (flag[i] && flag[i + 1]) count++;
    }

    /* Extract result */
    NDArray* result = ndarray_empty(1, (int32_t[]){(int32_t)count}, aux->dtype);
    char* result_data = (char*)result->data;

    size_t idx = 0;
    for (size_t i = 0; i < n; i++) {
        if (flag[i] && flag[i + 1]) {
            memcpy(result_data + idx * elem_size,
                   aux_data + i * elem_size,
                   elem_size);
            idx++;
        }
    }

    free(flag);
    ndarray_free(aux);

    return result;
}

EXPORT NDArray* ndarray_setdiff1d(NDArray* ar1, NDArray* ar2, bool assume_unique) {
    if (!ar1 || !ar2) return NULL;

    NDArray* u1;
    NDArray* u2;

    if (!assume_unique) {
        u1 = ndarray_unique_values(ar1);
        u2 = ndarray_unique_values(ar2);
    } else {
        u1 = ndarray_flatten(ar1);
        u2 = ndarray_flatten(ar2);
    }

    if (!u1 || !u2) {
        ndarray_free(u1);
        ndarray_free(u2);
        return NULL;
    }

    /* Use isin to find elements of u1 that are in u2 */
    NDArray* mask = ndarray_isin(u1, u2, true, true);  /* invert=true */
    if (!mask) {
        ndarray_free(u1);
        ndarray_free(u2);
        return NULL;
    }

    /* Extract elements where mask is true */
    NDArray* result = ndarray_compress(mask, u1, 0);

    ndarray_free(u1);
    ndarray_free(u2);
    ndarray_free(mask);

    return result;
}
```

---

### 8.4 Membership Testing (C)

**Add to:** `src/wasm/setops.h`

```c
/**
 * Method selection for isin
 */
typedef enum {
    ISIN_AUTO = 0,     /* Automatically select based on memory/performance */
    ISIN_SORT = 1,     /* Use sorting-based method */
    ISIN_TABLE = 2     /* Use lookup table (integers only) */
} IsinKind;

/**
 * Test whether elements of ar1 are in ar2.
 * Returns boolean array of same shape as ar1.
 *
 * Reference: numpy/lib/_arraysetops_impl.py (isin, _isin)
 *
 * @param ar1           Elements to test
 * @param ar2           Test elements (flattened)
 * @param assume_unique If true, both arrays assumed unique
 * @param invert        If true, return ~isin result
 * @param kind          Algorithm selection
 * @return              Boolean array or NULL on error
 */
NDArray* ndarray_isin(NDArray* ar1, NDArray* ar2,
                       bool assume_unique, bool invert,
                       IsinKind kind);

/**
 * Deprecated alias for isin with flattened output.
 */
NDArray* ndarray_in1d(NDArray* ar1, NDArray* ar2,
                       bool assume_unique, bool invert,
                       IsinKind kind);
```

**Add to:** `src/wasm/setops.c`

```c
/**
 * Table-based isin for integer arrays.
 * Creates a lookup table for O(1) membership testing.
 */
static NDArray* isin_table(NDArray* ar1, NDArray* ar2, bool invert) {
    /* Get min/max of ar2 */
    int64_t ar2_min = INT64_MAX;
    int64_t ar2_max = INT64_MIN;

    /* Type-specific min/max computation */
    switch (ar2->dtype) {
        case DTYPE_INT32: {
            int32_t* data = (int32_t*)ar2->data;
            for (size_t i = 0; i < ar2->size; i++) {
                if (data[i] < ar2_min) ar2_min = data[i];
                if (data[i] > ar2_max) ar2_max = data[i];
            }
            break;
        }
        /* ... other integer types ... */
        default:
            return NULL;
    }

    int64_t range = ar2_max - ar2_min;

    /* Check memory constraint */
    if (range > 6 * (int64_t)(ar1->size + ar2->size)) {
        return NULL;  /* Fall back to sort method */
    }

    /* Create lookup table */
    uint8_t* table = calloc((size_t)(range + 1), 1);
    if (!table) return NULL;

    /* Populate table */
    switch (ar2->dtype) {
        case DTYPE_INT32: {
            int32_t* data = (int32_t*)ar2->data;
            for (size_t i = 0; i < ar2->size; i++) {
                table[data[i] - ar2_min] = 1;
            }
            break;
        }
        /* ... other integer types ... */
    }

    /* Create output boolean array (same shape as ar1) */
    NDArray* result = ndarray_empty(ar1->ndim, ar1->shape, DTYPE_BOOL);
    if (!result) {
        free(table);
        return NULL;
    }

    uint8_t* result_data = (uint8_t*)result->data;
    uint8_t fill_val = invert ? 1 : 0;
    memset(result_data, fill_val, ar1->size);

    /* Check membership */
    switch (ar1->dtype) {
        case DTYPE_INT32: {
            int32_t* data = (int32_t*)ar1->data;
            for (size_t i = 0; i < ar1->size; i++) {
                if (data[i] >= ar2_min && data[i] <= ar2_max) {
                    uint8_t in_set = table[data[i] - ar2_min];
                    result_data[i] = invert ? !in_set : in_set;
                }
            }
            break;
        }
        /* ... other integer types ... */
    }

    free(table);
    return result;
}

/**
 * Sort-based isin implementation.
 */
static NDArray* isin_sort(NDArray* ar1, NDArray* ar2,
                          bool assume_unique, bool invert) {
    /* Flatten inputs */
    NDArray* flat1 = ndarray_flatten(ar1);
    NDArray* flat2 = ndarray_flatten(ar2);

    if (!flat1 || !flat2) {
        ndarray_free(flat1);
        ndarray_free(flat2);
        return NULL;
    }

    size_t n1 = flat1->size;
    size_t n2 = flat2->size;

    NDArray* u1 = flat1;
    NDArray* u2 = flat2;
    NDArray* rev_idx = NULL;

    /* Get unique if needed */
    if (!assume_unique) {
        UniqueResult* r1 = ndarray_unique(flat1, false, true, false, true);
        if (!r1) {
            ndarray_free(flat1);
            ndarray_free(flat2);
            return NULL;
        }
        u1 = r1->values; r1->values = NULL;
        rev_idx = r1->inverse; r1->inverse = NULL;
        unique_result_free(r1);

        u2 = ndarray_unique_values(flat2);
        ndarray_free(flat1);
        ndarray_free(flat2);
        flat1 = flat2 = NULL;
    }

    /* Concatenate u1, u2 */
    NDArray* ar = ndarray_concatenate(u1, u2, 0);
    if (!ar) {
        ndarray_free(u1);
        ndarray_free(u2);
        ndarray_free(rev_idx);
        return NULL;
    }

    size_t n = ar->size;

    /* Argsort with mergesort (stable) */
    NDArray* order = ndarray_argsort(ar, 0, SORT_MERGESORT);
    NDArray* sar = ndarray_take(ar, order, 0, CLIP_RAISE);

    /* Create comparison result */
    uint8_t* bool_ar = malloc(n);
    size_t elem_size = dtype_size(sar->dtype);
    char* sar_data = (char*)sar->data;

    for (size_t i = 0; i < n - 1; i++) {
        bool eq = elements_equal(sar_data + (i + 1) * elem_size,
                                  sar_data + i * elem_size,
                                  sar->dtype, true);
        bool_ar[i] = invert ? !eq : eq;
    }
    bool_ar[n - 1] = invert ? 1 : 0;

    /* Scatter back to original order */
    NDArray* ret = ndarray_empty(1, (int32_t[]){(int32_t)n}, DTYPE_BOOL);
    uint8_t* ret_data = (uint8_t*)ret->data;
    int32_t* order_data = (int32_t*)order->data;

    for (size_t i = 0; i < n; i++) {
        ret_data[order_data[i]] = bool_ar[i];
    }

    /* Extract ar1 portion */
    NDArray* result;
    if (assume_unique) {
        /* Take first n1 elements */
        result = ndarray_empty(ar1->ndim, ar1->shape, DTYPE_BOOL);
        memcpy(result->data, ret_data, u1->size);
    } else {
        /* Use reverse index to expand back */
        result = ndarray_empty(ar1->ndim, ar1->shape, DTYPE_BOOL);
        uint8_t* result_data = (uint8_t*)result->data;
        int32_t* rev_data = (int32_t*)rev_idx->data;

        for (size_t i = 0; i < ar1->size; i++) {
            result_data[i] = ret_data[rev_data[i]];
        }
    }

    /* Cleanup */
    free(bool_ar);
    ndarray_free(u1);
    ndarray_free(u2);
    ndarray_free(ar);
    ndarray_free(order);
    ndarray_free(sar);
    ndarray_free(ret);
    ndarray_free(rev_idx);

    return result;
}

EXPORT NDArray* ndarray_isin(NDArray* ar1, NDArray* ar2,
                              bool assume_unique, bool invert,
                              IsinKind kind) {
    if (!ar1 || !ar2) return NULL;

    /* Handle empty ar2 */
    if (ar2->size == 0) {
        NDArray* result = ndarray_empty(ar1->ndim, ar1->shape, DTYPE_BOOL);
        uint8_t fill = invert ? 1 : 0;
        memset(result->data, fill, ar1->size);
        return result;
    }

    /* Check if table method is applicable */
    bool is_int = (ar1->dtype == DTYPE_INT32 || ar1->dtype == DTYPE_INT64 ||
                   ar1->dtype == DTYPE_UINT32 || ar1->dtype == DTYPE_UINT64 ||
                   ar1->dtype == DTYPE_BOOL);

    if (kind == ISIN_TABLE || (kind == ISIN_AUTO && is_int)) {
        NDArray* result = isin_table(ar1, ar2, invert);
        if (result) return result;
        /* Fall through to sort method if table failed */
    }

    return isin_sort(ar1, ar2, assume_unique, invert);
}

EXPORT NDArray* ndarray_in1d(NDArray* ar1, NDArray* ar2,
                              bool assume_unique, bool invert,
                              IsinKind kind) {
    /* in1d is just isin with flattened result */
    NDArray* flat1 = ndarray_flatten(ar1);
    if (!flat1) return NULL;

    NDArray* result = ndarray_isin(flat1, ar2, assume_unique, invert, kind);
    ndarray_free(flat1);

    return result;
}
```

---

### 8.5 TypeScript Wrappers

**File:** `src/ts/setops.ts` (new file)

```typescript
/**
 * NumJS Set Operations
 *
 * TypeScript wrappers for WASM set operations.
 * Adapted from NumPy's lib/_arraysetops_impl.py
 */

import { NDArray } from './NDArray.js';
import { DType } from './types.js';

/* ============ Unique Functions ============ */

/**
 * Result type for unique with multiple return values.
 */
export interface UniqueResult {
  values: NDArray;
  indices?: NDArray;
  inverse?: NDArray;
  counts?: NDArray;
}

/**
 * Find the unique elements of an array.
 *
 * @param arr - Input array (will be flattened)
 * @param options - Optional parameters
 * @param options.returnIndex - If true, return indices of first occurrences
 * @param options.returnInverse - If true, return indices to reconstruct input
 * @param options.returnCounts - If true, return counts of each unique value
 * @param options.equalNan - If true, treat NaN values as equal (default: true)
 * @returns Unique values, or object with values and optional arrays
 *
 * @example
 * ```typescript
 * // Basic usage
 * const arr = await NDArray.fromArray([1, 2, 2, 3, 1]);
 * const uniq = await unique(arr);
 * // [1, 2, 3]
 *
 * // With all returns
 * const { values, indices, inverse, counts } = await unique(arr, {
 *   returnIndex: true,
 *   returnInverse: true,
 *   returnCounts: true
 * });
 * // values: [1, 2, 3]
 * // indices: [0, 1, 3]  (first occurrence of each)
 * // inverse: [0, 1, 1, 2, 0]  (arr = values[inverse])
 * // counts: [2, 2, 1]  (count of each unique)
 * ```
 */
export async function unique(
  arr: NDArray,
  options: {
    returnIndex?: boolean;
    returnInverse?: boolean;
    returnCounts?: boolean;
    equalNan?: boolean;
  } = {}
): Promise<NDArray | UniqueResult> {
  const module = arr._wasmModule;
  const {
    returnIndex = false,
    returnInverse = false,
    returnCounts = false,
    equalNan = true
  } = options;

  const resultPtr = module._ndarray_unique(
    arr._wasmPtr,
    returnIndex,
    returnInverse,
    returnCounts,
    equalNan
  );

  if (resultPtr === 0) {
    throw new Error('unique failed');
  }

  // Read UniqueResult structure from WASM memory
  const valuesPtr = module.HEAP32[resultPtr / 4];
  const indicesPtr = module.HEAP32[resultPtr / 4 + 1];
  const inversePtr = module.HEAP32[resultPtr / 4 + 2];
  const countsPtr = module.HEAP32[resultPtr / 4 + 3];

  const values = NDArray._fromPtr(valuesPtr, module);

  // If no optional returns requested, just return values
  if (!returnIndex && !returnInverse && !returnCounts) {
    module._unique_result_free(resultPtr);
    return values;
  }

  const result: UniqueResult = { values };

  if (returnIndex && indicesPtr !== 0) {
    result.indices = NDArray._fromPtr(indicesPtr, module);
  }
  if (returnInverse && inversePtr !== 0) {
    result.inverse = NDArray._fromPtr(inversePtr, module);
  }
  if (returnCounts && countsPtr !== 0) {
    result.counts = NDArray._fromPtr(countsPtr, module);
  }

  // Don't free the arrays we're returning
  module.HEAP32[resultPtr / 4] = 0;
  module.HEAP32[resultPtr / 4 + 1] = 0;
  module.HEAP32[resultPtr / 4 + 2] = 0;
  module.HEAP32[resultPtr / 4 + 3] = 0;
  module._unique_result_free(resultPtr);

  return result;
}

/**
 * Array API compatible: unique with all return values.
 */
export async function uniqueAll(arr: NDArray): Promise<{
  values: NDArray;
  indices: NDArray;
  inverseIndices: NDArray;
  counts: NDArray;
}> {
  const result = await unique(arr, {
    returnIndex: true,
    returnInverse: true,
    returnCounts: true,
    equalNan: false
  }) as UniqueResult;

  return {
    values: result.values,
    indices: result.indices!,
    inverseIndices: result.inverse!,
    counts: result.counts!
  };
}

/**
 * Array API compatible: unique values and counts.
 */
export async function uniqueCounts(arr: NDArray): Promise<{
  values: NDArray;
  counts: NDArray;
}> {
  const result = await unique(arr, {
    returnCounts: true,
    equalNan: false
  }) as UniqueResult;

  return {
    values: result.values,
    counts: result.counts!
  };
}

/**
 * Array API compatible: unique values and inverse indices.
 */
export async function uniqueInverse(arr: NDArray): Promise<{
  values: NDArray;
  inverseIndices: NDArray;
}> {
  const result = await unique(arr, {
    returnInverse: true,
    equalNan: false
  }) as UniqueResult;

  return {
    values: result.values,
    inverseIndices: result.inverse!
  };
}

/**
 * Array API compatible: just unique values.
 */
export async function uniqueValues(arr: NDArray): Promise<NDArray> {
  return unique(arr, { equalNan: false }) as Promise<NDArray>;
}

/* ============ Set Combination Functions ============ */

/**
 * Find the union of two arrays.
 * Returns unique, sorted values that are in either array.
 *
 * @param ar1 - First input array
 * @param ar2 - Second input array
 * @returns Sorted unique union
 *
 * @example
 * ```typescript
 * const a = await NDArray.fromArray([1, 2, 3]);
 * const b = await NDArray.fromArray([2, 3, 4]);
 * const u = await union1d(a, b);
 * // [1, 2, 3, 4]
 * ```
 */
export async function union1d(ar1: NDArray, ar2: NDArray): Promise<NDArray> {
  const module = ar1._wasmModule;

  const resultPtr = module._ndarray_union1d(ar1._wasmPtr, ar2._wasmPtr);

  if (resultPtr === 0) {
    throw new Error('union1d failed');
  }

  return NDArray._fromPtr(resultPtr, module);
}

/**
 * Find the intersection of two arrays.
 * Returns sorted, unique values that are in both arrays.
 *
 * @param ar1 - First input array
 * @param ar2 - Second input array
 * @param options - Optional parameters
 * @param options.assumeUnique - If true, inputs assumed unique (faster)
 * @param options.returnIndices - If true, also return indices in original arrays
 * @returns Intersection, or object with intersection and indices
 *
 * @example
 * ```typescript
 * const a = await NDArray.fromArray([1, 2, 3, 4]);
 * const b = await NDArray.fromArray([2, 4, 6]);
 * const inter = await intersect1d(a, b);
 * // [2, 4]
 *
 * // With indices
 * const { values, ar1Indices, ar2Indices } = await intersect1d(a, b, {
 *   returnIndices: true
 * });
 * // values: [2, 4]
 * // ar1Indices: [1, 3]  (indices in a)
 * // ar2Indices: [0, 1]  (indices in b)
 * ```
 */
export async function intersect1d(
  ar1: NDArray,
  ar2: NDArray,
  options: {
    assumeUnique?: boolean;
    returnIndices?: boolean;
  } = {}
): Promise<NDArray | { values: NDArray; ar1Indices: NDArray; ar2Indices: NDArray }> {
  const module = ar1._wasmModule;
  const { assumeUnique = false, returnIndices = false } = options;

  // Allocate output pointers for indices
  const outPtr1 = module._malloc(4);
  const outPtr2 = module._malloc(4);

  const resultPtr = module._ndarray_intersect1d(
    ar1._wasmPtr,
    ar2._wasmPtr,
    assumeUnique,
    returnIndices,
    outPtr1,
    outPtr2
  );

  if (resultPtr === 0) {
    module._free(outPtr1);
    module._free(outPtr2);
    throw new Error('intersect1d failed');
  }

  const values = NDArray._fromPtr(resultPtr, module);

  if (!returnIndices) {
    module._free(outPtr1);
    module._free(outPtr2);
    return values;
  }

  const idx1Ptr = module.HEAP32[outPtr1 / 4];
  const idx2Ptr = module.HEAP32[outPtr2 / 4];

  module._free(outPtr1);
  module._free(outPtr2);

  return {
    values,
    ar1Indices: NDArray._fromPtr(idx1Ptr, module),
    ar2Indices: NDArray._fromPtr(idx2Ptr, module)
  };
}

/**
 * Find the set difference (ar1 - ar2).
 * Returns unique values in ar1 that are not in ar2.
 *
 * @param ar1 - First input array
 * @param ar2 - Second input array
 * @param assumeUnique - If true, inputs assumed unique
 * @returns Set difference
 *
 * @example
 * ```typescript
 * const a = await NDArray.fromArray([1, 2, 3, 4]);
 * const b = await NDArray.fromArray([2, 4]);
 * const diff = await setdiff1d(a, b);
 * // [1, 3]
 * ```
 */
export async function setdiff1d(
  ar1: NDArray,
  ar2: NDArray,
  assumeUnique: boolean = false
): Promise<NDArray> {
  const module = ar1._wasmModule;

  const resultPtr = module._ndarray_setdiff1d(
    ar1._wasmPtr,
    ar2._wasmPtr,
    assumeUnique
  );

  if (resultPtr === 0) {
    throw new Error('setdiff1d failed');
  }

  return NDArray._fromPtr(resultPtr, module);
}

/**
 * Find the set exclusive-or (symmetric difference).
 * Returns unique values in exactly one of the arrays.
 *
 * @param ar1 - First input array
 * @param ar2 - Second input array
 * @param assumeUnique - If true, inputs assumed unique
 * @returns Symmetric difference
 *
 * @example
 * ```typescript
 * const a = await NDArray.fromArray([1, 2, 3]);
 * const b = await NDArray.fromArray([2, 3, 4]);
 * const xor = await setxor1d(a, b);
 * // [1, 4]
 * ```
 */
export async function setxor1d(
  ar1: NDArray,
  ar2: NDArray,
  assumeUnique: boolean = false
): Promise<NDArray> {
  const module = ar1._wasmModule;

  const resultPtr = module._ndarray_setxor1d(
    ar1._wasmPtr,
    ar2._wasmPtr,
    assumeUnique
  );

  if (resultPtr === 0) {
    throw new Error('setxor1d failed');
  }

  return NDArray._fromPtr(resultPtr, module);
}

/* ============ Membership Testing ============ */

/**
 * Algorithm selection for isin.
 */
export type IsinKind = 'sort' | 'table' | null;

/**
 * Test whether each element of ar1 is in ar2.
 * Returns boolean array of same shape as ar1.
 *
 * @param element - Input array to test
 * @param testElements - Values to test against (flattened)
 * @param options - Optional parameters
 * @param options.assumeUnique - If true, inputs assumed unique
 * @param options.invert - If true, return ~isin result
 * @param options.kind - Algorithm selection ('sort', 'table', or null for auto)
 * @returns Boolean array
 *
 * @example
 * ```typescript
 * const element = await NDArray.fromArray([[1, 2], [3, 4]]);
 * const test = await NDArray.fromArray([2, 4, 6]);
 * const result = await isin(element, test);
 * // [[false, true], [false, true]]
 *
 * // Inverted
 * const notIn = await isin(element, test, { invert: true });
 * // [[true, false], [true, false]]
 * ```
 */
export async function isin(
  element: NDArray,
  testElements: NDArray,
  options: {
    assumeUnique?: boolean;
    invert?: boolean;
    kind?: IsinKind;
  } = {}
): Promise<NDArray> {
  const module = element._wasmModule;
  const { assumeUnique = false, invert = false, kind = null } = options;

  const kindEnum = kind === 'sort' ? 1 : kind === 'table' ? 2 : 0;

  const resultPtr = module._ndarray_isin(
    element._wasmPtr,
    testElements._wasmPtr,
    assumeUnique,
    invert,
    kindEnum
  );

  if (resultPtr === 0) {
    throw new Error('isin failed');
  }

  return NDArray._fromPtr(resultPtr, module);
}

/**
 * @deprecated Use isin instead.
 * Test whether each element of ar1 is in ar2 (flattened output).
 */
export async function in1d(
  ar1: NDArray,
  ar2: NDArray,
  options: {
    assumeUnique?: boolean;
    invert?: boolean;
    kind?: IsinKind;
  } = {}
): Promise<NDArray> {
  const module = ar1._wasmModule;
  const { assumeUnique = false, invert = false, kind = null } = options;

  const kindEnum = kind === 'sort' ? 1 : kind === 'table' ? 2 : 0;

  const resultPtr = module._ndarray_in1d(
    ar1._wasmPtr,
    ar2._wasmPtr,
    assumeUnique,
    invert,
    kindEnum
  );

  if (resultPtr === 0) {
    throw new Error('in1d failed');
  }

  return NDArray._fromPtr(resultPtr, module);
}

/* ============ Difference Operations ============ */

/**
 * Compute differences between consecutive elements.
 *
 * @param arr - Input array (flattened if not 1D)
 * @param options - Optional parameters
 * @param options.toEnd - Values to append at the end
 * @param options.toBegin - Values to prepend at the beginning
 * @returns Differences array
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([1, 2, 4, 7, 0]);
 * const diff = await ediff1d(arr);
 * // [1, 2, 3, -7]
 *
 * // With prepend/append
 * const diff2 = await ediff1d(arr, {
 *   toBegin: await NDArray.fromArray([-99]),
 *   toEnd: await NDArray.fromArray([88, 99])
 * });
 * // [-99, 1, 2, 3, -7, 88, 99]
 * ```
 */
export async function ediff1d(
  arr: NDArray,
  options: {
    toEnd?: NDArray;
    toBegin?: NDArray;
  } = {}
): Promise<NDArray> {
  const { toEnd, toBegin } = options;

  // Flatten input
  const flat = arr.flatten();
  const n = flat.size;

  if (n === 0) {
    const result: number[] = [];
    if (toBegin) result.push(...toBegin.toArray());
    if (toEnd) result.push(...toEnd.toArray());
    flat.dispose();
    return NDArray.fromArray(result, [result.length], { dtype: arr.dtype });
  }

  // Fast path: no prepend/append
  if (!toBegin && !toEnd) {
    // diff = flat[1:] - flat[:-1]
    const data = flat.toArray();
    const diff = new Array(n - 1);
    for (let i = 0; i < n - 1; i++) {
      diff[i] = data[i + 1] - data[i];
    }
    flat.dispose();
    return NDArray.fromArray(diff, [diff.length], { dtype: arr.dtype });
  }

  // With prepend/append
  const data = flat.toArray();
  const beginData = toBegin ? toBegin.toArray() : [];
  const endData = toEnd ? toEnd.toArray() : [];

  const diffLen = Math.max(n - 1, 0);
  const resultLen = beginData.length + diffLen + endData.length;
  const result = new Array(resultLen);

  // Copy toBegin
  for (let i = 0; i < beginData.length; i++) {
    result[i] = beginData[i];
  }

  // Compute diff
  const offset = beginData.length;
  for (let i = 0; i < diffLen; i++) {
    result[offset + i] = data[i + 1] - data[i];
  }

  // Copy toEnd
  const endOffset = offset + diffLen;
  for (let i = 0; i < endData.length; i++) {
    result[endOffset + i] = endData[i];
  }

  flat.dispose();
  return NDArray.fromArray(result, [resultLen], { dtype: arr.dtype });
}
```

---

## File Changes Summary

### New Files to Create

```
src/wasm/
├── sorting.h          # Sort function declarations
├── sorting.c          # Quicksort, mergesort implementations
├── setops.h           # Set operation declarations
└── setops.c           # unique, union, intersect, setdiff, setxor, isin

src/ts/
└── setops.ts          # TypeScript wrappers for all set operations

tests/ts/
└── phase8.test.ts     # Set operations test suite
```

### Files to Modify

```
src/wasm/ndarray.h
├── Add ndarray_concatenate() if not present
└── Include sorting.h and setops.h

src/ts/index.ts
├── Export all from setops.ts
└── Export sort, argsort from indexing.ts updates

src/ts/indexing.ts
├── Update argsort() to use WASM implementation
└── Add sort() function

scripts/build-wasm.sh
├── Add sorting.c to compilation
├── Add setops.c to compilation
└── Add new EXPORTED_FUNCTIONS
```

---

## Build Script Updates

Add to `scripts/build-wasm.sh` EXPORTED_FUNCTIONS:

```bash
# Sorting
"_ndarray_sort",
"_ndarray_argsort",
"_ndarray_lexsort",

# Set operations
"_ndarray_unique",
"_ndarray_unique_values",
"_unique_result_free",
"_ndarray_union1d",
"_ndarray_intersect1d",
"_ndarray_setdiff1d",
"_ndarray_setxor1d",
"_ndarray_isin",
"_ndarray_in1d"
```

Add new source files:

```bash
emcc \
    "$SRC_DIR/ndarray.c" \
    "$SRC_DIR/dtype.c" \
    "$SRC_DIR/pairwise_sum.c" \
    "$SRC_DIR/broadcast.c" \
    "$SRC_DIR/indexing.c" \
    "$SRC_DIR/sorting.c" \
    "$SRC_DIR/setops.c" \
    -o "$OUT_DIR/numjs.cjs" \
    ...
```

---

## WasmModule Interface Updates

Add to `src/ts/types.ts`:

```typescript
// Sorting
_ndarray_sort(ptr: number, axis: number, kind: number): number;
_ndarray_argsort(ptr: number, axis: number, kind: number): number;
_ndarray_lexsort(keysPtr: number, numKeys: number, axis: number): number;

// Set operations
_ndarray_unique(ptr: number, returnIndex: boolean, returnInverse: boolean,
                returnCounts: boolean, equalNan: boolean): number;
_ndarray_unique_values(ptr: number): number;
_unique_result_free(ptr: number): void;
_ndarray_union1d(ptr1: number, ptr2: number): number;
_ndarray_intersect1d(ptr1: number, ptr2: number, assumeUnique: boolean,
                     returnIndices: boolean, outIdx1Ptr: number, outIdx2Ptr: number): number;
_ndarray_setdiff1d(ptr1: number, ptr2: number, assumeUnique: boolean): number;
_ndarray_setxor1d(ptr1: number, ptr2: number, assumeUnique: boolean): number;
_ndarray_isin(ptr1: number, ptr2: number, assumeUnique: boolean,
              invert: boolean, kind: number): number;
_ndarray_in1d(ptr1: number, ptr2: number, assumeUnique: boolean,
              invert: boolean, kind: number): number;
```

---

## Implementation Order

```
Week 1: Sorting Infrastructure
├── Day 1: C: Comparison functions per dtype
├── Day 2: C: Quicksort implementation
├── Day 3: C: Mergesort implementation (stable)
├── Day 4: C: ndarray_sort(), ndarray_argsort() for 1D
├── Day 5: C: Multi-dimensional axis support + TS wrappers

Week 2: Unique Function
├── Day 1: C: ndarray_unique() core algorithm
├── Day 2: C: return_index, return_inverse support
├── Day 3: C: return_counts, equal_nan support
├── Day 4: TS: unique() wrapper with options
└── Day 5: TS: Array API functions (uniqueAll, uniqueCounts, etc.)

Week 3: Set Combination Operations
├── Day 1: C: ndarray_union1d()
├── Day 2: C: ndarray_intersect1d() with return_indices
├── Day 3: C: ndarray_setdiff1d()
├── Day 4: C: ndarray_setxor1d()
└── Day 5: TS: All wrappers + integration tests

Week 4: Membership Testing & Polish
├── Day 1: C: isin_table() for integers
├── Day 2: C: isin_sort() general implementation
├── Day 3: TS: isin(), in1d(), ediff1d()
├── Day 4: Comprehensive test suite
└── Day 5: Documentation, examples, edge cases
```

---

## Verification Plan

After Phase 8 completion, verify:

```bash
# Build
npm run build

# Run all tests
npm test

# Specific Phase 8 tests should pass:

# Sorting
✓ sort() sorts array in place along axis
✓ argsort() returns indices that would sort
✓ Sorting handles NaN correctly (sorts to end)
✓ Mergesort is stable

# Unique
✓ unique() returns sorted unique elements
✓ unique() with returnIndex returns first occurrence indices
✓ unique() with returnInverse can reconstruct input
✓ unique() with returnCounts returns element counts
✓ unique() handles empty arrays
✓ unique() with equalNan=true collapses NaN values

# Set Combinations
✓ union1d() returns sorted unique union
✓ intersect1d() returns common elements
✓ intersect1d() with returnIndices returns source indices
✓ setdiff1d() returns elements in ar1 not in ar2
✓ setxor1d() returns elements in exactly one array

# Membership Testing
✓ isin() returns boolean array matching ar1 shape
✓ isin() with invert returns negated result
✓ isin() kind='table' uses lookup table for integers
✓ isin() kind='sort' uses sorting for all types
✓ in1d() returns flattened result

# Edge Cases
✓ All functions handle empty arrays
✓ All functions handle single-element arrays
✓ All functions work with all numeric dtypes
```

Generate NumPy comparison vectors:

```python
# tests/python/generate_phase8_tests.py
import numpy as np
import json

tests = {
    "unique": [
        {"data": [1, 2, 2, 3, 1], "expected": [1, 2, 3]},
        {"data": [3, 1, 2, 1, 3], "expected": [1, 2, 3]},
        {"data": [], "expected": []},
        {"data": [5], "expected": [5]},
    ],
    "unique_with_returns": [
        {
            "data": [1, 2, 6, 4, 2, 3, 2],
            "values": [1, 2, 3, 4, 6],
            "indices": [0, 1, 5, 3, 2],
            "inverse": [0, 1, 4, 3, 1, 2, 1],
            "counts": [1, 3, 1, 1, 1]
        }
    ],
    "union1d": [
        {"ar1": [1, 2, 3], "ar2": [2, 3, 4], "expected": [1, 2, 3, 4]},
        {"ar1": [-1, 0, 1], "ar2": [-2, 0, 2], "expected": [-2, -1, 0, 1, 2]},
    ],
    "intersect1d": [
        {"ar1": [1, 3, 4, 3], "ar2": [3, 1, 2, 1], "expected": [1, 3]},
        {"ar1": [1, 2, 3], "ar2": [4, 5, 6], "expected": []},
    ],
    "setdiff1d": [
        {"ar1": [1, 2, 3, 2, 4, 1], "ar2": [3, 4, 5, 6], "expected": [1, 2]},
    ],
    "setxor1d": [
        {"ar1": [1, 2, 3, 2, 4], "ar2": [2, 3, 5, 7, 5], "expected": [1, 4, 5, 7]},
    ],
    "isin": [
        {
            "element": [[0, 2], [4, 6]],
            "test_elements": [1, 2, 4, 8],
            "expected": [[False, True], [True, False]]
        }
    ],
    "ediff1d": [
        {"data": [1, 2, 4, 7, 0], "expected": [1, 2, 3, -7]},
    ]
}

with open("tests/fixtures/phase8_vectors.json", "w") as f:
    json.dump(tests, f, indent=2)
```

---

## Critical Dependencies for Later Phases

Phase 8 completion enables:

- **Phase 6 Statistics**: `unique` used for mode calculation, histograms
- **Phase 10 Functional**: `unique` useful for categorical data processing
- **General Use**: Set operations are fundamental for data manipulation

Phase 8 REQUIRES:

- **Phase 1-2**: Core NDArray, element access, slicing, broadcasting
- **Concatenate function**: Needed for union, intersect, setxor (implement in Phase 5 or as prerequisite)

---

## NumPy Reference Files

| Component | NumPy Reference File |
|-----------|---------------------|
| Set operations | `numpy/lib/_arraysetops_impl.py` |
| Sorting | `numpy/_core/src/npysort/quicksort.cpp` |
| Sorting | `numpy/_core/src/npysort/mergesort.cpp` |
| Argsort | `numpy/_core/src/multiarray/item_selection.c` |
| Set op tests | `numpy/lib/tests/test_arraysetops.py` |
