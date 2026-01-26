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
    (void)kth;  /* Suppress unused warning */
    return ndarray_argsort(arr, axis, SORT_QUICKSORT);
}
