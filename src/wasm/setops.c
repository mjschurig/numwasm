/**
 * NumJS-WASM Set Operations Implementation
 *
 * Adapted from NumPy's _arraysetops_impl.py.
 * Provides set operations on arrays: unique, union, intersection, etc.
 */

#include "setops.h"
#include "sorting.h"
#include "indexing.h"
#include "manipulation.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#define EXPORT EMSCRIPTEN_KEEPALIVE
#else
#define EXPORT
#endif

/* ============ Helper Functions ============ */

/**
 * Compare two elements for equality.
 * Handles NaN equality based on equal_nan flag.
 */
static bool elements_equal(const void* a, const void* b, DType dtype, bool equal_nan) {
    switch (dtype) {
        case DTYPE_FLOAT64: {
            double va = *(const double*)a;
            double vb = *(const double*)b;
            if (equal_nan && isnan(va) && isnan(vb)) return true;
            return va == vb;
        }
        case DTYPE_FLOAT32: {
            float va = *(const float*)a;
            float vb = *(const float*)b;
            if (equal_nan && isnan(va) && isnan(vb)) return true;
            return va == vb;
        }
        case DTYPE_INT64: {
            return *(const int64_t*)a == *(const int64_t*)b;
        }
        case DTYPE_INT32: {
            return *(const int32_t*)a == *(const int32_t*)b;
        }
        case DTYPE_INT16: {
            return *(const int16_t*)a == *(const int16_t*)b;
        }
        case DTYPE_INT8: {
            return *(const int8_t*)a == *(const int8_t*)b;
        }
        case DTYPE_UINT64: {
            return *(const uint64_t*)a == *(const uint64_t*)b;
        }
        case DTYPE_UINT32: {
            return *(const uint32_t*)a == *(const uint32_t*)b;
        }
        case DTYPE_UINT16: {
            return *(const uint16_t*)a == *(const uint16_t*)b;
        }
        case DTYPE_UINT8: {
            return *(const uint8_t*)a == *(const uint8_t*)b;
        }
        case DTYPE_BOOL: {
            return *(const uint8_t*)a == *(const uint8_t*)b;
        }
        default:
            return memcmp(a, b, dtype_size(dtype)) == 0;
    }
}

/**
 * Get pointer to element at index in array data.
 */
static inline void* get_element_ptr(void* data, size_t index, size_t elem_size) {
    return (uint8_t*)data + index * elem_size;
}

/* ============ UniqueResult Memory Management ============ */

EXPORT void unique_result_free(UniqueResult* result) {
    if (!result) return;
    if (result->values) ndarray_free(result->values);
    if (result->indices) ndarray_free(result->indices);
    if (result->inverse) ndarray_free(result->inverse);
    if (result->counts) ndarray_free(result->counts);
    free(result);
}

/* ============ Unique Implementation ============ */

/**
 * Find unique elements of an array.
 *
 * Algorithm (based on NumPy _unique1d):
 * 1. Flatten input
 * 2. argsort to get permutation (mergesort for stability if return_index)
 * 3. Create sorted array via permutation
 * 4. Build mask where adjacent elements differ
 * 5. Extract unique values, indices, inverse, counts as needed
 */
EXPORT UniqueResult* ndarray_unique(NDArray* arr, bool return_index, bool return_inverse,
                                    bool return_counts, bool equal_nan) {
    if (!arr) return NULL;

    /* Flatten input */
    NDArray* flat = ndarray_flatten(arr);
    if (!flat) return NULL;

    size_t n = flat->size;

    /* Handle empty array */
    if (n == 0) {
        UniqueResult* result = (UniqueResult*)calloc(1, sizeof(UniqueResult));
        if (!result) {
            ndarray_free(flat);
            return NULL;
        }
        int32_t empty_shape[1] = {0};
        result->values = ndarray_empty(1, empty_shape, flat->dtype);
        if (return_index) {
            result->indices = ndarray_empty(1, empty_shape, DTYPE_INT32);
        }
        if (return_inverse) {
            result->inverse = ndarray_empty(1, empty_shape, DTYPE_INT32);
        }
        if (return_counts) {
            result->counts = ndarray_empty(1, empty_shape, DTYPE_INT32);
        }
        ndarray_free(flat);
        return result;
    }

    /* Use mergesort for stability when return_index is needed */
    int32_t sort_kind = return_index ? SORT_MERGESORT : SORT_QUICKSORT;

    /* Get argsort permutation */
    NDArray* perm = ndarray_argsort(flat, 0, sort_kind);
    if (!perm) {
        ndarray_free(flat);
        return NULL;
    }

    /* Create sorted array via take */
    NDArray* aux = ndarray_take(flat, perm, 0, CLIP_RAISE);
    if (!aux) {
        ndarray_free(flat);
        ndarray_free(perm);
        return NULL;
    }

    /* Build mask: mask[i] = (aux[i] != aux[i-1]) */
    /* mask[0] is always true (first element is unique) */
    uint8_t* mask = (uint8_t*)malloc(n);
    if (!mask) {
        ndarray_free(flat);
        ndarray_free(perm);
        ndarray_free(aux);
        return NULL;
    }

    size_t elem_size = dtype_size(aux->dtype);
    mask[0] = 1;
    size_t unique_count = 1;
    for (size_t i = 1; i < n; i++) {
        void* prev = get_element_ptr(aux->data, i - 1, elem_size);
        void* curr = get_element_ptr(aux->data, i, elem_size);
        mask[i] = !elements_equal(prev, curr, aux->dtype, equal_nan);
        if (mask[i]) unique_count++;
    }

    /* Allocate result structure */
    UniqueResult* result = (UniqueResult*)calloc(1, sizeof(UniqueResult));
    if (!result) {
        free(mask);
        ndarray_free(flat);
        ndarray_free(perm);
        ndarray_free(aux);
        return NULL;
    }

    /* Extract unique values */
    int32_t values_shape[1] = {(int32_t)unique_count};
    result->values = ndarray_empty(1, values_shape, aux->dtype);
    if (!result->values) {
        free(mask);
        ndarray_free(flat);
        ndarray_free(perm);
        ndarray_free(aux);
        free(result);
        return NULL;
    }

    size_t out_idx = 0;
    for (size_t i = 0; i < n; i++) {
        if (mask[i]) {
            double val = ndarray_get_flat(aux, i);
            ndarray_set_flat(result->values, out_idx++, val);
        }
    }

    /* return_index: indices of first occurrence in original array */
    if (return_index) {
        result->indices = ndarray_empty(1, values_shape, DTYPE_INT32);
        if (!result->indices) {
            free(mask);
            ndarray_free(flat);
            ndarray_free(perm);
            ndarray_free(aux);
            unique_result_free(result);
            return NULL;
        }
        int32_t* perm_data = (int32_t*)perm->data;
        int32_t* idx_data = (int32_t*)result->indices->data;
        out_idx = 0;
        for (size_t i = 0; i < n; i++) {
            if (mask[i]) {
                idx_data[out_idx++] = perm_data[i];
            }
        }
    }

    /* return_inverse: indices to reconstruct original from unique */
    if (return_inverse) {
        int32_t inverse_shape[1] = {(int32_t)n};
        result->inverse = ndarray_empty(1, inverse_shape, DTYPE_INT32);
        if (!result->inverse) {
            free(mask);
            ndarray_free(flat);
            ndarray_free(perm);
            ndarray_free(aux);
            unique_result_free(result);
            return NULL;
        }

        /* Compute cumulative sum of mask (0-indexed unique ids) */
        int32_t* cumsum = (int32_t*)malloc(n * sizeof(int32_t));
        if (!cumsum) {
            free(mask);
            ndarray_free(flat);
            ndarray_free(perm);
            ndarray_free(aux);
            unique_result_free(result);
            return NULL;
        }
        cumsum[0] = 0;
        for (size_t i = 1; i < n; i++) {
            cumsum[i] = cumsum[i - 1] + (mask[i] ? 1 : 0);
        }

        /* Scatter back to original positions using perm */
        int32_t* perm_data = (int32_t*)perm->data;
        int32_t* inv_data = (int32_t*)result->inverse->data;
        for (size_t i = 0; i < n; i++) {
            inv_data[perm_data[i]] = cumsum[i];
        }
        free(cumsum);
    }

    /* return_counts: count of each unique element */
    if (return_counts) {
        result->counts = ndarray_empty(1, values_shape, DTYPE_INT32);
        if (!result->counts) {
            free(mask);
            ndarray_free(flat);
            ndarray_free(perm);
            ndarray_free(aux);
            unique_result_free(result);
            return NULL;
        }

        int32_t* cnt_data = (int32_t*)result->counts->data;
        out_idx = 0;

        /* Find positions where mask is true */
        size_t* boundaries = (size_t*)malloc((unique_count + 1) * sizeof(size_t));
        if (!boundaries) {
            free(mask);
            ndarray_free(flat);
            ndarray_free(perm);
            ndarray_free(aux);
            unique_result_free(result);
            return NULL;
        }

        size_t b_idx = 0;
        for (size_t i = 0; i < n; i++) {
            if (mask[i]) boundaries[b_idx++] = i;
        }
        boundaries[b_idx] = n;  /* End boundary */

        /* Count = diff of boundaries */
        for (size_t i = 0; i < unique_count; i++) {
            cnt_data[i] = (int32_t)(boundaries[i + 1] - boundaries[i]);
        }
        free(boundaries);
    }

    free(mask);
    ndarray_free(flat);
    ndarray_free(perm);
    ndarray_free(aux);

    return result;
}

/**
 * Return only the unique values (simplified interface).
 */
EXPORT NDArray* ndarray_unique_values(NDArray* arr) {
    UniqueResult* result = ndarray_unique(arr, false, false, false, false);
    if (!result) return NULL;

    NDArray* values = result->values;
    result->values = NULL;  /* Prevent free */
    unique_result_free(result);
    return values;
}

/* ============ Set Combination Operations ============ */

/**
 * Find the union of two arrays.
 *
 * Algorithm: concatenate flattened arrays, then unique.
 */
EXPORT NDArray* ndarray_union1d(NDArray* ar1, NDArray* ar2) {
    if (!ar1 || !ar2) return NULL;

    NDArray* flat1 = ndarray_flatten(ar1);
    if (!flat1) return NULL;

    NDArray* flat2 = ndarray_flatten(ar2);
    if (!flat2) {
        ndarray_free(flat1);
        return NULL;
    }

    NDArray* arrays[2] = {flat1, flat2};
    NDArray* concat = ndarray_concatenate(arrays, 2, 0);
    ndarray_free(flat1);
    ndarray_free(flat2);

    if (!concat) return NULL;

    NDArray* result = ndarray_unique_values(concat);
    ndarray_free(concat);

    return result;
}

/**
 * Find the intersection of two arrays.
 *
 * Algorithm (based on NumPy intersect1d):
 * 1. unique each array (with indices if return_indices)
 * 2. concatenate
 * 3. argsort (mergesort if return_indices for stability)
 * 4. find adjacent equal elements -> intersection
 */
EXPORT NDArray* ndarray_intersect1d(NDArray* ar1, NDArray* ar2, bool assume_unique,
                                    bool return_indices, NDArray** indices1_out,
                                    NDArray** indices2_out) {
    if (!ar1 || !ar2) return NULL;
    if (return_indices && (!indices1_out || !indices2_out)) return NULL;

    NDArray* u1 = NULL;
    NDArray* u2 = NULL;
    NDArray* u1_indices = NULL;
    NDArray* u2_indices = NULL;

    if (assume_unique) {
        u1 = ndarray_flatten(ar1);
        u2 = ndarray_flatten(ar2);
    } else {
        if (return_indices) {
            UniqueResult* res1 = ndarray_unique(ar1, true, false, false, false);
            UniqueResult* res2 = ndarray_unique(ar2, true, false, false, false);
            if (!res1 || !res2) {
                if (res1) unique_result_free(res1);
                if (res2) unique_result_free(res2);
                return NULL;
            }
            u1 = res1->values;
            u1_indices = res1->indices;
            res1->values = NULL;
            res1->indices = NULL;
            unique_result_free(res1);

            u2 = res2->values;
            u2_indices = res2->indices;
            res2->values = NULL;
            res2->indices = NULL;
            unique_result_free(res2);
        } else {
            u1 = ndarray_unique_values(ar1);
            u2 = ndarray_unique_values(ar2);
        }
    }

    if (!u1 || !u2) {
        if (u1) ndarray_free(u1);
        if (u2) ndarray_free(u2);
        if (u1_indices) ndarray_free(u1_indices);
        if (u2_indices) ndarray_free(u2_indices);
        return NULL;
    }

    size_t n1 = u1->size;
    size_t n2 = u2->size;

    /* Concatenate */
    NDArray* arrays[2] = {u1, u2};
    NDArray* aux = ndarray_concatenate(arrays, 2, 0);
    if (!aux) {
        ndarray_free(u1);
        ndarray_free(u2);
        if (u1_indices) ndarray_free(u1_indices);
        if (u2_indices) ndarray_free(u2_indices);
        return NULL;
    }

    /* argsort */
    int32_t sort_kind = return_indices ? SORT_MERGESORT : SORT_QUICKSORT;
    NDArray* order = ndarray_argsort(aux, 0, sort_kind);
    if (!order) {
        ndarray_free(u1);
        ndarray_free(u2);
        ndarray_free(aux);
        if (u1_indices) ndarray_free(u1_indices);
        if (u2_indices) ndarray_free(u2_indices);
        return NULL;
    }

    /* Sort aux by order */
    NDArray* sorted = ndarray_take(aux, order, 0, CLIP_RAISE);
    if (!sorted) {
        ndarray_free(u1);
        ndarray_free(u2);
        ndarray_free(aux);
        ndarray_free(order);
        if (u1_indices) ndarray_free(u1_indices);
        if (u2_indices) ndarray_free(u2_indices);
        return NULL;
    }

    /* Find where adjacent elements are equal */
    size_t total = n1 + n2;
    size_t elem_size = dtype_size(sorted->dtype);
    size_t intersect_count = 0;

    uint8_t* mask = (uint8_t*)calloc(total, 1);
    if (!mask) {
        ndarray_free(u1);
        ndarray_free(u2);
        ndarray_free(aux);
        ndarray_free(order);
        ndarray_free(sorted);
        if (u1_indices) ndarray_free(u1_indices);
        if (u2_indices) ndarray_free(u2_indices);
        return NULL;
    }

    for (size_t i = 0; i < total - 1; i++) {
        void* curr = get_element_ptr(sorted->data, i, elem_size);
        void* next = get_element_ptr(sorted->data, i + 1, elem_size);
        if (elements_equal(curr, next, sorted->dtype, false)) {
            mask[i] = 1;
            intersect_count++;
        }
    }

    /* Extract intersection values */
    int32_t result_shape[1] = {(int32_t)intersect_count};
    NDArray* result = ndarray_empty(1, result_shape, sorted->dtype);
    if (!result) {
        free(mask);
        ndarray_free(u1);
        ndarray_free(u2);
        ndarray_free(aux);
        ndarray_free(order);
        ndarray_free(sorted);
        if (u1_indices) ndarray_free(u1_indices);
        if (u2_indices) ndarray_free(u2_indices);
        return NULL;
    }

    size_t out_idx = 0;
    for (size_t i = 0; i < total; i++) {
        if (mask[i]) {
            double val = ndarray_get_flat(sorted, i);
            ndarray_set_flat(result, out_idx++, val);
        }
    }

    /* Return indices if requested */
    if (return_indices) {
        *indices1_out = ndarray_empty(1, result_shape, DTYPE_INT32);
        *indices2_out = ndarray_empty(1, result_shape, DTYPE_INT32);

        if (!*indices1_out || !*indices2_out) {
            free(mask);
            ndarray_free(u1);
            ndarray_free(u2);
            ndarray_free(aux);
            ndarray_free(order);
            ndarray_free(sorted);
            ndarray_free(result);
            if (u1_indices) ndarray_free(u1_indices);
            if (u2_indices) ndarray_free(u2_indices);
            if (*indices1_out) ndarray_free(*indices1_out);
            if (*indices2_out) ndarray_free(*indices2_out);
            return NULL;
        }

        int32_t* order_data = (int32_t*)order->data;
        int32_t* idx1_data = (int32_t*)(*indices1_out)->data;
        int32_t* idx2_data = (int32_t*)(*indices2_out)->data;
        int32_t* u1_idx_data = u1_indices ? (int32_t*)u1_indices->data : NULL;
        int32_t* u2_idx_data = u2_indices ? (int32_t*)u2_indices->data : NULL;

        out_idx = 0;
        for (size_t i = 0; i < total; i++) {
            if (mask[i]) {
                /* order[i] tells which original position */
                int32_t pos1 = order_data[i];
                int32_t pos2 = order_data[i + 1];

                /* Ensure pos1 is from ar1 and pos2 is from ar2 */
                if (pos1 >= (int32_t)n1) {
                    int32_t tmp = pos1;
                    pos1 = pos2;
                    pos2 = tmp;
                }

                /* Map back to original array indices */
                if (u1_idx_data) {
                    idx1_data[out_idx] = u1_idx_data[pos1];
                } else {
                    idx1_data[out_idx] = pos1;
                }

                if (u2_idx_data) {
                    idx2_data[out_idx] = u2_idx_data[pos2 - (int32_t)n1];
                } else {
                    idx2_data[out_idx] = pos2 - (int32_t)n1;
                }
                out_idx++;
            }
        }
    }

    free(mask);
    ndarray_free(u1);
    ndarray_free(u2);
    ndarray_free(aux);
    ndarray_free(order);
    ndarray_free(sorted);
    if (u1_indices) ndarray_free(u1_indices);
    if (u2_indices) ndarray_free(u2_indices);

    return result;
}

/**
 * Find the set exclusive-or of two arrays.
 *
 * Algorithm (based on NumPy setxor1d):
 * 1. unique each array
 * 2. concatenate + sort
 * 3. flag = [True, aux[1:] != aux[:-1], True]
 * 4. return aux[flag[1:] & flag[:-1]]
 */
EXPORT NDArray* ndarray_setxor1d(NDArray* ar1, NDArray* ar2, bool assume_unique) {
    if (!ar1 || !ar2) return NULL;

    NDArray* u1 = assume_unique ? ndarray_flatten(ar1) : ndarray_unique_values(ar1);
    NDArray* u2 = assume_unique ? ndarray_flatten(ar2) : ndarray_unique_values(ar2);

    if (!u1 || !u2) {
        if (u1) ndarray_free(u1);
        if (u2) ndarray_free(u2);
        return NULL;
    }

    /* Concatenate */
    NDArray* arrays[2] = {u1, u2};
    NDArray* aux = ndarray_concatenate(arrays, 2, 0);
    ndarray_free(u1);
    ndarray_free(u2);

    if (!aux) return NULL;

    size_t n = aux->size;

    /* Handle small cases */
    if (n == 0) return aux;
    if (n == 1) return aux;

    /* Sort in place */
    if (ndarray_sort(aux, 0, SORT_QUICKSORT) != 0) {
        ndarray_free(aux);
        return NULL;
    }

    /* Build flag: flag[i] = (aux[i] != aux[i+1]) for i in [0,n-1]
     * Extended with true at boundaries
     * Result includes aux[i] where flag[i-1] && flag[i] both true
     */
    size_t elem_size = dtype_size(aux->dtype);

    /* flag has n+1 elements: flag[0]=true, flag[i]=diff[i-1] for i in [1,n-1], flag[n]=true */
    uint8_t* flag = (uint8_t*)malloc(n + 1);
    if (!flag) {
        ndarray_free(aux);
        return NULL;
    }

    flag[0] = 1;
    flag[n] = 1;
    for (size_t i = 0; i < n - 1; i++) {
        void* curr = get_element_ptr(aux->data, i, elem_size);
        void* next = get_element_ptr(aux->data, i + 1, elem_size);
        flag[i + 1] = !elements_equal(curr, next, aux->dtype, false);
    }

    /* Count elements where flag[i] && flag[i+1] */
    size_t count = 0;
    for (size_t i = 0; i < n; i++) {
        if (flag[i] && flag[i + 1]) count++;
    }

    /* Extract result */
    int32_t result_shape[1] = {(int32_t)count};
    NDArray* result = ndarray_empty(1, result_shape, aux->dtype);
    if (!result) {
        free(flag);
        ndarray_free(aux);
        return NULL;
    }

    size_t out_idx = 0;
    for (size_t i = 0; i < n; i++) {
        if (flag[i] && flag[i + 1]) {
            double val = ndarray_get_flat(aux, i);
            ndarray_set_flat(result, out_idx++, val);
        }
    }

    free(flag);
    ndarray_free(aux);

    return result;
}

/* ============ Membership Testing ============ */

/**
 * Check if ar2 values fit in a lookup table.
 * Returns true if table method is feasible.
 */
static bool can_use_table_method(NDArray* ar1, NDArray* ar2, int64_t* min_out, int64_t* max_out) {
    /* Only for integer types */
    if (!dtype_is_integer(ar2->dtype) || ar2->dtype == DTYPE_INT64 || ar2->dtype == DTYPE_UINT64) {
        return false;
    }

    if (ar2->size == 0) return false;

    /* Find min/max of ar2 */
    int64_t min_val = INT64_MAX;
    int64_t max_val = INT64_MIN;

    for (size_t i = 0; i < ar2->size; i++) {
        int64_t val = (int64_t)ndarray_get_flat(ar2, i);
        if (val < min_val) min_val = val;
        if (val > max_val) max_val = val;
    }

    /* Check if range is reasonable (6 * (n1 + n2) memory limit from NumPy) */
    int64_t range = max_val - min_val + 1;
    if (range < 0 || range > (int64_t)(6 * (ar1->size + ar2->size))) {
        return false;
    }

    *min_out = min_val;
    *max_out = max_val;
    return true;
}

/**
 * isin using table lookup (O(n1 + n2) for bounded integers).
 */
static NDArray* isin_table(NDArray* ar1, NDArray* ar2, bool invert, int64_t min_val, int64_t max_val) {
    int64_t range = max_val - min_val + 1;

    /* Build lookup table */
    uint8_t* table = (uint8_t*)calloc((size_t)range, 1);
    if (!table) return NULL;

    for (size_t i = 0; i < ar2->size; i++) {
        int64_t val = (int64_t)ndarray_get_flat(ar2, i);
        table[val - min_val] = 1;
    }

    /* Create boolean result with same shape as ar1 */
    NDArray* result = ndarray_empty(ar1->ndim, ar1->shape, DTYPE_BOOL);
    if (!result) {
        free(table);
        return NULL;
    }

    uint8_t* res_data = (uint8_t*)result->data;
    for (size_t i = 0; i < ar1->size; i++) {
        int64_t val = (int64_t)ndarray_get_flat(ar1, i);
        bool found;
        if (val < min_val || val > max_val) {
            found = false;
        } else {
            found = table[val - min_val] != 0;
        }
        res_data[i] = invert ? !found : found;
    }

    free(table);
    return result;
}

/**
 * isin using sorting (general case).
 *
 * Algorithm:
 * 1. Flatten ar1, get unique with inverse
 * 2. Get unique ar2
 * 3. Concatenate uniques
 * 4. argsort (mergesort)
 * 5. Find where adjacent elements equal
 * 6. Map back via inverse
 */
static NDArray* isin_sort(NDArray* ar1, NDArray* ar2, bool invert, bool assume_unique) {
    NDArray* flat1 = ndarray_flatten(ar1);
    if (!flat1) return NULL;

    size_t n1 = flat1->size;

    /* Handle empty cases */
    if (n1 == 0) {
        ndarray_free(flat1);
        return ndarray_empty(ar1->ndim, ar1->shape, DTYPE_BOOL);
    }

    NDArray* flat2 = ndarray_flatten(ar2);
    if (!flat2) {
        ndarray_free(flat1);
        return NULL;
    }

    /* Get unique ar1 with inverse to map back */
    UniqueResult* res1;
    if (assume_unique) {
        /* Create dummy result */
        res1 = (UniqueResult*)calloc(1, sizeof(UniqueResult));
        if (!res1) {
            ndarray_free(flat1);
            ndarray_free(flat2);
            return NULL;
        }
        res1->values = flat1;
        flat1 = NULL;  /* Prevent double free */

        /* Create identity inverse */
        int32_t inv_shape[1] = {(int32_t)n1};
        res1->inverse = ndarray_empty(1, inv_shape, DTYPE_INT32);
        if (!res1->inverse) {
            unique_result_free(res1);
            ndarray_free(flat2);
            return NULL;
        }
        int32_t* inv_data = (int32_t*)res1->inverse->data;
        for (size_t i = 0; i < n1; i++) {
            inv_data[i] = (int32_t)i;
        }
    } else {
        res1 = ndarray_unique(ar1, false, true, false, false);
        if (!res1) {
            if (flat1) ndarray_free(flat1);
            ndarray_free(flat2);
            return NULL;
        }
        if (flat1) ndarray_free(flat1);
    }

    NDArray* u1 = res1->values;
    NDArray* rev_idx = res1->inverse;
    size_t nu1 = u1->size;

    /* Get unique ar2 */
    NDArray* u2 = assume_unique ? flat2 : ndarray_unique_values(flat2);
    if (!assume_unique) ndarray_free(flat2);
    if (!u2) {
        unique_result_free(res1);
        return NULL;
    }

    size_t nu2 = u2->size;

    /* Handle empty ar2 */
    if (nu2 == 0) {
        ndarray_free(u2);
        NDArray* result = ndarray_create(ar1->ndim, ar1->shape, DTYPE_BOOL);
        if (result && invert) {
            /* Fill with true */
            ndarray_fill(result, 1.0);
        }
        unique_result_free(res1);
        return result;
    }

    /* Concatenate */
    NDArray* arrays[2] = {u1, u2};
    NDArray* ar = ndarray_concatenate(arrays, 2, 0);
    ndarray_free(u2);
    if (!ar) {
        unique_result_free(res1);
        return NULL;
    }

    /* argsort with mergesort */
    NDArray* order = ndarray_argsort(ar, 0, SORT_MERGESORT);
    if (!order) {
        ndarray_free(ar);
        unique_result_free(res1);
        return NULL;
    }

    /* Sort ar by order */
    NDArray* sar = ndarray_take(ar, order, 0, CLIP_RAISE);
    ndarray_free(ar);
    if (!sar) {
        ndarray_free(order);
        unique_result_free(res1);
        return NULL;
    }

    /* Find where adjacent equal and first came from ar1 */
    size_t total = nu1 + nu2;
    size_t elem_size = dtype_size(sar->dtype);
    int32_t* order_data = (int32_t*)order->data;

    /* bool_ar: which positions in u1 have a match in u2 */
    uint8_t* bool_ar = (uint8_t*)calloc(nu1, 1);
    if (!bool_ar) {
        ndarray_free(sar);
        ndarray_free(order);
        unique_result_free(res1);
        return NULL;
    }

    for (size_t i = 0; i < total - 1; i++) {
        void* curr = get_element_ptr(sar->data, i, elem_size);
        void* next = get_element_ptr(sar->data, i + 1, elem_size);
        if (elements_equal(curr, next, sar->dtype, false)) {
            /* Both order[i] and order[i+1] point to same value */
            /* Mark if either is from u1 (index < nu1) */
            if (order_data[i] < (int32_t)nu1) {
                bool_ar[order_data[i]] = 1;
            }
            if (order_data[i + 1] < (int32_t)nu1) {
                bool_ar[order_data[i + 1]] = 1;
            }
        }
    }

    ndarray_free(sar);
    ndarray_free(order);

    /* Map back to original shape using reverse index */
    NDArray* result = ndarray_empty(ar1->ndim, ar1->shape, DTYPE_BOOL);
    if (!result) {
        free(bool_ar);
        unique_result_free(res1);
        return NULL;
    }

    uint8_t* res_data = (uint8_t*)result->data;
    int32_t* rev_data = (int32_t*)rev_idx->data;
    for (size_t i = 0; i < n1; i++) {
        bool found = bool_ar[rev_data[i]] != 0;
        res_data[i] = invert ? !found : found;
    }

    free(bool_ar);
    unique_result_free(res1);

    return result;
}

/**
 * Test whether each element of ar1 is also present in ar2.
 */
EXPORT NDArray* ndarray_isin(NDArray* ar1, NDArray* ar2, bool invert,
                             bool assume_unique, int32_t kind) {
    if (!ar1 || !ar2) return NULL;

    /* Empty ar2 means nothing is in it */
    if (ar2->size == 0) {
        NDArray* result = ndarray_create(ar1->ndim, ar1->shape, DTYPE_BOOL);
        if (result && invert) {
            ndarray_fill(result, 1.0);
        }
        return result;
    }

    /* Choose algorithm */
    int64_t min_val, max_val;
    bool use_table = (kind != ISIN_SORT) && can_use_table_method(ar1, ar2, &min_val, &max_val);

    if (kind == ISIN_TABLE && !use_table) {
        /* User requested table but it's not feasible */
        return NULL;
    }

    if (use_table || kind == ISIN_TABLE) {
        return isin_table(ar1, ar2, invert, min_val, max_val);
    }

    return isin_sort(ar1, ar2, invert, assume_unique);
}

/**
 * Test whether each element of ar1 is also present in ar2 (flat version).
 */
EXPORT NDArray* ndarray_in1d(NDArray* ar1, NDArray* ar2, bool invert,
                             bool assume_unique, int32_t kind) {
    NDArray* result = ndarray_isin(ar1, ar2, invert, assume_unique, kind);
    if (!result) return NULL;

    /* Flatten result if not already 1D */
    if (result->ndim != 1) {
        NDArray* flat = ndarray_flatten(result);
        ndarray_free(result);
        return flat;
    }
    return result;
}

/* ============ Set Difference ============ */

/**
 * Find the set difference of two arrays.
 *
 * Uses isin with invert=true.
 */
EXPORT NDArray* ndarray_setdiff1d(NDArray* ar1, NDArray* ar2, bool assume_unique) {
    if (!ar1 || !ar2) return NULL;

    NDArray* u1 = assume_unique ? ndarray_flatten(ar1) : ndarray_unique_values(ar1);
    if (!u1) return NULL;

    /* If ar2 is empty, return u1 */
    if (ar2->size == 0) return u1;

    /* Get mask: elements of u1 NOT in ar2 */
    NDArray* mask = ndarray_isin(u1, ar2, true, assume_unique, ISIN_AUTO);
    if (!mask) {
        ndarray_free(u1);
        return NULL;
    }

    /* Compress u1 by mask */
    NDArray* result = ndarray_compress(mask, u1, 0);
    ndarray_free(u1);
    ndarray_free(mask);

    return result;
}
