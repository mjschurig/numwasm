# Level 1: Core Operations Implementation Plan

Level 1 extends the NumJS-WASM foundation with three major components: Element Access (get/set by index), Iteration Infrastructure, and Shape Manipulation. This layer enables direct array element manipulation and view-based transformations.

---

## Current State (Level 0 Complete)

```
src/wasm/
├── ndarray.h          # NDArray struct, DType enum (4 types)
├── ndarray.c          # create, from_data, free, sum, fill, accessors
├── pairwise_sum.h     # Pairwise summation declarations
└── pairwise_sum.c     # Pairwise sum for float32/float64

src/ts/
├── types.ts           # DType enum, WasmModule interface
├── NDArray.ts         # zeros, ones, fromArray, arange, sum, fill, toArray
├── wasm-loader.ts     # WASM module loading
└── index.ts           # Public exports
```

**Existing DTypes:** Float32, Float64, Int32, Int64
**Existing Operations:** zeros, ones, fromArray, arange, sum, fill, toArray, dispose

---

## Level 1 Implementation Tree

```
LEVEL 1: CORE OPERATIONS
│
├── 1.1 View Infrastructure
│   │
│   ├── 1.1.1 C: Update NDArray Struct
│   │   ├── Add: struct NDArray* base (NULL if owns data)
│   │   ├── Add: NDARRAY_C_CONTIGUOUS = 0x0004
│   │   └── Add: NDARRAY_F_CONTIGUOUS = 0x0008
│   │
│   ├── 1.1.2 C: View Creation
│   │   ├── ndarray_view(src, ndim, shape, strides) → NDArray*
│   │   │   └── Creates view sharing data with src
│   │   └── ndarray_view_with_offset(src, ndim, shape, strides, offset) → NDArray*
│   │       └── View with byte offset into data
│   │
│   ├── 1.1.3 C: Contiguity Checks
│   │   ├── ndarray_is_c_contiguous(arr) → int
│   │   ├── ndarray_is_f_contiguous(arr) → int
│   │   └── ndarray_update_contiguity_flags(arr)
│   │
│   ├── 1.1.4 C: Property Accessors
│   │   ├── ndarray_get_strides(arr) → int32_t*
│   │   ├── ndarray_get_flags(arr) → int32_t
│   │   └── ndarray_get_base(arr) → NDArray*
│   │
│   ├── 1.1.5 C: Update ndarray_free
│   │   └── Only free data if base == NULL
│   │
│   └── 1.1.6 TypeScript: View Tracking
│       ├── _base: NDArray | null (reference to base)
│       ├── isView property
│       └── _wrapView(ptr) helper method
│
├── 1.2 Element Access
│   │
│   ├── 1.2.1 C: Byte Offset Calculation
│   │   └── calc_offset(arr, indices) → size_t
│   │       └── Returns: sum(indices[i] * strides[i])
│   │
│   ├── 1.2.2 C: Bounds Checking
│   │   └── ndarray_check_bounds(arr, indices, ndim) → int
│   │       └── Validates all indices within shape
│   │
│   ├── 1.2.3 C: Element Getter
│   │   └── ndarray_get_item(arr, indices, ndim) → double
│   │       ├── Calculates offset using strides
│   │       ├── Reads value at (data + offset)
│   │       └── Converts to double based on dtype
│   │
│   ├── 1.2.4 C: Element Setter
│   │   └── ndarray_set_item(arr, indices, ndim, value)
│   │       ├── Checks NDARRAY_WRITEABLE flag
│   │       ├── Calculates offset using strides
│   │       └── Writes value converted to dtype
│   │
│   ├── 1.2.5 TypeScript: get() Method
│   │   └── get(...indices: number[]) → number
│   │       ├── Validates index count matches ndim
│   │       ├── Allocates indices in WASM memory
│   │       ├── Calls _ndarray_get_item
│   │       └── Frees WASM memory
│   │
│   ├── 1.2.6 TypeScript: set() Method
│   │   └── set(value: number, ...indices: number[])
│   │       ├── Validates index count matches ndim
│   │       ├── Allocates indices in WASM memory
│   │       ├── Calls _ndarray_set_item
│   │       └── Frees WASM memory
│   │
│   ├── 1.2.7 TypeScript: item() Method
│   │   └── item() → number
│   │       ├── Validates size === 1
│   │       └── Returns single element value
│   │
│   └── 1.2.8 TypeScript: strides Property
│       └── get strides() → number[]
│           └── Reads strides array from WASM
│
├── 1.3 Shape Manipulation
│   │
│   ├── 1.3.1 C: Reshape
│   │   └── ndarray_reshape(arr, new_ndim, new_shape) → NDArray*
│   │       ├── Handle -1 dimension (auto-calculate)
│   │       ├── Verify total size unchanged
│   │       ├── Check if contiguous (view possible)
│   │       ├── Compute new C-contiguous strides
│   │       └── Return view or NULL if copy needed
│   │
│   ├── 1.3.2 C: Transpose
│   │   └── ndarray_transpose(arr, axes) → NDArray*
│   │       ├── Default: reverse all axes
│   │       ├── Custom: permute by axes array
│   │       ├── Permute shape and strides
│   │       └── Always returns view
│   │
│   ├── 1.3.3 C: Ravel (flatten to 1D)
│   │   └── ndarray_ravel(arr) → NDArray*
│   │       ├── If contiguous: return view with shape=[size]
│   │       └── If non-contiguous: return NULL (copy needed)
│   │
│   ├── 1.3.4 C: Flatten (always copy)
│   │   └── ndarray_flatten(arr) → NDArray*
│   │       ├── Allocate new contiguous array
│   │       └── Copy elements in row-major order
│   │
│   ├── 1.3.5 C: Squeeze
│   │   └── ndarray_squeeze(arr, axis) → NDArray*
│   │       ├── axis=-1: remove all size-1 dimensions
│   │       ├── axis>=0: remove specific axis if size=1
│   │       └── Returns view with reduced shape
│   │
│   ├── 1.3.6 C: Expand Dims
│   │   └── ndarray_expand_dims(arr, axis) → NDArray*
│   │       ├── Insert size-1 dimension at axis
│   │       ├── Handle negative axis (from end)
│   │       └── Returns view with expanded shape
│   │
│   ├── 1.3.7 C: Swap Axes
│   │   └── ndarray_swapaxes(arr, axis1, axis2) → NDArray*
│   │       ├── Swap shape[axis1] with shape[axis2]
│   │       ├── Swap strides[axis1] with strides[axis2]
│   │       └── Returns view
│   │
│   ├── 1.3.8 C: Copy
│   │   └── ndarray_copy(arr) → NDArray*
│   │       ├── Allocate new array with same shape
│   │       ├── Copy all elements
│   │       └── Always returns owned data
│   │
│   ├── 1.3.9 TypeScript: reshape()
│   │   └── reshape(newShape: number[]) → NDArray
│   │
│   ├── 1.3.10 TypeScript: T Property
│   │   └── get T() → NDArray (transpose)
│   │
│   ├── 1.3.11 TypeScript: transpose()
│   │   └── transpose(axes?: number[]) → NDArray
│   │
│   ├── 1.3.12 TypeScript: ravel()
│   │   └── ravel() → NDArray
│   │
│   ├── 1.3.13 TypeScript: flatten()
│   │   └── flatten() → NDArray
│   │
│   ├── 1.3.14 TypeScript: squeeze()
│   │   └── squeeze(axis?: number) → NDArray
│   │
│   ├── 1.3.15 TypeScript: expandDims()
│   │   └── expandDims(axis: number) → NDArray
│   │
│   ├── 1.3.16 TypeScript: swapaxes()
│   │   └── swapaxes(axis1: number, axis2: number) → NDArray
│   │
│   └── 1.3.17 TypeScript: copy()
│       └── copy() → NDArray
│
└── 1.4 Iteration Infrastructure
    │
    ├── 1.4.1 TypeScript: FlatIterator Class
    │   └── class FlatIterator implements IterableIterator<number>
    │       ├── constructor(arr: NDArray)
    │       ├── [Symbol.iterator]() → this
    │       ├── next() → IteratorResult<number>
    │       └── _unravelIndex(flat) → number[]
    │
    ├── 1.4.2 TypeScript: flat Property
    │   └── NDArray.flat → FlatIterator
    │       └── Iterates elements in row-major order
    │
    ├── 1.4.3 TypeScript: nditer Generator
    │   └── function* nditer(arr) → Generator<number>
    │       └── Yields each element value
    │
    ├── 1.4.4 TypeScript: ndenumerate Generator
    │   └── function* ndenumerate(arr) → Generator<[number[], number]>
    │       └── Yields [indices, value] pairs
    │
    └── 1.4.5 TypeScript: ndindex Generator
        └── function* ndindex(...shape) → Generator<number[]>
            └── Yields all index combinations for shape
```

---

## Detailed Implementation Specifications

### 1.1.1 Update NDArray Struct (C)

```c
// In ndarray.h

// Additional flags
#define NDARRAY_C_CONTIGUOUS  0x0004
#define NDARRAY_F_CONTIGUOUS  0x0008

// Updated struct
typedef struct NDArray {
    void* data;              // Pointer to raw data buffer
    int32_t ndim;            // Number of dimensions
    int32_t* shape;          // Size in each dimension
    int32_t* strides;        // Bytes to jump in each dimension
    DType dtype;             // Data type
    int32_t flags;           // Memory ownership and layout flags
    size_t size;             // Total number of elements
    struct NDArray* base;    // NEW: Base array for views (NULL if owns data)
} NDArray;
```

### 1.1.2 View Creation (C)

```c
// In ndarray.c
EXPORT NDArray* ndarray_view(NDArray* src, int32_t ndim,
                              int32_t* shape, int32_t* strides)
{
    return ndarray_view_with_offset(src, ndim, shape, strides, 0);
}

EXPORT NDArray* ndarray_view_with_offset(NDArray* src, int32_t ndim,
                                          int32_t* shape, int32_t* strides,
                                          size_t byte_offset)
{
    if (!src) return NULL;

    NDArray* view = (NDArray*)malloc(sizeof(NDArray));
    if (!view) return NULL;

    view->ndim = ndim;
    view->dtype = src->dtype;
    view->flags = NDARRAY_WRITEABLE;  // No OWNDATA flag
    view->size = compute_size(ndim, shape);
    view->data = (char*)src->data + byte_offset;

    // Find ultimate base (for chained views)
    view->base = (src->base != NULL) ? src->base : src;

    // Allocate and copy shape/strides
    view->shape = (int32_t*)malloc(ndim * sizeof(int32_t));
    view->strides = (int32_t*)malloc(ndim * sizeof(int32_t));
    if (!view->shape || !view->strides) {
        free(view->shape);
        free(view->strides);
        free(view);
        return NULL;
    }

    memcpy(view->shape, shape, ndim * sizeof(int32_t));
    memcpy(view->strides, strides, ndim * sizeof(int32_t));

    // Update contiguity flags
    ndarray_update_contiguity_flags(view);

    return view;
}
```

### 1.1.3 Contiguity Checks (C)

```c
// In ndarray.c
EXPORT int ndarray_is_c_contiguous(NDArray* arr)
{
    if (!arr || arr->ndim == 0) return 1;

    size_t expected = dtype_size(arr->dtype);
    for (int i = arr->ndim - 1; i >= 0; i--) {
        if (arr->shape[i] == 0) return 1;  // Empty array
        if (arr->shape[i] != 1 && (size_t)arr->strides[i] != expected) {
            return 0;
        }
        expected *= arr->shape[i];
    }
    return 1;
}

EXPORT int ndarray_is_f_contiguous(NDArray* arr)
{
    if (!arr || arr->ndim == 0) return 1;

    size_t expected = dtype_size(arr->dtype);
    for (int i = 0; i < arr->ndim; i++) {
        if (arr->shape[i] == 0) return 1;
        if (arr->shape[i] != 1 && (size_t)arr->strides[i] != expected) {
            return 0;
        }
        expected *= arr->shape[i];
    }
    return 1;
}

EXPORT void ndarray_update_contiguity_flags(NDArray* arr)
{
    if (!arr) return;

    arr->flags &= ~(NDARRAY_C_CONTIGUOUS | NDARRAY_F_CONTIGUOUS);

    if (ndarray_is_c_contiguous(arr)) {
        arr->flags |= NDARRAY_C_CONTIGUOUS;
    }
    if (ndarray_is_f_contiguous(arr)) {
        arr->flags |= NDARRAY_F_CONTIGUOUS;
    }
}
```

### 1.1.5 Update ndarray_free (C)

```c
// In ndarray.c
EXPORT void ndarray_free(NDArray* arr)
{
    if (!arr) return;

    // Only free data if we own it (base == NULL and OWNDATA flag set)
    if (arr->base == NULL && (arr->flags & NDARRAY_OWNDATA)) {
        free(arr->data);
    }

    free(arr->shape);
    free(arr->strides);
    free(arr);
}
```

### 1.2.1-1.2.4 Element Access (C)

```c
// In ndarray.c

// Calculate byte offset from N-dimensional indices
static size_t calc_offset(NDArray* arr, int32_t* indices)
{
    size_t offset = 0;
    for (int i = 0; i < arr->ndim; i++) {
        offset += (size_t)indices[i] * (size_t)arr->strides[i];
    }
    return offset;
}

EXPORT int ndarray_check_bounds(NDArray* arr, int32_t* indices, int32_t ndim)
{
    if (!arr || ndim != arr->ndim) return 0;

    for (int i = 0; i < ndim; i++) {
        if (indices[i] < 0 || indices[i] >= arr->shape[i]) {
            return 0;
        }
    }
    return 1;
}

EXPORT double ndarray_get_item(NDArray* arr, int32_t* indices, int32_t ndim)
{
    if (!arr || !arr->data || ndim != arr->ndim) return 0.0;
    if (!ndarray_check_bounds(arr, indices, ndim)) return 0.0;

    size_t byte_offset = calc_offset(arr, indices);
    void* ptr = (char*)arr->data + byte_offset;

    switch (arr->dtype) {
        case DTYPE_FLOAT64: return *((double*)ptr);
        case DTYPE_FLOAT32: return (double)*((float*)ptr);
        case DTYPE_INT32:   return (double)*((int32_t*)ptr);
        case DTYPE_INT64:   return (double)*((int64_t*)ptr);
        default: return 0.0;
    }
}

EXPORT void ndarray_set_item(NDArray* arr, int32_t* indices,
                              int32_t ndim, double value)
{
    if (!arr || !arr->data || ndim != arr->ndim) return;
    if (!(arr->flags & NDARRAY_WRITEABLE)) return;
    if (!ndarray_check_bounds(arr, indices, ndim)) return;

    size_t byte_offset = calc_offset(arr, indices);
    void* ptr = (char*)arr->data + byte_offset;

    switch (arr->dtype) {
        case DTYPE_FLOAT64: *((double*)ptr) = value; break;
        case DTYPE_FLOAT32: *((float*)ptr) = (float)value; break;
        case DTYPE_INT32:   *((int32_t*)ptr) = (int32_t)value; break;
        case DTYPE_INT64:   *((int64_t*)ptr) = (int64_t)value; break;
    }
}
```

### 1.2.5-1.2.8 Element Access (TypeScript)

```typescript
// In NDArray.ts

/**
 * Get element value at specified indices.
 * @param indices - N-dimensional indices
 * @returns Element value as number
 * @example
 * arr.get(0)       // 1D array
 * arr.get(1, 2)    // 2D array, row 1, col 2
 * arr.get(0, 1, 2) // 3D array
 */
get(...indices: number[]): number {
    this.ensureNotDisposed();

    if (indices.length !== this.ndim) {
        throw new Error(`Expected ${this.ndim} indices, got ${indices.length}`);
    }

    // Bounds checking
    const shape = this.shape;
    for (let i = 0; i < indices.length; i++) {
        if (indices[i] < 0 || indices[i] >= shape[i]) {
            throw new RangeError(
                `Index ${indices[i]} out of bounds for axis ${i} with size ${shape[i]}`
            );
        }
    }

    // Allocate indices array in WASM memory
    const indicesPtr = this._module._malloc(indices.length * 4);
    for (let i = 0; i < indices.length; i++) {
        this._module.setValue(indicesPtr + i * 4, indices[i], 'i32');
    }

    const value = this._module._ndarray_get_item(
        this._ptr, indicesPtr, indices.length
    );
    this._module._free(indicesPtr);

    return value;
}

/**
 * Set element value at specified indices.
 * @param value - Value to set
 * @param indices - N-dimensional indices
 */
set(value: number, ...indices: number[]): void {
    this.ensureNotDisposed();

    if (indices.length !== this.ndim) {
        throw new Error(`Expected ${this.ndim} indices, got ${indices.length}`);
    }

    // Bounds checking
    const shape = this.shape;
    for (let i = 0; i < indices.length; i++) {
        if (indices[i] < 0 || indices[i] >= shape[i]) {
            throw new RangeError(
                `Index ${indices[i]} out of bounds for axis ${i} with size ${shape[i]}`
            );
        }
    }

    const indicesPtr = this._module._malloc(indices.length * 4);
    for (let i = 0; i < indices.length; i++) {
        this._module.setValue(indicesPtr + i * 4, indices[i], 'i32');
    }

    this._module._ndarray_set_item(
        this._ptr, indicesPtr, indices.length, value
    );
    this._module._free(indicesPtr);
}

/**
 * Return the array as a scalar (for 0-d or single-element arrays).
 * @returns Single element value
 */
item(): number {
    this.ensureNotDisposed();

    if (this.size !== 1) {
        throw new Error(
            `item() only works for single-element arrays, got size ${this.size}`
        );
    }

    const indices = new Array(this.ndim).fill(0);
    return this.get(...indices);
}

/**
 * Array strides in bytes.
 */
get strides(): number[] {
    this.ensureNotDisposed();

    const stridesPtr = this._module._ndarray_get_strides(this._ptr);
    const strides: number[] = [];
    for (let i = 0; i < this.ndim; i++) {
        strides.push(this._module.getValue(stridesPtr + i * 4, 'i32'));
    }
    return strides;
}

/**
 * Whether this array is a view (shares data with another array).
 */
get isView(): boolean {
    this.ensureNotDisposed();
    return this._module._ndarray_get_base(this._ptr) !== 0;
}
```

### 1.3.1 Reshape (C)

```c
// In ndarray.c
EXPORT NDArray* ndarray_reshape(NDArray* arr, int32_t new_ndim, int32_t* new_shape)
{
    if (!arr) return NULL;

    // Resolve -1 dimension
    int32_t resolved[32];
    int unknown_idx = -1;
    size_t known_product = 1;

    for (int i = 0; i < new_ndim; i++) {
        if (new_shape[i] == -1) {
            if (unknown_idx != -1) return NULL;  // Multiple -1 not allowed
            unknown_idx = i;
        } else if (new_shape[i] < 0) {
            return NULL;  // Invalid negative dimension
        } else {
            resolved[i] = new_shape[i];
            known_product *= new_shape[i];
        }
    }

    if (unknown_idx >= 0) {
        if (known_product == 0 || arr->size % known_product != 0) {
            return NULL;  // Cannot determine unknown dimension
        }
        resolved[unknown_idx] = arr->size / known_product;
    }

    // Verify total size matches
    size_t new_size = 1;
    for (int i = 0; i < new_ndim; i++) {
        new_size *= resolved[i];
    }
    if (new_size != arr->size) return NULL;

    // Only contiguous arrays can be reshaped as views
    if (!ndarray_is_c_contiguous(arr)) return NULL;

    // Compute new C-contiguous strides
    int32_t new_strides[32];
    compute_strides(new_ndim, resolved, arr->dtype, new_strides);

    return ndarray_view(arr, new_ndim, resolved, new_strides);
}
```

### 1.3.2 Transpose (C)

```c
// In ndarray.c
EXPORT NDArray* ndarray_transpose(NDArray* arr, int32_t* axes)
{
    if (!arr) return NULL;

    int32_t new_shape[32];
    int32_t new_strides[32];
    int32_t default_axes[32];

    // Default: reverse all axes
    if (axes == NULL) {
        for (int i = 0; i < arr->ndim; i++) {
            default_axes[i] = arr->ndim - 1 - i;
        }
        axes = default_axes;
    }

    // Validate axes (check for valid range and no duplicates)
    uint32_t seen = 0;
    for (int i = 0; i < arr->ndim; i++) {
        int ax = axes[i];
        if (ax < 0) ax += arr->ndim;
        if (ax < 0 || ax >= arr->ndim) return NULL;
        if (seen & (1u << ax)) return NULL;  // Duplicate axis
        seen |= (1u << ax);
    }

    // Permute shape and strides
    for (int i = 0; i < arr->ndim; i++) {
        int ax = axes[i];
        if (ax < 0) ax += arr->ndim;
        new_shape[i] = arr->shape[ax];
        new_strides[i] = arr->strides[ax];
    }

    return ndarray_view(arr, arr->ndim, new_shape, new_strides);
}
```

### 1.3.3-1.3.8 Other Shape Operations (C)

```c
// In ndarray.c

EXPORT NDArray* ndarray_ravel(NDArray* arr)
{
    if (!arr) return NULL;

    // If contiguous, can return view
    if (ndarray_is_c_contiguous(arr)) {
        int32_t shape[1] = { (int32_t)arr->size };
        int32_t strides[1] = { (int32_t)dtype_size(arr->dtype) };
        return ndarray_view(arr, 1, shape, strides);
    }

    // Non-contiguous: return NULL to signal copy needed
    return NULL;
}

EXPORT NDArray* ndarray_flatten(NDArray* arr)
{
    if (!arr) return NULL;

    // Always create a copy
    int32_t shape[1] = { (int32_t)arr->size };
    NDArray* result = ndarray_create(1, shape, arr->dtype);
    if (!result) return NULL;

    // Copy elements in row-major order
    // Use iterator pattern for non-contiguous arrays
    size_t elem_size = dtype_size(arr->dtype);

    if (ndarray_is_c_contiguous(arr)) {
        // Fast path: memcpy
        memcpy(result->data, arr->data, arr->size * elem_size);
    } else {
        // Slow path: element by element
        int32_t indices[32] = {0};
        for (size_t flat = 0; flat < arr->size; flat++) {
            double val = ndarray_get_item(arr, indices, arr->ndim);

            // Set in result (which is contiguous)
            void* dst = (char*)result->data + flat * elem_size;
            switch (arr->dtype) {
                case DTYPE_FLOAT64: *((double*)dst) = val; break;
                case DTYPE_FLOAT32: *((float*)dst) = (float)val; break;
                case DTYPE_INT32:   *((int32_t*)dst) = (int32_t)val; break;
                case DTYPE_INT64:   *((int64_t*)dst) = (int64_t)val; break;
            }

            // Increment indices (row-major order)
            for (int d = arr->ndim - 1; d >= 0; d--) {
                indices[d]++;
                if (indices[d] < arr->shape[d]) break;
                indices[d] = 0;
            }
        }
    }

    return result;
}

EXPORT NDArray* ndarray_squeeze(NDArray* arr, int32_t axis)
{
    if (!arr) return NULL;

    int32_t new_shape[32];
    int32_t new_strides[32];
    int new_ndim = 0;

    if (axis == -1) {
        // Remove all size-1 dimensions
        for (int i = 0; i < arr->ndim; i++) {
            if (arr->shape[i] != 1) {
                new_shape[new_ndim] = arr->shape[i];
                new_strides[new_ndim] = arr->strides[i];
                new_ndim++;
            }
        }
    } else {
        // Remove specific axis
        if (axis < 0) axis += arr->ndim;
        if (axis < 0 || axis >= arr->ndim) return NULL;
        if (arr->shape[axis] != 1) return NULL;  // Can only squeeze size-1

        for (int i = 0; i < arr->ndim; i++) {
            if (i != axis) {
                new_shape[new_ndim] = arr->shape[i];
                new_strides[new_ndim] = arr->strides[i];
                new_ndim++;
            }
        }
    }

    // Handle case where all dimensions are squeezed (scalar)
    if (new_ndim == 0) {
        new_ndim = 1;
        new_shape[0] = 1;
        new_strides[0] = dtype_size(arr->dtype);
    }

    return ndarray_view(arr, new_ndim, new_shape, new_strides);
}

EXPORT NDArray* ndarray_expand_dims(NDArray* arr, int32_t axis)
{
    if (!arr) return NULL;

    // Handle negative axis
    if (axis < 0) axis += arr->ndim + 1;
    if (axis < 0 || axis > arr->ndim) return NULL;

    int32_t new_shape[33];
    int32_t new_strides[33];
    int new_ndim = arr->ndim + 1;

    int j = 0;
    for (int i = 0; i < new_ndim; i++) {
        if (i == axis) {
            new_shape[i] = 1;
            // Stride for size-1 dimension doesn't matter for indexing
            // Use next stride or element size
            if (j < arr->ndim) {
                new_strides[i] = arr->strides[j];
            } else {
                new_strides[i] = dtype_size(arr->dtype);
            }
        } else {
            new_shape[i] = arr->shape[j];
            new_strides[i] = arr->strides[j];
            j++;
        }
    }

    return ndarray_view(arr, new_ndim, new_shape, new_strides);
}

EXPORT NDArray* ndarray_swapaxes(NDArray* arr, int32_t axis1, int32_t axis2)
{
    if (!arr) return NULL;

    // Handle negative axes
    if (axis1 < 0) axis1 += arr->ndim;
    if (axis2 < 0) axis2 += arr->ndim;
    if (axis1 < 0 || axis1 >= arr->ndim) return NULL;
    if (axis2 < 0 || axis2 >= arr->ndim) return NULL;

    int32_t new_shape[32];
    int32_t new_strides[32];

    memcpy(new_shape, arr->shape, arr->ndim * sizeof(int32_t));
    memcpy(new_strides, arr->strides, arr->ndim * sizeof(int32_t));

    // Swap
    int32_t tmp = new_shape[axis1];
    new_shape[axis1] = new_shape[axis2];
    new_shape[axis2] = tmp;

    tmp = new_strides[axis1];
    new_strides[axis1] = new_strides[axis2];
    new_strides[axis2] = tmp;

    return ndarray_view(arr, arr->ndim, new_shape, new_strides);
}

EXPORT NDArray* ndarray_copy(NDArray* arr)
{
    if (!arr) return NULL;

    NDArray* result = ndarray_create(arr->ndim, arr->shape, arr->dtype);
    if (!result) return NULL;

    if (ndarray_is_c_contiguous(arr)) {
        // Fast path
        memcpy(result->data, arr->data, arr->size * dtype_size(arr->dtype));
    } else {
        // Slow path: use flatten logic
        int32_t indices[32] = {0};
        size_t elem_size = dtype_size(arr->dtype);

        for (size_t flat = 0; flat < arr->size; flat++) {
            double val = ndarray_get_item(arr, indices, arr->ndim);
            void* dst = (char*)result->data + flat * elem_size;

            switch (arr->dtype) {
                case DTYPE_FLOAT64: *((double*)dst) = val; break;
                case DTYPE_FLOAT32: *((float*)dst) = (float)val; break;
                case DTYPE_INT32:   *((int32_t*)dst) = (int32_t)val; break;
                case DTYPE_INT64:   *((int64_t*)dst) = (int64_t)val; break;
            }

            for (int d = arr->ndim - 1; d >= 0; d--) {
                indices[d]++;
                if (indices[d] < arr->shape[d]) break;
                indices[d] = 0;
            }
        }
    }

    return result;
}
```

### 1.3.9-1.3.17 Shape Manipulation (TypeScript)

```typescript
// In NDArray.ts

/**
 * Private helper to wrap a WASM pointer as a view.
 */
private _wrapView(ptr: number): NDArray {
    return new NDArray(ptr, this._module);
}

/**
 * Returns array reshaped to new shape.
 * Returns a view when possible (contiguous arrays).
 * @param newShape - New shape (can include one -1 for auto-calculation)
 */
reshape(newShape: number[]): NDArray {
    this.ensureNotDisposed();

    const shapePtr = this._module._malloc(newShape.length * 4);
    for (let i = 0; i < newShape.length; i++) {
        this._module.setValue(shapePtr + i * 4, newShape[i], 'i32');
    }

    const resultPtr = this._module._ndarray_reshape(
        this._ptr, newShape.length, shapePtr
    );
    this._module._free(shapePtr);

    if (resultPtr === 0) {
        throw new Error(
            'Cannot reshape: incompatible shape or non-contiguous array'
        );
    }

    return this._wrapView(resultPtr);
}

/**
 * Transpose (reverse axes).
 */
get T(): NDArray {
    this.ensureNotDisposed();

    const resultPtr = this._module._ndarray_transpose(this._ptr, 0);
    if (resultPtr === 0) {
        throw new Error('Transpose failed');
    }

    return this._wrapView(resultPtr);
}

/**
 * Transpose with custom axes permutation.
 * @param axes - Optional axes permutation
 */
transpose(axes?: number[]): NDArray {
    this.ensureNotDisposed();

    let axesPtr = 0;
    if (axes) {
        if (axes.length !== this.ndim) {
            throw new Error(
                `Axes must have ${this.ndim} elements, got ${axes.length}`
            );
        }
        axesPtr = this._module._malloc(axes.length * 4);
        for (let i = 0; i < axes.length; i++) {
            this._module.setValue(axesPtr + i * 4, axes[i], 'i32');
        }
    }

    const resultPtr = this._module._ndarray_transpose(this._ptr, axesPtr);

    if (axesPtr) this._module._free(axesPtr);

    if (resultPtr === 0) {
        throw new Error('Transpose failed: invalid axes');
    }

    return this._wrapView(resultPtr);
}

/**
 * Return flattened array as view (if contiguous) or copy.
 */
ravel(): NDArray {
    this.ensureNotDisposed();

    const resultPtr = this._module._ndarray_ravel(this._ptr);
    if (resultPtr === 0) {
        // Non-contiguous: fall back to flatten (copy)
        return this.flatten();
    }

    return this._wrapView(resultPtr);
}

/**
 * Return flattened array (always a copy).
 */
flatten(): NDArray {
    this.ensureNotDisposed();

    const resultPtr = this._module._ndarray_flatten(this._ptr);
    if (resultPtr === 0) {
        throw new Error('Flatten failed');
    }

    return this._wrapView(resultPtr);
}

/**
 * Remove size-1 dimensions.
 * @param axis - Specific axis to squeeze, or undefined for all
 */
squeeze(axis?: number): NDArray {
    this.ensureNotDisposed();

    const ax = axis ?? -1;
    const resultPtr = this._module._ndarray_squeeze(this._ptr, ax);
    if (resultPtr === 0) {
        throw new Error(
            axis !== undefined
                ? `Cannot squeeze axis ${axis}: size is not 1`
                : 'Squeeze failed'
        );
    }

    return this._wrapView(resultPtr);
}

/**
 * Add a size-1 dimension at the specified position.
 * @param axis - Position for new axis (supports negative indexing)
 */
expandDims(axis: number): NDArray {
    this.ensureNotDisposed();

    const resultPtr = this._module._ndarray_expand_dims(this._ptr, axis);
    if (resultPtr === 0) {
        throw new Error(`Invalid axis ${axis} for array with ${this.ndim} dimensions`);
    }

    return this._wrapView(resultPtr);
}

/**
 * Swap two axes.
 */
swapaxes(axis1: number, axis2: number): NDArray {
    this.ensureNotDisposed();

    const resultPtr = this._module._ndarray_swapaxes(this._ptr, axis1, axis2);
    if (resultPtr === 0) {
        throw new Error(`Invalid axes: ${axis1}, ${axis2}`);
    }

    return this._wrapView(resultPtr);
}

/**
 * Return a deep copy of the array.
 */
copy(): NDArray {
    this.ensureNotDisposed();

    const resultPtr = this._module._ndarray_copy(this._ptr);
    if (resultPtr === 0) {
        throw new Error('Copy failed');
    }

    return this._wrapView(resultPtr);
}
```

### 1.4.1-1.4.5 Iteration Infrastructure (TypeScript)

```typescript
// In src/ts/iterators.ts (new file)

import type { NDArray } from './NDArray';

/**
 * Iterator over array elements in row-major (C) order.
 * Implements JavaScript Iterator protocol.
 */
export class FlatIterator implements IterableIterator<number> {
    private _arr: NDArray;
    private _index: number = 0;
    private _size: number;
    private _shape: number[];

    constructor(arr: NDArray) {
        this._arr = arr;
        this._size = arr.size;
        this._shape = arr.shape;
    }

    [Symbol.iterator](): IterableIterator<number> {
        return this;
    }

    next(): IteratorResult<number> {
        if (this._index >= this._size) {
            return { done: true, value: undefined };
        }

        const indices = this._unravelIndex(this._index);
        const value = this._arr.get(...indices);
        this._index++;

        return { done: false, value };
    }

    /**
     * Convert flat index to N-dimensional indices (row-major order).
     */
    private _unravelIndex(flatIndex: number): number[] {
        const indices: number[] = new Array(this._shape.length);
        let remaining = flatIndex;

        for (let i = this._shape.length - 1; i >= 0; i--) {
            indices[i] = remaining % this._shape[i];
            remaining = Math.floor(remaining / this._shape[i]);
        }

        return indices;
    }

    /**
     * Current flat index.
     */
    get index(): number {
        return this._index;
    }

    /**
     * Reset iterator to beginning.
     */
    reset(): void {
        this._index = 0;
    }
}

/**
 * Iterate over array elements.
 * @param arr - Array to iterate
 * @yields Element values in row-major order
 */
export function* nditer(arr: NDArray): Generator<number> {
    for (const value of arr.flat) {
        yield value;
    }
}

/**
 * Iterate over array with indices and values.
 * @param arr - Array to iterate
 * @yields [indices, value] pairs
 * @example
 * for (const [idx, val] of ndenumerate(arr)) {
 *     console.log(`arr[${idx}] = ${val}`);
 * }
 */
export function* ndenumerate(arr: NDArray): Generator<[number[], number]> {
    const shape = arr.shape;
    const ndim = shape.length;
    const indices = new Array(ndim).fill(0);

    for (let flat = 0; flat < arr.size; flat++) {
        yield [indices.slice(), arr.get(...indices)];

        // Increment indices (row-major order)
        for (let d = ndim - 1; d >= 0; d--) {
            indices[d]++;
            if (indices[d] < shape[d]) break;
            indices[d] = 0;
        }
    }
}

/**
 * Generate all index combinations for a shape.
 * @param shape - Dimensions
 * @yields Index arrays
 * @example
 * for (const idx of ndindex(2, 3)) {
 *     console.log(idx);  // [0,0], [0,1], [0,2], [1,0], [1,1], [1,2]
 * }
 */
export function* ndindex(...shape: number[]): Generator<number[]> {
    const ndim = shape.length;
    if (ndim === 0) return;

    const size = shape.reduce((a, b) => a * b, 1);
    const indices = new Array(ndim).fill(0);

    for (let flat = 0; flat < size; flat++) {
        yield indices.slice();

        for (let d = ndim - 1; d >= 0; d--) {
            indices[d]++;
            if (indices[d] < shape[d]) break;
            indices[d] = 0;
        }
    }
}
```

```typescript
// Add to NDArray.ts
import { FlatIterator } from './iterators';

// Inside NDArray class:

/**
 * Flat iterator over array elements in row-major order.
 */
get flat(): FlatIterator {
    this.ensureNotDisposed();
    return new FlatIterator(this);
}
```

---

## File Changes Summary

### New Files to Create

```
src/ts/
└── iterators.ts       # FlatIterator class, nditer, ndenumerate, ndindex
```

### Files to Modify

```
src/wasm/ndarray.h
├── Add base pointer to NDArray struct
├── Add NDARRAY_C_CONTIGUOUS, NDARRAY_F_CONTIGUOUS flags
├── Declare view creation functions
├── Declare element access functions
├── Declare shape manipulation functions
└── Declare copy function

src/wasm/ndarray.c
├── Implement ndarray_view, ndarray_view_with_offset
├── Implement contiguity checks
├── Update ndarray_free for views
├── Implement element access (get_item, set_item, check_bounds)
├── Implement reshape, transpose, ravel, flatten
├── Implement squeeze, expand_dims, swapaxes
└── Implement copy

src/ts/types.ts
├── Add new WasmModule function declarations
└── (No new types needed)

src/ts/NDArray.ts
├── Import FlatIterator
├── Add get(), set(), item() methods
├── Add strides, isView properties
├── Add reshape(), T, transpose() methods
├── Add ravel(), flatten() methods
├── Add squeeze(), expandDims(), swapaxes() methods
├── Add copy() method
├── Add flat property
└── Add private _wrapView() helper

src/ts/index.ts
└── Export FlatIterator, nditer, ndenumerate, ndindex

scripts/build-wasm.sh
└── Add new EXPORTED_FUNCTIONS
```

---

## Build Script Updates

Add to `scripts/build-wasm.sh` EXPORTED_FUNCTIONS:

```bash
"_ndarray_view",
"_ndarray_view_with_offset",
"_ndarray_is_c_contiguous",
"_ndarray_is_f_contiguous",
"_ndarray_get_strides",
"_ndarray_get_flags",
"_ndarray_get_base",
"_ndarray_get_item",
"_ndarray_set_item",
"_ndarray_check_bounds",
"_ndarray_reshape",
"_ndarray_transpose",
"_ndarray_ravel",
"_ndarray_flatten",
"_ndarray_squeeze",
"_ndarray_expand_dims",
"_ndarray_swapaxes",
"_ndarray_copy"
```

---

## WasmModule Interface Updates

```typescript
// Add to src/ts/types.ts WasmModule interface

// View system
_ndarray_view(ptr: number, ndim: number, shapePtr: number, stridesPtr: number): number;
_ndarray_view_with_offset(ptr: number, ndim: number, shapePtr: number,
                           stridesPtr: number, byteOffset: number): number;
_ndarray_is_c_contiguous(ptr: number): number;
_ndarray_is_f_contiguous(ptr: number): number;
_ndarray_get_strides(ptr: number): number;
_ndarray_get_flags(ptr: number): number;
_ndarray_get_base(ptr: number): number;

// Element access
_ndarray_get_item(ptr: number, indicesPtr: number, ndim: number): number;
_ndarray_set_item(ptr: number, indicesPtr: number, ndim: number, value: number): void;
_ndarray_check_bounds(ptr: number, indicesPtr: number, ndim: number): number;

// Shape manipulation
_ndarray_reshape(ptr: number, ndim: number, shapePtr: number): number;
_ndarray_transpose(ptr: number, axesPtr: number): number;
_ndarray_ravel(ptr: number): number;
_ndarray_flatten(ptr: number): number;
_ndarray_squeeze(ptr: number, axis: number): number;
_ndarray_expand_dims(ptr: number, axis: number): number;
_ndarray_swapaxes(ptr: number, axis1: number, axis2: number): number;
_ndarray_copy(ptr: number): number;
```

---

## Implementation Order (Dependency-Driven)

```
Week 1: View Infrastructure
├── Day 1: Update NDArray struct with base pointer
├── Day 2: Implement ndarray_view functions
├── Day 3: Implement contiguity checks
├── Day 4: Update ndarray_free, add property accessors
└── Day 5: Test view creation and memory handling

Week 2: Element Access
├── Day 1: Implement calc_offset, check_bounds
├── Day 2: Implement ndarray_get_item, ndarray_set_item
├── Day 3: TypeScript get(), set() methods
├── Day 4: TypeScript item(), strides property
└── Day 5: Tests for element access

Week 3: Shape Manipulation
├── Day 1: Implement reshape (C)
├── Day 2: Implement transpose (C)
├── Day 3: Implement ravel, flatten (C)
├── Day 4: Implement squeeze, expand_dims, swapaxes (C)
├── Day 5: Implement copy (C)

Week 4: TypeScript + Iterators
├── Day 1: TypeScript reshape, T, transpose
├── Day 2: TypeScript ravel, flatten, copy
├── Day 3: TypeScript squeeze, expandDims, swapaxes
├── Day 4: FlatIterator class
├── Day 5: nditer, ndenumerate, ndindex generators

Week 5: Integration & Testing
├── Day 1-2: Comprehensive unit tests
├── Day 3: NumPy comparison test vectors
├── Day 4: Edge cases and error handling
└── Day 5: Documentation and examples
```

---

## Verification Plan

After Level 1 completion, verify:

```bash
# Build
npm run build

# Run all tests
npm test

# Specific Level 1 tests should pass:
✓ get() returns correct element for 1D, 2D, 3D arrays
✓ set() modifies correct element
✓ Bounds checking throws RangeError
✓ Wrong index count throws Error
✓ reshape() creates view for contiguous arrays
✓ reshape() handles -1 dimension
✓ reshape() fails on incompatible size
✓ T property returns transposed view
✓ transpose() with custom axes works
✓ Views share data (modifications visible)
✓ ravel() returns view for contiguous
✓ flatten() always returns copy
✓ squeeze() removes size-1 dimensions
✓ expandDims() adds dimension at position
✓ swapaxes() swaps two axes
✓ copy() creates independent copy
✓ FlatIterator traverses row-major order
✓ ndenumerate yields [indices, value] pairs
✓ ndindex generates correct index sequences
```

Generate NumPy comparison vectors:

```python
# tests/python/generate_level1_tests.py
import numpy as np
import json

tests = {
    "element_access": [
        {
            "shape": [2, 3],
            "data": [1, 2, 3, 4, 5, 6],
            "tests": [
                {"indices": [0, 0], "expected": 1},
                {"indices": [0, 2], "expected": 3},
                {"indices": [1, 0], "expected": 4},
                {"indices": [1, 2], "expected": 6},
            ]
        },
    ],
    "reshape": [
        {"input_shape": [6], "new_shape": [2, 3], "expected_shape": [2, 3]},
        {"input_shape": [6], "new_shape": [3, 2], "expected_shape": [3, 2]},
        {"input_shape": [2, 3], "new_shape": [-1], "expected_shape": [6]},
        {"input_shape": [2, 3], "new_shape": [3, -1], "expected_shape": [3, 2]},
    ],
    "transpose": [
        {"shape": [2, 3], "expected_shape": [3, 2]},
        {"shape": [2, 3, 4], "expected_shape": [4, 3, 2]},
        {"shape": [2, 3, 4], "axes": [2, 0, 1], "expected_shape": [4, 2, 3]},
    ],
    "iteration_order": [
        {
            "shape": [2, 3],
            "data": [1, 2, 3, 4, 5, 6],
            "expected_flat_order": [1, 2, 3, 4, 5, 6]
        },
    ],
}

with open("tests/fixtures/level1_vectors.json", "w") as f:
    json.dump(tests, f, indent=2)
```

---

## Critical Dependencies for Level 2

Level 1 completion enables:

- **Level 2.1 Basic Slicing**: Requires views (1.1), element access (1.2)
- **Level 2.2 Broadcasting**: Requires views (1.1), shape ops (1.3)
- **Level 2.3 Advanced Indexing**: Requires element access (1.2)

Level 1 MUST be complete before proceeding to Level 2.
