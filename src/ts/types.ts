/**
 * NumJS Type Definitions
 *
 * TypeScript types for the WASM-based NDArray implementation.
 */

/**
 * Data type enum - matches the C DType enum in dtype.h
 */
export enum DType {
  // Original types (maintain backward compatibility)
  Float32 = 0,
  Float64 = 1,
  Int32 = 2,
  Int64 = 3,

  // Boolean
  Bool = 4,

  // Additional integers
  Int8 = 5,
  Int16 = 6,
  Uint8 = 7,
  Uint16 = 8,
  Uint32 = 9,
  Uint64 = 10,

  // Half precision (limited support)
  Float16 = 11,

  // Complex types
  Complex64 = 12,
  Complex128 = 13,
}

/**
 * Maps DType to byte size
 */
export const DTYPE_SIZES: Record<DType, number> = {
  [DType.Float32]: 4,
  [DType.Float64]: 8,
  [DType.Int32]: 4,
  [DType.Int64]: 8,
  [DType.Bool]: 1,
  [DType.Int8]: 1,
  [DType.Int16]: 2,
  [DType.Uint8]: 1,
  [DType.Uint16]: 2,
  [DType.Uint32]: 4,
  [DType.Uint64]: 8,
  [DType.Float16]: 2,
  [DType.Complex64]: 8,
  [DType.Complex128]: 16,
};

/**
 * Maps DType to string name
 */
export const DTYPE_NAMES: Record<DType, string> = {
  [DType.Float32]: 'float32',
  [DType.Float64]: 'float64',
  [DType.Int32]: 'int32',
  [DType.Int64]: 'int64',
  [DType.Bool]: 'bool',
  [DType.Int8]: 'int8',
  [DType.Int16]: 'int16',
  [DType.Uint8]: 'uint8',
  [DType.Uint16]: 'uint16',
  [DType.Uint32]: 'uint32',
  [DType.Uint64]: 'uint64',
  [DType.Float16]: 'float16',
  [DType.Complex64]: 'complex64',
  [DType.Complex128]: 'complex128',
};

/**
 * Complex number representation
 */
export interface Complex {
  real: number;
  imag: number;
}

/**
 * Options for NDArray creation
 */
export interface NDArrayOptions {
  dtype?: DType;
}

/**
 * Supported TypedArray types for array creation
 */
export type TypedArrayType =
  | Float32Array
  | Float64Array
  | Int32Array
  | Int16Array
  | Int8Array
  | Uint32Array
  | Uint16Array
  | Uint8Array
  | BigInt64Array
  | BigUint64Array;

/**
 * Array flags interface
 */
export interface NDArrayFlags {
  owndata: boolean;
  writeable: boolean;
  c_contiguous: boolean;
  f_contiguous: boolean;
  aligned: boolean;
}

/**
 * Flag constants (must match C defines)
 */
export const NDARRAY_OWNDATA = 0x0001;
export const NDARRAY_WRITEABLE = 0x0002;
export const NDARRAY_C_CONTIGUOUS = 0x0004;
export const NDARRAY_F_CONTIGUOUS = 0x0008;
export const NDARRAY_ALIGNED = 0x0010;

/**
 * Clip mode constants for index functions (must match C defines)
 */
export const CLIP_RAISE = 0; // Raise error on out-of-bounds
export const CLIP_WRAP = 1; // Wrap around (modulo)
export const CLIP_CLIP = 2; // Clip to valid range

/**
 * Index specification types (must match C defines)
 */
export const INDEX_TYPE_INTEGER = 0;
export const INDEX_TYPE_SLICE = 1;
export const INDEX_TYPE_NEWAXIS = 2;
export const INDEX_TYPE_ELLIPSIS = 3;

/**
 * WASM Module interface - defines the functions exported from the C code
 */
export interface WasmModule {
  // NDArray creation and destruction
  _ndarray_create(ndim: number, shapePtr: number, dtype: number): number;
  _ndarray_from_data(
    dataPtr: number,
    ndim: number,
    shapePtr: number,
    dtype: number
  ): number;
  _ndarray_empty(ndim: number, shapePtr: number, dtype: number): number;
  _ndarray_full(
    ndim: number,
    shapePtr: number,
    dtype: number,
    value: number
  ): number;
  _ndarray_scalar(value: number, dtype: number): number;
  _ndarray_free(ptr: number): void;
  _ndarray_copy(ptr: number): number;
  _ndarray_astype(ptr: number, dtype: number): number;

  // Element access
  _ndarray_flat_index(ptr: number, indicesPtr: number, ndim: number): number;
  _ndarray_check_bounds(
    ptr: number,
    indicesPtr: number,
    ndim: number
  ): boolean;
  _ndarray_get_item(ptr: number, indicesPtr: number, ndim: number): number;
  _ndarray_set_item(
    ptr: number,
    indicesPtr: number,
    ndim: number,
    value: number
  ): void;
  _ndarray_get_flat(ptr: number, flatIdx: number): number;
  _ndarray_set_flat(ptr: number, flatIdx: number, value: number): void;
  _ndarray_get_complex_real(ptr: number, flatIdx: number): number;
  _ndarray_get_complex_imag(ptr: number, flatIdx: number): number;
  _ndarray_set_complex(
    ptr: number,
    flatIdx: number,
    real: number,
    imag: number
  ): void;

  // NDArray operations
  _ndarray_sum(ptr: number): number;
  _ndarray_fill(ptr: number, value: number): void;

  // NDArray property accessors
  _ndarray_get_ndim(ptr: number): number;
  _ndarray_get_shape(ptr: number): number;
  _ndarray_get_strides(ptr: number): number;
  _ndarray_get_data(ptr: number): number;
  _ndarray_get_size(ptr: number): number;
  _ndarray_get_dtype(ptr: number): number;
  _ndarray_get_flags(ptr: number): number;
  _ndarray_get_base(ptr: number): number;

  // Contiguity checks
  _ndarray_is_c_contiguous(ptr: number): boolean;
  _ndarray_is_f_contiguous(ptr: number): boolean;

  // Views
  _ndarray_view(
    ptr: number,
    ndim: number,
    shapePtr: number,
    stridesPtr: number
  ): number;
  _ndarray_view_with_offset(
    ptr: number,
    ndim: number,
    shapePtr: number,
    stridesPtr: number,
    byteOffset: number
  ): number;

  // Shape manipulation
  _ndarray_reshape(ptr: number, newNdim: number, newShapePtr: number): number;
  _ndarray_transpose(ptr: number, axesPtr: number): number;
  _ndarray_ravel(ptr: number): number;
  _ndarray_flatten(ptr: number): number;
  _ndarray_squeeze(ptr: number, axis: number): number;
  _ndarray_expand_dims(ptr: number, axis: number): number;
  _ndarray_swapaxes(ptr: number, axis1: number, axis2: number): number;

  // Slicing (Level 2)
  _ndarray_slice(ptr: number, indicesPtr: number, numIndices: number): number;
  _ndarray_get_subarray(ptr: number, index: number): number;

  // View extensions (Level 2)
  _ndarray_view_dtype(ptr: number, dtype: number): number;
  _ndarray_ascontiguousarray(ptr: number): number;
  _ndarray_asfortranarray(ptr: number): number;

  // Broadcasting (Level 2)
  _broadcast_shapes(
    shape1Ptr: number,
    ndim1: number,
    shape2Ptr: number,
    ndim2: number,
    outShapePtr: number,
    outNdimPtr: number
  ): number;
  _broadcast_shapes_multi(
    shapesPtr: number,
    ndimsPtr: number,
    numArrays: number,
    outShapePtr: number,
    outNdimPtr: number
  ): number;
  _broadcast_strides(
    ptr: number,
    targetShapePtr: number,
    targetNdim: number,
    outStridesPtr: number
  ): number;
  _ndarray_broadcast_to(
    ptr: number,
    targetShapePtr: number,
    targetNdim: number
  ): number;
  _shapes_are_broadcastable(
    shape1Ptr: number,
    ndim1: number,
    shape2Ptr: number,
    ndim2: number
  ): number;

  // Index functions (Level 2)
  _ndarray_take(
    ptr: number,
    indicesPtr: number,
    axis: number,
    clipmode: number
  ): number;
  _ndarray_take_flat(ptr: number, indicesPtr: number, clipmode: number): number;
  _ndarray_put(
    ptr: number,
    indicesPtr: number,
    valuesPtr: number,
    clipmode: number
  ): number;
  _ndarray_count_nonzero(ptr: number): number;
  _ndarray_nonzero(ptr: number): number;
  _ndarray_flatnonzero(ptr: number): number;
  _ndarray_where(conditionPtr: number, xPtr: number, yPtr: number): number;
  _ndarray_compress(conditionPtr: number, ptr: number, axis: number): number;
  _ndarray_extract(conditionPtr: number, ptr: number): number;
  _ndarray_choose(
    indicesPtr: number,
    choicesPtr: number,
    numChoices: number,
    clipmode: number
  ): number;
  _ndarray_diagonal(
    ptr: number,
    offset: number,
    axis1: number,
    axis2: number
  ): number;

  // DType functions
  _dtype_size(dtype: number): number;
  _dtype_is_integer(dtype: number): boolean;
  _dtype_is_floating(dtype: number): boolean;
  _dtype_is_complex(dtype: number): boolean;
  _dtype_is_signed(dtype: number): boolean;
  _dtype_is_bool(dtype: number): boolean;
  _dtype_promote(dtype1: number, dtype2: number): number;
  _dtype_can_cast(from: number, to: number, casting: number): boolean;

  // Logic functions (Phase 7) - Element-wise predicates
  _ndarray_isfinite(ptr: number): number;
  _ndarray_isinf(ptr: number): number;
  _ndarray_isnan(ptr: number): number;
  _ndarray_isneginf(ptr: number): number;
  _ndarray_isposinf(ptr: number): number;
  _ndarray_iscomplex_elem(ptr: number): number;
  _ndarray_isreal_elem(ptr: number): number;

  // Logic functions (Phase 7) - Reductions
  _ndarray_all(ptr: number): number;
  _ndarray_all_axis(ptr: number, axis: number, keepdims: number): number;
  _ndarray_any(ptr: number): number;
  _ndarray_any_axis(ptr: number, axis: number, keepdims: number): number;

  // Logic functions (Phase 7) - Comparison functions
  _ndarray_isclose(
    aPtr: number,
    bPtr: number,
    rtol: number,
    atol: number,
    equal_nan: number
  ): number;
  _ndarray_allclose(
    aPtr: number,
    bPtr: number,
    rtol: number,
    atol: number,
    equal_nan: number
  ): number;
  _ndarray_array_equal(
    a1Ptr: number,
    a2Ptr: number,
    equal_nan: number
  ): number;
  _ndarray_array_equiv(a1Ptr: number, a2Ptr: number): number;

  // Manipulation functions (Phase 5)
  _ndarray_concatenate(
    arraysPtr: number,
    nArrays: number,
    axis: number
  ): number;

  // Ufunc - Unary operations (Level 3)
  _ufunc_negative(ptr: number): number;
  _ufunc_positive(ptr: number): number;
  _ufunc_absolute(ptr: number): number;
  _ufunc_abs(ptr: number): number;
  _ufunc_sign(ptr: number): number;
  _ufunc_sqrt(ptr: number): number;
  _ufunc_square(ptr: number): number;
  _ufunc_cbrt(ptr: number): number;
  _ufunc_reciprocal(ptr: number): number;
  _ufunc_exp(ptr: number): number;
  _ufunc_exp2(ptr: number): number;
  _ufunc_expm1(ptr: number): number;
  _ufunc_log(ptr: number): number;
  _ufunc_log2(ptr: number): number;
  _ufunc_log10(ptr: number): number;
  _ufunc_log1p(ptr: number): number;
  _ufunc_sin(ptr: number): number;
  _ufunc_cos(ptr: number): number;
  _ufunc_tan(ptr: number): number;
  _ufunc_arcsin(ptr: number): number;
  _ufunc_arccos(ptr: number): number;
  _ufunc_arctan(ptr: number): number;
  _ufunc_sinh(ptr: number): number;
  _ufunc_cosh(ptr: number): number;
  _ufunc_tanh(ptr: number): number;
  _ufunc_arcsinh(ptr: number): number;
  _ufunc_arccosh(ptr: number): number;
  _ufunc_arctanh(ptr: number): number;
  _ufunc_floor(ptr: number): number;
  _ufunc_ceil(ptr: number): number;
  _ufunc_trunc(ptr: number): number;
  _ufunc_rint(ptr: number): number;
  _ufunc_round(ptr: number): number;
  _ufunc_degrees(ptr: number): number;
  _ufunc_radians(ptr: number): number;
  _ufunc_rad2deg(ptr: number): number;
  _ufunc_deg2rad(ptr: number): number;
  _ufunc_logical_not(ptr: number): number;
  _ufunc_invert(ptr: number): number;
  _ufunc_bitwise_not(ptr: number): number;

  // Ufunc - Binary operations (Level 3)
  _ufunc_add(ptr1: number, ptr2: number): number;
  _ufunc_subtract(ptr1: number, ptr2: number): number;
  _ufunc_multiply(ptr1: number, ptr2: number): number;
  _ufunc_divide(ptr1: number, ptr2: number): number;
  _ufunc_true_divide(ptr1: number, ptr2: number): number;
  _ufunc_floor_divide(ptr1: number, ptr2: number): number;
  _ufunc_remainder(ptr1: number, ptr2: number): number;
  _ufunc_mod(ptr1: number, ptr2: number): number;
  _ufunc_fmod(ptr1: number, ptr2: number): number;
  _ufunc_power(ptr1: number, ptr2: number): number;
  _ufunc_equal(ptr1: number, ptr2: number): number;
  _ufunc_not_equal(ptr1: number, ptr2: number): number;
  _ufunc_less(ptr1: number, ptr2: number): number;
  _ufunc_less_equal(ptr1: number, ptr2: number): number;
  _ufunc_greater(ptr1: number, ptr2: number): number;
  _ufunc_greater_equal(ptr1: number, ptr2: number): number;
  _ufunc_maximum(ptr1: number, ptr2: number): number;
  _ufunc_minimum(ptr1: number, ptr2: number): number;
  _ufunc_fmax(ptr1: number, ptr2: number): number;
  _ufunc_fmin(ptr1: number, ptr2: number): number;
  _ufunc_logical_and(ptr1: number, ptr2: number): number;
  _ufunc_logical_or(ptr1: number, ptr2: number): number;
  _ufunc_logical_xor(ptr1: number, ptr2: number): number;
  _ufunc_bitwise_and(ptr1: number, ptr2: number): number;
  _ufunc_bitwise_or(ptr1: number, ptr2: number): number;
  _ufunc_bitwise_xor(ptr1: number, ptr2: number): number;
  _ufunc_left_shift(ptr1: number, ptr2: number): number;
  _ufunc_right_shift(ptr1: number, ptr2: number): number;
  _ufunc_arctan2(ptr1: number, ptr2: number): number;
  _ufunc_hypot(ptr1: number, ptr2: number): number;
  _ufunc_copysign(ptr1: number, ptr2: number): number;
  _ufunc_signbit(ptr: number): number;
  _ufunc_logaddexp(ptr1: number, ptr2: number): number;
  _ufunc_logaddexp2(ptr1: number, ptr2: number): number;

  // Sorting functions (Phase 6)
  _ndarray_sort(ptr: number, axis: number, kind: number): number;
  _ndarray_sort_copy(ptr: number, axis: number, kind: number): number;
  _ndarray_argsort(ptr: number, axis: number, kind: number): number;
  _ndarray_partition(ptr: number, kth: number, axis: number): number;
  _ndarray_argpartition(ptr: number, kth: number, axis: number): number;

  // Searching functions (Phase 6)
  _ndarray_argmax(ptr: number, axis: number, keepdims: boolean): number;
  _ndarray_argmin(ptr: number, axis: number, keepdims: boolean): number;
  _ndarray_searchsorted(
    sortedPtr: number,
    valuesPtr: number,
    side: number,
    sorterPtr: number
  ): number;

  // Statistics functions (Phase 6)
  _ndarray_sum_axis(
    ptr: number,
    axis: number,
    keepdims: boolean,
    dtype: number
  ): number;
  _ndarray_mean_axis(
    ptr: number,
    axis: number,
    keepdims: boolean,
    dtype: number
  ): number;
  _ndarray_var_axis(
    ptr: number,
    axis: number,
    keepdims: boolean,
    dtype: number,
    ddof: number
  ): number;
  _ndarray_std_axis(
    ptr: number,
    axis: number,
    keepdims: boolean,
    dtype: number,
    ddof: number
  ): number;
  _ndarray_min_axis(ptr: number, axis: number, keepdims: boolean): number;
  _ndarray_max_axis(ptr: number, axis: number, keepdims: boolean): number;
  _ndarray_median(ptr: number, axis: number, keepdims: boolean): number;
  _ndarray_percentile(
    ptr: number,
    q: number,
    axis: number,
    keepdims: boolean
  ): number;
  _ndarray_quantile(
    ptr: number,
    q: number,
    axis: number,
    keepdims: boolean
  ): number;
  _ndarray_nansum(
    ptr: number,
    axis: number,
    keepdims: boolean,
    dtype: number
  ): number;
  _ndarray_nanmean(
    ptr: number,
    axis: number,
    keepdims: boolean,
    dtype: number
  ): number;
  _ndarray_nanvar(
    ptr: number,
    axis: number,
    keepdims: boolean,
    dtype: number,
    ddof: number
  ): number;
  _ndarray_nanstd(
    ptr: number,
    axis: number,
    keepdims: boolean,
    dtype: number,
    ddof: number
  ): number;

  // Set operations (Phase 8)
  _ndarray_unique(
    ptr: number,
    returnIndex: boolean,
    returnInverse: boolean,
    returnCounts: boolean,
    equalNan: boolean
  ): number;
  _ndarray_unique_values(ptr: number): number;
  _unique_result_free(ptr: number): void;
  _ndarray_union1d(ptr1: number, ptr2: number): number;
  _ndarray_intersect1d(
    ptr1: number,
    ptr2: number,
    assumeUnique: boolean,
    returnIndices: boolean,
    indices1Out: number,
    indices2Out: number
  ): number;
  _ndarray_setdiff1d(ptr1: number, ptr2: number, assumeUnique: boolean): number;
  _ndarray_setxor1d(ptr1: number, ptr2: number, assumeUnique: boolean): number;
  _ndarray_isin(
    ptr1: number,
    ptr2: number,
    invert: boolean,
    assumeUnique: boolean,
    kind: number
  ): number;
  _ndarray_in1d(
    ptr1: number,
    ptr2: number,
    invert: boolean,
    assumeUnique: boolean,
    kind: number
  ): number;

  // Memory management
  _wasm_malloc(size: number): number;
  _wasm_free(ptr: number): void;
  _malloc(size: number): number;
  _free(ptr: number): void;

  // Emscripten runtime methods
  getValue(ptr: number, type: string): number;
  setValue(ptr: number, value: number, type: string): void;

  // Heap views for direct memory access
  HEAPF64: Float64Array;
  HEAPF32: Float32Array;
  HEAP32: Int32Array;
  HEAP16: Int16Array;
  HEAP8: Int8Array;
  HEAPU32: Uint32Array;
  HEAPU16: Uint16Array;
  HEAPU8: Uint8Array;
}

/**
 * Factory function type for the WASM module
 */
export type WasmModuleFactory = () => Promise<WasmModule>;
