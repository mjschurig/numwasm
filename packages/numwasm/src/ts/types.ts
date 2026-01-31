/**
 * NumJS Type Definitions
 *
 * Re-exports from _core/types.js for backwards compatibility.
 */

export {
  DType,
  DTYPE_SIZES,
  DTYPE_NAMES,
  NDARRAY_OWNDATA,
  NDARRAY_WRITEABLE,
  NDARRAY_C_CONTIGUOUS,
  NDARRAY_F_CONTIGUOUS,
  NDARRAY_ALIGNED,
  CLIP_RAISE,
  CLIP_WRAP,
  CLIP_CLIP,
  INDEX_TYPE_INTEGER,
  INDEX_TYPE_SLICE,
  INDEX_TYPE_NEWAXIS,
  INDEX_TYPE_ELLIPSIS,
  isStructuredDType,
  dtypeSize,
  dtypeAlignment,
} from "./_core/types.js";
export type {
  Complex,
  NDArrayOptions,
  TypedArrayType,
  NDArrayFlags,
  WasmModule,
  WasmModuleOptions,
  WasmModuleFactory,
  FieldDescriptor,
  StructuredDType,
} from "./_core/types.js";
