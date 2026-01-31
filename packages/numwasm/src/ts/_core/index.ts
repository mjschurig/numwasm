/**
 * NumJS Core Module
 *
 * Core array infrastructure that everything depends on.
 * This module contains NDArray, types, dtype utilities, slicing, broadcasting, and WASM loading.
 */

// Main array class
export { NDArray } from "./NDArray.js";

// Type definitions
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
} from "./types.js";
export type {
  NDArrayOptions,
  TypedArrayType,
  WasmModule,
  Complex,
  NDArrayFlags,
  WasmModuleOptions,
  WasmModuleFactory,
  FieldDescriptor,
  StructuredDType,
} from "./types.js";

// DType utilities
export {
  dtypeFromString,
  dtypeToString,
  dtypeSize as getDtypeSize,
  isIntegerDType,
  isFloatDType,
  isComplexDType,
  isSignedDType,
  isBoolDType,
  isNumericDType,
  dtypeToTypedArrayConstructor,
  typedArrayToDType,
  promoteTypes,
  commonType,
  CastingKind,
} from "./dtype.js";

// WASM module management
export {
  loadWasmModule,
  isWasmLoaded,
  getWasmModule,
  configureWasm,
  resetWasmModule,
} from "./wasm-loader.js";
export type { WasmLoadConfig } from "./wasm-loader.js";

// Slice utilities
export {
  Slice,
  slice,
  ellipsis,
  newaxis,
  expandEllipsis,
  buildIndexSpecs,
  computeResultShape,
} from "./slice.js";
export type { IndexElement, IndexSpec, Ellipsis, Newaxis } from "./slice.js";

// Broadcasting functions
export {
  broadcastShapes,
  broadcastShapesMulti,
  shapesAreBroadcastable,
  broadcastTo,
  broadcastArrays,
  computeBroadcastStrides,
} from "./broadcast.js";
