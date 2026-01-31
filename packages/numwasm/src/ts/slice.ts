/**
 * Slice and indexing utilities for NumJS-WASM
 *
 * Re-exports from _core/slice.js for backwards compatibility.
 */

export {
  Slice,
  slice,
  ellipsis,
  newaxis,
  expandEllipsis,
  buildIndexSpecs,
  computeResultShape,
  INDEX_TYPE_INTEGER,
  INDEX_TYPE_SLICE,
  INDEX_TYPE_NEWAXIS,
  INDEX_TYPE_ELLIPSIS,
} from "./_core/slice.js";
export type { IndexElement, IndexSpec, Ellipsis, Newaxis } from "./_core/slice.js";
