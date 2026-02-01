/**
 * Type Information Utilities for NumJS
 *
 * Re-exports from _lib/typeinfo.js for backwards compatibility.
 */

export {
  finfo,
  iinfo,
  can_cast,
  issubdtype,
  isdtype,
  common_type,
  result_type,
  min_scalar_type,
} from "./_lib/typeinfo.js";
export type { CastingMode } from "./_lib/typeinfo.js";
