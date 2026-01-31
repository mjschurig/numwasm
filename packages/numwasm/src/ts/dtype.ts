/**
 * NumJS DType Utilities
 *
 * Re-exports from _core/dtype.js for backwards compatibility.
 */

export {
  dtypeFromString,
  dtypeToString,
  dtypeSize,
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
} from "./_core/dtype.js";
