/**
 * NumJS DType Utilities
 *
 * TypeScript utilities for working with data types.
 */

import { DType, DTYPE_NAMES, DTYPE_SIZES, type TypedArrayType } from './types.js';

/**
 * Map from string name to DType enum value
 */
const NAME_TO_DTYPE: Record<string, DType> = {
  float32: DType.Float32,
  float64: DType.Float64,
  int32: DType.Int32,
  int64: DType.Int64,
  bool: DType.Bool,
  int8: DType.Int8,
  int16: DType.Int16,
  uint8: DType.Uint8,
  uint16: DType.Uint16,
  uint32: DType.Uint32,
  uint64: DType.Uint64,
  float16: DType.Float16,
  complex64: DType.Complex64,
  complex128: DType.Complex128,
  // Aliases
  f4: DType.Float32,
  f8: DType.Float64,
  i4: DType.Int32,
  i8: DType.Int64,
  i1: DType.Int8,
  i2: DType.Int16,
  u1: DType.Uint8,
  u2: DType.Uint16,
  u4: DType.Uint32,
  u8: DType.Uint64,
  f2: DType.Float16,
  c8: DType.Complex64,
  c16: DType.Complex128,
  double: DType.Float64,
  float: DType.Float32,
  int: DType.Int32,
  boolean: DType.Bool,
};

/**
 * Convert string dtype name to DType enum.
 *
 * @param name - String name (e.g., 'float64', 'int32')
 * @returns DType enum value
 * @throws Error if name is not recognized
 */
export function dtypeFromString(name: string): DType {
  const normalized = name.toLowerCase().trim();
  const dtype = NAME_TO_DTYPE[normalized];
  if (dtype === undefined) {
    throw new Error(`Unknown dtype: ${name}`);
  }
  return dtype;
}

/**
 * Get dtype string name.
 *
 * @param dtype - DType enum value
 * @returns String name
 */
export function dtypeToString(dtype: DType): string {
  return DTYPE_NAMES[dtype] ?? 'unknown';
}

/**
 * Check if dtype is an integer type (signed or unsigned).
 */
export function isIntegerDType(dtype: DType): boolean {
  return (
    dtype === DType.Int8 ||
    dtype === DType.Int16 ||
    dtype === DType.Int32 ||
    dtype === DType.Int64 ||
    dtype === DType.Uint8 ||
    dtype === DType.Uint16 ||
    dtype === DType.Uint32 ||
    dtype === DType.Uint64
  );
}

/**
 * Check if dtype is a floating point type.
 */
export function isFloatDType(dtype: DType): boolean {
  return (
    dtype === DType.Float16 ||
    dtype === DType.Float32 ||
    dtype === DType.Float64
  );
}

/**
 * Check if dtype is a complex type.
 */
export function isComplexDType(dtype: DType): boolean {
  return dtype === DType.Complex64 || dtype === DType.Complex128;
}

/**
 * Check if dtype is a signed type.
 */
export function isSignedDType(dtype: DType): boolean {
  return (
    dtype === DType.Int8 ||
    dtype === DType.Int16 ||
    dtype === DType.Int32 ||
    dtype === DType.Int64 ||
    dtype === DType.Float16 ||
    dtype === DType.Float32 ||
    dtype === DType.Float64 ||
    dtype === DType.Complex64 ||
    dtype === DType.Complex128
  );
}

/**
 * Check if dtype is boolean.
 */
export function isBoolDType(dtype: DType): boolean {
  return dtype === DType.Bool;
}

/**
 * Check if dtype is numeric (not bool).
 */
export function isNumericDType(dtype: DType): boolean {
  return isIntegerDType(dtype) || isFloatDType(dtype) || isComplexDType(dtype);
}

/**
 * Get the TypedArray constructor for a given dtype.
 *
 * @param dtype - DType enum value
 * @returns TypedArray constructor
 */
export function dtypeToTypedArrayConstructor(
  dtype: DType
):
  | Float32ArrayConstructor
  | Float64ArrayConstructor
  | Int32ArrayConstructor
  | Int16ArrayConstructor
  | Int8ArrayConstructor
  | Uint32ArrayConstructor
  | Uint16ArrayConstructor
  | Uint8ArrayConstructor
  | BigInt64ArrayConstructor
  | BigUint64ArrayConstructor {
  switch (dtype) {
    case DType.Float32:
      return Float32Array;
    case DType.Float64:
      return Float64Array;
    case DType.Int32:
      return Int32Array;
    case DType.Int64:
      return BigInt64Array;
    case DType.Int16:
      return Int16Array;
    case DType.Int8:
      return Int8Array;
    case DType.Uint32:
      return Uint32Array;
    case DType.Uint64:
      return BigUint64Array;
    case DType.Uint16:
      return Uint16Array;
    case DType.Uint8:
    case DType.Bool:
      return Uint8Array;
    case DType.Float16:
      return Uint16Array; // Float16 stored as uint16
    case DType.Complex64:
      return Float32Array; // Complex stored as pairs
    case DType.Complex128:
      return Float64Array;
    default:
      return Float64Array;
  }
}

/**
 * Infer DType from a TypedArray.
 *
 * @param arr - TypedArray instance
 * @returns Inferred DType
 */
export function typedArrayToDType(arr: TypedArrayType): DType {
  if (arr instanceof Float64Array) return DType.Float64;
  if (arr instanceof Float32Array) return DType.Float32;
  if (arr instanceof Int32Array) return DType.Int32;
  if (arr instanceof Int16Array) return DType.Int16;
  if (arr instanceof Int8Array) return DType.Int8;
  if (arr instanceof Uint32Array) return DType.Uint32;
  if (arr instanceof Uint16Array) return DType.Uint16;
  if (arr instanceof Uint8Array) return DType.Uint8;
  if (arr instanceof BigInt64Array) return DType.Int64;
  if (arr instanceof BigUint64Array) return DType.Uint64;
  return DType.Float64; // Default
}

/**
 * Get the size of a dtype in bytes.
 */
export function dtypeSize(dtype: DType): number {
  return DTYPE_SIZES[dtype] ?? 0;
}

/**
 * Casting safety levels
 */
export enum CastingKind {
  No = 0,
  Equiv = 1,
  Safe = 2,
  SameKind = 3,
  Unsafe = 4,
}

/**
 * Priority values for type promotion (higher wins).
 */
const DTYPE_PRIORITY: Record<DType, number> = {
  [DType.Bool]: 0,
  [DType.Uint8]: 1,
  [DType.Uint16]: 2,
  [DType.Uint32]: 3,
  [DType.Uint64]: 4,
  [DType.Int8]: 5,
  [DType.Int16]: 6,
  [DType.Int32]: 7,
  [DType.Int64]: 8,
  [DType.Float16]: 9,
  [DType.Float32]: 10,
  [DType.Float64]: 11,
  [DType.Complex64]: 12,
  [DType.Complex128]: 13,
  [DType.String]: -1, // String does not participate in numeric promotion
};

/**
 * Get result type when combining two dtypes.
 * Follows NumPy promotion rules.
 *
 * @param dtype1 - First dtype
 * @param dtype2 - Second dtype
 * @returns Result dtype
 */
export function promoteTypes(dtype1: DType, dtype2: DType): DType {
  // Same type
  if (dtype1 === dtype2) return dtype1;

  // Complex promotion
  if (isComplexDType(dtype1) || isComplexDType(dtype2)) {
    if (
      dtype1 === DType.Complex128 ||
      dtype2 === DType.Complex128 ||
      dtype1 === DType.Float64 ||
      dtype2 === DType.Float64 ||
      dtype1 === DType.Int64 ||
      dtype2 === DType.Int64 ||
      dtype1 === DType.Uint64 ||
      dtype2 === DType.Uint64
    ) {
      return DType.Complex128;
    }
    return DType.Complex64;
  }

  // Float promotion
  if (isFloatDType(dtype1) || isFloatDType(dtype2)) {
    const hasLargeInt =
      dtype1 === DType.Int64 ||
      dtype1 === DType.Uint64 ||
      dtype1 === DType.Int32 ||
      dtype1 === DType.Uint32 ||
      dtype2 === DType.Int64 ||
      dtype2 === DType.Uint64 ||
      dtype2 === DType.Int32 ||
      dtype2 === DType.Uint32;
    const hasFloat32 = dtype1 === DType.Float32 || dtype2 === DType.Float32;

    if (hasLargeInt && hasFloat32) {
      return DType.Float64;
    }

    if (dtype1 === DType.Float64 || dtype2 === DType.Float64) {
      return DType.Float64;
    }
    if (dtype1 === DType.Float32 || dtype2 === DType.Float32) {
      return DType.Float32;
    }
    return DType.Float16;
  }

  // Bool promotion
  if (dtype1 === DType.Bool) return dtype2;
  if (dtype2 === DType.Bool) return dtype1;

  // Integer promotion
  const signed1 = isSignedDType(dtype1);
  const signed2 = isSignedDType(dtype2);
  const size1 = DTYPE_SIZES[dtype1];
  const size2 = DTYPE_SIZES[dtype2];

  if (signed1 !== signed2) {
    // Mixed signed/unsigned
    const unsignedSize = signed1 ? size2 : size1;
    const signedSize = signed1 ? size1 : size2;

    if (unsignedSize >= signedSize) {
      if (unsignedSize >= 8) return DType.Float64;
      if (unsignedSize >= 4) return DType.Int64;
      if (unsignedSize >= 2) return DType.Int32;
      return DType.Int16;
    }
    return signed1 ? dtype1 : dtype2;
  }

  // Same signedness: larger wins
  return DTYPE_PRIORITY[dtype1] > DTYPE_PRIORITY[dtype2] ? dtype1 : dtype2;
}

/**
 * Get common type for multiple dtypes.
 */
export function commonType(...dtypes: DType[]): DType {
  if (dtypes.length === 0) return DType.Float64;
  let result = dtypes[0];
  for (let i = 1; i < dtypes.length; i++) {
    result = promoteTypes(result, dtypes[i]);
  }
  return result;
}
