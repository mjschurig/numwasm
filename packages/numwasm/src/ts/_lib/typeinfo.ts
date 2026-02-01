/**
 * Type Information Utilities for NumJS
 *
 * Provides finfo and iinfo classes for inspecting floating-point
 * and integer type properties, similar to NumPy's numpy.finfo and numpy.iinfo.
 *
 * Reference: numpy/_core/getlimits.py
 */

import { DType } from '../_core/types.js';

/* ============ Floating-Point Type Information ============ */

/**
 * Machine limits for floating-point types.
 *
 * Provides properties describing the numerical limits and precision
 * of floating-point data types.
 *
 * @example
 * import { finfo, DType } from 'numjs';
 *
 * const f32 = new finfo(DType.Float32);
 * console.log(f32.eps);     // Machine epsilon
 * console.log(f32.min);     // Smallest positive normal
 * console.log(f32.max);     // Largest representable number
 * console.log(f32.bits);    // Number of bits
 *
 * // Also accepts string dtype
 * const f64 = new finfo('float64');
 *
 * @see https://numpy.org/doc/stable/reference/generated/numpy.finfo.html
 */
export class finfo {
  /** Number of bits in the type */
  readonly bits: number;

  /** Machine epsilon: smallest x such that 1.0 + x != 1.0 */
  readonly eps: number;

  /** The largest usable exponent for this type */
  readonly maxexp: number;

  /** The smallest usable exponent for this type */
  readonly minexp: number;

  /** The largest representable number */
  readonly max: number;

  /** The smallest positive usable number (normal) */
  readonly min: number;

  /** The smallest positive number (subnormal or normal) */
  readonly tiny: number;

  /** Number of decimal digits of precision */
  readonly precision: number;

  /** Number of bits in the exponent */
  readonly iexp: number;

  /** Number of bits in the mantissa (significand) */
  readonly nmant: number;

  /** Machine parameters for the dtype */
  readonly dtype: DType;

  constructor(dtype: DType | 'float16' | 'float32' | 'float64') {
    // Resolve string dtype to enum
    let resolvedDtype: DType;
    if (typeof dtype === 'string') {
      switch (dtype) {
        case 'float16':
          resolvedDtype = DType.Float16;
          break;
        case 'float32':
          resolvedDtype = DType.Float32;
          break;
        case 'float64':
          resolvedDtype = DType.Float64;
          break;
        default:
          throw new Error(`Unknown float dtype: ${dtype}`);
      }
    } else {
      resolvedDtype = dtype;
    }

    this.dtype = resolvedDtype;

    switch (resolvedDtype) {
      case DType.Float16:
        // IEEE 754 half precision
        this.bits = 16;
        this.eps = 0.0009765625; // 2^-10
        this.maxexp = 16;
        this.minexp = -14;
        this.max = 65504.0;
        this.min = 6.103515625e-5; // Smallest positive normal
        this.tiny = 6.103515625e-5;
        this.precision = 3;
        this.iexp = 5;
        this.nmant = 10;
        break;

      case DType.Float32:
        // IEEE 754 single precision
        this.bits = 32;
        this.eps = 1.1920928955078125e-7; // 2^-23
        this.maxexp = 128;
        this.minexp = -126;
        this.max = 3.4028234663852886e38;
        this.min = 1.1754943508222875e-38;
        this.tiny = 1.1754943508222875e-38;
        this.precision = 6;
        this.iexp = 8;
        this.nmant = 23;
        break;

      case DType.Float64:
        // IEEE 754 double precision
        this.bits = 64;
        this.eps = 2.220446049250313e-16; // 2^-52
        this.maxexp = 1024;
        this.minexp = -1022;
        this.max = 1.7976931348623157e308;
        this.min = 2.2250738585072014e-308;
        this.tiny = 2.2250738585072014e-308;
        this.precision = 15;
        this.iexp = 11;
        this.nmant = 52;
        break;

      default:
        throw new Error(
          `finfo not available for dtype: ${resolvedDtype}. Use float16, float32, or float64.`
        );
    }
  }

  /**
   * String representation
   */
  toString(): string {
    return (
      `finfo(dtype=${this.dtype}, bits=${this.bits}, eps=${this.eps}, ` +
      `max=${this.max}, min=${this.min})`
    );
  }
}

/* ============ Integer Type Information ============ */

/**
 * Machine limits for integer types.
 *
 * Provides properties describing the limits of integer data types.
 *
 * @example
 * import { iinfo, DType } from 'numjs';
 *
 * const i32 = new iinfo(DType.Int32);
 * console.log(i32.min);  // -2147483648
 * console.log(i32.max);  // 2147483647
 * console.log(i32.bits); // 32
 *
 * // Also accepts string dtype
 * const u8 = new iinfo('uint8');
 *
 * @see https://numpy.org/doc/stable/reference/generated/numpy.iinfo.html
 */
export class iinfo {
  /** Number of bits in the type */
  readonly bits: number;

  /** Minimum value of the type */
  readonly min: number;

  /** Maximum value of the type */
  readonly max: number;

  /** Integer type being queried */
  readonly dtype: DType;

  constructor(
    dtype:
      | DType
      | 'bool'
      | 'int8'
      | 'int16'
      | 'int32'
      | 'int64'
      | 'uint8'
      | 'uint16'
      | 'uint32'
      | 'uint64'
  ) {
    // Resolve string dtype to enum
    let resolvedDtype: DType;
    if (typeof dtype === 'string') {
      switch (dtype) {
        case 'bool':
          resolvedDtype = DType.Bool;
          break;
        case 'int8':
          resolvedDtype = DType.Int8;
          break;
        case 'int16':
          resolvedDtype = DType.Int16;
          break;
        case 'int32':
          resolvedDtype = DType.Int32;
          break;
        case 'int64':
          resolvedDtype = DType.Int64;
          break;
        case 'uint8':
          resolvedDtype = DType.Uint8;
          break;
        case 'uint16':
          resolvedDtype = DType.Uint16;
          break;
        case 'uint32':
          resolvedDtype = DType.Uint32;
          break;
        case 'uint64':
          resolvedDtype = DType.Uint64;
          break;
        default:
          throw new Error(`Unknown integer dtype: ${dtype}`);
      }
    } else {
      resolvedDtype = dtype;
    }

    this.dtype = resolvedDtype;

    switch (resolvedDtype) {
      case DType.Bool:
        this.bits = 8;
        this.min = 0;
        this.max = 1;
        break;

      case DType.Int8:
        this.bits = 8;
        this.min = -128;
        this.max = 127;
        break;

      case DType.Uint8:
        this.bits = 8;
        this.min = 0;
        this.max = 255;
        break;

      case DType.Int16:
        this.bits = 16;
        this.min = -32768;
        this.max = 32767;
        break;

      case DType.Uint16:
        this.bits = 16;
        this.min = 0;
        this.max = 65535;
        break;

      case DType.Int32:
        this.bits = 32;
        this.min = -2147483648;
        this.max = 2147483647;
        break;

      case DType.Uint32:
        this.bits = 32;
        this.min = 0;
        this.max = 4294967295;
        break;

      case DType.Int64:
        // Note: JavaScript number can't represent full int64 range precisely
        // Values beyond Number.MAX_SAFE_INTEGER may lose precision
        this.bits = 64;
        this.min = -9223372036854775808;
        this.max = 9223372036854775807;
        break;

      case DType.Uint64:
        // Note: JavaScript number can't represent full uint64 range precisely
        this.bits = 64;
        this.min = 0;
        this.max = 18446744073709551615;
        break;

      default:
        throw new Error(
          `iinfo not available for dtype: ${resolvedDtype}. Use an integer dtype.`
        );
    }
  }

  /**
   * String representation
   */
  toString(): string {
    return `iinfo(dtype=${this.dtype}, bits=${this.bits}, min=${this.min}, max=${this.max})`;
  }
}

/* ============ Casting Mode Constants ============ */

/** No casting at all */
export const CASTING_NO = 0;
/** Only equivalent dtypes allowed */
export const CASTING_EQUIV = 1;
/** Casting between dtypes of the same kind allowed */
export const CASTING_SAME_KIND = 2;
/** Safe casts (no precision loss) allowed */
export const CASTING_SAFE = 3;
/** Any cast is allowed */
export const CASTING_UNSAFE = 4;

/** Casting mode type */
export type CastingMode = 'no' | 'equiv' | 'same_kind' | 'safe' | 'unsafe';

/**
 * Convert casting mode string to integer constant.
 */
function castingToInt(casting: CastingMode): number {
  switch (casting) {
    case 'no':
      return CASTING_NO;
    case 'equiv':
      return CASTING_EQUIV;
    case 'same_kind':
      return CASTING_SAME_KIND;
    case 'safe':
      return CASTING_SAFE;
    case 'unsafe':
      return CASTING_UNSAFE;
    default:
      throw new Error(`Unknown casting mode: ${casting}`);
  }
}

/* ============ Dtype Categories ============ */

/** Integer dtypes */
const INTEGER_DTYPES = new Set([
  DType.Int8,
  DType.Int16,
  DType.Int32,
  DType.Int64,
  DType.Uint8,
  DType.Uint16,
  DType.Uint32,
  DType.Uint64,
]);

/** Signed integer dtypes */
const SIGNED_INTEGER_DTYPES = new Set([
  DType.Int8,
  DType.Int16,
  DType.Int32,
  DType.Int64,
]);

/** Unsigned integer dtypes */
const UNSIGNED_INTEGER_DTYPES = new Set([
  DType.Uint8,
  DType.Uint16,
  DType.Uint32,
  DType.Uint64,
]);

/** Floating-point dtypes */
const FLOATING_DTYPES = new Set([DType.Float16, DType.Float32, DType.Float64]);

/** Complex dtypes */
const COMPLEX_DTYPES = new Set([DType.Complex64, DType.Complex128]);

/** Numeric dtypes (all types that can participate in arithmetic) */
const NUMERIC_DTYPES = new Set([
  ...INTEGER_DTYPES,
  ...FLOATING_DTYPES,
  ...COMPLEX_DTYPES,
]);

/* ============ Type Utility Functions ============ */

/**
 * Returns true if cast between data types can occur according to the casting rule.
 *
 * @param from - Data type to cast from
 * @param to - Data type to cast to
 * @param casting - Controls what kind of data casting may occur:
 *   - 'no': No casting allowed
 *   - 'equiv': Only byte-order changes allowed
 *   - 'safe': Only casts that can preserve values allowed
 *   - 'same_kind': Only safe casts or casts within a kind allowed
 *   - 'unsafe': Any data conversions allowed
 * @returns True if the cast is allowed
 *
 * @example
 * ```typescript
 * can_cast(DType.Int32, DType.Float64, 'safe');  // true
 * can_cast(DType.Float64, DType.Int32, 'safe');  // false
 * can_cast(DType.Float64, DType.Int32, 'unsafe'); // true
 * ```
 */
export function can_cast(
  from: DType,
  to: DType,
  casting: CastingMode = 'safe',
): boolean {
  if (from === to) {
    return true;
  }

  const castingInt = castingToInt(casting);

  // 'no' casting - only same dtype
  if (castingInt === CASTING_NO) {
    return false;
  }

  // 'equiv' - same dtype or equivalent (e.g., same size, byte order changes)
  if (castingInt === CASTING_EQUIV) {
    return false; // In WASM, we don't have byte-order variations
  }

  // 'unsafe' - anything goes
  if (castingInt === CASTING_UNSAFE) {
    return true;
  }

  // 'safe' or 'same_kind'
  const fromInfo = getDtypeInfo(from);
  const toInfo = getDtypeInfo(to);

  // Bool can cast safely to any numeric
  if (from === DType.Bool && NUMERIC_DTYPES.has(to)) {
    return true;
  }

  // Check 'same_kind' first
  if (castingInt === CASTING_SAME_KIND) {
    // Allow casts within the same kind
    if (fromInfo.kind === toInfo.kind) {
      return true;
    }
    // Also allow safe casts
  }

  // Check for safe casting (no precision loss)
  // Integer to larger integer (signed or unsigned of same sign)
  if (INTEGER_DTYPES.has(from) && INTEGER_DTYPES.has(to)) {
    if (SIGNED_INTEGER_DTYPES.has(from) && SIGNED_INTEGER_DTYPES.has(to)) {
      return toInfo.bits >= fromInfo.bits;
    }
    if (UNSIGNED_INTEGER_DTYPES.has(from) && UNSIGNED_INTEGER_DTYPES.has(to)) {
      return toInfo.bits >= fromInfo.bits;
    }
    // Unsigned to signed: need extra bit
    if (UNSIGNED_INTEGER_DTYPES.has(from) && SIGNED_INTEGER_DTYPES.has(to)) {
      return toInfo.bits > fromInfo.bits;
    }
    // Signed to unsigned: never safe (negative values)
    return false;
  }

  // Integer to float: safe if float can represent all integer values
  if (INTEGER_DTYPES.has(from) && FLOATING_DTYPES.has(to)) {
    // Float64 can safely hold int32 and below
    // Float32 can safely hold int16 and below
    if (to === DType.Float64) {
      return fromInfo.bits <= 53; // 53 mantissa bits
    }
    if (to === DType.Float32) {
      return fromInfo.bits <= 24; // 24 mantissa bits (including implicit bit)
    }
    if (to === DType.Float16) {
      return fromInfo.bits <= 11; // 11 mantissa bits
    }
    return false;
  }

  // Float to larger float
  if (FLOATING_DTYPES.has(from) && FLOATING_DTYPES.has(to)) {
    return toInfo.bits >= fromInfo.bits;
  }

  // Float to integer: never safe
  if (FLOATING_DTYPES.has(from) && INTEGER_DTYPES.has(to)) {
    return false;
  }

  // Complex handling
  if (COMPLEX_DTYPES.has(from)) {
    // Complex to complex
    if (COMPLEX_DTYPES.has(to)) {
      return toInfo.bits >= fromInfo.bits;
    }
    // Complex to real: never safe
    return false;
  }

  // Real to complex: safe if underlying float is safe
  if (COMPLEX_DTYPES.has(to)) {
    const complexFloatBits = to === DType.Complex64 ? 32 : 64;
    if (FLOATING_DTYPES.has(from)) {
      return fromInfo.bits <= complexFloatBits;
    }
    if (INTEGER_DTYPES.has(from)) {
      // Check if integer fits in the complex's float component
      const mantissaBits = complexFloatBits === 64 ? 53 : 24;
      return fromInfo.bits <= mantissaBits;
    }
  }

  return false;
}

/**
 * Get information about a dtype.
 */
function getDtypeInfo(dtype: DType): { bits: number; kind: string } {
  switch (dtype) {
    case DType.Bool:
      return { bits: 1, kind: 'b' };
    case DType.Int8:
      return { bits: 8, kind: 'i' };
    case DType.Int16:
      return { bits: 16, kind: 'i' };
    case DType.Int32:
      return { bits: 32, kind: 'i' };
    case DType.Int64:
      return { bits: 64, kind: 'i' };
    case DType.Uint8:
      return { bits: 8, kind: 'u' };
    case DType.Uint16:
      return { bits: 16, kind: 'u' };
    case DType.Uint32:
      return { bits: 32, kind: 'u' };
    case DType.Uint64:
      return { bits: 64, kind: 'u' };
    case DType.Float16:
      return { bits: 16, kind: 'f' };
    case DType.Float32:
      return { bits: 32, kind: 'f' };
    case DType.Float64:
      return { bits: 64, kind: 'f' };
    case DType.Complex64:
      return { bits: 64, kind: 'c' };
    case DType.Complex128:
      return { bits: 128, kind: 'c' };
    default:
      return { bits: 0, kind: 'O' };
  }
}

/**
 * Return true if first argument is a typecode lower/equal in type hierarchy.
 *
 * @param arg1 - dtype to check
 * @param arg2 - dtype category or abstract type to compare against
 * @returns True if arg1 is a subtype of arg2
 *
 * @example
 * ```typescript
 * issubdtype(DType.Float32, 'floating');     // true
 * issubdtype(DType.Int32, 'integer');        // true
 * issubdtype(DType.Float32, 'integer');      // false
 * issubdtype(DType.Int32, 'signedinteger');  // true
 * issubdtype(DType.Uint32, 'signedinteger'); // false
 * ```
 */
export function issubdtype(
  arg1: DType,
  arg2: DType | 'integer' | 'signedinteger' | 'unsignedinteger' | 'floating' | 'complexfloating' | 'number' | 'bool',
): boolean {
  // If arg2 is a DType enum, check for exact match or subtype relation
  if (typeof arg2 === 'number') {
    return arg1 === arg2;
  }

  // Check against dtype category strings
  switch (arg2) {
    case 'bool':
      return arg1 === DType.Bool;
    case 'integer':
      return INTEGER_DTYPES.has(arg1);
    case 'signedinteger':
      return SIGNED_INTEGER_DTYPES.has(arg1);
    case 'unsignedinteger':
      return UNSIGNED_INTEGER_DTYPES.has(arg1);
    case 'floating':
      return FLOATING_DTYPES.has(arg1);
    case 'complexfloating':
      return COMPLEX_DTYPES.has(arg1);
    case 'number':
      return NUMERIC_DTYPES.has(arg1);
    default:
      return false;
  }
}

/**
 * Returns a boolean indicating whether the provided dtype is of the specified kind.
 *
 * This is the NumPy 2.0 version of issubdtype, with a more explicit API.
 *
 * @param dtype - The dtype to check
 * @param kind - The dtype kind to check against
 * @returns True if dtype matches the specified kind
 *
 * @example
 * ```typescript
 * isdtype(DType.Float32, 'real floating');   // true
 * isdtype(DType.Int32, 'integral');          // true
 * isdtype(DType.Complex128, 'complex floating'); // true
 * ```
 */
export function isdtype(
  dtype: DType,
  kind: 'bool' | 'integral' | 'real floating' | 'complex floating' | 'numeric',
): boolean {
  switch (kind) {
    case 'bool':
      return dtype === DType.Bool;
    case 'integral':
      return INTEGER_DTYPES.has(dtype);
    case 'real floating':
      return FLOATING_DTYPES.has(dtype);
    case 'complex floating':
      return COMPLEX_DTYPES.has(dtype);
    case 'numeric':
      return NUMERIC_DTYPES.has(dtype);
    default:
      return false;
  }
}

/**
 * Determine common type following standard coercion rules.
 *
 * @param arrays_or_dtypes - Array of dtypes or arrays to find common type for
 * @returns The common dtype that can represent all inputs
 *
 * @example
 * ```typescript
 * common_type(DType.Float32, DType.Float64);  // DType.Float64
 * common_type(DType.Int32, DType.Float32);    // DType.Float64
 * common_type(DType.Int32, DType.Int64);      // DType.Int64
 * ```
 */
export function common_type(...dtypes: DType[]): DType {
  if (dtypes.length === 0) {
    return DType.Float64; // Default
  }

  if (dtypes.length === 1) {
    return dtypes[0];
  }

  let result = dtypes[0];
  for (let i = 1; i < dtypes.length; i++) {
    result = promoteTypes(result, dtypes[i]);
  }

  return result;
}

/**
 * Promote two dtypes to a common type that can safely hold both.
 */
function promoteTypes(dtype1: DType, dtype2: DType): DType {
  if (dtype1 === dtype2) {
    return dtype1;
  }

  // Complex promotion
  if (COMPLEX_DTYPES.has(dtype1) || COMPLEX_DTYPES.has(dtype2)) {
    // If either is complex128, result is complex128
    if (dtype1 === DType.Complex128 || dtype2 === DType.Complex128) {
      return DType.Complex128;
    }
    // If both involve complex64 and no larger float, return complex64
    // If one is float64 or large integer, promote to complex128
    if (dtype1 === DType.Float64 || dtype2 === DType.Float64) {
      return DType.Complex128;
    }
    if (INTEGER_DTYPES.has(dtype1) && getDtypeInfo(dtype1).bits > 32) {
      return DType.Complex128;
    }
    if (INTEGER_DTYPES.has(dtype2) && getDtypeInfo(dtype2).bits > 32) {
      return DType.Complex128;
    }
    return DType.Complex64;
  }

  // Float promotion
  if (FLOATING_DTYPES.has(dtype1) || FLOATING_DTYPES.has(dtype2)) {
    const info1 = getDtypeInfo(dtype1);
    const info2 = getDtypeInfo(dtype2);

    // Both floats: take larger
    if (FLOATING_DTYPES.has(dtype1) && FLOATING_DTYPES.has(dtype2)) {
      return info1.bits >= info2.bits ? dtype1 : dtype2;
    }

    // Float + integer: promote based on integer size
    if (FLOATING_DTYPES.has(dtype1)) {
      if (info2.bits <= 16) return dtype1;
      if (info2.bits <= 32) return info1.bits >= 64 ? dtype1 : DType.Float64;
      return DType.Float64;
    }
    if (FLOATING_DTYPES.has(dtype2)) {
      if (info1.bits <= 16) return dtype2;
      if (info1.bits <= 32) return info2.bits >= 64 ? dtype2 : DType.Float64;
      return DType.Float64;
    }
  }

  // Integer promotion
  if (INTEGER_DTYPES.has(dtype1) && INTEGER_DTYPES.has(dtype2)) {
    const info1 = getDtypeInfo(dtype1);
    const info2 = getDtypeInfo(dtype2);

    // Same signedness: take larger
    if (
      (SIGNED_INTEGER_DTYPES.has(dtype1) &&
        SIGNED_INTEGER_DTYPES.has(dtype2)) ||
      (UNSIGNED_INTEGER_DTYPES.has(dtype1) &&
        UNSIGNED_INTEGER_DTYPES.has(dtype2))
    ) {
      return info1.bits >= info2.bits ? dtype1 : dtype2;
    }

    // Mixed signedness: promote to signed type that can hold both
    const unsignedBits = UNSIGNED_INTEGER_DTYPES.has(dtype1)
      ? info1.bits
      : info2.bits;
    if (unsignedBits < 64) {
      // Find next larger signed type
      if (unsignedBits <= 8) return DType.Int16;
      if (unsignedBits <= 16) return DType.Int32;
      if (unsignedBits <= 32) return DType.Int64;
    }
    // Fallback to float64 for uint64 + signed
    return DType.Float64;
  }

  // Bool promotion
  if (dtype1 === DType.Bool) return dtype2;
  if (dtype2 === DType.Bool) return dtype1;

  // Fallback
  return DType.Float64;
}

/**
 * Return the type that results from applying the type promotion rules to the arguments.
 *
 * @param dtypes - Input dtypes
 * @returns The result dtype from type promotion
 *
 * @example
 * ```typescript
 * result_type(DType.Int32, DType.Float32);  // DType.Float64
 * result_type(DType.Int16, DType.Int32);    // DType.Int32
 * ```
 */
export function result_type(...dtypes: DType[]): DType {
  return common_type(...dtypes);
}

/**
 * Return the minimum scalar type needed to represent the value.
 *
 * For scalars, returns the smallest dtype that can hold the value.
 * For integers, finds the smallest integer type.
 * For floats, returns float64 (to preserve precision).
 *
 * @param value - A scalar value
 * @returns The minimum dtype needed to represent the value
 *
 * @example
 * ```typescript
 * min_scalar_type(10);      // DType.Uint8 (fits in 0-255)
 * min_scalar_type(300);     // DType.Uint16 (doesn't fit in Uint8)
 * min_scalar_type(-1);      // DType.Int8 (negative, fits in -128 to 127)
 * min_scalar_type(3.14);    // DType.Float64
 * min_scalar_type(true);    // DType.Bool
 * ```
 */
export function min_scalar_type(value: number | boolean): DType {
  // Boolean
  if (typeof value === "boolean") {
    return DType.Bool;
  }

  // Non-finite or NaN: return Float64
  if (!Number.isFinite(value)) {
    return DType.Float64;
  }

  // Floating point
  if (!Number.isInteger(value)) {
    // Could potentially return Float32 if value fits, but for safety return Float64
    return DType.Float64;
  }

  // Integer handling
  const intVal = value;

  if (intVal >= 0) {
    // Unsigned integers
    if (intVal <= 255) {
      return DType.Uint8;
    }
    if (intVal <= 65535) {
      return DType.Uint16;
    }
    if (intVal <= 4294967295) {
      return DType.Uint32;
    }
    // For very large integers, use Int64 or Float64
    if (intVal <= Number.MAX_SAFE_INTEGER) {
      return DType.Uint64;
    }
    return DType.Float64;
  } else {
    // Signed integers (negative values)
    if (intVal >= -128) {
      return DType.Int8;
    }
    if (intVal >= -32768) {
      return DType.Int16;
    }
    if (intVal >= -2147483648) {
      return DType.Int32;
    }
    if (intVal >= Number.MIN_SAFE_INTEGER) {
      return DType.Int64;
    }
    return DType.Float64;
  }
}
