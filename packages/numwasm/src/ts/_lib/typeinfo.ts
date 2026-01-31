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
