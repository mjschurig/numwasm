/**
 * DType System Tests
 *
 * Tests for data type utilities: conversion, promotion, and predicates.
 */

import { describe, it, expect } from 'vitest';
import {
  DType,
  dtypeFromString,
  dtypeToString,
  isIntegerDType,
  isFloatDType,
  isComplexDType,
  isSignedDType,
  isBoolDType,
  isNumericDType,
  promoteTypes,
  commonType,
  dtypeSize,
  dtypeToTypedArrayConstructor,
  typedArrayToDType,
} from '../../dist/numjs.mjs';

describe('DType Enum', () => {
  it('should have all 14 types defined', () => {
    expect(DType.Float32).toBe(0);
    expect(DType.Float64).toBe(1);
    expect(DType.Int32).toBe(2);
    expect(DType.Int64).toBe(3);
    expect(DType.Bool).toBe(4);
    expect(DType.Int8).toBe(5);
    expect(DType.Int16).toBe(6);
    expect(DType.Uint8).toBe(7);
    expect(DType.Uint16).toBe(8);
    expect(DType.Uint32).toBe(9);
    expect(DType.Uint64).toBe(10);
    expect(DType.Float16).toBe(11);
    expect(DType.Complex64).toBe(12);
    expect(DType.Complex128).toBe(13);
  });
});

describe('dtypeFromString', () => {
  it('should convert standard type names', () => {
    expect(dtypeFromString('float32')).toBe(DType.Float32);
    expect(dtypeFromString('float64')).toBe(DType.Float64);
    expect(dtypeFromString('int32')).toBe(DType.Int32);
    expect(dtypeFromString('int64')).toBe(DType.Int64);
    expect(dtypeFromString('bool')).toBe(DType.Bool);
    expect(dtypeFromString('int8')).toBe(DType.Int8);
    expect(dtypeFromString('int16')).toBe(DType.Int16);
    expect(dtypeFromString('uint8')).toBe(DType.Uint8);
    expect(dtypeFromString('uint16')).toBe(DType.Uint16);
    expect(dtypeFromString('uint32')).toBe(DType.Uint32);
    expect(dtypeFromString('uint64')).toBe(DType.Uint64);
    expect(dtypeFromString('float16')).toBe(DType.Float16);
    expect(dtypeFromString('complex64')).toBe(DType.Complex64);
    expect(dtypeFromString('complex128')).toBe(DType.Complex128);
  });

  it('should handle aliases', () => {
    expect(dtypeFromString('double')).toBe(DType.Float64);
    expect(dtypeFromString('float')).toBe(DType.Float32);
    expect(dtypeFromString('int')).toBe(DType.Int32);
    expect(dtypeFromString('boolean')).toBe(DType.Bool);
    expect(dtypeFromString('f4')).toBe(DType.Float32);
    expect(dtypeFromString('f8')).toBe(DType.Float64);
    expect(dtypeFromString('i4')).toBe(DType.Int32);
    expect(dtypeFromString('i8')).toBe(DType.Int64);
  });

  it('should be case-insensitive', () => {
    expect(dtypeFromString('FLOAT32')).toBe(DType.Float32);
    expect(dtypeFromString('Float64')).toBe(DType.Float64);
    expect(dtypeFromString('INT32')).toBe(DType.Int32);
  });

  it('should throw for unknown types', () => {
    expect(() => dtypeFromString('unknown')).toThrow('Unknown dtype');
    expect(() => dtypeFromString('string')).toThrow('Unknown dtype');
  });
});

describe('dtypeToString', () => {
  it('should convert all dtypes to string', () => {
    expect(dtypeToString(DType.Float32)).toBe('float32');
    expect(dtypeToString(DType.Float64)).toBe('float64');
    expect(dtypeToString(DType.Int32)).toBe('int32');
    expect(dtypeToString(DType.Int64)).toBe('int64');
    expect(dtypeToString(DType.Bool)).toBe('bool');
    expect(dtypeToString(DType.Int8)).toBe('int8');
    expect(dtypeToString(DType.Int16)).toBe('int16');
    expect(dtypeToString(DType.Uint8)).toBe('uint8');
    expect(dtypeToString(DType.Uint16)).toBe('uint16');
    expect(dtypeToString(DType.Uint32)).toBe('uint32');
    expect(dtypeToString(DType.Uint64)).toBe('uint64');
    expect(dtypeToString(DType.Float16)).toBe('float16');
    expect(dtypeToString(DType.Complex64)).toBe('complex64');
    expect(dtypeToString(DType.Complex128)).toBe('complex128');
  });
});

describe('Type Predicates', () => {
  describe('isIntegerDType', () => {
    it('should identify integer types', () => {
      expect(isIntegerDType(DType.Int8)).toBe(true);
      expect(isIntegerDType(DType.Int16)).toBe(true);
      expect(isIntegerDType(DType.Int32)).toBe(true);
      expect(isIntegerDType(DType.Int64)).toBe(true);
      expect(isIntegerDType(DType.Uint8)).toBe(true);
      expect(isIntegerDType(DType.Uint16)).toBe(true);
      expect(isIntegerDType(DType.Uint32)).toBe(true);
      expect(isIntegerDType(DType.Uint64)).toBe(true);
    });

    it('should reject non-integer types', () => {
      expect(isIntegerDType(DType.Float32)).toBe(false);
      expect(isIntegerDType(DType.Float64)).toBe(false);
      expect(isIntegerDType(DType.Bool)).toBe(false);
      expect(isIntegerDType(DType.Complex64)).toBe(false);
    });
  });

  describe('isFloatDType', () => {
    it('should identify float types', () => {
      expect(isFloatDType(DType.Float16)).toBe(true);
      expect(isFloatDType(DType.Float32)).toBe(true);
      expect(isFloatDType(DType.Float64)).toBe(true);
    });

    it('should reject non-float types', () => {
      expect(isFloatDType(DType.Int32)).toBe(false);
      expect(isFloatDType(DType.Bool)).toBe(false);
      expect(isFloatDType(DType.Complex64)).toBe(false);
    });
  });

  describe('isComplexDType', () => {
    it('should identify complex types', () => {
      expect(isComplexDType(DType.Complex64)).toBe(true);
      expect(isComplexDType(DType.Complex128)).toBe(true);
    });

    it('should reject non-complex types', () => {
      expect(isComplexDType(DType.Float64)).toBe(false);
      expect(isComplexDType(DType.Int32)).toBe(false);
    });
  });

  describe('isSignedDType', () => {
    it('should identify signed types', () => {
      expect(isSignedDType(DType.Int8)).toBe(true);
      expect(isSignedDType(DType.Int16)).toBe(true);
      expect(isSignedDType(DType.Int32)).toBe(true);
      expect(isSignedDType(DType.Int64)).toBe(true);
      expect(isSignedDType(DType.Float32)).toBe(true);
      expect(isSignedDType(DType.Float64)).toBe(true);
      expect(isSignedDType(DType.Complex64)).toBe(true);
    });

    it('should reject unsigned types', () => {
      expect(isSignedDType(DType.Uint8)).toBe(false);
      expect(isSignedDType(DType.Uint16)).toBe(false);
      expect(isSignedDType(DType.Uint32)).toBe(false);
      expect(isSignedDType(DType.Uint64)).toBe(false);
      expect(isSignedDType(DType.Bool)).toBe(false);
    });
  });

  describe('isBoolDType', () => {
    it('should identify bool type', () => {
      expect(isBoolDType(DType.Bool)).toBe(true);
    });

    it('should reject non-bool types', () => {
      expect(isBoolDType(DType.Int32)).toBe(false);
      expect(isBoolDType(DType.Uint8)).toBe(false);
    });
  });

  describe('isNumericDType', () => {
    it('should identify numeric types', () => {
      expect(isNumericDType(DType.Int32)).toBe(true);
      expect(isNumericDType(DType.Float64)).toBe(true);
      expect(isNumericDType(DType.Complex128)).toBe(true);
    });

    it('should reject bool', () => {
      expect(isNumericDType(DType.Bool)).toBe(false);
    });
  });
});

describe('dtypeSize', () => {
  it('should return correct sizes', () => {
    expect(dtypeSize(DType.Bool)).toBe(1);
    expect(dtypeSize(DType.Int8)).toBe(1);
    expect(dtypeSize(DType.Uint8)).toBe(1);
    expect(dtypeSize(DType.Int16)).toBe(2);
    expect(dtypeSize(DType.Uint16)).toBe(2);
    expect(dtypeSize(DType.Float16)).toBe(2);
    expect(dtypeSize(DType.Int32)).toBe(4);
    expect(dtypeSize(DType.Uint32)).toBe(4);
    expect(dtypeSize(DType.Float32)).toBe(4);
    expect(dtypeSize(DType.Int64)).toBe(8);
    expect(dtypeSize(DType.Uint64)).toBe(8);
    expect(dtypeSize(DType.Float64)).toBe(8);
    expect(dtypeSize(DType.Complex64)).toBe(8);
    expect(dtypeSize(DType.Complex128)).toBe(16);
  });
});

describe('Type Promotion', () => {
  describe('promoteTypes', () => {
    it('should return same type when both are equal', () => {
      expect(promoteTypes(DType.Float64, DType.Float64)).toBe(DType.Float64);
      expect(promoteTypes(DType.Int32, DType.Int32)).toBe(DType.Int32);
    });

    it('should promote bool to other types', () => {
      expect(promoteTypes(DType.Bool, DType.Int32)).toBe(DType.Int32);
      expect(promoteTypes(DType.Bool, DType.Float64)).toBe(DType.Float64);
      expect(promoteTypes(DType.Int32, DType.Bool)).toBe(DType.Int32);
    });

    it('should promote integers to floats', () => {
      expect(promoteTypes(DType.Int32, DType.Float32)).toBe(DType.Float64);
      expect(promoteTypes(DType.Int16, DType.Float32)).toBe(DType.Float32);
      expect(promoteTypes(DType.Int64, DType.Float32)).toBe(DType.Float64);
    });

    it('should promote to complex when either is complex', () => {
      expect(promoteTypes(DType.Float64, DType.Complex64)).toBe(DType.Complex128);
      expect(promoteTypes(DType.Float32, DType.Complex64)).toBe(DType.Complex64);
      expect(promoteTypes(DType.Int32, DType.Complex64)).toBe(DType.Complex64);
      expect(promoteTypes(DType.Complex64, DType.Complex128)).toBe(DType.Complex128);
    });

    it('should handle mixed signed/unsigned integers', () => {
      expect(promoteTypes(DType.Int32, DType.Uint32)).toBe(DType.Int64);
      expect(promoteTypes(DType.Int16, DType.Uint8)).toBe(DType.Int16);
      expect(promoteTypes(DType.Int64, DType.Uint64)).toBe(DType.Float64);
    });

    it('should promote same-signedness integers by size', () => {
      expect(promoteTypes(DType.Int8, DType.Int32)).toBe(DType.Int32);
      expect(promoteTypes(DType.Uint8, DType.Uint32)).toBe(DType.Uint32);
    });
  });

  describe('commonType', () => {
    it('should handle empty input', () => {
      expect(commonType()).toBe(DType.Float64);
    });

    it('should handle single type', () => {
      expect(commonType(DType.Int32)).toBe(DType.Int32);
    });

    it('should find common type for multiple types', () => {
      expect(commonType(DType.Int32, DType.Float32)).toBe(DType.Float64);
      expect(commonType(DType.Int8, DType.Int16, DType.Int32)).toBe(DType.Int32);
      expect(commonType(DType.Int32, DType.Float32, DType.Complex64)).toBe(DType.Complex128);
    });
  });
});

describe('TypedArray Conversion', () => {
  describe('dtypeToTypedArrayConstructor', () => {
    it('should return correct constructors', () => {
      expect(dtypeToTypedArrayConstructor(DType.Float32)).toBe(Float32Array);
      expect(dtypeToTypedArrayConstructor(DType.Float64)).toBe(Float64Array);
      expect(dtypeToTypedArrayConstructor(DType.Int32)).toBe(Int32Array);
      expect(dtypeToTypedArrayConstructor(DType.Int16)).toBe(Int16Array);
      expect(dtypeToTypedArrayConstructor(DType.Int8)).toBe(Int8Array);
      expect(dtypeToTypedArrayConstructor(DType.Uint32)).toBe(Uint32Array);
      expect(dtypeToTypedArrayConstructor(DType.Uint16)).toBe(Uint16Array);
      expect(dtypeToTypedArrayConstructor(DType.Uint8)).toBe(Uint8Array);
      expect(dtypeToTypedArrayConstructor(DType.Bool)).toBe(Uint8Array);
      expect(dtypeToTypedArrayConstructor(DType.Int64)).toBe(BigInt64Array);
      expect(dtypeToTypedArrayConstructor(DType.Uint64)).toBe(BigUint64Array);
    });
  });

  describe('typedArrayToDType', () => {
    it('should infer dtype from typed arrays', () => {
      expect(typedArrayToDType(new Float64Array(1))).toBe(DType.Float64);
      expect(typedArrayToDType(new Float32Array(1))).toBe(DType.Float32);
      expect(typedArrayToDType(new Int32Array(1))).toBe(DType.Int32);
      expect(typedArrayToDType(new Int16Array(1))).toBe(DType.Int16);
      expect(typedArrayToDType(new Int8Array(1))).toBe(DType.Int8);
      expect(typedArrayToDType(new Uint32Array(1))).toBe(DType.Uint32);
      expect(typedArrayToDType(new Uint16Array(1))).toBe(DType.Uint16);
      expect(typedArrayToDType(new Uint8Array(1))).toBe(DType.Uint8);
      expect(typedArrayToDType(new BigInt64Array(1))).toBe(DType.Int64);
      expect(typedArrayToDType(new BigUint64Array(1))).toBe(DType.Uint64);
    });
  });
});
