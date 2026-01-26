/**
 * NumJS I/O Operations Tests
 *
 * Tests for NPY format, text I/O, binary I/O, array printing, and base conversion.
 */

import { describe, it, expect, beforeAll } from 'vitest';
import {
  NDArray,
  DType,
  loadWasmModule,
  // NPY format
  save,
  load,
  dtypeToDescr,
  descrToDtype,
  // Text I/O
  formatValue,
  // Binary I/O
  frombuffer,
  // Array printing
  setPrintoptions,
  getPrintoptions,
  resetPrintoptions,
  withPrintoptions,
  array2string,
  arrayRepr,
  formatFloatPositional,
  formatFloatScientific,
  // Base conversion
  binaryRepr,
  baseRepr,
  fromBinaryRepr,
  fromBaseRepr,
  hexRepr,
  octalRepr,
  binaryReprArray,
  baseReprArray,
} from 'numjs-wasm';

// Initialize WASM module
beforeAll(async () => {
  await loadWasmModule();
});

describe('NPY Format', () => {
  describe('dtypeToDescr', () => {
    it('should convert Float64 to <f8', () => {
      expect(dtypeToDescr(DType.Float64)).toBe('<f8');
    });

    it('should convert Float32 to <f4', () => {
      expect(dtypeToDescr(DType.Float32)).toBe('<f4');
    });

    it('should convert Int32 to <i4', () => {
      expect(dtypeToDescr(DType.Int32)).toBe('<i4');
    });

    it('should convert Bool to |b1', () => {
      expect(dtypeToDescr(DType.Bool)).toBe('|b1');
    });

    it('should convert Uint8 to |u1', () => {
      expect(dtypeToDescr(DType.Uint8)).toBe('|u1');
    });
  });

  describe('descrToDtype', () => {
    it('should parse <f8 to Float64', () => {
      expect(descrToDtype('<f8')).toBe(DType.Float64);
    });

    it('should parse <f4 to Float32', () => {
      expect(descrToDtype('<f4')).toBe(DType.Float32);
    });

    it('should parse <i4 to Int32', () => {
      expect(descrToDtype('<i4')).toBe(DType.Int32);
    });

    it('should parse |b1 to Bool', () => {
      expect(descrToDtype('|b1')).toBe(DType.Bool);
    });

    it('should handle endianness prefix', () => {
      expect(descrToDtype('>f8')).toBe(DType.Float64);
      expect(descrToDtype('=f8')).toBe(DType.Float64);
    });
  });

  describe('save and load', () => {
    it('should round-trip 1D Float64 array', async () => {
      const arr = await NDArray.fromArray([1.5, 2.5, 3.5, 4.5]);
      const buffer = await save(null, arr);

      expect(buffer).toBeInstanceOf(ArrayBuffer);

      const loaded = await load(buffer as ArrayBuffer);
      expect(loaded.shape).toEqual([4]);
      expect(loaded.dtype).toBe(DType.Float64);
      expect(loaded.toArray()).toEqual([1.5, 2.5, 3.5, 4.5]);

      arr.dispose();
      loaded.dispose();
    });

    it('should round-trip 2D Float64 array', async () => {
      const arr = await NDArray.fromArray([[1, 2, 3], [4, 5, 6]]);
      const buffer = await save(null, arr);

      const loaded = await load(buffer as ArrayBuffer);
      expect(loaded.shape).toEqual([2, 3]);
      expect(loaded.toArray()).toEqual([1, 2, 3, 4, 5, 6]);

      arr.dispose();
      loaded.dispose();
    });

    it('should round-trip Float32 array', async () => {
      const arr = await NDArray.fromArray([1.0, 2.0, 3.0], undefined, { dtype: DType.Float32 });
      const buffer = await save(null, arr);

      const loaded = await load(buffer as ArrayBuffer);
      expect(loaded.dtype).toBe(DType.Float32);
      expect(loaded.toArray()).toEqual([1, 2, 3]);

      arr.dispose();
      loaded.dispose();
    });

    it('should round-trip Int32 array', async () => {
      const arr = await NDArray.fromArray([10, 20, 30, 40], undefined, { dtype: DType.Int32 });
      const buffer = await save(null, arr);

      const loaded = await load(buffer as ArrayBuffer);
      expect(loaded.dtype).toBe(DType.Int32);
      expect(loaded.toArray()).toEqual([10, 20, 30, 40]);

      arr.dispose();
      loaded.dispose();
    });
  });
});

describe('Text I/O', () => {
  describe('formatValue', () => {
    it('should format integer with %d', () => {
      expect(formatValue(42, '%d')).toBe('42');
      expect(formatValue(42.9, '%d')).toBe('42');
    });

    it('should format float with %f', () => {
      expect(formatValue(3.14159, '%.2f')).toBe('3.14');
      expect(formatValue(3.14159, '%.4f')).toBe('3.1416');
    });

    it('should format scientific with %e', () => {
      expect(formatValue(1234.5, '%.2e')).toBe('1.23e+3');
      expect(formatValue(0.00123, '%.2e')).toBe('1.23e-3');
    });

    it('should format with width padding', () => {
      expect(formatValue(42, '%5d')).toBe('   42');
      expect(formatValue(42, '%-5d')).toBe('42   ');
      expect(formatValue(42, '%05d')).toBe('00042');
    });

    it('should format with sign', () => {
      expect(formatValue(42, '%+d')).toBe('+42');
      expect(formatValue(-42, '%+d')).toBe('-42');
      expect(formatValue(42, '% d')).toBe(' 42');
    });
  });
});

describe('Binary I/O', () => {
  describe('frombuffer', () => {
    it('should create array from Float64 buffer', async () => {
      const floats = new Float64Array([1.5, 2.5, 3.5, 4.5]);
      const arr = await frombuffer(floats.buffer, { dtype: DType.Float64 });

      expect(arr.shape).toEqual([4]);
      expect(arr.toArray()).toEqual([1.5, 2.5, 3.5, 4.5]);

      arr.dispose();
    });

    it('should create array from Float32 buffer', async () => {
      const floats = new Float32Array([1, 2, 3, 4]);
      const arr = await frombuffer(floats.buffer, { dtype: DType.Float32 });

      expect(arr.shape).toEqual([4]);
      expect(arr.toArray()).toEqual([1, 2, 3, 4]);

      arr.dispose();
    });

    it('should respect count option', async () => {
      const floats = new Float64Array([1, 2, 3, 4, 5]);
      const arr = await frombuffer(floats.buffer, { dtype: DType.Float64, count: 3 });

      expect(arr.shape).toEqual([3]);
      expect(arr.toArray()).toEqual([1, 2, 3]);

      arr.dispose();
    });

    it('should respect offset option', async () => {
      const floats = new Float64Array([1, 2, 3, 4, 5]);
      // Offset by 2 Float64 values = 16 bytes
      const arr = await frombuffer(floats.buffer, { dtype: DType.Float64, offset: 16 });

      expect(arr.shape).toEqual([3]);
      expect(arr.toArray()).toEqual([3, 4, 5]);

      arr.dispose();
    });

    it('should handle empty result', async () => {
      const floats = new Float64Array([1, 2, 3]);
      // Offset past all data
      const arr = await frombuffer(floats.buffer, { dtype: DType.Float64, offset: 24 });

      expect(arr.shape).toEqual([0]);
      expect(arr.size).toBe(0);

      arr.dispose();
    });
  });
});

describe('Array Printing', () => {
  describe('print options', () => {
    it('should get default options', () => {
      resetPrintoptions();
      const opts = getPrintoptions();

      expect(opts.precision).toBe(8);
      expect(opts.threshold).toBe(1000);
      expect(opts.edgeitems).toBe(3);
    });

    it('should set options', () => {
      resetPrintoptions();
      setPrintoptions({ precision: 4 });

      const opts = getPrintoptions();
      expect(opts.precision).toBe(4);

      resetPrintoptions();
    });

    it('should restore options with withPrintoptions', () => {
      resetPrintoptions();
      setPrintoptions({ precision: 8 });

      withPrintoptions({ precision: 2 }, () => {
        const opts = getPrintoptions();
        expect(opts.precision).toBe(2);
      });

      const opts = getPrintoptions();
      expect(opts.precision).toBe(8);
    });
  });

  describe('array2string', () => {
    it('should format 1D array', async () => {
      const arr = await NDArray.fromArray([1, 2, 3]);
      const str = array2string(arr);

      expect(str).toContain('1');
      expect(str).toContain('2');
      expect(str).toContain('3');
      expect(str.startsWith('[')).toBe(true);
      expect(str.endsWith(']')).toBe(true);

      arr.dispose();
    });

    it('should format 2D array', async () => {
      const arr = await NDArray.fromArray([[1, 2], [3, 4]]);
      const str = array2string(arr);

      expect(str).toContain('[');
      expect(str).toContain(']');

      arr.dispose();
    });
  });

  describe('arrayRepr', () => {
    it('should include array() wrapper', async () => {
      const arr = await NDArray.fromArray([1, 2, 3]);
      const repr = arrayRepr(arr);

      expect(repr.startsWith('array(')).toBe(true);
      expect(repr.endsWith(')')).toBe(true);

      arr.dispose();
    });

    it('should include dtype for non-float64', async () => {
      const arr = await NDArray.fromArray([1, 2, 3], undefined, { dtype: DType.Int32 });
      const repr = arrayRepr(arr);

      expect(repr).toContain('dtype=int32');

      arr.dispose();
    });
  });

  describe('formatFloatPositional', () => {
    it('should format with precision', () => {
      expect(formatFloatPositional(3.14159, 2)).toBe('3.14');
      expect(formatFloatPositional(3.14159, 4)).toBe('3.1416');
    });

    it('should trim trailing zeros', () => {
      expect(formatFloatPositional(3.10000, 5, true, true, '0')).toBe('3.1');
    });

    it('should add sign', () => {
      expect(formatFloatPositional(3.14, 2, true, true, 'k', true)).toBe('+3.14');
    });
  });

  describe('formatFloatScientific', () => {
    it('should format in scientific notation', () => {
      const result = formatFloatScientific(1234.5, 2);
      expect(result).toMatch(/1\.23e\+0*3/);
    });

    it('should handle small numbers', () => {
      const result = formatFloatScientific(0.00123, 2);
      expect(result).toMatch(/1\.23e-0*3/);
    });
  });
});

describe('Base Conversion', () => {
  describe('binaryRepr', () => {
    it('should convert positive integers', () => {
      expect(binaryRepr(5)).toBe('101');
      expect(binaryRepr(8)).toBe('1000');
      expect(binaryRepr(255)).toBe('11111111');
    });

    it('should pad with width', () => {
      expect(binaryRepr(5, 8)).toBe('00000101');
      expect(binaryRepr(1, 4)).toBe('0001');
    });

    it('should handle negative without width', () => {
      expect(binaryRepr(-5)).toBe('-101');
    });

    it('should use twos complement with width', () => {
      expect(binaryRepr(-1, 8)).toBe('11111111');
      expect(binaryRepr(-5, 8)).toBe('11111011');
      expect(binaryRepr(-128, 8)).toBe('10000000');
    });

    it('should throw on width too small', () => {
      expect(() => binaryRepr(256, 8)).toThrow();
      expect(() => binaryRepr(-129, 8)).toThrow();
    });
  });

  describe('baseRepr', () => {
    it('should convert to various bases', () => {
      expect(baseRepr(10, 2)).toBe('1010');
      expect(baseRepr(10, 8)).toBe('12');
      expect(baseRepr(10, 16)).toBe('A');
      expect(baseRepr(255, 16)).toBe('FF');
    });

    it('should handle negative numbers', () => {
      expect(baseRepr(-10, 16)).toBe('-A');
    });

    it('should pad with zeros', () => {
      expect(baseRepr(10, 2, 8)).toBe('00001010');
      expect(baseRepr(15, 16, 4)).toBe('000F');
    });

    it('should throw on invalid base', () => {
      expect(() => baseRepr(10, 1)).toThrow();
      expect(() => baseRepr(10, 37)).toThrow();
    });
  });

  describe('fromBinaryRepr', () => {
    it('should parse binary strings', () => {
      expect(fromBinaryRepr('101')).toBe(5);
      expect(fromBinaryRepr('1000')).toBe(8);
      expect(fromBinaryRepr('0b101')).toBe(5);
    });

    it('should handle negative sign', () => {
      expect(fromBinaryRepr('-101')).toBe(-5);
    });

    it('should parse twos complement', () => {
      expect(fromBinaryRepr('11111111', true, 8)).toBe(-1);
      expect(fromBinaryRepr('11111011', true, 8)).toBe(-5);
    });
  });

  describe('fromBaseRepr', () => {
    it('should parse base-N strings', () => {
      expect(fromBaseRepr('1010', 2)).toBe(10);
      expect(fromBaseRepr('12', 8)).toBe(10);
      expect(fromBaseRepr('A', 16)).toBe(10);
      expect(fromBaseRepr('FF', 16)).toBe(255);
    });

    it('should handle negative', () => {
      expect(fromBaseRepr('-A', 16)).toBe(-10);
    });
  });

  describe('hexRepr', () => {
    it('should convert to hex', () => {
      expect(hexRepr(255)).toBe('FF');
      expect(hexRepr(16)).toBe('10');
    });

    it('should pad with zeros', () => {
      expect(hexRepr(255, 4)).toBe('00FF');
    });
  });

  describe('octalRepr', () => {
    it('should convert to octal', () => {
      expect(octalRepr(8)).toBe('10');
      expect(octalRepr(64)).toBe('100');
    });
  });

  describe('binaryReprArray', () => {
    it('should convert array of numbers', () => {
      expect(binaryReprArray([1, 2, 3, 4])).toEqual(['1', '10', '11', '100']);
      expect(binaryReprArray([1, 2, 3, 4], 4)).toEqual(['0001', '0010', '0011', '0100']);
    });
  });

  describe('baseReprArray', () => {
    it('should convert array of numbers', () => {
      expect(baseReprArray([10, 20, 30], 16)).toEqual(['A', '14', '1E']);
    });
  });
});
