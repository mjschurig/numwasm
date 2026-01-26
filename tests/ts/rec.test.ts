/**
 * NumJS Record Arrays Module Tests
 *
 * Tests for numpy.rec compatible record arrays.
 */

import { describe, it, expect } from 'vitest';
import {
  rec,
  recarray,
  record,
  format_parser,
  find_duplicate,
  NDArray,
  DType,
} from '../../dist/numjs.mjs';
import type { StructuredDType } from '../../dist/numjs.mjs';

/* ============ format_parser Tests ============ */

describe('format_parser', () => {
  describe('format string parsing', () => {
    it('parses comma-separated format string', () => {
      const parser = new format_parser('f8, i4, S10');
      expect(parser.dtype.names).toEqual(['f0', 'f1', 'f2']);
      expect(parser.dtype.fieldList.length).toBe(3);
    });

    it('parses format list', () => {
      const parser = new format_parser(['f8', 'i4', 'S10'], ['x', 'y', 'name']);
      expect(parser.dtype.names).toEqual(['x', 'y', 'name']);
    });

    it('parses format with explicit names', () => {
      const parser = new format_parser('f8, i4', 'x, y');
      expect(parser.dtype.names).toEqual(['x', 'y']);
    });

    it('generates default names when not provided', () => {
      const parser = new format_parser(['i4', 'f8', 'U10']);
      expect(parser.dtype.names).toEqual(['f0', 'f1', 'f2']);
    });
  });

  describe('dtype interpretation', () => {
    it('parses integer formats correctly', () => {
      const parser = new format_parser('i1, i2, i4, i8');
      const fields = parser.dtype.fieldList;

      expect(fields[0].dtype).toBe(DType.Int8);
      expect(fields[1].dtype).toBe(DType.Int16);
      expect(fields[2].dtype).toBe(DType.Int32);
      expect(fields[3].dtype).toBe(DType.Int64);
    });

    it('parses unsigned integer formats', () => {
      const parser = new format_parser('u1, u2, u4');
      const fields = parser.dtype.fieldList;

      expect(fields[0].dtype).toBe(DType.Uint8);
      expect(fields[1].dtype).toBe(DType.Uint16);
      expect(fields[2].dtype).toBe(DType.Uint32);
    });

    it('parses float formats', () => {
      const parser = new format_parser('f4, f8');
      const fields = parser.dtype.fieldList;

      expect(fields[0].dtype).toBe(DType.Float32);
      expect(fields[1].dtype).toBe(DType.Float64);
    });

    it('parses string formats', () => {
      const parser = new format_parser('S10, U20');
      const fields = parser.dtype.fieldList;

      expect(fields[0].dtype).toBe(DType.String);
      expect(fields[0].charType).toBe('S');
      expect(fields[0].itemsize).toBe(10);

      expect(fields[1].dtype).toBe(DType.String);
      expect(fields[1].charType).toBe('U');
      expect(fields[1].itemsize).toBe(80); // 20 * 4 bytes per Unicode char
    });

    it('parses boolean format', () => {
      const parser = new format_parser('?');
      expect(parser.dtype.fieldList[0].dtype).toBe(DType.Bool);
    });
  });

  describe('offset calculation', () => {
    it('calculates packed offsets', () => {
      const parser = new format_parser('i4, f8, i2', null, null, false);
      const fields = parser.dtype.fieldList;

      expect(fields[0].offset).toBe(0);
      expect(fields[1].offset).toBe(4);
      expect(fields[2].offset).toBe(12);
      expect(parser.dtype.itemsize).toBe(14);
    });

    it('calculates aligned offsets', () => {
      const parser = new format_parser('i1, i4', null, null, true);
      const fields = parser.dtype.fieldList;

      expect(fields[0].offset).toBe(0);
      expect(fields[1].offset).toBe(4); // Padded for 4-byte alignment
      expect(parser.dtype.itemsize).toBe(8);
    });
  });

  describe('error handling', () => {
    it('throws on duplicate field names', () => {
      expect(() => new format_parser('f8, i4', 'x, x')).toThrow();
    });

    it('throws on invalid format string', () => {
      expect(() => new format_parser('invalid')).toThrow();
    });

    it('throws on mismatched names count', () => {
      expect(() => new format_parser('f8, i4, f8', 'x, y')).toThrow();
    });

    it('throws on empty formats', () => {
      expect(() => new format_parser('')).toThrow();
    });

    it('throws on invalid field names', () => {
      expect(() => new format_parser('f8', '123invalid')).toThrow();
    });
  });

  describe('titles', () => {
    it('parses titles', () => {
      const parser = new format_parser('f8, i4', 'x, y', 'X Value, Y Value');
      const fields = parser.dtype.fieldList;

      expect(fields[0].title).toBe('X Value');
      expect(fields[1].title).toBe('Y Value');
    });

    it('creates title aliases in fields map', () => {
      const parser = new format_parser('f8', 'x', 'X Value');

      expect(parser.dtype.fields.has('x')).toBe(true);
      expect(parser.dtype.fields.has('X Value')).toBe(true);
    });
  });
});

/* ============ find_duplicate Tests ============ */

describe('find_duplicate', () => {
  it('finds duplicate strings', () => {
    const dupes = find_duplicate(['a', 'b', 'a', 'c', 'b']);
    expect(dupes.sort()).toEqual(['a', 'b'].sort());
  });

  it('returns empty array for no duplicates', () => {
    const dupes = find_duplicate(['a', 'b', 'c']);
    expect(dupes).toEqual([]);
  });

  it('handles empty array', () => {
    const dupes = find_duplicate([]);
    expect(dupes).toEqual([]);
  });
});

/* ============ record Tests ============ */

describe('record', () => {
  it('provides field access via methods', async () => {
    const arr = await rec.fromrecords(
      [[1, 'Alice'], [2, 'Bob']],
      { formats: ['i4', 'U10'], names: ['age', 'name'] }
    );

    const r = arr.getRecord(0);
    expect(r.getField('age')).toBe(1);
    expect(r.getField('name')).toBe('Alice');
  });

  it('provides field access via attributes', async () => {
    const arr = await rec.fromrecords(
      [[25, 'Alice']],
      { formats: ['i4', 'U10'], names: ['age', 'name'] }
    );

    const r = arr.getRecord(0);
    // Proxy enables attribute access
    expect((r as any).age).toBe(25);
    expect((r as any).name).toBe('Alice');
  });

  it('supports field modification', async () => {
    const arr = await rec.fromrecords(
      [[1, 'Alice']],
      { formats: ['i4', 'U10'], names: ['age', 'name'] }
    );

    const r = arr.getRecord(0);
    r.setField('age', 30);
    expect(r.getField('age')).toBe(30);
  });

  it('converts to array', async () => {
    const arr = await rec.fromrecords(
      [[1, 'Alice']],
      { formats: ['i4', 'U10'], names: ['age', 'name'] }
    );

    expect(arr.getRecord(0).toArray()).toEqual([1, 'Alice']);
  });

  it('converts to object', async () => {
    const arr = await rec.fromrecords(
      [[1, 'Alice']],
      { formats: ['i4', 'U10'], names: ['age', 'name'] }
    );

    expect(arr.getRecord(0).toObject()).toEqual({ age: 1, name: 'Alice' });
  });

  it('has correct string representation', async () => {
    const arr = await rec.fromrecords(
      [[25, 'Alice']],
      { formats: ['i4', 'U10'], names: ['age', 'name'] }
    );

    const str = arr.getRecord(0).toString();
    expect(str).toContain('age=25');
    expect(str).toContain("name='Alice'");
  });

  it('supports iteration over field values', async () => {
    const arr = await rec.fromrecords(
      [[1, 'Alice']],
      { formats: ['i4', 'U10'], names: ['age', 'name'] }
    );

    const values = [...arr.getRecord(0)];
    expect(values).toEqual([1, 'Alice']);
  });

  it('exposes dtype and names', async () => {
    const arr = await rec.fromrecords(
      [[1, 'Alice']],
      { formats: ['i4', 'U10'], names: ['age', 'name'] }
    );

    const r = arr.getRecord(0);
    expect(r.names).toEqual(['age', 'name']);
    expect(r.dtype.names).toEqual(['age', 'name']);
    expect(r.length).toBe(2);
  });

  it('throws on invalid field access', async () => {
    const arr = await rec.fromrecords(
      [[1, 'Alice']],
      { formats: ['i4', 'U10'], names: ['age', 'name'] }
    );

    const r = arr.getRecord(0);
    expect(() => r.getField('invalid')).toThrow();
  });
});

/* ============ recarray Tests ============ */

describe('recarray', () => {
  describe('creation', () => {
    it('creates from fromarrays', async () => {
      const arr = await rec.fromarrays(
        [[1, 2, 3], ['a', 'b', 'c']],
        { names: ['x', 'y'] }
      );

      expect(arr.size).toBe(3);
      expect(arr.names).toEqual(['x', 'y']);
    });

    it('creates from fromrecords', async () => {
      const arr = await rec.fromrecords(
        [[1, 'a'], [2, 'b'], [3, 'c']],
        { formats: ['i4', 'U1'], names: ['x', 'y'] }
      );

      expect(arr.size).toBe(3);
    });

    it('creates from array function', async () => {
      const arr = await rec.array(
        [[1, 'a'], [2, 'b']],
        { formats: ['i4', 'U1'], names: ['x', 'y'] }
      );

      expect(arr.size).toBe(2);
    });

    it('infers formats when not provided', async () => {
      const arr = await rec.fromarrays(
        [[1.5, 2.5, 3.5], [1, 2, 3], ['a', 'bb', 'ccc']],
        { names: ['floats', 'ints', 'strings'] }
      );

      expect(arr.dtype.fieldList[0].dtype).toBe(DType.Float64);
      expect(arr.dtype.fieldList[1].dtype).toBe(DType.Int32);
      expect(arr.dtype.fieldList[2].dtype).toBe(DType.String);
    });
  });

  describe('field access', () => {
    it('returns field as NDArray', async () => {
      const arr = await rec.fromarrays(
        [[1, 2, 3], ['a', 'b', 'c']],
        { names: ['x', 'y'] }
      );

      const xField = arr.field('x');
      expect(xField).toBeInstanceOf(NDArray);
      expect(xField.toArray()).toEqual([1, 2, 3]);
    });

    it('returns string field as NDArray', async () => {
      const arr = await rec.fromarrays(
        [[1, 2], ['Alice', 'Bob']],
        { names: ['id', 'name'] }
      );

      const nameField = arr.field('name');
      expect(nameField.isStringArray).toBe(true);
      expect(nameField.toArray()).toEqual(['Alice', 'Bob']);
    });

    it('supports field access by index', async () => {
      const arr = await rec.fromarrays(
        [[1, 2], ['a', 'b']],
        { names: ['x', 'y'] }
      );

      const field = arr.field(0);
      expect(field.toArray()).toEqual([1, 2]);
    });

    it('throws on invalid field name', async () => {
      const arr = await rec.fromarrays(
        [[1, 2]],
        { names: ['x'] }
      );

      expect(() => arr.field('invalid')).toThrow();
    });

    it('supports attribute-style field access via proxy', async () => {
      const arr = await rec.fromarrays(
        [[1, 2, 3]],
        { names: ['values'] }
      );

      // Proxy enables this
      const values = (arr as any).values;
      expect(values.toArray()).toEqual([1, 2, 3]);
    });
  });

  describe('record access', () => {
    it('returns record at index', async () => {
      const arr = await rec.fromarrays(
        [[1, 2, 3], ['a', 'b', 'c']],
        { names: ['x', 'y'] }
      );

      const r = arr.getRecord(1);
      expect(r.getField('x')).toBe(2);
      expect(r.getField('y')).toBe('b');
    });

    it('supports negative indexing', async () => {
      const arr = await rec.fromarrays(
        [[1, 2, 3]],
        { names: ['x'] }
      );

      expect(arr.getRecord(-1).getField('x')).toBe(3);
      expect(arr.getRecord(-2).getField('x')).toBe(2);
    });

    it('throws on out-of-bounds index', async () => {
      const arr = await rec.fromarrays(
        [[1, 2, 3]],
        { names: ['x'] }
      );

      expect(() => arr.getRecord(10)).toThrow();
      expect(() => arr.getRecord(-10)).toThrow();
    });

    it('setRecord updates values by array', async () => {
      const arr = await rec.fromarrays(
        [[1, 2], ['a', 'b']],
        { names: ['x', 'y'] }
      );

      arr.setRecord(0, [99, 'z']);
      expect(arr.getRecord(0).toArray()).toEqual([99, 'z']);
    });

    it('setRecord updates values by object', async () => {
      const arr = await rec.fromarrays(
        [[1, 2], ['a', 'b']],
        { names: ['x', 'y'] }
      );

      arr.setRecord(1, { x: 42, y: 'test' });
      expect(arr.getRecord(1).toObject()).toEqual({ x: 42, y: 'test' });
    });
  });

  describe('iteration', () => {
    it('iterates over records', async () => {
      const arr = await rec.fromarrays(
        [[1, 2, 3], ['a', 'b', 'c']],
        { names: ['x', 'y'] }
      );

      const records = [...arr];
      expect(records.length).toBe(3);
      expect(records[0].getField('x')).toBe(1);
      expect(records[2].getField('y')).toBe('c');
    });
  });

  describe('conversion', () => {
    it('tolist returns array of tuples', async () => {
      const arr = await rec.fromarrays(
        [[1, 2], ['a', 'b']],
        { names: ['x', 'y'] }
      );

      expect(arr.tolist()).toEqual([[1, 'a'], [2, 'b']]);
    });

    it('toObjects returns array of objects', async () => {
      const arr = await rec.fromarrays(
        [[1, 2], ['a', 'b']],
        { names: ['x', 'y'] }
      );

      expect(arr.toObjects()).toEqual([
        { x: 1, y: 'a' },
        { x: 2, y: 'b' }
      ]);
    });

    it('toString produces readable output', async () => {
      const arr = await rec.fromarrays(
        [[1, 2]],
        { names: ['x'] }
      );

      const str = arr.toString();
      expect(str).toContain('recarray');
      expect(str).toContain('x=1');
    });
  });

  describe('properties', () => {
    it('has correct shape', async () => {
      const arr = await rec.fromarrays(
        [[1, 2, 3, 4, 5]],
        { names: ['x'] }
      );

      expect(arr.shape).toEqual([5]);
      expect(arr.ndim).toBe(1);
      expect(arr.size).toBe(5);
    });

    it('exposes dtype', async () => {
      const arr = await rec.fromarrays(
        [[1, 2], ['a', 'b']],
        { formats: ['i4', 'U10'], names: ['x', 'y'] }
      );

      expect(arr.dtype.names).toEqual(['x', 'y']);
      expect(arr.dtype.fieldList[0].dtype).toBe(DType.Int32);
    });
  });

  describe('memory management', () => {
    it('copy creates independent array', async () => {
      const arr = await rec.fromarrays(
        [[1, 2, 3]],
        { names: ['x'] }
      );

      const copy = arr.copy();
      arr.setRecord(0, [99]);

      expect(arr.getRecord(0).getField('x')).toBe(99);
      expect(copy.getRecord(0).getField('x')).toBe(1);

      copy.dispose();
    });

    it('dispose frees resources', async () => {
      const arr = await rec.fromarrays(
        [[1, 2, 3]],
        { names: ['x'] }
      );

      arr.dispose();
      expect(arr.isDisposed).toBe(true);
      expect(() => arr.field('x')).toThrow();
    });
  });
});

/* ============ Creation Functions Tests ============ */

describe('Creation Functions', () => {
  describe('rec.array', () => {
    it('creates from null (empty)', async () => {
      const arr = await rec.array(null, {
        formats: ['i4', 'f8'],
        names: ['x', 'y']
      });

      expect(arr.size).toBe(0);
    });

    it('creates from list of tuples', async () => {
      const arr = await rec.array(
        [[1, 'a'], [2, 'b'], [3, 'c']],
        { formats: ['i4', 'U1'], names: ['x', 'y'] }
      );

      expect(arr.size).toBe(3);
      expect(arr.getRecord(0).toArray()).toEqual([1, 'a']);
    });

    it('infers formats from data', async () => {
      const arr = await rec.array(
        [[1, 'hello'], [2, 'world']],
        { names: ['num', 'text'] }
      );

      expect(arr.dtype.fieldList[0].dtype).toBe(DType.Int32);
      expect(arr.dtype.fieldList[1].dtype).toBe(DType.String);
    });

    it('copies existing recarray when copy=true', async () => {
      const original = await rec.fromarrays(
        [[1, 2]],
        { names: ['x'] }
      );

      const copied = await rec.array(original, { copy: true });
      original.setRecord(0, [99]);

      expect(copied.getRecord(0).getField('x')).toBe(1);

      original.dispose();
      copied.dispose();
    });
  });

  describe('rec.fromarrays', () => {
    it('creates from column arrays', async () => {
      const ages = [25, 30, 35];
      const names = ['Alice', 'Bob', 'Carol'];

      const arr = await rec.fromarrays([ages, names], {
        names: ['age', 'name']
      });

      expect(arr.field('age').toArray()).toEqual([25, 30, 35]);
      expect(arr.field('name').toArray()).toEqual(['Alice', 'Bob', 'Carol']);
    });

    it('creates from NDArray columns', async () => {
      const xArr = await NDArray.fromArray([1, 2, 3]);
      const yArr = await NDArray.fromArray([4, 5, 6]);

      const arr = await rec.fromarrays([xArr, yArr], {
        names: ['x', 'y']
      });

      expect(arr.field('x').toArray()).toEqual([1, 2, 3]);

      xArr.dispose();
      yArr.dispose();
    });

    it('validates array lengths match', async () => {
      await expect(rec.fromarrays(
        [[1, 2, 3], ['a', 'b']],
        { names: ['x', 'y'] }
      )).rejects.toThrow();
    });

    it('validates non-empty arrayList', async () => {
      await expect(rec.fromarrays([], { names: [] })).rejects.toThrow();
    });
  });

  describe('rec.fromrecords', () => {
    it('creates from row records', async () => {
      const arr = await rec.fromrecords(
        [[25, 'Alice'], [30, 'Bob']],
        { formats: ['i4', 'U10'], names: ['age', 'name'] }
      );

      expect(arr.getRecord(0).toArray()).toEqual([25, 'Alice']);
      expect(arr.getRecord(1).toArray()).toEqual([30, 'Bob']);
    });
  });

  describe('rec.fromstring', () => {
    it('parses binary data', async () => {
      const parser = new format_parser('i4, f8', 'x, y');
      const dtype = parser.dtype;

      // Create binary buffer: int32 + float64 = 12 bytes per record
      const buffer = new ArrayBuffer(24); // 2 records
      const view = new DataView(buffer);

      // Record 0
      view.setInt32(0, 42, true);
      view.setFloat64(4, 3.14, true);

      // Record 1
      view.setInt32(12, 100, true);
      view.setFloat64(16, 2.71, true);

      const arr = await rec.fromstring(buffer, dtype);

      expect(arr.size).toBe(2);
      expect(arr.getRecord(0).getField('x')).toBe(42);
      expect(arr.getRecord(0).getField('y')).toBeCloseTo(3.14, 5);
      expect(arr.getRecord(1).getField('x')).toBe(100);
    });

    it('respects offset parameter', async () => {
      const parser = new format_parser('i4', 'x');
      const dtype = parser.dtype;

      const buffer = new ArrayBuffer(12);
      const view = new DataView(buffer);
      view.setInt32(0, 1, true);
      view.setInt32(4, 2, true);
      view.setInt32(8, 3, true);

      const arr = await rec.fromstring(buffer, dtype, [], 4); // Skip first int

      expect(arr.size).toBe(2);
      expect(arr.getRecord(0).getField('x')).toBe(2);
      expect(arr.getRecord(1).getField('x')).toBe(3);
    });

    it('handles Uint8Array input', async () => {
      const parser = new format_parser('i4', 'x');
      const dtype = parser.dtype;

      const buffer = new ArrayBuffer(8);
      const view = new DataView(buffer);
      view.setInt32(0, 10, true);
      view.setInt32(4, 20, true);

      const uint8 = new Uint8Array(buffer);
      const arr = await rec.fromstring(uint8, dtype);

      expect(arr.size).toBe(2);
      expect(arr.getRecord(0).getField('x')).toBe(10);
    });
  });
});

/* ============ Integration Tests ============ */

describe('Integration', () => {
  it('round-trip: create -> modify -> read', async () => {
    // Create
    const arr = await rec.fromarrays(
      [[1, 2, 3], ['a', 'b', 'c']],
      { formats: ['i4', 'U10'], names: ['id', 'letter'] }
    );

    // Modify via record
    arr.getRecord(1).setField('id', 99);
    arr.getRecord(1).setField('letter', 'X');

    // Read back
    expect(arr.getRecord(1).getField('id')).toBe(99);
    expect(arr.getRecord(1).getField('letter')).toBe('X');

    // Other records unchanged
    expect(arr.getRecord(0).toArray()).toEqual([1, 'a']);
    expect(arr.getRecord(2).toArray()).toEqual([3, 'c']);

    arr.dispose();
  });

  it('mixed numeric types', async () => {
    const arr = await rec.fromarrays(
      [
        [1, 2, 3],           // i4
        [1.5, 2.5, 3.5],     // f8
        [true, false, true]  // bool
      ],
      { formats: ['i4', 'f8', '?'], names: ['int', 'float', 'bool'] }
    );

    expect(arr.getRecord(0).toObject()).toEqual({
      int: 1,
      float: 1.5,
      bool: 1  // Stored as 0/1 in WASM
    });

    expect(arr.getRecord(1).getField('bool')).toBe(0);  // false is 0

    arr.dispose();
  });

  it('large array performance', async () => {
    const n = 10000;
    const ids = Array.from({ length: n }, (_, i) => i);
    const values = Array.from({ length: n }, (_, i) => i * 0.1);

    const arr = await rec.fromarrays([ids, values], {
      formats: ['i4', 'f8'],
      names: ['id', 'value']
    });

    expect(arr.size).toBe(n);
    expect(arr.getRecord(9999).getField('id')).toBe(9999);
    expect(arr.getRecord(9999).getField('value')).toBeCloseTo(999.9, 4);

    arr.dispose();
  });
});

/* ============ rec Namespace Tests ============ */

describe('rec namespace', () => {
  it('exports all expected members', () => {
    expect(rec.array).toBeDefined();
    expect(rec.fromarrays).toBeDefined();
    expect(rec.fromrecords).toBeDefined();
    expect(rec.fromstring).toBeDefined();
    expect(rec.fromfile).toBeDefined();
    expect(rec.format_parser).toBeDefined();
    expect(rec.recarray).toBeDefined();
    expect(rec.record).toBeDefined();
    expect(rec.find_duplicate).toBeDefined();
    expect(rec.ValueError).toBeDefined();
    expect(rec.KeyError).toBeDefined();
    expect(rec.IndexError).toBeDefined();
  });
});
