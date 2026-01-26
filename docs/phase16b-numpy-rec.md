# Phase 16b: numpy.rec Implementation Plan

Complete implementation roadmap for the NumJS-WASM record arrays module, providing NumPy-compatible structured arrays with named field access.

---

## ⚠️ Implementation Guideline

**The original NumPy source code is available in `/numpy`.**

Key reference files:
- `numpy/_core/records.py` - Main implementation (1,086 lines)
- `numpy/rec/__init__.py` - Public exports
- `numpy/_core/numerictypes.py` - Type system integration

Implementation should follow NumPy's algorithms and error handling for consistency.

---

## Current State (Pre-Phase 16b)

```
src/ts/
├── NDArray.ts         # Core array class
├── types.ts           # Type definitions
├── dtype.ts           # Type utilities
└── index.ts           # Public exports
```

**Required Infrastructure:**
- NDArray core
- DType system with structured dtype support (needs to be added)

---

## Phase 16b Dependency Tree

```
PHASE 16b: NUMPY.REC
│
├── 16b.1 Structured DType Infrastructure
│   ├── StructuredDType interface
│   ├── FieldDescriptor interface
│   ├── Field offset calculation
│   └── Structured array memory layout
│
├── 16b.2 Format Parser
│   ├── format_parser class
│   │   ├── Parse format strings ('f8, i4, S10')
│   │   ├── Parse format lists (['f8', 'i4', 'S10'])
│   │   ├── Parse field names
│   │   ├── Parse field titles
│   │   ├── Handle alignment
│   │   └── Handle byte order
│   └── dtype property
│
├── 16b.3 Record Class
│   ├── Single record access
│   ├── Field access via attributes
│   ├── Field access via indexing
│   └── String representation
│
├── 16b.4 RecArray Class
│   ├── Array of records
│   ├── Field access as attributes
│   ├── field(name) method
│   ├── tolist() conversion
│   └── Iterator support
│
└── 16b.5 Creation Functions
    ├── array(data, dtype, formats, names)
    ├── fromarrays(arrayList, names)
    ├── fromrecords(recList, names)
    ├── fromstring(datastring, dtype)
    └── fromfile(fd, dtype)

Dependencies: NDArray core, DType system
```

---

## Detailed Implementation Specifications

### 16b.1 Structured DType Infrastructure

**File:** `src/ts/types.ts` (additions)

```typescript
/**
 * Descriptor for a single field in a structured dtype.
 */
export interface FieldDescriptor {
  /** Field name */
  name: string;

  /** Field data type */
  dtype: DType;

  /** Byte offset from start of record */
  offset: number;

  /** Optional display title */
  title?: string | null;

  /** Size in bytes (for fixed-width strings) */
  itemsize?: number;
}

/**
 * Structured data type with named fields.
 */
export interface StructuredDType {
  /** Ordered list of field names */
  names: string[];

  /** Field descriptors */
  fields: FieldDescriptor[];

  /** Total size of one record in bytes */
  itemsize: number;

  /** Alignment requirements */
  alignment?: number;
}

/**
 * Type guard for structured dtypes.
 */
export function isStructuredDType(dtype: any): dtype is StructuredDType {
  return typeof dtype === 'object' &&
         dtype !== null &&
         'names' in dtype &&
         'fields' in dtype &&
         'itemsize' in dtype;
}

/**
 * Get byte size for a given dtype.
 */
export function dtypeSize(dtype: DType): number {
  const sizes: Record<DType, number> = {
    [DType.Bool]: 1,
    [DType.Int8]: 1,
    [DType.Int16]: 2,
    [DType.Int32]: 4,
    [DType.Int64]: 8,
    [DType.UInt8]: 1,
    [DType.UInt16]: 2,
    [DType.UInt32]: 4,
    [DType.UInt64]: 8,
    [DType.Float16]: 2,
    [DType.Float32]: 4,
    [DType.Float64]: 8,
    [DType.Complex64]: 8,
    [DType.Complex128]: 16,
    [DType.String]: 0, // Variable
  };
  return sizes[dtype] ?? 0;
}

/**
 * Get alignment requirement for a dtype.
 */
export function dtypeAlignment(dtype: DType): number {
  // Most types align to their size, up to 8 bytes
  return Math.min(dtypeSize(dtype), 8);
}
```

---

### 16b.2 Format Parser

**File:** `src/ts/rec/format_parser.ts`

```typescript
import { DType, StructuredDType, FieldDescriptor, dtypeSize, dtypeAlignment } from '../types.js';

/**
 * Class to convert format strings, names, and titles to a structured dtype.
 *
 * Parses format specifications like 'f8, i4, S10' into proper field descriptors.
 *
 * @example
 * const parser = new format_parser('f8, i4, S10', 'x, y, name');
 * console.log(parser.dtype);
 * // { names: ['x', 'y', 'name'], fields: [...], itemsize: 22 }
 */
export class format_parser {
  private _dtype: StructuredDType;

  /**
   * Create a format parser.
   *
   * @param formats - Format specification (string like 'f8, i4' or list like ['f8', 'i4'])
   * @param names - Field names (string like 'x, y' or list like ['x', 'y'])
   * @param titles - Optional field titles for display
   * @param aligned - If true, pad fields for C struct alignment
   * @param byteorder - Byte order: '<' little, '>' big, '=' native
   */
  constructor(
    formats: string | string[],
    names: string | string[] | null = null,
    titles: string | string[] | null = null,
    aligned: boolean = false,
    byteorder: string | null = null
  ) {
    // Parse format specification
    const formatList = this._parseFormats(formats);

    // Parse names (generate defaults if not provided)
    const nameList = this._parseNames(names, formatList.length);

    // Parse titles
    const titleList = this._parseTitles(titles, formatList.length);

    // Check for duplicate names
    this._checkDuplicates(nameList);

    // Create the structured dtype
    this._dtype = this._createDtype(formatList, nameList, titleList, aligned, byteorder);
  }

  /**
   * The resulting structured dtype.
   */
  get dtype(): StructuredDType {
    return this._dtype;
  }

  /**
   * Parse format specification into list of format strings.
   */
  private _parseFormats(formats: string | string[]): string[] {
    let formatList: string[];

    if (typeof formats === 'string') {
      // Handle comma-separated string: 'f8, i4, S10'
      formatList = formats.split(',').map(f => f.trim()).filter(f => f.length > 0);
    } else {
      formatList = [...formats];
    }

    // Validate each format
    for (const fmt of formatList) {
      if (!this._isValidFormat(fmt)) {
        throw new TypeError(`Invalid format specification: '${fmt}'`);
      }
    }

    if (formatList.length === 0) {
      throw new ValueError('At least one format must be specified');
    }

    return formatList;
  }

  /**
   * Parse names into list of field names.
   */
  private _parseNames(names: string | string[] | null, count: number): string[] {
    if (names === null) {
      // Generate default names: f0, f1, f2, ...
      return Array.from({ length: count }, (_, i) => `f${i}`);
    }

    let nameList: string[];

    if (typeof names === 'string') {
      // Handle comma-separated string: 'x, y, z'
      nameList = names.split(',').map(n => n.trim()).filter(n => n.length > 0);
    } else {
      nameList = [...names];
    }

    if (nameList.length !== count) {
      throw new ValueError(
        `Mismatch between number of formats (${count}) and names (${nameList.length})`
      );
    }

    // Validate names
    for (const name of nameList) {
      if (!/^[a-zA-Z_][a-zA-Z0-9_]*$/.test(name)) {
        throw new ValueError(`Invalid field name: '${name}'`);
      }
    }

    return nameList;
  }

  /**
   * Parse titles into list of field titles.
   */
  private _parseTitles(titles: string | string[] | null, count: number): (string | null)[] {
    if (titles === null) {
      return Array(count).fill(null);
    }

    let titleList: (string | null)[];

    if (typeof titles === 'string') {
      titleList = titles.split(',').map(t => {
        const trimmed = t.trim();
        return trimmed.length > 0 ? trimmed : null;
      });
    } else {
      titleList = titles.map(t => t || null);
    }

    // Pad with nulls if needed
    while (titleList.length < count) {
      titleList.push(null);
    }

    if (titleList.length > count) {
      titleList = titleList.slice(0, count);
    }

    return titleList;
  }

  /**
   * Check for duplicate field names.
   */
  private _checkDuplicates(names: string[]): void {
    const seen = new Set<string>();
    const duplicates: string[] = [];

    for (const name of names) {
      if (seen.has(name)) {
        if (!duplicates.includes(name)) {
          duplicates.push(name);
        }
      }
      seen.add(name);
    }

    if (duplicates.length > 0) {
      throw new ValueError(`Duplicate field names: ${duplicates.join(', ')}`);
    }
  }

  /**
   * Validate a format string.
   *
   * Valid formats:
   * - Type codes: b, i, u, f, c, S, U, V, O, m, M
   * - With size: i4, f8, S10, U20
   * - With byte order: <f8, >i4, =f4
   */
  private _isValidFormat(fmt: string): boolean {
    // Pattern: optional byteorder + type char + optional size
    return /^[<>=|]?[biufcmMOSUV]\d*$/.test(fmt);
  }

  /**
   * Create the structured dtype from parsed components.
   */
  private _createDtype(
    formats: string[],
    names: string[],
    titles: (string | null)[],
    aligned: boolean,
    byteorder: string | null
  ): StructuredDType {
    const fields: FieldDescriptor[] = [];
    let offset = 0;

    for (let i = 0; i < formats.length; i++) {
      const { dtype, itemsize } = this._parseFormat(formats[i], byteorder);

      // Handle alignment padding
      if (aligned && offset > 0) {
        const alignment = dtypeAlignment(dtype);
        const padding = (alignment - (offset % alignment)) % alignment;
        offset += padding;
      }

      fields.push({
        name: names[i],
        dtype: dtype,
        offset: offset,
        title: titles[i],
        itemsize: itemsize,
      });

      offset += itemsize || dtypeSize(dtype);
    }

    return {
      names: names,
      fields: fields,
      itemsize: offset,
      alignment: aligned ? this._maxAlignment(fields) : 1,
    };
  }

  /**
   * Parse a single format string into dtype and size.
   */
  private _parseFormat(
    fmt: string,
    defaultByteorder: string | null
  ): { dtype: DType; itemsize: number } {
    // Extract components
    const match = fmt.match(/^([<>=|])?([biufcmMOSUV])(\d*)$/);
    if (!match) {
      throw new TypeError(`Invalid format: ${fmt}`);
    }

    const [, byteorder, typeChar, sizeStr] = match;
    const size = sizeStr ? parseInt(sizeStr, 10) : null;

    // Convert type character to DType
    const dtype = this._charToDtype(typeChar, size);

    // Calculate itemsize
    let itemsize: number;
    if (typeChar === 'S' || typeChar === 'U') {
      // String types: size is number of characters
      itemsize = size ?? 1;
      if (typeChar === 'U') {
        itemsize *= 4; // Unicode uses 4 bytes per char
      }
    } else if (typeChar === 'V') {
      // Void type: size is number of bytes
      itemsize = size ?? 1;
    } else {
      // Numeric types: use dtype size
      itemsize = dtypeSize(dtype);
    }

    return { dtype, itemsize };
  }

  /**
   * Convert type character to DType.
   */
  private _charToDtype(typeChar: string, size: number | null): DType {
    switch (typeChar) {
      case 'b':
        return DType.Int8;
      case 'i':
        switch (size) {
          case 1: return DType.Int8;
          case 2: return DType.Int16;
          case 4: return DType.Int32;
          case 8: return DType.Int64;
          default: return DType.Int32; // Default
        }
      case 'u':
        switch (size) {
          case 1: return DType.UInt8;
          case 2: return DType.UInt16;
          case 4: return DType.UInt32;
          case 8: return DType.UInt64;
          default: return DType.UInt32; // Default
        }
      case 'f':
        switch (size) {
          case 2: return DType.Float16;
          case 4: return DType.Float32;
          case 8: return DType.Float64;
          default: return DType.Float64; // Default
        }
      case 'c':
        switch (size) {
          case 8: return DType.Complex64;
          case 16: return DType.Complex128;
          default: return DType.Complex128; // Default
        }
      case 'S':
      case 'U':
        return DType.String;
      case 'V':
        return DType.UInt8; // Void treated as bytes
      case 'O':
        throw new TypeError('Object dtype not supported');
      case 'm':
      case 'M':
        throw new TypeError('Datetime dtypes not yet supported');
      default:
        throw new TypeError(`Unknown type character: ${typeChar}`);
    }
  }

  /**
   * Get maximum alignment requirement from fields.
   */
  private _maxAlignment(fields: FieldDescriptor[]): number {
    return Math.max(1, ...fields.map(f => dtypeAlignment(f.dtype)));
  }
}

/* ============ Error Classes ============ */

export class ValueError extends Error {
  constructor(message: string) {
    super(message);
    this.name = 'ValueError';
  }
}
```

---

### 16b.3 Record Class

**File:** `src/ts/rec/record.ts`

```typescript
import { NDArray } from '../NDArray.js';
import { StructuredDType, FieldDescriptor } from '../types.js';

/**
 * A single record from a record array.
 *
 * Allows field access via attributes or indexing.
 *
 * @example
 * const rec = recarray.getRecord(0);
 * console.log(rec.name);  // Field access via attribute
 * console.log(rec['age']); // Field access via indexing
 */
export class record {
  /** Reference to parent array */
  private _data: NDArray;

  /** Index of this record in parent array */
  private _index: number;

  /** Structured dtype describing fields */
  private _dtype: StructuredDType;

  /**
   * Create a record view.
   *
   * @param data - Parent recarray
   * @param index - Index of this record
   */
  constructor(data: NDArray, index: number) {
    this._data = data;
    this._index = index;
    this._dtype = data.structuredDtype!;

    // Create property getters/setters for each field
    for (const field of this._dtype.fields) {
      Object.defineProperty(this, field.name, {
        get: () => this._getField(field.name),
        set: (value) => this._setField(field.name, value),
        enumerable: true,
        configurable: true,
      });
    }
  }

  /**
   * Get field value by name.
   */
  private _getField(name: string): any {
    const fieldArr = this._data.getField(name);
    return fieldArr.getFlat(this._index);
  }

  /**
   * Set field value by name.
   */
  private _setField(name: string, value: any): void {
    const fieldArr = this._data.getField(name);
    fieldArr.setFlat(this._index, value);
  }

  /**
   * Number of fields in this record.
   */
  get length(): number {
    return this._dtype.fields.length;
  }

  /**
   * Get field names.
   */
  get dtype(): StructuredDType {
    return this._dtype;
  }

  /**
   * Index-based field access.
   */
  [Symbol.iterator](): Iterator<any> {
    let i = 0;
    const fields = this._dtype.fields;
    return {
      next: () => {
        if (i < fields.length) {
          return { value: this._getField(fields[i++].name), done: false };
        }
        return { value: undefined, done: true };
      }
    };
  }

  /**
   * Convert to plain JavaScript array.
   */
  toArray(): any[] {
    return this._dtype.fields.map(f => this._getField(f.name));
  }

  /**
   * Convert to plain JavaScript object.
   */
  toObject(): Record<string, any> {
    const result: Record<string, any> = {};
    for (const field of this._dtype.fields) {
      result[field.name] = this._getField(field.name);
    }
    return result;
  }

  /**
   * String representation.
   */
  toString(): string {
    const values = this._dtype.fields.map(f => {
      const val = this._getField(f.name);
      return `${f.name}=${typeof val === 'string' ? `'${val}'` : val}`;
    });
    return `(${values.join(', ')})`;
  }

  /**
   * Pretty print the record.
   */
  pprint(): void {
    console.log(this.toString());
  }

  /**
   * Support bracket notation for field access.
   */
  static {
    // Use Proxy for bracket notation support
  }
}

/**
 * Create a record with Proxy for bracket notation.
 */
export function createRecord(data: NDArray, index: number): record {
  const rec = new record(data, index);

  return new Proxy(rec, {
    get(target, prop) {
      if (typeof prop === 'string' && target._dtype.names.includes(prop)) {
        return target['_getField'](prop);
      }
      return (target as any)[prop];
    },
    set(target, prop, value) {
      if (typeof prop === 'string' && target._dtype.names.includes(prop)) {
        target['_setField'](prop, value);
        return true;
      }
      (target as any)[prop] = value;
      return true;
    }
  });
}
```

---

### 16b.4 RecArray Class

**File:** `src/ts/rec/recarray.ts`

```typescript
import { NDArray } from '../NDArray.js';
import { DType, StructuredDType, FieldDescriptor } from '../types.js';
import { format_parser } from './format_parser.js';
import { record, createRecord } from './record.js';

/**
 * Record array with named fields.
 *
 * A recarray allows field access as attributes in addition to
 * standard array indexing.
 *
 * @example
 * const arr = rec.fromarrays(
 *   [[25, 35, 45], ['Alice', 'Bob', 'Carol']],
 *   { names: ['age', 'name'] }
 * );
 *
 * arr.age;      // NDArray([25, 35, 45])
 * arr['name'];  // NDArray(['Alice', 'Bob', 'Carol'])
 * arr[0];       // record(age=25, name='Alice')
 */
export class recarray extends NDArray {
  /** Cached field arrays */
  private _fieldCache: Map<string, NDArray> = new Map();

  /**
   * Create a new record array.
   *
   * Use the factory methods (array, fromarrays, etc.) instead of
   * calling this constructor directly.
   */
  constructor() {
    super();
  }

  /**
   * Create a recarray from shape and dtype.
   */
  static create(
    shape: number | number[],
    dtype: StructuredDType | null = null,
    formats: string | string[] | null = null,
    names: string | string[] | null = null,
    titles: string | string[] | null = null,
    aligned: boolean = false,
    byteorder: string | null = null
  ): recarray {
    // Determine structured dtype
    let structDtype: StructuredDType;

    if (dtype !== null) {
      structDtype = dtype;
    } else if (formats !== null) {
      const parser = new format_parser(formats, names, titles, aligned, byteorder);
      structDtype = parser.dtype;
    } else {
      throw new TypeError('Must specify either dtype or formats');
    }

    // Create the recarray
    const shapeArr = typeof shape === 'number' ? [shape] : shape;
    const arr = new recarray();

    arr._initStructured(shapeArr, structDtype);

    return arr;
  }

  /**
   * Initialize as structured array.
   */
  private _initStructured(shape: number[], dtype: StructuredDType): void {
    this._shape = shape;
    this._size = shape.reduce((a, b) => a * b, 1);
    this._structuredDtype = dtype;

    // Allocate storage for each field
    this._fieldData = new Map();
    for (const field of dtype.fields) {
      const fieldShape = [...shape];
      const fieldArr = NDArray.empty(fieldShape, field.dtype);
      this._fieldData.set(field.name, fieldArr);
    }

    // Create property accessors for each field
    for (const field of dtype.fields) {
      if (!(field.name in this)) {
        Object.defineProperty(this, field.name, {
          get: () => this.field(field.name),
          set: (value) => this._setFieldArray(field.name, value),
          enumerable: true,
          configurable: true,
        });
      }
    }
  }

  /**
   * Get structured dtype.
   */
  get structuredDtype(): StructuredDType | undefined {
    return this._structuredDtype;
  }

  /**
   * Get a single field as an array.
   *
   * @param name - Field name
   * @returns NDArray containing field values
   */
  field(name: string): NDArray {
    if (!this._fieldData?.has(name)) {
      throw new KeyError(`Field '${name}' not found`);
    }

    // Return cached or create new view
    if (!this._fieldCache.has(name)) {
      this._fieldCache.set(name, this._fieldData.get(name)!);
    }

    return this._fieldCache.get(name)!;
  }

  /**
   * Set a field array.
   */
  private _setFieldArray(name: string, value: NDArray | number[] | string[]): void {
    const fieldArr = this.field(name);

    if (Array.isArray(value)) {
      for (let i = 0; i < value.length && i < fieldArr.size; i++) {
        fieldArr.setFlat(i, value[i]);
      }
    } else {
      for (let i = 0; i < value.size && i < fieldArr.size; i++) {
        fieldArr.setFlat(i, value.getFlat(i));
      }
    }

    // Clear cache
    this._fieldCache.delete(name);
  }

  /**
   * Get a single record by index.
   */
  getRecord(index: number): record {
    if (index < 0) index += this._size;
    if (index < 0 || index >= this._size) {
      throw new IndexError(`Index ${index} out of bounds for size ${this._size}`);
    }
    return createRecord(this, index);
  }

  /**
   * Get field data for a specific index.
   */
  getField(name: string): NDArray {
    return this.field(name);
  }

  /**
   * Iterate over records.
   */
  *[Symbol.iterator](): Iterator<record> {
    for (let i = 0; i < this._size; i++) {
      yield this.getRecord(i);
    }
  }

  /**
   * Convert to list of tuples.
   */
  tolist(): any[][] {
    const result: any[][] = [];
    const dtype = this._structuredDtype!;

    for (let i = 0; i < this._size; i++) {
      const row: any[] = [];
      for (const field of dtype.fields) {
        row.push(this.field(field.name).getFlat(i));
      }
      result.push(row);
    }

    return result;
  }

  /**
   * Convert to list of objects.
   */
  toObjects(): Record<string, any>[] {
    const result: Record<string, any>[] = [];
    const dtype = this._structuredDtype!;

    for (let i = 0; i < this._size; i++) {
      const obj: Record<string, any> = {};
      for (const field of dtype.fields) {
        obj[field.name] = this.field(field.name).getFlat(i);
      }
      result.push(obj);
    }

    return result;
  }

  /**
   * List field names (for tab completion support).
   */
  __dir__(): string[] {
    return this._structuredDtype?.names ?? [];
  }

  /**
   * String representation.
   */
  toString(): string {
    const dtype = this._structuredDtype;
    if (!dtype) return 'recarray([])';

    const lines: string[] = ['recarray(['];

    const maxShow = 10;
    const numToShow = Math.min(this._size, maxShow);

    for (let i = 0; i < numToShow; i++) {
      const rec = this.getRecord(i);
      lines.push(`  ${rec.toString()},`);
    }

    if (this._size > maxShow) {
      lines.push('  ...');
    }

    lines.push(`], dtype=[${dtype.names.map(n => `'${n}'`).join(', ')}])`);

    return lines.join('\n');
  }
}

/* ============ Error Classes ============ */

export class KeyError extends Error {
  constructor(message: string) {
    super(message);
    this.name = 'KeyError';
  }
}

export class IndexError extends Error {
  constructor(message: string) {
    super(message);
    this.name = 'IndexError';
  }
}
```

---

### 16b.5 Creation Functions

**File:** `src/ts/rec/index.ts`

```typescript
/**
 * NumJS Record Arrays Module
 *
 * Provides record array functionality compatible with NumPy's numpy.rec module.
 */

import { NDArray } from '../NDArray.js';
import { DType, StructuredDType, isStructuredDType } from '../types.js';
import { format_parser, ValueError } from './format_parser.js';
import { recarray, KeyError, IndexError } from './recarray.js';
import { record, createRecord } from './record.js';

export { format_parser, recarray, record, createRecord };
export { ValueError, KeyError, IndexError };

/**
 * Options for creating record arrays.
 */
export interface RecArrayOptions {
  dtype?: StructuredDType;
  formats?: string | string[];
  names?: string | string[];
  titles?: string | string[];
  aligned?: boolean;
  byteorder?: string;
}

/**
 * Create a record array.
 *
 * @param data - Input data (list of tuples or existing array)
 * @param options - Creation options
 * @returns A new recarray
 *
 * @example
 * const arr = rec.array([
 *   [25, 'Alice'],
 *   [30, 'Bob'],
 *   [35, 'Carol']
 * ], {
 *   formats: ['i4', 'U10'],
 *   names: ['age', 'name']
 * });
 */
export function array(
  data: any[][] | NDArray | null = null,
  options: RecArrayOptions = {}
): recarray {
  const {
    dtype = null,
    formats = null,
    names = null,
    titles = null,
    aligned = false,
    byteorder = null,
  } = options;

  // Determine structured dtype
  let structDtype: StructuredDType;

  if (dtype !== null) {
    structDtype = dtype;
  } else if (formats !== null) {
    const parser = new format_parser(formats, names, titles, aligned, byteorder);
    structDtype = parser.dtype;
  } else {
    throw new TypeError('Must specify either dtype or formats');
  }

  // Determine shape from data
  let shape: number[];
  if (data === null) {
    shape = [0];
  } else if (Array.isArray(data)) {
    shape = [data.length];
  } else {
    shape = data.shape;
  }

  // Create recarray
  const result = recarray.create(shape, structDtype);

  // Fill with data
  if (data !== null && Array.isArray(data)) {
    for (let i = 0; i < data.length; i++) {
      const row = data[i];
      for (let j = 0; j < structDtype.fields.length; j++) {
        const field = structDtype.fields[j];
        const value = row[j];
        result.field(field.name).setFlat(i, value);
      }
    }
  }

  return result;
}

/**
 * Create a record array from a list of arrays (columns).
 *
 * @param arrayList - List of arrays, one per field
 * @param options - Creation options
 * @returns A new recarray
 *
 * @example
 * const ages = [25, 30, 35];
 * const names = ['Alice', 'Bob', 'Carol'];
 * const arr = rec.fromarrays([ages, names], {
 *   names: ['age', 'name']
 * });
 */
export function fromarrays(
  arrayList: (NDArray | number[] | string[])[],
  options: RecArrayOptions = {}
): recarray {
  if (arrayList.length === 0) {
    throw new ValueError('arrayList must contain at least one array');
  }

  let {
    dtype = null,
    formats = null,
    names = null,
    titles = null,
    aligned = false,
    byteorder = null,
  } = options;

  // Infer formats from arrays if not provided
  if (formats === null && dtype === null) {
    formats = arrayList.map(arr => _inferFormat(arr));
  }

  // Determine shape from first array
  const firstArr = arrayList[0];
  const shape = Array.isArray(firstArr) ? [firstArr.length] : firstArr.shape;

  // Verify all arrays have compatible shapes
  for (let i = 1; i < arrayList.length; i++) {
    const arr = arrayList[i];
    const arrLen = Array.isArray(arr) ? arr.length : arr.size;
    if (arrLen !== shape[0]) {
      throw new ValueError(
        `Array ${i} has length ${arrLen}, expected ${shape[0]}`
      );
    }
  }

  // Create structured dtype
  let structDtype: StructuredDType;
  if (dtype !== null) {
    structDtype = dtype;
  } else {
    const parser = new format_parser(formats!, names, titles, aligned, byteorder);
    structDtype = parser.dtype;
  }

  // Create recarray and fill
  const result = recarray.create(shape, structDtype);

  for (let j = 0; j < arrayList.length; j++) {
    const field = structDtype.fields[j];
    const srcArr = arrayList[j];
    const dstArr = result.field(field.name);

    if (Array.isArray(srcArr)) {
      for (let i = 0; i < srcArr.length; i++) {
        dstArr.setFlat(i, srcArr[i]);
      }
    } else {
      for (let i = 0; i < srcArr.size; i++) {
        dstArr.setFlat(i, srcArr.getFlat(i));
      }
    }
  }

  return result;
}

/**
 * Create a record array from a list of records (rows).
 *
 * @param recList - List of records (tuples or arrays)
 * @param options - Creation options
 * @returns A new recarray
 *
 * @example
 * const arr = rec.fromrecords([
 *   [25, 'Alice'],
 *   [30, 'Bob']
 * ], {
 *   names: ['age', 'name']
 * });
 */
export function fromrecords(
  recList: any[][],
  options: RecArrayOptions = {}
): recarray {
  return array(recList, options);
}

/**
 * Create a record array from a binary string/buffer.
 *
 * @param datastring - Binary data (ArrayBuffer or typed array)
 * @param dtype - Structured dtype
 * @param shape - Shape of output array
 * @param offset - Byte offset into data
 * @returns A new recarray
 */
export function fromstring(
  datastring: ArrayBuffer | Uint8Array,
  dtype: StructuredDType,
  shape: number[] = [],
  offset: number = 0
): recarray {
  // Get buffer
  let buffer: ArrayBuffer;
  if (datastring instanceof ArrayBuffer) {
    buffer = datastring;
  } else {
    buffer = datastring.buffer;
    offset += datastring.byteOffset;
  }

  // Calculate number of records
  const itemsize = dtype.itemsize;
  const availableBytes = buffer.byteLength - offset;
  const numRecords = shape.length > 0
    ? shape.reduce((a, b) => a * b, 1)
    : Math.floor(availableBytes / itemsize);

  const finalShape = shape.length > 0 ? shape : [numRecords];

  // Create recarray
  const result = recarray.create(finalShape, dtype);

  // Parse binary data
  const view = new DataView(buffer, offset);
  let byteOffset = 0;

  for (let i = 0; i < numRecords; i++) {
    for (const field of dtype.fields) {
      const fieldOffset = byteOffset + field.offset;
      const value = _readFieldValue(view, fieldOffset, field.dtype, field.itemsize);
      result.field(field.name).setFlat(i, value);
    }
    byteOffset += itemsize;
  }

  return result;
}

/**
 * Create a record array from a file.
 *
 * @param file - File path (Node.js) or File object (browser)
 * @param dtype - Structured dtype
 * @param shape - Shape of output array
 * @param offset - Byte offset into file
 * @returns Promise resolving to a new recarray
 */
export async function fromfile(
  file: File | string,
  dtype: StructuredDType,
  shape: number[] = [],
  offset: number = 0
): Promise<recarray> {
  let buffer: ArrayBuffer;

  if (typeof file === 'string') {
    // Node.js: read from file path
    if (typeof process !== 'undefined' && process.versions?.node) {
      const fs = await import('fs/promises');
      const data = await fs.readFile(file);
      buffer = data.buffer.slice(data.byteOffset, data.byteOffset + data.byteLength);
    } else {
      throw new Error('File paths only supported in Node.js');
    }
  } else {
    // Browser: read from File object
    buffer = await file.arrayBuffer();
  }

  return fromstring(buffer, dtype, shape, offset);
}

/**
 * Find duplicate elements in a list.
 * Utility function for checking field names.
 */
export function find_duplicate(list: string[]): string[] {
  const counts = new Map<string, number>();

  for (const item of list) {
    counts.set(item, (counts.get(item) || 0) + 1);
  }

  return Array.from(counts.entries())
    .filter(([_, count]) => count > 1)
    .map(([item, _]) => item);
}

/* ============ Helper Functions ============ */

/**
 * Infer format string from array.
 */
function _inferFormat(arr: NDArray | number[] | string[]): string {
  if (arr instanceof NDArray) {
    return _dtypeToFormat(arr.dtype);
  }

  if (arr.length === 0) {
    return 'f8'; // Default to float64
  }

  const sample = arr[0];
  if (typeof sample === 'string') {
    const maxLen = Math.max(...(arr as string[]).map(s => s.length));
    return `U${maxLen || 1}`;
  }
  if (typeof sample === 'number') {
    if (Number.isInteger(sample)) {
      return 'i4';
    }
    return 'f8';
  }
  if (typeof sample === 'boolean') {
    return 'b1';
  }

  return 'O'; // Object (not supported)
}

/**
 * Convert DType to format string.
 */
function _dtypeToFormat(dtype: DType): string {
  const formats: Record<DType, string> = {
    [DType.Bool]: 'b1',
    [DType.Int8]: 'i1',
    [DType.Int16]: 'i2',
    [DType.Int32]: 'i4',
    [DType.Int64]: 'i8',
    [DType.UInt8]: 'u1',
    [DType.UInt16]: 'u2',
    [DType.UInt32]: 'u4',
    [DType.UInt64]: 'u8',
    [DType.Float16]: 'f2',
    [DType.Float32]: 'f4',
    [DType.Float64]: 'f8',
    [DType.Complex64]: 'c8',
    [DType.Complex128]: 'c16',
    [DType.String]: 'U',
  };
  return formats[dtype] || 'V';
}

/**
 * Read a field value from binary data.
 */
function _readFieldValue(
  view: DataView,
  offset: number,
  dtype: DType,
  itemsize?: number
): any {
  switch (dtype) {
    case DType.Bool:
      return view.getUint8(offset) !== 0;
    case DType.Int8:
      return view.getInt8(offset);
    case DType.Int16:
      return view.getInt16(offset, true);
    case DType.Int32:
      return view.getInt32(offset, true);
    case DType.UInt8:
      return view.getUint8(offset);
    case DType.UInt16:
      return view.getUint16(offset, true);
    case DType.UInt32:
      return view.getUint32(offset, true);
    case DType.Float32:
      return view.getFloat32(offset, true);
    case DType.Float64:
      return view.getFloat64(offset, true);
    case DType.String:
      // Read fixed-width string
      if (itemsize) {
        const bytes = new Uint8Array(view.buffer, view.byteOffset + offset, itemsize);
        const decoder = new TextDecoder('utf-8');
        return decoder.decode(bytes).replace(/\0+$/, ''); // Trim null bytes
      }
      return '';
    default:
      return 0;
  }
}
```

---

## File Changes Summary

### New Files to Create

```
src/ts/rec/
├── index.ts           # Public exports and creation functions (~350 lines)
├── format_parser.ts   # Format string parsing (~250 lines)
├── record.ts          # Single record class (~150 lines)
└── recarray.ts        # Record array class (~300 lines)
```

### Files to Modify

```
src/ts/types.ts
├── Add StructuredDType interface
├── Add FieldDescriptor interface
├── Add isStructuredDType type guard
├── Add dtypeSize function
└── Add dtypeAlignment function

src/ts/NDArray.ts
├── Add _structuredDtype property
├── Add _fieldData property (Map<string, NDArray>)
├── Add structuredDtype getter
├── Add getField method
└── Add emptyStructured factory method

src/ts/index.ts
└── Export rec module
```

---

## Implementation Order

```
Week 3: numpy.rec Module
├── Day 1: Type infrastructure
│   ├── Add StructuredDType interface to types.ts
│   ├── Add FieldDescriptor interface
│   ├── Add dtypeSize, dtypeAlignment utilities
│   └── Update NDArray with structured support
│
├── Day 2: Format parser
│   ├── Implement format_parser class
│   ├── Parse format strings ('f8, i4, S10')
│   ├── Parse format lists
│   └── Handle alignment and byte order
│
├── Day 3: Record class
│   ├── Implement record class
│   ├── Add field access via attributes
│   ├── Add iteration support
│   └── Add toString, toArray, toObject
│
├── Day 4: RecArray class
│   ├── Implement recarray class
│   ├── Add field() method
│   ├── Add getRecord() method
│   ├── Add iteration over records
│   └── Add tolist(), toObjects()
│
└── Day 5: Creation functions & tests
    ├── Implement array()
    ├── Implement fromarrays()
    ├── Implement fromrecords()
    ├── Implement fromstring(), fromfile()
    └── Write comprehensive tests
```

---

## Verification Plan

```bash
# Build
npm run build

# Run tests
npm test -- --grep "rec"

# Test cases to verify:

# Format Parser
✓ format_parser('f8, i4') creates correct dtype
✓ format_parser with names creates named fields
✓ format_parser detects duplicate names
✓ format_parser handles alignment

# Record
✓ record.fieldname returns field value
✓ record['fieldname'] returns field value
✓ record.toArray() returns tuple
✓ record.toObject() returns dict-like object

# RecArray
✓ recarray.fieldname returns field array
✓ recarray.field('name') returns field array
✓ recarray.getRecord(i) returns record
✓ recarray iteration yields records
✓ recarray.tolist() returns list of tuples

# Creation Functions
✓ array([[1,'a'],[2,'b']], formats=['i4','U1']) works
✓ fromarrays([[1,2], ['a','b']], names=['x','y']) works
✓ fromrecords([[1,'a'],[2,'b']], names=['x','y']) works
✓ fromstring(buffer, dtype) parses binary data
```

---

## API Compatibility Notes

### NumPy API Mapping

```python
# NumPy
rec = np.rec.array(data, dtype=[('x', 'f8'), ('y', 'i4')])
rec.x  # Field access
rec[0]  # Record access

# NumJS
const rec = rec.array(data, { formats: ['f8', 'i4'], names: ['x', 'y'] });
rec.x;  // Field access
rec.getRecord(0);  // Record access
```

### Differences from NumPy

1. **Dtype specification**: NumJS uses `formats` + `names` separately instead of combined dtype list
2. **Record indexing**: Use `getRecord(i)` instead of `arr[i]` for single record access
3. **No view semantics**: Field arrays are stored separately, not as views into contiguous memory

---

## Estimated Lines of Code

| File | Lines |
|------|-------|
| format_parser.ts | ~250 |
| record.ts | ~150 |
| recarray.ts | ~300 |
| index.ts | ~350 |
| **Total** | **~1,050** |

Plus modifications to types.ts (~50 lines) and NDArray.ts (~80 lines).
