/**
 * NumJS Record Arrays Module
 *
 * Provides record array functionality compatible with NumPy's numpy.rec module.
 * Record arrays allow field access as attributes in addition to standard indexing.
 *
 * @example
 * import { rec } from 'numjs';
 *
 * // Create from column arrays
 * const arr = rec.fromarrays(
 *   [[25, 30, 35], ['Alice', 'Bob', 'Carol']],
 *   { names: ['age', 'name'] }
 * );
 *
 * // Field access
 * console.log(arr.age);     // NDArray([25, 30, 35])
 * console.log(arr.field('name')); // NDArray(['Alice', 'Bob', 'Carol'])
 *
 * // Record access
 * console.log(arr.getRecord(0)); // record(age=25, name='Alice')
 */

import { NDArray } from '../NDArray.js';
import { DType, type StructuredDType } from '../types.js';
import { format_parser, ValueError, find_duplicate, inferFormat, dtypeToFormat } from './format_parser.js';
import { recarray, IndexError, createRecarray } from './recarray.js';
import { record, KeyError } from './record.js';

/* ============ Re-exports ============ */

export { format_parser, ValueError, find_duplicate, inferFormat, dtypeToFormat } from './format_parser.js';
export { recarray, IndexError, createRecarray } from './recarray.js';
export { record, createRecord, KeyError } from './record.js';

/* ============ Types ============ */

/**
 * Options for creating record arrays.
 */
export interface RecArrayOptions {
  /** Pre-built structured dtype */
  dtype?: StructuredDType;

  /** Format specifications (e.g., 'f8, i4' or ['f8', 'i4']) */
  formats?: string | string[];

  /** Field names (e.g., 'x, y' or ['x', 'y']) */
  names?: string | string[];

  /** Optional field titles */
  titles?: string | string[];

  /** Whether to align fields for C struct compatibility */
  aligned?: boolean;

  /** Byte order: '<' little, '>' big, '=' native */
  byteorder?: '<' | '>' | '=' | '|';

  /** Whether to copy input data */
  copy?: boolean;
}

/* ============ Creation Functions ============ */

/**
 * Create a record array.
 *
 * @param data - Input data (list of tuples, existing array, or null for empty)
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
export async function array(
  data: unknown[][] | NDArray | recarray | null = null,
  options: RecArrayOptions = {}
): Promise<recarray> {
  const {
    dtype = undefined,
    formats = undefined,
    names = undefined,
    titles = undefined,
    aligned = false,
    byteorder = undefined,
    copy = false,
  } = options;

  // Handle recarray input
  if (data instanceof recarray) {
    return copy ? data.copy() : data;
  }

  // Determine structured dtype
  let structDtype: StructuredDType;

  if (dtype !== undefined) {
    structDtype = dtype;
  } else if (formats !== undefined) {
    const parser = new format_parser(formats, names ?? null, titles ?? null, aligned, byteorder ?? null);
    structDtype = parser.dtype;
  } else if (data !== null && Array.isArray(data) && data.length > 0) {
    // Try to infer from data
    const firstRow = data[0];
    if (Array.isArray(firstRow)) {
      const inferredFormats = firstRow.map((val) => {
        if (typeof val === 'string') return 'U' + Math.max(1, val.length);
        if (typeof val === 'number') return Number.isInteger(val) ? 'i4' : 'f8';
        if (typeof val === 'boolean') return '?';
        return 'f8';
      });
      const parser = new format_parser(inferredFormats, names ?? null, titles ?? null, aligned, byteorder ?? null);
      structDtype = parser.dtype;
    } else {
      throw new TypeError('Cannot infer dtype from data');
    }
  } else {
    throw new TypeError('Must specify either dtype or formats');
  }

  // Determine shape from data
  let shape: number[];
  if (data === null) {
    shape = [0];
  } else if (Array.isArray(data)) {
    shape = [data.length];
  } else if (data instanceof NDArray) {
    shape = data.shape;
  } else {
    shape = [0];
  }

  // Create field arrays
  const fields = new Map<string, NDArray>();

  for (const field of structDtype.fieldList) {
    if (field.dtype === DType.String) {
      fields.set(field.name, NDArray.emptyString(shape));
    } else {
      fields.set(field.name, await NDArray.zeros(shape, { dtype: field.dtype }));
    }
  }

  // Fill with data
  if (data !== null && Array.isArray(data)) {
    for (let i = 0; i < data.length; i++) {
      const row = data[i] as unknown[];
      for (let j = 0; j < structDtype.fieldList.length && j < row.length; j++) {
        const field = structDtype.fieldList[j];
        const fieldArr = fields.get(field.name)!;
        const value = row[j];

        if (field.dtype === DType.String) {
          fieldArr.setStringFlat(i, String(value));
        } else {
          fieldArr.setFlat(i, value as number);
        }
      }
    }
  }

  // Create and return proxied recarray
  return createRecarray(structDtype, shape, fields);
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
 * const arr = await rec.fromarrays([ages, names], {
 *   names: ['age', 'name']
 * });
 */
export async function fromarrays(
  arrayList: (NDArray | number[] | string[])[],
  options: RecArrayOptions = {}
): Promise<recarray> {
  if (arrayList.length === 0) {
    throw new ValueError('arrayList must contain at least one array');
  }

  const {
    dtype = undefined,
    formats = undefined,
    names = undefined,
    titles = undefined,
    aligned = false,
    byteorder = undefined,
  } = options;

  // Infer formats from arrays if not provided
  let structDtype: StructuredDType;

  if (dtype !== undefined) {
    structDtype = dtype;
  } else {
    let formatList: string[];

    if (formats !== undefined) {
      formatList = typeof formats === 'string'
        ? formats.split(',').map(f => f.trim())
        : formats;
    } else {
      // Infer formats from arrays
      formatList = arrayList.map(arr => {
        if (arr instanceof NDArray) {
          return dtypeToFormat(arr.dtype);
        }
        return inferFormat(arr as unknown[]);
      });
    }

    const parser = new format_parser(formatList, names ?? null, titles ?? null, aligned, byteorder ?? null);
    structDtype = parser.dtype;
  }

  if (arrayList.length !== structDtype.fieldList.length) {
    throw new ValueError(
      `Number of arrays (${arrayList.length}) must match number of fields (${structDtype.fieldList.length})`
    );
  }

  // Determine shape from first array
  const firstArr = arrayList[0];
  const size = Array.isArray(firstArr) ? firstArr.length : firstArr.size;
  const shape = [size];

  // Verify all arrays have compatible lengths
  for (let i = 1; i < arrayList.length; i++) {
    const arr = arrayList[i];
    const arrLen = Array.isArray(arr) ? arr.length : arr.size;
    if (arrLen !== size) {
      throw new ValueError(
        `Array ${i} has length ${arrLen}, expected ${size}`
      );
    }
  }

  // Create field arrays
  const fields = new Map<string, NDArray>();

  for (let j = 0; j < arrayList.length; j++) {
    const field = structDtype.fieldList[j];
    const srcArr = arrayList[j];

    if (srcArr instanceof NDArray) {
      // Copy the NDArray
      fields.set(field.name, srcArr.copy());
    } else if (field.dtype === DType.String) {
      // String array from JavaScript array
      fields.set(field.name, NDArray.fromStringArray(srcArr as string[]));
    } else {
      // Numeric array from JavaScript array
      fields.set(field.name, await NDArray.fromArray(srcArr as number[], shape, { dtype: field.dtype }));
    }
  }

  return createRecarray(structDtype, shape, fields);
}

/**
 * Create a record array from a list of records (rows).
 *
 * @param recList - List of records (tuples or arrays)
 * @param options - Creation options
 * @returns A new recarray
 *
 * @example
 * const arr = await rec.fromrecords([
 *   [25, 'Alice'],
 *   [30, 'Bob']
 * ], {
 *   names: ['age', 'name']
 * });
 */
export async function fromrecords(
  recList: unknown[][],
  options: RecArrayOptions = {}
): Promise<recarray> {
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
export async function fromstring(
  datastring: ArrayBuffer | Uint8Array,
  dtype: StructuredDType,
  shape: number[] = [],
  offset: number = 0
): Promise<recarray> {
  // Get buffer
  let buffer: ArrayBuffer;
  let bufferOffset = offset;

  if (datastring instanceof ArrayBuffer) {
    buffer = datastring;
  } else {
    buffer = datastring.buffer as ArrayBuffer;
    bufferOffset += datastring.byteOffset;
  }

  // Calculate number of records
  const itemsize = dtype.itemsize;
  const availableBytes = buffer.byteLength - bufferOffset;
  const numRecords = shape.length > 0
    ? shape.reduce((a, b) => a * b, 1)
    : Math.floor(availableBytes / itemsize);

  const finalShape = shape.length > 0 ? shape : [numRecords];

  // Create field arrays
  const fields = new Map<string, NDArray>();

  for (const field of dtype.fieldList) {
    if (field.dtype === DType.String) {
      fields.set(field.name, NDArray.emptyString(finalShape));
    } else {
      fields.set(field.name, await NDArray.zeros(finalShape, { dtype: field.dtype }));
    }
  }

  // Parse binary data
  const view = new DataView(buffer, bufferOffset);
  let byteOffset = 0;

  for (let i = 0; i < numRecords; i++) {
    for (const field of dtype.fieldList) {
      const fieldOffset = byteOffset + field.offset;
      const value = readFieldValue(view, fieldOffset, field.dtype, field.itemsize, field.charType);
      const fieldArr = fields.get(field.name)!;

      if (field.dtype === DType.String) {
        fieldArr.setStringFlat(i, value as string);
      } else {
        fieldArr.setFlat(i, value as number);
      }
    }
    byteOffset += itemsize;
  }

  return createRecarray(dtype, finalShape, fields);
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

/* ============ Helper Functions ============ */

/**
 * Read a field value from binary data.
 */
function readFieldValue(
  view: DataView,
  offset: number,
  dtype: DType,
  itemsize: number,
  charType?: 'S' | 'U'
): unknown {
  switch (dtype) {
    case DType.Bool:
      return view.getUint8(offset) !== 0;
    case DType.Int8:
      return view.getInt8(offset);
    case DType.Int16:
      return view.getInt16(offset, true);
    case DType.Int32:
      return view.getInt32(offset, true);
    case DType.Int64:
      // Read as two 32-bit values (limited precision)
      return view.getInt32(offset, true);
    case DType.Uint8:
      return view.getUint8(offset);
    case DType.Uint16:
      return view.getUint16(offset, true);
    case DType.Uint32:
      return view.getUint32(offset, true);
    case DType.Uint64:
      return view.getUint32(offset, true);
    case DType.Float16:
      // Float16 not directly supported, read as bytes
      return view.getFloat32(offset, true);
    case DType.Float32:
      return view.getFloat32(offset, true);
    case DType.Float64:
      return view.getFloat64(offset, true);
    case DType.String:
      // Read fixed-width string
      if (charType === 'U') {
        // Unicode: 4 bytes per character
        const chars: string[] = [];
        const numChars = itemsize / 4;
        for (let i = 0; i < numChars; i++) {
          const code = view.getUint32(offset + i * 4, true);
          if (code === 0) break;
          chars.push(String.fromCodePoint(code));
        }
        return chars.join('');
      } else {
        // ASCII: 1 byte per character
        const bytes = new Uint8Array(view.buffer, view.byteOffset + offset, itemsize);
        const decoder = new TextDecoder('utf-8');
        return decoder.decode(bytes).replace(/\0+$/, ''); // Trim null bytes
      }
    default:
      return 0;
  }
}

/* ============ Namespace Object ============ */

/**
 * rec namespace object for convenient grouped import.
 *
 * @example
 * import { rec } from 'numjs';
 * const arr = await rec.fromarrays([[1,2,3], ['a','b','c']], {names: ['x', 'y']});
 */
export const rec = {
  array,
  fromarrays,
  fromrecords,
  fromstring,
  fromfile,
  format_parser,
  recarray,
  record,
  find_duplicate,
  ValueError,
  KeyError,
  IndexError,
};

export default rec;
