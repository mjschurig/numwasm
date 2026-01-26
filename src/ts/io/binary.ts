/**
 * Binary File I/O
 *
 * Functions for reading arrays from raw binary data.
 *
 * Reference: numpy/core/records.py, numpy/lib/_npyio_impl.py
 */

import { NDArray } from '../NDArray.js';
import { DType, DTYPE_SIZES } from '../types.js';
import { isNode } from './format.js';

/**
 * Options for fromfile()
 */
export interface FromfileOptions {
  /** Data type of the file */
  dtype?: DType;
  /** Number of items to read (-1 for all) */
  count?: number;
  /** Separator string (empty for binary) */
  sep?: string;
  /** Byte offset from start of file */
  offset?: number;
}

/**
 * Options for frombuffer()
 */
export interface FrombufferOptions {
  /** Data type to interpret buffer as */
  dtype?: DType;
  /** Number of items to read (-1 for all) */
  count?: number;
  /** Byte offset into buffer */
  offset?: number;
}

/**
 * Construct an array from data in a binary file.
 *
 * @param file - File source (path, File object, URL, or ArrayBuffer)
 * @param options - Reading options
 * @returns Loaded array
 *
 * @example
 * // Read as float64 (default)
 * const arr = await fromfile('data.bin');
 *
 * // Read specific dtype
 * const arr = await fromfile('data.bin', { dtype: DType.Int32 });
 *
 * // Read first 100 elements
 * const arr = await fromfile('data.bin', { dtype: DType.Float32, count: 100 });
 *
 * // Read with offset
 * const arr = await fromfile('data.bin', { dtype: DType.Float64, offset: 1024 });
 *
 * // Read text file with separator
 * const arr = await fromfile('data.txt', { sep: ',' });
 */
export async function fromfile(
  file: string | File | ArrayBuffer | URL,
  options: FromfileOptions = {}
): Promise<NDArray> {
  const {
    dtype = DType.Float64,
    count = -1,
    sep = '',
    offset = 0,
  } = options;

  // Text mode if separator is provided
  if (sep !== '') {
    return fromfileText(file, dtype, count, sep, offset);
  }

  // Binary mode
  const buffer = await readFileBuffer(file);
  const bytes = new Uint8Array(buffer);

  // Apply offset
  const startOffset = offset;
  if (startOffset >= bytes.length) {
    throw new Error(`Offset ${startOffset} exceeds file size ${bytes.length}`);
  }

  const availableBytes = bytes.length - startOffset;
  const itemSize = DTYPE_SIZES[dtype];
  const maxItems = Math.floor(availableBytes / itemSize);

  // Determine number of items to read
  const numItems = count === -1 ? maxItems : Math.min(count, maxItems);

  if (numItems === 0) {
    return NDArray.zeros([0], { dtype });
  }

  // Create typed array from buffer
  const dataBytes = numItems * itemSize;
  const slice = bytes.slice(startOffset, startOffset + dataBytes);
  const typedArray = createTypedArrayFromBuffer(dtype, slice.buffer, slice.byteOffset, numItems);

  return NDArray.fromTypedArray(typedArray, [numItems], dtype);
}

/**
 * Read text file with separator.
 */
async function fromfileText(
  file: string | File | ArrayBuffer | URL,
  dtype: DType,
  count: number,
  sep: string,
  offset: number
): Promise<NDArray> {
  // Read as text
  let text: string;

  if (file instanceof ArrayBuffer) {
    const decoder = new TextDecoder('utf-8');
    text = decoder.decode(file);
  } else if (typeof File !== 'undefined' && file instanceof File) {
    text = await file.text();
  } else if (file instanceof URL || (typeof file === 'string' && file.startsWith('http'))) {
    const response = await fetch(file.toString());
    text = await response.text();
  } else if (typeof file === 'string') {
    if (isNode()) {
      const fs = await import('fs/promises');
      text = await fs.readFile(file, { encoding: 'utf-8' });
    } else {
      throw new Error('String file paths only supported in Node.js');
    }
  } else {
    throw new Error('Unsupported file source');
  }

  // Skip offset bytes (as characters in text mode)
  text = text.slice(offset);

  // Split by separator and parse
  const parts = text.split(sep).filter(s => s.trim().length > 0);
  const maxItems = parts.length;
  const numItems = count === -1 ? maxItems : Math.min(count, maxItems);

  const data = parts.slice(0, numItems).map(s => parseFloat(s.trim()));

  return NDArray.fromArray(data, [data.length], { dtype });
}

/**
 * Interpret a buffer as a 1-dimensional array.
 *
 * @param buffer - Buffer containing array data
 * @param options - Interpretation options
 * @returns Promise resolving to array from buffer data
 *
 * @example
 * // From ArrayBuffer
 * const buffer = new ArrayBuffer(32);
 * const arr = await frombuffer(buffer, { dtype: DType.Float64 });
 *
 * // From TypedArray
 * const floats = new Float32Array([1, 2, 3, 4]);
 * const arr = await frombuffer(floats.buffer, { dtype: DType.Float32 });
 *
 * // With offset
 * const arr = await frombuffer(buffer, { dtype: DType.Int32, offset: 8 });
 *
 * // Read specific count
 * const arr = await frombuffer(buffer, { dtype: DType.Float64, count: 2 });
 */
export async function frombuffer(
  buffer: ArrayBuffer | ArrayBufferView,
  options: FrombufferOptions = {}
): Promise<NDArray> {
  const {
    dtype = DType.Float64,
    count = -1,
    offset = 0,
  } = options;

  // Get underlying ArrayBuffer
  let arrayBuffer: ArrayBuffer;
  let baseOffset: number;

  if (buffer instanceof ArrayBuffer) {
    arrayBuffer = buffer;
    baseOffset = 0;
  } else {
    // ArrayBufferView (TypedArray, DataView)
    arrayBuffer = buffer.buffer as ArrayBuffer;
    baseOffset = buffer.byteOffset;
  }

  const totalOffset = baseOffset + offset;
  const availableBytes = arrayBuffer.byteLength - totalOffset;
  const itemSize = DTYPE_SIZES[dtype];

  if (totalOffset > arrayBuffer.byteLength) {
    throw new Error(`Offset ${offset} exceeds buffer size`);
  }

  const maxItems = Math.floor(availableBytes / itemSize);
  const numItems = count === -1 ? maxItems : Math.min(count, maxItems);

  if (numItems === 0) {
    return NDArray.zeros([0], { dtype });
  }

  // Create typed array view
  const typedArray = createTypedArrayFromBuffer(dtype, arrayBuffer, totalOffset, numItems);

  // Create NDArray from typed array
  return NDArray.fromTypedArray(typedArray, [numItems], dtype);
}

/**
 * Read binary data from various file sources.
 */
async function readFileBuffer(
  file: string | File | ArrayBuffer | URL
): Promise<ArrayBuffer> {
  if (file instanceof ArrayBuffer) {
    return file;
  }

  // Check for File in browser
  if (typeof File !== 'undefined' && file instanceof File) {
    return file.arrayBuffer();
  }

  if (file instanceof URL || (typeof file === 'string' && file.startsWith('http'))) {
    const response = await fetch(file.toString());
    return response.arrayBuffer();
  }

  if (typeof file === 'string') {
    // Node.js file path
    if (isNode()) {
      const fs = await import('fs/promises');
      const buffer = await fs.readFile(file);
      return buffer.buffer.slice(
        buffer.byteOffset,
        buffer.byteOffset + buffer.byteLength
      );
    }
    throw new Error('String file paths only supported in Node.js');
  }

  throw new Error('Unsupported file source');
}

/**
 * Create a typed array from a buffer based on dtype.
 */
function createTypedArrayFromBuffer(
  dtype: DType,
  buffer: ArrayBuffer,
  byteOffset: number,
  length: number
): Float64Array | Float32Array | Int32Array | Int16Array | Int8Array | Uint32Array | Uint16Array | Uint8Array | BigInt64Array | BigUint64Array {
  switch (dtype) {
    case DType.Float64:
      return new Float64Array(buffer, byteOffset, length);
    case DType.Float32:
      return new Float32Array(buffer, byteOffset, length);
    case DType.Int64:
      return new BigInt64Array(buffer, byteOffset, length);
    case DType.Int32:
      return new Int32Array(buffer, byteOffset, length);
    case DType.Int16:
      return new Int16Array(buffer, byteOffset, length);
    case DType.Int8:
      return new Int8Array(buffer, byteOffset, length);
    case DType.Uint64:
      return new BigUint64Array(buffer, byteOffset, length);
    case DType.Uint32:
      return new Uint32Array(buffer, byteOffset, length);
    case DType.Uint16:
      return new Uint16Array(buffer, byteOffset, length);
    case DType.Uint8:
    case DType.Bool:
      return new Uint8Array(buffer, byteOffset, length);
    case DType.Float16:
      // Float16 not natively supported, store as Uint16
      return new Uint16Array(buffer, byteOffset, length);
    case DType.Complex64:
      // Complex64 = 2 x Float32
      return new Float32Array(buffer, byteOffset, length * 2);
    case DType.Complex128:
      // Complex128 = 2 x Float64
      return new Float64Array(buffer, byteOffset, length * 2);
    default:
      return new Float64Array(buffer, byteOffset, length);
  }
}
