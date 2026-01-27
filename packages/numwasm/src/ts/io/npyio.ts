/**
 * NPY/NPZ File I/O
 *
 * Functions for saving and loading NDArrays in NumPy's binary format.
 *
 * Reference: numpy/lib/_npyio_impl.py
 */

import { NDArray } from '../NDArray.js';
import { DType } from '../types.js';
import {
  createHeader,
  parseHeader,
  dtypeToDescr,
  descrToDtype,
  dtypeSize,
  getMinVersion,
  isNode,
  type NpyHeader,
} from './format.js';

/**
 * Options for save()
 */
export interface SaveOptions {
  /** Allow saving object arrays (not supported in NumJS) */
  allow_pickle?: boolean;
}

/**
 * Options for load()
 */
export interface LoadOptions {
  /** Memory mapping mode (not yet supported) */
  mmap_mode?: 'r+' | 'r' | 'w+' | 'c' | null;
  /** Allow loading pickled objects (not supported) */
  allow_pickle?: boolean;
  /** Maximum header size */
  max_header_size?: number;
}

/**
 * Save an array to a binary file in NPY format.
 *
 * @param file - File path (Node.js), FileSystemFileHandle (browser), or null to return ArrayBuffer
 * @param arr - Array to save
 * @param options - Save options
 * @returns ArrayBuffer if file is null, void otherwise
 *
 * @example
 * // In Node.js
 * await save('data.npy', arr);
 *
 * // In browser (returns ArrayBuffer)
 * const buffer = await save(null, arr);
 *
 * // With File System Access API
 * const handle = await showSaveFilePicker();
 * await save(handle, arr);
 */
export async function save(
  file: string | FileSystemFileHandle | null,
  arr: NDArray,
  _options: SaveOptions = {}
): Promise<ArrayBuffer | void> {
  // Ensure array is contiguous for efficient writing
  const contiguous = arr.flags.c_contiguous ? arr : arr.copy();

  // Build header
  const header: NpyHeader = {
    descr: dtypeToDescr(contiguous.dtype),
    fortran_order: false,
    shape: [...contiguous.shape],
  };

  // Calculate header size to determine version
  const shapeStr = header.shape.length === 0
    ? '()'
    : header.shape.length === 1
    ? `(${header.shape[0]},)`
    : `(${header.shape.join(', ')})`;
  const headerStrLen = `{'descr': '${header.descr}', 'fortran_order': False, 'shape': ${shapeStr}, }`.length;
  const version = getMinVersion(headerStrLen);

  // Create header bytes
  const headerBytes = createHeader(header, version);

  // Get array data as bytes
  const dataBytes = contiguous.toTypedArray();
  const dataView = new Uint8Array(dataBytes.buffer, dataBytes.byteOffset, contiguous.size * dtypeSize(contiguous.dtype));

  // Combine header and data
  const totalSize = headerBytes.length + dataView.length;
  const result = new Uint8Array(totalSize);
  result.set(headerBytes, 0);
  result.set(dataView, headerBytes.length);

  // Write to destination
  if (file === null) {
    return result.buffer;
  }

  if (typeof file === 'string') {
    // Node.js file path
    if (isNode()) {
      const fs = await import('fs/promises');
      await fs.writeFile(file, result);
      return;
    }
    throw new Error('String file paths only supported in Node.js');
  }

  // FileSystemFileHandle (File System Access API)
  if ('createWritable' in file) {
    const writable = await (file as FileSystemFileHandle).createWritable();
    await writable.write(result);
    await writable.close();
    return;
  }

  throw new Error('Unsupported file target');
}

/**
 * Load an array from a NPY file.
 *
 * @param file - File path (Node.js), File object (browser), ArrayBuffer, or URL
 * @param options - Load options
 * @returns Loaded NDArray
 *
 * @example
 * // From ArrayBuffer
 * const arr = await load(buffer);
 *
 * // From File object (browser)
 * const arr = await load(fileInput.files[0]);
 *
 * // From URL
 * const arr = await load('https://example.com/data.npy');
 *
 * // In Node.js
 * const arr = await load('data.npy');
 */
export async function load(
  file: string | File | ArrayBuffer | URL,
  _options: LoadOptions = {}
): Promise<NDArray> {
  const data = await readFileData(file);
  const bytes = new Uint8Array(data);

  // Parse header
  const { header, dataOffset } = parseHeader(bytes);

  // Convert dtype
  const dtype = descrToDtype(header.descr);

  // Calculate expected size
  const size = header.shape.length === 0 ? 1 : header.shape.reduce((a, b) => a * b, 1);
  const expectedBytes = size * dtypeSize(dtype);

  if (bytes.length - dataOffset < expectedBytes) {
    throw new Error(
      `NPY file truncated: expected ${expectedBytes} data bytes, ` +
      `got ${bytes.length - dataOffset}`
    );
  }

  // Create array from data
  const dataSlice = bytes.slice(dataOffset, dataOffset + expectedBytes);
  const typedArray = createTypedArrayFromBuffer(dtype, dataSlice.buffer, dataSlice.byteOffset, size);

  // Handle 0-d arrays (scalars)
  if (header.shape.length === 0) {
    const arr = await NDArray.zeros([1], { dtype });
    arr.setFlat(0, typedArray[0] as number);
    return arr;
  }

  const arr = await NDArray.fromTypedArray(typedArray, header.shape, dtype);

  // Handle Fortran order by setting appropriate strides
  if (header.fortran_order) {
    // Data is in Fortran order, transpose to get C-order view
    return arr.T;
  }

  return arr;
}

/**
 * Read file data from various sources.
 */
async function readFileData(
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

// Type declaration for FileSystemFileHandle (not always available)
interface FileSystemFileHandle {
  createWritable(): Promise<FileSystemWritableFileStream>;
  getFile(): Promise<File>;
}

interface FileSystemWritableFileStream {
  write(data: ArrayBuffer | Uint8Array): Promise<void>;
  seek?(position: number): Promise<void>;
  close(): Promise<void>;
}
