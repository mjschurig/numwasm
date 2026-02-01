/**
 * NPY/NPZ File I/O
 *
 * Functions for saving and loading NDArrays in NumPy's binary format.
 *
 * Reference: numpy/lib/_npyio_impl.py
 */

import { NDArray } from "../_core/NDArray.js";
import { DType } from "../types.js";
import {
  createHeader,
  parseHeader,
  dtypeToDescr,
  descrToDtype,
  dtypeSize,
  getMinVersion,
  isNode,
  type NpyHeader,
} from "./format.js";

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
  mmap_mode?: "r+" | "r" | "w+" | "c" | null;
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
  _options: SaveOptions = {},
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
  const shapeStr =
    header.shape.length === 0
      ? "()"
      : header.shape.length === 1
        ? `(${header.shape[0]},)`
        : `(${header.shape.join(", ")})`;
  const headerStrLen =
    `{'descr': '${header.descr}', 'fortran_order': False, 'shape': ${shapeStr}, }`
      .length;
  const version = getMinVersion(headerStrLen);

  // Create header bytes
  const headerBytes = createHeader(header, version);

  // Get array data as bytes
  const dataBytes = contiguous.toTypedArray();
  const dataView = new Uint8Array(
    dataBytes.buffer,
    dataBytes.byteOffset,
    contiguous.size * dtypeSize(contiguous.dtype),
  );

  // Combine header and data
  const totalSize = headerBytes.length + dataView.length;
  const result = new Uint8Array(totalSize);
  result.set(headerBytes, 0);
  result.set(dataView, headerBytes.length);

  // Write to destination
  if (file === null) {
    return result.buffer;
  }

  if (typeof file === "string") {
    // Node.js file path
    if (isNode()) {
      const fs = await import("fs/promises");
      await fs.writeFile(file, result);
      return;
    }
    throw new Error("String file paths only supported in Node.js");
  }

  // FileSystemFileHandle (File System Access API)
  if ("createWritable" in file) {
    const writable = await (file as FileSystemFileHandle).createWritable();
    await writable.write(result);
    await writable.close();
    return;
  }

  throw new Error("Unsupported file target");
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
  _options: LoadOptions = {},
): Promise<NDArray> {
  const data = await readFileData(file);
  const bytes = new Uint8Array(data);

  // Parse header
  const { header, dataOffset } = parseHeader(bytes);

  // Convert dtype
  const dtype = descrToDtype(header.descr);

  // Calculate expected size
  const size =
    header.shape.length === 0 ? 1 : header.shape.reduce((a, b) => a * b, 1);
  const expectedBytes = size * dtypeSize(dtype);

  if (bytes.length - dataOffset < expectedBytes) {
    throw new Error(
      `NPY file truncated: expected ${expectedBytes} data bytes, ` +
        `got ${bytes.length - dataOffset}`,
    );
  }

  // Create array from data
  const dataSlice = bytes.slice(dataOffset, dataOffset + expectedBytes);
  const typedArray = createTypedArrayFromBuffer(
    dtype,
    dataSlice.buffer,
    dataSlice.byteOffset,
    size,
  );

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
  file: string | File | ArrayBuffer | URL,
): Promise<ArrayBuffer> {
  if (file instanceof ArrayBuffer) {
    return file;
  }

  // Check for File in browser
  if (typeof File !== "undefined" && file instanceof File) {
    return file.arrayBuffer();
  }

  if (
    file instanceof URL ||
    (typeof file === "string" && file.startsWith("http"))
  ) {
    const response = await fetch(file.toString());
    return response.arrayBuffer();
  }

  if (typeof file === "string") {
    // Node.js file path
    if (isNode()) {
      const fs = await import("fs/promises");
      const buffer = await fs.readFile(file);
      return buffer.buffer.slice(
        buffer.byteOffset,
        buffer.byteOffset + buffer.byteLength,
      );
    }
    throw new Error("String file paths only supported in Node.js");
  }

  throw new Error("Unsupported file source");
}

/**
 * Create a typed array from a buffer based on dtype.
 */
function createTypedArrayFromBuffer(
  dtype: DType,
  buffer: ArrayBuffer,
  byteOffset: number,
  length: number,
):
  | Float64Array
  | Float32Array
  | Int32Array
  | Int16Array
  | Int8Array
  | Uint32Array
  | Uint16Array
  | Uint8Array
  | BigInt64Array
  | BigUint64Array {
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

/**
 * Options for savez()
 */
export interface SavezOptions {
  /** Compression level (0-9, only used for savez_compressed) */
  compressLevel?: number;
}

/**
 * Result type for loadz() - a dictionary-like object of arrays
 */
export interface NpzFile {
  /** Get array by name */
  [key: string]: NDArray | string[];
  /** Get all array names */
  files: string[];
}

/**
 * Save multiple arrays to a single file in uncompressed `.npz` format.
 *
 * The `.npz` file format is a ZIP archive containing multiple `.npy` files,
 * one for each array. Arrays are named 'arr_0', 'arr_1', etc. for positional
 * arguments, or by their keyword argument names.
 *
 * @param file - File path (Node.js), FileSystemFileHandle (browser), or null to return ArrayBuffer
 * @param args - Arrays to save (positional) or object mapping names to arrays
 * @returns ArrayBuffer if file is null, void otherwise
 *
 * @example
 * // Save with automatic names
 * await savez('data.npz', arr1, arr2);
 * // Creates data.npz with arr_0.npy, arr_1.npy
 *
 * // Save with custom names
 * await savez('data.npz', { x: arr1, y: arr2 });
 * // Creates data.npz with x.npy, y.npy
 *
 * // Get as ArrayBuffer
 * const buffer = await savez(null, arr1, arr2);
 */
export async function savez(
  file: string | FileSystemFileHandle | null,
  ...args: (NDArray | Record<string, NDArray>)[]
): Promise<ArrayBuffer | void> {
  return saveNpz(file, args, false);
}

/**
 * Save multiple arrays to a single file in compressed `.npz` format.
 *
 * Similar to `savez`, but each array is individually compressed using DEFLATE.
 * This can significantly reduce file size for sparse or repetitive data.
 *
 * @param file - File path (Node.js), FileSystemFileHandle (browser), or null to return ArrayBuffer
 * @param args - Arrays to save (positional) or object mapping names to arrays
 * @returns ArrayBuffer if file is null, void otherwise
 *
 * @example
 * // Save compressed with automatic names
 * await savez_compressed('data.npz', arr1, arr2);
 *
 * // Save compressed with custom names
 * await savez_compressed('data.npz', { weights: weights, biases: biases });
 */
export async function savez_compressed(
  file: string | FileSystemFileHandle | null,
  ...args: (NDArray | Record<string, NDArray>)[]
): Promise<ArrayBuffer | void> {
  return saveNpz(file, args, true);
}

/**
 * Internal function to save NPZ files.
 */
async function saveNpz(
  file: string | FileSystemFileHandle | null,
  args: (NDArray | Record<string, NDArray>)[],
  compress: boolean,
): Promise<ArrayBuffer | void> {
  // Collect named arrays
  const arrays: Map<string, NDArray> = new Map();
  let autoIndex = 0;

  for (const arg of args) {
    if (arg instanceof NDArray) {
      arrays.set(`arr_${autoIndex}`, arg);
      autoIndex++;
    } else {
      // It's a Record<string, NDArray>
      for (const [name, arr] of Object.entries(arg)) {
        if (arr instanceof NDArray) {
          arrays.set(name, arr);
        }
      }
    }
  }

  // Create ZIP archive
  const zipData = await createZipArchive(arrays, compress);

  // Write to destination
  if (file === null) {
    return zipData;
  }

  if (typeof file === "string") {
    if (isNode()) {
      const fs = await import("fs/promises");
      await fs.writeFile(file, new Uint8Array(zipData));
      return;
    }
    throw new Error("String file paths only supported in Node.js");
  }

  if ("createWritable" in file) {
    const writable = await (file as FileSystemFileHandle).createWritable();
    await writable.write(zipData);
    await writable.close();
    return;
  }

  throw new Error("Unsupported file target");
}

/**
 * Load arrays from a `.npz` file.
 *
 * Returns an object with array names as keys and NDArrays as values.
 * The `files` property contains the list of array names.
 *
 * @param file - File path (Node.js), File object (browser), ArrayBuffer, or URL
 * @returns Object containing named arrays
 *
 * @example
 * // Load all arrays
 * const data = await loadz('data.npz');
 * console.log(data.files);  // ['arr_0', 'arr_1']
 * const arr = data['arr_0'];
 *
 * // Load from URL
 * const data = await loadz('https://example.com/data.npz');
 */
export async function loadz(
  file: string | File | ArrayBuffer | URL,
): Promise<NpzFile> {
  const data = await readNpzData(file);
  const bytes = new Uint8Array(data);

  // Parse ZIP archive
  const entries = parseZipArchive(bytes);

  // Load each .npy file
  const result: NpzFile = { files: [] };

  for (const [name, content] of entries) {
    // Remove .npy extension from name
    const arrayName = name.endsWith('.npy') ? name.slice(0, -4) : name;
    result.files.push(arrayName);
    result[arrayName] = await loadNpyFromBuffer(content);
  }

  return result;
}

/**
 * Read NPZ file data from various sources.
 */
async function readNpzData(
  file: string | File | ArrayBuffer | URL,
): Promise<ArrayBuffer> {
  if (file instanceof ArrayBuffer) {
    return file;
  }

  if (typeof File !== "undefined" && file instanceof File) {
    return file.arrayBuffer();
  }

  if (
    file instanceof URL ||
    (typeof file === "string" && file.startsWith("http"))
  ) {
    const response = await fetch(file.toString());
    return response.arrayBuffer();
  }

  if (typeof file === "string") {
    if (isNode()) {
      const fs = await import("fs/promises");
      const buffer = await fs.readFile(file);
      return buffer.buffer.slice(
        buffer.byteOffset,
        buffer.byteOffset + buffer.byteLength,
      );
    }
    throw new Error("String file paths only supported in Node.js");
  }

  throw new Error("Unsupported file source");
}

/**
 * Load an NDArray from a buffer containing NPY data.
 */
async function loadNpyFromBuffer(bytes: Uint8Array): Promise<NDArray> {
  const { header, dataOffset } = parseHeader(bytes);
  const dtype = descrToDtype(header.descr);
  const size = header.shape.length === 0 ? 1 : header.shape.reduce((a, b) => a * b, 1);
  const expectedBytes = size * dtypeSize(dtype);

  if (bytes.length - dataOffset < expectedBytes) {
    throw new Error(
      `NPY data truncated: expected ${expectedBytes} data bytes, ` +
      `got ${bytes.length - dataOffset}`,
    );
  }

  const dataSlice = bytes.slice(dataOffset, dataOffset + expectedBytes);
  const typedArray = createTypedArrayFromNpyBuffer(
    dtype,
    dataSlice.buffer,
    dataSlice.byteOffset,
    size,
  );

  if (header.shape.length === 0) {
    const arr = await NDArray.zeros([1], { dtype });
    arr.setFlat(0, typedArray[0] as number);
    return arr;
  }

  const arr = await NDArray.fromTypedArray(typedArray, header.shape, dtype);
  return header.fortran_order ? arr.T : arr;
}

/**
 * Create a typed array from NPY buffer.
 */
function createTypedArrayFromNpyBuffer(
  dtype: DType,
  buffer: ArrayBuffer,
  byteOffset: number,
  length: number,
):
  | Float64Array
  | Float32Array
  | Int32Array
  | Int16Array
  | Int8Array
  | Uint32Array
  | Uint16Array
  | Uint8Array
  | BigInt64Array
  | BigUint64Array {
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
      return new Uint16Array(buffer, byteOffset, length);
    case DType.Complex64:
      return new Float32Array(buffer, byteOffset, length * 2);
    case DType.Complex128:
      return new Float64Array(buffer, byteOffset, length * 2);
    default:
      return new Float64Array(buffer, byteOffset, length);
  }
}

// =============================================================================
// ZIP Archive Implementation (minimal, for NPZ files)
// =============================================================================

// ZIP format constants
const ZIP_LOCAL_HEADER_SIG = 0x04034b50;
const ZIP_CENTRAL_HEADER_SIG = 0x02014b50;
const ZIP_END_OF_CENTRAL_SIG = 0x06054b50;
const ZIP_COMPRESSION_STORED = 0;
const ZIP_COMPRESSION_DEFLATE = 8;

/**
 * Create a ZIP archive from named arrays.
 */
async function createZipArchive(
  arrays: Map<string, NDArray>,
  compress: boolean,
): Promise<ArrayBuffer> {
  const localFiles: Uint8Array[] = [];
  const centralDirs: Uint8Array[] = [];
  let offset = 0;

  for (const [name, arr] of arrays) {
    const fileName = name + '.npy';
    const fileNameBytes = new TextEncoder().encode(fileName);

    // Get NPY data for this array
    const npyBuffer = await save(null, arr) as ArrayBuffer;
    let fileData = new Uint8Array(npyBuffer);

    let compressionMethod = ZIP_COMPRESSION_STORED;
    const uncompressedSize = fileData.length;
    let compressedSize = uncompressedSize;

    if (compress) {
      // Use DEFLATE compression if available
      if (typeof CompressionStream !== 'undefined') {
        fileData = await deflateData(fileData);
        compressionMethod = ZIP_COMPRESSION_DEFLATE;
        compressedSize = fileData.length;
      }
    }

    // Calculate CRC32
    const crc = crc32(new Uint8Array(npyBuffer));

    // Local file header (30 bytes + filename + data)
    const localHeader = new Uint8Array(30 + fileNameBytes.length);
    const localView = new DataView(localHeader.buffer);

    localView.setUint32(0, ZIP_LOCAL_HEADER_SIG, true);  // Signature
    localView.setUint16(4, 20, true);                     // Version needed
    localView.setUint16(6, 0, true);                      // Flags
    localView.setUint16(8, compressionMethod, true);      // Compression
    localView.setUint16(10, 0, true);                     // Mod time
    localView.setUint16(12, 0, true);                     // Mod date
    localView.setUint32(14, crc, true);                   // CRC32
    localView.setUint32(18, compressedSize, true);        // Compressed size
    localView.setUint32(22, uncompressedSize, true);      // Uncompressed size
    localView.setUint16(26, fileNameBytes.length, true);  // Filename length
    localView.setUint16(28, 0, true);                     // Extra field length
    localHeader.set(fileNameBytes, 30);

    // Central directory entry (46 bytes + filename)
    const centralDir = new Uint8Array(46 + fileNameBytes.length);
    const centralView = new DataView(centralDir.buffer);

    centralView.setUint32(0, ZIP_CENTRAL_HEADER_SIG, true);  // Signature
    centralView.setUint16(4, 20, true);                       // Version made by
    centralView.setUint16(6, 20, true);                       // Version needed
    centralView.setUint16(8, 0, true);                        // Flags
    centralView.setUint16(10, compressionMethod, true);       // Compression
    centralView.setUint16(12, 0, true);                       // Mod time
    centralView.setUint16(14, 0, true);                       // Mod date
    centralView.setUint32(16, crc, true);                     // CRC32
    centralView.setUint32(20, compressedSize, true);          // Compressed size
    centralView.setUint32(24, uncompressedSize, true);        // Uncompressed size
    centralView.setUint16(28, fileNameBytes.length, true);    // Filename length
    centralView.setUint16(30, 0, true);                       // Extra field length
    centralView.setUint16(32, 0, true);                       // Comment length
    centralView.setUint16(34, 0, true);                       // Disk number
    centralView.setUint16(36, 0, true);                       // Internal attrs
    centralView.setUint32(38, 0, true);                       // External attrs
    centralView.setUint32(42, offset, true);                  // Relative offset
    centralDir.set(fileNameBytes, 46);

    localFiles.push(localHeader);
    localFiles.push(fileData);
    centralDirs.push(centralDir);

    offset += localHeader.length + fileData.length;
  }

  // End of central directory (22 bytes)
  const centralDirOffset = offset;
  const centralDirSize = centralDirs.reduce((sum, cd) => sum + cd.length, 0);

  const endOfCentral = new Uint8Array(22);
  const endView = new DataView(endOfCentral.buffer);

  endView.setUint32(0, ZIP_END_OF_CENTRAL_SIG, true);    // Signature
  endView.setUint16(4, 0, true);                          // Disk number
  endView.setUint16(6, 0, true);                          // Central dir disk
  endView.setUint16(8, arrays.size, true);                // Entries on disk
  endView.setUint16(10, arrays.size, true);               // Total entries
  endView.setUint32(12, centralDirSize, true);            // Central dir size
  endView.setUint32(16, centralDirOffset, true);          // Central dir offset
  endView.setUint16(20, 0, true);                         // Comment length

  // Combine all parts
  const totalSize = offset + centralDirSize + 22;
  const result = new Uint8Array(totalSize);
  let pos = 0;

  for (const part of localFiles) {
    result.set(part, pos);
    pos += part.length;
  }
  for (const part of centralDirs) {
    result.set(part, pos);
    pos += part.length;
  }
  result.set(endOfCentral, pos);

  return result.buffer;
}

/**
 * Parse a ZIP archive and extract entries.
 */
function parseZipArchive(data: Uint8Array): Map<string, Uint8Array> {
  const entries = new Map<string, Uint8Array>();
  const view = new DataView(data.buffer, data.byteOffset, data.byteLength);

  // Find end of central directory
  let endOfCentralPos = -1;
  for (let i = data.length - 22; i >= 0; i--) {
    if (view.getUint32(i, true) === ZIP_END_OF_CENTRAL_SIG) {
      endOfCentralPos = i;
      break;
    }
  }

  if (endOfCentralPos < 0) {
    throw new Error('Invalid ZIP file: end of central directory not found');
  }

  const centralDirOffset = view.getUint32(endOfCentralPos + 16, true);
  const numEntries = view.getUint16(endOfCentralPos + 10, true);

  // Parse central directory entries
  let pos = centralDirOffset;
  for (let i = 0; i < numEntries; i++) {
    if (view.getUint32(pos, true) !== ZIP_CENTRAL_HEADER_SIG) {
      throw new Error('Invalid ZIP file: bad central directory signature');
    }

    const compressionMethod = view.getUint16(pos + 10, true);
    const compressedSize = view.getUint32(pos + 20, true);
    const uncompressedSize = view.getUint32(pos + 24, true);
    const fileNameLen = view.getUint16(pos + 28, true);
    const extraLen = view.getUint16(pos + 30, true);
    const commentLen = view.getUint16(pos + 32, true);
    const localHeaderOffset = view.getUint32(pos + 42, true);

    const fileName = new TextDecoder().decode(
      data.slice(pos + 46, pos + 46 + fileNameLen)
    );

    // Read file data from local header
    const localPos = localHeaderOffset;
    const localFileNameLen = view.getUint16(localPos + 26, true);
    const localExtraLen = view.getUint16(localPos + 28, true);
    const dataStart = localPos + 30 + localFileNameLen + localExtraLen;

    let fileData = data.slice(dataStart, dataStart + compressedSize);

    // Decompress if needed
    if (compressionMethod === ZIP_COMPRESSION_DEFLATE) {
      fileData = inflateData(fileData, uncompressedSize);
    } else if (compressionMethod !== ZIP_COMPRESSION_STORED) {
      throw new Error(`Unsupported compression method: ${compressionMethod}`);
    }

    entries.set(fileName, fileData);

    pos += 46 + fileNameLen + extraLen + commentLen;
  }

  return entries;
}

/**
 * Compress data using DEFLATE.
 */
async function deflateData(data: Uint8Array): Promise<Uint8Array> {
  if (typeof CompressionStream !== 'undefined') {
    const cs = new CompressionStream('deflate-raw');
    const writer = cs.writable.getWriter();
    writer.write(data);
    writer.close();

    const chunks: Uint8Array[] = [];
    const reader = cs.readable.getReader();

    while (true) {
      const { done, value } = await reader.read();
      if (done) break;
      chunks.push(value);
    }

    const totalLen = chunks.reduce((sum, c) => sum + c.length, 0);
    const result = new Uint8Array(totalLen);
    let offset = 0;
    for (const chunk of chunks) {
      result.set(chunk, offset);
      offset += chunk.length;
    }
    return result;
  }

  // Fallback: return uncompressed (caller should use STORED method)
  return data;
}

/**
 * Decompress DEFLATE data.
 */
function inflateData(data: Uint8Array, uncompressedSize: number): Uint8Array {
  if (typeof DecompressionStream !== 'undefined') {
    // Use synchronous pako-style inflate if available in browser
    // For now, use a simple synchronous DEFLATE implementation
    return inflateSync(data, uncompressedSize);
  }
  return inflateSync(data, uncompressedSize);
}

/**
 * Synchronous DEFLATE inflation.
 * This is a minimal implementation for NPZ compatibility.
 */
function inflateSync(data: Uint8Array, uncompressedSize: number): Uint8Array {
  // Use DecompressionStream if available (async, but we need sync)
  // For browser environments without pako, implement basic inflate

  // Try using Node.js zlib if available
  if (isNode()) {
    try {
      // Dynamic import won't work synchronously, so we use a workaround
      const zlibModule = require('zlib');
      const result = zlibModule.inflateRawSync(data);
      return new Uint8Array(result.buffer, result.byteOffset, result.byteLength);
    } catch {
      // Fall through to pure JS implementation
    }
  }

  // Pure JavaScript DEFLATE implementation
  return pureJsInflate(data, uncompressedSize);
}

/**
 * Pure JavaScript DEFLATE inflate implementation.
 * Handles both fixed and dynamic Huffman codes.
 */
function pureJsInflate(data: Uint8Array, expectedSize: number): Uint8Array {
  const output = new Uint8Array(expectedSize);
  let outPos = 0;
  let bitPos = 0;
  let bytePos = 0;

  function readBits(n: number): number {
    let result = 0;
    for (let i = 0; i < n; i++) {
      if (bytePos >= data.length) throw new Error('Unexpected end of data');
      result |= ((data[bytePos] >> bitPos) & 1) << i;
      bitPos++;
      if (bitPos === 8) {
        bitPos = 0;
        bytePos++;
      }
    }
    return result;
  }

  function readByte(): number {
    if (bytePos >= data.length) throw new Error('Unexpected end of data');
    return data[bytePos++];
  }

  // Fixed Huffman code tables
  const fixedLitLenLengths = new Uint8Array(288);
  for (let i = 0; i <= 143; i++) fixedLitLenLengths[i] = 8;
  for (let i = 144; i <= 255; i++) fixedLitLenLengths[i] = 9;
  for (let i = 256; i <= 279; i++) fixedLitLenLengths[i] = 7;
  for (let i = 280; i <= 287; i++) fixedLitLenLengths[i] = 8;

  const fixedDistLengths = new Uint8Array(32);
  for (let i = 0; i < 32; i++) fixedDistLengths[i] = 5;

  function buildHuffmanTable(lengths: Uint8Array): { codes: Uint16Array; lens: Uint8Array } {
    const maxLen = Math.max(...lengths);
    const count = new Uint16Array(maxLen + 1);
    for (let i = 0; i < lengths.length; i++) {
      if (lengths[i]) count[lengths[i]]++;
    }

    const nextCode = new Uint16Array(maxLen + 1);
    let code = 0;
    for (let bits = 1; bits <= maxLen; bits++) {
      code = (code + count[bits - 1]) << 1;
      nextCode[bits] = code;
    }

    const codes = new Uint16Array(lengths.length);
    const lens = new Uint8Array(lengths.length);
    for (let i = 0; i < lengths.length; i++) {
      if (lengths[i]) {
        codes[i] = nextCode[lengths[i]]++;
        lens[i] = lengths[i];
      }
    }

    return { codes, lens };
  }

  function decodeSymbol(table: { codes: Uint16Array; lens: Uint8Array }): number {
    let code = 0;
    let len = 0;

    while (len < 16) {
      code = (code << 1) | readBits(1);
      len++;

      for (let i = 0; i < table.codes.length; i++) {
        if (table.lens[i] === len && table.codes[i] === code) {
          return i;
        }
      }
    }

    throw new Error('Invalid Huffman code');
  }

  // Length and distance extra bits tables
  const lenBase = [3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 17, 19, 23, 27, 31, 35, 43, 51, 59, 67, 83, 99, 115, 131, 163, 195, 227, 258];
  const lenExtra = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 0];
  const distBase = [1, 2, 3, 4, 5, 7, 9, 13, 17, 25, 33, 49, 65, 97, 129, 193, 257, 385, 513, 769, 1025, 1537, 2049, 3073, 4097, 6145, 8193, 12289, 16385, 24577];
  const distExtra = [0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13];

  let bfinal = 0;

  while (!bfinal) {
    bfinal = readBits(1);
    const btype = readBits(2);

    if (btype === 0) {
      // Stored block
      bitPos = 0;
      bytePos++;
      const len = readByte() | (readByte() << 8);
      readByte(); readByte(); // nlen (complement, ignored)
      for (let i = 0; i < len; i++) {
        output[outPos++] = readByte();
      }
    } else if (btype === 1 || btype === 2) {
      let litLenTable: { codes: Uint16Array; lens: Uint8Array };
      let distTable: { codes: Uint16Array; lens: Uint8Array };

      if (btype === 1) {
        // Fixed Huffman
        litLenTable = buildHuffmanTable(fixedLitLenLengths);
        distTable = buildHuffmanTable(fixedDistLengths);
      } else {
        // Dynamic Huffman
        const hlit = readBits(5) + 257;
        const hdist = readBits(5) + 1;
        const hclen = readBits(4) + 4;

        const codeLenOrder = [16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15];
        const codeLenLengths = new Uint8Array(19);
        for (let i = 0; i < hclen; i++) {
          codeLenLengths[codeLenOrder[i]] = readBits(3);
        }

        const codeLenTable = buildHuffmanTable(codeLenLengths);
        const allLengths = new Uint8Array(hlit + hdist);
        let i = 0;

        while (i < hlit + hdist) {
          const sym = decodeSymbol(codeLenTable);
          if (sym < 16) {
            allLengths[i++] = sym;
          } else if (sym === 16) {
            const repeat = readBits(2) + 3;
            const val = allLengths[i - 1];
            for (let j = 0; j < repeat; j++) allLengths[i++] = val;
          } else if (sym === 17) {
            const repeat = readBits(3) + 3;
            for (let j = 0; j < repeat; j++) allLengths[i++] = 0;
          } else {
            const repeat = readBits(7) + 11;
            for (let j = 0; j < repeat; j++) allLengths[i++] = 0;
          }
        }

        litLenTable = buildHuffmanTable(allLengths.slice(0, hlit));
        distTable = buildHuffmanTable(allLengths.slice(hlit));
      }

      // Decode symbols
      while (true) {
        const sym = decodeSymbol(litLenTable);
        if (sym < 256) {
          output[outPos++] = sym;
        } else if (sym === 256) {
          break;
        } else {
          const lenIdx = sym - 257;
          const length = lenBase[lenIdx] + readBits(lenExtra[lenIdx]);

          const distSym = decodeSymbol(distTable);
          const distance = distBase[distSym] + readBits(distExtra[distSym]);

          for (let i = 0; i < length; i++) {
            output[outPos] = output[outPos - distance];
            outPos++;
          }
        }
      }
    } else {
      throw new Error('Invalid DEFLATE block type');
    }
  }

  return output.slice(0, outPos);
}

/**
 * CRC32 calculation for ZIP files.
 */
function crc32(data: Uint8Array): number {
  let crc = 0xffffffff;

  // Generate CRC table
  const table = new Uint32Array(256);
  for (let i = 0; i < 256; i++) {
    let c = i;
    for (let j = 0; j < 8; j++) {
      c = (c & 1) ? (0xedb88320 ^ (c >>> 1)) : (c >>> 1);
    }
    table[i] = c;
  }

  for (let i = 0; i < data.length; i++) {
    crc = table[(crc ^ data[i]) & 0xff] ^ (crc >>> 8);
  }

  return (crc ^ 0xffffffff) >>> 0;
}

/* ============ Bit Packing Functions ============ */

/**
 * Order of bits in packbits/unpackbits operations.
 */
export type BitOrder = "big" | "little";

/**
 * Pack the elements of a binary array into bits in a uint8 array.
 *
 * The result is padded to full bytes with zeros.
 *
 * @param a - Input array of uint8 values (each treated as a single bit: 0 or non-zero)
 * @param axis - The axis along which to pack. If null, the input is flattened
 * @param bitorder - Order of bits in the packed representation. 'big' (default) means
 *                   the first bit goes into the highest bit of the byte; 'little' means
 *                   the first bit goes into the lowest bit of the byte
 * @returns Array of uint8 with packed bits
 *
 * @example
 * // Pack 8 bits into a single byte
 * const a = await NDArray.fromArray([1, 0, 1, 1, 0, 0, 0, 1]);
 * const packed = await packbits(a);
 * // packed.toArray() -> [0b10110001] = [177] for big-endian
 *
 * @example
 * // Pack with little-endian bit order
 * const a = await NDArray.fromArray([1, 0, 1, 1, 0, 0, 0, 1]);
 * const packed = await packbits(a, null, 'little');
 * // First bit (1) goes to bit 0, so result differs
 */
export async function packbits(
  a: NDArray | number[],
  axis: number | null = null,
  bitorder: BitOrder = "big",
): Promise<NDArray> {
  // Convert to NDArray if needed
  let arr: NDArray;
  if (Array.isArray(a)) {
    arr = await NDArray.fromArray(a);
  } else {
    arr = a;
  }

  // Validate bitorder
  if (bitorder !== "big" && bitorder !== "little") {
    throw new Error(`bitorder must be 'big' or 'little', got: '${bitorder}'`);
  }

  // If axis is null, flatten the array
  if (axis === null) {
    const flatSize = arr.size;
    const nBytes = Math.ceil(flatSize / 8);
    const result = new Uint8Array(nBytes);

    for (let byteIdx = 0; byteIdx < nBytes; byteIdx++) {
      let byte = 0;
      for (let bitIdx = 0; bitIdx < 8; bitIdx++) {
        const flatIdx = byteIdx * 8 + bitIdx;
        if (flatIdx < flatSize) {
          const value = arr.getFlat(flatIdx);
          const bit = value !== 0 ? 1 : 0;
          if (bitorder === "big") {
            byte |= bit << (7 - bitIdx);
          } else {
            byte |= bit << bitIdx;
          }
        }
      }
      result[byteIdx] = byte;
    }

    return NDArray.fromArray(Array.from(result), undefined, { dtype: DType.Uint8 });
  }

  // Handle axis parameter
  const ndim = arr.ndim;
  if (axis < -ndim || axis >= ndim) {
    throw new Error(
      `axis ${axis} is out of bounds for array of dimension ${ndim}`,
    );
  }
  // Normalize negative axis
  const normalizedAxis = axis < 0 ? ndim + axis : axis;

  // Calculate output shape
  const inputShape = arr.shape;
  const axisSize = inputShape[normalizedAxis];
  const outputAxisSize = Math.ceil(axisSize / 8);

  const outputShape = [...inputShape];
  outputShape[normalizedAxis] = outputAxisSize;

  // Calculate strides for iteration
  const outputSize = outputShape.reduce((a, b) => a * b, 1);
  const result = new Uint8Array(outputSize);

  // Calculate the number of "slices" before and after the axis
  let preAxisSize = 1;
  let postAxisSize = 1;
  for (let i = 0; i < normalizedAxis; i++) {
    preAxisSize *= inputShape[i];
  }
  for (let i = normalizedAxis + 1; i < ndim; i++) {
    postAxisSize *= inputShape[i];
  }

  // Pack along the specified axis
  for (let pre = 0; pre < preAxisSize; pre++) {
    for (let post = 0; post < postAxisSize; post++) {
      for (let outByte = 0; outByte < outputAxisSize; outByte++) {
        let byte = 0;
        for (let bitIdx = 0; bitIdx < 8; bitIdx++) {
          const inAxisIdx = outByte * 8 + bitIdx;
          if (inAxisIdx < axisSize) {
            // Calculate flat index in input array
            const inputFlatIdx =
              pre * axisSize * postAxisSize + inAxisIdx * postAxisSize + post;
            const value = arr.getFlat(inputFlatIdx);
            const bit = value !== 0 ? 1 : 0;
            if (bitorder === "big") {
              byte |= bit << (7 - bitIdx);
            } else {
              byte |= bit << bitIdx;
            }
          }
        }
        // Calculate flat index in output array
        const outputFlatIdx =
          pre * outputAxisSize * postAxisSize +
          outByte * postAxisSize +
          post;
        result[outputFlatIdx] = byte;
      }
    }
  }

  return NDArray.fromArray(Array.from(result), outputShape, { dtype: DType.Uint8 });
}

/**
 * Unpack elements of a uint8 array into a binary array.
 *
 * Each element is unpacked into 8 bits (or fewer for the last element
 * if count is specified).
 *
 * @param a - Input array of uint8 values
 * @param axis - The axis along which to unpack. If null, the input is flattened
 * @param count - Number of bits to unpack from the last dimension. If null,
 *                all bits are unpacked. Negative count removes bits from the end.
 * @param bitorder - Order of bits in the packed representation. 'big' (default) means
 *                   the highest bit is unpacked first; 'little' means the lowest
 *                   bit is unpacked first
 * @returns Array of uint8 with unpacked bits (0 or 1)
 *
 * @example
 * // Unpack a single byte
 * const packed = await NDArray.fromArray([177]);  // 0b10110001
 * const unpacked = await unpackbits(packed);
 * // unpacked.toArray() -> [1, 0, 1, 1, 0, 0, 0, 1] for big-endian
 *
 * @example
 * // Unpack with count limit
 * const packed = await NDArray.fromArray([177]);
 * const unpacked = await unpackbits(packed, null, 4);
 * // unpacked.toArray() -> [1, 0, 1, 1] (only first 4 bits)
 */
export async function unpackbits(
  a: NDArray | number[],
  axis: number | null = null,
  count: number | null = null,
  bitorder: BitOrder = "big",
): Promise<NDArray> {
  // Convert to NDArray if needed
  let arr: NDArray;
  if (Array.isArray(a)) {
    arr = await NDArray.fromArray(a, undefined, { dtype: DType.Uint8 });
  } else {
    arr = a;
  }

  // Validate bitorder
  if (bitorder !== "big" && bitorder !== "little") {
    throw new Error(`bitorder must be 'big' or 'little', got: '${bitorder}'`);
  }

  // If axis is null, flatten the array
  if (axis === null) {
    const nBytes = arr.size;
    const totalBits = nBytes * 8;
    let nBits = totalBits;

    if (count !== null) {
      if (count < 0) {
        nBits = Math.max(0, totalBits + count);
      } else {
        nBits = Math.min(count, totalBits);
      }
    }

    const result = new Uint8Array(nBits);

    for (let bitIdx = 0; bitIdx < nBits; bitIdx++) {
      const byteIdx = Math.floor(bitIdx / 8);
      const bitInByte = bitIdx % 8;
      const byte = arr.getFlat(byteIdx);

      let bit: number;
      if (bitorder === "big") {
        bit = (byte >> (7 - bitInByte)) & 1;
      } else {
        bit = (byte >> bitInByte) & 1;
      }
      result[bitIdx] = bit;
    }

    return NDArray.fromArray(Array.from(result), undefined, { dtype: DType.Uint8 });
  }

  // Handle axis parameter
  const ndim = arr.ndim;
  if (axis < -ndim || axis >= ndim) {
    throw new Error(
      `axis ${axis} is out of bounds for array of dimension ${ndim}`,
    );
  }
  // Normalize negative axis
  const normalizedAxis = axis < 0 ? ndim + axis : axis;

  // Calculate output shape
  const inputShape = arr.shape;
  const axisSize = inputShape[normalizedAxis];
  let outputAxisSize = axisSize * 8;

  if (count !== null) {
    if (count < 0) {
      outputAxisSize = Math.max(0, outputAxisSize + count);
    } else {
      outputAxisSize = Math.min(count, outputAxisSize);
    }
  }

  const outputShape = [...inputShape];
  outputShape[normalizedAxis] = outputAxisSize;

  // Calculate strides for iteration
  const outputSize = outputShape.reduce((a, b) => a * b, 1);
  const result = new Uint8Array(outputSize);

  // Calculate the number of "slices" before and after the axis
  let preAxisSize = 1;
  let postAxisSize = 1;
  for (let i = 0; i < normalizedAxis; i++) {
    preAxisSize *= inputShape[i];
  }
  for (let i = normalizedAxis + 1; i < ndim; i++) {
    postAxisSize *= inputShape[i];
  }

  // Unpack along the specified axis
  for (let pre = 0; pre < preAxisSize; pre++) {
    for (let post = 0; post < postAxisSize; post++) {
      for (let outBit = 0; outBit < outputAxisSize; outBit++) {
        const inByteIdx = Math.floor(outBit / 8);
        const bitInByte = outBit % 8;

        // Calculate flat index in input array
        const inputFlatIdx =
          pre * axisSize * postAxisSize + inByteIdx * postAxisSize + post;
        const byte = arr.getFlat(inputFlatIdx);

        let bit: number;
        if (bitorder === "big") {
          bit = (byte >> (7 - bitInByte)) & 1;
        } else {
          bit = (byte >> bitInByte) & 1;
        }

        // Calculate flat index in output array
        const outputFlatIdx =
          pre * outputAxisSize * postAxisSize + outBit * postAxisSize + post;
        result[outputFlatIdx] = bit;
      }
    }
  }

  return NDArray.fromArray(Array.from(result), outputShape, { dtype: DType.Uint8 });
}
