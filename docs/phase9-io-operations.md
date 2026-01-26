# Phase 9: I/O Operations Implementation Plan

Phase 9 implements NumJS-WASM's input/output system for array persistence, text parsing, and string formatting. This enables saving/loading arrays to files, reading/writing text data, and customizing array display.

---

## Current State (Phases 1-8)

```
src/wasm/
├── ndarray.h/c        # Core array operations
├── dtype.h/c          # DType utilities
├── broadcast.h/c      # Broadcasting
├── indexing.h/c       # Index operations
└── pairwise_sum.h/c   # Accurate summation

src/ts/
├── types.ts           # DType enum, WasmModule interface
├── NDArray.ts         # Core NDArray class
├── dtype.ts           # Type promotion utilities
├── broadcast.ts       # Broadcasting functions
├── indexing.ts        # Index operations
├── slice.ts           # Slicing utilities
├── iterators.ts       # Iterator implementations
├── wasm-loader.ts     # WASM module loading
└── index.ts           # Public exports
```

**Prerequisites for Phase 9:**
- ✅ NDArray core with all dtypes
- ✅ Shape manipulation (reshape, transpose)
- ✅ Type promotion and casting
- ✅ Element access (get/set)
- ✅ Contiguity flags (C/F order)

---

## Phase 9 Implementation Tree

```
PHASE 9: I/O OPERATIONS
│
├── 9.1 Binary File Format (NPY/NPZ)
│   ├── 9.1.1 NPY Format Constants & Header
│   ├── 9.1.2 save(file, arr) → write .npy file
│   ├── 9.1.3 load(file) → read .npy file
│   ├── 9.1.4 savez(file, *arrays) → write .npz archive
│   ├── 9.1.5 savez_compressed(file, *arrays)
│   └── 9.1.6 NpzFile class (lazy dict-like)
│
│   Dependencies: NDArray core, dtype system
│
├── 9.2 Text File I/O
│   ├── 9.2.1 loadtxt(file, dtype, delimiter, ...)
│   ├── 9.2.2 savetxt(file, arr, fmt, delimiter, ...)
│   ├── 9.2.3 genfromtxt(file, ...) with missing values
│   └── 9.2.4 fromregex(file, regexp, dtype)
│
│   Dependencies: NDArray creation, dtype conversion
│
├── 9.3 Raw Binary I/O
│   ├── 9.3.1 fromfile(file, dtype, count, sep, offset)
│   ├── 9.3.2 NDArray.tofile(file, sep, format)
│   └── 9.3.3 frombuffer(buffer, dtype, count, offset)
│
│   Dependencies: TypedArray interop
│
├── 9.4 String Formatting
│   ├── 9.4.1 Print Options System
│   │   ├── set_printoptions(precision, threshold, ...)
│   │   ├── get_printoptions()
│   │   └── printoptions() context manager
│   ├── 9.4.2 array2string(arr, max_line_width, precision, ...)
│   ├── 9.4.3 array_repr(arr, ...) → full representation
│   ├── 9.4.4 array_str(arr, ...) → string representation
│   ├── 9.4.5 format_float_positional(x, precision, ...)
│   └── 9.4.6 format_float_scientific(x, precision, ...)
│
│   Dependencies: NDArray iteration
│
├── 9.5 Memory Mapping
│   ├── 9.5.1 memmap class (subclass of NDArray)
│   └── 9.5.2 open_memmap(filename, mode, dtype, shape)
│
│   Dependencies: NPY format, File API
│
└── 9.6 Base Conversion
    ├── 9.6.1 binary_repr(num, width)
    └── 9.6.2 base_repr(number, base, padding)

    Dependencies: None (standalone utilities)
```

---

## Detailed Implementation Specifications

### 9.1 Binary File Format (NPY/NPZ)

The NPY format is NumPy's native binary format for efficient array storage.

#### 9.1.1 NPY Format Constants & Header

**File:** `src/ts/io/format.ts` (new file)

```typescript
/**
 * NPY file format implementation.
 *
 * Format Structure:
 * - Magic string: \x93NUMPY (6 bytes)
 * - Version: major (1 byte) + minor (1 byte)
 * - Header length: 2 bytes (v1.0) or 4 bytes (v2.0+)
 * - Header: ASCII dict literal with 'descr', 'fortran_order', 'shape'
 * - Padding: spaces/newline to align to 64 bytes
 * - Data: raw array bytes
 */

/** Magic prefix for NPY files */
export const MAGIC_PREFIX = new Uint8Array([0x93, 0x4e, 0x55, 0x4d, 0x50, 0x59]); // \x93NUMPY

/** Total magic length including version bytes */
export const MAGIC_LEN = 8;

/** Alignment boundary for data section (enables memory mapping) */
export const ARRAY_ALIGN = 64;

/** Maximum header size for security (prevents eval attacks) */
export const MAX_HEADER_SIZE = 10000;

/** NPY format version */
export interface NpyVersion {
  major: number;
  minor: number;
}

/** NPY header dictionary */
export interface NpyHeader {
  descr: string;           // dtype descriptor (e.g., '<f8', '|b1')
  fortran_order: boolean;  // True if Fortran (column-major) order
  shape: number[];         // Array dimensions
}

/**
 * Determine minimum NPY version required for a dtype.
 * - v1.0: Basic types, header < 65536 bytes
 * - v2.0: Large headers (> 65535 bytes)
 * - v3.0: Unicode field names in structured dtypes
 */
export function getMinVersion(dtype: DType, headerSize: number): NpyVersion {
  if (headerSize > 65535) {
    return { major: 2, minor: 0 };
  }
  return { major: 1, minor: 0 };
}

/**
 * Convert DType to NPY descriptor string.
 *
 * Format: <endian><type><size>
 * - Endian: '<' (little), '>' (big), '|' (not applicable)
 * - Type: 'b' (int8), 'i' (int), 'u' (uint), 'f' (float), 'c' (complex), '?' (bool)
 * - Size: bytes (1, 2, 4, 8, 16)
 *
 * @example
 * dtypeToDescr(DType.Float64) // '<f8'
 * dtypeToDescr(DType.Int32)   // '<i4'
 * dtypeToDescr(DType.Bool)    // '|b1'
 */
export function dtypeToDescr(dtype: DType): string {
  const descriptors: Record<DType, string> = {
    [DType.Bool]: '|b1',
    [DType.Int8]: '|i1',
    [DType.Int16]: '<i2',
    [DType.Int32]: '<i4',
    [DType.Int64]: '<i8',
    [DType.UInt8]: '|u1',
    [DType.UInt16]: '<u2',
    [DType.UInt32]: '<u4',
    [DType.UInt64]: '<u8',
    [DType.Float16]: '<f2',
    [DType.Float32]: '<f4',
    [DType.Float64]: '<f8',
    [DType.Complex64]: '<c8',
    [DType.Complex128]: '<c16',
  };
  return descriptors[dtype] ?? '<f8';
}

/**
 * Parse NPY descriptor string to DType.
 */
export function descrToDtype(descr: string): DType {
  // Remove endianness for matching (all our types are native/little-endian)
  const normalized = descr.replace(/^[<>|=]/, '');

  const mapping: Record<string, DType> = {
    'b1': DType.Bool,
    '?': DType.Bool,
    'i1': DType.Int8,
    'i2': DType.Int16,
    'i4': DType.Int32,
    'i8': DType.Int64,
    'u1': DType.UInt8,
    'u2': DType.UInt16,
    'u4': DType.UInt32,
    'u8': DType.UInt64,
    'f2': DType.Float16,
    'f4': DType.Float32,
    'f8': DType.Float64,
    'c8': DType.Complex64,
    'c16': DType.Complex128,
  };

  const dtype = mapping[normalized];
  if (dtype === undefined) {
    throw new Error(`Unsupported dtype descriptor: ${descr}`);
  }
  return dtype;
}

/**
 * Create NPY header bytes.
 */
export function createHeader(header: NpyHeader, version: NpyVersion): Uint8Array {
  // Format header as Python dict literal
  const shapeStr = header.shape.length === 1
    ? `(${header.shape[0]},)`
    : `(${header.shape.join(', ')})`;

  const headerStr = `{'descr': '${header.descr}', 'fortran_order': ${header.fortran_order ? 'True' : 'False'}, 'shape': ${shapeStr}, }`;

  // Calculate padding for alignment
  const headerLenSize = version.major >= 2 ? 4 : 2;
  const baseLen = MAGIC_LEN + headerLenSize + headerStr.length + 1; // +1 for newline
  const padding = (ARRAY_ALIGN - (baseLen % ARRAY_ALIGN)) % ARRAY_ALIGN;
  const paddedHeader = headerStr + ' '.repeat(padding) + '\n';

  // Build complete header
  const headerLen = paddedHeader.length;
  const totalSize = MAGIC_LEN + headerLenSize + headerLen;
  const result = new Uint8Array(totalSize);

  // Magic prefix
  result.set(MAGIC_PREFIX, 0);

  // Version
  result[6] = version.major;
  result[7] = version.minor;

  // Header length (little-endian)
  const view = new DataView(result.buffer);
  if (version.major >= 2) {
    view.setUint32(8, headerLen, true);
  } else {
    view.setUint16(8, headerLen, true);
  }

  // Header string
  const encoder = new TextEncoder();
  const headerBytes = encoder.encode(paddedHeader);
  result.set(headerBytes, MAGIC_LEN + headerLenSize);

  return result;
}

/**
 * Parse NPY header from bytes.
 */
export function parseHeader(data: Uint8Array): { header: NpyHeader; dataOffset: number } {
  // Verify magic
  for (let i = 0; i < MAGIC_PREFIX.length; i++) {
    if (data[i] !== MAGIC_PREFIX[i]) {
      throw new Error('Invalid NPY file: bad magic number');
    }
  }

  const major = data[6];
  const minor = data[7];

  // Read header length
  const view = new DataView(data.buffer, data.byteOffset);
  let headerLen: number;
  let headerStart: number;

  if (major >= 2) {
    headerLen = view.getUint32(8, true);
    headerStart = 12;
  } else {
    headerLen = view.getUint16(8, true);
    headerStart = 10;
  }

  if (headerLen > MAX_HEADER_SIZE) {
    throw new Error(`NPY header too large: ${headerLen} > ${MAX_HEADER_SIZE}`);
  }

  // Parse header string
  const decoder = new TextDecoder('latin1');
  const headerStr = decoder.decode(data.slice(headerStart, headerStart + headerLen)).trim();

  // Parse Python dict literal (simplified parser for our subset)
  const header = parsePythonDict(headerStr);

  return {
    header,
    dataOffset: headerStart + headerLen,
  };
}

/**
 * Simple Python dict literal parser for NPY headers.
 * Only handles the specific format used in NPY files.
 */
function parsePythonDict(s: string): NpyHeader {
  // Extract 'descr' value
  const descrMatch = s.match(/'descr'\s*:\s*'([^']+)'/);
  if (!descrMatch) throw new Error('Invalid NPY header: missing descr');

  // Extract 'fortran_order' value
  const fortranMatch = s.match(/'fortran_order'\s*:\s*(True|False)/);
  if (!fortranMatch) throw new Error('Invalid NPY header: missing fortran_order');

  // Extract 'shape' value
  const shapeMatch = s.match(/'shape'\s*:\s*\(([^)]*)\)/);
  if (!shapeMatch) throw new Error('Invalid NPY header: missing shape');

  const shape = shapeMatch[1]
    .split(',')
    .map(x => x.trim())
    .filter(x => x.length > 0)
    .map(x => parseInt(x, 10));

  return {
    descr: descrMatch[1],
    fortran_order: fortranMatch[1] === 'True',
    shape,
  };
}
```

---

#### 9.1.2 save() - Write Single Array

**File:** `src/ts/io/npyio.ts` (new file)

```typescript
import { NDArray } from '../NDArray.js';
import { DType, dtypeSize } from '../types.js';
import {
  createHeader,
  dtypeToDescr,
  getMinVersion,
  NpyHeader
} from './format.js';

/**
 * Save an array to a binary file in NPY format.
 *
 * @param file - File path, File object, or ArrayBuffer to write to
 * @param arr - Array to save
 * @param options - Save options
 *
 * @example
 * // In Node.js
 * await save('data.npy', arr);
 *
 * // In browser (returns ArrayBuffer)
 * const buffer = save(null, arr);
 *
 * // With File System Access API
 * const handle = await showSaveFilePicker();
 * await save(handle, arr);
 */
export async function save(
  file: string | FileSystemFileHandle | null,
  arr: NDArray,
  options: SaveOptions = {}
): Promise<ArrayBuffer | void> {
  // Ensure array is contiguous for efficient writing
  const contiguous = arr.flags.c_contiguous ? arr : arr.copy();

  // Build header
  const header: NpyHeader = {
    descr: dtypeToDescr(contiguous.dtype),
    fortran_order: false,
    shape: contiguous.shape,
  };

  // Determine version
  const headerStr = JSON.stringify(header);
  const version = getMinVersion(contiguous.dtype, headerStr.length);

  // Create header bytes
  const headerBytes = createHeader(header, version);

  // Get array data as bytes
  const dataBytes = new Uint8Array(
    contiguous.data.buffer,
    contiguous.data.byteOffset,
    contiguous.size * dtypeSize(contiguous.dtype)
  );

  // Combine header and data
  const totalSize = headerBytes.length + dataBytes.length;
  const result = new Uint8Array(totalSize);
  result.set(headerBytes, 0);
  result.set(dataBytes, headerBytes.length);

  // Write to destination
  if (file === null) {
    return result.buffer;
  }

  if (typeof file === 'string') {
    // Node.js file path
    if (typeof process !== 'undefined' && process.versions?.node) {
      const fs = await import('fs/promises');
      await fs.writeFile(file, result);
      return;
    }
    throw new Error('String file paths only supported in Node.js');
  }

  // FileSystemFileHandle (File System Access API)
  if ('createWritable' in file) {
    const writable = await file.createWritable();
    await writable.write(result);
    await writable.close();
    return;
  }

  throw new Error('Unsupported file target');
}

export interface SaveOptions {
  /** Allow saving object arrays (not supported in NumJS) */
  allow_pickle?: boolean;
}
```

---

#### 9.1.3 load() - Read Array

```typescript
/**
 * Load an array from a NPY file.
 *
 * @param file - File path, File object, ArrayBuffer, or URL
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
  options: LoadOptions = {}
): Promise<NDArray> {
  const data = await readFileData(file);
  const bytes = new Uint8Array(data);

  // Parse header
  const { header, dataOffset } = parseHeader(bytes);

  // Convert dtype
  const dtype = descrToDtype(header.descr);

  // Calculate expected size
  const size = header.shape.reduce((a, b) => a * b, 1);
  const expectedBytes = size * dtypeSize(dtype);

  if (bytes.length - dataOffset < expectedBytes) {
    throw new Error(
      `NPY file truncated: expected ${expectedBytes} data bytes, ` +
      `got ${bytes.length - dataOffset}`
    );
  }

  // Create array from data
  const dataSlice = bytes.slice(dataOffset, dataOffset + expectedBytes);
  const arr = NDArray.fromTypedArray(
    createTypedArray(dtype, dataSlice.buffer, dataSlice.byteOffset, size),
    header.shape,
    dtype
  );

  // Handle Fortran order
  if (header.fortran_order) {
    // Data is in Fortran order, need to set strides appropriately
    return arr.T; // Transpose to get correct view
  }

  return arr;
}

export interface LoadOptions {
  /** Memory mapping mode (not yet supported) */
  mmap_mode?: 'r+' | 'r' | 'w+' | 'c' | null;
  /** Allow loading pickled objects (not supported) */
  allow_pickle?: boolean;
  /** Maximum header size */
  max_header_size?: number;
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

  if (file instanceof File) {
    return file.arrayBuffer();
  }

  if (file instanceof URL || (typeof file === 'string' && file.startsWith('http'))) {
    const response = await fetch(file.toString());
    return response.arrayBuffer();
  }

  if (typeof file === 'string') {
    // Node.js file path
    if (typeof process !== 'undefined' && process.versions?.node) {
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
```

---

#### 9.1.4-9.1.5 savez() and savez_compressed()

```typescript
/**
 * Save multiple arrays to an uncompressed NPZ archive.
 *
 * @param file - Output file
 * @param arrays - Arrays as positional args (arr_0, arr_1, ...) or named
 *
 * @example
 * // Positional arrays
 * await savez('data.npz', arr1, arr2, arr3);
 * // Creates: arr_0.npy, arr_1.npy, arr_2.npy
 *
 * // Named arrays
 * await savez('data.npz', { x: arr1, y: arr2 });
 * // Creates: x.npy, y.npy
 */
export async function savez(
  file: string | FileSystemFileHandle | null,
  ...arrays: (NDArray | Record<string, NDArray>)[]
): Promise<ArrayBuffer | void> {
  return _savez(file, arrays, false);
}

/**
 * Save multiple arrays to a compressed NPZ archive.
 */
export async function savez_compressed(
  file: string | FileSystemFileHandle | null,
  ...arrays: (NDArray | Record<string, NDArray>)[]
): Promise<ArrayBuffer | void> {
  return _savez(file, arrays, true);
}

async function _savez(
  file: string | FileSystemFileHandle | null,
  arrays: (NDArray | Record<string, NDArray>)[],
  compress: boolean
): Promise<ArrayBuffer | void> {
  // Normalize arrays to Record<string, NDArray>
  const namedArrays: Record<string, NDArray> = {};
  let autoIndex = 0;

  for (const item of arrays) {
    if (item instanceof NDArray) {
      namedArrays[`arr_${autoIndex++}`] = item;
    } else {
      Object.assign(namedArrays, item);
    }
  }

  // Use JSZip or similar library for ZIP creation
  // (dependency would need to be added)
  const zip = new JSZip();

  for (const [name, arr] of Object.entries(namedArrays)) {
    const npyData = await save(null, arr) as ArrayBuffer;
    zip.file(`${name}.npy`, npyData, {
      compression: compress ? 'DEFLATE' : 'STORE',
    });
  }

  const zipData = await zip.generateAsync({ type: 'arraybuffer' });

  if (file === null) {
    return zipData;
  }

  // Write to file (same logic as save())
  // ...
}
```

---

#### 9.1.6 NpzFile Class

```typescript
/**
 * Lazy-loading dictionary-like container for NPZ files.
 * Arrays are loaded on-demand when accessed.
 *
 * @example
 * const npz = await NpzFile.open('data.npz');
 * console.log(npz.files); // ['arr_0', 'arr_1', 'x', 'y']
 *
 * const x = npz.get('x');  // Loads x.npy
 * const y = npz.get('y');  // Loads y.npy
 *
 * npz.close();
 */
export class NpzFile implements Iterable<[string, NDArray]> {
  private _zip: JSZip;
  private _cache: Map<string, NDArray> = new Map();
  private _files: string[];

  private constructor(zip: JSZip) {
    this._zip = zip;
    this._files = Object.keys(zip.files)
      .filter(name => name.endsWith('.npy'))
      .map(name => name.slice(0, -4)); // Remove .npy extension
  }

  /**
   * Open an NPZ file.
   */
  static async open(
    file: string | File | ArrayBuffer | URL
  ): Promise<NpzFile> {
    const data = await readFileData(file);
    const zip = await JSZip.loadAsync(data);
    return new NpzFile(zip);
  }

  /**
   * List of array names in the archive.
   */
  get files(): string[] {
    return [...this._files];
  }

  /**
   * Get an array by name.
   * Loads from archive on first access, then cached.
   */
  async get(name: string): Promise<NDArray> {
    if (this._cache.has(name)) {
      return this._cache.get(name)!;
    }

    const npyFile = this._zip.file(`${name}.npy`);
    if (!npyFile) {
      throw new Error(`Array '${name}' not found in NPZ archive`);
    }

    const data = await npyFile.async('arraybuffer');
    const arr = await load(data);
    this._cache.set(name, arr);
    return arr;
  }

  /**
   * Check if array exists.
   */
  has(name: string): boolean {
    return this._files.includes(name);
  }

  /**
   * Get array names.
   */
  keys(): IterableIterator<string> {
    return this._files.values();
  }

  /**
   * Iterate over all arrays (loads all into memory).
   */
  async *[Symbol.asyncIterator](): AsyncIterableIterator<[string, NDArray]> {
    for (const name of this._files) {
      yield [name, await this.get(name)];
    }
  }

  /**
   * Close and release resources.
   */
  close(): void {
    this._cache.clear();
    // @ts-ignore - clear internal state
    this._zip = null;
  }
}
```

---

### 9.2 Text File I/O

#### 9.2.1 loadtxt()

**File:** `src/ts/io/text.ts` (new file)

```typescript
import { NDArray } from '../NDArray.js';
import { DType } from '../types.js';

/**
 * Load data from a text file.
 *
 * Each row in the text file must have the same number of values.
 *
 * @param file - File source
 * @param options - Parsing options
 * @returns Loaded array
 *
 * @example
 * // Simple CSV
 * const arr = await loadtxt('data.csv', { delimiter: ',' });
 *
 * // Skip header row
 * const arr = await loadtxt('data.txt', { skiprows: 1 });
 *
 * // Select specific columns
 * const arr = await loadtxt('data.txt', { usecols: [0, 2, 3] });
 *
 * // With custom converter
 * const arr = await loadtxt('data.txt', {
 *   converters: { 0: (s) => s === 'yes' ? 1 : 0 }
 * });
 */
export async function loadtxt(
  file: string | File | URL,
  options: LoadtxtOptions = {}
): Promise<NDArray> {
  const {
    dtype = DType.Float64,
    comments = '#',
    delimiter,  // undefined = whitespace
    converters = {},
    skiprows = 0,
    usecols,
    unpack = false,
    ndmin = 0,
    encoding = 'utf-8',
    max_rows,
  } = options;

  // Read file content
  const text = await readTextFile(file, encoding);
  const lines = text.split(/\r?\n/);

  // Parse data
  const data: number[][] = [];
  let rowCount = 0;

  for (let i = skiprows; i < lines.length; i++) {
    if (max_rows !== undefined && rowCount >= max_rows) break;

    let line = lines[i].trim();

    // Skip empty lines and comments
    if (line.length === 0) continue;
    if (comments && line.startsWith(comments)) continue;

    // Remove inline comments
    if (comments) {
      const commentIdx = line.indexOf(comments);
      if (commentIdx !== -1) {
        line = line.slice(0, commentIdx).trim();
      }
    }

    // Split by delimiter
    const parts = delimiter
      ? line.split(delimiter)
      : line.split(/\s+/);

    // Select columns
    const selectedParts = usecols
      ? usecols.map(i => parts[i])
      : parts;

    // Convert values
    const row = selectedParts.map((val, colIdx) => {
      const actualCol = usecols ? usecols[colIdx] : colIdx;
      if (converters[actualCol]) {
        return converters[actualCol](val);
      }
      return parseFloat(val);
    });

    data.push(row);
    rowCount++;
  }

  if (data.length === 0) {
    throw new Error('No data found in file');
  }

  // Create array
  let arr = NDArray.fromArray(data, dtype);

  // Handle ndmin
  while (arr.ndim < ndmin) {
    arr = arr.reshape([1, ...arr.shape]);
  }

  // Handle unpack (transpose)
  if (unpack && arr.ndim >= 2) {
    arr = arr.T;
  }

  return arr;
}

export interface LoadtxtOptions {
  /** Output data type */
  dtype?: DType;
  /** Characters indicating start of comment */
  comments?: string | null;
  /** Column delimiter (whitespace if undefined) */
  delimiter?: string;
  /** Column converters: { colIndex: (str) => number } */
  converters?: Record<number, (val: string) => number>;
  /** Number of rows to skip at beginning */
  skiprows?: number;
  /** Columns to read (0-indexed) */
  usecols?: number[];
  /** Transpose result (for unpacking columns) */
  unpack?: boolean;
  /** Minimum number of dimensions */
  ndmin?: 0 | 1 | 2;
  /** File encoding */
  encoding?: string;
  /** Maximum rows to read */
  max_rows?: number;
}
```

---

#### 9.2.2 savetxt()

```typescript
/**
 * Save an array to a text file.
 *
 * @param file - Output file
 * @param arr - 1D or 2D array to save
 * @param options - Formatting options
 *
 * @example
 * // Default scientific notation
 * await savetxt('data.txt', arr);
 *
 * // CSV with header
 * await savetxt('data.csv', arr, {
 *   delimiter: ',',
 *   header: 'x,y,z',
 * });
 *
 * // Custom format per column
 * await savetxt('data.txt', arr, {
 *   fmt: ['%.3f', '%.6e', '%d'],
 * });
 */
export async function savetxt(
  file: string | FileSystemFileHandle,
  arr: NDArray,
  options: SavetxtOptions = {}
): Promise<void> {
  const {
    fmt = '%.18e',
    delimiter = ' ',
    newline = '\n',
    header = '',
    footer = '',
    comments = '# ',
    encoding = 'utf-8',
  } = options;

  if (arr.ndim > 2) {
    throw new Error('savetxt only supports 1D and 2D arrays');
  }

  // Ensure 2D
  const data = arr.ndim === 1 ? arr.reshape([arr.size, 1]) : arr;
  const [nrows, ncols] = data.shape;

  // Parse format specifiers
  const formats = Array.isArray(fmt) ? fmt : Array(ncols).fill(fmt);
  if (formats.length !== ncols) {
    throw new Error(`Number of formats (${formats.length}) doesn't match columns (${ncols})`);
  }

  // Build output
  const lines: string[] = [];

  // Header
  if (header) {
    const headerLines = header.split('\n');
    for (const line of headerLines) {
      lines.push(comments + line);
    }
  }

  // Data rows
  for (let i = 0; i < nrows; i++) {
    const row: string[] = [];
    for (let j = 0; j < ncols; j++) {
      const value = data.get(i, j);
      row.push(formatValue(value, formats[j]));
    }
    lines.push(row.join(delimiter));
  }

  // Footer
  if (footer) {
    const footerLines = footer.split('\n');
    for (const line of footerLines) {
      lines.push(comments + line);
    }
  }

  const text = lines.join(newline) + newline;

  // Write to file
  await writeTextFile(file, text, encoding);
}

export interface SavetxtOptions {
  /** Format string(s) for values */
  fmt?: string | string[];
  /** Column delimiter */
  delimiter?: string;
  /** Line separator */
  newline?: string;
  /** Header string (prepended with comments) */
  header?: string;
  /** Footer string (prepended with comments) */
  footer?: string;
  /** Comment prefix for header/footer */
  comments?: string;
  /** Output encoding */
  encoding?: string;
}

/**
 * Format a number according to printf-style format string.
 */
function formatValue(value: number, fmt: string): string {
  // Parse format: %[flags][width][.precision]specifier
  const match = fmt.match(/^%([+\- 0#]*)(\d*)(?:\.(\d+))?([diouxXeEfFgGaAcs%])$/);
  if (!match) {
    throw new Error(`Invalid format string: ${fmt}`);
  }

  const [, flags, width, precision, specifier] = match;
  const prec = precision ? parseInt(precision, 10) : 6;

  let result: string;

  switch (specifier) {
    case 'd':
    case 'i':
      result = Math.trunc(value).toString();
      break;
    case 'e':
      result = value.toExponential(prec);
      break;
    case 'E':
      result = value.toExponential(prec).toUpperCase();
      break;
    case 'f':
    case 'F':
      result = value.toFixed(prec);
      break;
    case 'g':
      result = value.toPrecision(prec || 1);
      break;
    case 'G':
      result = value.toPrecision(prec || 1).toUpperCase();
      break;
    default:
      result = value.toString();
  }

  // Apply width padding
  if (width) {
    const w = parseInt(width, 10);
    const pad = flags.includes('-') ? 'end' : 'start';
    const char = flags.includes('0') ? '0' : ' ';
    result = pad === 'start'
      ? result.padStart(w, char)
      : result.padEnd(w, char);
  }

  // Apply sign
  if (flags.includes('+') && value >= 0 && !result.startsWith('-')) {
    result = '+' + result;
  } else if (flags.includes(' ') && value >= 0 && !result.startsWith('-')) {
    result = ' ' + result;
  }

  return result;
}
```

---

#### 9.2.3 genfromtxt()

```typescript
/**
 * Load data from a text file with handling for missing values.
 * More sophisticated than loadtxt, supports structured arrays.
 *
 * @param file - File source
 * @param options - Parsing options
 * @returns Loaded array
 *
 * @example
 * // Handle missing values
 * const arr = await genfromtxt('data.csv', {
 *   delimiter: ',',
 *   missing_values: ['', 'NA', 'N/A'],
 *   filling_values: NaN,
 * });
 *
 * // Auto-detect column names from header
 * const arr = await genfromtxt('data.csv', {
 *   delimiter: ',',
 *   names: true,
 * });
 */
export async function genfromtxt(
  file: string | File | URL,
  options: GenfromtxtOptions = {}
): Promise<NDArray> {
  const {
    dtype = DType.Float64,
    comments = '#',
    delimiter,
    skip_header = 0,
    skip_footer = 0,
    converters = {},
    missing_values = [''],
    filling_values = NaN,
    usecols,
    names,
    encoding = 'utf-8',
    max_rows,
    invalid_raise = true,
  } = options;

  // Read and parse similar to loadtxt but with:
  // - Missing value detection and replacement
  // - Auto dtype detection if dtype=null
  // - Column name extraction from header
  // - Skip footer support

  const text = await readTextFile(file, encoding);
  const allLines = text.split(/\r?\n/);

  // Remove footer lines
  const lines = skip_footer > 0
    ? allLines.slice(0, -skip_footer)
    : allLines;

  // Handle names
  let columnNames: string[] | null = null;
  let dataStart = skip_header;

  if (names === true) {
    // Use first non-comment line as names
    for (let i = skip_header; i < lines.length; i++) {
      const line = lines[i].trim();
      if (line && !line.startsWith(comments || '')) {
        columnNames = delimiter
          ? line.split(delimiter).map(s => s.trim())
          : line.split(/\s+/);
        dataStart = i + 1;
        break;
      }
    }
  } else if (Array.isArray(names)) {
    columnNames = names;
  }

  // Parse data with missing value handling
  const missingSet = new Set(missing_values);
  const data: number[][] = [];
  let rowCount = 0;

  for (let i = dataStart; i < lines.length; i++) {
    if (max_rows !== undefined && rowCount >= max_rows) break;

    let line = lines[i].trim();
    if (!line || (comments && line.startsWith(comments))) continue;

    const parts = delimiter
      ? line.split(delimiter)
      : line.split(/\s+/);

    const selectedParts = usecols
      ? usecols.map(i => parts[i])
      : parts;

    const row = selectedParts.map((val, colIdx) => {
      const trimmed = val.trim();

      // Check for missing
      if (missingSet.has(trimmed)) {
        return typeof filling_values === 'number'
          ? filling_values
          : (filling_values[colIdx] ?? NaN);
      }

      // Apply converter
      const actualCol = usecols ? usecols[colIdx] : colIdx;
      if (converters[actualCol]) {
        return converters[actualCol](trimmed);
      }

      const num = parseFloat(trimmed);
      if (isNaN(num) && invalid_raise) {
        throw new Error(`Cannot convert '${trimmed}' to float at row ${i}, column ${colIdx}`);
      }
      return num;
    });

    data.push(row);
    rowCount++;
  }

  // TODO: Structured array support when names are provided
  // For now, return regular array
  return NDArray.fromArray(data, dtype);
}

export interface GenfromtxtOptions extends LoadtxtOptions {
  /** Rows to skip at end of file */
  skip_footer?: number;
  /** Skip header rows (replaces skiprows) */
  skip_header?: number;
  /** Strings treated as missing values */
  missing_values?: string[];
  /** Values to use for missing data */
  filling_values?: number | number[];
  /** Column names: true=from first row, array=explicit names */
  names?: boolean | string[];
  /** Raise exception on invalid values */
  invalid_raise?: boolean;
}
```

---

#### 9.2.4 fromregex()

```typescript
/**
 * Construct an array from a text file using regular expression parsing.
 *
 * @param file - File source
 * @param regexp - Regular expression with groups to extract
 * @param dtype - Output data type
 * @param options - Additional options
 *
 * @example
 * // Extract numbers from log file
 * const arr = await fromregex('log.txt', /value=(\d+\.?\d*)/g, DType.Float64);
 *
 * // Multiple groups
 * const arr = await fromregex('data.txt', /(\d+)\s+(\d+\.?\d*)/g, DType.Float64);
 */
export async function fromregex(
  file: string | File | URL,
  regexp: RegExp,
  dtype: DType = DType.Float64,
  options: FromregexOptions = {}
): Promise<NDArray> {
  const { encoding = 'utf-8' } = options;

  const text = await readTextFile(file, encoding);

  // Ensure global flag for matchAll
  const globalRegexp = regexp.global
    ? regexp
    : new RegExp(regexp.source, regexp.flags + 'g');

  const matches = [...text.matchAll(globalRegexp)];

  if (matches.length === 0) {
    throw new Error('No matches found');
  }

  // Extract captured groups
  const numGroups = matches[0].length - 1;

  if (numGroups === 0) {
    // No groups, use full match
    const data = matches.map(m => parseFloat(m[0]));
    return NDArray.fromArray(data, dtype);
  }

  if (numGroups === 1) {
    // Single group, return 1D array
    const data = matches.map(m => parseFloat(m[1]));
    return NDArray.fromArray(data, dtype);
  }

  // Multiple groups, return 2D array
  const data = matches.map(m => {
    const row: number[] = [];
    for (let i = 1; i <= numGroups; i++) {
      row.push(parseFloat(m[i]));
    }
    return row;
  });

  return NDArray.fromArray(data, dtype);
}

export interface FromregexOptions {
  /** File encoding */
  encoding?: string;
}
```

---

### 9.3 Raw Binary I/O

**File:** `src/ts/io/binary.ts` (new file)

```typescript
import { NDArray } from '../NDArray.js';
import { DType, dtypeSize, createTypedArray } from '../types.js';

/**
 * Construct an array from raw binary data in a file.
 *
 * @param file - File source
 * @param dtype - Data type to interpret as
 * @param options - Read options
 *
 * @example
 * // Read entire file as float64
 * const arr = await fromfile('data.bin', DType.Float64);
 *
 * // Read specific count with offset
 * const arr = await fromfile('data.bin', DType.Int32, {
 *   count: 100,
 *   offset: 256,
 * });
 *
 * // Read text file with separator
 * const arr = await fromfile('numbers.txt', DType.Float64, {
 *   sep: ' ',
 * });
 */
export async function fromfile(
  file: string | File | ArrayBuffer,
  dtype: DType = DType.Float64,
  options: FromfileOptions = {}
): Promise<NDArray> {
  const {
    count = -1,
    sep = '',
    offset = 0,
  } = options;

  if (sep) {
    // Text mode: parse separated values
    const text = await readTextFile(file as string | File, 'utf-8');
    const values = text.split(sep)
      .map(s => s.trim())
      .filter(s => s.length > 0)
      .map(s => parseFloat(s));

    const finalValues = count > 0 ? values.slice(0, count) : values;
    return NDArray.fromArray(finalValues, dtype);
  }

  // Binary mode
  const buffer = await readFileBuffer(file);
  const itemSize = dtypeSize(dtype);
  const startByte = offset;
  const availableItems = Math.floor((buffer.byteLength - startByte) / itemSize);
  const numItems = count > 0 ? Math.min(count, availableItems) : availableItems;

  const typedArray = createTypedArray(dtype, buffer, startByte, numItems);
  return NDArray.fromTypedArray(typedArray, [numItems], dtype);
}

export interface FromfileOptions {
  /** Number of items to read (-1 for all) */
  count?: number;
  /** Separator for text mode (empty string for binary) */
  sep?: string;
  /** Byte offset in file */
  offset?: number;
}

/**
 * Construct an array from a buffer object.
 *
 * @param buffer - Buffer-like object
 * @param dtype - Data type
 * @param options - Read options
 *
 * @example
 * const buf = new ArrayBuffer(32);
 * new Float64Array(buf).set([1, 2, 3, 4]);
 * const arr = frombuffer(buf, DType.Float64);
 */
export function frombuffer(
  buffer: ArrayBuffer | ArrayBufferView,
  dtype: DType = DType.Float64,
  options: FrombufferOptions = {}
): NDArray {
  const {
    count = -1,
    offset = 0,
  } = options;

  const arrayBuffer = buffer instanceof ArrayBuffer
    ? buffer
    : buffer.buffer;
  const bufferOffset = buffer instanceof ArrayBuffer
    ? 0
    : buffer.byteOffset;

  const itemSize = dtypeSize(dtype);
  const startByte = bufferOffset + offset;
  const availableItems = Math.floor((arrayBuffer.byteLength - startByte) / itemSize);
  const numItems = count > 0 ? Math.min(count, availableItems) : availableItems;

  const typedArray = createTypedArray(dtype, arrayBuffer, startByte, numItems);
  return NDArray.fromTypedArray(typedArray, [numItems], dtype);
}

export interface FrombufferOptions {
  /** Number of items to read */
  count?: number;
  /** Byte offset */
  offset?: number;
}
```

**NDArray.tofile() method:**

Add to `src/ts/NDArray.ts`:

```typescript
/**
 * Write array to a file as binary or text.
 *
 * @param file - Output file
 * @param options - Write options
 *
 * @example
 * // Binary output
 * await arr.tofile('data.bin');
 *
 * // Text output with separator
 * await arr.tofile('data.txt', { sep: ' ' });
 *
 * // Formatted text
 * await arr.tofile('data.txt', { sep: '\n', format: '%.6f' });
 */
async tofile(
  file: string | FileSystemFileHandle,
  options: TofileOptions = {}
): Promise<void> {
  const { sep = '', format = '' } = options;

  if (sep) {
    // Text mode
    const values = this.toArray();
    const formatted = format
      ? values.map(v => formatValue(v, format))
      : values.map(v => v.toString());
    const text = formatted.join(sep);
    await writeTextFile(file, text, 'utf-8');
  } else {
    // Binary mode
    const contiguous = this.flags.c_contiguous ? this : this.copy();
    const bytes = new Uint8Array(
      contiguous.data.buffer,
      contiguous.data.byteOffset,
      contiguous.size * dtypeSize(this.dtype)
    );
    await writeBinaryFile(file, bytes);
  }
}

interface TofileOptions {
  /** Separator for text mode */
  sep?: string;
  /** Format string for text values */
  format?: string;
}
```

---

### 9.4 String Formatting

**File:** `src/ts/io/arrayprint.ts` (new file)

```typescript
import { NDArray } from '../NDArray.js';
import { DType } from '../types.js';

/**
 * Print options for array formatting.
 */
export interface PrintOptions {
  /** Number of edge items to show per dimension */
  edgeitems: number;
  /** Total elements before summarization */
  threshold: number;
  /** Float formatting mode */
  floatmode: 'fixed' | 'unique' | 'maxprec' | 'maxprec_equal';
  /** Decimal precision for floats */
  precision: number;
  /** Suppress small values (show as 0) */
  suppress: boolean;
  /** Maximum line width */
  linewidth: number;
  /** String for NaN values */
  nanstr: string;
  /** String for Inf values */
  infstr: string;
  /** Sign display: '-' (default), '+' (always), ' ' (space for positive) */
  sign: '-' | '+' | ' ';
  /** Custom formatter functions by dtype */
  formatter: Record<string, (x: number) => string> | null;
}

/** Default print options */
const defaultOptions: PrintOptions = {
  edgeitems: 3,
  threshold: 1000,
  floatmode: 'maxprec',
  precision: 8,
  suppress: false,
  linewidth: 75,
  nanstr: 'nan',
  infstr: 'inf',
  sign: '-',
  formatter: null,
};

/** Current print options (mutable) */
let currentOptions: PrintOptions = { ...defaultOptions };

/**
 * Set printing options.
 *
 * @param options - Options to set (partial)
 *
 * @example
 * setPrintoptions({ precision: 4, linewidth: 120 });
 */
export function setPrintoptions(options: Partial<PrintOptions>): void {
  currentOptions = { ...currentOptions, ...options };
}

/**
 * Get current print options.
 */
export function getPrintoptions(): PrintOptions {
  return { ...currentOptions };
}

/**
 * Context manager for temporary print options.
 *
 * @example
 * await withPrintoptions({ precision: 2 }, async () => {
 *   console.log(array2string(arr));
 * });
 */
export async function withPrintoptions<T>(
  options: Partial<PrintOptions>,
  fn: () => T | Promise<T>
): Promise<T> {
  const saved = { ...currentOptions };
  try {
    currentOptions = { ...currentOptions, ...options };
    return await fn();
  } finally {
    currentOptions = saved;
  }
}

/**
 * Return a string representation of an array.
 *
 * @param arr - Input array
 * @param options - Formatting options (overrides global settings)
 *
 * @example
 * const str = array2string(arr, { precision: 4, separator: ', ' });
 * console.log(str);
 * // [[1.0000, 2.0000, 3.0000],
 * //  [4.0000, 5.0000, 6.0000]]
 */
export function array2string(
  arr: NDArray,
  options: Partial<Array2StringOptions> = {}
): string {
  const opts = { ...currentOptions, ...options };

  // Handle 0-d array
  if (arr.ndim === 0) {
    return formatScalar(arr.item(), arr.dtype, opts);
  }

  // Determine if summarization is needed
  const summarize = arr.size > opts.threshold;

  // Build string recursively
  return formatArray(arr, 0, opts, summarize);
}

interface Array2StringOptions extends PrintOptions {
  /** Element separator */
  separator?: string;
  /** Prefix for alignment */
  prefix?: string;
  /** Suffix for alignment */
  suffix?: string;
}

/**
 * Format a scalar value.
 */
function formatScalar(
  value: number,
  dtype: DType,
  opts: PrintOptions
): string {
  // Handle special values
  if (Number.isNaN(value)) return opts.nanstr;
  if (!Number.isFinite(value)) {
    return value > 0 ? opts.infstr : `-${opts.infstr}`;
  }

  // Check for custom formatter
  if (opts.formatter) {
    const key = getDtypeFormatterKey(dtype);
    if (opts.formatter[key]) {
      return opts.formatter[key](value);
    }
    if (opts.formatter['all']) {
      return opts.formatter['all'](value);
    }
  }

  // Format based on dtype
  if (isIntegerDtype(dtype)) {
    return formatInteger(value, opts);
  } else if (isFloatDtype(dtype)) {
    return formatFloat(value, opts);
  } else if (isComplexDtype(dtype)) {
    return formatComplex(value, opts);
  }

  return value.toString();
}

/**
 * Format an integer value.
 */
function formatInteger(value: number, opts: PrintOptions): string {
  let str = Math.trunc(value).toString();

  if (opts.sign === '+' && value >= 0) {
    str = '+' + str;
  } else if (opts.sign === ' ' && value >= 0) {
    str = ' ' + str;
  }

  return str;
}

/**
 * Format a floating-point value.
 */
function formatFloat(value: number, opts: PrintOptions): string {
  let str: string;

  switch (opts.floatmode) {
    case 'fixed':
      str = value.toFixed(opts.precision);
      break;
    case 'unique':
      // Minimum digits for uniqueness (using toPrecision with high value)
      str = value.toPrecision(17).replace(/\.?0+$/, '');
      break;
    case 'maxprec':
      str = value.toPrecision(opts.precision);
      // Remove trailing zeros
      if (str.includes('.')) {
        str = str.replace(/\.?0+$/, '');
      }
      break;
    case 'maxprec_equal':
    default:
      str = value.toFixed(opts.precision);
      break;
  }

  // Handle suppress (show very small as 0)
  if (opts.suppress && Math.abs(value) < 1e-10) {
    str = '0';
  }

  // Add sign
  if (opts.sign === '+' && value >= 0 && !str.startsWith('-')) {
    str = '+' + str;
  } else if (opts.sign === ' ' && value >= 0 && !str.startsWith('-')) {
    str = ' ' + str;
  }

  return str;
}

/**
 * Format array recursively.
 */
function formatArray(
  arr: NDArray,
  depth: number,
  opts: PrintOptions & { separator?: string },
  summarize: boolean
): string {
  const sep = opts.separator ?? ' ';
  const indent = ' '.repeat(depth + 1);

  if (arr.ndim === 1) {
    return formatRow(arr, opts, summarize);
  }

  // Multi-dimensional: format each slice
  const parts: string[] = [];
  const n = arr.shape[0];

  if (summarize && n > 2 * opts.edgeitems) {
    // Show edge items with ellipsis
    for (let i = 0; i < opts.edgeitems; i++) {
      const slice = arr.at(i) as NDArray;
      parts.push(indent + formatArray(slice, depth + 1, opts, summarize));
    }
    parts.push(indent + '...');
    for (let i = n - opts.edgeitems; i < n; i++) {
      const slice = arr.at(i) as NDArray;
      parts.push(indent + formatArray(slice, depth + 1, opts, summarize));
    }
  } else {
    for (let i = 0; i < n; i++) {
      const slice = arr.at(i) as NDArray;
      parts.push(indent + formatArray(slice, depth + 1, opts, summarize));
    }
  }

  // Join with appropriate separators
  const innerSep = arr.ndim > 2 ? ',\n\n' : ',\n';
  return '[' + parts.join(innerSep).slice(depth + 1) + ']';
}

/**
 * Format a 1D array row.
 */
function formatRow(
  arr: NDArray,
  opts: PrintOptions & { separator?: string },
  summarize: boolean
): string {
  const sep = opts.separator ?? ' ';
  const values: string[] = [];
  const n = arr.size;

  if (summarize && n > 2 * opts.edgeitems) {
    for (let i = 0; i < opts.edgeitems; i++) {
      values.push(formatScalar(arr.getFlat(i), arr.dtype, opts));
    }
    values.push('...');
    for (let i = n - opts.edgeitems; i < n; i++) {
      values.push(formatScalar(arr.getFlat(i), arr.dtype, opts));
    }
  } else {
    for (let i = 0; i < n; i++) {
      values.push(formatScalar(arr.getFlat(i), arr.dtype, opts));
    }
  }

  return '[' + values.join(',' + sep) + ']';
}

/**
 * Return the official string representation of an array.
 * Includes dtype information.
 *
 * @example
 * console.log(arrayRepr(arr));
 * // array([[1., 2., 3.],
 * //        [4., 5., 6.]], dtype=float64)
 */
export function arrayRepr(
  arr: NDArray,
  options: Partial<PrintOptions> = {}
): string {
  const str = array2string(arr, options);
  const dtypeName = getDtypeName(arr.dtype);

  if (arr.ndim === 0) {
    return `array(${str}, dtype=${dtypeName})`;
  }

  // Indent continuation lines
  const prefix = 'array(';
  const lines = str.split('\n');
  const indentedLines = lines.map((line, i) =>
    i === 0 ? prefix + line : ' '.repeat(prefix.length) + line
  );

  return indentedLines.join('\n') + `, dtype=${dtypeName})`;
}

/**
 * Return a string representation of an array.
 * Without dtype information (like str()).
 */
export function arrayStr(
  arr: NDArray,
  options: Partial<PrintOptions> = {}
): string {
  return array2string(arr, options);
}

/**
 * Format a floating point number in positional notation.
 */
export function formatFloatPositional(
  x: number,
  options: FormatFloatOptions = {}
): string {
  const {
    precision,
    unique = true,
    fractional = true,
    trim = 'k',
    sign = false,
    pad_left,
    pad_right,
  } = options;

  let str: string;

  if (unique) {
    str = x.toPrecision(17);
    // Trim to unique representation
    str = str.replace(/\.?0+$/, '');
  } else if (precision !== undefined) {
    str = fractional
      ? x.toFixed(precision)
      : x.toPrecision(precision);
  } else {
    str = x.toString();
  }

  // Trim options: 'k'=keep, '.'=trim trailing '.', '0'=trim trailing zeros
  if (trim === '.' && str.endsWith('.')) {
    str = str.slice(0, -1);
  } else if (trim === '0') {
    str = str.replace(/\.?0+$/, '');
  }

  // Sign
  if (sign && x >= 0) {
    str = '+' + str;
  }

  // Padding
  if (pad_left !== undefined) {
    const [intPart] = str.split('.');
    if (intPart.length < pad_left) {
      str = ' '.repeat(pad_left - intPart.length) + str;
    }
  }

  return str;
}

interface FormatFloatOptions {
  precision?: number;
  unique?: boolean;
  fractional?: boolean;
  trim?: 'k' | '.' | '0' | '-';
  sign?: boolean;
  pad_left?: number;
  pad_right?: number;
}

/**
 * Format a floating point number in scientific notation.
 */
export function formatFloatScientific(
  x: number,
  options: FormatFloatOptions = {}
): string {
  const {
    precision = 8,
    unique = true,
    trim = 'k',
    sign = false,
    pad_left,
    exp_digits,
  } = options;

  let str = unique
    ? x.toExponential()
    : x.toExponential(precision);

  // Normalize exponent format
  const [mantissa, exp] = str.split('e');
  let expNum = parseInt(exp, 10);
  let expStr = expNum >= 0 ? `+${expNum}` : `${expNum}`;

  if (exp_digits !== undefined) {
    const digits = Math.abs(expNum).toString();
    expStr = (expNum >= 0 ? '+' : '-') + digits.padStart(exp_digits, '0');
  }

  str = mantissa + 'e' + expStr;

  // Sign
  if (sign && x >= 0) {
    str = '+' + str;
  }

  return str;
}

// Helper functions
function getDtypeFormatterKey(dtype: DType): string {
  if (isIntegerDtype(dtype)) return 'int';
  if (isFloatDtype(dtype)) return 'float';
  if (isComplexDtype(dtype)) return 'complexfloat';
  if (dtype === DType.Bool) return 'bool';
  return 'all';
}

function isIntegerDtype(dtype: DType): boolean {
  return [
    DType.Int8, DType.Int16, DType.Int32, DType.Int64,
    DType.UInt8, DType.UInt16, DType.UInt32, DType.UInt64,
  ].includes(dtype);
}

function isFloatDtype(dtype: DType): boolean {
  return [DType.Float16, DType.Float32, DType.Float64].includes(dtype);
}

function isComplexDtype(dtype: DType): boolean {
  return [DType.Complex64, DType.Complex128].includes(dtype);
}

function getDtypeName(dtype: DType): string {
  const names: Record<DType, string> = {
    [DType.Bool]: 'bool',
    [DType.Int8]: 'int8',
    [DType.Int16]: 'int16',
    [DType.Int32]: 'int32',
    [DType.Int64]: 'int64',
    [DType.UInt8]: 'uint8',
    [DType.UInt16]: 'uint16',
    [DType.UInt32]: 'uint32',
    [DType.UInt64]: 'uint64',
    [DType.Float16]: 'float16',
    [DType.Float32]: 'float32',
    [DType.Float64]: 'float64',
    [DType.Complex64]: 'complex64',
    [DType.Complex128]: 'complex128',
  };
  return names[dtype] ?? 'unknown';
}

function formatComplex(value: number, opts: PrintOptions): string {
  // Complex numbers are stored as real, imag pairs
  // This would need special handling - placeholder for now
  return value.toString();
}
```

---

### 9.5 Memory Mapping

**File:** `src/ts/io/memmap.ts` (new file)

```typescript
import { NDArray } from '../NDArray.js';
import { DType, dtypeSize } from '../types.js';
import { parseHeader, createHeader, dtypeToDescr, descrToDtype, NpyHeader } from './format.js';

/**
 * Memory-mapped array stored in a binary file on disk.
 *
 * Memory mapping allows working with arrays larger than available RAM
 * by mapping file contents directly into virtual memory.
 *
 * NOTE: Full memory mapping requires Node.js with mmap support or
 * File System Access API with proper streaming. This implementation
 * provides a compatible API that loads data lazily.
 *
 * @example
 * // Create new memory-mapped file
 * const fp = await Memmap.create('data.npy', DType.Float64, [1000, 1000]);
 * fp.set(0, 0, 42);
 * await fp.flush();
 *
 * // Open existing file
 * const fp = await Memmap.open('data.npy', 'r');
 * console.log(fp.get(0, 0));
 */
export class Memmap extends NDArray {
  private _fileHandle: FileSystemFileHandle | null = null;
  private _filename: string | null = null;
  private _mode: MemmapMode;
  private _offset: number;
  private _dirty: boolean = false;

  private constructor(
    ptr: number,
    module: WasmModule,
    mode: MemmapMode,
    offset: number
  ) {
    super(ptr, module);
    this._mode = mode;
    this._offset = offset;
  }

  /**
   * Create a new memory-mapped file.
   */
  static async create(
    filename: string | FileSystemFileHandle,
    dtype: DType,
    shape: number[],
    options: MemmapCreateOptions = {}
  ): Promise<Memmap> {
    const { order = 'C' } = options;

    // Calculate size
    const size = shape.reduce((a, b) => a * b, 1);
    const dataBytes = size * dtypeSize(dtype);

    // Create header
    const header: NpyHeader = {
      descr: dtypeToDescr(dtype),
      fortran_order: order === 'F',
      shape,
    };
    const headerBytes = createHeader(header, { major: 1, minor: 0 });

    // Create file with header + zeroed data
    const totalSize = headerBytes.length + dataBytes;
    const fileData = new Uint8Array(totalSize);
    fileData.set(headerBytes, 0);
    // Data section is already zeroed

    // Write to file
    if (typeof filename === 'string') {
      if (typeof process !== 'undefined' && process.versions?.node) {
        const fs = await import('fs/promises');
        await fs.writeFile(filename, fileData);
      } else {
        throw new Error('String paths only supported in Node.js');
      }
    } else {
      const writable = await filename.createWritable();
      await writable.write(fileData);
      await writable.close();
    }

    // Open the created file
    return Memmap.open(filename, 'r+');
  }

  /**
   * Open an existing memory-mapped file.
   */
  static async open(
    filename: string | FileSystemFileHandle,
    mode: MemmapMode = 'r'
  ): Promise<Memmap> {
    // Read file
    let data: ArrayBuffer;

    if (typeof filename === 'string') {
      if (typeof process !== 'undefined' && process.versions?.node) {
        const fs = await import('fs/promises');
        const buffer = await fs.readFile(filename);
        data = buffer.buffer.slice(
          buffer.byteOffset,
          buffer.byteOffset + buffer.byteLength
        );
      } else {
        throw new Error('String paths only supported in Node.js');
      }
    } else {
      const file = await filename.getFile();
      data = await file.arrayBuffer();
    }

    // Parse header
    const bytes = new Uint8Array(data);
    const { header, dataOffset } = parseHeader(bytes);

    const dtype = descrToDtype(header.descr);
    const size = header.shape.reduce((a, b) => a * b, 1);

    // Create NDArray from data
    // In a full implementation, this would use actual mmap
    const dataSlice = new Uint8Array(data, dataOffset);
    const arr = NDArray.fromTypedArray(
      createTypedArray(dtype, dataSlice.buffer, dataSlice.byteOffset, size),
      header.shape,
      dtype
    );

    // Create Memmap wrapper
    const memmap = Object.create(Memmap.prototype) as Memmap;
    Object.assign(memmap, arr);
    memmap._mode = mode;
    memmap._offset = dataOffset;
    memmap._filename = typeof filename === 'string' ? filename : null;
    memmap._fileHandle = typeof filename === 'string' ? null : filename;

    // Set writeable flag based on mode
    if (mode === 'r' || mode === 'c') {
      memmap.flags.writeable = false;
    }

    return memmap;
  }

  /**
   * File name.
   */
  get filename(): string | null {
    return this._filename;
  }

  /**
   * File offset where data begins.
   */
  get offset(): number {
    return this._offset;
  }

  /**
   * File mode.
   */
  get mode(): MemmapMode {
    return this._mode;
  }

  /**
   * Write any changes to disk.
   */
  async flush(): Promise<void> {
    if (!this._dirty) return;
    if (this._mode === 'r' || this._mode === 'c') {
      throw new Error('Cannot flush read-only memmap');
    }

    // Get current data
    const dataBytes = new Uint8Array(
      this.data.buffer,
      this.data.byteOffset,
      this.size * dtypeSize(this.dtype)
    );

    // Write to file
    if (this._filename) {
      const fs = await import('fs/promises');
      const handle = await fs.open(this._filename, 'r+');
      await handle.write(dataBytes, 0, dataBytes.length, this._offset);
      await handle.close();
    } else if (this._fileHandle) {
      const writable = await this._fileHandle.createWritable({ keepExistingData: true });
      await writable.seek(this._offset);
      await writable.write(dataBytes);
      await writable.close();
    }

    this._dirty = false;
  }

  /**
   * Override set to track modifications.
   */
  set(...args: number[]): void {
    if (this._mode === 'r') {
      throw new Error('Cannot modify read-only memmap');
    }
    super.set(...args);
    this._dirty = true;
  }
}

/** Memory map file modes */
export type MemmapMode = 'r' | 'r+' | 'w+' | 'c';

interface MemmapCreateOptions {
  /** Memory order: 'C' (row-major) or 'F' (column-major) */
  order?: 'C' | 'F';
}

/**
 * Open a .npy file as a memory-mapped array.
 * Convenience function wrapping Memmap.open().
 */
export async function openMemmap(
  filename: string | FileSystemFileHandle,
  mode: MemmapMode = 'r',
  dtype?: DType,
  shape?: number[],
  options: OpenMemmapOptions = {}
): Promise<Memmap> {
  if (mode === 'w+' && dtype !== undefined && shape !== undefined) {
    return Memmap.create(filename, dtype, shape, options);
  }
  return Memmap.open(filename, mode);
}

interface OpenMemmapOptions extends MemmapCreateOptions {
  /** Fortran (column-major) order */
  fortran_order?: boolean;
}

// Helper function (would be imported from types.ts)
function createTypedArray(dtype: DType, buffer: ArrayBuffer, byteOffset: number, length: number): ArrayBufferView {
  // Implementation would create appropriate typed array based on dtype
  switch (dtype) {
    case DType.Float64: return new Float64Array(buffer, byteOffset, length);
    case DType.Float32: return new Float32Array(buffer, byteOffset, length);
    case DType.Int32: return new Int32Array(buffer, byteOffset, length);
    case DType.Int16: return new Int16Array(buffer, byteOffset, length);
    case DType.Int8: return new Int8Array(buffer, byteOffset, length);
    case DType.UInt32: return new Uint32Array(buffer, byteOffset, length);
    case DType.UInt16: return new Uint16Array(buffer, byteOffset, length);
    case DType.UInt8: return new Uint8Array(buffer, byteOffset, length);
    default: return new Float64Array(buffer, byteOffset, length);
  }
}
```

---

### 9.6 Base Conversion

**File:** `src/ts/io/base.ts` (new file)

```typescript
/**
 * Return the binary representation of the input number as a string.
 *
 * For negative numbers, if width is not given, a minus sign is added.
 * If width is given, returns the two's complement representation.
 *
 * @param num - Integer to convert
 * @param width - Number of bits (enables two's complement for negatives)
 *
 * @example
 * binaryRepr(3)           // '11'
 * binaryRepr(-3)          // '-11'
 * binaryRepr(3, 4)        // '0011'
 * binaryRepr(-3, 4)       // '1101' (two's complement)
 */
export function binaryRepr(num: number, width?: number): string {
  const n = Math.trunc(num);

  if (width === undefined) {
    // Simple case: just convert to binary string
    if (n >= 0) {
      return n.toString(2);
    }
    return '-' + Math.abs(n).toString(2);
  }

  // With width: use two's complement for negatives
  if (n >= 0) {
    return n.toString(2).padStart(width, '0');
  }

  // Two's complement: 2^width + n
  const twosComplement = (1 << width) + n;
  if (twosComplement < 0) {
    throw new Error(
      `Width ${width} is not sufficient to represent ${num} in two's complement`
    );
  }
  return twosComplement.toString(2).padStart(width, '0');
}

/**
 * Return a string representation of a number in the given base.
 *
 * @param number - Integer to convert
 * @param base - Base for conversion (2-36)
 * @param padding - Number of zeros to pad on left
 *
 * @example
 * baseRepr(5)              // '101' (binary, base 2)
 * baseRepr(5, 5)           // '10' (base 5)
 * baseRepr(5, 16)          // '5' (hex)
 * baseRepr(15, 16)         // 'F'
 * baseRepr(5, 2, 4)        // '0101' (4-digit binary)
 */
export function baseRepr(
  number: number,
  base: number = 2,
  padding: number = 0
): string {
  const n = Math.trunc(number);

  if (base < 2 || base > 36) {
    throw new Error('Base must be between 2 and 36');
  }

  if (n === 0) {
    return '0'.repeat(Math.max(1, padding));
  }

  const negative = n < 0;
  let absN = Math.abs(n);
  const digits = '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ';
  let result = '';

  while (absN > 0) {
    result = digits[absN % base] + result;
    absN = Math.floor(absN / base);
  }

  // Apply padding
  if (result.length < padding) {
    result = '0'.repeat(padding - result.length) + result;
  }

  return negative ? '-' + result : result;
}
```

---

## File Structure Summary

### New Files to Create

```
src/ts/io/
├── index.ts           # Re-exports all I/O functions
├── format.ts          # NPY format constants and header parsing
├── npyio.ts           # save, load, savez, savez_compressed, NpzFile
├── text.ts            # loadtxt, savetxt, genfromtxt, fromregex
├── binary.ts          # fromfile, frombuffer, tofile
├── arrayprint.ts      # Print options, array2string, arrayRepr, arrayStr
├── memmap.ts          # Memmap class, openMemmap
└── base.ts            # binaryRepr, baseRepr
```

### Files to Modify

```
src/ts/index.ts
├── Export I/O functions
└── Export print options

src/ts/NDArray.ts
├── Add tofile() method
├── Add toString() using array2string
└── Add [Symbol.toStringTag]

src/ts/types.ts
├── Add I/O-related type definitions
└── Add file helper type utilities

package.json
├── Add optional JSZip dependency for NPZ support
└── Add file system polyfill for browser support (optional)
```

---

## Dependencies

### Required
- None (core TypeScript/JavaScript only)

### Optional
- `jszip`: For NPZ archive support (savez, savez_compressed, NpzFile)
- `pako`: Alternative for compression if avoiding JSZip

### Platform Considerations

| Feature | Node.js | Browser |
|---------|---------|---------|
| save/load to file path | ✅ fs module | ❌ Not supported |
| save/load to ArrayBuffer | ✅ | ✅ |
| save/load to File object | ✅ | ✅ |
| save/load via URL | ✅ fetch | ✅ fetch |
| FileSystemFileHandle | ❌ | ✅ File System Access API |
| Memory mapping | ✅ mmap module | ⚠️ Simulated via ArrayBuffer |
| NPZ archives | ✅ with jszip | ✅ with jszip |

---

## Implementation Order

```
Week 1: Core Binary Format
├── Day 1: format.ts - NPY constants, header creation/parsing
├── Day 2: npyio.ts - save() function
├── Day 3: npyio.ts - load() function
├── Day 4: Unit tests for save/load round-trip
└── Day 5: Edge cases (dtypes, shapes, endianness)

Week 2: NPZ and Text I/O
├── Day 1: npyio.ts - savez(), savez_compressed()
├── Day 2: npyio.ts - NpzFile class
├── Day 3: text.ts - loadtxt()
├── Day 4: text.ts - savetxt()
└── Day 5: text.ts - genfromtxt(), fromregex()

Week 3: Binary I/O and Formatting
├── Day 1: binary.ts - fromfile(), frombuffer()
├── Day 2: NDArray.tofile() method
├── Day 3: arrayprint.ts - PrintOptions, setPrintoptions
├── Day 4: arrayprint.ts - array2string(), arrayRepr()
├── Day 5: arrayprint.ts - formatFloat functions

Week 4: Memory Mapping and Polish
├── Day 1: memmap.ts - Memmap class structure
├── Day 2: memmap.ts - create/open methods
├── Day 3: base.ts - binaryRepr(), baseRepr()
├── Day 4: Integration tests, browser compatibility
└── Day 5: Documentation, examples, cleanup
```

---

## Verification Plan

After Phase 9 completion, verify:

```bash
# Build
npm run build

# Run all tests
npm test

# Specific Phase 9 tests should pass:

# NPY Format
✓ save() creates valid NPY file with correct header
✓ save() handles all supported dtypes
✓ load() reads NPY files created by NumPy
✓ load() reads NPY files created by save()
✓ Round-trip preserves data exactly
✓ Handles C and Fortran order

# NPZ Archives
✓ savez() creates valid ZIP with multiple arrays
✓ savez_compressed() creates compressed archive
✓ NpzFile lazily loads arrays on access
✓ NpzFile supports iteration over all arrays

# Text I/O
✓ loadtxt() parses CSV correctly
✓ loadtxt() handles comments and headers
✓ loadtxt() supports custom converters
✓ savetxt() produces NumPy-compatible output
✓ genfromtxt() handles missing values
✓ fromregex() extracts data correctly

# Binary I/O
✓ fromfile() reads raw binary
✓ fromfile() reads text with separator
✓ frombuffer() creates array from buffer
✓ tofile() writes binary and text

# String Formatting
✓ array2string() produces correct output
✓ Print options affect formatting
✓ Summarization works for large arrays
✓ Custom formatters work
✓ formatFloatPositional/Scientific work

# Memory Mapping
✓ Memmap.create() creates valid file
✓ Memmap.open() reads existing file
✓ Changes persist after flush()
✓ Read-only mode prevents writes

# Base Conversion
✓ binaryRepr() handles positive/negative
✓ binaryRepr() with width uses two's complement
✓ baseRepr() works for bases 2-36
```

Generate NumPy comparison vectors:

```python
# tests/python/generate_phase9_tests.py
import numpy as np
import json
import io

tests = {
    "npy_format": [],
    "text_io": [],
    "formatting": [],
}

# NPY format tests
for dtype in [np.float32, np.float64, np.int32, np.int64, np.bool_]:
    arr = np.arange(12, dtype=dtype).reshape(3, 4)
    buf = io.BytesIO()
    np.save(buf, arr)
    tests["npy_format"].append({
        "dtype": str(dtype),
        "shape": list(arr.shape),
        "data": arr.tolist(),
        "npy_bytes": list(buf.getvalue()),
    })

# Text I/O tests
arr = np.array([[1.5, 2.5, 3.5], [4.5, 5.5, 6.5]])
buf = io.StringIO()
np.savetxt(buf, arr, delimiter=',', fmt='%.2f')
tests["text_io"].append({
    "data": arr.tolist(),
    "csv_output": buf.getvalue(),
})

# Formatting tests
arr = np.array([1.23456789, 0.000001, 1e10, np.nan, np.inf])
tests["formatting"].append({
    "data": arr.tolist(),
    "precision_4": np.array2string(arr, precision=4),
    "precision_8": np.array2string(arr, precision=8),
})

with open("tests/fixtures/phase9_vectors.json", "w") as f:
    json.dump(tests, f, indent=2)
```

---

## NumPy Reference Files

| Component | NumPy Reference File |
|-----------|---------------------|
| NPY format | `numpy/lib/_format_impl.py` |
| save/load/savez | `numpy/lib/_npyio_impl.py` |
| loadtxt/savetxt | `numpy/lib/_npyio_impl.py` |
| genfromtxt | `numpy/lib/_npyio_impl.py` |
| NpzFile | `numpy/lib/_npyio_impl.py` |
| memmap | `numpy/_core/memmap.py` |
| array2string | `numpy/_core/arrayprint.py` |
| printoptions | `numpy/_core/printoptions.py` |
| binary_repr | `numpy/_core/numeric.py` |
| base_repr | `numpy/_core/numeric.py` |

---

## Critical Notes

### Security Considerations
- **MAX_HEADER_SIZE**: Limit header parsing to prevent attacks via malformed files
- **allow_pickle**: Disabled by default in load() to prevent arbitrary code execution
- **Header validation**: Strict parsing of NPY headers to prevent injection

### Browser Compatibility
- File paths only work in Node.js
- Use FileSystemFileHandle or ArrayBuffer in browsers
- Memory mapping is simulated (no true mmap in browsers)

### Performance Considerations
- Large arrays should use chunked writing (16 MiB buffers)
- NPZ lazy-loading prevents memory issues with large archives
- Consider streaming for very large text files

---

## Phase 9 Enables

After Phase 9 completion:
- Arrays can be persisted and shared between sessions
- Data can be imported from CSV, text files, and raw binary
- Arrays can be displayed with customizable formatting
- Memory-efficient handling of large datasets via memory mapping
- Interoperability with NumPy's native file formats
