/**
 * NPY File Format Implementation
 *
 * Implements NumPy's .npy binary file format for array serialization.
 *
 * Format Structure:
 * - Magic string: \x93NUMPY (6 bytes)
 * - Version: major (1 byte) + minor (1 byte)
 * - Header length: 2 bytes (v1.0) or 4 bytes (v2.0+)
 * - Header: ASCII dict literal with 'descr', 'fortran_order', 'shape'
 * - Padding: spaces/newline to align to 64 bytes
 * - Data: raw array bytes
 *
 * Reference: numpy/lib/_format_impl.py
 */

import { DType, DTYPE_SIZES } from '../types.js';

/** Magic prefix for NPY files: \x93NUMPY */
export const MAGIC_PREFIX = new Uint8Array([0x93, 0x4e, 0x55, 0x4d, 0x50, 0x59]);

/** Total magic length including version bytes */
export const MAGIC_LEN = 8;

/** Alignment boundary for data section (enables memory mapping) */
export const ARRAY_ALIGN = 64;

/** Maximum header size for security (prevents parsing attacks) */
export const MAX_HEADER_SIZE = 10000;

/** NPY format version */
export interface NpyVersion {
  major: number;
  minor: number;
}

/** NPY header dictionary */
export interface NpyHeader {
  /** dtype descriptor (e.g., '<f8', '|b1') */
  descr: string;
  /** True if Fortran (column-major) order */
  fortran_order: boolean;
  /** Array dimensions */
  shape: number[];
}

/**
 * Determine minimum NPY version required for a dtype.
 * - v1.0: Basic types, header < 65536 bytes
 * - v2.0: Large headers (> 65535 bytes)
 * - v3.0: Unicode field names in structured dtypes
 */
export function getMinVersion(headerSize: number): NpyVersion {
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
 * - Type: 'b' (bool), 'i' (int), 'u' (uint), 'f' (float), 'c' (complex)
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
    [DType.Uint8]: '|u1',
    [DType.Uint16]: '<u2',
    [DType.Uint32]: '<u4',
    [DType.Uint64]: '<u8',
    [DType.Float16]: '<f2',
    [DType.Float32]: '<f4',
    [DType.Float64]: '<f8',
    [DType.Complex64]: '<c8',
    [DType.Complex128]: '<c16',
    [DType.String]: '|U0', // Variable-length Unicode string
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
    'u1': DType.Uint8,
    'u2': DType.Uint16,
    'u4': DType.Uint32,
    'u8': DType.Uint64,
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
 * Get the byte size for a dtype.
 */
export function dtypeSize(dtype: DType): number {
  return DTYPE_SIZES[dtype];
}

/**
 * Create NPY header bytes.
 */
export function createHeader(header: NpyHeader, version: NpyVersion): Uint8Array {
  // Format header as Python dict literal
  const shapeStr = header.shape.length === 0
    ? '()'
    : header.shape.length === 1
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
  // data[7] is minor version, currently unused

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

/**
 * Check if running in Node.js environment.
 */
export function isNode(): boolean {
  return typeof process !== 'undefined' && process.versions?.node !== undefined;
}
