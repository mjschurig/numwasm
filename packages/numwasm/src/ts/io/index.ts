/**
 * I/O Operations Module
 *
 * Re-exports all I/O functions for array persistence, text parsing,
 * and string formatting.
 */

// NPY format constants and utilities
export {
  MAGIC_PREFIX,
  MAGIC_LEN,
  ARRAY_ALIGN,
  MAX_HEADER_SIZE,
  type NpyVersion,
  type NpyHeader,
  getMinVersion,
  dtypeToDescr,
  descrToDtype,
  dtypeSize,
  createHeader,
  parseHeader,
  isNode,
} from './format.js';

// NPY file I/O
export {
  save,
  load,
  type SaveOptions,
  type LoadOptions,
  // NPZ file I/O
  savez,
  savez_compressed,
  loadz,
  type SavezOptions,
  type NpzFile,
  // Bit packing
  packbits,
  unpackbits,
  type BitOrder,
} from './npyio.js';

// Text file I/O
export {
  loadtxt,
  savetxt,
  genfromtxt,
  fromregex,
  formatValue,
  type LoadtxtOptions,
  type SavetxtOptions,
  type GenfromtxtOptions,
  type FromregexOptions,
} from './text.js';

// Binary file I/O
export {
  fromfile,
  frombuffer,
  type FromfileOptions,
  type FrombufferOptions,
} from './binary.js';

// Array printing and formatting
export {
  setPrintoptions,
  getPrintoptions,
  resetPrintoptions,
  withPrintoptions,
  array2string,
  arrayRepr,
  arrayStr,
  formatFloatPositional,
  formatFloatScientific,
  type PrintOptions,
} from './arrayprint.js';

// Memory-mapped files
export {
  Memmap,
  openMemmap,
  type MemmapMode,
  type MemmapOptions,
} from './memmap.js';

// Base conversion
export {
  binaryRepr,
  baseRepr,
  fromBinaryRepr,
  fromBaseRepr,
  hexRepr,
  octalRepr,
  binaryReprArray,
  baseReprArray,
} from './base.js';
