/**
 * Text File I/O
 *
 * Functions for reading and writing arrays from/to text files.
 *
 * Reference: numpy/lib/_npyio_impl.py
 */

import { NDArray } from '../NDArray.js';
import { DType } from '../types.js';
import { isNode } from './format.js';

/**
 * Options for loadtxt()
 */
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

/**
 * Options for savetxt()
 */
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
 * Options for genfromtxt()
 */
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

/**
 * Options for fromregex()
 */
export interface FromregexOptions {
  /** File encoding */
  encoding?: string;
}

/**
 * Load data from a text file.
 *
 * Each row in the text file must have the same number of values.
 *
 * @param file - File source (path, File object, or URL)
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
    delimiter,
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

    if (line.length === 0) continue;

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
  let arr: NDArray;
  if (data[0].length === 1) {
    // Single column - create 1D array
    arr = await NDArray.fromArray(data.map(row => row[0]), dtype);
  } else {
    arr = await NDArray.fromArray(data, dtype);
  }

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

/**
 * Save an array to a text file.
 *
 * @param file - Output file (path or FileSystemFileHandle)
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

/**
 * Load data from a text file with handling for missing values.
 * More sophisticated than loadtxt, supports missing value replacement.
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
      ? usecols.map(idx => parts[idx])
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

  if (data.length === 0) {
    throw new Error('No data found in file');
  }

  return NDArray.fromArray(data, dtype);
}

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

/**
 * Format a number according to printf-style format string.
 */
export function formatValue(value: number, fmt: string): string {
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

  // Apply sign
  if (flags.includes('+') && value >= 0 && !result.startsWith('-')) {
    result = '+' + result;
  } else if (flags.includes(' ') && value >= 0 && !result.startsWith('-')) {
    result = ' ' + result;
  }

  // Apply width padding
  if (width) {
    const w = parseInt(width, 10);
    const pad = flags.includes('-') ? 'end' : 'start';
    const char = flags.includes('0') && !flags.includes('-') ? '0' : ' ';
    result = pad === 'start'
      ? result.padStart(w, char)
      : result.padEnd(w, char);
  }

  return result;
}

/**
 * Read text from various file sources.
 */
async function readTextFile(
  file: string | File | URL,
  encoding: string = 'utf-8'
): Promise<string> {
  // Check for File in browser
  if (typeof File !== 'undefined' && file instanceof File) {
    return file.text();
  }

  if (file instanceof URL || (typeof file === 'string' && file.startsWith('http'))) {
    const response = await fetch(file.toString());
    return response.text();
  }

  if (typeof file === 'string') {
    // Node.js file path
    if (isNode()) {
      const fs = await import('fs/promises');
      return fs.readFile(file, { encoding: encoding as BufferEncoding });
    }
    throw new Error('String file paths only supported in Node.js');
  }

  throw new Error('Unsupported file source');
}

/**
 * Write text to file.
 */
async function writeTextFile(
  file: string | FileSystemFileHandle,
  text: string,
  encoding: string = 'utf-8'
): Promise<void> {
  if (typeof file === 'string') {
    // Node.js file path
    if (isNode()) {
      const fs = await import('fs/promises');
      await fs.writeFile(file, text, { encoding: encoding as BufferEncoding });
      return;
    }
    throw new Error('String file paths only supported in Node.js');
  }

  // FileSystemFileHandle (File System Access API)
  if ('createWritable' in file) {
    const writable = await (file as FileSystemFileHandle).createWritable();
    await writable.write(text);
    await writable.close();
    return;
  }

  throw new Error('Unsupported file target');
}

// Type declaration for FileSystemFileHandle
interface FileSystemFileHandle {
  createWritable(): Promise<FileSystemWritableFileStream>;
}

interface FileSystemWritableFileStream {
  write(data: string | ArrayBuffer | Uint8Array): Promise<void>;
  close(): Promise<void>;
}
