/**
 * Array Printing and String Formatting
 *
 * Functions for converting arrays to human-readable string representations.
 *
 * Reference: numpy/_core/arrayprint.py
 */

import { NDArray } from '../NDArray.js';
import { DType } from '../types.js';

/**
 * Print formatting options.
 */
export interface PrintOptions {
  /** Total digits of precision for floating point (default: 8) */
  precision?: number;
  /** Maximum total array elements to print (default: 1000) */
  threshold?: number;
  /** Number of items at beginning/end when summarizing (default: 3) */
  edgeitems?: number;
  /** Characters per line for insertion of line breaks (default: 75) */
  linewidth?: number;
  /** Whether to suppress small floating point values (default: false) */
  suppress?: boolean;
  /** String inserted between elements (default: ' ') */
  separator?: string;
  /** Prefix for array output */
  prefix?: string;
  /** Suffix for array output */
  suffix?: string;
  /** Custom formatter functions per dtype */
  formatter?: Partial<Record<DType, (x: number) => string>>;
  /** Float formatting mode */
  floatmode?: 'fixed' | 'unique' | 'maxprec' | 'maxprec_equal';
  /** Show sign for positive numbers */
  sign?: '-' | '+' | ' ';
  /** Representation of infinity */
  infstr?: string;
  /** Representation of NaN */
  nanstr?: string;
  /** Legacy mode for compatibility */
  legacy?: string | false;
}

/**
 * Default print options matching NumPy defaults.
 */
const defaultOptions: Required<PrintOptions> = {
  precision: 8,
  threshold: 1000,
  edgeitems: 3,
  linewidth: 75,
  suppress: false,
  separator: ' ',
  prefix: '',
  suffix: '',
  formatter: {},
  floatmode: 'maxprec',
  sign: '-',
  infstr: 'inf',
  nanstr: 'nan',
  legacy: false,
};

/**
 * Current global print options.
 */
let currentOptions: Required<PrintOptions> = { ...defaultOptions };

/**
 * Set printing options globally.
 *
 * @param options - New print options (partial, merged with current)
 *
 * @example
 * setPrintoptions({ precision: 4 });
 * setPrintoptions({ threshold: 100, edgeitems: 5 });
 */
export function setPrintoptions(options: PrintOptions): void {
  currentOptions = { ...currentOptions, ...options };
}

/**
 * Get current print options.
 *
 * @returns Current print options
 */
export function getPrintoptions(): Required<PrintOptions> {
  return { ...currentOptions };
}

/**
 * Reset print options to defaults.
 */
export function resetPrintoptions(): void {
  currentOptions = { ...defaultOptions };
}

/**
 * Execute a function with temporary print options.
 *
 * @param options - Temporary print options
 * @param fn - Function to execute
 * @returns Result of function
 *
 * @example
 * withPrintoptions({ precision: 2 }, () => {
 *   console.log(array2string(arr));
 * });
 */
export function withPrintoptions<T>(options: PrintOptions, fn: () => T): T {
  const saved = currentOptions;
  currentOptions = { ...currentOptions, ...options };
  try {
    return fn();
  } finally {
    currentOptions = saved;
  }
}

/**
 * Return a string representation of an array.
 *
 * @param arr - Array to convert
 * @param options - Override print options
 * @returns String representation
 *
 * @example
 * const arr = NDArray.fromArray([[1, 2, 3], [4, 5, 6]], DType.Float64);
 * console.log(array2string(arr));
 * // [[1. 2. 3.]
 * //  [4. 5. 6.]]
 */
export function array2string(
  arr: NDArray,
  options: PrintOptions = {}
): string {
  const opts = { ...currentOptions, ...options };

  if (arr.ndim === 0) {
    // Scalar
    const value = arr.get();
    return formatScalar(value, arr.dtype, opts);
  }

  const shouldSummarize = arr.size > opts.threshold;

  return formatArray(arr, opts, shouldSummarize, 0);
}

/**
 * Return the string representation of an array (with dtype info).
 *
 * @param arr - Array to represent
 * @param options - Override print options
 * @returns Full representation with dtype
 *
 * @example
 * console.log(arrayRepr(arr));
 * // array([[1., 2., 3.],
 * //        [4., 5., 6.]], dtype=float64)
 */
export function arrayRepr(
  arr: NDArray,
  options: PrintOptions = {}
): string {
  const opts = { ...currentOptions, ...options };
  const prefix = 'array(';
  const suffix = ')';

  if (arr.ndim === 0) {
    // Scalar
    const value = arr.get();
    return `${prefix}${formatScalar(value, arr.dtype, opts)}, dtype=${dtypeName(arr.dtype)}${suffix}`;
  }

  const shouldSummarize = arr.size > opts.threshold;
  const body = formatArray(arr, opts, shouldSummarize, prefix.length);

  // Determine if we need dtype suffix
  const dtypeSuffix = needsDtypeSuffix(arr.dtype)
    ? `, dtype=${dtypeName(arr.dtype)}`
    : '';

  return `${prefix}${body}${dtypeSuffix}${suffix}`;
}

/**
 * Return a string representation of an array (no dtype).
 *
 * @param arr - Array to represent
 * @param options - Override print options
 * @returns String representation without dtype info
 */
export function arrayStr(
  arr: NDArray,
  options: PrintOptions = {}
): string {
  return array2string(arr, options);
}

/**
 * Format a single scalar value.
 */
function formatScalar(
  value: number | bigint,
  dtype: DType,
  opts: Required<PrintOptions>
): string {
  // Check for custom formatter
  if (opts.formatter[dtype]) {
    return opts.formatter[dtype]!(Number(value));
  }

  const numValue = typeof value === 'bigint' ? Number(value) : value;

  // Handle special values
  if (Number.isNaN(numValue)) {
    return opts.nanstr;
  }
  if (!Number.isFinite(numValue)) {
    return numValue > 0 ? opts.infstr : `-${opts.infstr}`;
  }

  // Integer types
  if (isIntegerDtype(dtype)) {
    return formatInteger(numValue, opts);
  }

  // Float types
  return formatFloat(numValue, opts);
}

/**
 * Format an integer value.
 */
function formatInteger(value: number, opts: Required<PrintOptions>): string {
  const intVal = Math.trunc(value);
  let result = intVal.toString();

  if (opts.sign === '+' && intVal >= 0) {
    result = '+' + result;
  } else if (opts.sign === ' ' && intVal >= 0) {
    result = ' ' + result;
  }

  return result;
}

/**
 * Format a floating point value.
 */
function formatFloat(value: number, opts: Required<PrintOptions>): string {
  let result: string;

  if (opts.suppress && Math.abs(value) < Math.pow(10, -opts.precision)) {
    result = '0.';
  } else if (opts.floatmode === 'fixed') {
    result = value.toFixed(opts.precision);
  } else {
    // Default: use reasonable precision
    result = formatFloatSmart(value, opts.precision);
  }

  // Add sign
  if (opts.sign === '+' && value >= 0 && !result.startsWith('-')) {
    result = '+' + result;
  } else if (opts.sign === ' ' && value >= 0 && !result.startsWith('-')) {
    result = ' ' + result;
  }

  return result;
}

/**
 * Smart float formatting that avoids unnecessary trailing zeros.
 */
function formatFloatSmart(value: number, precision: number): string {
  // Check if value is effectively an integer
  if (Number.isInteger(value) && Math.abs(value) < 1e15) {
    return value.toString() + '.';
  }

  // Use exponential for very large or very small numbers
  const absValue = Math.abs(value);
  if (absValue !== 0 && (absValue >= 1e16 || absValue < 1e-4)) {
    return value.toExponential(precision - 1);
  }

  // Standard decimal notation
  const fixed = value.toFixed(precision);

  // Remove trailing zeros but keep at least one decimal digit
  const parts = fixed.split('.');
  if (parts.length === 2) {
    let decimal = parts[1].replace(/0+$/, '');
    if (decimal.length === 0) {
      decimal = '';
    }
    return parts[0] + '.' + decimal;
  }

  return fixed;
}

/**
 * Format a float in positional notation.
 */
export function formatFloatPositional(
  value: number,
  precision: number = 8,
  unique: boolean = true,
  fractional: boolean = true,
  trim: 'k' | '.' | '0' | '-' = 'k',
  sign: boolean = false,
  padLeft: number = 0,
  padRight: number = 0
): string {
  let result: string;

  if (fractional) {
    result = value.toFixed(precision);
  } else {
    result = value.toPrecision(precision);
  }

  // Trim trailing zeros
  if (trim !== 'k') {
    if (trim === '0' || trim === '.') {
      result = result.replace(/(\.\d*?)0+$/, '$1');
    }
    if (trim === '.') {
      result = result.replace(/\.$/, '');
    }
    if (trim === '-') {
      result = result.replace(/\.?0+$/, '');
    }
  }

  // Add sign
  if (sign && value >= 0 && !result.startsWith('-')) {
    result = '+' + result;
  }

  // Padding
  if (padLeft > 0) {
    result = result.padStart(padLeft, ' ');
  }
  if (padRight > 0) {
    result = result.padEnd(padRight, ' ');
  }

  return result;
}

/**
 * Format a float in scientific notation.
 */
export function formatFloatScientific(
  value: number,
  precision: number = 8,
  unique: boolean = true,
  trim: 'k' | '.' | '0' | '-' = 'k',
  sign: boolean = false,
  padLeft: number = 0,
  expDigits: number = 2
): string {
  let result = value.toExponential(precision);

  // Normalize exponent digits
  const match = result.match(/e([+-])(\d+)$/);
  if (match) {
    const expSign = match[1];
    const exp = match[2].padStart(expDigits, '0');
    result = result.replace(/e[+-]\d+$/, `e${expSign}${exp}`);
  }

  // Trim trailing zeros in mantissa
  if (trim !== 'k') {
    const [mantissa, exp] = result.split('e');
    let trimmed = mantissa;
    if (trim === '0' || trim === '.') {
      trimmed = mantissa.replace(/(\.\d*?)0+$/, '$1');
    }
    if (trim === '.') {
      trimmed = trimmed.replace(/\.$/, '');
    }
    result = trimmed + 'e' + exp;
  }

  // Add sign
  if (sign && value >= 0 && !result.startsWith('-')) {
    result = '+' + result;
  }

  // Padding
  if (padLeft > 0) {
    result = result.padStart(padLeft, ' ');
  }

  return result;
}

/**
 * A lightweight view into a subarray for formatting purposes.
 * Implements the minimal interface needed for formatArray/formatRow.
 */
class SubarrayView {
  private _parent: NDArray | SubarrayView;
  private _shape: number[];
  private _offset: number;
  private _dtype: DType;

  constructor(parent: NDArray | SubarrayView, shape: number[], offset: number, dtype: DType) {
    this._parent = parent;
    this._shape = shape;
    this._offset = offset;
    this._dtype = dtype;
  }

  get ndim(): number {
    return this._shape.length;
  }

  get shape(): number[] {
    return this._shape;
  }

  get size(): number {
    return this._shape.length === 0 ? 1 : this._shape.reduce((a, b) => a * b, 1);
  }

  get dtype(): DType {
    return this._dtype;
  }

  get(i: number): number {
    // Get the root NDArray and calculate absolute offset
    let root: NDArray | SubarrayView = this._parent;
    let absoluteOffset = this._offset + i;

    while (root instanceof SubarrayView) {
      absoluteOffset += root._offset;
      root = root._parent;
    }

    return (root as NDArray).getFlat(absoluteOffset);
  }

  getFlat(i: number): number {
    return this.get(i);
  }
}

/**
 * Get a subarray at a given index along the first axis.
 * This creates a simulated view for formatting purposes.
 */
function getSubarray(arr: NDArray | SubarrayView, rowIndex: number): SubarrayView {
  // Get the shape without the first dimension
  const parentShape = arr.shape;
  const newShape = parentShape.slice(1);

  // Calculate the stride for the first dimension
  const firstStride = newShape.length > 0 ? newShape.reduce((a, b) => a * b, 1) : 1;
  const startOffset = rowIndex * firstStride;

  return new SubarrayView(arr, newShape, startOffset, arr.dtype);
}

/**
 * Array-like type for formatting (supports both NDArray and SubarrayView).
 */
type ArrayLike = NDArray | SubarrayView;

/**
 * Format an n-dimensional array recursively.
 */
function formatArray(
  arr: ArrayLike,
  opts: Required<PrintOptions>,
  summarize: boolean,
  indent: number
): string {
  if (arr.ndim === 1) {
    return formatRow(arr, opts, summarize);
  }

  // Multi-dimensional: format as nested arrays
  const lines: string[] = [];
  const n = arr.shape[0];
  const edge = opts.edgeitems;

  lines.push('[');

  for (let i = 0; i < n; i++) {
    // Summarize middle elements
    if (summarize && n > 2 * edge && i === edge) {
      const spaces = ' '.repeat(indent + 1);
      lines.push(`${spaces}...`);
      i = n - edge - 1;
      continue;
    }

    // Get subarray at index i along first axis
    const subArr = getSubarray(arr, i);
    const formatted = formatArray(subArr, opts, summarize, indent + 1);
    const prefix = i === 0 ? '' : ' '.repeat(indent + 1);
    const suffix = i < n - 1 ? ',' : '';

    if (arr.ndim === 2) {
      lines.push(`${prefix}${formatted}${suffix}`);
    } else {
      // For 3D+, add extra newlines between sub-arrays
      if (i > 0 && (i < edge || i >= n - edge)) {
        lines.push('');
      }
      lines.push(`${prefix}${formatted}${suffix}`);
    }
  }

  lines.push(']');

  return lines.join('\n');
}

/**
 * Format a 1D array row.
 */
function formatRow(
  arr: ArrayLike,
  opts: Required<PrintOptions>,
  summarize: boolean
): string {
  const n = arr.shape[0];
  const edge = opts.edgeitems;
  const items: string[] = [];

  for (let i = 0; i < n; i++) {
    // Summarize middle elements
    if (summarize && n > 2 * edge && i === edge) {
      items.push('...');
      i = n - edge - 1;
      continue;
    }

    const value = arr.get(i);
    items.push(formatScalar(value, arr.dtype, opts));
  }

  return '[' + items.join(opts.separator) + ']';
}

/**
 * Check if dtype is an integer type.
 */
function isIntegerDtype(dtype: DType): boolean {
  return [
    DType.Int8, DType.Int16, DType.Int32, DType.Int64,
    DType.Uint8, DType.Uint16, DType.Uint32, DType.Uint64,
    DType.Bool,
  ].includes(dtype);
}

/**
 * Check if dtype needs explicit suffix in repr.
 */
function needsDtypeSuffix(dtype: DType): boolean {
  // Float64 is the default, so no suffix needed
  return dtype !== DType.Float64;
}

/**
 * Get human-readable dtype name.
 */
function dtypeName(dtype: DType): string {
  const names: Record<DType, string> = {
    [DType.Bool]: 'bool',
    [DType.Int8]: 'int8',
    [DType.Int16]: 'int16',
    [DType.Int32]: 'int32',
    [DType.Int64]: 'int64',
    [DType.Uint8]: 'uint8',
    [DType.Uint16]: 'uint16',
    [DType.Uint32]: 'uint32',
    [DType.Uint64]: 'uint64',
    [DType.Float16]: 'float16',
    [DType.Float32]: 'float32',
    [DType.Float64]: 'float64',
    [DType.Complex64]: 'complex64',
    [DType.Complex128]: 'complex128',
    [DType.String]: 'string',
  };
  return names[dtype] ?? 'unknown';
}
