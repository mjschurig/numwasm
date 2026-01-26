/**
 * NumJS Record Arrays - Format Parser
 *
 * Parses format specifications into structured dtypes.
 * Compatible with NumPy's format string syntax.
 */

import {
  DType,
  type StructuredDType,
  type FieldDescriptor,
  dtypeSize,
  dtypeAlignment,
} from '../types.js';

/* ============ Error Classes ============ */

/**
 * Error for invalid values or arguments.
 */
export class ValueError extends Error {
  constructor(message: string) {
    super(message);
    this.name = 'ValueError';
  }
}

/* ============ Internal Types ============ */

interface ParsedFormat {
  dtype: DType;
  itemsize: number;
  byteorder?: string;
  charType?: 'S' | 'U';
}

/* ============ Format Parser Class ============ */

/**
 * Class to convert format strings, names, and titles to a structured dtype.
 *
 * Parses format specifications like 'f8, i4, S10' into proper field descriptors.
 *
 * @example
 * const parser = new format_parser('f8, i4, S10', 'x, y, name');
 * console.log(parser.dtype);
 * // { names: ['x', 'y', 'name'], fields: Map {...}, itemsize: 22 }
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
   * @param byteorder - Byte order: '<' little, '>' big, '=' native, '|' not applicable
   */
  constructor(
    formats: string | string[],
    names: string | string[] | null = null,
    titles: string | string[] | null = null,
    aligned: boolean = false,
    byteorder: '<' | '>' | '=' | '|' | string | null = null
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
   * Parse format specification into list of ParsedFormat objects.
   */
  private _parseFormats(formats: string | string[]): ParsedFormat[] {
    let formatList: string[];

    if (typeof formats === 'string') {
      // Handle comma-separated string: 'f8, i4, S10'
      formatList = formats.split(',').map(f => f.trim()).filter(f => f.length > 0);
    } else {
      formatList = [...formats];
    }

    if (formatList.length === 0) {
      throw new ValueError('At least one format must be specified');
    }

    // Parse each format string
    return formatList.map(fmt => this._parseSingleFormat(fmt));
  }

  /**
   * Parse a single format string into dtype and size.
   *
   * Valid formats:
   * - Type codes: b, i, u, f, c, S, U, V
   * - With size: i4, f8, S10, U20
   * - With byte order: <f8, >i4, =f4
   */
  private _parseSingleFormat(fmt: string): ParsedFormat {
    // Pattern: optional byteorder + type char + optional size
    const match = fmt.match(/^([<>=|])?([biufcmMOSUV?])(\d*)$/);
    if (!match) {
      throw new TypeError(`Invalid format specification: '${fmt}'`);
    }

    const [, byteorder, typeChar, sizeStr] = match;
    const size = sizeStr ? parseInt(sizeStr, 10) : null;

    // Convert type character to DType
    const { dtype, charType } = this._charToDtype(typeChar, size);

    // Calculate itemsize
    let itemsize: number;
    if (typeChar === 'S') {
      // ASCII string: size is number of characters (1 byte each)
      itemsize = size ?? 1;
    } else if (typeChar === 'U') {
      // Unicode string: size is number of characters (4 bytes each)
      itemsize = (size ?? 1) * 4;
    } else if (typeChar === 'V') {
      // Void type: size is number of bytes
      itemsize = size ?? 1;
    } else {
      // Numeric types: use dtype size
      itemsize = dtypeSize(dtype);
    }

    return { dtype, itemsize, byteorder, charType };
  }

  /**
   * Convert type character to DType.
   */
  private _charToDtype(
    typeChar: string,
    size: number | null
  ): { dtype: DType; charType?: 'S' | 'U' } {
    switch (typeChar) {
      case 'b':
        return { dtype: DType.Int8 };
      case '?':
        return { dtype: DType.Bool };
      case 'i':
        switch (size) {
          case 1: return { dtype: DType.Int8 };
          case 2: return { dtype: DType.Int16 };
          case 4: return { dtype: DType.Int32 };
          case 8: return { dtype: DType.Int64 };
          default: return { dtype: DType.Int32 }; // Default
        }
      case 'u':
        switch (size) {
          case 1: return { dtype: DType.Uint8 };
          case 2: return { dtype: DType.Uint16 };
          case 4: return { dtype: DType.Uint32 };
          case 8: return { dtype: DType.Uint64 };
          default: return { dtype: DType.Uint32 }; // Default
        }
      case 'f':
        switch (size) {
          case 2: return { dtype: DType.Float16 };
          case 4: return { dtype: DType.Float32 };
          case 8: return { dtype: DType.Float64 };
          default: return { dtype: DType.Float64 }; // Default
        }
      case 'c':
        switch (size) {
          case 8: return { dtype: DType.Complex64 };
          case 16: return { dtype: DType.Complex128 };
          default: return { dtype: DType.Complex128 }; // Default
        }
      case 'S':
        return { dtype: DType.String, charType: 'S' };
      case 'U':
        return { dtype: DType.String, charType: 'U' };
      case 'V':
        return { dtype: DType.Uint8 }; // Void treated as bytes
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
    const duplicates = find_duplicate(names);

    if (duplicates.length > 0) {
      throw new ValueError(`Duplicate field names: ${duplicates.join(', ')}`);
    }
  }

  /**
   * Create the structured dtype from parsed components.
   */
  private _createDtype(
    formats: ParsedFormat[],
    names: string[],
    titles: (string | null)[],
    aligned: boolean,
    _byteorder: string | null
  ): StructuredDType {
    const fieldList: FieldDescriptor[] = [];
    const fields = new Map<string, FieldDescriptor>();
    let offset = 0;
    let maxAlignment = 1;

    for (let i = 0; i < formats.length; i++) {
      const { dtype, itemsize, charType } = formats[i];

      // Handle alignment padding
      if (aligned && offset > 0) {
        const alignment = dtypeAlignment(dtype);
        maxAlignment = Math.max(maxAlignment, alignment);
        const padding = (alignment - (offset % alignment)) % alignment;
        offset += padding;
      }

      const field: FieldDescriptor = {
        name: names[i],
        dtype: dtype,
        offset: offset,
        title: titles[i],
        itemsize: itemsize,
      };

      if (charType) {
        field.charType = charType;
      }

      fieldList.push(field);
      fields.set(names[i], field);

      // Add title as alias if provided
      if (titles[i] && titles[i] !== names[i]) {
        fields.set(titles[i]!, field);
      }

      offset += itemsize;
    }

    // Add final padding for struct alignment
    if (aligned && maxAlignment > 1) {
      const padding = (maxAlignment - (offset % maxAlignment)) % maxAlignment;
      offset += padding;
    }

    return {
      names: names,
      fields: fields,
      fieldList: fieldList,
      itemsize: offset,
      alignment: aligned ? maxAlignment : 1,
      isAligned: aligned,
    };
  }
}

/* ============ Utility Functions ============ */

/**
 * Find duplicate elements in a list.
 *
 * @param list - Array of strings to check
 * @returns Array of duplicate strings
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

/**
 * Convert DType to format string.
 */
export function dtypeToFormat(dtype: DType): string {
  const formats: Record<DType, string> = {
    [DType.Bool]: '?',
    [DType.Int8]: 'i1',
    [DType.Int16]: 'i2',
    [DType.Int32]: 'i4',
    [DType.Int64]: 'i8',
    [DType.Uint8]: 'u1',
    [DType.Uint16]: 'u2',
    [DType.Uint32]: 'u4',
    [DType.Uint64]: 'u8',
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
 * Infer format string from a JavaScript array.
 */
export function inferFormat(arr: unknown[]): string {
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
    return '?';
  }

  // Default
  return 'f8';
}
