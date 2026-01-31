/**
 * Memory-Mapped File Arrays
 *
 * Provides memory-mapped array access for large files.
 * In browsers, this is simulated using ArrayBuffer (no true mmap).
 * In Node.js, it can use fs.promises for file access.
 *
 * Reference: numpy/core/memmap.py
 */

import { NDArray } from "../_core/NDArray.js";
import { DType, DTYPE_SIZES } from "../types.js";
import { isNode } from "./format.js";

/**
 * Memory mapping modes.
 */
export type MemmapMode = "r" | "r+" | "w+" | "c";

/**
 * Options for creating/opening a memmap.
 */
export interface MemmapOptions {
  /** Data type */
  dtype?: DType;
  /** File access mode */
  mode?: MemmapMode;
  /** Byte offset in file where array data starts */
  offset?: number;
  /** Array shape (required for new files) */
  shape?: number[];
  /** Storage order: 'C' (row-major) or 'F' (column-major) */
  order?: "C" | "F";
}

/**
 * Memory-mapped array that reads/writes directly to a file.
 *
 * This class extends NDArray to provide file-backed storage.
 * Changes to the array can be flushed to disk.
 *
 * @example
 * // Create a new memory-mapped file
 * const mmap = await Memmap.create('large_data.dat', {
 *   dtype: DType.Float64,
 *   shape: [1000, 1000],
 *   mode: 'w+',
 * });
 *
 * // Modify data
 * mmap.set(0, 0, 3.14);
 *
 * // Flush changes to disk
 * await mmap.flush();
 *
 * // Open existing memory-mapped file
 * const existing = await Memmap.open('large_data.dat', {
 *   dtype: DType.Float64,
 *   shape: [1000, 1000],
 *   mode: 'r+',
 * });
 */
export class Memmap extends NDArray {
  /** File path (Node.js only) */
  private _filePath: string | null = null;

  /** Original file buffer */
  private _fileBuffer: ArrayBuffer | null = null;

  /** Byte offset where array data starts */
  private _offset: number = 0;

  /** Access mode */
  private _mode: MemmapMode = "r";

  /** Whether data has been modified */
  private _dirty: boolean = false;

  /**
   * Create a new memory-mapped file.
   *
   * @param file - File path (Node.js) or filename for reference
   * @param options - Creation options
   * @returns New Memmap instance
   */
  static async create(
    file: string,
    options: MemmapOptions = {},
  ): Promise<Memmap> {
    const {
      dtype = DType.Float64,
      mode = "w+",
      offset = 0,
      shape = [],
      order: _order = "C",
    } = options;

    if (shape.length === 0) {
      throw new Error("Shape must be specified for new memmap");
    }

    const size = shape.reduce((a, b) => a * b, 1);
    const itemSize = DTYPE_SIZES[dtype];
    const dataSize = size * itemSize;

    // Create buffer for the data
    const buffer = new ArrayBuffer(offset + dataSize);

    // Create base array from buffer
    const typedArray = createTypedArrayForDtype(dtype, buffer, offset, size);
    const baseArr = await NDArray.fromTypedArray(typedArray, shape, dtype);

    // Create Memmap instance
    const mmap = Object.create(Memmap.prototype) as Memmap;
    Object.assign(mmap, baseArr);

    mmap._filePath = file;
    mmap._fileBuffer = buffer;
    mmap._offset = offset;
    mmap._mode = mode;
    mmap._dirty = false;

    // Write initial file if in Node.js
    if (isNode() && (mode === "w+" || mode === "r+")) {
      await mmap.flush();
    }

    return mmap;
  }

  /**
   * Open an existing memory-mapped file.
   *
   * @param file - File path (Node.js) or ArrayBuffer (browser)
   * @param options - Open options
   * @returns Memmap instance
   */
  static async open(
    file: string | ArrayBuffer,
    options: MemmapOptions = {},
  ): Promise<Memmap> {
    const {
      dtype = DType.Float64,
      mode = "r",
      offset = 0,
      shape,
      order: _order = "C",
    } = options;

    let buffer: ArrayBuffer;
    let filePath: string | null = null;

    if (file instanceof ArrayBuffer) {
      buffer = file;
    } else if (typeof file === "string") {
      filePath = file;

      if (isNode()) {
        const fs = await import("fs/promises");
        const fileBuffer = await fs.readFile(file);
        buffer = fileBuffer.buffer.slice(
          fileBuffer.byteOffset,
          fileBuffer.byteOffset + fileBuffer.byteLength,
        );
      } else {
        throw new Error("String file paths only supported in Node.js");
      }
    } else {
      throw new Error("Unsupported file source");
    }

    // Determine shape from file size if not provided
    const itemSize = DTYPE_SIZES[dtype];
    const availableBytes = buffer.byteLength - offset;
    const maxElements = Math.floor(availableBytes / itemSize);

    const actualShape = shape ?? [maxElements];
    const size = actualShape.reduce((a, b) => a * b, 1);

    if (size > maxElements) {
      throw new Error(
        `Shape ${actualShape} requires ${size} elements, ` +
          `but file only has ${maxElements} available`,
      );
    }

    // Create typed array view
    const typedArray = createTypedArrayForDtype(dtype, buffer, offset, size);
    const baseArr = await NDArray.fromTypedArray(
      typedArray,
      actualShape,
      dtype,
    );

    // Create Memmap instance
    const mmap = Object.create(Memmap.prototype) as Memmap;
    Object.assign(mmap, baseArr);

    mmap._filePath = filePath;
    mmap._fileBuffer = buffer;
    mmap._offset = offset;
    mmap._mode = mode;
    mmap._dirty = false;

    return mmap;
  }

  /**
   * Flush changes to disk.
   *
   * Only works in Node.js with 'r+' or 'w+' mode.
   */
  async flush(): Promise<void> {
    if (this._mode === "r" || this._mode === "c") {
      return; // Read-only or copy-on-write, nothing to flush
    }

    if (!this._filePath) {
      throw new Error("Cannot flush: no file path associated");
    }

    if (!isNode()) {
      throw new Error("File flushing only supported in Node.js");
    }

    const fs = await import("fs/promises");

    // Get the current data
    const typedArray = this.toTypedArray();
    const dataBuffer = new Uint8Array(
      typedArray.buffer,
      typedArray.byteOffset,
      typedArray.byteLength,
    );

    // If offset is 0, write entire buffer
    if (this._offset === 0) {
      await fs.writeFile(this._filePath, dataBuffer);
    } else {
      // Need to preserve header/prefix data
      const fullBuffer = new Uint8Array(this._offset + dataBuffer.length);

      if (this._fileBuffer) {
        // Copy existing prefix
        const prefix = new Uint8Array(this._fileBuffer, 0, this._offset);
        fullBuffer.set(prefix, 0);
      }

      fullBuffer.set(dataBuffer, this._offset);
      await fs.writeFile(this._filePath, fullBuffer);
    }

    this._dirty = false;
  }

  /**
   * Get the file path associated with this memmap.
   */
  get filename(): string | null {
    return this._filePath;
  }

  /**
   * Get the access mode.
   */
  get mode(): MemmapMode {
    return this._mode;
  }

  /**
   * Get the byte offset.
   */
  get offset(): number {
    return this._offset;
  }

  /**
   * Check if data has been modified since last flush.
   */
  get dirty(): boolean {
    return this._dirty;
  }

  /**
   * Override set to track modifications.
   */
  set(value: number, ...indices: number[]): void {
    if (this._mode === "r") {
      throw new Error("Cannot modify read-only memmap");
    }
    super.set(value, ...indices);
    this._dirty = true;
  }

  /**
   * Override setFlat to track modifications.
   */
  setFlat(index: number, value: number): void {
    if (this._mode === "r") {
      throw new Error("Cannot modify read-only memmap");
    }
    super.setFlat(index, value);
    this._dirty = true;
  }
}

/**
 * Open a memory-mapped file.
 *
 * Convenience function wrapping Memmap.open().
 *
 * @param file - File path or ArrayBuffer
 * @param options - Memmap options
 * @returns Memmap instance
 *
 * @example
 * // Open existing file
 * const arr = await openMemmap('data.bin', {
 *   dtype: DType.Float32,
 *   shape: [100, 100],
 *   mode: 'r+',
 * });
 *
 * // Modify and flush
 * arr.set(0, 0, 42);
 * await arr.flush();
 */
export async function openMemmap(
  file: string | ArrayBuffer,
  options: MemmapOptions = {},
): Promise<Memmap> {
  return Memmap.open(file, options);
}

/**
 * Create a typed array view for a specific dtype.
 */
function createTypedArrayForDtype(
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
