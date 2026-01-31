/**
 * NumJS Record Arrays - RecArray Class
 *
 * A record array with named fields, compatible with NumPy's recarray.
 * Uses columnar storage where each field is stored as a separate NDArray.
 */

import { NDArray } from "../_core/NDArray.js";
import { DType, type StructuredDType } from "../types.js";
import { format_parser, ValueError } from "./format_parser.js";
import { record, createRecord, type RecArrayLike, KeyError } from "./record.js";

/* ============ Error Classes ============ */

/**
 * Error for out-of-bounds index access.
 */
export class IndexError extends Error {
  constructor(message: string) {
    super(message);
    this.name = "IndexError";
  }
}

/* ============ RecArray Class ============ */

/**
 * Record array with named fields.
 *
 * A recarray allows field access as attributes in addition to
 * standard array indexing. Uses columnar storage internally.
 *
 * @example
 * const arr = rec.fromarrays(
 *   [[25, 35, 45], ['Alice', 'Bob', 'Carol']],
 *   { names: ['age', 'name'] }
 * );
 *
 * arr.age;      // NDArray([25, 35, 45])
 * arr.field('name');  // NDArray(['Alice', 'Bob', 'Carol'])
 * arr.getRecord(0);   // record(age=25, name='Alice')
 */
export class recarray implements RecArrayLike {
  /** Structured dtype describing the fields */
  private _dtype: StructuredDType;

  /** Shape of the array */
  private _shape: number[];

  /** Columnar storage: field name -> NDArray */
  private _fields: Map<string, NDArray>;

  /** Whether this array has been disposed */
  private _disposed: boolean = false;

  /**
   * Private constructor - use static factory methods.
   */
  private constructor(
    dtype: StructuredDType,
    shape: number[],
    fields: Map<string, NDArray>,
  ) {
    this._dtype = dtype;
    this._shape = shape;
    this._fields = fields;
  }

  /* ============ Properties ============ */

  /**
   * Shape of the array (number of records).
   */
  get shape(): number[] {
    this.ensureNotDisposed();
    return [...this._shape];
  }

  /**
   * Number of dimensions.
   */
  get ndim(): number {
    this.ensureNotDisposed();
    return this._shape.length;
  }

  /**
   * Total number of records.
   */
  get size(): number {
    this.ensureNotDisposed();
    return this._shape.reduce((a, b) => a * b, 1);
  }

  /**
   * Structured dtype.
   */
  get dtype(): StructuredDType {
    this.ensureNotDisposed();
    return this._dtype;
  }

  /**
   * Field names.
   */
  get names(): string[] {
    this.ensureNotDisposed();
    return this._dtype.names;
  }

  /**
   * Whether this array has been disposed.
   */
  get isDisposed(): boolean {
    return this._disposed;
  }

  /* ============ Field Access ============ */

  /**
   * Get a field as an NDArray.
   *
   * @param attr - Field name or index
   * @param val - Optional value to set (if provided, sets the field)
   * @returns NDArray containing field values
   */
  field(attr: string | number, val?: NDArray | number[] | string[]): NDArray {
    this.ensureNotDisposed();

    // Handle numeric index
    let fieldName: string;
    if (typeof attr === "number") {
      if (attr < 0 || attr >= this._dtype.fieldList.length) {
        throw new IndexError(`Field index ${attr} out of bounds`);
      }
      fieldName = this._dtype.fieldList[attr].name;
    } else {
      fieldName = attr;
    }

    if (!this._fields.has(fieldName)) {
      throw new KeyError(`Field '${fieldName}' not found`);
    }

    // Set value if provided
    if (val !== undefined) {
      this._setFieldArray(fieldName, val);
    }

    return this._fields.get(fieldName)!;
  }

  /**
   * Set a field array.
   */
  private _setFieldArray(
    name: string,
    value: NDArray | number[] | string[],
  ): void {
    const fieldArr = this._fields.get(name)!;
    const descriptor = this._dtype.fields.get(name)!;

    if (Array.isArray(value)) {
      if (descriptor.dtype === DType.String) {
        // String field
        for (let i = 0; i < value.length && i < fieldArr.size; i++) {
          fieldArr.setStringFlat(i, String(value[i]));
        }
      } else {
        // Numeric field
        for (let i = 0; i < value.length && i < fieldArr.size; i++) {
          fieldArr.setFlat(i, value[i] as number);
        }
      }
    } else {
      // NDArray value
      if (descriptor.dtype === DType.String) {
        for (let i = 0; i < value.size && i < fieldArr.size; i++) {
          fieldArr.setStringFlat(i, value.getStringFlat(i));
        }
      } else {
        for (let i = 0; i < value.size && i < fieldArr.size; i++) {
          fieldArr.setFlat(i, value.getFlat(i));
        }
      }
    }
  }

  /**
   * Get a field value at a specific index.
   * Used by record class.
   */
  getFieldValue(fieldName: string, index: number): unknown {
    this.ensureNotDisposed();

    const fieldArr = this._fields.get(fieldName);
    if (!fieldArr) {
      throw new KeyError(`Field '${fieldName}' not found`);
    }

    const descriptor = this._dtype.fields.get(fieldName)!;

    if (descriptor.dtype === DType.String) {
      return fieldArr.getStringFlat(index);
    } else {
      return fieldArr.getFlat(index);
    }
  }

  /**
   * Set a field value at a specific index.
   * Used by record class.
   */
  setFieldValue(fieldName: string, index: number, value: unknown): void {
    this.ensureNotDisposed();

    const fieldArr = this._fields.get(fieldName);
    if (!fieldArr) {
      throw new KeyError(`Field '${fieldName}' not found`);
    }

    const descriptor = this._dtype.fields.get(fieldName)!;

    if (descriptor.dtype === DType.String) {
      fieldArr.setStringFlat(index, String(value));
    } else {
      fieldArr.setFlat(index, value as number);
    }
  }

  /* ============ Record Access ============ */

  /**
   * Get a single record by index.
   *
   * @param index - Record index (supports negative indexing)
   * @returns A record object
   */
  getRecord(index: number): record {
    this.ensureNotDisposed();

    // Handle negative indexing
    if (index < 0) index += this.size;
    if (index < 0 || index >= this.size) {
      throw new IndexError(
        `Index ${index} out of bounds for size ${this.size}`,
      );
    }

    return createRecord(this, index);
  }

  /**
   * Set a single record by index.
   *
   * @param index - Record index (supports negative indexing)
   * @param values - Array or object of field values
   */
  setRecord(index: number, values: unknown[] | Record<string, unknown>): void {
    this.ensureNotDisposed();

    // Handle negative indexing
    if (index < 0) index += this.size;
    if (index < 0 || index >= this.size) {
      throw new IndexError(
        `Index ${index} out of bounds for size ${this.size}`,
      );
    }

    if (Array.isArray(values)) {
      // Set by position
      for (
        let i = 0;
        i < values.length && i < this._dtype.fieldList.length;
        i++
      ) {
        this.setFieldValue(this._dtype.fieldList[i].name, index, values[i]);
      }
    } else {
      // Set by name
      for (const [name, value] of Object.entries(values)) {
        if (this._dtype.fields.has(name)) {
          this.setFieldValue(name, index, value);
        }
      }
    }
  }

  /* ============ Iteration ============ */

  /**
   * Iterate over records.
   */
  *[Symbol.iterator](): Iterator<record> {
    for (let i = 0; i < this.size; i++) {
      yield this.getRecord(i);
    }
  }

  /* ============ Conversion ============ */

  /**
   * Convert to list of tuples.
   */
  tolist(): unknown[][] {
    this.ensureNotDisposed();

    const result: unknown[][] = [];

    for (let i = 0; i < this.size; i++) {
      const row: unknown[] = [];
      for (const field of this._dtype.fieldList) {
        row.push(this.getFieldValue(field.name, i));
      }
      result.push(row);
    }

    return result;
  }

  /**
   * Convert to list of objects.
   */
  toObjects(): Record<string, unknown>[] {
    this.ensureNotDisposed();

    const result: Record<string, unknown>[] = [];

    for (let i = 0; i < this.size; i++) {
      const obj: Record<string, unknown> = {};
      for (const field of this._dtype.fieldList) {
        obj[field.name] = this.getFieldValue(field.name, i);
      }
      result.push(obj);
    }

    return result;
  }

  /* ============ Copy & Memory ============ */

  /**
   * Create a copy of this recarray.
   */
  copy(): recarray {
    this.ensureNotDisposed();

    const newFields = new Map<string, NDArray>();

    for (const [name, arr] of this._fields) {
      newFields.set(name, arr.copy());
    }

    return new recarray(this._dtype, [...this._shape], newFields);
  }

  /**
   * Free memory for all field arrays.
   */
  dispose(): void {
    if (!this._disposed) {
      for (const arr of this._fields.values()) {
        arr.dispose();
      }
      this._fields.clear();
      this._disposed = true;
    }
  }

  /**
   * String representation.
   */
  toString(): string {
    if (this._disposed) return "recarray(disposed)";

    const lines: string[] = ["recarray(["];

    const maxShow = 10;
    const numToShow = Math.min(this.size, maxShow);

    for (let i = 0; i < numToShow; i++) {
      const rec = this.getRecord(i);
      lines.push(`  ${rec.toString()},`);
    }

    if (this.size > maxShow) {
      lines.push("  ...");
    }

    lines.push(
      `], dtype=[${this._dtype.names.map((n) => `'${n}'`).join(", ")}])`,
    );

    return lines.join("\n");
  }

  /* ============ Private Methods ============ */

  private ensureNotDisposed(): void {
    if (this._disposed) {
      throw new Error("recarray has been disposed");
    }
  }

  /* ============ Static Factory Methods ============ */

  /**
   * Create a recarray from shape and dtype.
   */
  static create(
    shape: number | number[],
    dtype: StructuredDType | null = null,
    formats: string | string[] | null = null,
    names: string | string[] | null = null,
    titles: string | string[] | null = null,
    aligned: boolean = false,
    byteorder: string | null = null,
  ): recarray {
    // Determine structured dtype
    let structDtype: StructuredDType;

    if (dtype !== null) {
      structDtype = dtype;
    } else if (formats !== null) {
      const parser = new format_parser(
        formats,
        names,
        titles,
        aligned,
        byteorder,
      );
      structDtype = parser.dtype;
    } else {
      throw new TypeError("Must specify either dtype or formats");
    }

    // Normalize shape
    const shapeArr = typeof shape === "number" ? [shape] : [...shape];
    const size = shapeArr.reduce((a, b) => a * b, 1);

    // Create field arrays
    const fields = new Map<string, NDArray>();

    for (const field of structDtype.fieldList) {
      if (field.dtype === DType.String) {
        // String field
        fields.set(field.name, NDArray.emptyString(shapeArr));
      } else {
        // Numeric field - use zeros for initialization
        // Note: This is synchronous using the pre-loaded module
        const arr = NDArray.fromStringArray(Array(size).fill(""));
        // Actually we need async zeros, but for now use string placeholder
        // We'll fix this in the index.ts with proper async factory
        fields.set(field.name, arr);
      }
    }

    return new recarray(structDtype, shapeArr, fields);
  }

  /**
   * Create a recarray from column arrays (internal sync version).
   */
  static _fromArraysSync(
    arrayList: NDArray[],
    dtype: StructuredDType,
  ): recarray {
    if (arrayList.length === 0) {
      throw new ValueError("arrayList must contain at least one array");
    }

    if (arrayList.length !== dtype.fieldList.length) {
      throw new ValueError(
        `Number of arrays (${arrayList.length}) must match number of fields (${dtype.fieldList.length})`,
      );
    }

    // Determine shape from first array
    const shape = arrayList[0].shape;

    // Verify all arrays have compatible shapes
    for (let i = 1; i < arrayList.length; i++) {
      const arr = arrayList[i];
      if (arr.size !== arrayList[0].size) {
        throw new ValueError(
          `Array ${i} has size ${arr.size}, expected ${arrayList[0].size}`,
        );
      }
    }

    // Build fields map
    const fields = new Map<string, NDArray>();

    for (let i = 0; i < arrayList.length; i++) {
      const fieldName = dtype.fieldList[i].name;
      fields.set(fieldName, arrayList[i]);
    }

    return new recarray(dtype, shape, fields);
  }
}

/* ============ Proxy Factory ============ */

/**
 * Create a recarray with Proxy for attribute-style field access.
 */
export function createRecarray(
  dtype: StructuredDType,
  shape: number[],
  fields: Map<string, NDArray>,
): recarray {
  // We need to use a factory that can access the private constructor
  // For now, return a Proxy-wrapped recarray using the static method
  const arr = (recarray as any)._createInternal(dtype, shape, fields);

  return new Proxy(arr, {
    get(target, prop) {
      // Handle symbol properties
      if (typeof prop === "symbol") {
        return (target as any)[prop];
      }

      // Allow access to isDisposed even after disposal
      if (prop === "isDisposed") {
        return target.isDisposed;
      }

      // Check if it's a field name (and not a method/property)
      if (
        typeof prop === "string" &&
        !target.isDisposed &&
        target.dtype.fields.has(prop) &&
        !(prop in target)
      ) {
        return target.field(prop);
      }

      // Fall back to actual property
      return (target as any)[prop];
    },

    set(target, prop, value) {
      // Handle symbol properties
      if (typeof prop === "symbol") {
        (target as any)[prop] = value;
        return true;
      }

      // Check if it's a field name
      if (typeof prop === "string" && target.dtype.fields.has(prop)) {
        target.field(prop, value);
        return true;
      }

      // Fall back to actual property
      (target as any)[prop] = value;
      return true;
    },

    has(target, prop) {
      if (typeof prop === "string" && target.dtype.fields.has(prop)) {
        return true;
      }
      return prop in target;
    },
  });
}

// Add internal factory to recarray
(recarray as any)._createInternal = function (
  dtype: StructuredDType,
  shape: number[],
  fields: Map<string, NDArray>,
): recarray {
  return new (recarray as any)(dtype, shape, fields);
};

export { KeyError } from "./record.js";
