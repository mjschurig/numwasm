/**
 * NumJS Record Arrays - Record Class
 *
 * Represents a single record (row) from a record array.
 * Provides field access via attributes and indexing.
 */

import type { StructuredDType } from '../types.js';

/* ============ Types ============ */

/**
 * Interface for the parent recarray (avoids circular import).
 */
export interface RecArrayLike {
  readonly dtype: StructuredDType;
  readonly size: number;
  getFieldValue(fieldName: string, index: number): unknown;
  setFieldValue(fieldName: string, index: number, value: unknown): void;
}

/* ============ Record Class ============ */

/**
 * A single record from a record array.
 *
 * Allows field access via attributes or indexing.
 *
 * @example
 * const rec = recarray.getRecord(0);
 * console.log(rec.name);  // Field access via attribute
 * console.log(rec['age']); // Field access via indexing
 */
export class record {
  /** Reference to parent array */
  private _parent: RecArrayLike;

  /** Index of this record in parent array */
  private _index: number;

  /** Structured dtype describing fields */
  private _dtype: StructuredDType;

  /**
   * Create a record view.
   *
   * @param parent - Parent recarray
   * @param index - Index of this record
   */
  constructor(parent: RecArrayLike, index: number) {
    this._parent = parent;
    this._index = index;
    this._dtype = parent.dtype;

    // Create property getters/setters for each field
    for (const field of this._dtype.fieldList) {
      if (!(field.name in this)) {
        Object.defineProperty(this, field.name, {
          get: () => this.getField(field.name),
          set: (value) => this.setField(field.name, value),
          enumerable: true,
          configurable: true,
        });
      }
    }
  }

  /**
   * Get field value by name.
   */
  getField(name: string): unknown {
    if (!this._dtype.fields.has(name)) {
      throw new KeyError(`Field '${name}' not found`);
    }
    return this._parent.getFieldValue(name, this._index);
  }

  /**
   * Set field value by name.
   */
  setField(name: string, value: unknown): void {
    if (!this._dtype.fields.has(name)) {
      throw new KeyError(`Field '${name}' not found`);
    }
    this._parent.setFieldValue(name, this._index, value);
  }

  /**
   * Number of fields in this record.
   */
  get length(): number {
    return this._dtype.fieldList.length;
  }

  /**
   * Get structured dtype.
   */
  get dtype(): StructuredDType {
    return this._dtype;
  }

  /**
   * Get field names.
   */
  get names(): string[] {
    return this._dtype.names;
  }

  /**
   * Get the record index in the parent array.
   */
  get index(): number {
    return this._index;
  }

  /**
   * Iteration over field values.
   */
  [Symbol.iterator](): Iterator<unknown> {
    let i = 0;
    const fields = this._dtype.fieldList;
    return {
      next: () => {
        if (i < fields.length) {
          return { value: this.getField(fields[i++].name), done: false };
        }
        return { value: undefined, done: true };
      }
    };
  }

  /**
   * Convert to plain JavaScript array (tuple-like).
   */
  toArray(): unknown[] {
    return this._dtype.fieldList.map(f => this.getField(f.name));
  }

  /**
   * Convert to plain JavaScript object.
   */
  toObject(): Record<string, unknown> {
    const result: Record<string, unknown> = {};
    for (const field of this._dtype.fieldList) {
      result[field.name] = this.getField(field.name);
    }
    return result;
  }

  /**
   * String representation.
   */
  toString(): string {
    const values = this._dtype.fieldList.map(f => {
      const val = this.getField(f.name);
      return `${f.name}=${typeof val === 'string' ? `'${val}'` : val}`;
    });
    return `(${values.join(', ')})`;
  }

  /**
   * Pretty print the record.
   */
  pprint(): string {
    const lines: string[] = [];
    const maxNameLen = Math.max(...this._dtype.names.map(n => n.length));

    for (const field of this._dtype.fieldList) {
      const val = this.getField(field.name);
      const paddedName = field.name.padEnd(maxNameLen);
      lines.push(`${paddedName}: ${typeof val === 'string' ? `'${val}'` : val}`);
    }

    return lines.join('\n');
  }
}

/* ============ Proxy Factory ============ */

/**
 * Create a record with Proxy for bracket notation support.
 *
 * This enables both `rec.fieldname` and `rec['fieldname']` access patterns.
 */
export function createRecord(parent: RecArrayLike, index: number): record {
  const rec = new record(parent, index);

  return new Proxy(rec, {
    get(target, prop) {
      // Handle symbol properties
      if (typeof prop === 'symbol') {
        return (target as any)[prop];
      }

      // Check if it's a field name
      if (typeof prop === 'string' && target.dtype.fields.has(prop)) {
        return target.getField(prop);
      }

      // Fall back to actual property
      return (target as any)[prop];
    },

    set(target, prop, value) {
      // Handle symbol properties
      if (typeof prop === 'symbol') {
        (target as any)[prop] = value;
        return true;
      }

      // Check if it's a field name
      if (typeof prop === 'string' && target.dtype.fields.has(prop)) {
        target.setField(prop, value);
        return true;
      }

      // Fall back to actual property
      (target as any)[prop] = value;
      return true;
    },

    has(target, prop) {
      if (typeof prop === 'string' && target.dtype.fields.has(prop)) {
        return true;
      }
      return prop in target;
    },

    ownKeys(target) {
      return [...target.dtype.names, ...Object.keys(target)];
    },

    getOwnPropertyDescriptor(target, prop) {
      if (typeof prop === 'string' && target.dtype.fields.has(prop)) {
        return {
          value: target.getField(prop),
          writable: true,
          enumerable: true,
          configurable: true,
        };
      }
      return Object.getOwnPropertyDescriptor(target, prop);
    }
  });
}

/* ============ Error Classes ============ */

/**
 * Error for missing keys/fields.
 */
export class KeyError extends Error {
  constructor(message: string) {
    super(message);
    this.name = 'KeyError';
  }
}
