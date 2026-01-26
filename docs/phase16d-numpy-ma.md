# Phase 16d: numpy.ma Implementation Plan

Complete implementation roadmap for the NumJS-WASM masked arrays module, providing NumPy-compatible arrays with support for missing or invalid data.

---

## ⚠️ Implementation Guideline

**The original NumPy source code is available in `/numpy`.**

Key reference files:
- `numpy/ma/core.py` - MaskedArray class and core functions (8,940 lines)
- `numpy/ma/extras.py` - Additional utilities (2,266 lines)
- `numpy/ma/mrecords.py` - Masked record arrays (762 lines)
- `numpy/ma/testutils.py` - Testing utilities (294 lines)

Implementation should follow NumPy's algorithms and error handling for consistency.

---

## Current State (Pre-Phase 16d)

```
src/ts/
├── NDArray.ts         # Core array class
├── types.ts           # Type definitions
├── broadcast.ts       # Broadcasting functions
└── index.ts           # Public exports
```

**Dependencies Required:**
- NDArray core
- Broadcasting
- All math operations (Phase 4)
- Statistics functions (Phase 6)

---

## Phase 16d Dependency Tree

```
PHASE 16d: NUMPY.MA
│
├── 16d.1 Core Types and Constants
│   ├── nomask → sentinel for no mask
│   ├── masked → masked scalar constant
│   ├── MaskType → boolean array or nomask
│   └── MaskedArrayError, MaskError
│
├── 16d.2 MaskedArray Class
│   ├── Attributes
│   │   ├── _data → underlying NDArray
│   │   ├── _mask → boolean mask (true=masked)
│   │   ├── _fill_value → value for masked elements
│   │   ├── _hardmask → if true, cannot unmask
│   │   └── _sharedmask → if true, mask is shared
│   │
│   ├── Properties
│   │   ├── data, mask, fill_value
│   │   ├── shape, ndim, size, dtype
│   │   ├── flat (FlatIterator)
│   │   └── T (transpose)
│   │
│   ├── Arithmetic (with mask propagation)
│   │   ├── add, subtract, multiply, divide
│   │   ├── power, mod, floor_divide
│   │   ├── Comparison operators
│   │   └── Bitwise operators
│   │
│   ├── Math Functions (masked versions)
│   │   ├── Trigonometric: sin, cos, tan, etc.
│   │   ├── Hyperbolic: sinh, cosh, tanh, etc.
│   │   ├── Exponential: exp, log, log10, etc.
│   │   └── Rounding: floor, ceil, round, etc.
│   │
│   ├── Reductions (respecting masks)
│   │   ├── sum, prod, mean, std, var
│   │   ├── min, max, ptp, argmin, argmax
│   │   ├── count (non-masked)
│   │   ├── cumsum, cumprod
│   │   └── all, any
│   │
│   ├── Shape Manipulation
│   │   ├── reshape, ravel, flatten
│   │   ├── transpose, swapaxes
│   │   └── squeeze, expand_dims
│   │
│   └── Mask Operations
│       ├── filled(fill_value) → NDArray
│       ├── compressed() → 1D non-masked
│       ├── harden_mask(), soften_mask()
│       └── shrink_mask()
│
├── 16d.3 Array Creation Functions
│   ├── masked_array(data, mask)
│   ├── array(data, mask, dtype)
│   ├── masked_equal(x, value)
│   ├── masked_greater(x, value)
│   ├── masked_less(x, value)
│   ├── masked_inside(x, v1, v2)
│   ├── masked_outside(x, v1, v2)
│   ├── masked_where(condition, x)
│   ├── masked_invalid(x) → mask NaN/Inf
│   ├── masked_values(x, value, rtol, atol)
│   └── zeros, ones, empty (masked versions)
│
├── 16d.4 Mask Operations
│   ├── make_mask(m, copy, shrink)
│   ├── make_mask_none(shape)
│   ├── make_mask_descr(dtype)
│   ├── getmask(a), getmaskarray(a)
│   ├── getdata(a)
│   ├── is_mask(m), is_masked(x)
│   ├── mask_or(m1, m2)
│   └── flatten_mask(mask)
│
├── 16d.5 Fill Value Operations
│   ├── common_fill_value(a, b)
│   ├── default_fill_value(obj)
│   ├── maximum_fill_value(dtype)
│   ├── minimum_fill_value(dtype)
│   └── set_fill_value(a, fill_value)
│
├── 16d.6 Extras Module
│   ├── Statistical Functions
│   │   ├── average(a, axis, weights)
│   │   ├── median(a, axis)
│   │   ├── cov(x, y, rowvar)
│   │   └── corrcoef(x, y, rowvar)
│   │
│   ├── Array Manipulation
│   │   ├── atleast_1d, atleast_2d, atleast_3d
│   │   ├── apply_along_axis(func, axis, arr)
│   │   └── apply_over_axes(func, arr, axes)
│   │
│   ├── Set Operations
│   │   ├── unique(ar, return_index, return_inverse)
│   │   ├── intersect1d, union1d
│   │   ├── setdiff1d, setxor1d
│   │   └── in1d, isin
│   │
│   ├── Edge Detection
│   │   ├── notmasked_edges(a, axis)
│   │   ├── notmasked_contiguous(a, axis)
│   │   ├── flatnotmasked_edges(a)
│   │   └── flatnotmasked_contiguous(a)
│   │
│   └── Clumping
│       ├── clump_masked(a)
│       └── clump_unmasked(a)
│
└── 16d.7 Masked Records (mrecords)
    ├── mrecarray class
    ├── MaskedRecord class
    └── fromrecords, fromarrays

Dependencies: NDArray, Broadcasting, Math ops, Statistics
```

---

## Detailed Implementation Specifications

### 16d.1 Core Types and Constants

**File:** `src/ts/ma/types.ts`

```typescript
import { NDArray } from '../NDArray.js';
import { DType } from '../types.js';

/**
 * Sentinel value indicating no mask.
 * When mask is nomask, all elements are valid.
 */
export const nomask: unique symbol = Symbol('nomask');

/**
 * Type for mask: either a boolean array or nomask sentinel.
 */
export type MaskType = NDArray | typeof nomask;

/**
 * Masked scalar constant.
 * Represents a single masked value.
 */
export const masked = Object.freeze({
  __masked__: true,
  toString: () => '--',
  valueOf: () => NaN,
});

/**
 * Check if value is the masked constant.
 */
export function isMaskedConstant(x: any): boolean {
  return x !== null && typeof x === 'object' && x.__masked__ === true;
}

/**
 * Error for masked array operations.
 */
export class MaskedArrayError extends Error {
  constructor(message: string) {
    super(message);
    this.name = 'MaskedArrayError';
  }
}

/**
 * Error for mask operations.
 */
export class MaskError extends Error {
  constructor(message: string) {
    super(message);
    this.name = 'MaskError';
  }
}

/**
 * Default fill values for each dtype.
 */
export const defaultFillValues: Record<DType, number | boolean | string> = {
  [DType.Bool]: true,
  [DType.Int8]: 999999,
  [DType.Int16]: 999999,
  [DType.Int32]: 999999,
  [DType.Int64]: 999999,
  [DType.UInt8]: 999999,
  [DType.UInt16]: 999999,
  [DType.UInt32]: 999999,
  [DType.UInt64]: 999999,
  [DType.Float16]: 1e20,
  [DType.Float32]: 1e20,
  [DType.Float64]: 1e20,
  [DType.Complex64]: 1e20,
  [DType.Complex128]: 1e20,
  [DType.String]: 'N/A',
};

/**
 * Maximum fill values for dtype (for min operations).
 */
export function maximumFillValue(dtype: DType): number {
  switch (dtype) {
    case DType.Int8: return 127;
    case DType.Int16: return 32767;
    case DType.Int32: return 2147483647;
    case DType.UInt8: return 255;
    case DType.UInt16: return 65535;
    case DType.UInt32: return 4294967295;
    case DType.Float16:
    case DType.Float32:
    case DType.Float64:
      return Number.POSITIVE_INFINITY;
    default:
      return 1e20;
  }
}

/**
 * Minimum fill values for dtype (for max operations).
 */
export function minimumFillValue(dtype: DType): number {
  switch (dtype) {
    case DType.Int8: return -128;
    case DType.Int16: return -32768;
    case DType.Int32: return -2147483648;
    case DType.UInt8:
    case DType.UInt16:
    case DType.UInt32:
      return 0;
    case DType.Float16:
    case DType.Float32:
    case DType.Float64:
      return Number.NEGATIVE_INFINITY;
    default:
      return -1e20;
  }
}
```

---

### 16d.2 MaskedArray Class

**File:** `src/ts/ma/core.ts`

```typescript
import { NDArray } from '../NDArray.js';
import { DType } from '../types.js';
import { broadcastArrays, broadcastShapes } from '../broadcast.js';
import {
  nomask,
  masked,
  MaskType,
  MaskedArrayError,
  MaskError,
  defaultFillValues,
  maximumFillValue,
  minimumFillValue,
} from './types.js';

/**
 * Array with masked (invalid) elements.
 *
 * A MaskedArray is a combination of a standard NDArray with a boolean
 * mask that indicates which elements should be excluded from operations.
 *
 * @example
 * const data = NDArray.fromArray([1, 2, 3, 4, 5]);
 * const mask = [false, false, true, false, false];
 * const ma = new MaskedArray(data, mask);
 *
 * ma.mean(); // (1+2+4+5)/4 = 3, excludes masked element
 * ma.filled(0); // [1, 2, 0, 4, 5]
 */
export class MaskedArray {
  /** Underlying data array */
  protected _data: NDArray;

  /** Boolean mask (true = masked/invalid) */
  protected _mask: MaskType;

  /** Value to use for masked elements in filled() */
  protected _fill_value: number | string | boolean;

  /** If true, mask cannot be changed (unmasked) */
  protected _hardmask: boolean;

  /** If true, mask is shared with parent view */
  protected _sharedmask: boolean;

  /** Base data offset for views */
  protected _baseclass: typeof MaskedArray = MaskedArray;

  /**
   * Create a masked array.
   *
   * @param data - Input data (array or MaskedArray)
   * @param mask - Boolean mask or nomask
   * @param dtype - Data type
   * @param copy - If true, copy data
   * @param fill_value - Fill value for masked elements
   * @param hard_mask - If true, mask cannot be unset
   * @param shrink - If true, shrink mask to nomask if all false
   */
  constructor(
    data: NDArray | number[] | MaskedArray,
    mask: MaskType | boolean[] | boolean = nomask,
    dtype: DType | null = null,
    copy: boolean = false,
    fill_value: number | string | boolean | null = null,
    hard_mask: boolean = false,
    shrink: boolean = true
  ) {
    // Extract data and existing mask from MaskedArray input
    if (data instanceof MaskedArray) {
      this._data = copy ? data._data.copy() : data._data;
      const existingMask = data._mask;

      // Combine masks if both exist
      if (mask !== nomask && existingMask !== nomask) {
        this._mask = mask_or(existingMask, this._processMask(mask, data.shape, shrink));
      } else if (mask !== nomask) {
        this._mask = this._processMask(mask, data.shape, shrink);
      } else {
        this._mask = existingMask;
      }

      this._fill_value = fill_value ?? data._fill_value;
      this._hardmask = hard_mask || data._hardmask;
    } else {
      // Convert to NDArray
      if (Array.isArray(data)) {
        this._data = NDArray.fromArray(data, dtype ?? DType.Float64);
      } else if (copy) {
        this._data = data.copy();
      } else {
        this._data = data;
      }

      // Apply dtype conversion
      if (dtype !== null && this._data.dtype !== dtype) {
        this._data = this._data.astype(dtype);
      }

      this._mask = this._processMask(mask, this._data.shape, shrink);
      this._fill_value = fill_value ?? defaultFillValues[this._data.dtype];
      this._hardmask = hard_mask;
    }

    this._sharedmask = false;
  }

  /* ============ Properties ============ */

  /**
   * The underlying data array.
   */
  get data(): NDArray {
    return this._data;
  }

  /**
   * The mask array (or nomask if no masking).
   */
  get mask(): MaskType {
    return this._mask;
  }

  /**
   * Set the mask.
   */
  set mask(value: MaskType | boolean[] | boolean) {
    if (this._hardmask && this._mask !== nomask) {
      // Hard mask: can only add to mask, not remove
      const newMask = this._processMask(value, this._data.shape, false);
      if (newMask === nomask) {
        throw new MaskError('Cannot unmask hard mask');
      }
      this._mask = mask_or(this._mask, newMask);
    } else {
      this._mask = this._processMask(value, this._data.shape, true);
    }
  }

  /**
   * Fill value for masked elements.
   */
  get fill_value(): number | string | boolean {
    return this._fill_value;
  }

  set fill_value(value: number | string | boolean) {
    this._fill_value = value;
  }

  /**
   * Shape of the array.
   */
  get shape(): number[] {
    return this._data.shape;
  }

  /**
   * Number of dimensions.
   */
  get ndim(): number {
    return this._data.ndim;
  }

  /**
   * Total number of elements.
   */
  get size(): number {
    return this._data.size;
  }

  /**
   * Data type.
   */
  get dtype(): DType {
    return this._data.dtype;
  }

  /**
   * Transpose.
   */
  get T(): MaskedArray {
    return this.transpose();
  }

  /* ============ Element Access ============ */

  /**
   * Get element at flat index.
   */
  getFlat(index: number): number | typeof masked {
    if (this._mask !== nomask && (this._mask as NDArray).getFlat(index)) {
      return masked;
    }
    return this._data.getFlat(index);
  }

  /**
   * Set element at flat index.
   */
  setFlat(index: number, value: number | typeof masked): void {
    if (value === masked || isMaskedConstant(value)) {
      if (this._mask === nomask) {
        this._mask = NDArray.zeros(this._data.shape, DType.Bool);
      }
      (this._mask as NDArray).setFlat(index, 1);
    } else {
      this._data.setFlat(index, value as number);
      if (this._mask !== nomask && !this._hardmask) {
        (this._mask as NDArray).setFlat(index, 0);
      }
    }
  }

  /* ============ Mask Operations ============ */

  /**
   * Count of non-masked elements.
   *
   * @param axis - Axis along which to count
   * @returns Count as number or NDArray
   */
  count(axis: number | null = null): number | NDArray {
    if (this._mask === nomask) {
      if (axis === null) {
        return this._data.size;
      }
      // Return array of sizes along axis
      const shape = [...this._data.shape];
      const axisSize = shape[axis];
      shape.splice(axis, 1);
      return NDArray.full(shape, axisSize, DType.Int32);
    }

    // Count non-masked (mask == false)
    const notMask = this._logicalNot(this._mask as NDArray);

    if (axis === null) {
      return notMask.sum() as number;
    }

    return notMask.sum(axis);
  }

  /**
   * Return data with masked values replaced by fill_value.
   *
   * @param fill_value - Value to use (default: array's fill_value)
   * @returns Regular NDArray with filled values
   */
  filled(fill_value: number | string | boolean | null = null): NDArray {
    const fv = fill_value ?? this._fill_value;

    if (this._mask === nomask) {
      return this._data.copy();
    }

    const result = this._data.copy();
    const mask = this._mask as NDArray;

    for (let i = 0; i < result.size; i++) {
      if (mask.getFlat(i)) {
        result.setFlat(i, fv as number);
      }
    }

    return result;
  }

  /**
   * Return 1D array of non-masked values.
   */
  compressed(): NDArray {
    if (this._mask === nomask) {
      return this._data.ravel();
    }

    const mask = this._mask as NDArray;
    const values: number[] = [];

    for (let i = 0; i < this._data.size; i++) {
      if (!mask.getFlat(i)) {
        values.push(this._data.getFlat(i));
      }
    }

    return NDArray.fromArray(values, this._data.dtype);
  }

  /**
   * Force the mask to hard (cannot unmask).
   */
  harden_mask(): this {
    this._hardmask = true;
    return this;
  }

  /**
   * Force the mask to soft (can unmask).
   */
  soften_mask(): this {
    this._hardmask = false;
    return this;
  }

  /**
   * Reduce mask to nomask if all elements are False.
   */
  shrink_mask(): this {
    if (this._mask !== nomask) {
      const mask = this._mask as NDArray;
      let anyMasked = false;
      for (let i = 0; i < mask.size; i++) {
        if (mask.getFlat(i)) {
          anyMasked = true;
          break;
        }
      }
      if (!anyMasked) {
        this._mask = nomask;
      }
    }
    return this;
  }

  /* ============ Arithmetic Operations ============ */

  /**
   * Add another array or scalar.
   */
  add(other: MaskedArray | NDArray | number): MaskedArray {
    return this._binaryOp(other, (a, b) => a + b);
  }

  /**
   * Subtract another array or scalar.
   */
  subtract(other: MaskedArray | NDArray | number): MaskedArray {
    return this._binaryOp(other, (a, b) => a - b);
  }

  /**
   * Multiply by another array or scalar.
   */
  multiply(other: MaskedArray | NDArray | number): MaskedArray {
    return this._binaryOp(other, (a, b) => a * b);
  }

  /**
   * Divide by another array or scalar.
   */
  divide(other: MaskedArray | NDArray | number): MaskedArray {
    // Also mask division by zero
    return this._binaryOp(other, (a, b) => a / b, (a, b) => b === 0);
  }

  /**
   * True divide.
   */
  true_divide(other: MaskedArray | NDArray | number): MaskedArray {
    return this.divide(other);
  }

  /**
   * Floor divide.
   */
  floor_divide(other: MaskedArray | NDArray | number): MaskedArray {
    return this._binaryOp(other, (a, b) => Math.floor(a / b), (a, b) => b === 0);
  }

  /**
   * Power.
   */
  power(other: MaskedArray | NDArray | number): MaskedArray {
    return this._binaryOp(other, (a, b) => Math.pow(a, b));
  }

  /**
   * Modulo.
   */
  mod(other: MaskedArray | NDArray | number): MaskedArray {
    return this._binaryOp(other, (a, b) => a % b, (a, b) => b === 0);
  }

  /**
   * Negate.
   */
  neg(): MaskedArray {
    return this._unaryOp(x => -x);
  }

  /**
   * Absolute value.
   */
  abs(): MaskedArray {
    return this._unaryOp(Math.abs);
  }

  /* ============ Comparison Operations ============ */

  /**
   * Equal comparison.
   */
  equal(other: MaskedArray | NDArray | number): MaskedArray {
    return this._binaryOp(other, (a, b) => a === b ? 1 : 0);
  }

  /**
   * Not equal comparison.
   */
  not_equal(other: MaskedArray | NDArray | number): MaskedArray {
    return this._binaryOp(other, (a, b) => a !== b ? 1 : 0);
  }

  /**
   * Less than comparison.
   */
  less(other: MaskedArray | NDArray | number): MaskedArray {
    return this._binaryOp(other, (a, b) => a < b ? 1 : 0);
  }

  /**
   * Less than or equal comparison.
   */
  less_equal(other: MaskedArray | NDArray | number): MaskedArray {
    return this._binaryOp(other, (a, b) => a <= b ? 1 : 0);
  }

  /**
   * Greater than comparison.
   */
  greater(other: MaskedArray | NDArray | number): MaskedArray {
    return this._binaryOp(other, (a, b) => a > b ? 1 : 0);
  }

  /**
   * Greater than or equal comparison.
   */
  greater_equal(other: MaskedArray | NDArray | number): MaskedArray {
    return this._binaryOp(other, (a, b) => a >= b ? 1 : 0);
  }

  /* ============ Math Functions ============ */

  sin(): MaskedArray { return this._unaryOp(Math.sin); }
  cos(): MaskedArray { return this._unaryOp(Math.cos); }
  tan(): MaskedArray { return this._unaryOp(Math.tan); }
  arcsin(): MaskedArray { return this._unaryOp(Math.asin, x => x < -1 || x > 1); }
  arccos(): MaskedArray { return this._unaryOp(Math.acos, x => x < -1 || x > 1); }
  arctan(): MaskedArray { return this._unaryOp(Math.atan); }

  sinh(): MaskedArray { return this._unaryOp(Math.sinh); }
  cosh(): MaskedArray { return this._unaryOp(Math.cosh); }
  tanh(): MaskedArray { return this._unaryOp(Math.tanh); }

  exp(): MaskedArray { return this._unaryOp(Math.exp); }
  log(): MaskedArray { return this._unaryOp(Math.log, x => x <= 0); }
  log10(): MaskedArray { return this._unaryOp(Math.log10, x => x <= 0); }
  log2(): MaskedArray { return this._unaryOp(Math.log2, x => x <= 0); }

  sqrt(): MaskedArray { return this._unaryOp(Math.sqrt, x => x < 0); }
  square(): MaskedArray { return this._unaryOp(x => x * x); }

  floor(): MaskedArray { return this._unaryOp(Math.floor); }
  ceil(): MaskedArray { return this._unaryOp(Math.ceil); }
  round(): MaskedArray { return this._unaryOp(Math.round); }
  trunc(): MaskedArray { return this._unaryOp(Math.trunc); }

  /* ============ Reductions ============ */

  /**
   * Sum of non-masked elements.
   */
  sum(axis: number | null = null, keepdims: boolean = false): number | MaskedArray {
    return this._reduction('sum', axis, keepdims);
  }

  /**
   * Product of non-masked elements.
   */
  prod(axis: number | null = null, keepdims: boolean = false): number | MaskedArray {
    return this._reduction('prod', axis, keepdims);
  }

  /**
   * Mean of non-masked elements.
   */
  mean(axis: number | null = null, keepdims: boolean = false): number | MaskedArray {
    const sum = this.sum(axis, keepdims);
    const count = this.count(axis);

    if (typeof sum === 'number' && typeof count === 'number') {
      return count === 0 ? NaN : sum / count;
    }

    // Element-wise division
    return this._divideResults(sum as MaskedArray, count as NDArray);
  }

  /**
   * Variance of non-masked elements.
   */
  var(axis: number | null = null, ddof: number = 0, keepdims: boolean = false): number | MaskedArray {
    const mean = this.mean(axis, true);
    const centered = this.subtract(mean as MaskedArray);
    const squared = centered.multiply(centered);
    const sumSquared = squared.sum(axis, keepdims);
    const count = this.count(axis);

    if (typeof sumSquared === 'number' && typeof count === 'number') {
      const n = count - ddof;
      return n <= 0 ? NaN : sumSquared / n;
    }

    return this._divideWithDdof(sumSquared as MaskedArray, count as NDArray, ddof);
  }

  /**
   * Standard deviation of non-masked elements.
   */
  std(axis: number | null = null, ddof: number = 0, keepdims: boolean = false): number | MaskedArray {
    const variance = this.var(axis, ddof, keepdims);

    if (typeof variance === 'number') {
      return Math.sqrt(variance);
    }

    return this._sqrtResult(variance);
  }

  /**
   * Minimum of non-masked elements.
   */
  min(axis: number | null = null, keepdims: boolean = false): number | MaskedArray {
    return this._reduction('min', axis, keepdims);
  }

  /**
   * Maximum of non-masked elements.
   */
  max(axis: number | null = null, keepdims: boolean = false): number | MaskedArray {
    return this._reduction('max', axis, keepdims);
  }

  /**
   * Peak-to-peak (max - min) of non-masked elements.
   */
  ptp(axis: number | null = null, keepdims: boolean = false): number | MaskedArray {
    const minVal = this.min(axis, keepdims);
    const maxVal = this.max(axis, keepdims);

    if (typeof minVal === 'number' && typeof maxVal === 'number') {
      return maxVal - minVal;
    }

    return (maxVal as MaskedArray).subtract(minVal as MaskedArray);
  }

  /**
   * Index of minimum non-masked element.
   */
  argmin(axis: number | null = null): number | NDArray {
    return this._argReduction('min', axis);
  }

  /**
   * Index of maximum non-masked element.
   */
  argmax(axis: number | null = null): number | NDArray {
    return this._argReduction('max', axis);
  }

  /**
   * Cumulative sum.
   */
  cumsum(axis: number | null = null): MaskedArray {
    return this._cumulative('sum', axis);
  }

  /**
   * Cumulative product.
   */
  cumprod(axis: number | null = null): MaskedArray {
    return this._cumulative('prod', axis);
  }

  /**
   * Test if all non-masked elements are true.
   */
  all(axis: number | null = null): boolean | NDArray {
    return this._booleanReduction('all', axis);
  }

  /**
   * Test if any non-masked element is true.
   */
  any(axis: number | null = null): boolean | NDArray {
    return this._booleanReduction('any', axis);
  }

  /* ============ Shape Manipulation ============ */

  /**
   * Reshape the array and mask.
   */
  reshape(newshape: number[]): MaskedArray {
    const newData = this._data.reshape(newshape);
    const newMask = this._mask === nomask
      ? nomask
      : (this._mask as NDArray).reshape(newshape);

    return new MaskedArray(newData, newMask, null, false, this._fill_value, this._hardmask);
  }

  /**
   * Flatten to 1D (view if possible).
   */
  ravel(): MaskedArray {
    return this.reshape([this._data.size]);
  }

  /**
   * Flatten to 1D (always copy).
   */
  flatten(): MaskedArray {
    return new MaskedArray(
      this._data.flatten(),
      this._mask === nomask ? nomask : (this._mask as NDArray).flatten(),
      null,
      false,
      this._fill_value,
      this._hardmask
    );
  }

  /**
   * Transpose the array.
   */
  transpose(axes: number[] | null = null): MaskedArray {
    const newData = this._data.transpose(axes);
    const newMask = this._mask === nomask
      ? nomask
      : (this._mask as NDArray).transpose(axes);

    return new MaskedArray(newData, newMask, null, false, this._fill_value, this._hardmask);
  }

  /**
   * Swap two axes.
   */
  swapaxes(axis1: number, axis2: number): MaskedArray {
    const newData = this._data.swapaxes(axis1, axis2);
    const newMask = this._mask === nomask
      ? nomask
      : (this._mask as NDArray).swapaxes(axis1, axis2);

    return new MaskedArray(newData, newMask, null, false, this._fill_value, this._hardmask);
  }

  /**
   * Remove axes of length 1.
   */
  squeeze(axis: number | null = null): MaskedArray {
    const newData = this._data.squeeze(axis);
    const newMask = this._mask === nomask
      ? nomask
      : (this._mask as NDArray).squeeze(axis);

    return new MaskedArray(newData, newMask, null, false, this._fill_value, this._hardmask);
  }

  /**
   * Add axis of length 1.
   */
  expand_dims(axis: number): MaskedArray {
    const newData = this._data.expand_dims(axis);
    const newMask = this._mask === nomask
      ? nomask
      : (this._mask as NDArray).expand_dims(axis);

    return new MaskedArray(newData, newMask, null, false, this._fill_value, this._hardmask);
  }

  /**
   * Create a copy.
   */
  copy(): MaskedArray {
    return new MaskedArray(
      this._data.copy(),
      this._mask === nomask ? nomask : (this._mask as NDArray).copy(),
      null,
      false,
      this._fill_value,
      this._hardmask
    );
  }

  /**
   * Convert to specified dtype.
   */
  astype(dtype: DType): MaskedArray {
    return new MaskedArray(
      this._data.astype(dtype),
      this._mask,
      null,
      false,
      this._fill_value,
      this._hardmask
    );
  }

  /* ============ String Representation ============ */

  toString(): string {
    const lines: string[] = [];
    const maxShow = 10;
    const numToShow = Math.min(this._data.size, maxShow);

    for (let i = 0; i < numToShow; i++) {
      const val = this.getFlat(i);
      lines.push(val === masked ? '--' : String(val));
    }

    if (this._data.size > maxShow) {
      lines.push('...');
    }

    return `masked_array([${lines.join(', ')}], fill_value=${this._fill_value})`;
  }

  /* ============ Private Helpers ============ */

  /**
   * Process mask input into standard form.
   */
  protected _processMask(
    mask: MaskType | boolean[] | boolean,
    shape: number[],
    shrink: boolean
  ): MaskType {
    if (mask === nomask) {
      return nomask;
    }

    if (typeof mask === 'boolean') {
      if (!mask) return nomask;
      return NDArray.full(shape, 1, DType.Bool);
    }

    let maskArr: NDArray;
    if (Array.isArray(mask)) {
      maskArr = NDArray.fromArray(mask.map(m => m ? 1 : 0), DType.Bool);
    } else {
      maskArr = mask.astype(DType.Bool);
    }

    // Broadcast mask to shape if needed
    if (!this._shapesEqual(maskArr.shape, shape)) {
      const [broadcastedMask] = broadcastArrays(maskArr, NDArray.zeros(shape, DType.Bool));
      maskArr = broadcastedMask;
    }

    // Shrink to nomask if all false
    if (shrink && !this._anyTrue(maskArr)) {
      return nomask;
    }

    return maskArr;
  }

  /**
   * Apply binary operation with mask propagation.
   */
  protected _binaryOp(
    other: MaskedArray | NDArray | number,
    op: (a: number, b: number) => number,
    extraMask?: (a: number, b: number) => boolean
  ): MaskedArray {
    let otherData: NDArray;
    let otherMask: MaskType;

    if (typeof other === 'number') {
      otherData = NDArray.full(this._data.shape, other, this._data.dtype);
      otherMask = nomask;
    } else if (other instanceof MaskedArray) {
      otherData = other._data;
      otherMask = other._mask;
    } else {
      otherData = other;
      otherMask = nomask;
    }

    // Broadcast
    const [bData1, bData2] = broadcastArrays(this._data, otherData);

    // Compute result
    const resultData = NDArray.empty(bData1.shape, this._data.dtype);
    const resultMaskArr = NDArray.zeros(bData1.shape, DType.Bool);

    for (let i = 0; i < resultData.size; i++) {
      const a = bData1.getFlat(i);
      const b = bData2.getFlat(i);
      resultData.setFlat(i, op(a, b));

      // Extra condition for masking (e.g., division by zero)
      if (extraMask && extraMask(a, b)) {
        resultMaskArr.setFlat(i, 1);
      }
    }

    // Combine masks
    const resultMask = this._combineMasks(this._mask, otherMask, resultMaskArr, bData1.shape);

    return new MaskedArray(resultData, resultMask, null, false, this._fill_value);
  }

  /**
   * Apply unary operation.
   */
  protected _unaryOp(
    op: (x: number) => number,
    shouldMask?: (x: number) => boolean
  ): MaskedArray {
    const resultData = NDArray.empty(this._data.shape, this._data.dtype);
    let resultMask = this._mask === nomask ? nomask : (this._mask as NDArray).copy();

    for (let i = 0; i < resultData.size; i++) {
      const x = this._data.getFlat(i);
      resultData.setFlat(i, op(x));

      if (shouldMask && shouldMask(x)) {
        if (resultMask === nomask) {
          resultMask = NDArray.zeros(this._data.shape, DType.Bool);
        }
        (resultMask as NDArray).setFlat(i, 1);
      }
    }

    return new MaskedArray(resultData, resultMask, null, false, this._fill_value);
  }

  /**
   * Perform reduction operation.
   */
  protected _reduction(
    type: 'sum' | 'prod' | 'min' | 'max',
    axis: number | null,
    keepdims: boolean
  ): number | MaskedArray {
    // Get fill value for this reduction type
    const fillValue = type === 'sum' ? 0 :
                      type === 'prod' ? 1 :
                      type === 'min' ? maximumFillValue(this._data.dtype) :
                      minimumFillValue(this._data.dtype);

    // Fill masked values
    const filledData = this.filled(fillValue);

    // Perform reduction
    let result: number | NDArray;
    switch (type) {
      case 'sum':
        result = filledData.sum(axis, keepdims);
        break;
      case 'prod':
        result = filledData.prod(axis, keepdims);
        break;
      case 'min':
        result = filledData.min(axis, keepdims);
        break;
      case 'max':
        result = filledData.max(axis, keepdims);
        break;
    }

    if (typeof result === 'number') {
      return result;
    }

    // For axis reductions, we need to check if all values along axis were masked
    // and mask those results
    return new MaskedArray(result as NDArray);
  }

  // ... additional helper methods ...

  protected _combineMasks(
    mask1: MaskType,
    mask2: MaskType,
    extraMask: NDArray | null,
    shape: number[]
  ): MaskType {
    const masks: NDArray[] = [];

    if (mask1 !== nomask) {
      const [bm1] = broadcastArrays(mask1 as NDArray, NDArray.zeros(shape, DType.Bool));
      masks.push(bm1);
    }
    if (mask2 !== nomask) {
      const [bm2] = broadcastArrays(mask2 as NDArray, NDArray.zeros(shape, DType.Bool));
      masks.push(bm2);
    }
    if (extraMask !== null) {
      masks.push(extraMask);
    }

    if (masks.length === 0) return nomask;

    // OR all masks
    const result = NDArray.zeros(shape, DType.Bool);
    for (const m of masks) {
      for (let i = 0; i < result.size; i++) {
        if (m.getFlat(i)) result.setFlat(i, 1);
      }
    }

    // Shrink if all false
    if (!this._anyTrue(result)) return nomask;
    return result;
  }

  protected _anyTrue(arr: NDArray): boolean {
    for (let i = 0; i < arr.size; i++) {
      if (arr.getFlat(i)) return true;
    }
    return false;
  }

  protected _logicalNot(arr: NDArray): NDArray {
    const result = NDArray.empty(arr.shape, DType.Bool);
    for (let i = 0; i < arr.size; i++) {
      result.setFlat(i, arr.getFlat(i) ? 0 : 1);
    }
    return result;
  }

  protected _shapesEqual(s1: number[], s2: number[]): boolean {
    if (s1.length !== s2.length) return false;
    return s1.every((v, i) => v === s2[i]);
  }

  protected _divideResults(a: MaskedArray, b: NDArray): MaskedArray {
    return a.divide(b);
  }

  protected _divideWithDdof(a: MaskedArray, b: NDArray, ddof: number): MaskedArray {
    const adjusted = NDArray.empty(b.shape, b.dtype);
    for (let i = 0; i < b.size; i++) {
      adjusted.setFlat(i, b.getFlat(i) - ddof);
    }
    return a.divide(adjusted);
  }

  protected _sqrtResult(a: MaskedArray): MaskedArray {
    return a.sqrt();
  }

  protected _argReduction(type: 'min' | 'max', axis: number | null): number | NDArray {
    // Simplified implementation
    if (axis !== null) {
      throw new MaskedArrayError('argmin/argmax with axis not yet implemented');
    }

    let bestIdx = -1;
    let bestVal = type === 'min' ? Infinity : -Infinity;

    for (let i = 0; i < this._data.size; i++) {
      if (this._mask !== nomask && (this._mask as NDArray).getFlat(i)) {
        continue;
      }
      const val = this._data.getFlat(i);
      if (type === 'min' ? val < bestVal : val > bestVal) {
        bestVal = val;
        bestIdx = i;
      }
    }

    return bestIdx;
  }

  protected _cumulative(type: 'sum' | 'prod', axis: number | null): MaskedArray {
    if (axis !== null) {
      throw new MaskedArrayError('cumsum/cumprod with axis not yet implemented');
    }

    const result = NDArray.empty([this._data.size], this._data.dtype);
    const resultMask = this._mask === nomask
      ? nomask
      : NDArray.zeros([this._data.size], DType.Bool);

    let acc = type === 'sum' ? 0 : 1;

    for (let i = 0; i < this._data.size; i++) {
      if (this._mask !== nomask && (this._mask as NDArray).getFlat(i)) {
        if (resultMask !== nomask) {
          (resultMask as NDArray).setFlat(i, 1);
        }
        result.setFlat(i, acc);
      } else {
        const val = this._data.getFlat(i);
        acc = type === 'sum' ? acc + val : acc * val;
        result.setFlat(i, acc);
      }
    }

    return new MaskedArray(result, resultMask);
  }

  protected _booleanReduction(type: 'all' | 'any', axis: number | null): boolean | NDArray {
    if (axis !== null) {
      throw new MaskedArrayError('all/any with axis not yet implemented');
    }

    for (let i = 0; i < this._data.size; i++) {
      if (this._mask !== nomask && (this._mask as NDArray).getFlat(i)) {
        continue;
      }
      const val = this._data.getFlat(i);
      if (type === 'all' && !val) return false;
      if (type === 'any' && val) return true;
    }

    return type === 'all';
  }
}

/* ============ Utility Functions ============ */

/**
 * Logical OR of two masks.
 */
export function mask_or(m1: MaskType, m2: MaskType): MaskType {
  if (m1 === nomask && m2 === nomask) return nomask;
  if (m1 === nomask) return m2;
  if (m2 === nomask) return m1;

  const arr1 = m1 as NDArray;
  const arr2 = m2 as NDArray;
  const [b1, b2] = broadcastArrays(arr1, arr2);

  const result = NDArray.zeros(b1.shape, DType.Bool);
  for (let i = 0; i < result.size; i++) {
    result.setFlat(i, b1.getFlat(i) || b2.getFlat(i) ? 1 : 0);
  }

  return result;
}

function isMaskedConstant(x: any): boolean {
  return x !== null && typeof x === 'object' && x.__masked__ === true;
}
```

---

### 16d.3 Array Creation Functions

**File:** `src/ts/ma/creation.ts`

```typescript
import { NDArray } from '../NDArray.js';
import { DType } from '../types.js';
import { MaskedArray, mask_or } from './core.js';
import { nomask, MaskType } from './types.js';

/**
 * Create a masked array.
 *
 * @param data - Input data
 * @param mask - Boolean mask
 * @param dtype - Data type
 * @param copy - If true, copy data
 * @param fill_value - Fill value for masked elements
 */
export function masked_array(
  data: NDArray | number[] | MaskedArray,
  mask: MaskType | boolean[] | boolean = nomask,
  dtype: DType | null = null,
  copy: boolean = false,
  fill_value: number | null = null
): MaskedArray {
  return new MaskedArray(data, mask, dtype, copy, fill_value);
}

/**
 * Alias for masked_array.
 */
export const array = masked_array;

/**
 * Mask elements equal to a value.
 */
export function masked_equal(x: NDArray | number[], value: number): MaskedArray {
  const arr = Array.isArray(x) ? NDArray.fromArray(x) : x;
  const mask = NDArray.empty(arr.shape, DType.Bool);

  for (let i = 0; i < arr.size; i++) {
    mask.setFlat(i, arr.getFlat(i) === value ? 1 : 0);
  }

  return new MaskedArray(arr, mask);
}

/**
 * Mask elements not equal to a value.
 */
export function masked_not_equal(x: NDArray | number[], value: number): MaskedArray {
  const arr = Array.isArray(x) ? NDArray.fromArray(x) : x;
  const mask = NDArray.empty(arr.shape, DType.Bool);

  for (let i = 0; i < arr.size; i++) {
    mask.setFlat(i, arr.getFlat(i) !== value ? 1 : 0);
  }

  return new MaskedArray(arr, mask);
}

/**
 * Mask elements greater than a value.
 */
export function masked_greater(x: NDArray | number[], value: number): MaskedArray {
  const arr = Array.isArray(x) ? NDArray.fromArray(x) : x;
  const mask = NDArray.empty(arr.shape, DType.Bool);

  for (let i = 0; i < arr.size; i++) {
    mask.setFlat(i, arr.getFlat(i) > value ? 1 : 0);
  }

  return new MaskedArray(arr, mask);
}

/**
 * Mask elements greater than or equal to a value.
 */
export function masked_greater_equal(x: NDArray | number[], value: number): MaskedArray {
  const arr = Array.isArray(x) ? NDArray.fromArray(x) : x;
  const mask = NDArray.empty(arr.shape, DType.Bool);

  for (let i = 0; i < arr.size; i++) {
    mask.setFlat(i, arr.getFlat(i) >= value ? 1 : 0);
  }

  return new MaskedArray(arr, mask);
}

/**
 * Mask elements less than a value.
 */
export function masked_less(x: NDArray | number[], value: number): MaskedArray {
  const arr = Array.isArray(x) ? NDArray.fromArray(x) : x;
  const mask = NDArray.empty(arr.shape, DType.Bool);

  for (let i = 0; i < arr.size; i++) {
    mask.setFlat(i, arr.getFlat(i) < value ? 1 : 0);
  }

  return new MaskedArray(arr, mask);
}

/**
 * Mask elements less than or equal to a value.
 */
export function masked_less_equal(x: NDArray | number[], value: number): MaskedArray {
  const arr = Array.isArray(x) ? NDArray.fromArray(x) : x;
  const mask = NDArray.empty(arr.shape, DType.Bool);

  for (let i = 0; i < arr.size; i++) {
    mask.setFlat(i, arr.getFlat(i) <= value ? 1 : 0);
  }

  return new MaskedArray(arr, mask);
}

/**
 * Mask elements inside an interval.
 */
export function masked_inside(
  x: NDArray | number[],
  v1: number,
  v2: number
): MaskedArray {
  const arr = Array.isArray(x) ? NDArray.fromArray(x) : x;
  const mask = NDArray.empty(arr.shape, DType.Bool);
  const lo = Math.min(v1, v2);
  const hi = Math.max(v1, v2);

  for (let i = 0; i < arr.size; i++) {
    const val = arr.getFlat(i);
    mask.setFlat(i, val >= lo && val <= hi ? 1 : 0);
  }

  return new MaskedArray(arr, mask);
}

/**
 * Mask elements outside an interval.
 */
export function masked_outside(
  x: NDArray | number[],
  v1: number,
  v2: number
): MaskedArray {
  const arr = Array.isArray(x) ? NDArray.fromArray(x) : x;
  const mask = NDArray.empty(arr.shape, DType.Bool);
  const lo = Math.min(v1, v2);
  const hi = Math.max(v1, v2);

  for (let i = 0; i < arr.size; i++) {
    const val = arr.getFlat(i);
    mask.setFlat(i, val < lo || val > hi ? 1 : 0);
  }

  return new MaskedArray(arr, mask);
}

/**
 * Mask elements where condition is true.
 */
export function masked_where(
  condition: NDArray | boolean[],
  x: NDArray | number[]
): MaskedArray {
  const arr = Array.isArray(x) ? NDArray.fromArray(x) : x;
  const cond = Array.isArray(condition)
    ? NDArray.fromArray(condition.map(c => c ? 1 : 0), DType.Bool)
    : condition;

  return new MaskedArray(arr, cond);
}

/**
 * Mask NaN and Inf values.
 */
export function masked_invalid(x: NDArray | number[]): MaskedArray {
  const arr = Array.isArray(x) ? NDArray.fromArray(x) : x;
  const mask = NDArray.empty(arr.shape, DType.Bool);

  for (let i = 0; i < arr.size; i++) {
    const val = arr.getFlat(i);
    mask.setFlat(i, !Number.isFinite(val) ? 1 : 0);
  }

  return new MaskedArray(arr, mask);
}

/**
 * Mask values close to a given value.
 */
export function masked_values(
  x: NDArray | number[],
  value: number,
  rtol: number = 1e-5,
  atol: number = 1e-8
): MaskedArray {
  const arr = Array.isArray(x) ? NDArray.fromArray(x) : x;
  const mask = NDArray.empty(arr.shape, DType.Bool);

  for (let i = 0; i < arr.size; i++) {
    const val = arr.getFlat(i);
    const close = Math.abs(val - value) <= atol + rtol * Math.abs(value);
    mask.setFlat(i, close ? 1 : 0);
  }

  return new MaskedArray(arr, mask);
}

/**
 * Create masked array of zeros.
 */
export function zeros(
  shape: number | number[],
  dtype: DType = DType.Float64
): MaskedArray {
  const shapeArr = typeof shape === 'number' ? [shape] : shape;
  return new MaskedArray(NDArray.zeros(shapeArr, dtype));
}

/**
 * Create masked array of ones.
 */
export function ones(
  shape: number | number[],
  dtype: DType = DType.Float64
): MaskedArray {
  const shapeArr = typeof shape === 'number' ? [shape] : shape;
  return new MaskedArray(NDArray.ones(shapeArr, dtype));
}

/**
 * Create empty masked array.
 */
export function empty(
  shape: number | number[],
  dtype: DType = DType.Float64
): MaskedArray {
  const shapeArr = typeof shape === 'number' ? [shape] : shape;
  return new MaskedArray(NDArray.empty(shapeArr, dtype));
}

/**
 * Create fully masked array.
 */
export function masked_all(
  shape: number | number[],
  dtype: DType = DType.Float64
): MaskedArray {
  const shapeArr = typeof shape === 'number' ? [shape] : shape;
  return new MaskedArray(
    NDArray.zeros(shapeArr, dtype),
    true // All masked
  );
}

/**
 * Create fully masked array like another.
 */
export function masked_all_like(prototype: MaskedArray | NDArray): MaskedArray {
  return new MaskedArray(
    NDArray.zeros(prototype.shape, prototype.dtype),
    true
  );
}
```

---

### Module Index

**File:** `src/ts/ma/index.ts`

```typescript
/**
 * NumJS Masked Arrays Module
 *
 * Provides masked array functionality compatible with NumPy's numpy.ma module.
 * Masked arrays support missing or invalid data by maintaining a boolean mask
 * alongside the data array.
 */

// Types and constants
export {
  nomask,
  masked,
  MaskType,
  MaskedArrayError,
  MaskError,
  defaultFillValues,
  maximumFillValue,
  minimumFillValue,
  isMaskedConstant,
} from './types.js';

// Core class
export { MaskedArray, mask_or } from './core.js';

// Creation functions
export {
  masked_array,
  array,
  masked_equal,
  masked_not_equal,
  masked_greater,
  masked_greater_equal,
  masked_less,
  masked_less_equal,
  masked_inside,
  masked_outside,
  masked_where,
  masked_invalid,
  masked_values,
  zeros,
  ones,
  empty,
  masked_all,
  masked_all_like,
} from './creation.js';

// Mask operations
export {
  make_mask,
  make_mask_none,
  getmask,
  getmaskarray,
  getdata,
  is_mask,
  is_masked,
  flatten_mask,
} from './mask_ops.js';

// Fill value operations
export {
  common_fill_value,
  default_fill_value,
  set_fill_value,
} from './fill_value.js';

// Extras (optional advanced functions)
export {
  average,
  median,
  cov,
  corrcoef,
  notmasked_edges,
  notmasked_contiguous,
  clump_masked,
  clump_unmasked,
} from './extras.js';
```

---

## Implementation Order

```
Weeks 7-9: numpy.ma Module

Week 7: Core MaskedArray
├── Day 1: types.ts
│   ├── nomask, masked constants
│   ├── MaskType definition
│   ├── Error classes
│   └── Fill value utilities
│
├── Day 2-3: core.ts (MaskedArray class)
│   ├── Constructor, properties
│   ├── Element access
│   ├── Mask operations (filled, compressed, etc.)
│   └── harden_mask, soften_mask
│
├── Day 4: Arithmetic operations
│   ├── add, subtract, multiply, divide
│   ├── Binary operation with mask propagation
│   └── Comparison operators
│
└── Day 5: Math functions
    ├── Trigonometric
    ├── Exponential/log
    └── Rounding

Week 8: Reductions & Shape
├── Day 1-2: Reductions
│   ├── sum, prod, mean, std, var
│   ├── min, max, ptp
│   ├── argmin, argmax
│   └── cumsum, cumprod, all, any
│
├── Day 3: Shape manipulation
│   ├── reshape, ravel, flatten
│   ├── transpose, swapaxes
│   └── squeeze, expand_dims
│
├── Day 4: Creation functions
│   ├── masked_array, array
│   ├── masked_equal/greater/less/etc.
│   ├── masked_where, masked_invalid
│   └── zeros, ones, empty
│
└── Day 5: Mask operations
    ├── make_mask, getmask, getdata
    ├── is_mask, is_masked
    └── mask_or

Week 9: Extras & Testing
├── Day 1-2: Extras module
│   ├── average, median
│   ├── cov, corrcoef
│   └── Set operations
│
├── Day 3: Edge detection
│   ├── notmasked_edges
│   ├── notmasked_contiguous
│   └── clump_masked/unmasked
│
├── Day 4: mrecords (optional)
│   └── Masked record arrays
│
└── Day 5: Integration tests
    └── NumPy compatibility tests
```

---

## Verification Plan

```bash
# Build
npm run build

# Run tests
npm test -- --grep "ma"

# Test cases:

# MaskedArray creation
✓ masked_array([1,2,3], [F,T,F]) masks element 2
✓ masked_invalid([1,NaN,3]) masks NaN
✓ masked_where([T,F,T], [1,2,3]) masks 1 and 3

# Reductions
✓ ma.sum() excludes masked values
✓ ma.mean() computes mean of non-masked only
✓ ma.count() returns number of non-masked

# filled/compressed
✓ ma.filled(0) replaces masked with 0
✓ ma.compressed() returns 1D non-masked values

# Operations
✓ ma1 + ma2 propagates masks (OR)
✓ ma / 0 masks division by zero
✓ ma.log() masks negative inputs

# Hard mask
✓ Hard mask prevents unmasking
✓ Soft mask allows unmasking
```

---

## Estimated Lines of Code

| File | Lines |
|------|-------|
| types.ts | ~100 |
| core.ts | ~800 |
| creation.ts | ~250 |
| mask_ops.ts | ~150 |
| fill_value.ts | ~80 |
| extras.ts | ~400 |
| index.ts | ~60 |
| **Total** | **~1,840** |

Note: This is a simplified implementation. A full NumPy-compatible implementation would be ~3,500+ lines.
