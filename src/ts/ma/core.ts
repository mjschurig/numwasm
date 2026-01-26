/**
 * NumJS Masked Arrays - Core MaskedArray Class
 *
 * The main MaskedArray class providing arrays with masked (invalid) elements.
 * Compatible with NumPy's numpy.ma.MaskedArray.
 */

import { NDArray } from '../NDArray.js';
import { DType } from '../types.js';
import {
  nomask,
  masked,
  MaskType,
  MaskedConstant,
  MaskedArrayError,
  MaskError,
  isMaskedConstant,
  getDefaultFillValue,
} from './types.js';

/**
 * Array with masked (invalid) elements.
 *
 * A MaskedArray is a combination of a standard NDArray with a boolean
 * mask that indicates which elements should be excluded from operations.
 * Elements where mask is True are considered invalid/masked.
 *
 * @example
 * ```typescript
 * // Create using async factory method
 * const ma = await MaskedArray.create([1, 2, 3, 4, 5], [false, true, false, false, true]);
 *
 * // Or from existing NDArray (sync)
 * const data = await NDArray.fromArray([1, 2, 3, 4, 5]);
 * const mask = await NDArray.fromArray([0, 1, 0, 0, 1], undefined, { dtype: DType.Bool });
 * const ma = new MaskedArray(data, mask);
 *
 * // Operations exclude masked elements
 * ma.mean();  // (1 + 3 + 4) / 3 = 2.667
 * ma.filled(0);  // [1, 0, 3, 4, 0]
 * ma.compressed();  // [1, 3, 4]
 * ```
 */
export class MaskedArray {
  /** Underlying data array */
  public _data: NDArray;

  /** Boolean mask (true = masked/invalid) */
  public _mask: MaskType;

  /** Value to use for masked elements in filled() */
  public _fill_value: number | string | boolean;

  /** If true, mask cannot be changed (unmasked) */
  protected _hardmask: boolean;

  /** If true, mask is shared with parent view */
  protected _sharedmask: boolean;

  /**
   * Create a masked array from NDArray data and optional mask.
   * For creating from raw arrays, use the static `create()` method.
   *
   * @param data - Input NDArray or MaskedArray
   * @param mask - Boolean mask NDArray or nomask
   * @param fill_value - Fill value for masked elements
   * @param hard_mask - If true, mask cannot be unset
   */
  constructor(
    data: NDArray | MaskedArray,
    mask: MaskType = nomask,
    fill_value: number | string | boolean | null = null,
    hard_mask: boolean = false
  ) {
    if (data instanceof MaskedArray) {
      this._data = data._data;
      // Combine masks if both exist
      if (mask !== nomask && data._mask !== nomask) {
        this._mask = this._combineMasksSync(data._mask as NDArray, mask as NDArray);
      } else if (mask !== nomask) {
        this._mask = mask;
      } else {
        this._mask = data._mask;
      }
      this._fill_value = fill_value ?? data._fill_value;
      this._hardmask = hard_mask || data._hardmask;
    } else {
      this._data = data;
      this._mask = mask;
      this._fill_value = fill_value ?? getDefaultFillValue(data.dtype);
      this._hardmask = hard_mask;
    }
    this._sharedmask = false;
  }

  /**
   * Create a masked array from raw data (async factory method).
   *
   * @param data - Input array data
   * @param mask - Boolean mask array or nomask
   * @param dtype - Data type
   * @param fill_value - Fill value for masked elements
   * @param hard_mask - If true, mask cannot be unset
   * @param shrink - If true, shrink mask to nomask if all false
   * @returns Promise resolving to MaskedArray
   */
  static async create(
    data: number[] | NDArray | MaskedArray,
    mask: boolean[] | MaskType = nomask,
    dtype: DType = DType.Float64,
    fill_value: number | string | boolean | null = null,
    hard_mask: boolean = false,
    shrink: boolean = true
  ): Promise<MaskedArray> {
    // Handle MaskedArray input
    if (data instanceof MaskedArray) {
      let finalMask = data._mask;
      if (mask !== nomask) {
        const newMask = await MaskedArray._processMaskAsync(mask, data.shape, shrink);
        if (newMask !== nomask && finalMask !== nomask) {
          finalMask = await MaskedArray._combineMasksAsync(finalMask as NDArray, newMask as NDArray);
        } else if (newMask !== nomask) {
          finalMask = newMask;
        }
      }
      return new MaskedArray(data._data, finalMask, fill_value ?? data._fill_value, hard_mask || data._hardmask);
    }

    // Convert data to NDArray
    let ndData: NDArray;
    if (data instanceof NDArray) {
      ndData = data;
    } else {
      ndData = await NDArray.fromArray(data, undefined, { dtype });
    }

    // Process mask
    const finalMask = await MaskedArray._processMaskAsync(mask, ndData.shape, shrink);

    return new MaskedArray(ndData, finalMask, fill_value, hard_mask);
  }

  /**
   * Process mask input into standard form (async).
   */
  private static async _processMaskAsync(
    mask: MaskType | boolean[],
    shape: number[],
    shrink: boolean
  ): Promise<MaskType> {
    if (mask === nomask) {
      return nomask;
    }

    if (typeof mask === 'boolean') {
      if (!mask && shrink) {
        return nomask;
      }
      if (mask) {
        const fullMask = await NDArray.full(shape, 1, { dtype: DType.Bool });
        return fullMask;
      }
      return nomask;
    }

    let maskArr: NDArray;
    if (Array.isArray(mask)) {
      const boolData = mask.map((m) => (m ? 1 : 0));
      maskArr = await NDArray.fromArray(boolData, undefined, { dtype: DType.Bool });
    } else {
      maskArr = mask.dtype === DType.Bool ? mask : mask.astype(DType.Bool);
    }

    // Broadcast mask to shape if needed
    if (!MaskedArray._shapesEqual(maskArr.shape, shape)) {
      maskArr = await MaskedArray._broadcastMaskAsync(maskArr, shape);
    }

    // Shrink to nomask if all false
    if (shrink && !MaskedArray._anyTrue(maskArr)) {
      return nomask;
    }

    return maskArr;
  }

  /**
   * Broadcast mask to target shape (async).
   */
  private static async _broadcastMaskAsync(mask: NDArray, targetShape: number[]): Promise<NDArray> {
    const result = await NDArray.zeros(targetShape, { dtype: DType.Bool });
    const maskShape = mask.shape;
    const maskNdim = maskShape.length;
    const targetNdim = targetShape.length;

    const paddedMaskShape = new Array(targetNdim).fill(1);
    for (let i = 0; i < maskNdim; i++) {
      paddedMaskShape[targetNdim - maskNdim + i] = maskShape[i];
    }

    for (let i = 0; i < targetNdim; i++) {
      if (paddedMaskShape[i] !== 1 && paddedMaskShape[i] !== targetShape[i]) {
        throw new MaskError(`Cannot broadcast mask shape ${maskShape} to ${targetShape}`);
      }
    }

    for (let i = 0; i < result.size; i++) {
      const indices = MaskedArray._flatToMultiIndex(i, targetShape);
      const maskIndices = new Array(maskNdim);
      for (let j = 0; j < maskNdim; j++) {
        const targetIdx = targetNdim - maskNdim + j;
        maskIndices[j] = maskShape[j] === 1 ? 0 : indices[targetIdx];
      }
      const maskFlatIdx = MaskedArray._multiToFlatIndex(maskIndices, maskShape);
      result.setFlat(i, mask.getFlat(maskFlatIdx));
    }

    return result;
  }

  /**
   * Combine two masks with OR (async).
   */
  private static async _combineMasksAsync(m1: NDArray, m2: NDArray): Promise<NDArray> {
    if (!MaskedArray._shapesEqual(m1.shape, m2.shape)) {
      throw new MaskError(`Mask shapes must match: ${m1.shape} vs ${m2.shape}`);
    }
    const result = await NDArray.zeros(m1.shape, { dtype: DType.Bool });
    for (let i = 0; i < result.size; i++) {
      result.setFlat(i, m1.getFlat(i) || m2.getFlat(i) ? 1 : 0);
    }
    return result;
  }

  /**
   * Combine two masks with OR (sync, uses pre-existing NDArrays).
   */
  private _combineMasksSync(m1: NDArray, m2: NDArray): NDArray {
    if (!MaskedArray._shapesEqual(m1.shape, m2.shape)) {
      throw new MaskError(`Mask shapes must match: ${m1.shape} vs ${m2.shape}`);
    }
    // Use m1.copy() and modify in place
    const result = m1.copy();
    for (let i = 0; i < result.size; i++) {
      result.setFlat(i, m1.getFlat(i) || m2.getFlat(i) ? 1 : 0);
    }
    return result;
  }

  /* ============ Properties ============ */

  /** The underlying data array. */
  get data(): NDArray {
    return this._data;
  }

  /** The mask array (or nomask if no masking). */
  get mask(): MaskType {
    return this._mask;
  }

  /** Fill value for masked elements. */
  get fill_value(): number | string | boolean {
    return this._fill_value;
  }

  set fill_value(value: number | string | boolean) {
    this._fill_value = value;
  }

  /** Shape of the array. */
  get shape(): number[] {
    return this._data.shape;
  }

  /** Number of dimensions. */
  get ndim(): number {
    return this._data.ndim;
  }

  /** Total number of elements. */
  get size(): number {
    return this._data.size;
  }

  /** Data type. */
  get dtype(): DType {
    return this._data.dtype;
  }

  /** Transpose. */
  get T(): MaskedArray {
    return this.transpose();
  }

  /** Whether the mask is hard (cannot be unset). */
  get hardmask(): boolean {
    return this._hardmask;
  }

  /* ============ Element Access ============ */

  /**
   * Get element at flat index.
   * @returns Value or masked constant if element is masked
   */
  getFlat(index: number): number | MaskedConstant {
    if (this._mask !== nomask && (this._mask as NDArray).getFlat(index)) {
      return masked;
    }
    return this._data.getFlat(index);
  }

  /**
   * Set element at flat index (async if creating new mask).
   * @param index - Flat index
   * @param value - Value to set (or masked to mask the element)
   */
  async setFlatAsync(index: number, value: number | MaskedConstant): Promise<void> {
    if (value === masked || isMaskedConstant(value)) {
      if (this._mask === nomask) {
        this._mask = await NDArray.zeros(this._data.shape, { dtype: DType.Bool });
      }
      (this._mask as NDArray).setFlat(index, 1);
    } else {
      this._data.setFlat(index, value as number);
      if (this._mask !== nomask && !this._hardmask) {
        (this._mask as NDArray).setFlat(index, 0);
      }
    }
  }

  /**
   * Set element at flat index (sync, assumes mask exists if needed).
   */
  setFlat(index: number, value: number | MaskedConstant): void {
    if (value === masked || isMaskedConstant(value)) {
      if (this._mask === nomask) {
        throw new MaskedArrayError('Cannot set masked value when mask is nomask. Use setFlatAsync instead.');
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
   */
  count(axis: number | null = null): number {
    if (axis !== null) {
      throw new MaskedArrayError('count with axis not yet implemented');
    }

    if (this._mask === nomask) {
      return this._data.size;
    }

    let count = 0;
    const mask = this._mask as NDArray;
    for (let i = 0; i < mask.size; i++) {
      if (!mask.getFlat(i)) {
        count++;
      }
    }
    return count;
  }

  /**
   * Return data with masked values replaced by fill_value.
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
   * Return 1D array of non-masked values (async).
   */
  async compressed(): Promise<NDArray> {
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

    return await NDArray.fromArray(values, undefined, { dtype: this._data.dtype });
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
      if (!MaskedArray._anyTrue(this._mask as NDArray)) {
        this._mask = nomask;
      }
    }
    return this;
  }

  /* ============ Arithmetic Operations ============ */

  /**
   * Add another array or scalar (async).
   */
  async add(other: MaskedArray | NDArray | number): Promise<MaskedArray> {
    return this._binaryOpAsync(other, (a, b) => a + b);
  }

  /**
   * Subtract another array or scalar (async).
   */
  async subtract(other: MaskedArray | NDArray | number): Promise<MaskedArray> {
    return this._binaryOpAsync(other, (a, b) => a - b);
  }

  /**
   * Multiply by another array or scalar (async).
   */
  async multiply(other: MaskedArray | NDArray | number): Promise<MaskedArray> {
    return this._binaryOpAsync(other, (a, b) => a * b);
  }

  /**
   * Divide by another array or scalar (async).
   */
  async divide(other: MaskedArray | NDArray | number): Promise<MaskedArray> {
    return this._binaryOpAsync(other, (a, b) => a / b, (_a, b) => b === 0);
  }

  /**
   * True divide (async).
   */
  async true_divide(other: MaskedArray | NDArray | number): Promise<MaskedArray> {
    return this.divide(other);
  }

  /**
   * Floor divide (async).
   */
  async floor_divide(other: MaskedArray | NDArray | number): Promise<MaskedArray> {
    return this._binaryOpAsync(other, (a, b) => Math.floor(a / b), (_a, b) => b === 0);
  }

  /**
   * Power (async).
   */
  async power(other: MaskedArray | NDArray | number): Promise<MaskedArray> {
    return this._binaryOpAsync(other, (a, b) => Math.pow(a, b));
  }

  /**
   * Modulo (async).
   */
  async mod(other: MaskedArray | NDArray | number): Promise<MaskedArray> {
    return this._binaryOpAsync(other, (a, b) => a % b, (_a, b) => b === 0);
  }

  /**
   * Negate (async).
   */
  async neg(): Promise<MaskedArray> {
    return this._unaryOpAsync((x) => -x);
  }

  /**
   * Absolute value (async).
   */
  async abs(): Promise<MaskedArray> {
    return this._unaryOpAsync(Math.abs);
  }

  /* ============ Comparison Operations ============ */

  async equal(other: MaskedArray | NDArray | number): Promise<MaskedArray> {
    return this._binaryOpAsync(other, (a, b) => (a === b ? 1 : 0));
  }

  async not_equal(other: MaskedArray | NDArray | number): Promise<MaskedArray> {
    return this._binaryOpAsync(other, (a, b) => (a !== b ? 1 : 0));
  }

  async less(other: MaskedArray | NDArray | number): Promise<MaskedArray> {
    return this._binaryOpAsync(other, (a, b) => (a < b ? 1 : 0));
  }

  async less_equal(other: MaskedArray | NDArray | number): Promise<MaskedArray> {
    return this._binaryOpAsync(other, (a, b) => (a <= b ? 1 : 0));
  }

  async greater(other: MaskedArray | NDArray | number): Promise<MaskedArray> {
    return this._binaryOpAsync(other, (a, b) => (a > b ? 1 : 0));
  }

  async greater_equal(other: MaskedArray | NDArray | number): Promise<MaskedArray> {
    return this._binaryOpAsync(other, (a, b) => (a >= b ? 1 : 0));
  }

  /* ============ Math Functions ============ */

  async sin(): Promise<MaskedArray> {
    return this._unaryOpAsync(Math.sin);
  }
  async cos(): Promise<MaskedArray> {
    return this._unaryOpAsync(Math.cos);
  }
  async tan(): Promise<MaskedArray> {
    return this._unaryOpAsync(Math.tan);
  }
  async arcsin(): Promise<MaskedArray> {
    return this._unaryOpAsync(Math.asin, (x) => x < -1 || x > 1);
  }
  async arccos(): Promise<MaskedArray> {
    return this._unaryOpAsync(Math.acos, (x) => x < -1 || x > 1);
  }
  async arctan(): Promise<MaskedArray> {
    return this._unaryOpAsync(Math.atan);
  }

  async sinh(): Promise<MaskedArray> {
    return this._unaryOpAsync(Math.sinh);
  }
  async cosh(): Promise<MaskedArray> {
    return this._unaryOpAsync(Math.cosh);
  }
  async tanh(): Promise<MaskedArray> {
    return this._unaryOpAsync(Math.tanh);
  }
  async arcsinh(): Promise<MaskedArray> {
    return this._unaryOpAsync(Math.asinh);
  }
  async arccosh(): Promise<MaskedArray> {
    return this._unaryOpAsync(Math.acosh, (x) => x < 1);
  }
  async arctanh(): Promise<MaskedArray> {
    return this._unaryOpAsync(Math.atanh, (x) => x <= -1 || x >= 1);
  }

  async exp(): Promise<MaskedArray> {
    return this._unaryOpAsync(Math.exp);
  }
  async log(): Promise<MaskedArray> {
    return this._unaryOpAsync(Math.log, (x) => x <= 0);
  }
  async log10(): Promise<MaskedArray> {
    return this._unaryOpAsync(Math.log10, (x) => x <= 0);
  }
  async log2(): Promise<MaskedArray> {
    return this._unaryOpAsync(Math.log2, (x) => x <= 0);
  }

  async sqrt(): Promise<MaskedArray> {
    return this._unaryOpAsync(Math.sqrt, (x) => x < 0);
  }
  async square(): Promise<MaskedArray> {
    return this._unaryOpAsync((x) => x * x);
  }

  async floor(): Promise<MaskedArray> {
    return this._unaryOpAsync(Math.floor);
  }
  async ceil(): Promise<MaskedArray> {
    return this._unaryOpAsync(Math.ceil);
  }
  async round(): Promise<MaskedArray> {
    return this._unaryOpAsync(Math.round);
  }
  async trunc(): Promise<MaskedArray> {
    return this._unaryOpAsync(Math.trunc);
  }

  /* ============ Reductions ============ */

  /**
   * Sum of non-masked elements.
   */
  sum(axis: number | null = null): number {
    if (axis !== null) {
      throw new MaskedArrayError('sum with axis not yet implemented');
    }

    let result = 0;

    for (let i = 0; i < this._data.size; i++) {
      if (this._mask === nomask || !(this._mask as NDArray).getFlat(i)) {
        result += this._data.getFlat(i);
      }
    }

    return result;
  }

  /**
   * Product of non-masked elements.
   */
  prod(axis: number | null = null): number {
    if (axis !== null) {
      throw new MaskedArrayError('prod with axis not yet implemented');
    }

    let result = 1;

    for (let i = 0; i < this._data.size; i++) {
      if (this._mask === nomask || !(this._mask as NDArray).getFlat(i)) {
        result *= this._data.getFlat(i);
      }
    }

    return result;
  }

  /**
   * Mean of non-masked elements.
   */
  mean(axis: number | null = null): number {
    if (axis !== null) {
      throw new MaskedArrayError('mean with axis not yet implemented');
    }

    const s = this.sum();
    const c = this.count();
    return c === 0 ? NaN : s / c;
  }

  /**
   * Variance of non-masked elements.
   */
  var(axis: number | null = null, ddof: number = 0): number {
    if (axis !== null) {
      throw new MaskedArrayError('var with axis not yet implemented');
    }

    const m = this.mean();
    const c = this.count();
    if (c - ddof <= 0) return NaN;

    let sumSq = 0;
    for (let i = 0; i < this._data.size; i++) {
      if (this._mask === nomask || !(this._mask as NDArray).getFlat(i)) {
        const diff = this._data.getFlat(i) - m;
        sumSq += diff * diff;
      }
    }

    return sumSq / (c - ddof);
  }

  /**
   * Standard deviation of non-masked elements.
   */
  std(axis: number | null = null, ddof: number = 0): number {
    return Math.sqrt(this.var(axis, ddof));
  }

  /**
   * Minimum of non-masked elements.
   */
  min(axis: number | null = null): number {
    if (axis !== null) {
      throw new MaskedArrayError('min with axis not yet implemented');
    }

    let result = Infinity;
    for (let i = 0; i < this._data.size; i++) {
      if (this._mask === nomask || !(this._mask as NDArray).getFlat(i)) {
        const v = this._data.getFlat(i);
        if (v < result) result = v;
      }
    }

    return result;
  }

  /**
   * Maximum of non-masked elements.
   */
  max(axis: number | null = null): number {
    if (axis !== null) {
      throw new MaskedArrayError('max with axis not yet implemented');
    }

    let result = -Infinity;
    for (let i = 0; i < this._data.size; i++) {
      if (this._mask === nomask || !(this._mask as NDArray).getFlat(i)) {
        const v = this._data.getFlat(i);
        if (v > result) result = v;
      }
    }

    return result;
  }

  /**
   * Peak-to-peak (max - min) of non-masked elements.
   */
  ptp(axis: number | null = null): number {
    return this.max(axis) - this.min(axis);
  }

  /**
   * Index of minimum non-masked element.
   */
  argmin(axis: number | null = null): number {
    if (axis !== null) {
      throw new MaskedArrayError('argmin with axis not yet implemented');
    }

    let bestIdx = -1;
    let bestVal = Infinity;

    for (let i = 0; i < this._data.size; i++) {
      if (this._mask !== nomask && (this._mask as NDArray).getFlat(i)) {
        continue;
      }
      const val = this._data.getFlat(i);
      if (val < bestVal) {
        bestVal = val;
        bestIdx = i;
      }
    }

    return bestIdx;
  }

  /**
   * Index of maximum non-masked element.
   */
  argmax(axis: number | null = null): number {
    if (axis !== null) {
      throw new MaskedArrayError('argmax with axis not yet implemented');
    }

    let bestIdx = -1;
    let bestVal = -Infinity;

    for (let i = 0; i < this._data.size; i++) {
      if (this._mask !== nomask && (this._mask as NDArray).getFlat(i)) {
        continue;
      }
      const val = this._data.getFlat(i);
      if (val > bestVal) {
        bestVal = val;
        bestIdx = i;
      }
    }

    return bestIdx;
  }

  /**
   * Test if all non-masked elements are true.
   */
  all(axis: number | null = null): boolean {
    if (axis !== null) {
      throw new MaskedArrayError('all with axis not yet implemented');
    }

    for (let i = 0; i < this._data.size; i++) {
      if (this._mask !== nomask && (this._mask as NDArray).getFlat(i)) {
        continue;
      }
      if (!this._data.getFlat(i)) {
        return false;
      }
    }

    return true;
  }

  /**
   * Test if any non-masked element is true.
   */
  any(axis: number | null = null): boolean {
    if (axis !== null) {
      throw new MaskedArrayError('any with axis not yet implemented');
    }

    for (let i = 0; i < this._data.size; i++) {
      if (this._mask !== nomask && (this._mask as NDArray).getFlat(i)) {
        continue;
      }
      if (this._data.getFlat(i)) {
        return true;
      }
    }

    return false;
  }

  /* ============ Shape Manipulation ============ */

  /**
   * Reshape the array and mask.
   */
  reshape(newshape: number[]): MaskedArray {
    const newData = this._data.reshape(newshape);
    const newMask =
      this._mask === nomask ? nomask : (this._mask as NDArray).reshape(newshape);

    return new MaskedArray(newData, newMask, this._fill_value, this._hardmask);
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
      this._fill_value,
      this._hardmask
    );
  }

  /**
   * Transpose the array.
   */
  transpose(axes?: number[]): MaskedArray {
    const newData = this._data.transpose(axes);
    const newMask =
      this._mask === nomask ? nomask : (this._mask as NDArray).transpose(axes);

    return new MaskedArray(newData, newMask, this._fill_value, this._hardmask);
  }

  /**
   * Swap two axes.
   */
  swapaxes(axis1: number, axis2: number): MaskedArray {
    const newData = this._data.swapaxes(axis1, axis2);
    const newMask =
      this._mask === nomask
        ? nomask
        : (this._mask as NDArray).swapaxes(axis1, axis2);

    return new MaskedArray(newData, newMask, this._fill_value, this._hardmask);
  }

  /**
   * Remove axes of length 1.
   */
  squeeze(axis?: number): MaskedArray {
    const newData = this._data.squeeze(axis);
    const newMask =
      this._mask === nomask ? nomask : (this._mask as NDArray).squeeze(axis);

    return new MaskedArray(newData, newMask, this._fill_value, this._hardmask);
  }

  /**
   * Add axis of length 1.
   */
  expand_dims(axis: number): MaskedArray {
    const newData = this._data.expandDims(axis);
    const newMask =
      this._mask === nomask ? nomask : (this._mask as NDArray).expandDims(axis);

    return new MaskedArray(newData, newMask, this._fill_value, this._hardmask);
  }

  /**
   * Create a copy.
   */
  copy(): MaskedArray {
    return new MaskedArray(
      this._data.copy(),
      this._mask === nomask ? nomask : (this._mask as NDArray).copy(),
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
      this._fill_value,
      this._hardmask
    );
  }

  /**
   * Convert to regular array (unmasked data).
   */
  toArray(): (number | null)[] {
    const result: (number | null)[] = [];
    for (let i = 0; i < this._data.size; i++) {
      if (this._mask !== nomask && (this._mask as NDArray).getFlat(i)) {
        result.push(null);
      } else {
        result.push(this._data.getFlat(i));
      }
    }
    return result;
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
   * Apply binary operation with mask propagation (async).
   */
  protected async _binaryOpAsync(
    other: MaskedArray | NDArray | number,
    op: (a: number, b: number) => number,
    extraMask?: (a: number, b: number) => boolean
  ): Promise<MaskedArray> {
    let otherData: NDArray;
    let otherMask: MaskType;

    if (typeof other === 'number') {
      otherData = await NDArray.full(this._data.shape, other, { dtype: this._data.dtype });
      otherMask = nomask;
    } else if (other instanceof MaskedArray) {
      otherData = other._data;
      otherMask = other._mask;
    } else {
      otherData = other;
      otherMask = nomask;
    }

    // Check shapes match
    if (!MaskedArray._shapesEqual(this._data.shape, otherData.shape)) {
      throw new MaskedArrayError(
        `Shape mismatch: ${this._data.shape} vs ${otherData.shape}`
      );
    }

    // Compute result
    const resultData = await NDArray.empty(this._data.shape, { dtype: this._data.dtype });
    let extraMaskArr: NDArray | null = null;

    for (let i = 0; i < resultData.size; i++) {
      const a = this._data.getFlat(i);
      const b = otherData.getFlat(i);
      resultData.setFlat(i, op(a, b));

      if (extraMask && extraMask(a, b)) {
        if (extraMaskArr === null) {
          extraMaskArr = await NDArray.zeros(this._data.shape, { dtype: DType.Bool });
        }
        extraMaskArr.setFlat(i, 1);
      }
    }

    // Combine masks
    let resultMask = await this._combineMasksOptional(this._mask, otherMask);
    if (extraMaskArr !== null) {
      resultMask = await this._combineMasksOptional(resultMask, extraMaskArr);
    }

    return new MaskedArray(resultData, resultMask, this._fill_value);
  }

  /**
   * Combine two optional masks.
   */
  private async _combineMasksOptional(m1: MaskType, m2: MaskType): Promise<MaskType> {
    if (m1 === nomask && m2 === nomask) return nomask;
    if (m1 === nomask) return m2;
    if (m2 === nomask) return m1;
    return await MaskedArray._combineMasksAsync(m1 as NDArray, m2 as NDArray);
  }

  /**
   * Apply unary operation (async).
   */
  protected async _unaryOpAsync(
    op: (x: number) => number,
    shouldMask?: (x: number) => boolean
  ): Promise<MaskedArray> {
    const resultData = await NDArray.empty(this._data.shape, { dtype: this._data.dtype });
    let resultMask: MaskType =
      this._mask === nomask ? nomask : (this._mask as NDArray).copy();

    for (let i = 0; i < resultData.size; i++) {
      const x = this._data.getFlat(i);
      resultData.setFlat(i, op(x));

      if (shouldMask && shouldMask(x)) {
        if (resultMask === nomask) {
          resultMask = await NDArray.zeros(this._data.shape, { dtype: DType.Bool });
        }
        (resultMask as NDArray).setFlat(i, 1);
      }
    }

    return new MaskedArray(resultData, resultMask, this._fill_value);
  }

  /* ============ Static Helpers ============ */

  private static _anyTrue(arr: NDArray): boolean {
    for (let i = 0; i < arr.size; i++) {
      if (arr.getFlat(i)) return true;
    }
    return false;
  }

  private static _shapesEqual(s1: number[], s2: number[]): boolean {
    if (s1.length !== s2.length) return false;
    return s1.every((v, i) => v === s2[i]);
  }

  private static _flatToMultiIndex(flatIdx: number, shape: number[]): number[] {
    const indices = new Array(shape.length);
    let remaining = flatIdx;
    for (let i = shape.length - 1; i >= 0; i--) {
      indices[i] = remaining % shape[i];
      remaining = Math.floor(remaining / shape[i]);
    }
    return indices;
  }

  private static _multiToFlatIndex(indices: number[], shape: number[]): number {
    let flatIdx = 0;
    let multiplier = 1;
    for (let i = shape.length - 1; i >= 0; i--) {
      flatIdx += indices[i] * multiplier;
      multiplier *= shape[i];
    }
    return flatIdx;
  }
}
