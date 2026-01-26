/**
 * NumJS Random Module
 *
 * Provides NumPy-compatible random number generation with WebAssembly acceleration.
 * Implements the modern Generator API with PCG64 as the default BitGenerator.
 *
 * @example
 * ```typescript
 * import { random } from 'numjs';
 *
 * // Create a seeded generator for reproducibility
 * const rng = random.default_rng(12345);
 *
 * // Generate random values
 * const uniform = rng.random(10);           // [0, 1) uniform
 * const normal = rng.normal(0, 1, [5, 5]);  // 5x5 standard normal
 * const integers = rng.integers(0, 100, 20); // 20 integers in [0, 100)
 *
 * // Spawn independent generators for parallel work
 * const [rng1, rng2] = rng.spawn(2);
 * ```
 */

import { NDArray } from './NDArray.js';
import { DType } from './types.js';
import { getWasmModule, loadWasmModule } from './wasm-loader.js';
import { take } from './indexing.js';

/* ============ Type Definitions ============ */

/**
 * Size specification for random output arrays.
 */
export type SizeType = number | number[] | null;

/**
 * State object for PCG64 generator.
 */
export interface PCG64State {
  bit_generator: 'PCG64';
  state: {
    state: bigint;
    inc: bigint;
  };
  has_uint32: number;
  uinteger: number;
}

/* ============ BitGenerator Base Class ============ */

/**
 * Abstract base class for all bit generators.
 * BitGenerators provide the core random number generation.
 */
export abstract class BitGenerator {
  protected _wasmState: number = 0;

  /**
   * Generate a random 64-bit unsigned integer.
   */
  abstract next_uint64(): bigint;

  /**
   * Generate a random 32-bit unsigned integer.
   */
  abstract next_uint32(): number;

  /**
   * Generate a random double in [0, 1).
   */
  abstract next_double(): number;

  /**
   * Get the current state as a serializable object.
   */
  abstract getState(): object;

  /**
   * Restore state from a serializable object.
   */
  abstract setState(state: object): void;

  /**
   * Create independent child BitGenerators.
   */
  abstract spawn(n_children: number): BitGenerator[];

  /**
   * Get the WASM state pointer for C calls.
   * @internal
   */
  get wasmStatePtr(): number {
    return this._wasmState;
  }

  /**
   * Clean up WASM resources.
   */
  abstract dispose(): void;
}

/* ============ SeedSequence ============ */

/**
 * SeedSequence for reproducible entropy management.
 *
 * Mixes entropy from multiple sources to generate high-quality seeds
 * for any BitGenerator. Supports spawning for parallel streams.
 *
 * @example
 * ```typescript
 * // Create from integer seed
 * const ss = new SeedSequence(12345);
 *
 * // Create from entropy array
 * const ss2 = new SeedSequence([1, 2, 3, 4]);
 *
 * // Create with OS entropy (no seed)
 * const ss3 = new SeedSequence();
 *
 * // Spawn independent sequences
 * const [child1, child2] = ss.spawn(2);
 * ```
 */
export class SeedSequence {
  private _entropy: Uint32Array;
  private _spawnKey: number[];
  private _poolSize: number;
  private _nChildrenSpawned: number = 0;

  constructor(
    entropy?: number | bigint | number[] | Uint32Array | null,
    spawnKey: number[] = [],
    poolSize: number = 4
  ) {
    this._spawnKey = spawnKey;
    this._poolSize = poolSize;

    if (entropy === null || entropy === undefined) {
      // Use crypto.getRandomValues for OS entropy
      this._entropy = new Uint32Array(poolSize);
      if (typeof crypto !== 'undefined' && crypto.getRandomValues) {
        crypto.getRandomValues(this._entropy);
      } else {
        // Fallback: use Date.now() + Math.random()
        for (let i = 0; i < poolSize; i++) {
          this._entropy[i] = ((Date.now() * Math.random()) >>> 0) ^ (Math.random() * 0xffffffff >>> 0);
        }
      }
    } else if (typeof entropy === 'number') {
      // Single integer seed - split into low and high 32-bit parts
      this._entropy = new Uint32Array([
        entropy >>> 0,
        Math.floor(entropy / 0x100000000) >>> 0
      ]);
    } else if (typeof entropy === 'bigint') {
      // BigInt seed
      this._entropy = new Uint32Array([
        Number(entropy & 0xFFFFFFFFn),
        Number((entropy >> 32n) & 0xFFFFFFFFn),
        Number((entropy >> 64n) & 0xFFFFFFFFn),
        Number((entropy >> 96n) & 0xFFFFFFFFn)
      ]);
    } else if (entropy instanceof Uint32Array) {
      this._entropy = new Uint32Array(entropy);
    } else {
      // Array of numbers
      this._entropy = new Uint32Array(entropy.map(n => n >>> 0));
    }
  }

  /**
   * Generate state words for seeding a BitGenerator.
   *
   * @param nWords - Number of 32-bit words to generate
   * @param dtype - Output type ('uint32' or 'uint64')
   */
  generateState(nWords: number, dtype: 'uint32' | 'uint64' = 'uint32'): Uint32Array | BigUint64Array {
    const wasm = getWasmModule();

    // Allocate WASM memory
    const entropyPtr = wasm._malloc(this._entropy.length * 4);
    const spawnKeyPtr = wasm._malloc(Math.max(this._spawnKey.length, 1) * 4);
    const outPtr = wasm._malloc(nWords * 4);

    try {
      // Copy entropy to WASM
      for (let i = 0; i < this._entropy.length; i++) {
        wasm.HEAPU32[(entropyPtr >> 2) + i] = this._entropy[i];
      }

      // Copy spawn key to WASM
      for (let i = 0; i < this._spawnKey.length; i++) {
        wasm.HEAPU32[(spawnKeyPtr >> 2) + i] = this._spawnKey[i];
      }

      // Call WASM seed_seq_generate_words
      wasm._seed_seq_generate_words(
        entropyPtr, this._entropy.length,
        spawnKeyPtr, this._spawnKey.length,
        this._poolSize,
        outPtr, nWords
      );

      // Copy result
      const result = new Uint32Array(nWords);
      for (let i = 0; i < nWords; i++) {
        result[i] = wasm.HEAPU32[(outPtr >> 2) + i];
      }

      if (dtype === 'uint64') {
        const result64 = new BigUint64Array(Math.floor(nWords / 2));
        for (let i = 0; i < result64.length; i++) {
          result64[i] = BigInt(result[i * 2]) | (BigInt(result[i * 2 + 1]) << 32n);
        }
        return result64;
      }

      return result;
    } finally {
      wasm._free(entropyPtr);
      wasm._free(spawnKeyPtr);
      wasm._free(outPtr);
    }
  }

  /**
   * Spawn n_children independent SeedSequences.
   *
   * Each child sequence is independent and can be used to seed
   * a separate BitGenerator for parallel random streams.
   */
  spawn(nChildren: number): SeedSequence[] {
    const children: SeedSequence[] = [];
    for (let i = 0; i < nChildren; i++) {
      children.push(new SeedSequence(
        this._entropy,
        [...this._spawnKey, this._nChildrenSpawned + i],
        this._poolSize
      ));
    }
    this._nChildrenSpawned += nChildren;
    return children;
  }

  /**
   * Get the entropy array.
   */
  get entropy(): Uint32Array {
    return new Uint32Array(this._entropy);
  }

  /**
   * Get the spawn key.
   */
  get spawnKey(): number[] {
    return [...this._spawnKey];
  }
}

/* ============ PCG64 BitGenerator ============ */

/**
 * PCG64 BitGenerator - the default for NumJS.
 *
 * Uses a 128-bit LCG with XSL-RR output function.
 * Provides a period of 2^128 and excellent statistical properties.
 *
 * @example
 * ```typescript
 * // Create with integer seed
 * const pcg = new PCG64(12345);
 *
 * // Create with SeedSequence
 * const ss = new SeedSequence([1, 2, 3]);
 * const pcg2 = new PCG64(ss);
 *
 * // Generate values
 * const u64 = pcg.next_uint64();
 * const u32 = pcg.next_uint32();
 * const dbl = pcg.next_double();
 * ```
 */
export class PCG64 extends BitGenerator {
  private _seedSequence: SeedSequence | null = null;
  private _disposed: boolean = false;

  constructor(seed?: number | bigint | number[] | SeedSequence | null) {
    super();

    const wasm = getWasmModule();
    this._wasmState = wasm._pcg64_create();

    if (this._wasmState === 0) {
      throw new Error('Failed to create PCG64 state: memory allocation failed');
    }

    if (seed instanceof SeedSequence) {
      this._seedSequence = seed;
    } else {
      this._seedSequence = new SeedSequence(seed);
    }

    // Generate 8 uint32 values for 4 x 32-bit parts of two 128-bit numbers
    const state = this._seedSequence.generateState(8, 'uint32') as Uint32Array;
    this._initFromState(state);
  }

  private _initFromState(state: Uint32Array): void {
    const wasm = getWasmModule();
    // PCG64 needs 128-bit seed and 128-bit increment
    // Each 128-bit value is represented as (high64, low64) where each 64-bit is (high32, low32)
    // state[0..3] = seed: [hh, hl, lh, ll] forms seed_high | seed_low
    // state[4..7] = inc:  [hh, hl, lh, ll] forms inc_high | inc_low
    wasm._pcg64_seed_parts(
      this._wasmState,
      state[0], state[1], state[2], state[3],  // seed_hh, seed_hl, seed_lh, seed_ll
      state[4], state[5], state[6], state[7]   // inc_hh, inc_hl, inc_lh, inc_ll
    );
  }

  private ensureNotDisposed(): void {
    if (this._disposed) {
      throw new Error('PCG64 has been disposed');
    }
  }

  next_uint64(): bigint {
    this.ensureNotDisposed();
    const wasm = getWasmModule();

    // Get both parts of the 64-bit value
    const highPtr = wasm._malloc(4);
    try {
      const low = wasm._pcg64_next64_parts(this._wasmState, highPtr);
      const high = wasm.HEAPU32[highPtr >> 2];
      return BigInt(low >>> 0) | (BigInt(high) << 32n);
    } finally {
      wasm._free(highPtr);
    }
  }

  next_uint32(): number {
    this.ensureNotDisposed();
    const wasm = getWasmModule();
    return wasm._pcg64_next32(this._wasmState) >>> 0;
  }

  next_double(): number {
    this.ensureNotDisposed();
    const wasm = getWasmModule();
    return wasm._pcg64_next_double(this._wasmState);
  }

  /**
   * Advance the state by delta steps.
   *
   * This is equivalent to calling next_uint64() delta times,
   * but runs in O(log(delta)) time.
   */
  advance(delta: bigint): this {
    this.ensureNotDisposed();
    const wasm = getWasmModule();
    const low = Number(delta & 0xFFFFFFFFFFFFFFFFn);
    const high = Number((delta >> 64n) & 0xFFFFFFFFFFFFFFFFn);
    wasm._pcg64_advance(this._wasmState, high, low);
    return this;
  }

  /**
   * Return a jumped copy of the generator.
   *
   * The jumped generator is advanced by a large number of steps,
   * providing an independent stream.
   */
  jumped(jumps: number = 1): PCG64 {
    this.ensureNotDisposed();
    const newGen = new PCG64(this._seedSequence);
    newGen.setState(this.getState());
    // Jump by 2^64 * jumps (using half the period for safety)
    const jumpSize = 1n << 64n;
    for (let i = 0; i < jumps; i++) {
      newGen.advance(jumpSize);
    }
    return newGen;
  }

  getState(): PCG64State {
    this.ensureNotDisposed();
    const wasm = getWasmModule();
    const statePtr = wasm._malloc(6 * 8);

    try {
      wasm._pcg64_get_state(this._wasmState, statePtr);

      // Read the state values
      const stateLow = BigInt(wasm.HEAPU32[statePtr >> 2]) |
                       (BigInt(wasm.HEAPU32[(statePtr >> 2) + 1]) << 32n);
      const stateHigh = BigInt(wasm.HEAPU32[(statePtr >> 2) + 2]) |
                        (BigInt(wasm.HEAPU32[(statePtr >> 2) + 3]) << 32n);
      const incLow = BigInt(wasm.HEAPU32[(statePtr >> 2) + 4]) |
                     (BigInt(wasm.HEAPU32[(statePtr >> 2) + 5]) << 32n);
      const incHigh = BigInt(wasm.HEAPU32[(statePtr >> 2) + 6]) |
                      (BigInt(wasm.HEAPU32[(statePtr >> 2) + 7]) << 32n);

      return {
        bit_generator: 'PCG64',
        state: {
          state: stateLow | (stateHigh << 64n),
          inc: incLow | (incHigh << 64n),
        },
        has_uint32: wasm.HEAPU32[(statePtr >> 2) + 8],
        uinteger: wasm.HEAPU32[(statePtr >> 2) + 9],
      };
    } finally {
      wasm._free(statePtr);
    }
  }

  setState(state: object): void {
    this.ensureNotDisposed();
    const wasm = getWasmModule();
    const s = state as PCG64State;

    const statePtr = wasm._malloc(6 * 8);
    try {
      const st = BigInt(s.state.state);
      const inc = BigInt(s.state.inc);

      // Write state values to WASM memory
      wasm.HEAPU32[statePtr >> 2] = Number(st & 0xFFFFFFFFn);
      wasm.HEAPU32[(statePtr >> 2) + 1] = Number((st >> 32n) & 0xFFFFFFFFn);
      wasm.HEAPU32[(statePtr >> 2) + 2] = Number((st >> 64n) & 0xFFFFFFFFn);
      wasm.HEAPU32[(statePtr >> 2) + 3] = Number((st >> 96n) & 0xFFFFFFFFn);
      wasm.HEAPU32[(statePtr >> 2) + 4] = Number(inc & 0xFFFFFFFFn);
      wasm.HEAPU32[(statePtr >> 2) + 5] = Number((inc >> 32n) & 0xFFFFFFFFn);
      wasm.HEAPU32[(statePtr >> 2) + 6] = Number((inc >> 64n) & 0xFFFFFFFFn);
      wasm.HEAPU32[(statePtr >> 2) + 7] = Number((inc >> 96n) & 0xFFFFFFFFn);
      wasm.HEAPU32[(statePtr >> 2) + 8] = s.has_uint32;
      wasm.HEAPU32[(statePtr >> 2) + 9] = s.uinteger;

      wasm._pcg64_set_state(this._wasmState, statePtr);
    } finally {
      wasm._free(statePtr);
    }
  }

  spawn(nChildren: number): PCG64[] {
    this.ensureNotDisposed();
    if (!this._seedSequence) {
      throw new Error('Cannot spawn from a BitGenerator without a SeedSequence');
    }
    const childSeqs = this._seedSequence.spawn(nChildren);
    return childSeqs.map(seq => new PCG64(seq));
  }

  dispose(): void {
    if (!this._disposed && this._wasmState !== 0) {
      const wasm = getWasmModule();
      wasm._pcg64_free(this._wasmState);
      this._wasmState = 0;
      this._disposed = true;
    }
  }
}

/* ============ Generator Class ============ */

/**
 * Generator - main interface for random number generation.
 *
 * Provides methods for generating random samples from various
 * probability distributions. Uses a BitGenerator (default: PCG64)
 * for the underlying random bits.
 *
 * @example
 * ```typescript
 * // Create a generator
 * const rng = new Generator(new PCG64(12345));
 *
 * // Or use the convenience function
 * const rng = default_rng(12345);
 *
 * // Generate samples
 * const uniform = rng.random(100);        // 100 uniform [0, 1)
 * const normal = rng.normal(0, 1, [10, 10]); // 10x10 standard normal
 * const ints = rng.integers(0, 100, 50);  // 50 integers in [0, 100)
 * const choice = rng.choice([1, 2, 3, 4, 5], 3); // 3 random choices
 * ```
 */
export class Generator {
  private _bitGenerator: BitGenerator;
  private _wasmBitgen: number = 0;
  private _disposed: boolean = false;

  constructor(bitGenerator?: BitGenerator) {
    this._bitGenerator = bitGenerator ?? new PCG64();
    this._initWasmBitgen();
  }

  private _initWasmBitgen(): void {
    const wasm = getWasmModule();
    // Allocate bitgen_t structure and initialize it
    this._wasmBitgen = wasm._malloc(40); // sizeof(bitgen_t)
    if (this._wasmBitgen === 0) {
      throw new Error('Failed to allocate bitgen_t');
    }
    wasm._pcg64_init_bitgen(this._wasmBitgen, this._bitGenerator.wasmStatePtr);
  }

  private ensureNotDisposed(): void {
    if (this._disposed) {
      throw new Error('Generator has been disposed');
    }
  }

  /**
   * Access the underlying BitGenerator.
   */
  get bitGenerator(): BitGenerator {
    return this._bitGenerator;
  }

  /**
   * Spawn n_children independent Generators.
   *
   * Each child generator uses an independent BitGenerator
   * and produces statistically independent streams.
   */
  spawn(nChildren: number): Generator[] {
    this.ensureNotDisposed();
    const childBitGens = this._bitGenerator.spawn(nChildren);
    return childBitGens.map(bg => new Generator(bg));
  }

  /**
   * Clean up WASM resources.
   */
  dispose(): void {
    if (!this._disposed) {
      if (this._wasmBitgen !== 0) {
        const wasm = getWasmModule();
        wasm._free(this._wasmBitgen);
        this._wasmBitgen = 0;
      }
      this._bitGenerator.dispose();
      this._disposed = true;
    }
  }

  /* ============ Uniform Distributions ============ */

  /**
   * Return random floats in the half-open interval [0.0, 1.0).
   *
   * @param size - Output shape. If null, returns a scalar.
   * @param dtype - Output dtype (Float32 or Float64)
   */
  random(size?: SizeType, dtype: DType = DType.Float64): NDArray | number {
    this.ensureNotDisposed();

    if (size === null || size === undefined) {
      return this._bitGenerator.next_double();
    }

    const shape = typeof size === 'number' ? [size] : size;
    const totalSize = shape.reduce((a, b) => a * b, 1);
    const wasm = getWasmModule();

    // Allocate output buffer
    const itemSize = dtype === DType.Float32 ? 4 : 8;
    const outPtr = wasm._malloc(totalSize * itemSize);

    try {
      if (dtype === DType.Float64) {
        wasm._random_uniform_fill(this._wasmBitgen, totalSize, outPtr);
      } else if (dtype === DType.Float32) {
        wasm._random_uniform_fill_f(this._wasmBitgen, totalSize, outPtr);
      } else {
        throw new Error(`Unsupported dtype for random(): ${dtype}`);
      }

      // Create NDArray from the data
      return this._createArrayFromPtr(outPtr, shape, dtype);
    } finally {
      wasm._free(outPtr);
    }
  }

  /**
   * Return random bytes.
   *
   * @param length - Number of bytes to generate
   */
  bytes(length: number): Uint8Array {
    this.ensureNotDisposed();
    const result = new Uint8Array(length);
    for (let i = 0; i < length; i += 4) {
      const val = this._bitGenerator.next_uint32();
      for (let j = 0; j < 4 && i + j < length; j++) {
        result[i + j] = (val >> (j * 8)) & 0xFF;
      }
    }
    return result;
  }

  /* ============ Integer Methods ============ */

  /**
   * Return random integers from low (inclusive) to high (exclusive).
   *
   * @param low - Lowest integer (or high if high is null)
   * @param high - Upper bound (exclusive)
   * @param size - Output shape
   * @param dtype - Output dtype
   * @param endpoint - If true, include high in the range
   */
  integers(
    low: number,
    high?: number | null,
    size?: SizeType,
    _dtype: DType = DType.Int64,
    endpoint: boolean = false
  ): NDArray | number {
    this.ensureNotDisposed();

    // Handle single argument case: integers(high) -> [0, high)
    if (high === null || high === undefined) {
      high = low;
      low = 0;
    }

    if (endpoint) {
      high = high + 1;
    }

    const wasm = getWasmModule();

    if (size === null || size === undefined) {
      // Use 32-bit version for scalar - safer for JS interop
      return wasm._random_integers32(this._wasmBitgen, low, high - 1);
    }

    const shape = typeof size === 'number' ? [size] : size;
    const totalSize = shape.reduce((a, b) => a * b, 1);

    // Allocate output buffer - use int32 for better JS compatibility
    const outPtr = wasm._malloc(totalSize * 4); // int32_t = 4 bytes

    try {
      wasm._random_integers32_fill(this._wasmBitgen, totalSize, low, high - 1, outPtr);
      // Use Int32 dtype for the result
      return this._createArrayFromPtr(outPtr, shape, DType.Int32);
    } finally {
      wasm._free(outPtr);
    }
  }

  /**
   * Generates a random sample from a given array.
   *
   * @param a - If an NDArray/array, sample from its elements. If int, sample from arange(a).
   * @param size - Output shape
   * @param replace - Whether sampling is with replacement
   * @param p - Probability weights (optional)
   * @param axis - Axis along which to sample (for multi-dimensional a)
   * @param shuffle - Whether to shuffle when replace=false
   */
  async choice(
    a: number | NDArray | number[],
    size?: SizeType,
    replace: boolean = true,
    _p?: NDArray | number[] | null,  // TODO: Implement probability weights
    axis: number = 0,
    _shuffle: boolean = true  // TODO: Implement shuffle option
  ): Promise<NDArray | number> {
    this.ensureNotDisposed();

    let population: number;
    let sourceArray: NDArray | null = null;

    if (typeof a === 'number') {
      population = a;
    } else if (Array.isArray(a)) {
      sourceArray = await NDArray.fromArray(a);
      population = sourceArray.size;
    } else {
      sourceArray = a;
      population = a.shape[axis];
    }

    if (!replace && size !== null && size !== undefined) {
      const sizeNum = typeof size === 'number' ? size : size.reduce((x, y) => x * y, 1);
      if (sizeNum > population) {
        throw new Error('Cannot take a larger sample than population when replace=false');
      }
    }

    if (size === null || size === undefined) {
      // Return single element
      const idx = this.integers(0, population) as number;
      if (typeof a === 'number') {
        return idx;
      }
      return sourceArray!.get(idx);
    }

    const shape = typeof size === 'number' ? [size] : size;
    const totalSize = shape.reduce((x, y) => x * y, 1);

    if (replace) {
      // With replacement: simple random indices
      const indices = this.integers(0, population, totalSize) as NDArray;
      if (typeof a === 'number') {
        return indices.reshape(shape);
      }
      const result = await take(sourceArray!, indices);
      indices.dispose();
      return result.reshape(shape);
    } else {
      // Without replacement: Fisher-Yates shuffle subset
      const indexArray = new Int32Array(population);
      for (let i = 0; i < population; i++) {
        indexArray[i] = i;
      }

      for (let i = 0; i < totalSize; i++) {
        const j = i + (this.integers(0, population - i) as number);
        [indexArray[i], indexArray[j]] = [indexArray[j], indexArray[i]];
      }

      const selectedIndices = await NDArray.fromArray(Array.from(indexArray.slice(0, totalSize)));

      if (typeof a === 'number') {
        return selectedIndices.reshape(shape);
      }

      const result = await take(sourceArray!, selectedIndices);
      selectedIndices.dispose();
      return result.reshape(shape);
    }
  }

  /* ============ Continuous Distributions ============ */

  /**
   * Draw samples from a uniform distribution.
   *
   * @param low - Lower boundary
   * @param high - Upper boundary
   * @param size - Output shape
   */
  uniform(low: number = 0.0, high: number = 1.0, size?: SizeType): NDArray | number {
    this.ensureNotDisposed();

    if (size === null || size === undefined) {
      return low + (high - low) * this._bitGenerator.next_double();
    }

    const result = this.random(size, DType.Float64) as NDArray;
    // Scale: result = low + (high - low) * result
    this._scaleShift(result, high - low, low);
    return result;
  }

  /**
   * Draw samples from a standard Normal distribution (mean=0, stdev=1).
   *
   * @param size - Output shape
   * @param dtype - Output dtype
   */
  standard_normal(size?: SizeType, dtype: DType = DType.Float64): NDArray | number {
    this.ensureNotDisposed();
    const wasm = getWasmModule();

    if (size === null || size === undefined) {
      return wasm._random_standard_normal(this._wasmBitgen);
    }

    const shape = typeof size === 'number' ? [size] : size;
    const totalSize = shape.reduce((a, b) => a * b, 1);
    const itemSize = dtype === DType.Float32 ? 4 : 8;
    const outPtr = wasm._malloc(totalSize * itemSize);

    try {
      if (dtype === DType.Float64) {
        wasm._random_standard_normal_fill(this._wasmBitgen, totalSize, outPtr);
      } else {
        wasm._random_standard_normal_fill_f(this._wasmBitgen, totalSize, outPtr);
      }
      return this._createArrayFromPtr(outPtr, shape, dtype);
    } finally {
      wasm._free(outPtr);
    }
  }

  /**
   * Draw random samples from a normal (Gaussian) distribution.
   *
   * @param loc - Mean of the distribution
   * @param scale - Standard deviation
   * @param size - Output shape
   */
  normal(loc: number = 0.0, scale: number = 1.0, size?: SizeType): NDArray | number {
    this.ensureNotDisposed();

    if (scale < 0) {
      throw new Error('scale must be non-negative');
    }

    if (size === null || size === undefined) {
      return loc + scale * (this.standard_normal() as number);
    }

    const result = this.standard_normal(size) as NDArray;
    this._scaleShift(result, scale, loc);
    return result;
  }

  /**
   * Draw samples from a standard exponential distribution.
   *
   * @param size - Output shape
   * @param dtype - Output dtype
   * @param method - 'zig' for Ziggurat, 'inv' for inverse transform
   */
  standard_exponential(
    size?: SizeType,
    dtype: DType = DType.Float64,
    method: 'zig' | 'inv' = 'zig'
  ): NDArray | number {
    this.ensureNotDisposed();
    const wasm = getWasmModule();

    if (size === null || size === undefined) {
      return wasm._random_standard_exponential(this._wasmBitgen);
    }

    const shape = typeof size === 'number' ? [size] : size;
    const totalSize = shape.reduce((a, b) => a * b, 1);
    const outPtr = wasm._malloc(totalSize * 8);

    try {
      if (method === 'zig') {
        wasm._random_standard_exponential_fill(this._wasmBitgen, totalSize, outPtr);
      } else {
        wasm._random_standard_exponential_inv_fill(this._wasmBitgen, totalSize, outPtr);
      }
      return this._createArrayFromPtr(outPtr, shape, dtype);
    } finally {
      wasm._free(outPtr);
    }
  }

  /**
   * Draw samples from an exponential distribution.
   *
   * @param scale - Scale parameter (1/lambda)
   * @param size - Output shape
   */
  exponential(scale: number = 1.0, size?: SizeType): NDArray | number {
    this.ensureNotDisposed();

    if (scale < 0) {
      throw new Error('scale must be non-negative');
    }

    if (size === null || size === undefined) {
      return scale * (this.standard_exponential() as number);
    }

    const result = this.standard_exponential(size) as NDArray;
    this._scale(result, scale);
    return result;
  }

  /**
   * Draw samples from a standard Gamma distribution.
   *
   * @param shape - Shape parameter (k > 0)
   * @param size - Output shape
   * @param dtype - Output dtype
   */
  standard_gamma(shape: number, size?: SizeType, dtype: DType = DType.Float64): NDArray | number {
    this.ensureNotDisposed();

    if (shape < 0) {
      throw new Error('shape must be non-negative');
    }

    const wasm = getWasmModule();

    if (size === null || size === undefined) {
      return wasm._random_standard_gamma(this._wasmBitgen, shape);
    }

    const sizeArr = typeof size === 'number' ? [size] : size;
    const totalSize = sizeArr.reduce((a, b) => a * b, 1);
    const outPtr = wasm._malloc(totalSize * 8);

    try {
      wasm._random_standard_gamma_fill(this._wasmBitgen, totalSize, shape, outPtr);
      return this._createArrayFromPtr(outPtr, sizeArr, dtype);
    } finally {
      wasm._free(outPtr);
    }
  }

  /**
   * Draw samples from a Gamma distribution.
   *
   * @param shape - Shape parameter (k > 0)
   * @param scale - Scale parameter (theta > 0)
   * @param size - Output shape
   */
  gamma(shape: number, scale: number = 1.0, size?: SizeType): NDArray | number {
    this.ensureNotDisposed();

    if (shape < 0) {
      throw new Error('shape must be non-negative');
    }
    if (scale < 0) {
      throw new Error('scale must be non-negative');
    }

    if (size === null || size === undefined) {
      return scale * (this.standard_gamma(shape) as number);
    }

    const result = this.standard_gamma(shape, size) as NDArray;
    this._scale(result, scale);
    return result;
  }

  /**
   * Draw samples from a Beta distribution.
   *
   * @param a - Alpha parameter (> 0)
   * @param b - Beta parameter (> 0)
   * @param size - Output shape
   */
  beta(a: number, b: number, size?: SizeType): NDArray | number {
    this.ensureNotDisposed();

    if (a <= 0) throw new Error('a must be positive');
    if (b <= 0) throw new Error('b must be positive');

    const wasm = getWasmModule();

    if (size === null || size === undefined) {
      return wasm._random_beta(this._wasmBitgen, a, b);
    }

    const sizeArr = typeof size === 'number' ? [size] : size;
    const totalSize = sizeArr.reduce((x, y) => x * y, 1);
    const outPtr = wasm._malloc(totalSize * 8);

    try {
      wasm._random_beta_fill(this._wasmBitgen, totalSize, a, b, outPtr);
      return this._createArrayFromPtr(outPtr, sizeArr, DType.Float64);
    } finally {
      wasm._free(outPtr);
    }
  }

  /**
   * Draw samples from a chi-square distribution.
   *
   * @param df - Degrees of freedom
   * @param size - Output shape
   */
  chisquare(df: number, size?: SizeType): NDArray | number {
    this.ensureNotDisposed();

    if (df <= 0) throw new Error('df must be positive');
    return this.gamma(df / 2.0, 2.0, size);
  }

  /**
   * Draw samples from an F distribution.
   *
   * @param dfnum - Degrees of freedom for numerator
   * @param dfden - Degrees of freedom for denominator
   * @param size - Output shape
   */
  f(dfnum: number, dfden: number, size?: SizeType): NDArray | number {
    this.ensureNotDisposed();

    if (dfnum <= 0) throw new Error('dfnum must be positive');
    if (dfden <= 0) throw new Error('dfden must be positive');

    const wasm = getWasmModule();

    if (size === null || size === undefined) {
      return wasm._random_f(this._wasmBitgen, dfnum, dfden);
    }

    // Generate using ratio of chi-square variates
    const sizeArr = typeof size === 'number' ? [size] : size;
    const totalSize = sizeArr.reduce((a, b) => a * b, 1);

    const num = this.chisquare(dfnum, totalSize) as NDArray;
    const den = this.chisquare(dfden, totalSize) as NDArray;

    // result = (num / dfnum) / (den / dfden) = num * dfden / (den * dfnum)
    const result = _createFloat64Array(sizeArr);
    for (let i = 0; i < totalSize; i++) {
      const n = num.get(i) as number;
      const d = den.get(i) as number;
      result.set((n * dfden) / (d * dfnum), i);
    }

    num.dispose();
    den.dispose();

    return result;
  }

  /**
   * Draw samples from a Student's t distribution.
   *
   * @param df - Degrees of freedom
   * @param size - Output shape
   */
  standard_t(df: number, size?: SizeType): NDArray | number {
    this.ensureNotDisposed();

    if (df <= 0) throw new Error('df must be positive');

    const wasm = getWasmModule();

    if (size === null || size === undefined) {
      return wasm._random_standard_t(this._wasmBitgen, df);
    }

    // Generate using ratio of normal and chi-square
    const sizeArr = typeof size === 'number' ? [size] : size;
    const totalSize = sizeArr.reduce((a, b) => a * b, 1);

    const normal = this.standard_normal(totalSize) as NDArray;
    const chi2 = this.chisquare(df, totalSize) as NDArray;

    const result = _createFloat64Array(sizeArr);
    for (let i = 0; i < totalSize; i++) {
      const n = normal.get(i) as number;
      const c = chi2.get(i) as number;
      result.set(n / Math.sqrt(c / df), i);
    }

    normal.dispose();
    chi2.dispose();

    return result;
  }

  /**
   * Draw samples from a standard Cauchy distribution.
   *
   * @param size - Output shape
   */
  standard_cauchy(size?: SizeType): NDArray | number {
    this.ensureNotDisposed();
    const wasm = getWasmModule();

    if (size === null || size === undefined) {
      return wasm._random_standard_cauchy(this._wasmBitgen);
    }

    // Generate using ratio of two standard normals
    const sizeArr = typeof size === 'number' ? [size] : size;
    const totalSize = sizeArr.reduce((a, b) => a * b, 1);

    const n1 = this.standard_normal(totalSize) as NDArray;
    const n2 = this.standard_normal(totalSize) as NDArray;

    const result = _createFloat64Array(sizeArr);
    for (let i = 0; i < totalSize; i++) {
      result.set((n1.get(i) as number) / (n2.get(i) as number), i);
    }

    n1.dispose();
    n2.dispose();

    return result;
  }

  /**
   * Draw samples from a Pareto II (Lomax) distribution.
   *
   * @param a - Shape parameter
   * @param size - Output shape
   */
  pareto(a: number, size?: SizeType): NDArray | number {
    this.ensureNotDisposed();

    if (a <= 0) throw new Error('a must be positive');

    const wasm = getWasmModule();

    if (size === null || size === undefined) {
      return wasm._random_pareto(this._wasmBitgen, a);
    }

    const sizeArr = typeof size === 'number' ? [size] : size;
    const totalSize = sizeArr.reduce((x, y) => x * y, 1);

    const exp = this.standard_exponential(totalSize) as NDArray;
    const result = _createFloat64Array(sizeArr);

    for (let i = 0; i < totalSize; i++) {
      result.set(Math.exp((exp.get(i) as number) / a) - 1.0, i);
    }

    exp.dispose();
    return result;
  }

  /**
   * Draw samples from a Weibull distribution.
   *
   * @param a - Shape parameter
   * @param size - Output shape
   */
  weibull(a: number, size?: SizeType): NDArray | number {
    this.ensureNotDisposed();

    if (a < 0) throw new Error('a must be non-negative');

    const wasm = getWasmModule();

    if (size === null || size === undefined) {
      return wasm._random_weibull(this._wasmBitgen, a);
    }

    const sizeArr = typeof size === 'number' ? [size] : size;
    const totalSize = sizeArr.reduce((x, y) => x * y, 1);

    const exp = this.standard_exponential(totalSize) as NDArray;
    const result = _createFloat64Array(sizeArr);

    for (let i = 0; i < totalSize; i++) {
      result.set(Math.pow(exp.get(i) as number, 1.0 / a), i);
    }

    exp.dispose();
    return result;
  }

  /**
   * Draw samples from a Laplace (double exponential) distribution.
   *
   * @param loc - Location parameter
   * @param scale - Scale parameter
   * @param size - Output shape
   */
  laplace(loc: number = 0.0, scale: number = 1.0, size?: SizeType): NDArray | number {
    this.ensureNotDisposed();

    if (scale < 0) throw new Error('scale must be non-negative');

    const wasm = getWasmModule();

    if (size === null || size === undefined) {
      return wasm._random_laplace(this._wasmBitgen, loc, scale);
    }

    const sizeArr = typeof size === 'number' ? [size] : size;
    const totalSize = sizeArr.reduce((x, y) => x * y, 1);
    const result = _createFloat64Array(sizeArr);

    for (let i = 0; i < totalSize; i++) {
      const u = this._bitGenerator.next_double();
      if (u < 0.5) {
        result.set(loc + scale * Math.log(2.0 * u), i);
      } else {
        result.set(loc - scale * Math.log(2.0 * (1.0 - u)), i);
      }
    }

    return result;
  }

  /**
   * Draw samples from a log-normal distribution.
   *
   * @param mean - Mean of the underlying normal distribution
   * @param sigma - Standard deviation of the underlying normal
   * @param size - Output shape
   */
  lognormal(mean: number = 0.0, sigma: number = 1.0, size?: SizeType): NDArray | number {
    this.ensureNotDisposed();

    if (sigma < 0) throw new Error('sigma must be non-negative');

    if (size === null || size === undefined) {
      return Math.exp(mean + sigma * (this.standard_normal() as number));
    }

    const sizeArr = typeof size === 'number' ? [size] : size;
    const totalSize = sizeArr.reduce((x, y) => x * y, 1);
    const normal = this.normal(mean, sigma, totalSize) as NDArray;

    const result = _createFloat64Array(sizeArr);
    for (let i = 0; i < totalSize; i++) {
      result.set(Math.exp(normal.get(i) as number), i);
    }

    normal.dispose();
    return result;
  }

  /**
   * Draw samples from a Rayleigh distribution.
   *
   * @param scale - Scale parameter
   * @param size - Output shape
   */
  rayleigh(scale: number = 1.0, size?: SizeType): NDArray | number {
    this.ensureNotDisposed();

    if (scale < 0) throw new Error('scale must be non-negative');

    const wasm = getWasmModule();

    if (size === null || size === undefined) {
      return wasm._random_rayleigh(this._wasmBitgen, scale);
    }

    const sizeArr = typeof size === 'number' ? [size] : size;
    const totalSize = sizeArr.reduce((x, y) => x * y, 1);
    const exp = this.standard_exponential(totalSize) as NDArray;

    const result = _createFloat64Array(sizeArr);
    for (let i = 0; i < totalSize; i++) {
      result.set(scale * Math.sqrt(2.0 * (exp.get(i) as number)), i);
    }

    exp.dispose();
    return result;
  }

  /* ============ Discrete Distributions ============ */

  /**
   * Draw samples from a binomial distribution.
   *
   * @param n - Number of trials
   * @param p - Probability of success
   * @param size - Output shape
   */
  binomial(n: number, p: number, size?: SizeType): NDArray | number {
    this.ensureNotDisposed();

    if (n < 0) throw new Error('n must be non-negative');
    if (p < 0 || p > 1) throw new Error('p must be in [0, 1]');

    const wasm = getWasmModule();

    if (size === null || size === undefined) {
      return Number(wasm._random_binomial(this._wasmBitgen, p, n));
    }

    const sizeArr = typeof size === 'number' ? [size] : size;
    const totalSize = sizeArr.reduce((x, y) => x * y, 1);
    const result = _createInt64Array(sizeArr);

    for (let i = 0; i < totalSize; i++) {
      result.set(Number(wasm._random_binomial(this._wasmBitgen, p, n)), i);
    }

    return result;
  }

  /**
   * Draw samples from a Poisson distribution.
   *
   * @param lam - Expected number of events
   * @param size - Output shape
   */
  poisson(lam: number = 1.0, size?: SizeType): NDArray | number {
    this.ensureNotDisposed();

    if (lam < 0) throw new Error('lam must be non-negative');

    const wasm = getWasmModule();

    if (size === null || size === undefined) {
      return Number(wasm._random_poisson(this._wasmBitgen, lam));
    }

    const sizeArr = typeof size === 'number' ? [size] : size;
    const totalSize = sizeArr.reduce((x, y) => x * y, 1);
    const result = _createInt64Array(sizeArr);

    for (let i = 0; i < totalSize; i++) {
      result.set(Number(wasm._random_poisson(this._wasmBitgen, lam)), i);
    }

    return result;
  }

  /**
   * Draw samples from a geometric distribution.
   *
   * @param p - Probability of success
   * @param size - Output shape
   */
  geometric(p: number, size?: SizeType): NDArray | number {
    this.ensureNotDisposed();

    if (p <= 0 || p > 1) throw new Error('p must be in (0, 1]');

    const wasm = getWasmModule();

    if (size === null || size === undefined) {
      return Number(wasm._random_geometric(this._wasmBitgen, p));
    }

    const sizeArr = typeof size === 'number' ? [size] : size;
    const totalSize = sizeArr.reduce((x, y) => x * y, 1);
    const result = _createInt64Array(sizeArr);

    for (let i = 0; i < totalSize; i++) {
      result.set(Number(wasm._random_geometric(this._wasmBitgen, p)), i);
    }

    return result;
  }

  /**
   * Draw samples from a negative binomial distribution.
   *
   * @param n - Number of successes
   * @param p - Probability of success
   * @param size - Output shape
   */
  negative_binomial(n: number, p: number, size?: SizeType): NDArray | number {
    this.ensureNotDisposed();

    if (n <= 0) throw new Error('n must be positive');
    if (p <= 0 || p >= 1) throw new Error('p must be in (0, 1)');

    const wasm = getWasmModule();

    if (size === null || size === undefined) {
      return Number(wasm._random_negative_binomial(this._wasmBitgen, n, p));
    }

    const sizeArr = typeof size === 'number' ? [size] : size;
    const totalSize = sizeArr.reduce((x, y) => x * y, 1);
    const result = _createInt64Array(sizeArr);

    for (let i = 0; i < totalSize; i++) {
      result.set(Number(wasm._random_negative_binomial(this._wasmBitgen, n, p)), i);
    }

    return result;
  }

  /**
   * Draw samples from a hypergeometric distribution.
   *
   * @param ngood - Number of good items
   * @param nbad - Number of bad items
   * @param nsample - Number of items sampled
   * @param size - Output shape
   */
  hypergeometric(ngood: number, nbad: number, nsample: number, size?: SizeType): NDArray | number {
    this.ensureNotDisposed();

    if (ngood < 0) throw new Error('ngood must be non-negative');
    if (nbad < 0) throw new Error('nbad must be non-negative');
    if (nsample < 0) throw new Error('nsample must be non-negative');
    if (nsample > ngood + nbad) throw new Error('nsample must not exceed ngood + nbad');

    const wasm = getWasmModule();

    if (size === null || size === undefined) {
      return Number(wasm._random_hypergeometric(this._wasmBitgen, ngood, nbad, nsample));
    }

    const sizeArr = typeof size === 'number' ? [size] : size;
    const totalSize = sizeArr.reduce((x, y) => x * y, 1);
    const result = _createInt64Array(sizeArr);

    for (let i = 0; i < totalSize; i++) {
      result.set(Number(wasm._random_hypergeometric(this._wasmBitgen, ngood, nbad, nsample)), i);
    }

    return result;
  }

  /* ============ Permutation Methods ============ */

  /**
   * Randomly shuffle elements along an axis.
   * Modifies the array in-place.
   *
   * @param x - Array to shuffle
   * @param axis - Axis along which to shuffle
   */
  async shuffle(x: NDArray, axis: number = 0): Promise<void> {
    this.ensureNotDisposed();

    const n = x.shape[axis];
    if (n <= 1) return;

    // Fisher-Yates shuffle
    for (let i = n - 1; i > 0; i--) {
      const j = this.integers(0, i + 1) as number;
      if (i !== j) {
        // Swap elements at positions i and j along axis
        await this._swapAlongAxis(x, i, j, axis);
      }
    }
  }

  /**
   * Randomly permute a sequence, or return a permuted range.
   *
   * @param x - If int, permute arange(x). If array, permute its elements.
   * @param axis - Axis along which to permute
   */
  async permutation(x: number | NDArray, axis: number = 0): Promise<NDArray> {
    this.ensureNotDisposed();

    let arr: NDArray;
    if (typeof x === 'number') {
      arr = await NDArray.arange(0, x);
    } else {
      arr = await x.copy();
    }

    await this.shuffle(arr, axis);
    return arr;
  }

  /* ============ Helper Methods ============ */

  /**
   * Create an NDArray from a WASM pointer containing data.
   * @internal
   */
  private _createArrayFromPtr(
    dataPtr: number,
    shape: number[],
    dtype: DType
  ): NDArray {
    const wasm = getWasmModule();
    const totalSize = shape.reduce((a, b) => a * b, 1);

    // Create the array synchronously
    const shapePtr = wasm._malloc(shape.length * 4);
    for (let i = 0; i < shape.length; i++) {
      wasm.setValue(shapePtr + i * 4, shape[i], 'i32');
    }

    const ptr = wasm._ndarray_create(shape.length, shapePtr, dtype);
    wasm._free(shapePtr);

    if (ptr === 0) {
      throw new Error('Failed to create NDArray');
    }

    const result = NDArray._fromPtr(ptr, wasm);

    // Copy data from temporary buffer to the array's data
    const itemSize = dtype === DType.Float32 ? 4 : 8;
    const destPtr = wasm._ndarray_get_data(result._wasmPtr);

    // Use memcpy for efficiency
    const src = new Uint8Array(wasm.HEAPU8.buffer, dataPtr, totalSize * itemSize);
    const dest = new Uint8Array(wasm.HEAPU8.buffer, destPtr, totalSize * itemSize);
    dest.set(src);

    return result;
  }

  /**
   * Scale and shift array values in-place: x = scale * x + shift
   * @internal
   */
  private _scaleShift(arr: NDArray, scale: number, shift: number): void {
    const wasm = getWasmModule();
    const dataPtr = wasm._ndarray_get_data(arr._wasmPtr);
    const size = arr.size;

    for (let i = 0; i < size; i++) {
      const val = wasm.HEAPF64[(dataPtr >> 3) + i];
      wasm.HEAPF64[(dataPtr >> 3) + i] = scale * val + shift;
    }
  }

  /**
   * Scale array values in-place: x = scale * x
   * @internal
   */
  private _scale(arr: NDArray, scale: number): void {
    const wasm = getWasmModule();
    const dataPtr = wasm._ndarray_get_data(arr._wasmPtr);
    const size = arr.size;

    for (let i = 0; i < size; i++) {
      wasm.HEAPF64[(dataPtr >> 3) + i] *= scale;
    }
  }

  /**
   * Swap elements at positions i and j along an axis.
   * @internal
   */
  private async _swapAlongAxis(arr: NDArray, i: number, j: number, _axis: number): Promise<void> {
    // Note: _axis is unused for now - this is a simplified 1D implementation
    const temp = arr.get(i);
    arr.set(arr.get(j), i);
    arr.set(temp, j);
  }
}

/* ============ Module Functions ============ */

/**
 * Construct a new Generator with the given BitGenerator.
 *
 * If seed is:
 * - null/undefined: Use OS entropy
 * - number/bigint: Create SeedSequence from the integer
 * - array: Create SeedSequence from the array
 * - SeedSequence: Use directly
 * - BitGenerator: Use directly
 *
 * @param seed - Seed value or BitGenerator
 * @returns A new Generator instance
 *
 * @example
 * ```typescript
 * // Random seed from OS
 * const rng = default_rng();
 *
 * // Reproducible with integer seed
 * const rng = default_rng(12345);
 *
 * // Reproducible with array seed
 * const rng = default_rng([1, 2, 3, 4]);
 *
 * // Using existing BitGenerator
 * const pcg = new PCG64(12345);
 * const rng = default_rng(pcg);
 * ```
 */
export function default_rng(
  seed?: number | bigint | number[] | SeedSequence | BitGenerator | null
): Generator {
  if (seed instanceof BitGenerator) {
    return new Generator(seed);
  }

  if (seed instanceof SeedSequence) {
    return new Generator(new PCG64(seed));
  }

  return new Generator(new PCG64(seed));
}

/* ============ Convenience Functions ============ */

// Default generator for module-level functions
let _defaultGenerator: Generator | null = null;

function getDefaultGenerator(): Generator {
  if (_defaultGenerator === null) {
    _defaultGenerator = default_rng();
  }
  return _defaultGenerator;
}

/**
 * Seed the default random generator.
 *
 * @param seed - Seed value
 */
export function seed(seed: number | bigint | number[]): void {
  _defaultGenerator = default_rng(seed);
}

/**
 * Return random floats in the half-open interval [0.0, 1.0).
 *
 * @param size - Output shape
 */
export function random(size?: SizeType): NDArray | number {
  return getDefaultGenerator().random(size);
}

/**
 * Draw samples from a standard normal distribution.
 *
 * @param size - Output shape
 */
export function randn(...args: number[]): NDArray | number {
  if (args.length === 0) {
    return getDefaultGenerator().standard_normal();
  }
  return getDefaultGenerator().standard_normal(args);
}

/**
 * Return random integers from low (inclusive) to high (exclusive).
 *
 * @param low - Lower bound
 * @param high - Upper bound (optional)
 * @param size - Output shape
 */
export function randint(
  low: number,
  high?: number | null,
  size?: SizeType
): NDArray | number {
  return getDefaultGenerator().integers(low, high, size);
}

/* ============ Async Initialization Helper ============ */

/**
 * Initialize the random module.
 *
 * This must be called before using any random functions.
 * It loads the WASM module and sets up the default generator.
 */
export async function initRandom(): Promise<void> {
  await loadWasmModule();
}

/* ============ Internal Helpers ============ */

/**
 * Create a float64 array synchronously using the pre-loaded WASM module.
 * @internal
 */
function _createFloat64Array(shape: number[]): NDArray {
  const module = getWasmModule();

  const shapePtr = module._malloc(shape.length * 4);
  for (let i = 0; i < shape.length; i++) {
    module.setValue(shapePtr + i * 4, shape[i], 'i32');
  }

  const ptr = module._ndarray_create(shape.length, shapePtr, DType.Float64);
  module._free(shapePtr);

  if (ptr === 0) {
    throw new Error('Failed to create float64 array');
  }

  return NDArray._fromPtr(ptr, module);
}

/**
 * Create an int64 array synchronously using the pre-loaded WASM module.
 * @internal
 */
function _createInt64Array(shape: number[]): NDArray {
  const module = getWasmModule();

  const shapePtr = module._malloc(shape.length * 4);
  for (let i = 0; i < shape.length; i++) {
    module.setValue(shapePtr + i * 4, shape[i], 'i32');
  }

  const ptr = module._ndarray_create(shape.length, shapePtr, DType.Int64);
  module._free(shapePtr);

  if (ptr === 0) {
    throw new Error('Failed to create int64 array');
  }

  return NDArray._fromPtr(ptr, module);
}

