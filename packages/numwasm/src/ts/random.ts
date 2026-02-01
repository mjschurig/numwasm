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

import { NDArray } from "./_core/NDArray.js";
import { DType } from "./types.js";
import { getWasmModule, loadWasmModule } from "./wasm-loader.js";
import { take } from "./indexing.js";

/* ============ Type Definitions ============ */

/**
 * Size specification for random output arrays.
 */
export type SizeType = number | number[] | null;

/**
 * State object for PCG64 generator.
 */
export interface PCG64State {
  bit_generator: "PCG64";
  state: {
    state: bigint;
    inc: bigint;
  };
  has_uint32: number;
  uinteger: number;
}

/**
 * State object for MT19937 generator.
 */
export interface MT19937State {
  bit_generator: "MT19937";
  state: {
    key: Uint32Array;
    pos: number;
  };
  has_uint32: number;
  uinteger: number;
}

/**
 * State object for Philox generator.
 */
export interface PhiloxState {
  bit_generator: "Philox";
  state: {
    counter: BigUint64Array;
    key: BigUint64Array;
  };
  buffer: BigUint64Array;
  buffer_pos: number;
  has_uint32: number;
  uinteger: number;
}

/**
 * State object for SFC64 generator.
 */
export interface SFC64State {
  bit_generator: "SFC64";
  state: {
    a: bigint;
    b: bigint;
    c: bigint;
    w: bigint;
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
    poolSize: number = 4,
  ) {
    this._spawnKey = spawnKey;
    this._poolSize = poolSize;

    if (entropy === null || entropy === undefined) {
      // Use crypto.getRandomValues for OS entropy
      this._entropy = new Uint32Array(poolSize);
      if (typeof crypto !== "undefined" && crypto.getRandomValues) {
        crypto.getRandomValues(this._entropy);
      } else {
        // Fallback: use Date.now() + Math.random()
        for (let i = 0; i < poolSize; i++) {
          this._entropy[i] =
            ((Date.now() * Math.random()) >>> 0) ^
            ((Math.random() * 0xffffffff) >>> 0);
        }
      }
    } else if (typeof entropy === "number") {
      // Single integer seed - split into low and high 32-bit parts
      this._entropy = new Uint32Array([
        entropy >>> 0,
        Math.floor(entropy / 0x100000000) >>> 0,
      ]);
    } else if (typeof entropy === "bigint") {
      // BigInt seed
      this._entropy = new Uint32Array([
        Number(entropy & 0xffffffffn),
        Number((entropy >> 32n) & 0xffffffffn),
        Number((entropy >> 64n) & 0xffffffffn),
        Number((entropy >> 96n) & 0xffffffffn),
      ]);
    } else if (entropy instanceof Uint32Array) {
      this._entropy = new Uint32Array(entropy);
    } else {
      // Array of numbers
      this._entropy = new Uint32Array(entropy.map((n) => n >>> 0));
    }
  }

  /**
   * Generate state words for seeding a BitGenerator.
   *
   * @param nWords - Number of 32-bit words to generate
   * @param dtype - Output type ('uint32' or 'uint64')
   */
  generateState(
    nWords: number,
    dtype: "uint32" | "uint64" = "uint32",
  ): Uint32Array | BigUint64Array {
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
        entropyPtr,
        this._entropy.length,
        spawnKeyPtr,
        this._spawnKey.length,
        this._poolSize,
        outPtr,
        nWords,
      );

      // Copy result
      const result = new Uint32Array(nWords);
      for (let i = 0; i < nWords; i++) {
        result[i] = wasm.HEAPU32[(outPtr >> 2) + i];
      }

      if (dtype === "uint64") {
        const result64 = new BigUint64Array(Math.floor(nWords / 2));
        for (let i = 0; i < result64.length; i++) {
          result64[i] =
            BigInt(result[i * 2]) | (BigInt(result[i * 2 + 1]) << 32n);
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
      children.push(
        new SeedSequence(
          this._entropy,
          [...this._spawnKey, this._nChildrenSpawned + i],
          this._poolSize,
        ),
      );
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
      throw new Error("Failed to create PCG64 state: memory allocation failed");
    }

    if (seed instanceof SeedSequence) {
      this._seedSequence = seed;
    } else {
      this._seedSequence = new SeedSequence(seed);
    }

    // Generate 8 uint32 values for 4 x 32-bit parts of two 128-bit numbers
    const state = this._seedSequence.generateState(8, "uint32") as Uint32Array;
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
      state[0],
      state[1],
      state[2],
      state[3], // seed_hh, seed_hl, seed_lh, seed_ll
      state[4],
      state[5],
      state[6],
      state[7], // inc_hh, inc_hl, inc_lh, inc_ll
    );
  }

  private ensureNotDisposed(): void {
    if (this._disposed) {
      throw new Error("PCG64 has been disposed");
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
    const low = Number(delta & 0xffffffffffffffffn);
    const high = Number((delta >> 64n) & 0xffffffffffffffffn);
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
      const stateLow =
        BigInt(wasm.HEAPU32[statePtr >> 2]) |
        (BigInt(wasm.HEAPU32[(statePtr >> 2) + 1]) << 32n);
      const stateHigh =
        BigInt(wasm.HEAPU32[(statePtr >> 2) + 2]) |
        (BigInt(wasm.HEAPU32[(statePtr >> 2) + 3]) << 32n);
      const incLow =
        BigInt(wasm.HEAPU32[(statePtr >> 2) + 4]) |
        (BigInt(wasm.HEAPU32[(statePtr >> 2) + 5]) << 32n);
      const incHigh =
        BigInt(wasm.HEAPU32[(statePtr >> 2) + 6]) |
        (BigInt(wasm.HEAPU32[(statePtr >> 2) + 7]) << 32n);

      return {
        bit_generator: "PCG64",
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
      wasm.HEAPU32[statePtr >> 2] = Number(st & 0xffffffffn);
      wasm.HEAPU32[(statePtr >> 2) + 1] = Number((st >> 32n) & 0xffffffffn);
      wasm.HEAPU32[(statePtr >> 2) + 2] = Number((st >> 64n) & 0xffffffffn);
      wasm.HEAPU32[(statePtr >> 2) + 3] = Number((st >> 96n) & 0xffffffffn);
      wasm.HEAPU32[(statePtr >> 2) + 4] = Number(inc & 0xffffffffn);
      wasm.HEAPU32[(statePtr >> 2) + 5] = Number((inc >> 32n) & 0xffffffffn);
      wasm.HEAPU32[(statePtr >> 2) + 6] = Number((inc >> 64n) & 0xffffffffn);
      wasm.HEAPU32[(statePtr >> 2) + 7] = Number((inc >> 96n) & 0xffffffffn);
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
      throw new Error(
        "Cannot spawn from a BitGenerator without a SeedSequence",
      );
    }
    const childSeqs = this._seedSequence.spawn(nChildren);
    return childSeqs.map((seq) => new PCG64(seq));
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

/* ============ PCG64DXSM BitGenerator ============ */

/**
 * State object for PCG64DXSM generator.
 */
export interface PCG64DXSMState {
  bit_generator: "PCG64DXSM";
  state: {
    state: bigint;
    inc: bigint;
  };
  has_uint32: number;
  uinteger: number;
}

/**
 * PCG64DXSM BitGenerator - PCG with DXSM output function.
 *
 * PCG64DXSM is a variant of PCG64 that uses the "Double Xor Shift Multiply"
 * (DXSM) output function. This variant has better statistical properties
 * for some applications and is the default in NumPy 1.17+.
 *
 * Uses the same 128-bit state as PCG64 but with the DXSM output transform:
 * 1. XOR high and low halves
 * 2. Multiply by an odd constant
 * 3. Right rotation
 *
 * @example
 * ```typescript
 * // Create with integer seed
 * const rng = new PCG64DXSM(12345);
 *
 * // Create with OS entropy
 * const rng2 = new PCG64DXSM();
 *
 * // Use with Generator
 * const gen = new Generator(new PCG64DXSM(12345));
 * ```
 */
export class PCG64DXSM extends BitGenerator {
  private _seedSequence: SeedSequence | null = null;
  private _disposed: boolean = false;
  private _stateHigh: bigint = 0n;
  private _stateLow: bigint = 0n;
  private _incHigh: bigint = 0n;
  private _incLow: bigint = 0n;
  private _hasUint32: number = 0;
  private _uinteger: number = 0;

  // DXSM multiplier constant
  private static readonly MULT = 0xda942042e4dd58b5n;

  constructor(seed?: number | bigint | number[] | SeedSequence | null) {
    super();

    if (seed instanceof SeedSequence) {
      this._seedSequence = seed;
    } else {
      this._seedSequence = new SeedSequence(seed);
    }

    // Generate 8 uint32 values for 4 x 32-bit parts of two 128-bit numbers
    const state = this._seedSequence.generateState(8, "uint32") as Uint32Array;
    this._initFromState(state);
  }

  private _initFromState(state: Uint32Array): void {
    // PCG64DXSM needs 128-bit seed and 128-bit increment
    // Combine 32-bit parts into 64-bit values
    this._stateLow =
      BigInt(state[2] >>> 0) | (BigInt(state[3] >>> 0) << 32n);
    this._stateHigh =
      BigInt(state[0] >>> 0) | (BigInt(state[1] >>> 0) << 32n);
    this._incLow =
      BigInt(state[6] >>> 0) | (BigInt(state[7] >>> 0) << 32n);
    this._incHigh =
      BigInt(state[4] >>> 0) | (BigInt(state[5] >>> 0) << 32n);

    // Ensure increment is odd (required for full period)
    this._incLow |= 1n;
  }

  private ensureNotDisposed(): void {
    if (this._disposed) {
      throw new Error("PCG64DXSM has been disposed");
    }
  }

  /**
   * Apply the DXSM (Double Xor Shift Multiply) output function.
   * This transforms the state into a high-quality random output.
   */
  private _dxsmOutput(stateHigh: bigint, stateLow: bigint): bigint {
    // hi = (low >> 32n) with top bit set
    let hi = (stateLow >> 32n) | (1n << 63n);

    // XOR with stateHigh
    hi ^= stateHigh;

    // Multiply
    hi = this._mul64(hi, PCG64DXSM.MULT);

    // Final XOR
    hi ^= hi >> 32n;
    hi = this._mul64(hi, PCG64DXSM.MULT);
    hi ^= hi >> 32n;

    return hi & 0xffffffffffffffffn;
  }

  /**
   * 64-bit multiply keeping only lower 64 bits.
   */
  private _mul64(a: bigint, b: bigint): bigint {
    return (a * b) & 0xffffffffffffffffn;
  }

  /**
   * Advance the internal state by one step.
   */
  private _step(): void {
    // 128-bit multiply by the LCG multiplier
    // state = state * mult + inc (mod 2^128)
    // Since we're using BigInt, we can do this directly
    const stateLo = this._stateLow;
    const stateHi = this._stateHigh;
    const multLow = 0x4385df649fccf645n;
    const multHigh = 0x2360ed051fc65da4n;

    // Compute full 128-bit product: state * mult
    const fullProd = stateLo * multLow +
                     ((stateHi * multLow + stateLo * multHigh) << 64n);
    const fullWithInc = fullProd + this._incLow + (this._incHigh << 64n);

    this._stateLow = fullWithInc & 0xffffffffffffffffn;
    this._stateHigh = (fullWithInc >> 64n) & 0xffffffffffffffffn;
  }

  next_uint64(): bigint {
    this.ensureNotDisposed();

    // Get output from current state
    const output = this._dxsmOutput(this._stateHigh, this._stateLow);

    // Step the generator
    this._step();

    return output;
  }

  next_uint32(): number {
    this.ensureNotDisposed();

    if (this._hasUint32) {
      this._hasUint32 = 0;
      return this._uinteger >>> 0;
    }

    const val = this.next_uint64();
    this._uinteger = Number((val >> 32n) & 0xffffffffn);
    this._hasUint32 = 1;
    return Number(val & 0xffffffffn) >>> 0;
  }

  next_double(): number {
    this.ensureNotDisposed();
    const val = this.next_uint64();
    // Convert to double in [0, 1): divide by 2^64
    return Number(val) / 18446744073709551616.0;
  }

  /**
   * Advance the state by delta steps.
   * Uses the classic O(log n) PCG advance algorithm.
   */
  advance(delta: bigint): this {
    this.ensureNotDisposed();

    // Use 128-bit modular arithmetic to advance by delta steps
    // new_state = mult^delta * state + (mult^delta - 1) / (mult - 1) * inc

    if (delta === 0n) return this;

    let curMult = 0x2360ed051fc65da44385df649fccf645n; // Full 128-bit multiplier
    let curAdd = this._incLow | (this._incHigh << 64n);
    let accMult = 1n;
    let accAdd = 0n;
    const mod = 1n << 128n;

    while (delta > 0n) {
      if (delta & 1n) {
        accMult = (accMult * curMult) % mod;
        accAdd = (accAdd * curMult + curAdd) % mod;
      }
      curAdd = ((curMult + 1n) * curAdd) % mod;
      curMult = (curMult * curMult) % mod;
      delta >>= 1n;
    }

    const state = this._stateLow | (this._stateHigh << 64n);
    const newState = ((accMult * state) % mod + accAdd) % mod;

    this._stateLow = newState & 0xffffffffffffffffn;
    this._stateHigh = (newState >> 64n) & 0xffffffffffffffffn;

    return this;
  }

  /**
   * Return a jumped copy of the generator.
   */
  jumped(jumps: number = 1): PCG64DXSM {
    this.ensureNotDisposed();
    const newGen = new PCG64DXSM(this._seedSequence);
    newGen.setState(this.getState());
    const jumpSize = 1n << 64n;
    for (let i = 0; i < jumps; i++) {
      newGen.advance(jumpSize);
    }
    return newGen;
  }

  getState(): PCG64DXSMState {
    this.ensureNotDisposed();
    return {
      bit_generator: "PCG64DXSM",
      state: {
        state: this._stateLow | (this._stateHigh << 64n),
        inc: this._incLow | (this._incHigh << 64n),
      },
      has_uint32: this._hasUint32,
      uinteger: this._uinteger,
    };
  }

  setState(state: object): void {
    this.ensureNotDisposed();
    const s = state as PCG64DXSMState;

    const st = BigInt(s.state.state);
    const inc = BigInt(s.state.inc);

    this._stateLow = st & 0xffffffffffffffffn;
    this._stateHigh = (st >> 64n) & 0xffffffffffffffffn;
    this._incLow = inc & 0xffffffffffffffffn;
    this._incHigh = (inc >> 64n) & 0xffffffffffffffffn;
    this._hasUint32 = s.has_uint32;
    this._uinteger = s.uinteger;
  }

  spawn(nChildren: number): PCG64DXSM[] {
    this.ensureNotDisposed();
    if (!this._seedSequence) {
      throw new Error(
        "Cannot spawn from a BitGenerator without a SeedSequence",
      );
    }
    const childSeqs = this._seedSequence.spawn(nChildren);
    return childSeqs.map((seq) => new PCG64DXSM(seq));
  }

  dispose(): void {
    if (!this._disposed) {
      this._disposed = true;
    }
  }
}

/* ============ SFC64 BitGenerator ============ */

/**
 * SFC64 BitGenerator - Small Fast Chaotic 64-bit generator.
 *
 * Uses a 256-bit state (4 × 64-bit values). Fast with good statistical properties.
 *
 * @example
 * ```typescript
 * // Create with integer seed
 * const sfc = new SFC64(12345);
 *
 * // Create with SeedSequence
 * const ss = new SeedSequence([1, 2, 3]);
 * const sfc2 = new SFC64(ss);
 *
 * // Generate values
 * const u64 = sfc.next_uint64();
 * const dbl = sfc.next_double();
 * ```
 */
export class SFC64 extends BitGenerator {
  private _seedSequence: SeedSequence | null = null;
  private _disposed: boolean = false;

  constructor(seed?: number | bigint | number[] | SeedSequence | null) {
    super();

    const wasm = getWasmModule();
    this._wasmState = wasm._sfc64_create();

    if (this._wasmState === 0) {
      throw new Error("Failed to create SFC64 state: memory allocation failed");
    }

    if (seed instanceof SeedSequence) {
      this._seedSequence = seed;
    } else {
      this._seedSequence = new SeedSequence(seed);
    }

    // Generate 8 uint32 values for 4 × 64-bit state values
    const state = this._seedSequence.generateState(8, "uint32") as Uint32Array;
    this._initFromState(state);
  }

  private _initFromState(state: Uint32Array): void {
    const wasm = getWasmModule();
    // SFC64 needs 4 × 64-bit seed values
    // state[0..7] represents 4 uint64 values as 8 uint32 parts
    const partsPtr = wasm._malloc(8 * 4);
    try {
      for (let i = 0; i < 8; i++) {
        wasm.HEAPU32[(partsPtr >> 2) + i] = state[i];
      }
      wasm._sfc64_seed_parts(this._wasmState, partsPtr);
    } finally {
      wasm._free(partsPtr);
    }
  }

  private ensureNotDisposed(): void {
    if (this._disposed) {
      throw new Error("SFC64 has been disposed");
    }
  }

  next_uint64(): bigint {
    this.ensureNotDisposed();
    const wasm = getWasmModule();

    const highPtr = wasm._malloc(4);
    try {
      const low = wasm._sfc64_next64_parts(this._wasmState, highPtr);
      const high = wasm.HEAPU32[highPtr >> 2];
      return BigInt(low >>> 0) | (BigInt(high) << 32n);
    } finally {
      wasm._free(highPtr);
    }
  }

  next_uint32(): number {
    this.ensureNotDisposed();
    const wasm = getWasmModule();
    return wasm._sfc64_next32(this._wasmState) >>> 0;
  }

  next_double(): number {
    this.ensureNotDisposed();
    const wasm = getWasmModule();
    return wasm._sfc64_next_double(this._wasmState);
  }

  getState(): SFC64State {
    this.ensureNotDisposed();
    const wasm = getWasmModule();
    const statePtr = wasm._malloc(6 * 8);

    try {
      wasm._sfc64_get_state(this._wasmState, statePtr);

      // Read the state values (4 × 64-bit + has_uint32 + uinteger)
      const readUint64 = (offset: number): bigint => {
        const lo = wasm.HEAPU32[(statePtr >> 2) + offset * 2];
        const hi = wasm.HEAPU32[(statePtr >> 2) + offset * 2 + 1];
        return BigInt(lo >>> 0) | (BigInt(hi) << 32n);
      };

      return {
        bit_generator: "SFC64",
        state: {
          a: readUint64(0),
          b: readUint64(1),
          c: readUint64(2),
          w: readUint64(3),
        },
        has_uint32: Number(readUint64(4)),
        uinteger: Number(readUint64(5)),
      };
    } finally {
      wasm._free(statePtr);
    }
  }

  setState(state: object): void {
    this.ensureNotDisposed();
    const wasm = getWasmModule();
    const s = state as SFC64State;

    const statePtr = wasm._malloc(6 * 8);
    try {
      const writeUint64 = (offset: number, val: bigint): void => {
        wasm.HEAPU32[(statePtr >> 2) + offset * 2] = Number(val & 0xffffffffn);
        wasm.HEAPU32[(statePtr >> 2) + offset * 2 + 1] = Number(
          (val >> 32n) & 0xffffffffn,
        );
      };

      writeUint64(0, s.state.a);
      writeUint64(1, s.state.b);
      writeUint64(2, s.state.c);
      writeUint64(3, s.state.w);
      writeUint64(4, BigInt(s.has_uint32));
      writeUint64(5, BigInt(s.uinteger));

      wasm._sfc64_set_state(this._wasmState, statePtr);
    } finally {
      wasm._free(statePtr);
    }
  }

  spawn(nChildren: number): SFC64[] {
    this.ensureNotDisposed();
    if (!this._seedSequence) {
      throw new Error(
        "Cannot spawn from a BitGenerator without a SeedSequence",
      );
    }
    const childSeqs = this._seedSequence.spawn(nChildren);
    return childSeqs.map((seq) => new SFC64(seq));
  }

  dispose(): void {
    if (!this._disposed && this._wasmState !== 0) {
      const wasm = getWasmModule();
      wasm._sfc64_free(this._wasmState);
      this._wasmState = 0;
      this._disposed = true;
    }
  }
}

/* ============ MT19937 BitGenerator ============ */

/**
 * MT19937 BitGenerator - Classic Mersenne Twister.
 *
 * Uses a 624 × 32-bit state array. Period is 2^19937 - 1.
 * While historically popular, PCG64 or SFC64 are recommended for new applications.
 *
 * @example
 * ```typescript
 * // Create with integer seed
 * const mt = new MT19937(12345);
 *
 * // Create with SeedSequence
 * const ss = new SeedSequence([1, 2, 3]);
 * const mt2 = new MT19937(ss);
 *
 * // Generate values
 * const u32 = mt.next_uint32();  // Native 32-bit output
 * const dbl = mt.next_double();
 * ```
 */
export class MT19937 extends BitGenerator {
  private _seedSequence: SeedSequence | null = null;
  private _disposed: boolean = false;

  constructor(seed?: number | bigint | number[] | SeedSequence | null) {
    super();

    const wasm = getWasmModule();
    this._wasmState = wasm._mt19937_create();

    if (this._wasmState === 0) {
      throw new Error(
        "Failed to create MT19937 state: memory allocation failed",
      );
    }

    if (seed instanceof SeedSequence) {
      this._seedSequence = seed;
    } else {
      this._seedSequence = new SeedSequence(seed);
    }

    // Generate 624 uint32 values for the full state
    const state = this._seedSequence.generateState(
      624,
      "uint32",
    ) as Uint32Array;
    this._initFromState(state);
  }

  private _initFromState(state: Uint32Array): void {
    const wasm = getWasmModule();
    // Use array seeding for SeedSequence-generated values
    const statePtr = wasm._malloc(state.length * 4);
    try {
      for (let i = 0; i < state.length; i++) {
        wasm.HEAPU32[(statePtr >> 2) + i] = state[i];
      }
      wasm._mt19937_seed_array(this._wasmState, statePtr, state.length);
    } finally {
      wasm._free(statePtr);
    }
  }

  private ensureNotDisposed(): void {
    if (this._disposed) {
      throw new Error("MT19937 has been disposed");
    }
  }

  next_uint64(): bigint {
    this.ensureNotDisposed();
    const wasm = getWasmModule();

    const highPtr = wasm._malloc(4);
    try {
      const low = wasm._mt19937_next64_parts(this._wasmState, highPtr);
      const high = wasm.HEAPU32[highPtr >> 2];
      return BigInt(low >>> 0) | (BigInt(high) << 32n);
    } finally {
      wasm._free(highPtr);
    }
  }

  next_uint32(): number {
    this.ensureNotDisposed();
    const wasm = getWasmModule();
    return wasm._mt19937_next32(this._wasmState) >>> 0;
  }

  next_double(): number {
    this.ensureNotDisposed();
    const wasm = getWasmModule();
    return wasm._mt19937_next_double(this._wasmState);
  }

  getState(): MT19937State {
    this.ensureNotDisposed();
    const wasm = getWasmModule();

    const keyPtr = wasm._malloc(624 * 4);
    const posPtr = wasm._malloc(4);
    const hasUint32Ptr = wasm._malloc(4);
    const uintegerPtr = wasm._malloc(4);

    try {
      wasm._mt19937_get_state(
        this._wasmState,
        keyPtr,
        posPtr,
        hasUint32Ptr,
        uintegerPtr,
      );

      const key = new Uint32Array(624);
      for (let i = 0; i < 624; i++) {
        key[i] = wasm.HEAPU32[(keyPtr >> 2) + i];
      }

      return {
        bit_generator: "MT19937",
        state: {
          key: key,
          pos: wasm.HEAP32[posPtr >> 2],
        },
        has_uint32: wasm.HEAP32[hasUint32Ptr >> 2],
        uinteger: wasm.HEAPU32[uintegerPtr >> 2],
      };
    } finally {
      wasm._free(keyPtr);
      wasm._free(posPtr);
      wasm._free(hasUint32Ptr);
      wasm._free(uintegerPtr);
    }
  }

  setState(state: object): void {
    this.ensureNotDisposed();
    const wasm = getWasmModule();
    const s = state as MT19937State;

    const keyPtr = wasm._malloc(624 * 4);
    try {
      for (let i = 0; i < 624; i++) {
        wasm.HEAPU32[(keyPtr >> 2) + i] = s.state.key[i];
      }
      wasm._mt19937_set_state(
        this._wasmState,
        keyPtr,
        s.state.pos,
        s.has_uint32,
        s.uinteger,
      );
    } finally {
      wasm._free(keyPtr);
    }
  }

  spawn(nChildren: number): MT19937[] {
    this.ensureNotDisposed();
    if (!this._seedSequence) {
      throw new Error(
        "Cannot spawn from a BitGenerator without a SeedSequence",
      );
    }
    const childSeqs = this._seedSequence.spawn(nChildren);
    return childSeqs.map((seq) => new MT19937(seq));
  }

  dispose(): void {
    if (!this._disposed && this._wasmState !== 0) {
      const wasm = getWasmModule();
      wasm._mt19937_free(this._wasmState);
      this._wasmState = 0;
      this._disposed = true;
    }
  }
}

/* ============ Philox BitGenerator ============ */

/**
 * Philox BitGenerator - Counter-based RNG (Philox4x64-10).
 *
 * Uses a 4 × 64-bit counter and 2 × 64-bit key. Supports efficient jump and advance operations.
 * Well-suited for parallel and GPU applications.
 *
 * @example
 * ```typescript
 * // Create with integer seed
 * const philox = new Philox(12345);
 *
 * // Jump ahead by 2^128 draws
 * philox.jump();
 *
 * // Advance by arbitrary number of steps
 * philox.advance(1000000n);
 *
 * // Generate values
 * const u64 = philox.next_uint64();
 * const dbl = philox.next_double();
 * ```
 */
export class Philox extends BitGenerator {
  private _seedSequence: SeedSequence | null = null;
  private _disposed: boolean = false;

  constructor(seed?: number | bigint | number[] | SeedSequence | null) {
    super();

    const wasm = getWasmModule();
    this._wasmState = wasm._philox_create();

    if (this._wasmState === 0) {
      throw new Error(
        "Failed to create Philox state: memory allocation failed",
      );
    }

    if (seed instanceof SeedSequence) {
      this._seedSequence = seed;
    } else {
      this._seedSequence = new SeedSequence(seed);
    }

    // Generate 4 uint32 values for the 2 × 64-bit key
    const state = this._seedSequence.generateState(4, "uint32") as Uint32Array;
    this._initFromState(state);
  }

  private _initFromState(state: Uint32Array): void {
    const wasm = getWasmModule();
    // Philox needs 2 × 64-bit key values
    const partsPtr = wasm._malloc(4 * 4);
    try {
      for (let i = 0; i < 4; i++) {
        wasm.HEAPU32[(partsPtr >> 2) + i] = state[i];
      }
      wasm._philox_seed_parts(this._wasmState, partsPtr);
    } finally {
      wasm._free(partsPtr);
    }
  }

  private ensureNotDisposed(): void {
    if (this._disposed) {
      throw new Error("Philox has been disposed");
    }
  }

  next_uint64(): bigint {
    this.ensureNotDisposed();
    const wasm = getWasmModule();

    const highPtr = wasm._malloc(4);
    try {
      const low = wasm._philox_next64_parts(this._wasmState, highPtr);
      const high = wasm.HEAPU32[highPtr >> 2];
      return BigInt(low >>> 0) | (BigInt(high) << 32n);
    } finally {
      wasm._free(highPtr);
    }
  }

  next_uint32(): number {
    this.ensureNotDisposed();
    const wasm = getWasmModule();
    return wasm._philox_next32(this._wasmState) >>> 0;
  }

  next_double(): number {
    this.ensureNotDisposed();
    const wasm = getWasmModule();
    return wasm._philox_next_double(this._wasmState);
  }

  /**
   * Jump ahead by 2^128 draws.
   *
   * This is useful for creating independent streams without overlap.
   */
  jump(): this {
    this.ensureNotDisposed();
    const wasm = getWasmModule();
    wasm._philox_jump(this._wasmState);
    return this;
  }

  /**
   * Advance the generator by delta steps.
   *
   * @param delta - Number of steps to advance
   */
  advance(delta: bigint): this {
    this.ensureNotDisposed();
    const wasm = getWasmModule();

    // Use the full 256-bit advance function (works for any delta size)
    const stepPtr = wasm._malloc(4 * 8);
    try {
      // Write delta as 4 × 64-bit values (little-endian)
      const mask64 = 0xffffffffffffffffn;
      for (let i = 0; i < 4; i++) {
        const val = (delta >> BigInt(i * 64)) & mask64;
        wasm.HEAPU32[(stepPtr >> 2) + i * 2] = Number(val & 0xffffffffn);
        wasm.HEAPU32[(stepPtr >> 2) + i * 2 + 1] = Number(
          (val >> 32n) & 0xffffffffn,
        );
      }
      wasm._philox_advance(this._wasmState, stepPtr);
    } finally {
      wasm._free(stepPtr);
    }
    return this;
  }

  /**
   * Return a jumped copy of the generator.
   *
   * @param jumps - Number of jumps (each jump is 2^128 steps)
   */
  jumped(jumps: number = 1): Philox {
    this.ensureNotDisposed();
    const newGen = new Philox(this._seedSequence);
    newGen.setState(this.getState());
    for (let i = 0; i < jumps; i++) {
      newGen.jump();
    }
    return newGen;
  }

  getState(): PhiloxState {
    this.ensureNotDisposed();
    const wasm = getWasmModule();

    const ctrPtr = wasm._malloc(4 * 8);
    const keyPtr = wasm._malloc(2 * 8);
    const bufferPosPtr = wasm._malloc(4);
    const bufferPtr = wasm._malloc(4 * 8);
    const hasUint32Ptr = wasm._malloc(4);
    const uintegerPtr = wasm._malloc(4);

    try {
      wasm._philox_get_state(
        this._wasmState,
        ctrPtr,
        keyPtr,
        bufferPosPtr,
        bufferPtr,
        hasUint32Ptr,
        uintegerPtr,
      );

      const readUint64Array = (ptr: number, count: number): BigUint64Array => {
        const result = new BigUint64Array(count);
        for (let i = 0; i < count; i++) {
          const lo = wasm.HEAPU32[(ptr >> 2) + i * 2];
          const hi = wasm.HEAPU32[(ptr >> 2) + i * 2 + 1];
          result[i] = BigInt(lo >>> 0) | (BigInt(hi) << 32n);
        }
        return result;
      };

      return {
        bit_generator: "Philox",
        state: {
          counter: readUint64Array(ctrPtr, 4),
          key: readUint64Array(keyPtr, 2),
        },
        buffer: readUint64Array(bufferPtr, 4),
        buffer_pos: wasm.HEAP32[bufferPosPtr >> 2],
        has_uint32: wasm.HEAP32[hasUint32Ptr >> 2],
        uinteger: wasm.HEAPU32[uintegerPtr >> 2],
      };
    } finally {
      wasm._free(ctrPtr);
      wasm._free(keyPtr);
      wasm._free(bufferPosPtr);
      wasm._free(bufferPtr);
      wasm._free(hasUint32Ptr);
      wasm._free(uintegerPtr);
    }
  }

  setState(state: object): void {
    this.ensureNotDisposed();
    const wasm = getWasmModule();
    const s = state as PhiloxState;

    const ctrPtr = wasm._malloc(4 * 8);
    const keyPtr = wasm._malloc(2 * 8);
    const bufferPtr = wasm._malloc(4 * 8);

    try {
      const writeUint64Array = (ptr: number, arr: BigUint64Array): void => {
        for (let i = 0; i < arr.length; i++) {
          wasm.HEAPU32[(ptr >> 2) + i * 2] = Number(arr[i] & 0xffffffffn);
          wasm.HEAPU32[(ptr >> 2) + i * 2 + 1] = Number(
            (arr[i] >> 32n) & 0xffffffffn,
          );
        }
      };

      writeUint64Array(ctrPtr, s.state.counter);
      writeUint64Array(keyPtr, s.state.key);
      writeUint64Array(bufferPtr, s.buffer);

      wasm._philox_set_state(
        this._wasmState,
        ctrPtr,
        keyPtr,
        s.buffer_pos,
        bufferPtr,
        s.has_uint32,
        s.uinteger,
      );
    } finally {
      wasm._free(ctrPtr);
      wasm._free(keyPtr);
      wasm._free(bufferPtr);
    }
  }

  spawn(nChildren: number): Philox[] {
    this.ensureNotDisposed();
    if (!this._seedSequence) {
      throw new Error(
        "Cannot spawn from a BitGenerator without a SeedSequence",
      );
    }
    const childSeqs = this._seedSequence.spawn(nChildren);
    return childSeqs.map((seq) => new Philox(seq));
  }

  dispose(): void {
    if (!this._disposed && this._wasmState !== 0) {
      const wasm = getWasmModule();
      wasm._philox_free(this._wasmState);
      this._wasmState = 0;
      this._disposed = true;
    }
  }
}

/* ============ BitGenerator Registry ============ */

/**
 * Registry of available BitGenerator classes.
 */
const BIT_GENERATORS: Record<
  string,
  new (seed?: number | bigint | number[] | SeedSequence | null) => BitGenerator
> = {
  PCG64: PCG64,
  PCG64DXSM: PCG64DXSM,
  MT19937: MT19937,
  PHILOX: Philox,
  SFC64: SFC64,
};

/**
 * Get a BitGenerator class by name.
 *
 * @param name - Name of the BitGenerator (case-insensitive)
 * @returns The BitGenerator constructor
 *
 * @example
 * ```typescript
 * const MT = getBitGenerator('mt19937');
 * const mt = new MT(12345);
 * ```
 */
export function getBitGenerator(
  name: string,
): new (
  seed?: number | bigint | number[] | SeedSequence | null,
) => BitGenerator {
  const normalized = name.toUpperCase();
  const cls = BIT_GENERATORS[normalized];
  if (!cls) {
    const available = Object.keys(BIT_GENERATORS).join(", ");
    throw new Error(`Unknown BitGenerator: ${name}. Available: ${available}`);
  }
  return cls;
}

/**
 * List available BitGenerator names.
 *
 * @returns Array of BitGenerator names
 */
export function listBitGenerators(): string[] {
  return Object.keys(BIT_GENERATORS);
}

/* ============ Helper Functions ============ */

/**
 * Compute Cholesky decomposition of a symmetric positive semi-definite matrix.
 * Returns the lower triangular matrix L such that A = L * L^T.
 * @internal
 */
function choleskyDecomposition(
  A: number[][],
  checkValid: "warn" | "raise" | "ignore" = "warn",
  tol: number = 1e-8,
): number[][] {
  const n = A.length;
  const L: number[][] = Array.from({ length: n }, () => new Array(n).fill(0));

  for (let i = 0; i < n; i++) {
    for (let j = 0; j <= i; j++) {
      let sum = 0;
      for (let k = 0; k < j; k++) {
        sum += L[i][k] * L[j][k];
      }

      if (i === j) {
        const diag = A[i][i] - sum;
        if (diag < -tol) {
          if (checkValid === "raise") {
            throw new Error("covariance matrix is not positive semi-definite");
          } else if (checkValid === "warn") {
            console.warn(
              "covariance matrix is not positive semi-definite, " +
              "results may be unreliable",
            );
          }
          // Use small positive value to continue
          L[i][j] = Math.sqrt(Math.max(diag, tol));
        } else {
          L[i][j] = Math.sqrt(Math.max(diag, 0));
        }
      } else {
        if (L[j][j] !== 0) {
          L[i][j] = (A[i][j] - sum) / L[j][j];
        } else {
          L[i][j] = 0;
        }
      }
    }
  }

  return L;
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
      throw new Error("Failed to allocate bitgen_t");
    }

    // Initialize based on the BitGenerator type
    if (this._bitGenerator instanceof PCG64) {
      wasm._pcg64_init_bitgen(
        this._wasmBitgen,
        this._bitGenerator.wasmStatePtr,
      );
    } else if (this._bitGenerator instanceof SFC64) {
      wasm._sfc64_init_bitgen(
        this._wasmBitgen,
        this._bitGenerator.wasmStatePtr,
      );
    } else if (this._bitGenerator instanceof MT19937) {
      wasm._mt19937_init_bitgen(
        this._wasmBitgen,
        this._bitGenerator.wasmStatePtr,
      );
    } else if (this._bitGenerator instanceof Philox) {
      wasm._philox_init_bitgen(
        this._wasmBitgen,
        this._bitGenerator.wasmStatePtr,
      );
    } else {
      // Fallback for custom BitGenerators - use PCG64 style (may not work correctly)
      wasm._pcg64_init_bitgen(
        this._wasmBitgen,
        this._bitGenerator.wasmStatePtr,
      );
    }
  }

  private ensureNotDisposed(): void {
    if (this._disposed) {
      throw new Error("Generator has been disposed");
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
    return childBitGens.map((bg) => new Generator(bg));
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

    const shape = typeof size === "number" ? [size] : size;
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
        result[i + j] = (val >> (j * 8)) & 0xff;
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
    endpoint: boolean = false,
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

    const shape = typeof size === "number" ? [size] : size;
    const totalSize = shape.reduce((a, b) => a * b, 1);

    // Allocate output buffer - use int32 for better JS compatibility
    const outPtr = wasm._malloc(totalSize * 4); // int32_t = 4 bytes

    try {
      wasm._random_integers32_fill(
        this._wasmBitgen,
        totalSize,
        low,
        high - 1,
        outPtr,
      );
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
    _p?: NDArray | number[] | null, // TODO: Implement probability weights
    axis: number = 0,
    _shuffle: boolean = true, // TODO: Implement shuffle option
  ): Promise<NDArray | number> {
    this.ensureNotDisposed();

    let population: number;
    let sourceArray: NDArray | null = null;

    if (typeof a === "number") {
      population = a;
    } else if (Array.isArray(a)) {
      sourceArray = await NDArray.fromArray(a);
      population = sourceArray.size;
    } else {
      sourceArray = a;
      population = a.shape[axis];
    }

    if (!replace && size !== null && size !== undefined) {
      const sizeNum =
        typeof size === "number" ? size : size.reduce((x, y) => x * y, 1);
      if (sizeNum > population) {
        throw new Error(
          "Cannot take a larger sample than population when replace=false",
        );
      }
    }

    if (size === null || size === undefined) {
      // Return single element
      const idx = this.integers(0, population) as number;
      if (typeof a === "number") {
        return idx;
      }
      return sourceArray!.get(idx);
    }

    const shape = typeof size === "number" ? [size] : size;
    const totalSize = shape.reduce((x, y) => x * y, 1);

    if (replace) {
      // With replacement: simple random indices
      const indices = this.integers(0, population, totalSize) as NDArray;
      if (typeof a === "number") {
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

      const selectedIndices = await NDArray.fromArray(
        Array.from(indexArray.slice(0, totalSize)),
      );

      if (typeof a === "number") {
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
  uniform(
    low: number = 0.0,
    high: number = 1.0,
    size?: SizeType,
  ): NDArray | number {
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
  standard_normal(
    size?: SizeType,
    dtype: DType = DType.Float64,
  ): NDArray | number {
    this.ensureNotDisposed();
    const wasm = getWasmModule();

    if (size === null || size === undefined) {
      return wasm._random_standard_normal(this._wasmBitgen);
    }

    const shape = typeof size === "number" ? [size] : size;
    const totalSize = shape.reduce((a, b) => a * b, 1);
    const itemSize = dtype === DType.Float32 ? 4 : 8;
    const outPtr = wasm._malloc(totalSize * itemSize);

    try {
      if (dtype === DType.Float64) {
        wasm._random_standard_normal_fill(this._wasmBitgen, totalSize, outPtr);
      } else {
        wasm._random_standard_normal_fill_f(
          this._wasmBitgen,
          totalSize,
          outPtr,
        );
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
  normal(
    loc: number = 0.0,
    scale: number = 1.0,
    size?: SizeType,
  ): NDArray | number {
    this.ensureNotDisposed();

    if (scale < 0) {
      throw new Error("scale must be non-negative");
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
    method: "zig" | "inv" = "zig",
  ): NDArray | number {
    this.ensureNotDisposed();
    const wasm = getWasmModule();

    if (size === null || size === undefined) {
      return wasm._random_standard_exponential(this._wasmBitgen);
    }

    const shape = typeof size === "number" ? [size] : size;
    const totalSize = shape.reduce((a, b) => a * b, 1);
    const outPtr = wasm._malloc(totalSize * 8);

    try {
      if (method === "zig") {
        wasm._random_standard_exponential_fill(
          this._wasmBitgen,
          totalSize,
          outPtr,
        );
      } else {
        wasm._random_standard_exponential_inv_fill(
          this._wasmBitgen,
          totalSize,
          outPtr,
        );
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
      throw new Error("scale must be non-negative");
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
  standard_gamma(
    shape: number,
    size?: SizeType,
    dtype: DType = DType.Float64,
  ): NDArray | number {
    this.ensureNotDisposed();

    if (shape < 0) {
      throw new Error("shape must be non-negative");
    }

    const wasm = getWasmModule();

    if (size === null || size === undefined) {
      return wasm._random_standard_gamma(this._wasmBitgen, shape);
    }

    const sizeArr = typeof size === "number" ? [size] : size;
    const totalSize = sizeArr.reduce((a, b) => a * b, 1);
    const outPtr = wasm._malloc(totalSize * 8);

    try {
      wasm._random_standard_gamma_fill(
        this._wasmBitgen,
        totalSize,
        shape,
        outPtr,
      );
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
      throw new Error("shape must be non-negative");
    }
    if (scale < 0) {
      throw new Error("scale must be non-negative");
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

    if (a <= 0) throw new Error("a must be positive");
    if (b <= 0) throw new Error("b must be positive");

    const wasm = getWasmModule();

    if (size === null || size === undefined) {
      return wasm._random_beta(this._wasmBitgen, a, b);
    }

    const sizeArr = typeof size === "number" ? [size] : size;
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

    if (df <= 0) throw new Error("df must be positive");
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

    if (dfnum <= 0) throw new Error("dfnum must be positive");
    if (dfden <= 0) throw new Error("dfden must be positive");

    const wasm = getWasmModule();

    if (size === null || size === undefined) {
      return wasm._random_f(this._wasmBitgen, dfnum, dfden);
    }

    // Generate using ratio of chi-square variates
    const sizeArr = typeof size === "number" ? [size] : size;
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

    if (df <= 0) throw new Error("df must be positive");

    const wasm = getWasmModule();

    if (size === null || size === undefined) {
      return wasm._random_standard_t(this._wasmBitgen, df);
    }

    // Generate using ratio of normal and chi-square
    const sizeArr = typeof size === "number" ? [size] : size;
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
    const sizeArr = typeof size === "number" ? [size] : size;
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

    if (a <= 0) throw new Error("a must be positive");

    const wasm = getWasmModule();

    if (size === null || size === undefined) {
      return wasm._random_pareto(this._wasmBitgen, a);
    }

    const sizeArr = typeof size === "number" ? [size] : size;
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

    if (a < 0) throw new Error("a must be non-negative");

    const wasm = getWasmModule();

    if (size === null || size === undefined) {
      return wasm._random_weibull(this._wasmBitgen, a);
    }

    const sizeArr = typeof size === "number" ? [size] : size;
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
  laplace(
    loc: number = 0.0,
    scale: number = 1.0,
    size?: SizeType,
  ): NDArray | number {
    this.ensureNotDisposed();

    if (scale < 0) throw new Error("scale must be non-negative");

    const wasm = getWasmModule();

    if (size === null || size === undefined) {
      return wasm._random_laplace(this._wasmBitgen, loc, scale);
    }

    const sizeArr = typeof size === "number" ? [size] : size;
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
  lognormal(
    mean: number = 0.0,
    sigma: number = 1.0,
    size?: SizeType,
  ): NDArray | number {
    this.ensureNotDisposed();

    if (sigma < 0) throw new Error("sigma must be non-negative");

    if (size === null || size === undefined) {
      return Math.exp(mean + sigma * (this.standard_normal() as number));
    }

    const sizeArr = typeof size === "number" ? [size] : size;
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

    if (scale < 0) throw new Error("scale must be non-negative");

    const wasm = getWasmModule();

    if (size === null || size === undefined) {
      return wasm._random_rayleigh(this._wasmBitgen, scale);
    }

    const sizeArr = typeof size === "number" ? [size] : size;
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

    if (n < 0) throw new Error("n must be non-negative");
    if (p < 0 || p > 1) throw new Error("p must be in [0, 1]");

    const wasm = getWasmModule();

    if (size === null || size === undefined) {
      return wasm._random_binomial32(this._wasmBitgen, p, n);
    }

    const sizeArr = typeof size === "number" ? [size] : size;
    const totalSize = sizeArr.reduce((x, y) => x * y, 1);
    const result = _createFloat64Array(sizeArr);

    for (let i = 0; i < totalSize; i++) {
      result.set(wasm._random_binomial32(this._wasmBitgen, p, n), i);
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

    if (lam < 0) throw new Error("lam must be non-negative");

    const wasm = getWasmModule();

    if (size === null || size === undefined) {
      return wasm._random_poisson32(this._wasmBitgen, lam);
    }

    const sizeArr = typeof size === "number" ? [size] : size;
    const totalSize = sizeArr.reduce((x, y) => x * y, 1);
    const result = _createFloat64Array(sizeArr);

    for (let i = 0; i < totalSize; i++) {
      result.set(wasm._random_poisson32(this._wasmBitgen, lam), i);
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

    if (p <= 0 || p > 1) throw new Error("p must be in (0, 1]");

    const wasm = getWasmModule();

    if (size === null || size === undefined) {
      return wasm._random_geometric32(this._wasmBitgen, p);
    }

    const sizeArr = typeof size === "number" ? [size] : size;
    const totalSize = sizeArr.reduce((x, y) => x * y, 1);
    const result = _createFloat64Array(sizeArr);

    for (let i = 0; i < totalSize; i++) {
      result.set(wasm._random_geometric32(this._wasmBitgen, p), i);
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

    if (n <= 0) throw new Error("n must be positive");
    if (p <= 0 || p >= 1) throw new Error("p must be in (0, 1)");

    const wasm = getWasmModule();

    if (size === null || size === undefined) {
      return wasm._random_negative_binomial32(this._wasmBitgen, n, p);
    }

    const sizeArr = typeof size === "number" ? [size] : size;
    const totalSize = sizeArr.reduce((x, y) => x * y, 1);
    const result = _createFloat64Array(sizeArr);

    for (let i = 0; i < totalSize; i++) {
      result.set(wasm._random_negative_binomial32(this._wasmBitgen, n, p), i);
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
  hypergeometric(
    ngood: number,
    nbad: number,
    nsample: number,
    size?: SizeType,
  ): NDArray | number {
    this.ensureNotDisposed();

    if (ngood < 0) throw new Error("ngood must be non-negative");
    if (nbad < 0) throw new Error("nbad must be non-negative");
    if (nsample < 0) throw new Error("nsample must be non-negative");
    if (nsample > ngood + nbad)
      throw new Error("nsample must not exceed ngood + nbad");

    const wasm = getWasmModule();

    if (size === null || size === undefined) {
      return wasm._random_hypergeometric32(
        this._wasmBitgen,
        ngood,
        nbad,
        nsample,
      );
    }

    const sizeArr = typeof size === "number" ? [size] : size;
    const totalSize = sizeArr.reduce((x, y) => x * y, 1);
    const result = _createFloat64Array(sizeArr);

    for (let i = 0; i < totalSize; i++) {
      result.set(
        wasm._random_hypergeometric32(this._wasmBitgen, ngood, nbad, nsample),
        i,
      );
    }

    return result;
  }

  /**
   * Draw samples from a Gumbel distribution.
   *
   * The Gumbel distribution is used to model the distribution of the maximum
   * (or minimum) of a number of samples of various distributions.
   *
   * @param loc - Location parameter (mode of the distribution)
   * @param scale - Scale parameter (> 0)
   * @param size - Output shape
   */
  gumbel(loc: number = 0.0, scale: number = 1.0, size?: SizeType): NDArray | number {
    this.ensureNotDisposed();

    if (scale < 0) throw new Error("scale must be non-negative");

    const wasm = getWasmModule();

    if (size === null || size === undefined) {
      return wasm._random_gumbel(this._wasmBitgen, loc, scale);
    }

    const sizeArr = typeof size === "number" ? [size] : size;
    const totalSize = sizeArr.reduce((x, y) => x * y, 1);
    const result = _createFloat64Array(sizeArr);

    for (let i = 0; i < totalSize; i++) {
      result.set(wasm._random_gumbel(this._wasmBitgen, loc, scale), i);
    }

    return result;
  }

  /**
   * Draw samples from a logistic distribution.
   *
   * @param loc - Location parameter
   * @param scale - Scale parameter (> 0)
   * @param size - Output shape
   */
  logistic(loc: number = 0.0, scale: number = 1.0, size?: SizeType): NDArray | number {
    this.ensureNotDisposed();

    if (scale < 0) throw new Error("scale must be non-negative");

    const wasm = getWasmModule();

    if (size === null || size === undefined) {
      return wasm._random_logistic(this._wasmBitgen, loc, scale);
    }

    const sizeArr = typeof size === "number" ? [size] : size;
    const totalSize = sizeArr.reduce((x, y) => x * y, 1);
    const result = _createFloat64Array(sizeArr);

    for (let i = 0; i < totalSize; i++) {
      result.set(wasm._random_logistic(this._wasmBitgen, loc, scale), i);
    }

    return result;
  }

  /**
   * Draw samples from a logarithmic series distribution.
   *
   * @param p - Shape parameter (0 < p < 1)
   * @param size - Output shape
   */
  logseries(p: number, size?: SizeType): NDArray | number {
    this.ensureNotDisposed();

    if (p <= 0 || p >= 1) throw new Error("p must be in (0, 1)");

    const wasm = getWasmModule();

    if (size === null || size === undefined) {
      return wasm._random_logseries32(this._wasmBitgen, p);
    }

    const sizeArr = typeof size === "number" ? [size] : size;
    const totalSize = sizeArr.reduce((x, y) => x * y, 1);
    const result = _createFloat64Array(sizeArr);

    for (let i = 0; i < totalSize; i++) {
      result.set(wasm._random_logseries32(this._wasmBitgen, p), i);
    }

    return result;
  }

  /**
   * Draw samples from a noncentral chi-square distribution.
   *
   * @param df - Degrees of freedom (> 0)
   * @param nonc - Non-centrality parameter (>= 0)
   * @param size - Output shape
   */
  noncentral_chisquare(df: number, nonc: number, size?: SizeType): NDArray | number {
    this.ensureNotDisposed();

    if (df <= 0) throw new Error("df must be positive");
    if (nonc < 0) throw new Error("nonc must be non-negative");

    const wasm = getWasmModule();

    if (size === null || size === undefined) {
      return wasm._random_noncentral_chisquare(this._wasmBitgen, df, nonc);
    }

    const sizeArr = typeof size === "number" ? [size] : size;
    const totalSize = sizeArr.reduce((x, y) => x * y, 1);
    const result = _createFloat64Array(sizeArr);

    for (let i = 0; i < totalSize; i++) {
      result.set(wasm._random_noncentral_chisquare(this._wasmBitgen, df, nonc), i);
    }

    return result;
  }

  /**
   * Draw samples from a noncentral F distribution.
   *
   * @param dfnum - Degrees of freedom for numerator (> 0)
   * @param dfden - Degrees of freedom for denominator (> 0)
   * @param nonc - Non-centrality parameter (>= 0)
   * @param size - Output shape
   */
  noncentral_f(dfnum: number, dfden: number, nonc: number, size?: SizeType): NDArray | number {
    this.ensureNotDisposed();

    if (dfnum <= 0) throw new Error("dfnum must be positive");
    if (dfden <= 0) throw new Error("dfden must be positive");
    if (nonc < 0) throw new Error("nonc must be non-negative");

    const wasm = getWasmModule();

    if (size === null || size === undefined) {
      return wasm._random_noncentral_f(this._wasmBitgen, dfnum, dfden, nonc);
    }

    const sizeArr = typeof size === "number" ? [size] : size;
    const totalSize = sizeArr.reduce((x, y) => x * y, 1);
    const result = _createFloat64Array(sizeArr);

    for (let i = 0; i < totalSize; i++) {
      result.set(wasm._random_noncentral_f(this._wasmBitgen, dfnum, dfden, nonc), i);
    }

    return result;
  }

  /**
   * Draw samples from a power distribution.
   *
   * Samples are drawn from a power distribution with positive exponent a - 1
   * in the interval [0, 1].
   *
   * @param a - Shape parameter (> 0)
   * @param size - Output shape
   */
  power(a: number, size?: SizeType): NDArray | number {
    this.ensureNotDisposed();

    if (a <= 0) throw new Error("a must be positive");

    const wasm = getWasmModule();

    if (size === null || size === undefined) {
      return wasm._random_power(this._wasmBitgen, a);
    }

    const sizeArr = typeof size === "number" ? [size] : size;
    const totalSize = sizeArr.reduce((x, y) => x * y, 1);
    const result = _createFloat64Array(sizeArr);

    for (let i = 0; i < totalSize; i++) {
      result.set(wasm._random_power(this._wasmBitgen, a), i);
    }

    return result;
  }

  /**
   * Draw samples from a triangular distribution.
   *
   * @param left - Lower limit
   * @param mode - Mode (peak) of the distribution (left <= mode <= right)
   * @param right - Upper limit
   * @param size - Output shape
   */
  triangular(left: number, mode: number, right: number, size?: SizeType): NDArray | number {
    this.ensureNotDisposed();

    if (left > mode) throw new Error("left must be <= mode");
    if (mode > right) throw new Error("mode must be <= right");
    if (left === right) throw new Error("left must be < right");

    const wasm = getWasmModule();

    if (size === null || size === undefined) {
      return wasm._random_triangular(this._wasmBitgen, left, mode, right);
    }

    const sizeArr = typeof size === "number" ? [size] : size;
    const totalSize = sizeArr.reduce((x, y) => x * y, 1);
    const result = _createFloat64Array(sizeArr);

    for (let i = 0; i < totalSize; i++) {
      result.set(wasm._random_triangular(this._wasmBitgen, left, mode, right), i);
    }

    return result;
  }

  /**
   * Draw samples from a von Mises distribution.
   *
   * Also known as the circular normal distribution. The von Mises distribution
   * is a continuous probability distribution on the circle.
   *
   * @param mu - Mode (peak) of the distribution, in radians
   * @param kappa - Concentration parameter (>= 0)
   * @param size - Output shape
   */
  vonmises(mu: number, kappa: number, size?: SizeType): NDArray | number {
    this.ensureNotDisposed();

    if (kappa < 0) throw new Error("kappa must be non-negative");

    const wasm = getWasmModule();

    if (size === null || size === undefined) {
      return wasm._random_vonmises(this._wasmBitgen, mu, kappa);
    }

    const sizeArr = typeof size === "number" ? [size] : size;
    const totalSize = sizeArr.reduce((x, y) => x * y, 1);
    const result = _createFloat64Array(sizeArr);

    for (let i = 0; i < totalSize; i++) {
      result.set(wasm._random_vonmises(this._wasmBitgen, mu, kappa), i);
    }

    return result;
  }

  /**
   * Draw samples from a Wald (inverse Gaussian) distribution.
   *
   * @param mean - Mean of the distribution (> 0)
   * @param scale - Scale parameter (> 0)
   * @param size - Output shape
   */
  wald(mean: number, scale: number, size?: SizeType): NDArray | number {
    this.ensureNotDisposed();

    if (mean <= 0) throw new Error("mean must be positive");
    if (scale <= 0) throw new Error("scale must be positive");

    const wasm = getWasmModule();

    if (size === null || size === undefined) {
      return wasm._random_wald(this._wasmBitgen, mean, scale);
    }

    const sizeArr = typeof size === "number" ? [size] : size;
    const totalSize = sizeArr.reduce((x, y) => x * y, 1);
    const result = _createFloat64Array(sizeArr);

    for (let i = 0; i < totalSize; i++) {
      result.set(wasm._random_wald(this._wasmBitgen, mean, scale), i);
    }

    return result;
  }

  /**
   * Draw samples from a Zipf distribution.
   *
   * @param a - Distribution parameter (> 1)
   * @param size - Output shape
   */
  zipf(a: number, size?: SizeType): NDArray | number {
    this.ensureNotDisposed();

    if (a <= 1) throw new Error("a must be > 1");

    const wasm = getWasmModule();

    if (size === null || size === undefined) {
      return wasm._random_zipf32(this._wasmBitgen, a);
    }

    const sizeArr = typeof size === "number" ? [size] : size;
    const totalSize = sizeArr.reduce((x, y) => x * y, 1);
    const result = _createFloat64Array(sizeArr);

    for (let i = 0; i < totalSize; i++) {
      result.set(wasm._random_zipf32(this._wasmBitgen, a), i);
    }

    return result;
  }

  /* ============ Multivariate Distributions ============ */

  /**
   * Draw samples from a Dirichlet distribution.
   *
   * The Dirichlet distribution is a multivariate generalization of the
   * Beta distribution. It is the conjugate prior of the categorical
   * distribution and multinomial distribution.
   *
   * @param alpha - Concentration parameters. Must be > 0.
   * @param size - Number of samples to draw (output shape is (*size, len(alpha)))
   * @returns Array of shape (*size, len(alpha)) if size specified, else (len(alpha),)
   *
   * @example
   * ```typescript
   * const rng = default_rng(12345);
   *
   * // Draw from symmetric Dirichlet
   * const alpha = [1, 1, 1];
   * const sample = await rng.dirichlet(alpha);  // Shape: (3,)
   *
   * // Draw multiple samples
   * const samples = await rng.dirichlet(alpha, 10);  // Shape: (10, 3)
   * ```
   */
  async dirichlet(
    alpha: number[] | NDArray,
    size?: number | number[],
  ): Promise<NDArray> {
    this.ensureNotDisposed();

    // Get alpha values as array
    let alphaArr: number[];
    if (Array.isArray(alpha)) {
      alphaArr = alpha;
    } else {
      alphaArr = [];
      for (let i = 0; i < alpha.size; i++) {
        alphaArr.push(alpha.getFlat(i) as number);
      }
    }

    const k = alphaArr.length;
    if (k < 2) {
      throw new Error("alpha must have at least 2 elements");
    }

    for (let i = 0; i < k; i++) {
      if (alphaArr[i] <= 0) {
        throw new Error("all alpha values must be > 0");
      }
    }

    // Determine output shape
    let numSamples: number;
    let outputShape: number[];

    if (size === undefined || size === null) {
      numSamples = 1;
      outputShape = [k];
    } else if (typeof size === "number") {
      numSamples = size;
      outputShape = [size, k];
    } else {
      numSamples = size.reduce((a, b) => a * b, 1);
      outputShape = [...size, k];
    }

    // Generate samples using the Gamma distribution method:
    // For Dirichlet(alpha), draw X_i ~ Gamma(alpha_i, 1) and normalize
    const result = await NDArray.zeros(outputShape);
    const wasm = getWasmModule();

    for (let s = 0; s < numSamples; s++) {
      let sum = 0;
      const values: number[] = [];

      // Draw from gamma distributions
      for (let i = 0; i < k; i++) {
        const g = wasm._random_gamma(this._wasmBitgen, alphaArr[i], 1.0);
        values.push(g);
        sum += g;
      }

      // Normalize and store
      for (let i = 0; i < k; i++) {
        result.setFlat(s * k + i, values[i] / sum);
      }
    }

    return result;
  }

  /**
   * Draw samples from a multinomial distribution.
   *
   * The multinomial distribution is a generalization of the binomial
   * distribution, modeling the number of occurrences of each of K
   * different outcomes in n independent trials.
   *
   * @param n - Number of trials (experiments)
   * @param pvals - Probabilities of each outcome (must sum to 1)
   * @param size - Number of samples to draw
   * @returns Array of shape (*size, len(pvals)) if size specified, else (len(pvals),)
   *
   * @example
   * ```typescript
   * const rng = default_rng(12345);
   *
   * // Roll a fair die 20 times
   * const pvals = [1/6, 1/6, 1/6, 1/6, 1/6, 1/6];
   * const sample = await rng.multinomial(20, pvals);  // Shape: (6,)
   * // sample might be [4, 3, 2, 5, 3, 3] (counts summing to 20)
   *
   * // Draw multiple experiments
   * const samples = await rng.multinomial(20, pvals, 100);  // Shape: (100, 6)
   * ```
   */
  async multinomial(
    n: number,
    pvals: number[] | NDArray,
    size?: number | number[],
  ): Promise<NDArray> {
    this.ensureNotDisposed();

    if (n < 0 || !Number.isInteger(n)) {
      throw new Error("n must be a non-negative integer");
    }

    // Get probability values as array
    let pArr: number[];
    if (Array.isArray(pvals)) {
      pArr = pvals;
    } else {
      pArr = [];
      for (let i = 0; i < pvals.size; i++) {
        pArr.push(pvals.getFlat(i) as number);
      }
    }

    const k = pArr.length;
    if (k < 2) {
      throw new Error("pvals must have at least 2 elements");
    }

    // Validate probabilities
    let pSum = 0;
    for (let i = 0; i < k; i++) {
      if (pArr[i] < 0 || pArr[i] > 1) {
        throw new Error("all pvals must be in [0, 1]");
      }
      pSum += pArr[i];
    }
    if (Math.abs(pSum - 1.0) > 1e-10) {
      throw new Error("pvals must sum to 1");
    }

    // Determine output shape
    let numSamples: number;
    let outputShape: number[];

    if (size === undefined || size === null) {
      numSamples = 1;
      outputShape = [k];
    } else if (typeof size === "number") {
      numSamples = size;
      outputShape = [size, k];
    } else {
      numSamples = size.reduce((a, b) => a * b, 1);
      outputShape = [...size, k];
    }

    const result = await NDArray.zeros(outputShape, { dtype: DType.Int32 });
    const wasm = getWasmModule();

    for (let s = 0; s < numSamples; s++) {
      // Use the conditional binomial method:
      // For each category, draw from binomial with remaining trials and
      // conditional probability
      let remaining = n;
      let pRemaining = 1.0;

      for (let i = 0; i < k - 1 && remaining > 0; i++) {
        // Conditional probability for category i
        const pCond = pArr[i] / pRemaining;

        // Draw from binomial(remaining, pCond)
        const count = wasm._random_binomial(this._wasmBitgen, remaining, pCond);

        result.setFlat(s * k + i, count);
        remaining -= count;
        pRemaining -= pArr[i];
      }

      // Last category gets remaining trials
      result.setFlat(s * k + (k - 1), remaining);
    }

    return result;
  }

  /**
   * Draw samples from a multivariate normal distribution.
   *
   * The multivariate normal distribution is a generalization of the
   * one-dimensional normal distribution to higher dimensions.
   *
   * @param mean - Mean of the distribution (length N)
   * @param cov - Covariance matrix (NxN, symmetric positive semi-definite)
   * @param size - Number of samples to draw
   * @param checkValid - Method for checking positive semi-definiteness
   * @param tol - Tolerance for singular value check
   * @returns Array of shape (*size, N) if size specified, else (N,)
   *
   * @example
   * ```typescript
   * const rng = default_rng(12345);
   *
   * // 2D multivariate normal
   * const mean = [0, 0];
   * const cov = [[1, 0.5], [0.5, 1]];
   * const sample = await rng.multivariate_normal(mean, cov);  // Shape: (2,)
   *
   * // Draw multiple samples
   * const samples = await rng.multivariate_normal(mean, cov, 1000);  // Shape: (1000, 2)
   * ```
   */
  async multivariate_normal(
    mean: number[] | NDArray,
    cov: number[][] | NDArray,
    size?: number | number[],
    checkValid: "warn" | "raise" | "ignore" = "warn",
    tol: number = 1e-8,
  ): Promise<NDArray> {
    this.ensureNotDisposed();

    // Get mean values as array
    let meanArr: number[];
    if (Array.isArray(mean)) {
      meanArr = mean;
    } else {
      meanArr = [];
      for (let i = 0; i < mean.size; i++) {
        meanArr.push(mean.getFlat(i) as number);
      }
    }

    const n = meanArr.length;
    if (n < 1) {
      throw new Error("mean must have at least 1 element");
    }

    // Get covariance matrix as 2D array
    let covArr: number[][];
    if (Array.isArray(cov)) {
      covArr = cov;
    } else {
      // NDArray - reshape to 2D
      if (cov.ndim !== 2) {
        throw new Error("cov must be a 2D array");
      }
      covArr = [];
      for (let i = 0; i < cov.shape[0]; i++) {
        const row: number[] = [];
        for (let j = 0; j < cov.shape[1]; j++) {
          row.push(cov.get(i, j) as number);
        }
        covArr.push(row);
      }
    }

    // Validate covariance matrix dimensions
    if (covArr.length !== n || covArr[0].length !== n) {
      throw new Error(`cov must be ${n}x${n} to match mean length ${n}`);
    }

    // Check symmetry
    for (let i = 0; i < n; i++) {
      for (let j = i + 1; j < n; j++) {
        if (Math.abs(covArr[i][j] - covArr[j][i]) > tol) {
          if (checkValid === "raise") {
            throw new Error("cov matrix is not symmetric");
          } else if (checkValid === "warn") {
            console.warn("cov matrix is not symmetric, using (cov + cov.T) / 2");
          }
          // Symmetrize
          const avg = (covArr[i][j] + covArr[j][i]) / 2;
          covArr[i][j] = avg;
          covArr[j][i] = avg;
        }
      }
    }

    // Compute Cholesky decomposition for sampling
    // L such that cov = L * L^T
    const L = choleskyDecomposition(covArr, checkValid, tol);

    // Determine output shape
    let numSamples: number;
    let outputShape: number[];

    if (size === undefined || size === null) {
      numSamples = 1;
      outputShape = [n];
    } else if (typeof size === "number") {
      numSamples = size;
      outputShape = [size, n];
    } else {
      numSamples = size.reduce((a, b) => a * b, 1);
      outputShape = [...size, n];
    }

    const result = await NDArray.zeros(outputShape);
    const wasm = getWasmModule();

    for (let s = 0; s < numSamples; s++) {
      // Draw standard normal samples
      const z: number[] = [];
      for (let i = 0; i < n; i++) {
        z.push(wasm._random_normal(this._wasmBitgen, 0, 1));
      }

      // Transform: x = mean + L * z
      for (let i = 0; i < n; i++) {
        let val = meanArr[i];
        for (let j = 0; j <= i; j++) {
          val += L[i][j] * z[j];
        }
        result.setFlat(s * n + i, val);
      }
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
    if (typeof x === "number") {
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
    dtype: DType,
  ): NDArray {
    const wasm = getWasmModule();
    const totalSize = shape.reduce((a, b) => a * b, 1);

    // Create the array synchronously
    const shapePtr = wasm._malloc(shape.length * 4);
    for (let i = 0; i < shape.length; i++) {
      wasm.setValue(shapePtr + i * 4, shape[i], "i32");
    }

    const ptr = wasm._ndarray_create(shape.length, shapePtr, dtype);
    wasm._free(shapePtr);

    if (ptr === 0) {
      throw new Error("Failed to create NDArray");
    }

    const result = NDArray._fromPtr(ptr, wasm);

    // Copy data from temporary buffer to the array's data
    const itemSize = dtype === DType.Float32 ? 4 : 8;
    const destPtr = wasm._ndarray_get_data(result._wasmPtr);

    // Use memcpy for efficiency
    const src = new Uint8Array(
      wasm.HEAPU8.buffer,
      dataPtr,
      totalSize * itemSize,
    );
    const dest = new Uint8Array(
      wasm.HEAPU8.buffer,
      destPtr,
      totalSize * itemSize,
    );
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
  private async _swapAlongAxis(
    arr: NDArray,
    i: number,
    j: number,
    _axis: number,
  ): Promise<void> {
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
 * @param options - Options object with optional bitGenerator name
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
 *
 * // Using a specific BitGenerator type
 * const rng = default_rng(12345, { bitGenerator: 'MT19937' });
 * const rng = default_rng(12345, { bitGenerator: 'Philox' });
 * const rng = default_rng(12345, { bitGenerator: 'SFC64' });
 * ```
 */
export function default_rng(
  seed?: number | bigint | number[] | SeedSequence | BitGenerator | null,
  options?: { bitGenerator?: string },
): Generator {
  if (seed instanceof BitGenerator) {
    return new Generator(seed);
  }

  const bitGeneratorName = options?.bitGenerator ?? "PCG64";
  const BitGenClass = getBitGenerator(bitGeneratorName);

  if (seed instanceof SeedSequence) {
    return new Generator(new BitGenClass(seed));
  }

  return new Generator(new BitGenClass(seed));
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
  size?: SizeType,
): NDArray | number {
  return getDefaultGenerator().integers(low, high, size);
}

/**
 * Generates a random sample from a given array or range.
 *
 * @param a - If an array, a random sample is generated from its elements.
 *            If an int, the random sample is generated from arange(a).
 * @param size - Output shape. If null (default), a single value is returned.
 * @param replace - Whether the sample is with or without replacement. Default is true.
 * @param p - Probabilities associated with each entry in a. If not given, the sample
 *            assumes a uniform distribution over all entries.
 * @returns Random samples from the array
 *
 * @example
 * // Choose 3 random elements from [1, 2, 3, 4, 5]
 * const c = await choice([1, 2, 3, 4, 5], 3);
 *
 * @example
 * // Choose 5 elements from range(10) with replacement
 * const c = await choice(10, 5, true);
 */
export async function choice(
  a: number | number[] | NDArray,
  size?: SizeType,
  replace: boolean = true,
  p?: number[] | NDArray | null,
): Promise<NDArray | number> {
  return getDefaultGenerator().choice(a, size, replace, p);
}

/**
 * Randomly shuffle elements of an array in-place along an axis.
 *
 * @param x - Array to shuffle (modified in-place)
 * @param axis - Axis along which to shuffle (default 0)
 *
 * @example
 * const arr = await NDArray.fromArray([1, 2, 3, 4, 5]);
 * await shuffle(arr);
 * // arr is now shuffled in-place
 */
export async function shuffle(x: NDArray, axis: number = 0): Promise<void> {
  return getDefaultGenerator().shuffle(x, axis);
}

/**
 * Randomly permute a sequence, or return a permuted range.
 *
 * @param x - If an integer, randomly permute arange(x).
 *            If an array, make a copy and shuffle the elements randomly.
 * @param axis - Axis along which to permute (default 0)
 * @returns Permuted array (new array, not in-place)
 *
 * @example
 * // Permute range(5)
 * const p = await permutation(5);
 * // e.g., [3, 1, 4, 0, 2]
 *
 * @example
 * // Permute an array
 * const arr = await NDArray.fromArray([10, 20, 30, 40]);
 * const p = await permutation(arr);
 */
export async function permutation(
  x: number | NDArray,
  axis: number = 0,
): Promise<NDArray> {
  return getDefaultGenerator().permutation(x, axis);
}

/**
 * Return random bytes.
 *
 * @param length - Number of random bytes to generate
 * @returns A Uint8Array of random bytes
 *
 * @example
 * // Generate 16 random bytes
 * const b = bytes(16);
 */
export function bytes(length: number): Uint8Array {
  return getDefaultGenerator().bytes(length);
}

/* ============ Legacy Functions ============ */

/**
 * Random values in a given shape.
 *
 * Create an array of the given shape and populate it with random samples
 * from a uniform distribution over [0, 1).
 *
 * @deprecated This is a legacy function from numpy.random. Use random.random()
 *   or Generator.random() instead. For new code, prefer default_rng().
 *
 * @param d0 - First dimension
 * @param d1, d2, ... - Additional dimensions
 * @returns Random values in the specified shape
 *
 * @example
 * // Random float
 * rand();
 *
 * // Random 1D array
 * rand(5);
 *
 * // Random 2D array
 * rand(3, 4);
 */
export function rand(...shape: number[]): NDArray | number {
  if (shape.length === 0) {
    return getDefaultGenerator().random(null) as number;
  }
  return getDefaultGenerator().random(shape);
}

/**
 * Return random floats in the half-open interval [0.0, 1.0).
 *
 * @deprecated This is a legacy alias for random_sample. Use random.random()
 *   or Generator.random() instead. For new code, prefer default_rng().
 *
 * @param size - Output shape
 * @returns Random values in [0, 1)
 */
export function ranf(size?: SizeType): NDArray | number {
  return getDefaultGenerator().random(size);
}

/**
 * Return random floats in the half-open interval [0.0, 1.0).
 *
 * @deprecated This is a legacy function. Use random.random()
 *   or Generator.random() instead. For new code, prefer default_rng().
 *
 * @param size - Output shape
 * @returns Random values in [0, 1)
 */
export function random_sample(size?: SizeType): NDArray | number {
  return getDefaultGenerator().random(size);
}

/**
 * Return random floats in the half-open interval [0.0, 1.0).
 *
 * @deprecated This is a legacy alias for random_sample. Use random.random()
 *   or Generator.random() instead. For new code, prefer default_rng().
 *
 * @param size - Output shape
 * @returns Random values in [0, 1)
 */
export function sample(size?: SizeType): NDArray | number {
  return getDefaultGenerator().random(size);
}

/**
 * Return random integers from low (inclusive) to high (exclusive).
 *
 * @deprecated This is a legacy function. Use random.randint() or
 *   Generator.integers() instead. For new code, prefer default_rng().
 *
 * @param low - Lowest integer (unless high is None, then range is [0, low))
 * @param high - Upper bound (exclusive). If None, range is [0, low)
 * @param size - Output shape
 * @returns Random integers in [low, high)
 */
export function random_integers(
  low: number,
  high?: number | null,
  size?: SizeType,
): NDArray | number {
  // In legacy random_integers, the range is inclusive: [low, high]
  // If high is None, range is [1, low] (inclusive)
  if (high === null || high === undefined) {
    return getDefaultGenerator().integers(1, low + 1, size);
  }
  return getDefaultGenerator().integers(low, high + 1, size);
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
    module.setValue(shapePtr + i * 4, shape[i], "i32");
  }

  const ptr = module._ndarray_create(shape.length, shapePtr, DType.Float64);
  module._free(shapePtr);

  if (ptr === 0) {
    throw new Error("Failed to create float64 array");
  }

  return NDArray._fromPtr(ptr, module);
}

/**
 * Create an int64 array synchronously using the pre-loaded WASM module.
 * @internal
 * @unused - Kept for future int64 distribution implementations
 */
export function _createInt64Array(shape: number[]): NDArray {
  const module = getWasmModule();

  const shapePtr = module._malloc(shape.length * 4);
  for (let i = 0; i < shape.length; i++) {
    module.setValue(shapePtr + i * 4, shape[i], "i32");
  }

  const ptr = module._ndarray_create(shape.length, shapePtr, DType.Int64);
  module._free(shapePtr);

  if (ptr === 0) {
    throw new Error("Failed to create int64 array");
  }

  return NDArray._fromPtr(ptr, module);
}
