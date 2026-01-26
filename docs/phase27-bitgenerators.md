# Phase 27: Additional BitGenerators Implementation Plan

Complete implementation roadmap for additional pseudo-random number generator engines: MT19937, Philox, and SFC64.

---

## ⚠️ Implementation Guideline

**The original NumPy source code is available in `/numpy`.**

Key reference files:
- `numpy/random/bit_generator.pyx` - BitGenerator base class
- `numpy/random/mt19937.pyx` - Mersenne Twister
- `numpy/random/philox.pyx` - Philox counter-based RNG
- `numpy/random/sfc64.pyx` - SFC64 (Small Fast Chaotic)
- `numpy/random/src/` - C implementations

Implementation should follow NumPy's algorithms and error handling for consistency.

---

## Current State (Pre-Phase 27)

```
Already Implemented:
├── Generator class
│   ├── All distribution methods (30+)
│   └── Permutation methods (shuffle, permutation, choice)
├── PCG64 BitGenerator (default)
│   ├── 128-bit state
│   ├── Jump/advance methods
│   └── WASM implementation
└── SeedSequence
    ├── Entropy processing
    └── Spawning for parallel streams

Missing BitGenerators:
├── MT19937 (Mersenne Twister)
├── Philox (counter-based)
└── SFC64 (Small Fast Chaotic)
```

---

## Phase 27 Dependency Tree

```
PHASE 27: ADDITIONAL BITGENERATORS
│
├── 27.1 MT19937 (Mersenne Twister) (TypeScript + WASM)
│   ├── 27.1.1 State management (624 × 32-bit words)
│   ├── 27.1.2 Twist operation
│   ├── 27.1.3 Tempering transformation
│   ├── 27.1.4 Seeding (array or scalar)
│   └── 27.1.5 State serialization
│
│   Dependencies: BitGenerator base, SeedSequence
│
├── 27.2 Philox (Counter-based) (TypeScript + WASM)
│   ├── 27.2.1 4×64-bit counter state
│   ├── 27.2.2 2×64-bit key
│   ├── 27.2.3 Philox round function
│   ├── 27.2.4 Jump/advance methods
│   └── 27.2.5 State serialization
│
│   Dependencies: BitGenerator base, SeedSequence
│
├── 27.3 SFC64 (Small Fast Chaotic) (TypeScript + WASM)
│   ├── 27.3.1 256-bit state (4×64-bit)
│   ├── 27.3.2 Update function
│   ├── 27.3.3 Seeding
│   └── 27.3.4 State serialization
│
│   Dependencies: BitGenerator base, SeedSequence
│
└── 27.4 BitGenerator Registry (TypeScript)
    ├── 27.4.1 Registration system
    ├── 27.4.2 Factory function
    └── 27.4.3 default_rng updates

    Dependencies: 27.1-27.3
```

---

## Detailed Implementation Specifications

### 27.1 MT19937 (Mersenne Twister)

#### 27.1.1 Overview

MT19937 is the classic Mersenne Twister algorithm with period 2^19937 - 1. It uses a 624×32-bit state array.

**File:** `src/ts/random/mt19937.ts` (new file)

```typescript
import { BitGenerator, SeedType } from './bitgenerator.js';
import { SeedSequence } from './seedsequence.js';

/**
 * MT19937 (Mersenne Twister) pseudo-random number generator.
 *
 * This is the classic Mersenne Twister algorithm with period 2^19937-1.
 * While widely used historically, PCG64 or SFC64 are recommended for
 * new applications due to better statistical properties and performance.
 *
 * @example
 * const rng = new Generator(new MT19937(12345));
 * rng.random();  // Random float in [0, 1)
 */
export class MT19937 extends BitGenerator {
  /** State array size */
  private static readonly N = 624;

  /** Middle point for twist */
  private static readonly M = 397;

  /** Constant vector a */
  private static readonly MATRIX_A = 0x9908b0df;

  /** Most significant w-r bits */
  private static readonly UPPER_MASK = 0x80000000;

  /** Least significant r bits */
  private static readonly LOWER_MASK = 0x7fffffff;

  /** State array */
  private mt: Uint32Array;

  /** Current position in state array */
  private mti: number;

  constructor(seed?: SeedType) {
    super();
    this.mt = new Uint32Array(MT19937.N);
    this.mti = MT19937.N + 1;  // Uninitialized

    if (seed !== undefined) {
      this.seed(seed);
    }
  }

  /**
   * Seed the generator.
   */
  seed(seed: SeedType): void {
    if (seed instanceof SeedSequence) {
      this._seedFromSequence(seed);
    } else if (Array.isArray(seed)) {
      this._seedFromArray(seed);
    } else {
      this._seedFromScalar(seed);
    }
  }

  /**
   * Seed from a scalar value.
   */
  private _seedFromScalar(seed: number): void {
    this.mt[0] = seed >>> 0;

    for (let i = 1; i < MT19937.N; i++) {
      // mt[i] = 1812433253 * (mt[i-1] ^ (mt[i-1] >> 30)) + i
      const s = this.mt[i - 1] ^ (this.mt[i - 1] >>> 30);
      this.mt[i] = (((((s & 0xffff0000) >>> 16) * 1812433253) << 16) +
                    (s & 0x0000ffff) * 1812433253 + i) >>> 0;
    }

    this.mti = MT19937.N;
  }

  /**
   * Seed from an array of values.
   */
  private _seedFromArray(initKey: number[]): void {
    this._seedFromScalar(19650218);

    let i = 1;
    let j = 0;
    let k = Math.max(MT19937.N, initKey.length);

    for (; k > 0; k--) {
      const s = this.mt[i - 1] ^ (this.mt[i - 1] >>> 30);
      this.mt[i] = (this.mt[i] ^
                    (((((s & 0xffff0000) >>> 16) * 1664525) << 16) +
                     (s & 0x0000ffff) * 1664525) +
                    initKey[j] + j) >>> 0;
      i++;
      j++;

      if (i >= MT19937.N) {
        this.mt[0] = this.mt[MT19937.N - 1];
        i = 1;
      }
      if (j >= initKey.length) {
        j = 0;
      }
    }

    for (k = MT19937.N - 1; k > 0; k--) {
      const s = this.mt[i - 1] ^ (this.mt[i - 1] >>> 30);
      this.mt[i] = (this.mt[i] ^
                    (((((s & 0xffff0000) >>> 16) * 1566083941) << 16) +
                     (s & 0x0000ffff) * 1566083941) - i) >>> 0;
      i++;

      if (i >= MT19937.N) {
        this.mt[0] = this.mt[MT19937.N - 1];
        i = 1;
      }
    }

    this.mt[0] = 0x80000000;  // MSB is 1; assuring non-zero initial array
  }

  /**
   * Seed from a SeedSequence.
   */
  private _seedFromSequence(seq: SeedSequence): void {
    const initKey = seq.generate_state(MT19937.N);
    this._seedFromArray(Array.from(initKey));
  }

  /**
   * Generate a random 32-bit unsigned integer.
   */
  nextUint32(): number {
    if (this.mti >= MT19937.N) {
      this._twist();
    }

    let y = this.mt[this.mti++];

    // Tempering
    y ^= (y >>> 11);
    y ^= (y << 7) & 0x9d2c5680;
    y ^= (y << 15) & 0xefc60000;
    y ^= (y >>> 18);

    return y >>> 0;
  }

  /**
   * Generate a random 64-bit unsigned integer.
   * Combines two 32-bit values.
   */
  nextUint64(): bigint {
    const high = BigInt(this.nextUint32());
    const low = BigInt(this.nextUint32());
    return (high << 32n) | low;
  }

  /**
   * Generate a random double in [0, 1).
   */
  nextDouble(): number {
    // Use 53 bits for full double precision
    const a = this.nextUint32() >>> 5;  // 27 bits
    const b = this.nextUint32() >>> 6;  // 26 bits
    return (a * 67108864.0 + b) / 9007199254740992.0;
  }

  /**
   * Perform the twist operation to regenerate the state.
   */
  private _twist(): void {
    const N = MT19937.N;
    const M = MT19937.M;
    const MATRIX_A = MT19937.MATRIX_A;
    const UPPER_MASK = MT19937.UPPER_MASK;
    const LOWER_MASK = MT19937.LOWER_MASK;

    const mag01 = [0, MATRIX_A];

    for (let i = 0; i < N - M; i++) {
      const y = (this.mt[i] & UPPER_MASK) | (this.mt[i + 1] & LOWER_MASK);
      this.mt[i] = this.mt[i + M] ^ (y >>> 1) ^ mag01[y & 0x1];
    }

    for (let i = N - M; i < N - 1; i++) {
      const y = (this.mt[i] & UPPER_MASK) | (this.mt[i + 1] & LOWER_MASK);
      this.mt[i] = this.mt[i + (M - N)] ^ (y >>> 1) ^ mag01[y & 0x1];
    }

    const y = (this.mt[N - 1] & UPPER_MASK) | (this.mt[0] & LOWER_MASK);
    this.mt[N - 1] = this.mt[M - 1] ^ (y >>> 1) ^ mag01[y & 0x1];

    this.mti = 0;
  }

  /**
   * Get the current state for serialization.
   */
  get state(): { key: Uint32Array; pos: number } {
    return {
      key: this.mt.slice(),
      pos: this.mti
    };
  }

  /**
   * Set the state for deserialization.
   */
  set state(value: { key: Uint32Array; pos: number }) {
    if (value.key.length !== MT19937.N) {
      throw new ValueError(`State key must have ${MT19937.N} elements`);
    }

    this.mt.set(value.key);
    this.mti = value.pos;
  }
}
```

---

### 27.2 Philox (Counter-based RNG)

#### 27.2.1 Overview

Philox is a counter-based RNG that uses a simple, parallelizable algorithm. It's particularly suited for GPU and parallel applications.

**File:** `src/ts/random/philox.ts` (new file)

```typescript
import { BitGenerator, SeedType } from './bitgenerator.js';
import { SeedSequence } from './seedsequence.js';

/**
 * Philox (4×64) counter-based pseudo-random number generator.
 *
 * Philox is a counter-based RNG that uses a simple multiplication-based
 * mixing function. It's designed for parallel random number generation.
 *
 * Features:
 * - Period: 2^256
 * - Supports efficient jumping and advancing
 * - Deterministic parallel streams via counter manipulation
 *
 * @example
 * const rng = new Generator(new Philox(12345));
 * rng.random();
 *
 * // Jump ahead by 2^128 draws
 * const philox = new Philox(12345);
 * philox.jump();
 */
export class Philox extends BitGenerator {
  /** Number of rounds (10 for Philox4x64-10) */
  private static readonly ROUNDS = 10;

  /** Multiplier constants */
  private static readonly PHILOX_M0 = 0xD2E7470EE14C6C93n;
  private static readonly PHILOX_M1 = 0xCA5A826395121157n;

  /** Weyl sequence constants */
  private static readonly PHILOX_W0 = 0x9E3779B97F4A7C15n;
  private static readonly PHILOX_W1 = 0xBB67AE8584CAA73Bn;

  /** 4×64-bit counter */
  private counter: BigUint64Array;

  /** 2×64-bit key */
  private key: BigUint64Array;

  /** Buffer for generated values */
  private buffer: BigUint64Array;
  private bufferPos: number;

  constructor(seed?: SeedType, counter?: bigint | bigint[]) {
    super();

    this.counter = new BigUint64Array(4);
    this.key = new BigUint64Array(2);
    this.buffer = new BigUint64Array(4);
    this.bufferPos = 4;  // Empty buffer

    if (seed !== undefined) {
      this.seed(seed);
    }

    if (counter !== undefined) {
      this.setCounter(counter);
    }
  }

  /**
   * Seed the generator.
   */
  seed(seed: SeedType): void {
    if (seed instanceof SeedSequence) {
      const state = seed.generate_state(4, 'uint64');
      this.key[0] = state[0];
      this.key[1] = state[1];
    } else if (typeof seed === 'bigint') {
      this.key[0] = seed;
      this.key[1] = 0n;
    } else if (typeof seed === 'number') {
      this.key[0] = BigInt(seed >>> 0);
      this.key[1] = 0n;
    } else {
      throw new ValueError('Invalid seed type for Philox');
    }

    // Reset counter and buffer
    this.counter.fill(0n);
    this.bufferPos = 4;
  }

  /**
   * Set the counter directly.
   */
  setCounter(counter: bigint | bigint[]): void {
    if (typeof counter === 'bigint') {
      this.counter[0] = counter;
      this.counter[1] = 0n;
      this.counter[2] = 0n;
      this.counter[3] = 0n;
    } else if (Array.isArray(counter)) {
      for (let i = 0; i < Math.min(4, counter.length); i++) {
        this.counter[i] = counter[i];
      }
    }
    this.bufferPos = 4;  // Invalidate buffer
  }

  /**
   * Generate a random 64-bit unsigned integer.
   */
  nextUint64(): bigint {
    if (this.bufferPos >= 4) {
      this._generateBlock();
    }
    return this.buffer[this.bufferPos++];
  }

  /**
   * Generate a random 32-bit unsigned integer.
   */
  nextUint32(): number {
    const val = this.nextUint64();
    return Number(val & 0xFFFFFFFFn);
  }

  /**
   * Generate a random double in [0, 1).
   */
  nextDouble(): number {
    const val = this.nextUint64();
    return Number(val >> 11n) / 9007199254740992.0;
  }

  /**
   * Generate a block of 4 random values.
   */
  private _generateBlock(): void {
    this._philox4x64(this.counter, this.key, this.buffer);
    this._incrementCounter();
    this.bufferPos = 0;
  }

  /**
   * Increment the counter.
   */
  private _incrementCounter(): void {
    this.counter[0]++;
    if (this.counter[0] === 0n) {
      this.counter[1]++;
      if (this.counter[1] === 0n) {
        this.counter[2]++;
        if (this.counter[2] === 0n) {
          this.counter[3]++;
        }
      }
    }
  }

  /**
   * Philox4x64 round function.
   */
  private _philox4x64(
    counter: BigUint64Array,
    key: BigUint64Array,
    output: BigUint64Array
  ): void {
    let v0 = counter[0];
    let v1 = counter[1];
    let v2 = counter[2];
    let v3 = counter[3];

    let k0 = key[0];
    let k1 = key[1];

    for (let round = 0; round < Philox.ROUNDS; round++) {
      // Philox round function
      const [hi0, lo0] = this._mulhilo64(v0, Philox.PHILOX_M0);
      const [hi1, lo1] = this._mulhilo64(v2, Philox.PHILOX_M1);

      v0 = hi1 ^ v1 ^ k0;
      v1 = lo1;
      v2 = hi0 ^ v3 ^ k1;
      v3 = lo0;

      // Bump key
      k0 = (k0 + Philox.PHILOX_W0) & 0xFFFFFFFFFFFFFFFFn;
      k1 = (k1 + Philox.PHILOX_W1) & 0xFFFFFFFFFFFFFFFFn;
    }

    output[0] = v0;
    output[1] = v1;
    output[2] = v2;
    output[3] = v3;
  }

  /**
   * Multiply two 64-bit values, returning high and low 64-bit parts.
   */
  private _mulhilo64(a: bigint, b: bigint): [bigint, bigint] {
    const mask = 0xFFFFFFFFn;

    const aLo = a & mask;
    const aHi = a >> 32n;
    const bLo = b & mask;
    const bHi = b >> 32n;

    const p0 = aLo * bLo;
    const p1 = aLo * bHi;
    const p2 = aHi * bLo;
    const p3 = aHi * bHi;

    const cy = ((p0 >> 32n) + (p1 & mask) + (p2 & mask)) >> 32n;

    const lo = (p0 & mask) | (((p0 >> 32n) + (p1 & mask) + (p2 & mask)) << 32n);
    const hi = p3 + (p1 >> 32n) + (p2 >> 32n) + cy;

    return [hi & 0xFFFFFFFFFFFFFFFFn, lo & 0xFFFFFFFFFFFFFFFFn];
  }

  /**
   * Jump ahead by 2^128 draws.
   */
  jump(): void {
    // Increment counter[2] (equivalent to adding 2^128 to counter)
    this.counter[2]++;
    if (this.counter[2] === 0n) {
      this.counter[3]++;
    }
    this.bufferPos = 4;  // Invalidate buffer
  }

  /**
   * Advance the generator by n draws.
   */
  advance(n: bigint): void {
    // Add n to counter (each counter increment generates 4 values)
    const blocks = n / 4n;
    const remainder = n % 4n;

    // Add blocks to counter
    const carry0 = this.counter[0];
    this.counter[0] = (this.counter[0] + blocks) & 0xFFFFFFFFFFFFFFFFn;

    if (this.counter[0] < carry0) {
      this.counter[1]++;
      if (this.counter[1] === 0n) {
        this.counter[2]++;
        if (this.counter[2] === 0n) {
          this.counter[3]++;
        }
      }
    }

    // Handle remainder by generating and discarding
    this.bufferPos = 4;  // Invalidate buffer
    for (let i = 0n; i < remainder; i++) {
      this.nextUint64();
    }
  }

  /**
   * Get the current state.
   */
  get state(): { counter: bigint[]; key: bigint[]; buffer: bigint[]; bufferPos: number } {
    return {
      counter: Array.from(this.counter),
      key: Array.from(this.key),
      buffer: Array.from(this.buffer),
      bufferPos: this.bufferPos
    };
  }

  /**
   * Set the state.
   */
  set state(value: { counter: bigint[]; key: bigint[]; buffer: bigint[]; bufferPos: number }) {
    for (let i = 0; i < 4; i++) {
      this.counter[i] = value.counter[i];
      this.buffer[i] = value.buffer[i];
    }
    this.key[0] = value.key[0];
    this.key[1] = value.key[1];
    this.bufferPos = value.bufferPos;
  }
}
```

---

### 27.3 SFC64 (Small Fast Chaotic)

#### 27.3.1 Overview

SFC64 is a fast, small-state generator with good statistical properties. It uses a 256-bit state.

**File:** `src/ts/random/sfc64.ts` (new file)

```typescript
import { BitGenerator, SeedType } from './bitgenerator.js';
import { SeedSequence } from './seedsequence.js';

/**
 * SFC64 (Small Fast Chaotic 64) pseudo-random number generator.
 *
 * SFC64 is a fast generator with a 256-bit state. It has excellent
 * statistical properties and is one of the fastest high-quality
 * generators available.
 *
 * Features:
 * - Period: approximately 2^255
 * - State size: 256 bits (4×64-bit)
 * - Very fast generation
 *
 * @example
 * const rng = new Generator(new SFC64(12345));
 * rng.random();
 */
export class SFC64 extends BitGenerator {
  /** State variables */
  private a: bigint;
  private b: bigint;
  private c: bigint;
  private w: bigint;  // Weyl sequence counter

  constructor(seed?: SeedType) {
    super();

    this.a = 0n;
    this.b = 0n;
    this.c = 0n;
    this.w = 0n;

    if (seed !== undefined) {
      this.seed(seed);
    }
  }

  /**
   * Seed the generator.
   */
  seed(seed: SeedType): void {
    if (seed instanceof SeedSequence) {
      const state = seed.generate_state(4, 'uint64');
      this.a = state[0];
      this.b = state[1];
      this.c = state[2];
      this.w = state[3];
    } else if (typeof seed === 'bigint') {
      this._seedFromScalar(seed);
    } else if (typeof seed === 'number') {
      this._seedFromScalar(BigInt(seed >>> 0));
    } else {
      throw new ValueError('Invalid seed type for SFC64');
    }

    // Run a few iterations to mix the state
    for (let i = 0; i < 12; i++) {
      this.nextUint64();
    }
  }

  /**
   * Seed from a scalar value.
   */
  private _seedFromScalar(seed: bigint): void {
    // Initialize using SplitMix64-like approach
    let x = seed;

    x = (x ^ (x >> 30n)) * 0xBF58476D1CE4E5B9n & 0xFFFFFFFFFFFFFFFFn;
    x = (x ^ (x >> 27n)) * 0x94D049BB133111EBn & 0xFFFFFFFFFFFFFFFFn;
    x = x ^ (x >> 31n);
    this.a = x & 0xFFFFFFFFFFFFFFFFn;

    x = (seed + 0x9E3779B97F4A7C15n) & 0xFFFFFFFFFFFFFFFFn;
    x = (x ^ (x >> 30n)) * 0xBF58476D1CE4E5B9n & 0xFFFFFFFFFFFFFFFFn;
    x = (x ^ (x >> 27n)) * 0x94D049BB133111EBn & 0xFFFFFFFFFFFFFFFFn;
    x = x ^ (x >> 31n);
    this.b = x & 0xFFFFFFFFFFFFFFFFn;

    x = (seed + 0x9E3779B97F4A7C15n * 2n) & 0xFFFFFFFFFFFFFFFFn;
    x = (x ^ (x >> 30n)) * 0xBF58476D1CE4E5B9n & 0xFFFFFFFFFFFFFFFFn;
    x = (x ^ (x >> 27n)) * 0x94D049BB133111EBn & 0xFFFFFFFFFFFFFFFFn;
    x = x ^ (x >> 31n);
    this.c = x & 0xFFFFFFFFFFFFFFFFn;

    this.w = 1n;  // Weyl counter
  }

  /**
   * Generate a random 64-bit unsigned integer.
   */
  nextUint64(): bigint {
    const MASK64 = 0xFFFFFFFFFFFFFFFFn;

    // SFC64 algorithm
    const tmp = (this.a + this.b + this.w) & MASK64;
    this.w = (this.w + 1n) & MASK64;

    this.a = this.b ^ (this.b >> 11n);
    this.b = (this.c + (this.c << 3n)) & MASK64;
    this.c = this._rotl(this.c, 24n) + tmp;
    this.c = this.c & MASK64;

    return tmp;
  }

  /**
   * Rotate left.
   */
  private _rotl(x: bigint, k: bigint): bigint {
    const MASK64 = 0xFFFFFFFFFFFFFFFFn;
    return ((x << k) | (x >> (64n - k))) & MASK64;
  }

  /**
   * Generate a random 32-bit unsigned integer.
   */
  nextUint32(): number {
    return Number(this.nextUint64() & 0xFFFFFFFFn);
  }

  /**
   * Generate a random double in [0, 1).
   */
  nextDouble(): number {
    const val = this.nextUint64();
    return Number(val >> 11n) / 9007199254740992.0;
  }

  /**
   * Get the current state.
   */
  get state(): { a: bigint; b: bigint; c: bigint; w: bigint } {
    return {
      a: this.a,
      b: this.b,
      c: this.c,
      w: this.w
    };
  }

  /**
   * Set the state.
   */
  set state(value: { a: bigint; b: bigint; c: bigint; w: bigint }) {
    this.a = value.a;
    this.b = value.b;
    this.c = value.c;
    this.w = value.w;
  }
}
```

---

### 27.4 BitGenerator Registry

**File:** `src/ts/random/registry.ts` (new file)

```typescript
import { BitGenerator } from './bitgenerator.js';
import { PCG64 } from './pcg64.js';
import { MT19937 } from './mt19937.js';
import { Philox } from './philox.js';
import { SFC64 } from './sfc64.js';

/**
 * BitGenerator constructor type.
 */
type BitGeneratorConstructor = new (seed?: any) => BitGenerator;

/**
 * Registry of available BitGenerators.
 */
const registry = new Map<string, BitGeneratorConstructor>([
  ['pcg64', PCG64],
  ['mt19937', MT19937],
  ['philox', Philox],
  ['sfc64', SFC64],
]);

/**
 * Get a BitGenerator class by name.
 */
export function getBitGenerator(name: string): BitGeneratorConstructor {
  const normalized = name.toLowerCase();
  const BitGen = registry.get(normalized);

  if (!BitGen) {
    const available = Array.from(registry.keys()).join(', ');
    throw new ValueError(
      `Unknown BitGenerator '${name}'. Available: ${available}`
    );
  }

  return BitGen;
}

/**
 * Register a custom BitGenerator.
 */
export function registerBitGenerator(
  name: string,
  constructor: BitGeneratorConstructor
): void {
  registry.set(name.toLowerCase(), constructor);
}

/**
 * List available BitGenerators.
 */
export function listBitGenerators(): string[] {
  return Array.from(registry.keys());
}

/**
 * Default BitGenerator class (PCG64).
 */
export const DefaultBitGenerator = PCG64;
```

#### Updates to default_rng

**File:** `src/ts/random.ts` (modifications)

```typescript
import { getBitGenerator, DefaultBitGenerator } from './random/registry.js';

/**
 * Create a Generator with the specified BitGenerator.
 *
 * @param seed - Seed value (number, bigint, array, or SeedSequence)
 * @param bitGenerator - BitGenerator name or instance
 * @returns Generator instance
 *
 * @example
 * // Default (PCG64)
 * const rng = default_rng(12345);
 *
 * // Specify BitGenerator by name
 * const rng_mt = default_rng(12345, 'mt19937');
 *
 * // Specify BitGenerator instance
 * const philox = new Philox(12345);
 * const rng_philox = default_rng(philox);
 */
export function default_rng(
  seed?: SeedType | BitGenerator,
  bitGenerator: string | BitGenerator = 'pcg64'
): Generator {
  let bg: BitGenerator;

  if (seed instanceof BitGenerator) {
    bg = seed;
  } else if (bitGenerator instanceof BitGenerator) {
    bg = bitGenerator;
  } else {
    const BitGen = getBitGenerator(bitGenerator);
    bg = new BitGen(seed);
  }

  return new Generator(bg);
}
```

---

### 27.5 WASM Acceleration

**File:** `src/wasm/random/mt19937.c` (new file)

```c
#ifndef NUMJS_MT19937_H
#define NUMJS_MT19937_H

#include <stdint.h>

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#define EXPORT EMSCRIPTEN_KEEPALIVE
#else
#define EXPORT
#endif

#define MT_N 624
#define MT_M 397

/**
 * MT19937 state structure.
 */
typedef struct {
    uint32_t mt[MT_N];
    int mti;
} MT19937State;

/**
 * Initialize MT19937 from a scalar seed.
 */
EXPORT void mt19937_seed(MT19937State* state, uint32_t seed);

/**
 * Initialize MT19937 from an array.
 */
EXPORT void mt19937_seed_array(MT19937State* state,
                                const uint32_t* init_key,
                                int key_length);

/**
 * Generate a random 32-bit unsigned integer.
 */
EXPORT uint32_t mt19937_next_uint32(MT19937State* state);

/**
 * Generate a random double in [0, 1).
 */
EXPORT double mt19937_next_double(MT19937State* state);

/**
 * Fill array with random 32-bit integers.
 */
EXPORT void mt19937_fill_uint32(MT19937State* state,
                                 uint32_t* out,
                                 int n);

/**
 * Fill array with random doubles.
 */
EXPORT void mt19937_fill_double(MT19937State* state,
                                 double* out,
                                 int n);

#endif /* NUMJS_MT19937_H */
```

---

## File Changes Summary

### New Files to Create

```
src/ts/random/
├── mt19937.ts           # Mersenne Twister
├── philox.ts            # Philox counter-based RNG
├── sfc64.ts             # Small Fast Chaotic
└── registry.ts          # BitGenerator registry

src/wasm/random/
├── mt19937.c            # MT19937 WASM implementation
├── philox.c             # Philox WASM implementation
└── sfc64.c              # SFC64 WASM implementation

tests/ts/
└── bitgenerators.test.ts  # Test suite
```

### Files to Modify

```
src/ts/random.ts
├── Import new BitGenerators
├── Update default_rng
└── Export BitGenerator classes

src/ts/index.ts
├── Export MT19937
├── Export Philox
├── Export SFC64
└── Export registry functions

scripts/build-wasm.sh
├── Add random/*.c files
└── Add EXPORTED_FUNCTIONS
```

---

## Implementation Order

```
Phase 27.1: MT19937 (Day 1-3)
├── Day 1: Core algorithm
│   ├── State management
│   ├── Twist operation
│   └── Basic tests
│
├── Day 2: Seeding
│   ├── Scalar seeding
│   ├── Array seeding
│   ├── SeedSequence integration
│   └── Tests
│
└── Day 3: WASM optimization
    ├── C implementation
    ├── Fill methods
    └── Performance tests

Phase 27.2: Philox (Day 4-6)
├── Day 4: Core algorithm
│   ├── Counter management
│   ├── Round function
│   └── Basic tests
│
├── Day 5: Advanced features
│   ├── Jump method
│   ├── Advance method
│   └── Tests
│
└── Day 6: WASM optimization
    ├── C implementation
    └── Performance tests

Phase 27.3: SFC64 (Day 7-8)
├── Day 7: Core algorithm
│   ├── State update
│   ├── Seeding
│   └── Basic tests
│
└── Day 8: WASM optimization
    ├── C implementation
    └── Performance tests

Phase 27.4: Registry & Integration (Day 9)
├── Registry implementation
├── default_rng updates
└── Integration tests

Phase 27.5: Polish (Day 10)
├── Statistical tests (TestU01 comparison)
├── NumPy compatibility verification
└── Documentation
```

---

## Verification Plan

After Phase 27 completion, verify:

```bash
# Build
npm run build

# Run tests
npm test

# Phase 27 specific tests:

# MT19937
✓ MT19937 with seed 12345 matches NumPy output
✓ MT19937 state serialization/deserialization works
✓ MT19937 array seeding produces correct sequence

# Philox
✓ Philox with seed 12345 matches NumPy output
✓ Philox jump() advances correctly
✓ Philox advance(n) produces same result as n draws

# SFC64
✓ SFC64 with seed 12345 matches NumPy output
✓ SFC64 state serialization works

# Registry
✓ default_rng('mt19937') creates MT19937-backed generator
✓ default_rng('philox') creates Philox-backed generator
✓ default_rng('sfc64') creates SFC64-backed generator
✓ Unknown BitGenerator name throws error

# Statistical tests
✓ All generators pass basic uniformity tests
✓ Chi-square test for uniform distribution
✓ Correlation test for sequential values
```

Generate NumPy comparison vectors:

```python
import numpy as np
from numpy.random import Generator, MT19937, Philox, SFC64
import json

def get_samples(bg_class, seed, n=10):
    bg = bg_class(seed)
    gen = Generator(bg)
    return {
        "integers": gen.integers(0, 1000000, n).tolist(),
        "random": gen.random(n).tolist(),
        "normal": gen.standard_normal(n).tolist()
    }

tests = {
    "mt19937": get_samples(MT19937, 12345),
    "philox": get_samples(Philox, 12345),
    "sfc64": get_samples(SFC64, 12345),
}

# Test Philox jump
philox = Philox(12345)
gen1 = Generator(philox)
samples_before_jump = gen1.random(5).tolist()
philox.jumped()
samples_after_jump = gen1.random(5).tolist()
tests["philox_jump"] = {
    "before": samples_before_jump,
    "after": samples_after_jump
}

with open("tests/fixtures/bitgenerator_vectors.json", "w") as f:
    json.dump(tests, f, indent=2)
```

---

## API Compatibility Notes

### NumPy Signature Match

```typescript
// MT19937
new MT19937(seed)
mt19937.state  // { key: Uint32Array, pos: number }

// Philox
new Philox(seed, counter?)
philox.jump()
philox.advance(n)
philox.state  // { counter, key, buffer, bufferPos }

// SFC64
new SFC64(seed)
sfc64.state  // { a, b, c, w }

// default_rng
default_rng(seed?, bitGenerator?)
```

### Differences from NumPy

1. **WASM state**: State is kept in JavaScript for simplicity; WASM functions are stateless.

2. **BigInt usage**: 64-bit values use BigInt in JavaScript vs. native 64-bit in NumPy.

3. **Philox variant**: Implements Philox4×64-10 (10 rounds, 4×64-bit output).

4. **Performance**: May be slower than NumPy's C implementation but still fast due to WASM.
