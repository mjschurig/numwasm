/**
 * Bessel functions via WebAssembly.
 * Mirrors scipy.special Bessel function implementations.
 */

import { loadXSFModule } from 'xsfwasm';

// ============================================================
// Bessel functions of the first kind (J)
// ============================================================

/**
 * Bessel function of the first kind, order 0: J₀(x)
 *
 * J₀(x) is a solution to Bessel's differential equation and appears
 * in many physical problems involving cylindrical symmetry.
 *
 * Properties:
 *   - J₀(0) = 1
 *   - J₀(x) oscillates with decreasing amplitude as x → ∞
 *
 * @param x - Input value(s). Can be a number or array of numbers.
 * @returns J₀(x) value(s).
 *
 * @example
 * ```ts
 * import { special } from 'sciwasm';
 *
 * // J₀(0) = 1
 * console.log(await special.j0(0)); // 1
 *
 * // First zero of J₀ is at x ≈ 2.4048
 * console.log(await special.j0(2.4048)); // ≈ 0
 * ```
 */
export async function j0(x: number | number[]): Promise<number | number[]> {
  const Module = await loadXSFModule();

  if (Array.isArray(x)) {
    return x.map((val) => Module._wasm_j0(val));
  }

  return Module._wasm_j0(x);
}

/**
 * Bessel function of the first kind, order 1: J₁(x)
 *
 * Properties:
 *   - J₁(0) = 0
 *   - J₁(x) oscillates with decreasing amplitude as x → ∞
 *
 * @param x - Input value(s). Can be a number or array of numbers.
 * @returns J₁(x) value(s).
 *
 * @example
 * ```ts
 * import { special } from 'sciwasm';
 *
 * // J₁(0) = 0
 * console.log(await special.j1(0)); // 0
 *
 * // Maximum of J₁ is at x ≈ 1.84
 * console.log(await special.j1(1.84)); // ≈ 0.582
 * ```
 */
export async function j1(x: number | number[]): Promise<number | number[]> {
  const Module = await loadXSFModule();

  if (Array.isArray(x)) {
    return x.map((val) => Module._wasm_j1(val));
  }

  return Module._wasm_j1(x);
}

/**
 * Bessel function of the first kind, order v: Jᵥ(x)
 *
 * General Bessel function of the first kind for arbitrary order v.
 *
 * @param v - Order of the Bessel function.
 * @param x - Input value(s). Can be a number or array of numbers.
 * @returns Jᵥ(x) value(s).
 *
 * @example
 * ```ts
 * import { special } from 'sciwasm';
 *
 * // J₂(1) ≈ 0.1149
 * console.log(await special.jv(2, 1)); // 0.1149...
 *
 * // Half-integer order
 * console.log(await special.jv(0.5, 1)); // ≈ 0.6714
 * ```
 */
export async function jv(
  v: number,
  x: number | number[]
): Promise<number | number[]> {
  const Module = await loadXSFModule();

  if (Array.isArray(x)) {
    return x.map((val) => Module._wasm_jv(v, val));
  }

  return Module._wasm_jv(v, x);
}

// ============================================================
// Bessel functions of the second kind (Y)
// ============================================================

/**
 * Bessel function of the second kind, order 0: Y₀(x)
 *
 * Also known as the Neumann function or Weber function.
 * Y₀(x) is singular at x = 0 and is real for x > 0.
 *
 * @param x - Input value(s). Must be positive. Can be a number or array of numbers.
 * @returns Y₀(x) value(s).
 *
 * @example
 * ```ts
 * import { special } from 'sciwasm';
 *
 * // Y₀(1) ≈ 0.0883
 * console.log(await special.y0(1)); // 0.0883...
 *
 * // Y₀ → -∞ as x → 0
 * console.log(await special.y0(0.01)); // large negative
 * ```
 */
export async function y0(x: number | number[]): Promise<number | number[]> {
  const Module = await loadXSFModule();

  if (Array.isArray(x)) {
    return x.map((val) => Module._wasm_y0(val));
  }

  return Module._wasm_y0(x);
}

/**
 * Bessel function of the second kind, order 1: Y₁(x)
 *
 * @param x - Input value(s). Must be positive. Can be a number or array of numbers.
 * @returns Y₁(x) value(s).
 *
 * @example
 * ```ts
 * import { special } from 'sciwasm';
 *
 * // Y₁(1) ≈ -0.7812
 * console.log(await special.y1(1)); // -0.7812...
 * ```
 */
export async function y1(x: number | number[]): Promise<number | number[]> {
  const Module = await loadXSFModule();

  if (Array.isArray(x)) {
    return x.map((val) => Module._wasm_y1(val));
  }

  return Module._wasm_y1(x);
}

/**
 * Bessel function of the second kind, order v: Yᵥ(x)
 *
 * General Bessel function of the second kind for arbitrary order v.
 *
 * @param v - Order of the Bessel function.
 * @param x - Input value(s). Must be positive. Can be a number or array of numbers.
 * @returns Yᵥ(x) value(s).
 *
 * @example
 * ```ts
 * import { special } from 'sciwasm';
 *
 * // Y₂(1) ≈ -1.6507
 * console.log(await special.yv(2, 1)); // -1.6507...
 * ```
 */
export async function yv(
  v: number,
  x: number | number[]
): Promise<number | number[]> {
  const Module = await loadXSFModule();

  if (Array.isArray(x)) {
    return x.map((val) => Module._wasm_yv(v, val));
  }

  return Module._wasm_yv(v, x);
}

// ============================================================
// Modified Bessel functions of the first kind (I)
// ============================================================

/**
 * Modified Bessel function of the first kind, order 0: I₀(x)
 *
 * I₀(x) is related to J₀(ix) and grows exponentially for large x.
 *
 * Properties:
 *   - I₀(0) = 1
 *   - I₀(x) = I₀(-x) (even function)
 *   - I₀(x) ∼ e^x / √(2πx) as x → ∞
 *
 * @param x - Input value(s). Can be a number or array of numbers.
 * @returns I₀(x) value(s).
 *
 * @example
 * ```ts
 * import { special } from 'sciwasm';
 *
 * // I₀(0) = 1
 * console.log(await special.i0(0)); // 1
 *
 * // I₀(1) ≈ 1.2661
 * console.log(await special.i0(1)); // 1.2661...
 * ```
 */
export async function i0(x: number | number[]): Promise<number | number[]> {
  const Module = await loadXSFModule();

  if (Array.isArray(x)) {
    return x.map((val) => Module._wasm_i0(val));
  }

  return Module._wasm_i0(x);
}

/**
 * Modified Bessel function of the first kind, order 1: I₁(x)
 *
 * Properties:
 *   - I₁(0) = 0
 *   - I₁(-x) = -I₁(x) (odd function)
 *
 * @param x - Input value(s). Can be a number or array of numbers.
 * @returns I₁(x) value(s).
 *
 * @example
 * ```ts
 * import { special } from 'sciwasm';
 *
 * // I₁(0) = 0
 * console.log(await special.i1(0)); // 0
 *
 * // I₁(1) ≈ 0.5652
 * console.log(await special.i1(1)); // 0.5652...
 * ```
 */
export async function i1(x: number | number[]): Promise<number | number[]> {
  const Module = await loadXSFModule();

  if (Array.isArray(x)) {
    return x.map((val) => Module._wasm_i1(val));
  }

  return Module._wasm_i1(x);
}

/**
 * Modified Bessel function of the first kind, order v: Iᵥ(x)
 *
 * General modified Bessel function of the first kind for arbitrary order v.
 *
 * @param v - Order of the Bessel function.
 * @param x - Input value(s). Can be a number or array of numbers.
 * @returns Iᵥ(x) value(s).
 *
 * @example
 * ```ts
 * import { special } from 'sciwasm';
 *
 * // I₂(1) ≈ 0.1357
 * console.log(await special.iv(2, 1)); // 0.1357...
 * ```
 */
export async function iv(
  v: number,
  x: number | number[]
): Promise<number | number[]> {
  const Module = await loadXSFModule();

  if (Array.isArray(x)) {
    return x.map((val) => Module._wasm_iv(v, val));
  }

  return Module._wasm_iv(v, x);
}

// ============================================================
// Modified Bessel functions of the second kind (K)
// ============================================================

/**
 * Modified Bessel function of the second kind, order 0: K₀(x)
 *
 * K₀(x) decays exponentially for large x and is singular at x = 0.
 *
 * Properties:
 *   - K₀(x) → ∞ as x → 0
 *   - K₀(x) ∼ √(π/(2x)) * e^(-x) as x → ∞
 *
 * @param x - Input value(s). Must be positive. Can be a number or array of numbers.
 * @returns K₀(x) value(s).
 *
 * @example
 * ```ts
 * import { special } from 'sciwasm';
 *
 * // K₀(1) ≈ 0.4210
 * console.log(await special.k0(1)); // 0.4210...
 *
 * // K₀ decays exponentially
 * console.log(await special.k0(5)); // ≈ 0.00369
 * ```
 */
export async function k0(x: number | number[]): Promise<number | number[]> {
  const Module = await loadXSFModule();

  if (Array.isArray(x)) {
    return x.map((val) => Module._wasm_k0(val));
  }

  return Module._wasm_k0(x);
}

/**
 * Modified Bessel function of the second kind, order 1: K₁(x)
 *
 * @param x - Input value(s). Must be positive. Can be a number or array of numbers.
 * @returns K₁(x) value(s).
 *
 * @example
 * ```ts
 * import { special } from 'sciwasm';
 *
 * // K₁(1) ≈ 0.6019
 * console.log(await special.k1(1)); // 0.6019...
 * ```
 */
export async function k1(x: number | number[]): Promise<number | number[]> {
  const Module = await loadXSFModule();

  if (Array.isArray(x)) {
    return x.map((val) => Module._wasm_k1(val));
  }

  return Module._wasm_k1(x);
}
