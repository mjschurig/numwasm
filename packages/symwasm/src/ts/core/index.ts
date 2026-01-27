/**
 * Core symbolic expression types.
 * @module core
 */

import { NotImplementedError } from '../errors.js';

/**
 * Base class for all symbolic expressions.
 * All symbolic objects (symbols, numbers, operations) extend this class.
 */
export abstract class Expr {
  /** String representation of this expression. */
  abstract toString(): string;

  /** Substitute a symbol with another expression. */
  subs(_old: Expr, _new: Expr): Expr {
    throw new NotImplementedError('symwasm.core.Expr.subs');
  }

  /** Evaluate the expression numerically. */
  evalf(_precision?: number): number {
    throw new NotImplementedError('symwasm.core.Expr.evalf');
  }

  /** Check structural equality with another expression. */
  equals(_other: Expr): boolean {
    throw new NotImplementedError('symwasm.core.Expr.equals');
  }

  /** Return free symbols in this expression. */
  free_symbols(): Symbol[] {
    throw new NotImplementedError('symwasm.core.Expr.free_symbols');
  }
}

/**
 * A named symbolic variable.
 * Mirrors sympy.Symbol.
 */
export class Symbol extends Expr {
  readonly name: string;

  constructor(name: string, _assumptions?: Record<string, boolean>) {
    super();
    this.name = name;
    throw new NotImplementedError('symwasm.core.Symbol');
  }

  toString(): string {
    return this.name;
  }
}

/**
 * Create multiple symbols at once.
 * Mirrors sympy.symbols.
 */
export function symbols(_names: string, _assumptions?: Record<string, boolean>): Symbol[] {
  throw new NotImplementedError('symwasm.core.symbols');
}

/** Exact integer value. Mirrors sympy.Integer. */
export class Integer extends Expr {
  readonly value: number;

  constructor(value: number) {
    super();
    this.value = value;
    throw new NotImplementedError('symwasm.core.Integer');
  }

  toString(): string {
    return String(this.value);
  }
}

/** Exact rational number p/q. Mirrors sympy.Rational. */
export class Rational extends Expr {
  readonly p: number;
  readonly q: number;

  constructor(p: number, q: number) {
    super();
    this.p = p;
    this.q = q;
    throw new NotImplementedError('symwasm.core.Rational');
  }

  toString(): string {
    return this.p + "/" + this.q;
  }
}

/** Floating-point number with arbitrary precision. Mirrors sympy.Float. */
export class Float extends Expr {
  readonly value: number;

  constructor(value: number, _precision?: number) {
    super();
    this.value = value;
    throw new NotImplementedError('symwasm.core.Float');
  }

  toString(): string {
    return String(this.value);
  }
}

/** Symbolic addition. Mirrors sympy.Add. */
export class Add extends Expr {
  constructor(readonly args: Expr[]) {
    super();
    throw new NotImplementedError('symwasm.core.Add');
  }

  toString(): string {
    return this.args.map(a => a.toString()).join(' + ');
  }
}

/** Symbolic multiplication. Mirrors sympy.Mul. */
export class Mul extends Expr {
  constructor(readonly args: Expr[]) {
    super();
    throw new NotImplementedError('symwasm.core.Mul');
  }

  toString(): string {
    return this.args.map(a => a.toString()).join('*');
  }
}

/** Symbolic exponentiation. Mirrors sympy.Pow. */
export class Pow extends Expr {
  constructor(readonly base: Expr, readonly exp: Expr) {
    super();
    throw new NotImplementedError('symwasm.core.Pow');
  }

  toString(): string {
    return this.base + '**' + this.exp;
  }
}

/**
 * Singleton-like namespace for special symbolic constants.
 * Mirrors sympy.S.
 */
export const S = {
  /** The number zero. */
  get Zero(): Expr { throw new NotImplementedError('symwasm.core.S.Zero'); },
  /** The number one. */
  get One(): Expr { throw new NotImplementedError('symwasm.core.S.One'); },
  /** Negative one. */
  get NegativeOne(): Expr { throw new NotImplementedError('symwasm.core.S.NegativeOne'); },
  /** One half. */
  get Half(): Expr { throw new NotImplementedError('symwasm.core.S.Half'); },
  /** Positive infinity. */
  get Infinity(): Expr { throw new NotImplementedError('symwasm.core.S.Infinity'); },
  /** Negative infinity. */
  get NegativeInfinity(): Expr { throw new NotImplementedError('symwasm.core.S.NegativeInfinity'); },
};

/** Positive infinity. */
export const oo: Expr = Object.freeze({ toString: () => 'oo' }) as unknown as Expr;

/** Pi (π). */
export const pi: Expr = Object.freeze({ toString: () => 'pi' }) as unknown as Expr;

/** Euler's number (e). */
export const E: Expr = Object.freeze({ toString: () => 'E' }) as unknown as Expr;

/** Imaginary unit (i). */
export const I: Expr = Object.freeze({ toString: () => 'I' }) as unknown as Expr;
