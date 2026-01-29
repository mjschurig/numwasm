/**
 * Tests for symwasm arithmetic operations (Phase 1.5)
 *
 * Tests add, sub, mul, div, pow, neg helper functions
 * and Add, Mul, Pow expression classes
 */
import { describe, it, expect, beforeAll } from 'vitest';
import * as sym from '../../src/ts/index';

describe('arithmetic: add()', () => {
  beforeAll(async () => {
    await sym.core.loadWasmModule();
  });

  it('adds two integers', () => {
    const a = new sym.core.Integer(2);
    const b = new sym.core.Integer(3);
    const result = sym.core.add(a, b);
    expect(result.toString()).toBe('5');
  });

  it('adds integer and symbol', () => {
    const x = new sym.core.Symbol('x');
    const one = new sym.core.Integer(1);
    const result = sym.core.add(x, one);
    expect(result.toString()).toBe('1 + x');
  });

  it('simplifies x + x to 2*x', () => {
    const x = new sym.core.Symbol('x');
    const result = sym.core.add(x, x);
    expect(result.toString()).toBe('2*x');
  });

  it('adds two symbols', () => {
    const x = new sym.core.Symbol('x');
    const y = new sym.core.Symbol('y');
    const result = sym.core.add(x, y);
    expect(result.toString()).toBe('x + y');
  });

  it('returns correct type for symbolic result', () => {
    const x = new sym.core.Symbol('x');
    const y = new sym.core.Symbol('y');
    const result = sym.core.add(x, y);
    expect(result.get_type()).toBe(sym.core.SymEngineTypeID.SYMENGINE_ADD);
  });

  it('returns Integer type for numeric result', () => {
    const a = new sym.core.Integer(2);
    const b = new sym.core.Integer(3);
    const result = sym.core.add(a, b);
    expect(result.get_type()).toBe(sym.core.SymEngineTypeID.SYMENGINE_INTEGER);
  });
});

describe('arithmetic: sub()', () => {
  beforeAll(async () => {
    await sym.core.loadWasmModule();
  });

  it('subtracts two integers', () => {
    const a = new sym.core.Integer(5);
    const b = new sym.core.Integer(3);
    const result = sym.core.sub(a, b);
    expect(result.toString()).toBe('2');
  });

  it('subtracts symbol and integer', () => {
    const x = new sym.core.Symbol('x');
    const one = new sym.core.Integer(1);
    const result = sym.core.sub(x, one);
    expect(result.toString()).toBe('-1 + x');
  });

  it('simplifies x - x to 0', () => {
    const x = new sym.core.Symbol('x');
    const result = sym.core.sub(x, x);
    expect(result.toString()).toBe('0');
  });
});

describe('arithmetic: mul()', () => {
  beforeAll(async () => {
    await sym.core.loadWasmModule();
  });

  it('multiplies two integers', () => {
    const a = new sym.core.Integer(2);
    const b = new sym.core.Integer(3);
    const result = sym.core.mul(a, b);
    expect(result.toString()).toBe('6');
  });

  it('multiplies symbol and integer', () => {
    const x = new sym.core.Symbol('x');
    const two = new sym.core.Integer(2);
    const result = sym.core.mul(x, two);
    expect(result.toString()).toBe('2*x');
  });

  it('simplifies x * x to x**2', () => {
    const x = new sym.core.Symbol('x');
    const result = sym.core.mul(x, x);
    expect(result.toString()).toBe('x**2');
  });

  it('multiplies two symbols', () => {
    const x = new sym.core.Symbol('x');
    const y = new sym.core.Symbol('y');
    const result = sym.core.mul(x, y);
    expect(result.toString()).toBe('x*y');
  });

  it('returns correct type for symbolic result', () => {
    const x = new sym.core.Symbol('x');
    const y = new sym.core.Symbol('y');
    const result = sym.core.mul(x, y);
    expect(result.get_type()).toBe(sym.core.SymEngineTypeID.SYMENGINE_MUL);
  });

  it('multiplying by zero returns zero', () => {
    const x = new sym.core.Symbol('x');
    const zero = new sym.core.Integer(0);
    const result = sym.core.mul(x, zero);
    expect(result.toString()).toBe('0');
  });
});

describe('arithmetic: div()', () => {
  beforeAll(async () => {
    await sym.core.loadWasmModule();
  });

  it('divides two integers (exact)', () => {
    const a = new sym.core.Integer(6);
    const b = new sym.core.Integer(2);
    const result = sym.core.div(a, b);
    expect(result.toString()).toBe('3');
  });

  it('divides two integers (rational result)', () => {
    const a = new sym.core.Integer(1);
    const b = new sym.core.Integer(2);
    const result = sym.core.div(a, b);
    expect(result.toString()).toBe('1/2');
  });

  it('divides symbol by integer', () => {
    const x = new sym.core.Symbol('x');
    const two = new sym.core.Integer(2);
    const result = sym.core.div(x, two);
    expect(result.toString()).toBe('(1/2)*x');
  });

  it('simplifies x / x to 1', () => {
    const x = new sym.core.Symbol('x');
    const result = sym.core.div(x, x);
    expect(result.toString()).toBe('1');
  });
});

describe('arithmetic: pow()', () => {
  beforeAll(async () => {
    await sym.core.loadWasmModule();
  });

  it('raises integer to integer power', () => {
    const a = new sym.core.Integer(2);
    const b = new sym.core.Integer(3);
    const result = sym.core.pow(a, b);
    expect(result.toString()).toBe('8');
  });

  it('raises symbol to integer power', () => {
    const x = new sym.core.Symbol('x');
    const two = new sym.core.Integer(2);
    const result = sym.core.pow(x, two);
    expect(result.toString()).toBe('x**2');
  });

  it('returns correct type for symbolic power', () => {
    const x = new sym.core.Symbol('x');
    const two = new sym.core.Integer(2);
    const result = sym.core.pow(x, two);
    expect(result.get_type()).toBe(sym.core.SymEngineTypeID.SYMENGINE_POW);
  });

  it('x**0 equals 1', () => {
    const x = new sym.core.Symbol('x');
    const zero = new sym.core.Integer(0);
    const result = sym.core.pow(x, zero);
    expect(result.toString()).toBe('1');
  });

  it('x**1 equals x', () => {
    const x = new sym.core.Symbol('x');
    const one = new sym.core.Integer(1);
    const result = sym.core.pow(x, one);
    expect(result.toString()).toBe('x');
  });
});

describe('arithmetic: neg()', () => {
  beforeAll(async () => {
    await sym.core.loadWasmModule();
  });

  it('negates positive integer', () => {
    const a = new sym.core.Integer(5);
    const result = sym.core.neg(a);
    expect(result.toString()).toBe('-5');
  });

  it('negates negative integer', () => {
    const a = new sym.core.Integer(-3);
    const result = sym.core.neg(a);
    expect(result.toString()).toBe('3');
  });

  it('negates symbol', () => {
    const x = new sym.core.Symbol('x');
    const result = sym.core.neg(x);
    expect(result.toString()).toBe('-x');
  });

  it('double negation returns original', () => {
    const x = new sym.core.Symbol('x');
    const negX = sym.core.neg(x);
    const result = sym.core.neg(negX);
    expect(result.toString()).toBe('x');
  });
});

describe('arithmetic: Add class', () => {
  beforeAll(async () => {
    await sym.core.loadWasmModule();
  });

  it('can extract args from Add expression', () => {
    const x = new sym.core.Symbol('x');
    const y = new sym.core.Symbol('y');
    const result = sym.core.add(x, y);

    // Check it's an Add
    expect(result.get_type()).toBe(sym.core.SymEngineTypeID.SYMENGINE_ADD);

    // Get args
    const args = result.get_args();
    expect(args.length).toBe(2);

    // Args should be x and y (order may vary)
    const argStrings = args.map((a) => a.toString()).sort();
    expect(argStrings).toEqual(['x', 'y']);
  });
});

describe('arithmetic: Mul class', () => {
  beforeAll(async () => {
    await sym.core.loadWasmModule();
  });

  it('can extract args from Mul expression', () => {
    const x = new sym.core.Symbol('x');
    const y = new sym.core.Symbol('y');
    const result = sym.core.mul(x, y);

    // Check it's a Mul
    expect(result.get_type()).toBe(sym.core.SymEngineTypeID.SYMENGINE_MUL);

    // Get args
    const args = result.get_args();
    expect(args.length).toBe(2);

    // Args should be x and y (order may vary)
    const argStrings = args.map((a) => a.toString()).sort();
    expect(argStrings).toEqual(['x', 'y']);
  });
});

describe('arithmetic: Pow class', () => {
  beforeAll(async () => {
    await sym.core.loadWasmModule();
  });

  it('can extract args from Pow expression', () => {
    const x = new sym.core.Symbol('x');
    const two = new sym.core.Integer(2);
    const result = sym.core.pow(x, two);

    // Check it's a Pow
    expect(result.get_type()).toBe(sym.core.SymEngineTypeID.SYMENGINE_POW);

    // Get args - should be [base, exponent]
    const args = result.get_args();
    expect(args.length).toBe(2);
    expect(args[0].toString()).toBe('x');
    expect(args[1].toString()).toBe('2');
  });
});

describe('arithmetic: complex expressions', () => {
  beforeAll(async () => {
    await sym.core.loadWasmModule();
  });

  it('builds x^2 + 2*x + 1', () => {
    const x = new sym.core.Symbol('x');
    const one = new sym.core.Integer(1);
    const two = new sym.core.Integer(2);

    // x^2
    const x2 = sym.core.pow(x, two);
    // 2*x
    const twoX = sym.core.mul(two, x);
    // x^2 + 2*x
    const partial = sym.core.add(x2, twoX);
    // x^2 + 2*x + 1
    const result = sym.core.add(partial, one);

    expect(result.toString()).toBe('1 + 2*x + x**2');
  });

  it('free_symbols returns symbols in expression', () => {
    const x = new sym.core.Symbol('x');
    const y = new sym.core.Symbol('y');
    const result = sym.core.add(x, y);

    const freeSyms = result.free_symbols();
    expect(freeSyms.length).toBe(2);

    const symNames = freeSyms.map((s) => s.toString()).sort();
    expect(symNames).toEqual(['x', 'y']);
  });
});
