/**
 * Equation solvers.
 * @module solvers
 */

import type { Expr, Symbol } from '../core/index.js';
import { NotImplementedError } from '../errors.js';

/**
 * Solve algebraic equations.
 * Mirrors sympy.solve.
 */
export function solve(
  _expr: Expr | Expr[],
  _symbols?: Symbol | Symbol[],
  _options?: { dict?: boolean }
): Expr[] | Record<string, Expr>[] {
  throw new NotImplementedError('symwasm.solvers.solve');
}

/**
 * Solve equation returning a set.
 * Mirrors sympy.solveset.
 */
export function solveset(
  _expr: Expr,
  _symbol: Symbol,
  _domain?: 'Reals' | 'Complexes'
): Expr[] {
  throw new NotImplementedError('symwasm.solvers.solveset');
}

/**
 * Solve a system of linear equations.
 * Mirrors sympy.linsolve.
 */
export function linsolve(
  _system: Expr[],
  _symbols: Symbol[]
): Expr[][] {
  throw new NotImplementedError('symwasm.solvers.linsolve');
}

/**
 * Solve a system of nonlinear equations.
 * Mirrors sympy.nonlinsolve.
 */
export function nonlinsolve(
  _system: Expr[],
  _symbols: Symbol[]
): Expr[][] {
  throw new NotImplementedError('symwasm.solvers.nonlinsolve');
}

/**
 * Solve ordinary differential equations.
 * Mirrors sympy.dsolve.
 */
export function dsolve(
  _eq: Expr,
  _func?: Expr
): Expr {
  throw new NotImplementedError('symwasm.solvers.dsolve');
}
