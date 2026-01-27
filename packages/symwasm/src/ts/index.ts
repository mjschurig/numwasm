/**
 * symwasm - SymPy-inspired symbolic mathematics in TypeScript.
 *
 * Provides symbolic computation capabilities including symbolic expressions,
 * algebraic simplification, equation solving, calculus operations,
 * matrix algebra, and expression formatting/printing.
 */

export * as core from './core/index.js';
export * as simplify from './simplify/index.js';
export * as solvers from './solvers/index.js';
export * as calculus from './calculus/index.js';
export * as matrices from './matrices/index.js';
export * as printing from './printing/index.js';
export { NotImplementedError } from './errors.js';
