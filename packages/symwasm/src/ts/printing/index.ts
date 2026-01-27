/**
 * Expression formatting and printing.
 * @module printing
 */

import type { Expr } from '../core/index.js';
import { NotImplementedError } from '../errors.js';

/**
 * Convert an expression to LaTeX string.
 * Mirrors sympy.latex.
 */
export function latex(_expr: Expr, _options?: { mode?: 'plain' | 'inline' | 'equation' }): string {
  throw new NotImplementedError('symwasm.printing.latex');
}

/**
 * Convert an expression to MathML string.
 * Mirrors sympy.mathml.
 */
export function mathml(_expr: Expr, _printer?: 'content' | 'presentation'): string {
  throw new NotImplementedError('symwasm.printing.mathml');
}

/**
 * Pretty-print an expression (Unicode art).
 * Mirrors sympy.pretty.
 */
export function pretty(_expr: Expr, _options?: { use_unicode?: boolean }): string {
  throw new NotImplementedError('symwasm.printing.pretty');
}

/**
 * Simple string representation.
 * Mirrors sympy.sstr.
 */
export function sstr(_expr: Expr): string {
  throw new NotImplementedError('symwasm.printing.sstr');
}
