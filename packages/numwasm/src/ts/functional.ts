/**
 * NumJS Functional Programming Utilities
 *
 * Re-exports from _lib/functional.js for backwards compatibility.
 */

export {
  applyAlongAxis,
  applyOverAxes,
  Vectorize,
  vectorize,
  frompyfunc,
  piecewise,
} from "./_lib/functional.js";
export type { VectorizeOptions, UfuncLike } from "./_lib/functional.js";
