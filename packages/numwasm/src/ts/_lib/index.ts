/**
 * NumJS Library Utilities Module
 *
 * Internal utilities for functional programming, iteration, type checking, and type info.
 */

// Iterators
export { FlatIterator, nditer, ndenumerate, ndindex } from "./iterators.js";

// Functional programming utilities
export {
  applyAlongAxis,
  applyOverAxes,
  Vectorize,
  vectorize,
  frompyfunc,
  piecewise,
} from "./functional.js";
export type { VectorizeOptions, UfuncLike } from "./functional.js";

// Type checking and complex utilities
export { angle, real, imag } from "./type_check.js";

// Type information classes
export { finfo, iinfo } from "./typeinfo.js";
