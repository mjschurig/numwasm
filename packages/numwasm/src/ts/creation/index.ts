/**
 * Array Creation Module
 *
 * Functions for creating NDArrays with various patterns and values.
 */

// Basic array creation (zeros, ones, empty, full, array)
export {
  zeros,
  ones,
  empty,
  full,
  array,
  zeros_like,
  ones_like,
  empty_like,
  full_like,
  copy,
  asanyarray,
  asarray_chkfinite,
  ascontiguousarray,
  asfortranarray,
} from "./basic.js";

// Range-based creation (arange, linspace, logspace, geomspace)
export {
  arange,
  linspace,
  logspace,
  geomspace,
} from "./ranges.js";

// Matrix construction (eye, identity, diag, tri, tril, triu, vander)
export {
  eye,
  identity,
  diag,
  tri,
  tril,
  triu,
  diagflat,
  vander,
} from "./matrices.js";

// Data-based creation (fromfunction, fromiter)
export {
  fromfunction,
  fromiter,
} from "./from_data.js";
