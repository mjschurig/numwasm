/**
 * NumJS String Operations Module
 *
 * Provides vectorized string operations on arrays of strings,
 * compatible with NumPy's numpy.strings module.
 *
 * @example
 * ```typescript
 * import { strings, NDArray } from 'numjs';
 *
 * // Create a string array
 * const arr = NDArray.fromStringArray(['hello', 'world']);
 *
 * // Use string operations
 * const upper = strings.upper(arr);
 * const found = strings.find(arr, 'o');
 * ```
 */

// Error classes
export { ValueError } from './errors.js';

// Comparison functions
export {
  equal,
  not_equal,
  less,
  less_equal,
  greater,
  greater_equal,
  compare_chararrays,
} from './compare.js';

// Property testing
export {
  isalpha,
  isdigit,
  isalnum,
  isspace,
  islower,
  isupper,
  istitle,
  isdecimal,
  isnumeric,
  str_len,
} from './properties.js';

// Search and indexing
export {
  find,
  rfind,
  index,
  rindex,
  count,
  startswith,
  endswith,
} from './search.js';

// String manipulation
export {
  // Case conversion
  lower,
  upper,
  swapcase,
  capitalize,
  title,
  // Concatenation
  add,
  multiply,
  // Whitespace
  strip,
  lstrip,
  rstrip,
  expandtabs,
  // Replacement
  replace,
  // Alignment
  center,
  ljust,
  rjust,
  zfill,
  // Partitioning
  partition,
  rpartition,
  // Encoding
  encode,
  decode,
  // Formatting (NumPy 2.0)
  mod,
  translate,
  slice,
} from './manipulation.js';

// Re-import for namespace object
import {
  equal,
  not_equal,
  less,
  less_equal,
  greater,
  greater_equal,
  compare_chararrays,
} from './compare.js';

import {
  isalpha,
  isdigit,
  isalnum,
  isspace,
  islower,
  isupper,
  istitle,
  isdecimal,
  isnumeric,
  str_len,
} from './properties.js';

import {
  find,
  rfind,
  index,
  rindex,
  count,
  startswith,
  endswith,
} from './search.js';

import {
  lower,
  upper,
  swapcase,
  capitalize,
  title,
  add,
  multiply,
  strip,
  lstrip,
  rstrip,
  expandtabs,
  replace,
  center,
  ljust,
  rjust,
  zfill,
  partition,
  rpartition,
  encode,
  decode,
  mod,
  translate,
  slice,
} from './manipulation.js';

/**
 * Strings module namespace object for convenient grouped import.
 *
 * @example
 * ```typescript
 * import { strings } from 'numjs';
 *
 * strings.upper(['hello', 'world']);
 * strings.find(['hello'], 'l');
 * ```
 */
export const strings = {
  // Comparison
  equal,
  not_equal,
  less,
  less_equal,
  greater,
  greater_equal,
  compare_chararrays,

  // Property testing
  isalpha,
  isdigit,
  isalnum,
  isspace,
  islower,
  isupper,
  istitle,
  isdecimal,
  isnumeric,
  str_len,

  // Search
  find,
  rfind,
  index,
  rindex,
  count,
  startswith,
  endswith,

  // Case conversion
  lower,
  upper,
  swapcase,
  capitalize,
  title,

  // Concatenation
  add,
  multiply,

  // Whitespace
  strip,
  lstrip,
  rstrip,
  expandtabs,

  // Replacement
  replace,

  // Alignment
  center,
  ljust,
  rjust,
  zfill,

  // Partitioning
  partition,
  rpartition,

  // Encoding
  encode,
  decode,

  // Formatting (NumPy 2.0)
  mod,
  translate,
  slice,
};
