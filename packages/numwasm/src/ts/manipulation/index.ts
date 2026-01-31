/**
 * Array manipulation functions for NumJS-WASM
 *
 * Re-exports from submodules for convenience.
 */

// Shape manipulation
export { reshape, ravel, flatten, squeeze, expand_dims } from "./shape.js";

// Transpose and axis manipulation
export {
  transpose,
  permute_dims,
  swapaxes,
  moveaxis,
  rollaxis,
} from "./transpose.js";

// Joining arrays
export {
  concatenate,
  concat,
  stack,
  vstack,
  row_stack,
  hstack,
  dstack,
  column_stack,
  block,
  append,
} from "./join.js";

// Splitting arrays
export {
  array_split,
  split,
  vsplit,
  hsplit,
  dsplit,
  unstack,
} from "./split.js";

// Tiling and repetition
export { tile, repeat, pad } from "./tile.js";

// Rearranging arrays
export {
  flip,
  fliplr,
  flipud,
  roll,
  rot90,
  resize,
  trim_zeros,
} from "./rearrange.js";

// Insert and delete
export { insert, deleteArr, copyto, asarray } from "./insert_delete.js";
