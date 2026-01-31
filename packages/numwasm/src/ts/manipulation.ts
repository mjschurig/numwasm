/**
 * Array manipulation functions for NumJS-WASM
 *
 * Re-exports from manipulation/ for backwards compatibility.
 */

// Shape manipulation
export { reshape, ravel, flatten, squeeze, expand_dims } from "./manipulation/shape.js";

// Transpose and axis manipulation
export {
  transpose,
  permute_dims,
  swapaxes,
  moveaxis,
  rollaxis,
} from "./manipulation/transpose.js";

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
} from "./manipulation/join.js";

// Splitting arrays
export {
  array_split,
  split,
  vsplit,
  hsplit,
  dsplit,
  unstack,
} from "./manipulation/split.js";

// Tiling and repetition
export { tile, repeat, pad } from "./manipulation/tile.js";

// Rearranging arrays
export {
  flip,
  fliplr,
  flipud,
  roll,
  rot90,
  resize,
  trim_zeros,
} from "./manipulation/rearrange.js";

// Insert and delete
export { insert, deleteArr, copyto, asarray } from "./manipulation/insert_delete.js";
