/**
 * NumJS - NumPy-inspired array operations for TypeScript/WebAssembly
 *
 * This library provides efficient n-dimensional array operations using
 * WebAssembly for performance-critical computations.
 *
 * @packageDocumentation
 */

// Main array class
export { NDArray } from './NDArray.js';

// Type definitions
export {
  DType,
  DTYPE_SIZES,
  DTYPE_NAMES,
  NDARRAY_OWNDATA,
  NDARRAY_WRITEABLE,
  NDARRAY_C_CONTIGUOUS,
  NDARRAY_F_CONTIGUOUS,
  NDARRAY_ALIGNED,
  CLIP_RAISE,
  CLIP_WRAP,
  CLIP_CLIP,
  INDEX_TYPE_INTEGER,
  INDEX_TYPE_SLICE,
  INDEX_TYPE_NEWAXIS,
  INDEX_TYPE_ELLIPSIS,
} from './types.js';
export type {
  NDArrayOptions,
  TypedArrayType,
  WasmModule,
  Complex,
  NDArrayFlags,
} from './types.js';

// DType utilities
export {
  dtypeFromString,
  dtypeToString,
  dtypeSize,
  isIntegerDType,
  isFloatDType,
  isComplexDType,
  isSignedDType,
  isBoolDType,
  isNumericDType,
  dtypeToTypedArrayConstructor,
  typedArrayToDType,
  promoteTypes,
  commonType,
  CastingKind,
} from './dtype.js';

// Iterators
export { FlatIterator, nditer, ndenumerate, ndindex } from './iterators.js';

// WASM module management
export { loadWasmModule, isWasmLoaded, getWasmModule } from './wasm-loader.js';

// Slice utilities (Level 2)
export {
  Slice,
  slice,
  ellipsis,
  newaxis,
  expandEllipsis,
  buildIndexSpecs,
  computeResultShape,
} from './slice.js';
export type { IndexElement, IndexSpec, Ellipsis, Newaxis } from './slice.js';

// Broadcasting functions (Level 2)
export {
  broadcastShapes,
  broadcastShapesMulti,
  shapesAreBroadcastable,
  broadcastTo,
  broadcastArrays,
  computeBroadcastStrides,
} from './broadcast.js';

// Index functions (Level 2)
export {
  take,
  takeFlat,
  put,
  countNonzero,
  nonzero,
  flatnonzero,
  where,
  compress,
  extract,
  choose,
  diagonal,
  unravelIndex,
  ravelMultiIndex,
  meshgrid,
  // Shape manipulation
  atleast_1d,
  atleast_2d,
  atleast_3d,
  // Index generation
  argwhere,
  indices,
  ix_,
  diag_indices,
  tril_indices,
  triu_indices,
  // Advanced indexing
  take_along_axis,
  put_along_axis,
  putmask,
  place,
  select,
} from './indexing.js';
export type { ClipMode } from './indexing.js';

// Sorting functions (Phase 6)
export { sort, argsort, partition, argpartition } from './sorting.js';
export type { SortKind } from './sorting.js';

// Statistics functions (Phase 6)
export {
  sum,
  mean,
  variance,
  var_,
  std,
  min,
  max,
  median,
  argmax,
  argmin,
  searchsorted,
} from './statistics.js';

// Manipulation functions (Phase 5)
export {
  // Joining
  concatenate,
  stack,
  vstack,
  row_stack,
  hstack,
  dstack,
  column_stack,
  block,
  append,
  // Splitting
  split,
  array_split,
  vsplit,
  hsplit,
  dsplit,
  unstack,
  // Tiling
  tile,
  repeat,
  pad,
  // Rearranging
  flip,
  fliplr,
  flipud,
  roll,
  rot90,
  resize,
  trim_zeros,
  // Insert/Delete
  insert,
  deleteArr,
  // Copying
  copyto,
  asarray,
} from './manipulation.js';

// Logic & comparison functions (Phase 7)
export {
  // Element-wise predicates
  isfinite,
  isinf,
  isnan,
  isneginf,
  isposinf,
  iscomplex,
  isreal,
  // Type checking (TypeScript only)
  iscomplexobj,
  isrealobj,
  isfortran,
  isscalar,
  // Truth testing reductions
  all,
  any,
  // Comparison functions
  isclose,
  allclose,
  array_equal,
  array_equiv,
} from './logic.js';

// Universal functions (Ufuncs) - Level 3
export {
  // Arithmetic
  negative,
  positive,
  absolute,
  abs,
  sign,
  // Powers and roots
  sqrt,
  square,
  cbrt,
  reciprocal,
  // Exponential
  exp,
  exp2,
  expm1,
  // Logarithmic
  log,
  log2,
  log10,
  log1p,
  // Trigonometric
  sin,
  cos,
  tan,
  arcsin,
  arccos,
  arctan,
  // Hyperbolic
  sinh,
  cosh,
  tanh,
  arcsinh,
  arccosh,
  arctanh,
  // Rounding
  floor,
  ceil,
  trunc,
  rint,
  round,
  // Angle conversion
  degrees,
  rad2deg,
  radians,
  deg2rad,
  // Logical (unary)
  logical_not,
  // Bitwise (unary)
  invert,
  bitwise_not,
  // Binary arithmetic
  add,
  subtract,
  multiply,
  divide,
  true_divide,
  floor_divide,
  remainder,
  mod,
  fmod,
  power,
  // Comparison
  equal,
  not_equal,
  less,
  less_equal,
  greater,
  greater_equal,
  // Extrema
  maximum,
  minimum,
  fmax,
  fmin,
  // Logical (binary)
  logical_and,
  logical_or,
  logical_xor,
  // Bitwise (binary)
  bitwise_and,
  bitwise_or,
  bitwise_xor,
  left_shift,
  right_shift,
  // Special math
  arctan2,
  hypot,
  copysign,
  signbit,
  logaddexp,
  logaddexp2,
} from './ufunc.js';

// Functional programming (Level 10)
export {
  applyAlongAxis,
  applyOverAxes,
  Vectorize,
  vectorize,
  frompyfunc,
  piecewise,
} from './functional.js';
export type { VectorizeOptions, UfuncLike } from './functional.js';

// Set operations (Phase 8)
export {
  // Unique functions
  unique,
  uniqueValues,
  uniqueIndex,
  uniqueInverse,
  uniqueCounts,
  uniqueAll,
  // Set combination operations
  union1d,
  intersect1d,
  setdiff1d,
  setxor1d,
  // Membership testing
  isin,
  in1d,
  // Differences
  ediff1d,
  // Constants
  ISIN_AUTO,
  ISIN_SORT,
  ISIN_TABLE,
} from './setops.js';
export type {
  UniqueResult,
  UniqueOptions,
  IsinOptions,
  IsinKind,
  IntersectOptions,
  SetDiffOptions,
  SetXorOptions,
  Ediff1dOptions,
} from './setops.js';

// I/O operations (Phase 9)
export {
  // NPY format
  save,
  load,
  // Text I/O
  loadtxt,
  savetxt,
  genfromtxt,
  fromregex,
  formatValue,
  // Binary I/O
  fromfile,
  frombuffer,
  // Array printing
  setPrintoptions,
  getPrintoptions,
  resetPrintoptions,
  withPrintoptions,
  array2string,
  arrayRepr,
  arrayStr,
  formatFloatPositional,
  formatFloatScientific,
  // Memory mapping
  Memmap,
  openMemmap,
  // Base conversion
  binaryRepr,
  baseRepr,
  fromBinaryRepr,
  fromBaseRepr,
  hexRepr,
  octalRepr,
  binaryReprArray,
  baseReprArray,
  // NPY format utilities
  dtypeToDescr,
  descrToDtype,
  isNode,
} from './io/index.js';
export type {
  SaveOptions,
  LoadOptions,
  LoadtxtOptions,
  SavetxtOptions,
  GenfromtxtOptions,
  FromregexOptions,
  FromfileOptions,
  FrombufferOptions,
  PrintOptions,
  MemmapMode,
  MemmapOptions,
} from './io/index.js';

// Window functions (Phase 11)
export {
  blackman,
  bartlett,
  hanning,
  hamming,
  kaiser,
  i0,
} from './window.js';

// Constants (Phase 12)
export {
  // Mathematical constants
  e,
  pi,
  euler_gamma,
  // Special floating-point values
  inf,
  PINF,
  NINF,
  nan,
  NAN,
  PZERO,
  NZERO,
} from './constants.js';

// Type information classes (Phase 12)
export { finfo, iinfo } from './typeinfo.js';
