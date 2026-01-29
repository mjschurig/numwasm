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
export {
  loadWasmModule,
  isWasmLoaded,
  getWasmModule,
  configureWasm,
} from './wasm-loader.js';
export type { WasmLoadConfig } from './wasm-loader.js';

// Array creation functions
export {
  zeros,
  ones,
  empty,
  full,
  arange,
  linspace,
  logspace,
  geomspace,
  eye,
  identity,
  diag,
  array,
  zeros_like,
  ones_like,
  empty_like,
  full_like,
} from './creation.js';

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
export { sort, argsort, partition, argpartition, sort_complex } from './sorting.js';
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
  // Cumulative operations (Phase 22)
  cumsum,
  cumprod,
  nancumsum,
  nancumprod,
} from './statistics.js';

// NaN-handling functions (Phase 23)
export {
  nansum,
  nanprod,
  nanmin,
  nanmax,
  nanargmin,
  nanargmax,
  nanmean,
  nanvar,
  nanstd,
  nanmedian,
  nanquantile,
  nanpercentile,
  nan_to_num,
} from './nanfunctions.js';

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
  // Complex utilities
  real_if_close,
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
  // Phase 26: Miscellaneous Ufuncs
  frexp,
  ldexp,
  nextafter,
  spacing,
  modf,
  gcd,
  lcm,
  sinc,
  heaviside,
  divmod,
  bitwise_count,
  // Complex number operations
  conjugate,
  conj,
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

// Histogram functions (Phase 24)
export {
  bincount,
  digitize,
  histogram_bin_edges,
  histogram,
  histogram2d,
  histogramdd,
  HistogramError,
} from './histogram.js';
export type {
  BinMethod,
  HistogramResult,
  Histogram2DResult,
  HistogramDDResult,
} from './histogram.js';

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

// Linear algebra (Phase 13)
export {
  linalg,
  LinAlgError,
  matmul,
  dot,
  vdot,
  inner,
  outer,
  cholesky,
  qr,
  svd,
  svdvals,
  eig,
  eigvals,
  eigh,
  eigvalsh,
  norm,
  det,
  slogdet,
  matrix_rank,
  trace,
  cond,
  solve,
  lstsq,
  inv,
  pinv,
  matrix_power,
  // Phase 25: Advanced Linear Algebra
  tensordot,
  multi_dot,
  kron,
  cross,
  tensorsolve,
  tensorinv,
  matrix_norm,
  vector_norm,
} from './linalg.js';
export type {
  EigResult,
  SVDResult,
  QRResult,
  LstsqResult,
  SlogdetResult,
} from './linalg.js';

// Einstein summation (Phase 25)
export { einsum, einsum_path } from './einsum.js';

// FFT (Phase 14)
export {
  // 1D transforms
  fft,
  ifft,
  rfft,
  irfft,
  hfft,
  ihfft,
  // 2D transforms
  fft2,
  ifft2,
  rfft2,
  irfft2,
  // N-D transforms
  fftn,
  ifftn,
  rfftn,
  irfftn,
  // Helper functions
  fftfreq,
  rfftfreq,
  fftshift,
  ifftshift,
  // Module object
  fftModule,
} from './fft.js';
export type { FFTNorm } from './fft.js';

// Random module (Phase 15 + Phase 27 BitGenerators)
export {
  // Classes
  Generator,
  PCG64,
  SeedSequence,
  BitGenerator,
  // Phase 27: Additional BitGenerators
  MT19937,
  Philox,
  SFC64,
  // Factory function
  default_rng,
  // BitGenerator registry
  getBitGenerator,
  listBitGenerators,
  // Module-level functions
  seed,
  random,
  randn,
  randint,
  // Initialization
  initRandom,
} from './random.js';
export type { SizeType, PCG64State, MT19937State, PhiloxState, SFC64State } from './random.js';


// String operations (Phase 16a)
export {
  // Error class
  ValueError,
  // Comparison functions
  equal as strEqual,
  not_equal as strNotEqual,
  less as strLess,
  less_equal as strLessEqual,
  greater as strGreater,
  greater_equal as strGreaterEqual,
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
  // Search functions
  find as strFind,
  rfind,
  index as strIndex,
  rindex,
  count as strCount,
  startswith,
  endswith,
  // Manipulation functions
  lower,
  upper,
  swapcase,
  capitalize as strCapitalize,
  title as strTitle,
  add as strAdd,
  multiply as strMultiply,
  strip,
  lstrip,
  rstrip,
  expandtabs,
  replace as strReplace,
  center as strCenter,
  ljust,
  rjust,
  zfill,
  partition as strPartition,
  rpartition,
  encode,
  decode,
  // Namespace object
  strings,
} from './strings/index.js';

// Polynomial module (Phase 16c)
export {
  // Utilities
  PolyError,
  PolyDomainWarning,
  trimseq,
  trimcoef,
  as_series,
  getdomain,
  mapdomain,
  mapparms,
  // Base class
  ABCPolyBase,
  maxpower,
  // Polynomial (power series)
  Polynomial,
  polyval,
  polyval2d,
  polyval3d,
  polyvander,
  polyvander2d,
  polyder,
  polyint,
  polyfit,
  polyroots,
  polycompanion,
  polyfromroots,
  polyadd,
  polysub,
  polymul,
  polydiv,
  polypow,
  // Chebyshev
  Chebyshev,
  chebval,
  chebval2d,
  chebvander,
  chebder,
  chebint,
  chebfit,
  chebroots,
  chebcompanion,
  chebfromroots,
  chebinterpolate,
  poly2cheb,
  cheb2poly,
  chebadd,
  chebsub,
  chebmul,
  chebdiv,
  chebpow,
  // Legendre
  Legendre,
  legval,
  legvander,
  legder,
  legint,
  legfit,
  legroots,
  legcompanion,
  legfromroots,
  poly2leg,
  leg2poly,
  legadd,
  legsub,
  legmul,
  legdiv,
  legpow,
  // Hermite (Physicist's)
  Hermite,
  hermval,
  hermvander,
  hermder,
  hermint,
  hermfit,
  hermroots,
  hermcompanion,
  hermfromroots,
  poly2herm,
  herm2poly,
  hermadd,
  hermsub,
  hermmul,
  hermdiv,
  hermpow,
  // HermiteE (Probabilist's)
  HermiteE,
  hermeval,
  hermevander,
  hermeder,
  hermeint,
  hermefit,
  hermeroots,
  hermecompanion,
  hermefromroots,
  poly2herme,
  herme2poly,
  hermeadd,
  hermesub,
  hermemul,
  hermediv,
  hermepow,
  // Laguerre
  Laguerre,
  lagval,
  lagvander,
  lagder,
  lagint,
  lagfit,
  lagroots,
  lagcompanion,
  lagfromroots,
  poly2lag,
  lag2poly,
  lagadd,
  lagsub,
  lagmul,
  lagdiv,
  lagpow,
} from './polynomial/index.js';

// Testing module (Phase 16e)
export * as testing from './testing/index.js';
export {
  AssertionError,
  SkipTest,
  KnownFailureException,
} from './testing/index.js';

// Record arrays module (Phase 16b)
export {
  rec,
  recarray,
  record,
  format_parser,
  fromarrays as recFromarrays,
  fromrecords as recFromrecords,
  fromstring as recFromstring,
  fromfile as recFromfile,
  array as recArray,
  find_duplicate,
  KeyError,
  IndexError as RecIndexError,
} from './rec/index.js';
export type {
  RecArrayOptions,
} from './rec/index.js';
export type {
  StructuredDType,
  FieldDescriptor,
} from './types.js';
export {
  isStructuredDType,
  dtypeSize as structDtypeSize,
  dtypeAlignment,
} from './types.js';

// Masked arrays module (Phase 16d)
export { ma } from './ma/index.js';
export { MaskedArray } from './ma/core.js';
export type { MaskType, MaskedConstant, SliceInfo } from './ma/index.js';

// Type checking and complex utilities
export { angle, real, imag } from './type_check.js';
