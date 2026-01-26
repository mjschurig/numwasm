# Phase 17: Error Handling & Configuration Implementation Plan

Complete implementation roadmap for the NumJS-WASM error handling system and configuration module, providing NumPy-compatible exceptions, error states, and runtime configuration.

---

## ⚠️ Implementation Guideline

**The original NumPy source code is available in `/numpy`.**

Key reference files:
- `numpy/_core/_exceptions.py` - Core exception classes (~200 lines)
- `numpy/_core/numerictypes.py` - Type-related errors
- `numpy/_core/arrayprint.py` - Print options (~1,779 lines)
- `numpy/_core/getlimits.py` - Machine limits/finfo (~600 lines)
- `numpy/_globals.py` - Global state management
- `numpy/exceptions.py` - Public exception API

Implementation should follow NumPy's error hierarchy and configuration patterns for consistency.

---

## Current State (Pre-Phase 17)

```
src/wasm/
├── ndarray.h/c        # Core NDArray with views, slicing
├── dtype.h/c          # DType system
├── broadcast.h/c      # Broadcasting
├── indexing.h/c       # Index operations
├── pairwise_sum.h/c   # Accurate summation
└── logic.c            # Logical operations

src/ts/
├── NDArray.ts         # Core array class
├── types.ts           # Type definitions
├── dtype.ts           # Type utilities
├── broadcast.ts       # Broadcasting functions
├── indexing.ts        # Index operations
├── slice.ts           # Slicing utilities
├── iterators.ts       # Iterators
└── index.ts           # Public exports
```

**Current Error Handling:**
- Basic JavaScript Error throws
- No structured exception hierarchy
- No floating-point error state tracking
- No configurable print options
- No machine limits (finfo/iinfo)

---

## Phase 17 Dependency Tree

```
PHASE 17: ERROR HANDLING & CONFIGURATION
│
├── 17.1 Exception Hierarchy (TypeScript)
│   ├── 17.1.1 Base Exceptions
│   │   ├── NumJSError (base class)
│   │   ├── AxisError
│   │   ├── DTypePromotionError
│   │   └── TooHardError
│   │
│   ├── 17.1.2 Array Exceptions
│   │   ├── BroadcastError (from AxisError)
│   │   ├── ComplexWarning
│   │   ├── VisibleDeprecationWarning
│   │   └── RankWarning
│   │
│   ├── 17.1.3 Module-Specific Exceptions
│   │   ├── LinAlgError (numpy.linalg)
│   │   ├── FFTError (numpy.fft)
│   │   └── RandomStateError (numpy.random)
│   │
│   └── 17.1.4 Warning System
│       ├── NumJSWarning (base)
│       ├── warn() function
│       ├── filterwarnings()
│       └── Warning categories
│
│   Dependencies: None (foundation)
│
├── 17.2 Floating-Point Error Handling (C/WASM + TypeScript)
│   ├── 17.2.1 Error State Flags
│   │   ├── FPE_DIVIDEBYZERO
│   │   ├── FPE_OVERFLOW
│   │   ├── FPE_UNDERFLOW
│   │   ├── FPE_INVALID
│   │   └── FPE_INEXACT (optional)
│   │
│   ├── 17.2.2 Error State Functions
│   │   ├── seterr(all, divide, over, under, invalid)
│   │   ├── geterr() → current settings
│   │   ├── seterrcall(func)
│   │   ├── geterrcall() → current callback
│   │   ├── errstate() context manager
│   │   └── floating-point status inspection
│   │
│   └── 17.2.3 Error Response Modes
│       ├── 'ignore' - silently continue
│       ├── 'warn' - issue RuntimeWarning
│       ├── 'raise' - raise FloatingPointError
│       ├── 'call' - invoke callback function
│       └── 'print' - print warning message
│
│   Dependencies: 17.1.* (Exception classes)
│
├── 17.3 Print Options & Formatting (TypeScript)
│   ├── 17.3.1 Print Configuration
│   │   ├── set_printoptions(precision, threshold, ...)
│   │   ├── get_printoptions() → current options
│   │   ├── printoptions() context manager
│   │   └── Default settings
│   │
│   ├── 17.3.2 Array Formatting
│   │   ├── array2string(a, max_line_width, precision, ...)
│   │   ├── array_repr(arr)
│   │   ├── array_str(a)
│   │   └── set_string_function(f, repr)
│   │
│   └── 17.3.3 Number Formatting
│       ├── format_float_positional(x, precision, ...)
│       ├── format_float_scientific(x, precision, ...)
│       └── Custom formatter support
│
│   Dependencies: NDArray core
│
├── 17.4 Machine Limits & Type Info (TypeScript)
│   ├── 17.4.1 Floating-Point Info (finfo)
│   │   ├── finfo(dtype) class
│   │   ├── bits, dtype, eps, epsneg
│   │   ├── iexp, machep, max, maxexp
│   │   ├── min, minexp, negep, nexp
│   │   ├── nmant, precision, resolution
│   │   ├── smallest_normal, smallest_subnormal
│   │   └── tiny
│   │
│   ├── 17.4.2 Integer Info (iinfo)
│   │   ├── iinfo(dtype) class
│   │   ├── bits, dtype, kind
│   │   ├── max, min
│   │   └── Signed/unsigned handling
│   │
│   └── 17.4.3 Machine Parameters
│       ├── MachAr class (legacy)
│       └── Platform-specific constants
│
│   Dependencies: DType system
│
└── 17.5 Global Configuration (TypeScript)
    ├── 17.5.1 Module Configuration
    │   ├── Configuration singleton
    │   ├── Thread-local state (for workers)
    │   └── Reset to defaults
    │
    ├── 17.5.2 Debug & Diagnostics
    │   ├── show_config()
    │   ├── __version__
    │   └── Debug flags
    │
    └── 17.5.3 Context Managers
        ├── errstate()
        ├── printoptions()
        └── Nested context support

    Dependencies: 17.2.*, 17.3.*
```

---

## Detailed Implementation Specifications

### 17.1 Exception Hierarchy

#### 17.1.1 Base Exceptions

**File:** `src/ts/exceptions.ts` (new file)

```typescript
/**
 * NumJS Exception Hierarchy
 *
 * Mirrors NumPy's exception structure for compatibility.
 * All NumJS-specific exceptions inherit from NumJSError.
 */

/* ============ Base Exception ============ */

/**
 * Base class for all NumJS exceptions.
 * Extends JavaScript Error with additional context.
 */
export class NumJSError extends Error {
  /** The array or value that caused the error, if applicable */
  public readonly context?: unknown;

  /** Error code for programmatic handling */
  public readonly code?: string;

  constructor(message: string, options?: { context?: unknown; code?: string }) {
    super(message);
    this.name = 'NumJSError';
    this.context = options?.context;
    this.code = options?.code;

    // Maintains proper stack trace in V8 environments
    if (Error.captureStackTrace) {
      Error.captureStackTrace(this, this.constructor);
    }
  }
}

/* ============ Axis Errors ============ */

/**
 * Exception raised when an axis is out of bounds or invalid.
 *
 * @example
 * // Axis 3 is out of bounds for array with 2 dimensions
 * throw new AxisError(3, 2);
 */
export class AxisError extends NumJSError {
  /** The invalid axis value */
  public readonly axis: number;

  /** Number of dimensions in the array */
  public readonly ndim: number;

  /** Optional custom message override */
  public readonly msgPrefix?: string;

  constructor(axis: number, ndim: number, msgPrefix?: string) {
    const message = msgPrefix
      ? `${msgPrefix}: axis ${axis} is out of bounds for array of dimension ${ndim}`
      : `axis ${axis} is out of bounds for array of dimension ${ndim}`;

    super(message, { code: 'AXIS_ERROR' });
    this.name = 'AxisError';
    this.axis = axis;
    this.ndim = ndim;
    this.msgPrefix = msgPrefix;
  }

  /**
   * Normalize axis value to positive index.
   * Throws AxisError if out of bounds.
   */
  static normalizeAxis(axis: number, ndim: number, msgPrefix?: string): number {
    const normalizedAxis = axis < 0 ? axis + ndim : axis;
    if (normalizedAxis < 0 || normalizedAxis >= ndim) {
      throw new AxisError(axis, ndim, msgPrefix);
    }
    return normalizedAxis;
  }

  /**
   * Normalize multiple axes.
   */
  static normalizeAxes(axes: number[], ndim: number, msgPrefix?: string): number[] {
    return axes.map(axis => AxisError.normalizeAxis(axis, ndim, msgPrefix));
  }
}

/* ============ DType Errors ============ */

/**
 * Exception raised when dtype promotion fails.
 *
 * @example
 * // Cannot safely cast float64 to int32
 * throw new DTypePromotionError(DType.Float64, DType.Int32);
 */
export class DTypePromotionError extends NumJSError {
  /** Source dtype */
  public readonly fromDtype: string;

  /** Target dtype */
  public readonly toDtype: string;

  constructor(fromDtype: string, toDtype: string, message?: string) {
    const msg = message ??
      `Cannot safely cast from ${fromDtype} to ${toDtype}`;

    super(msg, { code: 'DTYPE_PROMOTION_ERROR' });
    this.name = 'DTypePromotionError';
    this.fromDtype = fromDtype;
    this.toDtype = toDtype;
  }
}

/**
 * Exception for incompatible dtype operations.
 */
export class DTypeError extends NumJSError {
  constructor(message: string) {
    super(message, { code: 'DTYPE_ERROR' });
    this.name = 'DTypeError';
  }
}

/* ============ Broadcasting Errors ============ */

/**
 * Exception raised when array shapes cannot be broadcast together.
 *
 * @example
 * // shapes (3, 4) and (5, 6) are not broadcastable
 * throw new BroadcastError([3, 4], [5, 6]);
 */
export class BroadcastError extends NumJSError {
  /** Shapes that failed to broadcast */
  public readonly shapes: number[][];

  constructor(shapes: number[][], message?: string) {
    const shapeStrs = shapes.map(s => `(${s.join(', ')})`).join(' and ');
    const msg = message ??
      `operands could not be broadcast together with shapes ${shapeStrs}`;

    super(msg, { code: 'BROADCAST_ERROR' });
    this.name = 'BroadcastError';
    this.shapes = shapes;
  }
}

/* ============ Value Errors ============ */

/**
 * Exception raised for invalid values.
 */
export class ValueError extends NumJSError {
  constructor(message: string) {
    super(message, { code: 'VALUE_ERROR' });
    this.name = 'ValueError';
  }
}

/**
 * Exception raised for invalid types.
 */
export class TypeError extends NumJSError {
  constructor(message: string) {
    super(message, { code: 'TYPE_ERROR' });
    this.name = 'TypeError';
  }
}

/* ============ Index Errors ============ */

/**
 * Exception raised for invalid indices.
 */
export class IndexError extends NumJSError {
  /** The invalid index */
  public readonly index?: number | number[];

  /** The array size or shape */
  public readonly size?: number | number[];

  constructor(message: string, index?: number | number[], size?: number | number[]) {
    super(message, { code: 'INDEX_ERROR' });
    this.name = 'IndexError';
    this.index = index;
    this.size = size;
  }
}

/* ============ Shape Errors ============ */

/**
 * Exception raised for shape mismatches.
 */
export class ShapeError extends NumJSError {
  /** Expected shape */
  public readonly expected?: number[];

  /** Actual shape */
  public readonly actual?: number[];

  constructor(message: string, expected?: number[], actual?: number[]) {
    super(message, { code: 'SHAPE_ERROR' });
    this.name = 'ShapeError';
    this.expected = expected;
    this.actual = actual;
  }
}

/* ============ Memory Errors ============ */

/**
 * Exception raised when memory allocation fails.
 */
export class MemoryError extends NumJSError {
  /** Requested size in bytes */
  public readonly requestedSize?: number;

  constructor(message: string, requestedSize?: number) {
    super(message, { code: 'MEMORY_ERROR' });
    this.name = 'MemoryError';
    this.requestedSize = requestedSize;
  }
}

/* ============ Computation Errors ============ */

/**
 * Exception raised when computation is too complex.
 * Used for operations that would require unreasonable resources.
 */
export class TooHardError extends NumJSError {
  constructor(message: string) {
    super(message, { code: 'TOO_HARD_ERROR' });
    this.name = 'TooHardError';
  }
}

/**
 * Exception raised for floating-point errors.
 */
export class FloatingPointError extends NumJSError {
  /** Type of FP error: 'divide', 'overflow', 'underflow', 'invalid' */
  public readonly errorType: string;

  constructor(message: string, errorType: string) {
    super(message, { code: 'FLOATING_POINT_ERROR' });
    this.name = 'FloatingPointError';
    this.errorType = errorType;
  }
}

/* ============ Module-Specific Errors ============ */

/**
 * Exception raised by numpy.linalg functions.
 * Indicates that a matrix operation failed (e.g., singular matrix).
 */
export class LinAlgError extends NumJSError {
  constructor(message: string) {
    super(message, { code: 'LINALG_ERROR' });
    this.name = 'LinAlgError';
  }
}

/**
 * Exception raised by numpy.fft functions.
 */
export class FFTError extends NumJSError {
  constructor(message: string) {
    super(message, { code: 'FFT_ERROR' });
    this.name = 'FFTError';
  }
}

/**
 * Exception raised by numpy.random functions.
 */
export class RandomStateError extends NumJSError {
  constructor(message: string) {
    super(message, { code: 'RANDOM_STATE_ERROR' });
    this.name = 'RandomStateError';
  }
}

/* ============ Export All ============ */

export const exceptions = {
  NumJSError,
  AxisError,
  DTypePromotionError,
  DTypeError,
  BroadcastError,
  ValueError,
  TypeError,
  IndexError,
  ShapeError,
  MemoryError,
  TooHardError,
  FloatingPointError,
  LinAlgError,
  FFTError,
  RandomStateError,
};
```

---

#### 17.1.4 Warning System

**File:** `src/ts/warnings.ts` (new file)

```typescript
/**
 * NumJS Warning System
 *
 * Provides NumPy-compatible warning infrastructure.
 */

/* ============ Warning Categories ============ */

/**
 * Base class for NumJS warnings.
 */
export class NumJSWarning {
  public readonly message: string;
  public readonly category: string;
  public readonly stackLevel: number;
  public readonly source?: string;

  constructor(
    message: string,
    category: string = 'NumJSWarning',
    stackLevel: number = 1,
    source?: string
  ) {
    this.message = message;
    this.category = category;
    this.stackLevel = stackLevel;
    this.source = source;
  }
}

/**
 * Warning for visible deprecations that users should address.
 */
export class VisibleDeprecationWarning extends NumJSWarning {
  constructor(message: string, stackLevel: number = 1) {
    super(message, 'VisibleDeprecationWarning', stackLevel);
  }
}

/**
 * Warning when complex values are cast to real, discarding imaginary part.
 */
export class ComplexWarning extends NumJSWarning {
  constructor(message: string = 'Casting complex values to real discards the imaginary part') {
    super(message, 'ComplexWarning');
  }
}

/**
 * Warning for rank deficiency in numerical computations.
 */
export class RankWarning extends NumJSWarning {
  constructor(message: string) {
    super(message, 'RankWarning');
  }
}

/**
 * Warning for runtime issues (e.g., floating-point anomalies).
 */
export class RuntimeWarning extends NumJSWarning {
  constructor(message: string) {
    super(message, 'RuntimeWarning');
  }
}

/**
 * Warning for user-facing issues.
 */
export class UserWarning extends NumJSWarning {
  constructor(message: string) {
    super(message, 'UserWarning');
  }
}

/* ============ Warning Filter ============ */

type WarningAction = 'error' | 'ignore' | 'always' | 'default' | 'module' | 'once';

interface WarningFilter {
  action: WarningAction;
  message?: string | RegExp;
  category?: string;
  module?: string | RegExp;
  lineno?: number;
}

/* ============ Warning State ============ */

interface WarningState {
  filters: WarningFilter[];
  onceRegistry: Set<string>;
  defaultAction: WarningAction;
}

const warningState: WarningState = {
  filters: [],
  onceRegistry: new Set(),
  defaultAction: 'default',
};

/* ============ Warning Functions ============ */

/**
 * Handlers for different output modes.
 */
type WarningHandler = (warning: NumJSWarning) => void;

let warningHandler: WarningHandler = defaultWarningHandler;

function defaultWarningHandler(warning: NumJSWarning): void {
  console.warn(`${warning.category}: ${warning.message}`);
}

/**
 * Issue a warning.
 *
 * @param message - Warning message
 * @param category - Warning category class
 * @param stackLevel - Stack frames to skip for source location
 */
export function warn(
  message: string | NumJSWarning,
  category: typeof NumJSWarning = NumJSWarning,
  stackLevel: number = 1
): void {
  const warning = message instanceof NumJSWarning
    ? message
    : new category(message, category.name, stackLevel);

  const action = getWarningAction(warning);

  switch (action) {
    case 'error':
      throw new Error(`${warning.category}: ${warning.message}`);

    case 'ignore':
      return;

    case 'once': {
      const key = `${warning.category}:${warning.message}`;
      if (warningState.onceRegistry.has(key)) {
        return;
      }
      warningState.onceRegistry.add(key);
      warningHandler(warning);
      break;
    }

    case 'always':
    case 'default':
    default:
      warningHandler(warning);
      break;
  }
}

/**
 * Determine action for a warning based on filters.
 */
function getWarningAction(warning: NumJSWarning): WarningAction {
  for (const filter of warningState.filters) {
    if (matchesFilter(warning, filter)) {
      return filter.action;
    }
  }
  return warningState.defaultAction;
}

/**
 * Check if warning matches a filter.
 */
function matchesFilter(warning: NumJSWarning, filter: WarningFilter): boolean {
  if (filter.category && warning.category !== filter.category) {
    return false;
  }

  if (filter.message) {
    if (filter.message instanceof RegExp) {
      if (!filter.message.test(warning.message)) return false;
    } else {
      if (!warning.message.includes(filter.message)) return false;
    }
  }

  if (filter.module && warning.source) {
    if (filter.module instanceof RegExp) {
      if (!filter.module.test(warning.source)) return false;
    } else {
      if (!warning.source.includes(filter.module)) return false;
    }
  }

  return true;
}

/**
 * Add a warning filter.
 *
 * @param action - Action to take: 'error', 'ignore', 'always', 'default', 'module', 'once'
 * @param message - Pattern to match against warning message
 * @param category - Warning category to match
 * @param module - Module pattern to match
 * @param lineno - Line number to match (0 matches all)
 * @param append - If true, append to filters; otherwise prepend
 */
export function filterwarnings(
  action: WarningAction,
  message?: string | RegExp,
  category?: string,
  module?: string | RegExp,
  lineno?: number,
  append: boolean = false
): void {
  const filter: WarningFilter = {
    action,
    message,
    category,
    module,
    lineno,
  };

  if (append) {
    warningState.filters.push(filter);
  } else {
    warningState.filters.unshift(filter);
  }
}

/**
 * Reset warning filters to default state.
 */
export function resetwarnings(): void {
  warningState.filters = [];
  warningState.onceRegistry.clear();
  warningState.defaultAction = 'default';
}

/**
 * Temporarily modify warning behavior.
 * Returns a restore function.
 */
export function simplefilter(action: WarningAction): () => void {
  const oldFilters = [...warningState.filters];
  const oldDefault = warningState.defaultAction;

  warningState.filters = [];
  warningState.defaultAction = action;

  return () => {
    warningState.filters = oldFilters;
    warningState.defaultAction = oldDefault;
  };
}

/**
 * Set custom warning handler.
 */
export function setWarningHandler(handler: WarningHandler): void {
  warningHandler = handler;
}

/**
 * Get current warning handler.
 */
export function getWarningHandler(): WarningHandler {
  return warningHandler;
}

/* ============ Exports ============ */

export const warnings = {
  // Warning classes
  NumJSWarning,
  VisibleDeprecationWarning,
  ComplexWarning,
  RankWarning,
  RuntimeWarning,
  UserWarning,

  // Functions
  warn,
  filterwarnings,
  resetwarnings,
  simplefilter,
  setWarningHandler,
  getWarningHandler,
};
```

---

### 17.2 Floating-Point Error Handling

#### 17.2.1 Error State (C/WASM)

**File:** `src/wasm/fperr.h` (new file)

```c
#ifndef NUMJS_FPERR_H
#define NUMJS_FPERR_H

#include <stdint.h>

/* ============ Floating-Point Error Flags ============ */

#define FPE_DIVIDEBYZERO  0x01
#define FPE_OVERFLOW      0x02
#define FPE_UNDERFLOW     0x04
#define FPE_INVALID       0x08
#define FPE_INEXACT       0x10
#define FPE_ALL           0x1F

/* ============ Error Response Modes ============ */

#define FPE_MODE_IGNORE   0
#define FPE_MODE_WARN     1
#define FPE_MODE_RAISE    2
#define FPE_MODE_CALL     3
#define FPE_MODE_PRINT    4

/* ============ Error State Structure ============ */

typedef struct {
    uint8_t flags;           /* Current error flags */
    uint8_t divide_mode;     /* Mode for divide-by-zero */
    uint8_t over_mode;       /* Mode for overflow */
    uint8_t under_mode;      /* Mode for underflow */
    uint8_t invalid_mode;    /* Mode for invalid operation */
} FPErrState;

/* ============ Global Error State ============ */

/**
 * Get pointer to current FP error state.
 */
FPErrState* fperr_get_state(void);

/**
 * Clear all error flags.
 */
void fperr_clear(void);

/**
 * Clear specific error flags.
 */
void fperr_clear_flags(uint8_t mask);

/**
 * Get current error flags.
 */
uint8_t fperr_get_flags(void);

/**
 * Set error flags (for testing/simulation).
 */
void fperr_set_flags(uint8_t flags);

/**
 * Check if specific error occurred.
 */
int fperr_check(uint8_t flag);

/* ============ Error Mode Configuration ============ */

/**
 * Set error handling mode for all error types.
 */
void fperr_seterr_all(uint8_t mode);

/**
 * Set error handling mode for divide-by-zero.
 */
void fperr_seterr_divide(uint8_t mode);

/**
 * Set error handling mode for overflow.
 */
void fperr_seterr_over(uint8_t mode);

/**
 * Set error handling mode for underflow.
 */
void fperr_seterr_under(uint8_t mode);

/**
 * Set error handling mode for invalid operations.
 */
void fperr_seterr_invalid(uint8_t mode);

/**
 * Get current error modes.
 * Returns packed modes: divide | (over << 4) | (under << 8) | (invalid << 12)
 */
uint32_t fperr_geterr(void);

/* ============ Error Detection Helpers ============ */

/**
 * Check if value is NaN and set invalid flag if so.
 */
int fperr_check_nan(double value);

/**
 * Check if value is infinite and set overflow flag if so.
 */
int fperr_check_inf(double value);

/**
 * Check for divide by zero.
 */
int fperr_check_divide(double divisor);

/**
 * Check for underflow (value very close to zero).
 */
int fperr_check_underflow(double value, double threshold);

/* ============ Safe Math Operations ============ */

/**
 * Safe division with error tracking.
 */
double fperr_divide(double a, double b);

/**
 * Safe log with error tracking.
 */
double fperr_log(double x);

/**
 * Safe sqrt with error tracking.
 */
double fperr_sqrt(double x);

/**
 * Safe power with error tracking.
 */
double fperr_pow(double base, double exp);

#endif /* NUMJS_FPERR_H */
```

**File:** `src/wasm/fperr.c` (new file)

```c
#include "fperr.h"
#include <math.h>
#include <float.h>

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#define EXPORT EMSCRIPTEN_KEEPALIVE
#else
#define EXPORT
#endif

/* ============ Global State ============ */

static FPErrState g_fperr_state = {
    .flags = 0,
    .divide_mode = FPE_MODE_WARN,
    .over_mode = FPE_MODE_WARN,
    .under_mode = FPE_MODE_IGNORE,
    .invalid_mode = FPE_MODE_WARN,
};

/* ============ State Access ============ */

EXPORT FPErrState* fperr_get_state(void)
{
    return &g_fperr_state;
}

EXPORT void fperr_clear(void)
{
    g_fperr_state.flags = 0;
}

EXPORT void fperr_clear_flags(uint8_t mask)
{
    g_fperr_state.flags &= ~mask;
}

EXPORT uint8_t fperr_get_flags(void)
{
    return g_fperr_state.flags;
}

EXPORT void fperr_set_flags(uint8_t flags)
{
    g_fperr_state.flags |= flags;
}

EXPORT int fperr_check(uint8_t flag)
{
    return (g_fperr_state.flags & flag) != 0;
}

/* ============ Mode Configuration ============ */

EXPORT void fperr_seterr_all(uint8_t mode)
{
    g_fperr_state.divide_mode = mode;
    g_fperr_state.over_mode = mode;
    g_fperr_state.under_mode = mode;
    g_fperr_state.invalid_mode = mode;
}

EXPORT void fperr_seterr_divide(uint8_t mode)
{
    g_fperr_state.divide_mode = mode;
}

EXPORT void fperr_seterr_over(uint8_t mode)
{
    g_fperr_state.over_mode = mode;
}

EXPORT void fperr_seterr_under(uint8_t mode)
{
    g_fperr_state.under_mode = mode;
}

EXPORT void fperr_seterr_invalid(uint8_t mode)
{
    g_fperr_state.invalid_mode = mode;
}

EXPORT uint32_t fperr_geterr(void)
{
    return (uint32_t)g_fperr_state.divide_mode |
           ((uint32_t)g_fperr_state.over_mode << 4) |
           ((uint32_t)g_fperr_state.under_mode << 8) |
           ((uint32_t)g_fperr_state.invalid_mode << 12);
}

/* ============ Error Detection ============ */

EXPORT int fperr_check_nan(double value)
{
    if (isnan(value)) {
        g_fperr_state.flags |= FPE_INVALID;
        return 1;
    }
    return 0;
}

EXPORT int fperr_check_inf(double value)
{
    if (isinf(value)) {
        g_fperr_state.flags |= FPE_OVERFLOW;
        return 1;
    }
    return 0;
}

EXPORT int fperr_check_divide(double divisor)
{
    if (divisor == 0.0) {
        g_fperr_state.flags |= FPE_DIVIDEBYZERO;
        return 1;
    }
    return 0;
}

EXPORT int fperr_check_underflow(double value, double threshold)
{
    if (value != 0.0 && fabs(value) < threshold) {
        g_fperr_state.flags |= FPE_UNDERFLOW;
        return 1;
    }
    return 0;
}

/* ============ Safe Math Operations ============ */

EXPORT double fperr_divide(double a, double b)
{
    if (b == 0.0) {
        g_fperr_state.flags |= FPE_DIVIDEBYZERO;
        if (a == 0.0) {
            g_fperr_state.flags |= FPE_INVALID;
            return NAN;
        }
        return a > 0 ? INFINITY : -INFINITY;
    }

    double result = a / b;

    if (isinf(result) && !isinf(a)) {
        g_fperr_state.flags |= FPE_OVERFLOW;
    }
    if (result == 0.0 && a != 0.0) {
        g_fperr_state.flags |= FPE_UNDERFLOW;
    }

    return result;
}

EXPORT double fperr_log(double x)
{
    if (x < 0.0) {
        g_fperr_state.flags |= FPE_INVALID;
        return NAN;
    }
    if (x == 0.0) {
        g_fperr_state.flags |= FPE_DIVIDEBYZERO;
        return -INFINITY;
    }
    return log(x);
}

EXPORT double fperr_sqrt(double x)
{
    if (x < 0.0) {
        g_fperr_state.flags |= FPE_INVALID;
        return NAN;
    }
    return sqrt(x);
}

EXPORT double fperr_pow(double base, double exp)
{
    /* Handle special cases */
    if (base < 0.0 && floor(exp) != exp) {
        g_fperr_state.flags |= FPE_INVALID;
        return NAN;
    }
    if (base == 0.0 && exp < 0.0) {
        g_fperr_state.flags |= FPE_DIVIDEBYZERO;
        return INFINITY;
    }

    double result = pow(base, exp);

    if (isinf(result) && !isinf(base) && !isinf(exp)) {
        g_fperr_state.flags |= FPE_OVERFLOW;
    }
    if (result == 0.0 && base != 0.0 && exp > 0.0) {
        g_fperr_state.flags |= FPE_UNDERFLOW;
    }

    return result;
}
```

---

#### 17.2.2 TypeScript Error State Interface

**File:** `src/ts/errstate.ts` (new file)

```typescript
/**
 * NumJS Floating-Point Error State Management
 *
 * Provides NumPy-compatible floating-point error handling.
 */

import { FloatingPointError } from './exceptions.js';
import { warn, RuntimeWarning } from './warnings.js';

/* ============ Types ============ */

export type ErrorMode = 'ignore' | 'warn' | 'raise' | 'call' | 'print';

export interface ErrorSettings {
  all?: ErrorMode;
  divide?: ErrorMode;
  over?: ErrorMode;
  under?: ErrorMode;
  invalid?: ErrorMode;
}

export interface ErrorState {
  divide: ErrorMode;
  over: ErrorMode;
  under: ErrorMode;
  invalid: ErrorMode;
}

export type ErrorCallback = (err: string, flag: number) => void;

/* ============ Constants ============ */

const FPE_DIVIDEBYZERO = 0x01;
const FPE_OVERFLOW = 0x02;
const FPE_UNDERFLOW = 0x04;
const FPE_INVALID = 0x08;

const MODE_IGNORE = 0;
const MODE_WARN = 1;
const MODE_RAISE = 2;
const MODE_CALL = 3;
const MODE_PRINT = 4;

/* ============ State ============ */

let currentState: ErrorState = {
  divide: 'warn',
  over: 'warn',
  under: 'ignore',
  invalid: 'warn',
};

let errorCallback: ErrorCallback | null = null;
let stateStack: ErrorState[] = [];

/* ============ Mode Conversion ============ */

function modeToNumber(mode: ErrorMode): number {
  switch (mode) {
    case 'ignore': return MODE_IGNORE;
    case 'warn': return MODE_WARN;
    case 'raise': return MODE_RAISE;
    case 'call': return MODE_CALL;
    case 'print': return MODE_PRINT;
    default: return MODE_WARN;
  }
}

function numberToMode(num: number): ErrorMode {
  switch (num) {
    case MODE_IGNORE: return 'ignore';
    case MODE_WARN: return 'warn';
    case MODE_RAISE: return 'raise';
    case MODE_CALL: return 'call';
    case MODE_PRINT: return 'print';
    default: return 'warn';
  }
}

/* ============ Public API ============ */

/**
 * Set how floating-point errors are handled.
 *
 * @param settings - Error handling settings
 * @returns Previous settings
 *
 * @example
 * // Ignore all floating-point errors
 * const old = seterr({ all: 'ignore' });
 *
 * // Raise on divide by zero, warn on overflow
 * seterr({ divide: 'raise', over: 'warn' });
 *
 * // Restore previous settings
 * seterr(old);
 */
export function seterr(settings: ErrorSettings): ErrorState {
  const old = { ...currentState };

  if (settings.all !== undefined) {
    currentState.divide = settings.all;
    currentState.over = settings.all;
    currentState.under = settings.all;
    currentState.invalid = settings.all;
  }

  if (settings.divide !== undefined) currentState.divide = settings.divide;
  if (settings.over !== undefined) currentState.over = settings.over;
  if (settings.under !== undefined) currentState.under = settings.under;
  if (settings.invalid !== undefined) currentState.invalid = settings.invalid;

  // Sync to WASM if available
  syncToWasm();

  return old;
}

/**
 * Get current floating-point error handling settings.
 *
 * @returns Current error settings
 *
 * @example
 * const settings = geterr();
 * console.log(settings.divide);  // 'warn'
 */
export function geterr(): ErrorState {
  return { ...currentState };
}

/**
 * Set the callback for 'call' mode.
 *
 * @param func - Callback function or null to clear
 * @returns Previous callback
 */
export function seterrcall(func: ErrorCallback | null): ErrorCallback | null {
  const old = errorCallback;
  errorCallback = func;
  return old;
}

/**
 * Get current error callback.
 */
export function geterrcall(): ErrorCallback | null {
  return errorCallback;
}

/**
 * Context manager for temporary error settings.
 * Returns a dispose function to restore previous settings.
 *
 * @param settings - Temporary error settings
 * @returns Dispose function
 *
 * @example
 * const restore = errstate({ divide: 'ignore' });
 * try {
 *   // Code that might divide by zero
 * } finally {
 *   restore();
 * }
 */
export function errstate(settings: ErrorSettings): () => void {
  stateStack.push({ ...currentState });
  seterr(settings);

  return () => {
    const previous = stateStack.pop();
    if (previous) {
      currentState = previous;
      syncToWasm();
    }
  };
}

/* ============ Error Checking ============ */

/**
 * Check for floating-point errors and handle according to current settings.
 * Called after operations that might produce FP errors.
 *
 * @param flags - Error flags from WASM
 */
export function checkFPErrors(flags: number): void {
  if (flags === 0) return;

  if (flags & FPE_DIVIDEBYZERO) {
    handleError('divide', 'divide by zero encountered');
  }
  if (flags & FPE_OVERFLOW) {
    handleError('over', 'overflow encountered');
  }
  if (flags & FPE_UNDERFLOW) {
    handleError('under', 'underflow encountered');
  }
  if (flags & FPE_INVALID) {
    handleError('invalid', 'invalid value encountered');
  }
}

/**
 * Handle a specific floating-point error.
 */
function handleError(type: keyof ErrorState, message: string): void {
  const mode = currentState[type];

  switch (mode) {
    case 'ignore':
      break;

    case 'warn':
      warn(new RuntimeWarning(message));
      break;

    case 'raise':
      throw new FloatingPointError(message, type);

    case 'call':
      if (errorCallback) {
        const flag = type === 'divide' ? FPE_DIVIDEBYZERO :
                     type === 'over' ? FPE_OVERFLOW :
                     type === 'under' ? FPE_UNDERFLOW : FPE_INVALID;
        errorCallback(message, flag);
      }
      break;

    case 'print':
      console.warn(`FloatingPointError: ${message}`);
      break;
  }
}

/* ============ WASM Sync ============ */

let wasmModule: any = null;

/**
 * Set WASM module reference for state synchronization.
 */
export function setWasmModule(module: any): void {
  wasmModule = module;
  syncToWasm();
}

/**
 * Sync current state to WASM module.
 */
function syncToWasm(): void {
  if (!wasmModule) return;

  try {
    wasmModule._fperr_seterr_divide(modeToNumber(currentState.divide));
    wasmModule._fperr_seterr_over(modeToNumber(currentState.over));
    wasmModule._fperr_seterr_under(modeToNumber(currentState.under));
    wasmModule._fperr_seterr_invalid(modeToNumber(currentState.invalid));
  } catch (e) {
    // WASM not available, ignore
  }
}

/**
 * Clear WASM error flags.
 */
export function clearFPErrors(): void {
  if (wasmModule) {
    try {
      wasmModule._fperr_clear();
    } catch (e) {
      // WASM not available, ignore
    }
  }
}

/**
 * Get WASM error flags.
 */
export function getFPErrors(): number {
  if (wasmModule) {
    try {
      return wasmModule._fperr_get_flags();
    } catch (e) {
      return 0;
    }
  }
  return 0;
}

/* ============ Exports ============ */

export const errhandling = {
  seterr,
  geterr,
  seterrcall,
  geterrcall,
  errstate,
  checkFPErrors,
  clearFPErrors,
  getFPErrors,
  setWasmModule,
};
```

---

### 17.3 Print Options & Formatting

**File:** `src/ts/printoptions.ts` (new file)

```typescript
/**
 * NumJS Print Options and Array Formatting
 *
 * Provides NumPy-compatible array printing configuration.
 */

/* ============ Types ============ */

export interface PrintOptions {
  /** Number of digits of precision for floating point output (default: 8) */
  precision: number;

  /** Total number of array elements which trigger summarization (default: 1000) */
  threshold: number;

  /** Number of array elements in summary at beginning and end (default: 3) */
  edgeitems: number;

  /** Number of characters per line (default: 75) */
  linewidth: number;

  /** If true, suppress printing of small floating point values (default: false) */
  suppress: boolean;

  /** String inserted between elements */
  separator: string;

  /** Prefix for the continuation line */
  prefix: string;

  /** Controls sign printing: '-', '+', or ' ' (default: '-') */
  sign: '-' | '+' | ' ';

  /** Controls float mode: 'fixed', 'unique', 'maxprec', or 'maxprec_equal' */
  floatmode: 'fixed' | 'unique' | 'maxprec' | 'maxprec_equal';

  /** String representation of infinity (default: 'inf') */
  infstr: string;

  /** String representation of NaN (default: 'nan') */
  nanstr: string;

  /** If true, use legacy printing mode (NumPy 1.13 style) */
  legacy: boolean | '1.13' | '1.21' | '1.25';
}

/* ============ Default Options ============ */

const defaultOptions: PrintOptions = {
  precision: 8,
  threshold: 1000,
  edgeitems: 3,
  linewidth: 75,
  suppress: false,
  separator: ' ',
  prefix: '',
  sign: '-',
  floatmode: 'maxprec',
  infstr: 'inf',
  nanstr: 'nan',
  legacy: false,
};

let currentOptions: PrintOptions = { ...defaultOptions };
const optionsStack: PrintOptions[] = [];

/* ============ Public API ============ */

/**
 * Set printing options.
 *
 * @param options - New print options (partial)
 *
 * @example
 * set_printoptions({ precision: 4, threshold: 100 });
 */
export function set_printoptions(options: Partial<PrintOptions>): void {
  currentOptions = { ...currentOptions, ...options };
}

/**
 * Get current printing options.
 *
 * @returns Current print options
 */
export function get_printoptions(): PrintOptions {
  return { ...currentOptions };
}

/**
 * Context manager for temporary print options.
 * Returns a dispose function to restore previous settings.
 *
 * @param options - Temporary print options
 * @returns Dispose function
 *
 * @example
 * const restore = printoptions({ precision: 2 });
 * try {
 *   console.log(array.toString());
 * } finally {
 *   restore();
 * }
 */
export function printoptions(options: Partial<PrintOptions>): () => void {
  optionsStack.push({ ...currentOptions });
  set_printoptions(options);

  return () => {
    const previous = optionsStack.pop();
    if (previous) {
      currentOptions = previous;
    }
  };
}

/**
 * Reset print options to defaults.
 */
export function reset_printoptions(): void {
  currentOptions = { ...defaultOptions };
}

/* ============ Number Formatting ============ */

/**
 * Format a floating-point number using positional notation.
 *
 * @param x - The value to format
 * @param precision - Number of digits after decimal point
 * @param unique - If true, use shortest representation
 * @param fractional - If true, always include fractional part
 * @param trim - Trim trailing zeros: 'k', '.', '0', or '-'
 * @param sign - Sign mode: '-', '+', or ' '
 * @param padLeft - Pad to this width on the left
 * @param padRight - Pad to this width on the right
 * @returns Formatted string
 */
export function format_float_positional(
  x: number,
  precision?: number,
  unique: boolean = true,
  fractional: boolean = true,
  trim: 'k' | '.' | '0' | '-' = 'k',
  sign: '-' | '+' | ' ' = '-',
  padLeft?: number,
  padRight?: number
): string {
  const opts = get_printoptions();
  const prec = precision ?? opts.precision;

  // Handle special values
  if (Number.isNaN(x)) return opts.nanstr;
  if (!Number.isFinite(x)) return x > 0 ? opts.infstr : `-${opts.infstr}`;

  // Format number
  let str: string;
  if (unique) {
    str = formatUnique(x, prec);
  } else {
    str = x.toFixed(prec);
  }

  // Apply sign
  if (x >= 0) {
    if (sign === '+') str = '+' + str;
    else if (sign === ' ') str = ' ' + str;
  }

  // Apply trim
  str = applyTrim(str, trim);

  // Apply padding
  if (padLeft && str.length < padLeft) {
    str = str.padStart(padLeft);
  }
  if (padRight && str.length < padRight) {
    str = str.padEnd(padRight);
  }

  return str;
}

/**
 * Format a floating-point number using scientific notation.
 *
 * @param x - The value to format
 * @param precision - Number of digits after decimal point
 * @param unique - If true, use shortest representation
 * @param trim - Trim trailing zeros
 * @param sign - Sign mode
 * @param padLeft - Pad to this width on the left
 * @param expDigits - Minimum digits in exponent
 * @returns Formatted string
 */
export function format_float_scientific(
  x: number,
  precision?: number,
  unique: boolean = true,
  trim: 'k' | '.' | '0' | '-' = 'k',
  sign: '-' | '+' | ' ' = '-',
  padLeft?: number,
  expDigits: number = 2
): string {
  const opts = get_printoptions();
  const prec = precision ?? opts.precision;

  // Handle special values
  if (Number.isNaN(x)) return opts.nanstr;
  if (!Number.isFinite(x)) return x > 0 ? opts.infstr : `-${opts.infstr}`;

  // Format in scientific notation
  let str = x.toExponential(prec);

  // Normalize exponent format
  str = normalizeExponent(str, expDigits);

  // Apply sign
  if (x >= 0) {
    if (sign === '+') str = '+' + str;
    else if (sign === ' ') str = ' ' + str;
  }

  // Apply trim to mantissa
  str = applyScientificTrim(str, trim);

  // Apply padding
  if (padLeft && str.length < padLeft) {
    str = str.padStart(padLeft);
  }

  return str;
}

/* ============ Array Formatting ============ */

/**
 * Return a string representation of an array.
 *
 * @param a - Input array
 * @param maxLineWidth - Maximum line width
 * @param precision - Floating point precision
 * @param suppressSmall - Suppress small values
 * @param separator - Element separator
 * @param prefix - Line prefix
 * @param suffix - Line suffix
 * @param threshold - Summarization threshold
 * @param edgeitems - Number of edge items in summary
 * @returns String representation
 */
export function array2string(
  a: { shape: number[]; toArray(): number[] | any[]; dtype: any },
  maxLineWidth?: number,
  precision?: number,
  suppressSmall?: boolean,
  separator?: string,
  prefix?: string,
  suffix?: string,
  threshold?: number,
  edgeitems?: number
): string {
  const opts = get_printoptions();

  const lineWidth = maxLineWidth ?? opts.linewidth;
  const prec = precision ?? opts.precision;
  const suppress = suppressSmall ?? opts.suppress;
  const sep = separator ?? opts.separator;
  const pfx = prefix ?? opts.prefix;
  const thresh = threshold ?? opts.threshold;
  const edge = edgeitems ?? opts.edgeitems;

  return formatArray(a, {
    lineWidth,
    precision: prec,
    suppress,
    separator: sep,
    prefix: pfx,
    suffix: suffix ?? '',
    threshold: thresh,
    edgeitems: edge,
  });
}

/**
 * Return the string representation of an array.
 * Uses repr-style formatting with full precision.
 */
export function array_repr(
  arr: { shape: number[]; dtype: any; toArray(): number[] | any[] },
  maxLineWidth?: number,
  precision?: number,
  suppressSmall?: boolean
): string {
  const opts = get_printoptions();
  const className = 'array';

  const content = array2string(
    arr,
    maxLineWidth,
    precision,
    suppressSmall
  );

  // Format as array(..., dtype=...)
  const dtypeName = getDTypeName(arr.dtype);
  return `${className}(${content}, dtype=${dtypeName})`;
}

/**
 * Return a string representation without dtype info.
 */
export function array_str(
  a: { shape: number[]; toArray(): number[] | any[] },
  maxLineWidth?: number,
  precision?: number,
  suppressSmall?: boolean
): string {
  return array2string(a, maxLineWidth, precision, suppressSmall);
}

/* ============ Helper Functions ============ */

function formatUnique(x: number, maxPrecision: number): string {
  // Find shortest representation that round-trips
  for (let p = 1; p <= maxPrecision; p++) {
    const str = x.toFixed(p);
    if (parseFloat(str) === x) {
      return str;
    }
  }
  return x.toFixed(maxPrecision);
}

function applyTrim(str: string, trim: 'k' | '.' | '0' | '-'): string {
  if (trim === 'k') return str;  // Keep as-is
  if (trim === '-') return str.replace(/\.?0+$/, '');  // Trim zeros and trailing dot
  if (trim === '0') return str.replace(/0+$/, '');  // Trim only trailing zeros
  if (trim === '.') {
    // Trim only trailing dot
    return str.endsWith('.') ? str.slice(0, -1) : str;
  }
  return str;
}

function applyScientificTrim(str: string, trim: 'k' | '.' | '0' | '-'): string {
  if (trim === 'k') return str;

  const eIndex = str.indexOf('e');
  if (eIndex === -1) return str;

  const mantissa = str.slice(0, eIndex);
  const exponent = str.slice(eIndex);

  return applyTrim(mantissa, trim) + exponent;
}

function normalizeExponent(str: string, minDigits: number): string {
  const match = str.match(/e([+-])(\d+)$/);
  if (!match) return str;

  const [, sign, digits] = match;
  const paddedDigits = digits.padStart(minDigits, '0');

  return str.replace(/e[+-]\d+$/, `e${sign}${paddedDigits}`);
}

interface FormatOptions {
  lineWidth: number;
  precision: number;
  suppress: boolean;
  separator: string;
  prefix: string;
  suffix: string;
  threshold: number;
  edgeitems: number;
}

function formatArray(
  arr: { shape: number[]; toArray(): number[] | any[] },
  opts: FormatOptions
): string {
  const data = arr.toArray();
  const shape = arr.shape;

  if (shape.length === 0) {
    // Scalar
    return formatElement(data as unknown as number, opts);
  }

  if (shape.length === 1) {
    // 1D array
    return formatArray1D(data as number[], opts);
  }

  // Multi-dimensional: recursively format
  return formatArrayND(data, shape, 0, opts);
}

function formatElement(x: number, opts: FormatOptions): string {
  if (opts.suppress && Math.abs(x) < 1e-10) {
    return '0.';
  }
  return format_float_positional(x, opts.precision);
}

function formatArray1D(arr: number[], opts: FormatOptions): string {
  const n = arr.length;

  // Check if summarization needed
  if (n > opts.threshold) {
    const head = arr.slice(0, opts.edgeitems);
    const tail = arr.slice(-opts.edgeitems);
    const headStr = head.map(x => formatElement(x, opts)).join(opts.separator);
    const tailStr = tail.map(x => formatElement(x, opts)).join(opts.separator);
    return `[${headStr}, ..., ${tailStr}]`;
  }

  const elements = arr.map(x => formatElement(x, opts)).join(opts.separator);
  return `[${elements}]`;
}

function formatArrayND(
  data: any[],
  shape: number[],
  depth: number,
  opts: FormatOptions
): string {
  if (depth === shape.length - 1) {
    return formatArray1D(data as number[], opts);
  }

  const n = shape[depth];
  const subShape = shape.slice(depth + 1);
  const subSize = subShape.reduce((a, b) => a * b, 1);

  const lines: string[] = [];
  const indent = ' '.repeat(depth + 1);

  // Check if summarization needed at this level
  if (n > opts.threshold / subSize) {
    // Summarize
    for (let i = 0; i < opts.edgeitems; i++) {
      const start = i * subSize;
      const subData = data.slice(start, start + subSize);
      lines.push(indent + formatArrayND(subData, shape, depth + 1, opts));
    }
    lines.push(indent + '...');
    for (let i = n - opts.edgeitems; i < n; i++) {
      const start = i * subSize;
      const subData = data.slice(start, start + subSize);
      lines.push(indent + formatArrayND(subData, shape, depth + 1, opts));
    }
  } else {
    for (let i = 0; i < n; i++) {
      const start = i * subSize;
      const subData = data.slice(start, start + subSize);
      lines.push(indent + formatArrayND(subData, shape, depth + 1, opts));
    }
  }

  return '[\n' + lines.join(',\n') + '\n' + ' '.repeat(depth) + ']';
}

function getDTypeName(dtype: any): string {
  // Map dtype to NumPy-compatible name
  const names: Record<number, string> = {
    0: 'float64',
    1: 'float32',
    2: 'int32',
    3: 'int64',
    4: 'bool',
    5: 'int8',
    6: 'int16',
    7: 'uint8',
    8: 'uint16',
    9: 'uint32',
    10: 'uint64',
    11: 'float16',
    12: 'complex64',
    13: 'complex128',
  };
  return names[dtype] ?? 'unknown';
}

/* ============ Exports ============ */

export const printopts = {
  set_printoptions,
  get_printoptions,
  printoptions,
  reset_printoptions,
  format_float_positional,
  format_float_scientific,
  array2string,
  array_repr,
  array_str,
};
```

---

### 17.4 Machine Limits & Type Info

**File:** `src/ts/getlimits.ts` (new file)

```typescript
/**
 * NumJS Machine Limits and Type Information
 *
 * Provides NumPy-compatible finfo and iinfo classes.
 */

import { DType } from './types.js';

/* ============ Floating-Point Info ============ */

/**
 * Machine limits for floating point types.
 *
 * @example
 * const f64 = finfo(DType.Float64);
 * console.log(f64.eps);       // Machine epsilon
 * console.log(f64.max);       // Largest representable number
 * console.log(f64.precision); // Number of decimal digits
 */
export class finfo {
  /** The dtype for which this info applies */
  public readonly dtype: DType;

  /** Number of bits in the type */
  public readonly bits: number;

  /** Machine epsilon: smallest representable positive number such that 1.0 + eps != 1.0 */
  public readonly eps: number;

  /** Negative machine epsilon: smallest representable negative number such that 1.0 - epsneg != 1.0 */
  public readonly epsneg: number;

  /** Number of bits in the exponent */
  public readonly iexp: number;

  /** The exponent that generates eps */
  public readonly machep: number;

  /** The largest representable number */
  public readonly max: number;

  /** The smallest positive power of the base (2) that causes overflow */
  public readonly maxexp: number;

  /** The most negative representable number */
  public readonly min: number;

  /** The most negative power of the base (2) consistent with no leading zeros */
  public readonly minexp: number;

  /** The exponent that generates epsneg */
  public readonly negep: number;

  /** Number of bits in the exponent including its sign and bias */
  public readonly nexp: number;

  /** Number of bits in the mantissa */
  public readonly nmant: number;

  /** Approximate number of decimal digits */
  public readonly precision: number;

  /** Approximate decimal resolution */
  public readonly resolution: number;

  /** The smallest positive normal number */
  public readonly smallest_normal: number;

  /** The smallest positive subnormal number */
  public readonly smallest_subnormal: number;

  /** Alias for smallest_normal (for NumPy compatibility) */
  public readonly tiny: number;

  constructor(dtype: DType | string) {
    const dt = typeof dtype === 'string' ? parseDType(dtype) : dtype;

    switch (dt) {
      case DType.Float16:
        this.dtype = DType.Float16;
        this.bits = 16;
        this.eps = 9.765625e-4;  // 2^-10
        this.epsneg = 4.8828125e-4;  // 2^-11
        this.iexp = 5;
        this.machep = -10;
        this.max = 65504;
        this.maxexp = 16;
        this.min = -65504;
        this.minexp = -14;
        this.negep = -11;
        this.nexp = 5;
        this.nmant = 10;
        this.precision = 3;
        this.resolution = 1e-3;
        this.smallest_normal = 6.103515625e-5;  // 2^-14
        this.smallest_subnormal = 5.960464477539063e-8;  // 2^-24
        this.tiny = this.smallest_normal;
        break;

      case DType.Float32:
        this.dtype = DType.Float32;
        this.bits = 32;
        this.eps = 1.1920929e-7;  // 2^-23
        this.epsneg = 5.9604645e-8;  // 2^-24
        this.iexp = 8;
        this.machep = -23;
        this.max = 3.4028235e+38;
        this.maxexp = 128;
        this.min = -3.4028235e+38;
        this.minexp = -126;
        this.negep = -24;
        this.nexp = 8;
        this.nmant = 23;
        this.precision = 6;
        this.resolution = 1e-6;
        this.smallest_normal = 1.1754944e-38;  // 2^-126
        this.smallest_subnormal = 1.4012985e-45;  // 2^-149
        this.tiny = this.smallest_normal;
        break;

      case DType.Float64:
      default:
        this.dtype = DType.Float64;
        this.bits = 64;
        this.eps = 2.220446049250313e-16;  // 2^-52
        this.epsneg = 1.1102230246251565e-16;  // 2^-53
        this.iexp = 11;
        this.machep = -52;
        this.max = 1.7976931348623157e+308;
        this.maxexp = 1024;
        this.min = -1.7976931348623157e+308;
        this.minexp = -1022;
        this.negep = -53;
        this.nexp = 11;
        this.nmant = 52;
        this.precision = 15;
        this.resolution = 1e-15;
        this.smallest_normal = 2.2250738585072014e-308;  // 2^-1022
        this.smallest_subnormal = 5e-324;  // 2^-1074
        this.tiny = this.smallest_normal;
        break;
    }
  }

  /**
   * String representation.
   */
  toString(): string {
    return `Machine parameters for ${getDTypeName(this.dtype)}
---------------------------------------------------------------
precision = ${this.precision}   resolution = ${this.resolution}
machep = ${this.machep}   eps = ${this.eps}
negep = ${this.negep}   epsneg = ${this.epsneg}
minexp = ${this.minexp}   tiny = ${this.tiny}
maxexp = ${this.maxexp}   max = ${this.max}
nexp = ${this.nexp}   min = ${this.min}
smallest_normal = ${this.smallest_normal}
smallest_subnormal = ${this.smallest_subnormal}
---------------------------------------------------------------`;
  }
}

/* ============ Integer Info ============ */

/**
 * Machine limits for integer types.
 *
 * @example
 * const i32 = iinfo(DType.Int32);
 * console.log(i32.min);  // -2147483648
 * console.log(i32.max);  // 2147483647
 */
export class iinfo {
  /** The dtype for which this info applies */
  public readonly dtype: DType;

  /** Number of bits in the type */
  public readonly bits: number;

  /** Character code for the type: 'i' for signed, 'u' for unsigned */
  public readonly kind: 'i' | 'u';

  /** Minimum value */
  public readonly min: number | bigint;

  /** Maximum value */
  public readonly max: number | bigint;

  constructor(dtype: DType | string) {
    const dt = typeof dtype === 'string' ? parseDType(dtype) : dtype;

    switch (dt) {
      case DType.Bool:
        this.dtype = DType.Bool;
        this.bits = 8;
        this.kind = 'u';
        this.min = 0;
        this.max = 1;
        break;

      case DType.Int8:
        this.dtype = DType.Int8;
        this.bits = 8;
        this.kind = 'i';
        this.min = -128;
        this.max = 127;
        break;

      case DType.UInt8:
        this.dtype = DType.UInt8;
        this.bits = 8;
        this.kind = 'u';
        this.min = 0;
        this.max = 255;
        break;

      case DType.Int16:
        this.dtype = DType.Int16;
        this.bits = 16;
        this.kind = 'i';
        this.min = -32768;
        this.max = 32767;
        break;

      case DType.UInt16:
        this.dtype = DType.UInt16;
        this.bits = 16;
        this.kind = 'u';
        this.min = 0;
        this.max = 65535;
        break;

      case DType.Int32:
        this.dtype = DType.Int32;
        this.bits = 32;
        this.kind = 'i';
        this.min = -2147483648;
        this.max = 2147483647;
        break;

      case DType.UInt32:
        this.dtype = DType.UInt32;
        this.bits = 32;
        this.kind = 'u';
        this.min = 0;
        this.max = 4294967295;
        break;

      case DType.Int64:
        this.dtype = DType.Int64;
        this.bits = 64;
        this.kind = 'i';
        this.min = BigInt('-9223372036854775808');
        this.max = BigInt('9223372036854775807');
        break;

      case DType.UInt64:
        this.dtype = DType.UInt64;
        this.bits = 64;
        this.kind = 'u';
        this.min = BigInt(0);
        this.max = BigInt('18446744073709551615');
        break;

      default:
        throw new Error(`Cannot create iinfo for non-integer dtype: ${dt}`);
    }
  }

  /**
   * String representation.
   */
  toString(): string {
    return `Machine parameters for ${getDTypeName(this.dtype)}
---------------------------------------------------------------
min = ${this.min}
max = ${this.max}
---------------------------------------------------------------`;
  }
}

/* ============ Helper Functions ============ */

function parseDType(name: string): DType {
  const mapping: Record<string, DType> = {
    'float16': DType.Float16,
    'float32': DType.Float32,
    'float64': DType.Float64,
    'int8': DType.Int8,
    'int16': DType.Int16,
    'int32': DType.Int32,
    'int64': DType.Int64,
    'uint8': DType.UInt8,
    'uint16': DType.UInt16,
    'uint32': DType.UInt32,
    'uint64': DType.UInt64,
    'bool': DType.Bool,
  };

  const dt = mapping[name.toLowerCase()];
  if (dt === undefined) {
    throw new Error(`Unknown dtype: ${name}`);
  }
  return dt;
}

function getDTypeName(dtype: DType): string {
  const names: Record<DType, string> = {
    [DType.Float64]: 'float64',
    [DType.Float32]: 'float32',
    [DType.Float16]: 'float16',
    [DType.Int64]: 'int64',
    [DType.Int32]: 'int32',
    [DType.Int16]: 'int16',
    [DType.Int8]: 'int8',
    [DType.UInt64]: 'uint64',
    [DType.UInt32]: 'uint32',
    [DType.UInt16]: 'uint16',
    [DType.UInt8]: 'uint8',
    [DType.Bool]: 'bool',
    [DType.Complex64]: 'complex64',
    [DType.Complex128]: 'complex128',
  };
  return names[dtype] ?? 'unknown';
}

/* ============ Exports ============ */

export const limits = {
  finfo,
  iinfo,
};
```

---

### 17.5 Global Configuration

**File:** `src/ts/config.ts` (new file)

```typescript
/**
 * NumJS Global Configuration
 *
 * Provides module-level configuration and diagnostics.
 */

import { get_printoptions, PrintOptions } from './printoptions.js';
import { geterr, ErrorState } from './errstate.js';

/* ============ Version Info ============ */

/** NumJS version */
export const __version__ = '0.1.0';

/** NumJS version as tuple */
export const version_info = [0, 1, 0] as const;

/* ============ Configuration ============ */

interface NumJSConfig {
  /** Enable debug mode */
  debug: boolean;

  /** Enable performance timing */
  timing: boolean;

  /** Prefer WASM acceleration when available */
  useWasm: boolean;

  /** Maximum memory for WASM heap (bytes) */
  maxWasmMemory: number;

  /** Default dtype for array creation */
  defaultDtype: string;

  /** Enable unsafe operations (skip bounds checks) */
  unsafe: boolean;
}

const defaultConfig: NumJSConfig = {
  debug: false,
  timing: false,
  useWasm: true,
  maxWasmMemory: 256 * 1024 * 1024,  // 256 MB
  defaultDtype: 'float64',
  unsafe: false,
};

let currentConfig: NumJSConfig = { ...defaultConfig };

/**
 * Get current configuration.
 */
export function getConfig(): NumJSConfig {
  return { ...currentConfig };
}

/**
 * Set configuration options.
 *
 * @param options - Configuration options to set
 * @returns Previous configuration
 */
export function setConfig(options: Partial<NumJSConfig>): NumJSConfig {
  const old = { ...currentConfig };
  currentConfig = { ...currentConfig, ...options };
  return old;
}

/**
 * Reset configuration to defaults.
 */
export function resetConfig(): void {
  currentConfig = { ...defaultConfig };
}

/* ============ Diagnostics ============ */

interface SystemInfo {
  version: string;
  platform: string;
  wasmSupported: boolean;
  simdSupported: boolean;
  threadsSupported: boolean;
  bigIntSupported: boolean;
  config: NumJSConfig;
  printOptions: PrintOptions;
  errorState: ErrorState;
}

/**
 * Show NumJS configuration and system information.
 */
export function show_config(): void {
  const info = getSystemInfo();

  console.log('NumJS Configuration');
  console.log('===================');
  console.log(`Version: ${info.version}`);
  console.log(`Platform: ${info.platform}`);
  console.log('');
  console.log('Capabilities:');
  console.log(`  WebAssembly: ${info.wasmSupported ? 'Yes' : 'No'}`);
  console.log(`  SIMD: ${info.simdSupported ? 'Yes' : 'No'}`);
  console.log(`  Threads: ${info.threadsSupported ? 'Yes' : 'No'}`);
  console.log(`  BigInt: ${info.bigIntSupported ? 'Yes' : 'No'}`);
  console.log('');
  console.log('Configuration:');
  console.log(`  Debug: ${info.config.debug}`);
  console.log(`  Use WASM: ${info.config.useWasm}`);
  console.log(`  Default dtype: ${info.config.defaultDtype}`);
  console.log(`  Unsafe mode: ${info.config.unsafe}`);
  console.log('');
  console.log('Print Options:');
  console.log(`  Precision: ${info.printOptions.precision}`);
  console.log(`  Threshold: ${info.printOptions.threshold}`);
  console.log(`  Line width: ${info.printOptions.linewidth}`);
  console.log('');
  console.log('Error Handling:');
  console.log(`  Divide: ${info.errorState.divide}`);
  console.log(`  Overflow: ${info.errorState.over}`);
  console.log(`  Underflow: ${info.errorState.under}`);
  console.log(`  Invalid: ${info.errorState.invalid}`);
}

/**
 * Get system information as object.
 */
export function getSystemInfo(): SystemInfo {
  return {
    version: __version__,
    platform: getPlatform(),
    wasmSupported: typeof WebAssembly !== 'undefined',
    simdSupported: checkSimdSupport(),
    threadsSupported: checkThreadsSupport(),
    bigIntSupported: typeof BigInt !== 'undefined',
    config: getConfig(),
    printOptions: get_printoptions(),
    errorState: geterr(),
  };
}

/**
 * Detect current platform.
 */
function getPlatform(): string {
  if (typeof window !== 'undefined') {
    return 'browser';
  }
  if (typeof process !== 'undefined' && process.versions?.node) {
    return `node ${process.versions.node}`;
  }
  if (typeof Deno !== 'undefined') {
    return 'deno';
  }
  if (typeof Bun !== 'undefined') {
    return 'bun';
  }
  return 'unknown';
}

/**
 * Check SIMD support.
 */
function checkSimdSupport(): boolean {
  try {
    return typeof WebAssembly !== 'undefined' &&
           typeof WebAssembly.validate === 'function' &&
           WebAssembly.validate(new Uint8Array([
             0x00, 0x61, 0x73, 0x6d, 0x01, 0x00, 0x00, 0x00,
             0x01, 0x05, 0x01, 0x60, 0x00, 0x01, 0x7b, 0x03,
             0x02, 0x01, 0x00, 0x0a, 0x0a, 0x01, 0x08, 0x00,
             0x41, 0x00, 0xfd, 0x0f, 0x00, 0x00, 0x0b
           ]));
  } catch {
    return false;
  }
}

/**
 * Check SharedArrayBuffer/threads support.
 */
function checkThreadsSupport(): boolean {
  try {
    return typeof SharedArrayBuffer !== 'undefined';
  } catch {
    return false;
  }
}

/* ============ Debug Utilities ============ */

/**
 * Log debug message if debug mode is enabled.
 */
export function debug(message: string, ...args: any[]): void {
  if (currentConfig.debug) {
    console.debug(`[NumJS] ${message}`, ...args);
  }
}

/**
 * Time a function execution if timing is enabled.
 */
export function timeIt<T>(name: string, fn: () => T): T {
  if (!currentConfig.timing) {
    return fn();
  }

  const start = performance.now();
  try {
    return fn();
  } finally {
    const elapsed = performance.now() - start;
    console.log(`[NumJS Timing] ${name}: ${elapsed.toFixed(3)}ms`);
  }
}

/**
 * Assert a condition in debug mode.
 */
export function assert(condition: boolean, message: string): asserts condition {
  if (currentConfig.debug && !condition) {
    throw new Error(`Assertion failed: ${message}`);
  }
}

/* ============ Exports ============ */

export const config = {
  __version__,
  version_info,
  getConfig,
  setConfig,
  resetConfig,
  show_config,
  getSystemInfo,
  debug,
  timeIt,
  assert,
};
```

---

## File Changes Summary

### New Files to Create

```
src/wasm/
└── fperr.h/c          # Floating-point error state (C/WASM)

src/ts/
├── exceptions.ts      # Exception hierarchy
├── warnings.ts        # Warning system
├── errstate.ts        # Floating-point error handling
├── printoptions.ts    # Print options and formatting
├── getlimits.ts       # finfo/iinfo classes
└── config.ts          # Global configuration
```

### Files to Modify

```
src/ts/types.ts
├── Import and re-export exceptions
└── Add error-related type definitions

src/ts/index.ts
├── Export exceptions module
├── Export warnings module
├── Export errstate functions (seterr, geterr, errstate)
├── Export printoptions functions
├── Export finfo, iinfo classes
├── Export config functions
└── Export __version__

scripts/build-wasm.sh
├── Add fperr.c to compilation
└── Add EXPORTED_FUNCTIONS for error handling
```

---

## Build Script Updates

Add to `scripts/build-wasm.sh` EXPORTED_FUNCTIONS:

```bash
# Floating-point error handling
"_fperr_get_state",
"_fperr_clear",
"_fperr_clear_flags",
"_fperr_get_flags",
"_fperr_set_flags",
"_fperr_check",
"_fperr_seterr_all",
"_fperr_seterr_divide",
"_fperr_seterr_over",
"_fperr_seterr_under",
"_fperr_seterr_invalid",
"_fperr_geterr",
"_fperr_check_nan",
"_fperr_check_inf",
"_fperr_check_divide",
"_fperr_check_underflow",
"_fperr_divide",
"_fperr_log",
"_fperr_sqrt",
"_fperr_pow"
```

Add new source file:

```bash
emcc \
    "$SRC_DIR/ndarray.c" \
    "$SRC_DIR/dtype.c" \
    "$SRC_DIR/pairwise_sum.c" \
    "$SRC_DIR/broadcast.c" \
    "$SRC_DIR/indexing.c" \
    "$SRC_DIR/fperr.c" \
    -o "$OUT_DIR/numjs.cjs" \
    ...
```

---

## WasmModule Interface Updates

Add to `src/ts/types.ts`:

```typescript
// Floating-point error handling
_fperr_clear(): void;
_fperr_get_flags(): number;
_fperr_set_flags(flags: number): void;
_fperr_seterr_all(mode: number): void;
_fperr_seterr_divide(mode: number): void;
_fperr_seterr_over(mode: number): void;
_fperr_seterr_under(mode: number): void;
_fperr_seterr_invalid(mode: number): void;
_fperr_geterr(): number;
_fperr_divide(a: number, b: number): number;
_fperr_log(x: number): number;
_fperr_sqrt(x: number): number;
_fperr_pow(base: number, exp: number): number;
```

---

## Implementation Order

```
Phase 17.1: Exception Hierarchy (Week 1)
├── Day 1: Base exceptions (NumJSError, AxisError, ValueError)
├── Day 2: DType and shape errors
├── Day 3: Broadcasting and index errors
├── Day 4: Module-specific errors (LinAlgError, FFTError)
└── Day 5: Warning system implementation

Phase 17.2: Floating-Point Error Handling (Week 2)
├── Day 1: C/WASM error state implementation
├── Day 2: TypeScript error state interface
├── Day 3: seterr/geterr/errstate functions
├── Day 4: Safe math operations with error tracking
└── Day 5: Integration with existing operations

Phase 17.3: Print Options (Week 3)
├── Day 1: PrintOptions interface and defaults
├── Day 2: set_printoptions/get_printoptions
├── Day 3: format_float_positional/scientific
├── Day 4: array2string implementation
├── Day 5: array_repr/array_str + context manager

Phase 17.4: Machine Limits (Week 4)
├── Day 1: finfo class for float types
├── Day 2: iinfo class for integer types
├── Day 3: Edge cases and precision
├── Day 4: Integration with dtype system
└── Day 5: Tests and documentation

Phase 17.5: Global Configuration (Week 5)
├── Day 1: Configuration singleton
├── Day 2: show_config and diagnostics
├── Day 3: Debug and timing utilities
├── Day 4: Integration and exports
└── Day 5: Comprehensive testing
```

---

## Verification Plan

After Phase 17 completion, verify:

```bash
# Build
npm run build

# Run all tests
npm test

# Phase 17 specific tests:

# Exceptions
✓ AxisError normalizes axis correctly
✓ AxisError throws for out-of-bounds axis
✓ BroadcastError formats shapes correctly
✓ LinAlgError is instanceof NumJSError
✓ Exception codes are set correctly

# Warnings
✓ warn() issues warnings to console
✓ filterwarnings('ignore') suppresses warnings
✓ filterwarnings('error') raises exceptions
✓ simplefilter returns restore function
✓ 'once' mode prevents duplicate warnings

# Floating-Point Errors
✓ seterr changes error handling mode
✓ geterr returns current settings
✓ errstate provides context manager
✓ divide by zero triggers error based on mode
✓ overflow detection works correctly
✓ 'call' mode invokes callback

# Print Options
✓ set_printoptions changes formatting
✓ get_printoptions returns current state
✓ printoptions provides context manager
✓ format_float_positional handles edge cases
✓ format_float_scientific formats correctly
✓ array2string respects threshold/edgeitems

# Machine Limits
✓ finfo(float64).eps matches IEEE 754
✓ finfo(float32).max is correct
✓ iinfo(int32).min == -2147483648
✓ iinfo(int32).max == 2147483647
✓ iinfo throws for float types

# Configuration
✓ show_config outputs info without errors
✓ getConfig returns current state
✓ setConfig changes options
✓ resetConfig restores defaults
✓ debug() only logs when enabled
```

---

## Critical Dependencies for Later Phases

Phase 17 completion enables:

- **All Phases**: Consistent error handling throughout codebase
- **Phase 4 (Ufuncs)**: Floating-point error tracking in operations
- **Phase 13 (numpy.linalg)**: LinAlgError for singular matrices
- **Phase 14 (numpy.fft)**: FFTError for transform errors
- **Phase 15 (numpy.random)**: RandomStateError for seed issues

Phase 17 can be implemented independently but should be integrated into existing phases as they are developed.

---

## Usage Examples

### Exception Handling

```typescript
import { AxisError, BroadcastError, LinAlgError } from 'numjs';

// Axis validation
function sumAlongAxis(arr, axis) {
  const normalizedAxis = AxisError.normalizeAxis(axis, arr.ndim);
  // ... perform sum
}

// Catch specific errors
try {
  const result = linalg.inv(singularMatrix);
} catch (e) {
  if (e instanceof LinAlgError) {
    console.log('Matrix is singular');
  }
}
```

### Floating-Point Error Control

```typescript
import { seterr, geterr, errstate } from 'numjs';

// Change error handling
const old = seterr({ divide: 'ignore', over: 'raise' });

// Use context manager
const restore = errstate({ all: 'ignore' });
try {
  // Code that might produce FP errors
} finally {
  restore();
}

// Check current settings
const settings = geterr();
console.log(settings.divide);  // 'ignore'
```

### Print Options

```typescript
import { set_printoptions, array2string, finfo } from 'numjs';

// Set precision
set_printoptions({ precision: 4, threshold: 50 });

// Format array
console.log(array2string(myArray));

// Get machine epsilon
const eps = finfo('float64').eps;
```
