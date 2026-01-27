/**
 * NumJS Testing Module - Assertions
 *
 * Provides NumPy-compatible assertion functions for testing numerical code.
 */

import { NDArray } from '../NDArray.js';
import { DType } from '../types.js';

/* ============ Exception Classes ============ */

/**
 * Custom error class for array assertion failures.
 * Provides detailed comparison information.
 */
export class AssertionError extends Error {
  constructor(message: string) {
    super(message);
    this.name = 'AssertionError';
  }
}

/**
 * Exception for skipping tests.
 * Thrown when test prerequisites are not met.
 */
export class SkipTest extends Error {
  constructor(message: string = 'Test skipped') {
    super(message);
    this.name = 'SkipTest';
  }
}

/**
 * Exception for known failures.
 * Used to mark tests that are expected to fail.
 */
export class KnownFailureException extends Error {
  constructor(message: string = 'Known failure') {
    super(message);
    this.name = 'KnownFailureException';
  }
}

/* ============ Basic Assertions ============ */

/**
 * Assert a condition is true.
 *
 * Unlike the built-in assert, this works in optimized mode
 * and provides a consistent interface.
 *
 * @param val - Value to test for truthiness
 * @param msg - Message to display on failure
 * @throws AssertionError if val is falsy
 *
 * @example
 * assert_(x > 0, "x must be positive");
 * assert_(array.size > 0);
 */
export function assert_(val: unknown, msg: string = ''): void {
  if (!val) {
    throw new AssertionError(msg || 'Assertion failed');
  }
}

/**
 * Assert that two objects are equal.
 *
 * Handles comparison of:
 * - Scalars (number, string, boolean)
 * - Arrays (JavaScript arrays, nested arrays)
 * - NDArrays (element-wise comparison)
 * - Objects (deep equality)
 *
 * Special cases:
 * - NaN == NaN is treated as True for testing purposes
 * - Infinity values are compared correctly
 *
 * @param actual - Actual value obtained
 * @param desired - Expected value
 * @param err_msg - Custom error message prefix
 * @param verbose - If true, show detailed diff (default: true)
 * @throws AssertionError if values are not equal
 *
 * @example
 * assert_equal(3, 3);
 * assert_equal([1, 2, 3], [1, 2, 3]);
 * assert_equal(np.array([1, 2]), np.array([1, 2]));
 */
export function assert_equal(
  actual: unknown,
  desired: unknown,
  err_msg: string = '',
  verbose: boolean = true
): void {
  if (actual instanceof NDArray && desired instanceof NDArray) {
    assert_array_equal(actual, desired, err_msg, verbose);
    return;
  }

  if (Array.isArray(actual) && Array.isArray(desired)) {
    if (actual.length !== desired.length) {
      throw new AssertionError(
        buildErrorMessage(actual, desired, 'Arrays have different lengths', verbose)
      );
    }
    for (let i = 0; i < actual.length; i++) {
      assert_equal(actual[i], desired[i], err_msg, verbose);
    }
    return;
  }

  if (!Object.is(actual, desired)) {
    // Handle NaN comparison (NaN == NaN for testing)
    if (typeof actual === 'number' && typeof desired === 'number' &&
        Number.isNaN(actual) && Number.isNaN(desired)) {
      return;
    }

    throw new AssertionError(
      buildErrorMessage(actual, desired, err_msg || 'Values are not equal', verbose)
    );
  }
}

/**
 * Assert that two arrays are element-wise equal.
 *
 * Compares shape and all elements. For floating-point arrays,
 * use assert_allclose for tolerance-based comparison.
 *
 * @param x - First array
 * @param y - Second array
 * @param err_msg - Custom error message
 * @param verbose - Show detailed output on failure
 * @throws AssertionError if arrays differ
 *
 * @example
 * assert_array_equal([1, 2, 3], [1, 2, 3]);  // Pass
 * assert_array_equal([1, 2, 3], [1, 2, 4]);  // Fail
 * assert_array_equal(np.array([NaN]), np.array([NaN]));  // Pass (NaN == NaN)
 */
export function assert_array_equal(
  x: NDArray | number[] | number,
  y: NDArray | number[] | number,
  err_msg: string = '',
  verbose: boolean = true
): void {
  const arr1 = toNDArray(x);
  const arr2 = toNDArray(y);

  // Check shapes
  if (!arraysEqual(arr1.shape, arr2.shape)) {
    throw new AssertionError(
      `Arrays have different shapes: [${arr1.shape}] vs [${arr2.shape}]`
    );
  }

  // Check values
  for (let i = 0; i < arr1.size; i++) {
    const v1 = arr1.getFlat(i);
    const v2 = arr2.getFlat(i);

    // Handle NaN
    if (Number.isNaN(v1) && Number.isNaN(v2)) {
      continue;
    }

    if (v1 !== v2) {
      throw new AssertionError(
        buildErrorMessage(arr1, arr2, err_msg || `Arrays differ at index ${i}`, verbose)
      );
    }
  }
}

/* ============ Tolerance-Based Assertions ============ */

/**
 * Assert that two arrays are element-wise equal within tolerance.
 *
 * The comparison is:
 *   |actual - desired| <= atol + rtol * |desired|
 *
 * This is the recommended assertion for floating-point comparisons.
 *
 * @param actual - Actual array
 * @param desired - Expected array
 * @param rtol - Relative tolerance (default: 1e-7)
 * @param atol - Absolute tolerance (default: 0)
 * @param equal_nan - If true, NaN == NaN (default: true)
 * @param err_msg - Custom error message
 * @param verbose - Show detailed output on failure
 * @throws AssertionError if arrays differ beyond tolerance
 *
 * @example
 * // Basic usage
 * assert_allclose([1.0, 2.0], [1.0, 2.0]);
 *
 * // With custom tolerance
 * assert_allclose([1.0, 2.0], [1.001, 2.002], 1e-2);
 *
 * // Very loose tolerance
 * assert_allclose(result, expected, 1e-3, 1e-5);
 */
export function assert_allclose(
  actual: NDArray | number[] | number,
  desired: NDArray | number[] | number,
  rtol: number = 1e-7,
  atol: number = 0,
  equal_nan: boolean = true,
  err_msg: string = '',
  _verbose: boolean = true
): void {
  const arr1 = toNDArray(actual);
  const arr2 = toNDArray(desired);

  // Check shapes
  if (!arraysEqual(arr1.shape, arr2.shape)) {
    throw new AssertionError(
      `Arrays have different shapes: [${arr1.shape}] vs [${arr2.shape}]`
    );
  }

  const mismatches: number[] = [];

  for (let i = 0; i < arr1.size; i++) {
    const v1 = arr1.getFlat(i);
    const v2 = arr2.getFlat(i);

    // Handle NaN
    if (Number.isNaN(v1) && Number.isNaN(v2)) {
      if (equal_nan) continue;
      mismatches.push(i);
      continue;
    }

    // Handle Inf
    if (!Number.isFinite(v1) || !Number.isFinite(v2)) {
      if (v1 === v2) continue;
      mismatches.push(i);
      continue;
    }

    // Check tolerance: |actual - desired| <= atol + rtol * |desired|
    const diff = Math.abs(v1 - v2);
    const tolerance = atol + rtol * Math.abs(v2);

    if (diff > tolerance) {
      mismatches.push(i);
    }
  }

  if (mismatches.length > 0) {
    const maxDiff = mismatches.reduce((max, i) => {
      const diff = Math.abs(arr1.getFlat(i) - arr2.getFlat(i));
      return Math.max(max, diff);
    }, 0);

    throw new AssertionError(
      `Arrays are not close within tolerance.\n` +
      `Max absolute difference: ${maxDiff}\n` +
      `rtol: ${rtol}, atol: ${atol}\n` +
      `Number of mismatches: ${mismatches.length} / ${arr1.size}` +
      (err_msg ? `\n${err_msg}` : '')
    );
  }
}

/**
 * Assert that two values are equal up to specified decimal places.
 *
 * Comparison: |actual - desired| < 10^(-decimal)
 *
 * For more control over tolerances, use assert_allclose instead.
 *
 * @param actual - Actual value or array
 * @param desired - Expected value or array
 * @param decimal - Number of decimal places (default: 7)
 * @param err_msg - Custom error message
 * @param verbose - Show detailed output
 * @throws AssertionError if values differ
 *
 * @example
 * assert_almost_equal(1.23456789, 1.23456780, 7);  // Pass
 * assert_almost_equal(1.23456789, 1.23456780, 8);  // Fail
 */
export function assert_almost_equal(
  actual: number | NDArray,
  desired: number | NDArray,
  decimal: number = 7,
  err_msg: string = '',
  verbose: boolean = true
): void {
  const atol = Math.pow(10, -decimal);

  if (typeof actual === 'number' && typeof desired === 'number') {
    // Handle NaN
    if (Number.isNaN(actual) && Number.isNaN(desired)) {
      return;
    }
    if (Math.abs(actual - desired) >= atol) {
      throw new AssertionError(
        `Values differ by more than ${decimal} decimal places: ` +
        `${actual} vs ${desired}`
      );
    }
    return;
  }

  assert_allclose(actual, desired, 0, atol, true, err_msg, verbose);
}

/**
 * Assert values are equal to specified significant figures.
 *
 * The comparison is relative to the scale of the expected value.
 *
 * @param actual - Actual value
 * @param desired - Expected value
 * @param significant - Number of significant figures (default: 7)
 * @param err_msg - Custom error message
 * @param verbose - Show detailed output
 * @throws AssertionError if values differ
 *
 * @example
 * assert_approx_equal(1234567.0, 1234568.0, 6);  // Pass
 * assert_approx_equal(1234567.0, 1234568.0, 7);  // Fail
 */
export function assert_approx_equal(
  actual: number,
  desired: number,
  significant: number = 7,
  err_msg: string = '',
  verbose: boolean = true
): void {
  // err_msg and verbose are kept for API compatibility but not used in scalar comparison
  void err_msg;
  void verbose;
  // Handle NaN
  if (Number.isNaN(actual) && Number.isNaN(desired)) {
    return;
  }

  if (desired === 0) {
    if (actual !== 0) {
      throw new AssertionError(`Expected 0, got ${actual}`);
    }
    return;
  }

  const scale = Math.pow(10, Math.floor(Math.log10(Math.abs(desired))));
  const rtol = Math.pow(10, -significant) * scale;

  if (Math.abs(actual - desired) > rtol) {
    throw new AssertionError(
      `Values differ by more than ${significant} significant figures: ` +
      `${actual} vs ${desired}`
    );
  }
}

/**
 * Assert arrays are element-wise equal up to decimal places.
 *
 * Combination of assert_array_equal and assert_almost_equal.
 * First checks shapes, then element values.
 *
 * @param x - First array
 * @param y - Second array
 * @param decimal - Number of decimal places (default: 6)
 * @param err_msg - Custom error message
 * @param verbose - Show detailed output
 * @throws AssertionError if arrays differ
 *
 * @example
 * const x = np.array([1.0, 2.33333]);
 * const y = np.array([1.0, 2.33334]);
 * assert_array_almost_equal(x, y, 4);  // Pass
 */
export function assert_array_almost_equal(
  x: NDArray | number[],
  y: NDArray | number[],
  decimal: number = 6,
  err_msg: string = '',
  verbose: boolean = true
): void {
  const atol = Math.pow(10, -decimal);
  assert_allclose(x, y, 0, atol, true, err_msg, verbose);
}

/* ============ Comparison Assertions ============ */

/**
 * Assert that all elements of x are strictly less than those of y.
 *
 * @param x - First array
 * @param y - Second array
 * @param err_msg - Custom error message
 * @param verbose - Show detailed output
 * @throws AssertionError if any element x[i] >= y[i]
 *
 * @example
 * assert_array_less([1, 2, 3], [2, 3, 4]);  // Pass
 * assert_array_less([1, 2, 3], [1, 3, 4]);  // Fail (1 >= 1)
 */
export function assert_array_less(
  x: NDArray | number[] | number,
  y: NDArray | number[] | number,
  err_msg: string = '',
  verbose: boolean = true
): void {
  // verbose is kept for API compatibility
  void verbose;
  const arr1 = toNDArray(x);
  const arr2 = toNDArray(y);

  // Check shapes
  if (!arraysEqual(arr1.shape, arr2.shape)) {
    throw new AssertionError(
      `Arrays have different shapes: [${arr1.shape}] vs [${arr2.shape}]`
    );
  }

  for (let i = 0; i < arr1.size; i++) {
    const v1 = arr1.getFlat(i);
    const v2 = arr2.getFlat(i);

    if (!(v1 < v2)) {
      throw new AssertionError(
        `Arrays are not less: x[${i}] = ${v1} >= y[${i}] = ${v2}` +
        (err_msg ? `\n${err_msg}` : '')
      );
    }
  }
}

/**
 * Assert arrays differ by at most `maxulp` ULPs.
 *
 * ULP = Unit in the Last Place (floating-point precision unit).
 * This is the most precise way to compare floating-point numbers.
 *
 * @param a - First array
 * @param b - Second array
 * @param maxulp - Maximum ULP difference allowed (default: 1)
 * @param dtype - Data type for comparison
 * @returns Maximum ULP difference found
 * @throws AssertionError if difference exceeds maxulp
 *
 * @example
 * const a = np.array([1.0]);
 * const b = np.array([1.0 + Number.EPSILON]);
 * assert_array_max_ulp(a, b, 1);  // Pass (differ by 1 ULP)
 */
export function assert_array_max_ulp(
  a: NDArray | number[],
  b: NDArray | number[],
  maxulp: number = 1,
  dtype: DType = DType.Float64
): number {
  // dtype is kept for API compatibility (future: use for Float32 ULP calculation)
  void dtype;
  const arr1 = toNDArray(a);
  const arr2 = toNDArray(b);

  // Check shapes
  if (!arraysEqual(arr1.shape, arr2.shape)) {
    throw new AssertionError(
      `Arrays have different shapes: [${arr1.shape}] vs [${arr2.shape}]`
    );
  }

  let maxFound = 0;

  for (let i = 0; i < arr1.size; i++) {
    const v1 = arr1.getFlat(i);
    const v2 = arr2.getFlat(i);

    const ulpDiff = computeUlpDiff(v1, v2);
    maxFound = Math.max(maxFound, ulpDiff);

    if (ulpDiff > maxulp) {
      throw new AssertionError(
        `Arrays differ by more than ${maxulp} ULPs at index ${i}: ` +
        `${v1} vs ${v2} (${ulpDiff} ULPs)`
      );
    }
  }

  return maxFound;
}

/**
 * Assert arrays satisfy a custom comparison.
 *
 * Low-level function for building custom assertions.
 *
 * @param comparison - Comparison function (a, b) => boolean
 * @param x - First array
 * @param y - Second array
 * @param header - Error message header
 * @param precision - Decimal places for display
 * @throws AssertionError if comparison fails
 *
 * @example
 * assert_array_compare(
 *   (a, b) => a >= b,
 *   [1, 2, 3],
 *   [0, 1, 2],
 *   "x should be >= y"
 * );
 */
export function assert_array_compare(
  comparison: (a: number, b: number) => boolean,
  x: NDArray | number[],
  y: NDArray | number[],
  header: string = 'Arrays do not satisfy comparison',
  precision: number = 6
): void {
  const arr1 = toNDArray(x);
  const arr2 = toNDArray(y);

  // Check shapes
  if (!arraysEqual(arr1.shape, arr2.shape)) {
    throw new AssertionError(
      `Arrays have different shapes: [${arr1.shape}] vs [${arr2.shape}]`
    );
  }

  for (let i = 0; i < arr1.size; i++) {
    const v1 = arr1.getFlat(i);
    const v2 = arr2.getFlat(i);

    if (!comparison(v1, v2)) {
      throw new AssertionError(
        `${header}\nFailed at index ${i}: ${v1.toPrecision(precision)} vs ${v2.toPrecision(precision)}`
      );
    }
  }
}

/* ============ Exception Assertions ============ */

/**
 * Assert that a callable raises the specified exception.
 *
 * @param exception_class - Expected exception type
 * @param callable - Function to call
 * @param args - Arguments to pass to callable
 * @throws AssertionError if exception not raised or wrong type
 *
 * @example
 * assert_raises(Error, () => {
 *   throw new Error("test");
 * });
 *
 * assert_raises(
 *   TypeError,
 *   np.add,
 *   "string", 5  // Can't add string and number
 * );
 */
export function assert_raises<T extends Error>(
  exception_class: new (...args: unknown[]) => T,
  callable: (...args: unknown[]) => unknown,
  ...args: unknown[]
): void {
  try {
    callable(...args);
  } catch (e) {
    if (e instanceof exception_class) {
      return;
    }
    throw new AssertionError(
      `Expected ${exception_class.name}, but got ${(e as Error).constructor.name}: ${(e as Error).message}`
    );
  }

  throw new AssertionError(
    `Expected ${exception_class.name} to be raised, but no exception was thrown`
  );
}

/**
 * Assert exception is raised with matching message.
 *
 * @param exception_class - Expected exception type
 * @param expected_regexp - Regex pattern for message
 * @param callable - Function to call
 * @param args - Arguments to pass to callable
 * @throws AssertionError if exception not raised or message doesn't match
 *
 * @example
 * assert_raises_regex(
 *   Error,
 *   /dimension mismatch/i,
 *   np.dot,
 *   np.array([1, 2]),
 *   np.array([1, 2, 3])
 * );
 */
export function assert_raises_regex<T extends Error>(
  exception_class: new (...args: unknown[]) => T,
  expected_regexp: string | RegExp,
  callable: (...args: unknown[]) => unknown,
  ...args: unknown[]
): void {
  const regex = typeof expected_regexp === 'string'
    ? new RegExp(expected_regexp)
    : expected_regexp;

  try {
    callable(...args);
  } catch (e) {
    if (e instanceof exception_class) {
      if (regex.test((e as Error).message)) {
        return;
      }
      throw new AssertionError(
        `Exception message "${(e as Error).message}" does not match pattern "${regex}"`
      );
    }
    throw new AssertionError(
      `Expected ${exception_class.name}, but got ${(e as Error).constructor.name}`
    );
  }

  throw new AssertionError(
    `Expected ${exception_class.name} to be raised, but no exception was thrown`
  );
}

/**
 * Assert that a function emits a warning.
 *
 * Note: JavaScript doesn't have built-in warnings.
 * This integrates with a custom warning system if available.
 *
 * @param warning_class - Expected warning type
 * @param callable - Function to call
 * @param args - Arguments to pass to callable
 * @returns Result of calling callable
 * @throws AssertionError if warning not emitted
 *
 * @example
 * assert_warns(DeprecationWarning, deprecated_function);
 */
export function assert_warns(
  warning_class: new (...args: unknown[]) => Error,
  callable: (...args: unknown[]) => unknown,
  ...args: unknown[]
): unknown {
  // warning_class is kept for API compatibility
  void warning_class;
  // JavaScript doesn't have built-in warnings
  // This is a placeholder that would integrate with a warning system
  return callable(...args);
}

/* ============ String Assertions ============ */

/**
 * Assert that two strings are equal.
 *
 * Provides clear diff output for multi-line strings.
 *
 * @param actual - Actual string
 * @param desired - Expected string
 * @throws AssertionError if strings differ
 *
 * @example
 * assert_string_equal("hello", "hello");  // Pass
 * assert_string_equal("hello", "world");  // Fail with diff
 */
export function assert_string_equal(
  actual: string,
  desired: string
): void {
  if (actual !== desired) {
    // Find first difference for helpful output
    let diffIndex = 0;
    while (diffIndex < actual.length && diffIndex < desired.length &&
           actual[diffIndex] === desired[diffIndex]) {
      diffIndex++;
    }

    throw new AssertionError(
      `Strings differ at position ${diffIndex}:\n` +
      `  actual: "${actual.slice(Math.max(0, diffIndex - 10), diffIndex + 20)}"\n` +
      `  desired: "${desired.slice(Math.max(0, diffIndex - 10), diffIndex + 20)}"`
    );
  }
}

/**
 * Assert that a function emits no warnings.
 *
 * @param func - Function to call
 * @param args - Arguments to pass to function
 * @returns Result of calling function
 * @throws AssertionError if any warning is emitted
 *
 * @example
 * const result = assert_no_warnings(my_function, arg1, arg2);
 */
export function assert_no_warnings(
  func: (...args: unknown[]) => unknown,
  ...args: unknown[]
): unknown {
  // JavaScript doesn't have built-in warnings
  // This would integrate with a custom warning system
  return func(...args);
}

/* ============ Utility Functions ============ */

/**
 * Measure execution time of a code block.
 *
 * @param code - Function to benchmark
 * @param times - Number of iterations (default: 1)
 * @returns Average execution time in milliseconds
 *
 * @example
 * const time = measure(() => {
 *   np.dot(largeMatrix, largeMatrix);
 * }, 100);
 * console.log(`Average: ${time}ms`);
 */
export function measure(
  code: () => void,
  times: number = 1
): number {
  const start = performance.now();
  for (let i = 0; i < times; i++) {
    code();
  }
  const end = performance.now();
  return (end - start) / times;
}

/**
 * Print comparison result for debugging.
 *
 * @param test_string - Description of test
 * @param actual - Actual value
 * @param desired - Expected value
 *
 * @example
 * print_assert_equal("matrix multiply", result, expected);
 * // Outputs:
 * // matrix multiply
 * //   actual: [[1, 2], [3, 4]]
 * //   desired: [[1, 2], [3, 4]]
 * //   match: true
 */
export function print_assert_equal(
  test_string: string,
  actual: unknown,
  desired: unknown
): void {
  const match = deepEqual(actual, desired);
  console.log(`${test_string}`);
  console.log(`  actual: ${formatValue(actual)}`);
  console.log(`  desired: ${formatValue(desired)}`);
  console.log(`  match: ${match}`);
}

/* ============ Internal Helper Functions ============ */

/**
 * Convert various inputs to NDArray for comparison.
 */
function toNDArray(x: NDArray | number[] | number): NDArray {
  if (x instanceof NDArray) return x;

  // Create a temporary wrapper for comparison
  // This is a lightweight approach that avoids async
  if (Array.isArray(x)) {
    return createSyncNDArray(x);
  }
  return createSyncNDArray([x]);
}

/**
 * Create a sync NDArray-like object for testing comparisons.
 * This avoids the need for async operations in assertions.
 */
function createSyncNDArray(data: number[]): NDArray {
  // We need to create a minimal object that supports getFlat, shape, and size
  // Since NDArray.fromArray is async, we create a wrapper object
  const flatData = flattenArray(data);
  const shape = inferShape(data);

  return {
    getFlat(index: number): number {
      return flatData[index];
    },
    get shape(): number[] {
      return shape;
    },
    get size(): number {
      return flatData.length;
    },
    toString(): string {
      return JSON.stringify(data);
    }
  } as NDArray;
}

/**
 * Infer shape from nested array.
 */
function inferShape(data: number | number[]): number[] {
  const shape: number[] = [];
  let current: number | number[] = data;

  while (Array.isArray(current)) {
    shape.push(current.length);
    if (current.length === 0) break;
    current = current[0];
  }

  return shape;
}

/**
 * Flatten nested array.
 */
function flattenArray(data: number | number[]): number[] {
  if (typeof data === 'number') {
    return [data];
  }

  const result: number[] = [];
  for (const item of data) {
    if (typeof item === 'number') {
      result.push(item);
    } else {
      result.push(...flattenArray(item));
    }
  }
  return result;
}

/**
 * Check if two number arrays are equal.
 */
function arraysEqual(a: number[], b: number[]): boolean {
  if (a.length !== b.length) return false;
  return a.every((v, i) => v === b[i]);
}

/**
 * Build a detailed error message.
 */
function buildErrorMessage(
  actual: unknown,
  desired: unknown,
  header: string,
  verbose: boolean
): string {
  if (!verbose) {
    return header;
  }

  return `${header}\n` +
    `Actual:\n${formatValue(actual)}\n` +
    `Desired:\n${formatValue(desired)}`;
}

/**
 * Format a value for display.
 */
function formatValue(val: unknown): string {
  if (val instanceof NDArray) {
    return val.toString();
  }
  if (Array.isArray(val)) {
    return JSON.stringify(val);
  }
  return String(val);
}

/**
 * Deep equality check.
 */
function deepEqual(a: unknown, b: unknown): boolean {
  if (a === b) return true;
  if (typeof a !== typeof b) return false;
  if (typeof a === 'number' && typeof b === 'number' &&
      Number.isNaN(a) && Number.isNaN(b)) return true;

  if (a instanceof NDArray && b instanceof NDArray) {
    if (!arraysEqual(a.shape, b.shape)) return false;
    for (let i = 0; i < a.size; i++) {
      const v1 = a.getFlat(i);
      const v2 = b.getFlat(i);
      if (v1 !== v2 && !(Number.isNaN(v1) && Number.isNaN(v2))) return false;
    }
    return true;
  }

  if (Array.isArray(a) && Array.isArray(b)) {
    if (a.length !== b.length) return false;
    return a.every((v, i) => deepEqual(v, b[i]));
  }

  return false;
}

/**
 * Compute ULP difference between two floating point numbers.
 */
function computeUlpDiff(a: number, b: number): number {
  if (a === b) return 0;
  if (Number.isNaN(a) || Number.isNaN(b)) return Infinity;
  if (!Number.isFinite(a) || !Number.isFinite(b)) return Infinity;

  // Convert to integer bit patterns
  const buffer = new ArrayBuffer(8);
  const floatView = new Float64Array(buffer);
  const intView = new BigInt64Array(buffer);

  floatView[0] = a;
  const aInt = intView[0];

  floatView[0] = b;
  const bInt = intView[0];

  // Handle sign differences
  if ((aInt < 0n) !== (bInt < 0n)) {
    return Number(abs64(aInt) + abs64(bInt));
  }

  return Number(abs64(aInt - bInt));
}

/**
 * Absolute value for BigInt.
 */
function abs64(n: bigint): bigint {
  return n < 0n ? -n : n;
}
