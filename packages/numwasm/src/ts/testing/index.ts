/**
 * NumJS Testing Module
 *
 * Provides NumPy-compatible assertion functions for testing numerical code.
 *
 * @module testing
 *
 * @example
 * import { testing } from 'numjs';
 *
 * // Basic assertions
 * testing.assert_equal(a, b);
 * testing.assert_array_equal(arr1, arr2);
 *
 * // Tolerance-based assertions
 * testing.assert_allclose(actual, expected, 1e-5);
 * testing.assert_almost_equal(result, 3.14159, 5);
 *
 * // Exception assertions
 * testing.assert_raises(Error, () => { throw new Error(); });
 */
export {
  // Exception classes
  AssertionError,
  SkipTest,
  KnownFailureException,

  // Basic assertions
  assert_,
  assert_equal,
  assert_array_equal,

  // Tolerance-based assertions
  assert_allclose,
  assert_almost_equal,
  assert_approx_equal,
  assert_array_almost_equal,

  // Comparison assertions
  assert_array_less,
  assert_array_max_ulp,
  assert_array_compare,

  // Exception assertions
  assert_raises,
  assert_raises_regex,
  assert_warns,

  // String assertions
  assert_string_equal,
  assert_no_warnings,

  // Utilities
  measure,
  print_assert_equal,
} from './assertions.js';
