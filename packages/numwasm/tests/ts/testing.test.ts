/**
 * NumJS Testing Module Tests
 *
 * Tests for assertion functions, exception classes, and utilities.
 */

import { describe, it, expect, beforeAll } from 'vitest';
import {
  NDArray,
  loadWasmModule,
  testing,
  AssertionError,
  SkipTest,
  KnownFailureException,
} from 'numjs-wasm';

// Initialize WASM module
beforeAll(async () => {
  await loadWasmModule();
});

describe('Exception Classes', () => {
  it('should create AssertionError with message', () => {
    const error = new AssertionError('test message');
    expect(error.name).toBe('AssertionError');
    expect(error.message).toBe('test message');
    expect(error instanceof Error).toBe(true);
  });

  it('should create SkipTest with default message', () => {
    const error = new SkipTest();
    expect(error.name).toBe('SkipTest');
    expect(error.message).toBe('Test skipped');
  });

  it('should create SkipTest with custom message', () => {
    const error = new SkipTest('prerequisite not met');
    expect(error.message).toBe('prerequisite not met');
  });

  it('should create KnownFailureException with default message', () => {
    const error = new KnownFailureException();
    expect(error.name).toBe('KnownFailureException');
    expect(error.message).toBe('Known failure');
  });

  it('should create KnownFailureException with custom message', () => {
    const error = new KnownFailureException('expected to fail');
    expect(error.message).toBe('expected to fail');
  });
});

describe('assert_', () => {
  it('should pass for truthy values', () => {
    expect(() => testing.assert_(true)).not.toThrow();
    expect(() => testing.assert_(1)).not.toThrow();
    expect(() => testing.assert_('string')).not.toThrow();
    expect(() => testing.assert_([])).not.toThrow(); // empty array is truthy
    expect(() => testing.assert_({})).not.toThrow();
  });

  it('should throw for falsy values', () => {
    expect(() => testing.assert_(false)).toThrow(AssertionError);
    expect(() => testing.assert_(0)).toThrow(AssertionError);
    expect(() => testing.assert_('')).toThrow(AssertionError);
    expect(() => testing.assert_(null)).toThrow(AssertionError);
    expect(() => testing.assert_(undefined)).toThrow(AssertionError);
  });

  it('should include custom message on failure', () => {
    expect(() => testing.assert_(false, 'custom message')).toThrow('custom message');
  });
});

describe('assert_equal', () => {
  it('should pass for equal scalars', () => {
    expect(() => testing.assert_equal(3, 3)).not.toThrow();
    expect(() => testing.assert_equal(-5.5, -5.5)).not.toThrow();
    expect(() => testing.assert_equal(0, 0)).not.toThrow();
  });

  it('should fail for unequal scalars', () => {
    expect(() => testing.assert_equal(3, 4)).toThrow(AssertionError);
    expect(() => testing.assert_equal(0, 1)).toThrow(AssertionError);
  });

  it('should treat NaN as equal to NaN', () => {
    expect(() => testing.assert_equal(NaN, NaN)).not.toThrow();
  });

  it('should pass for equal arrays', () => {
    expect(() => testing.assert_equal([1, 2, 3], [1, 2, 3])).not.toThrow();
    expect(() => testing.assert_equal([[1, 2], [3, 4]], [[1, 2], [3, 4]])).not.toThrow();
  });

  it('should fail for unequal arrays', () => {
    expect(() => testing.assert_equal([1, 2, 3], [1, 2, 4])).toThrow(AssertionError);
    expect(() => testing.assert_equal([1, 2, 3], [1, 2])).toThrow(AssertionError);
  });

  it('should pass for equal strings', () => {
    expect(() => testing.assert_equal('hello', 'hello')).not.toThrow();
  });

  it('should fail for unequal strings', () => {
    expect(() => testing.assert_equal('hello', 'world')).toThrow(AssertionError);
  });
});

describe('assert_array_equal', () => {
  it('should pass for equal arrays', () => {
    expect(() => testing.assert_array_equal([1, 2, 3], [1, 2, 3])).not.toThrow();
  });

  it('should fail for arrays with different shapes', () => {
    expect(() => testing.assert_array_equal([1, 2], [1, 2, 3])).toThrow(AssertionError);
  });

  it('should fail for arrays with different values', () => {
    expect(() => testing.assert_array_equal([1, 2, 3], [1, 2, 4])).toThrow(AssertionError);
  });

  it('should treat NaN as equal to NaN', () => {
    expect(() => testing.assert_array_equal([NaN], [NaN])).not.toThrow();
    expect(() => testing.assert_array_equal([1, NaN, 3], [1, NaN, 3])).not.toThrow();
  });

  it('should handle single element arrays', () => {
    expect(() => testing.assert_array_equal(5, 5)).not.toThrow();
    expect(() => testing.assert_array_equal([5], 5)).not.toThrow();
  });

  it('should work with NDArray', async () => {
    const arr1 = await NDArray.fromArray([1, 2, 3]);
    const arr2 = await NDArray.fromArray([1, 2, 3]);

    expect(() => testing.assert_array_equal(arr1, arr2)).not.toThrow();

    arr1.dispose();
    arr2.dispose();
  });
});

describe('assert_allclose', () => {
  it('should pass for exactly equal arrays', () => {
    expect(() => testing.assert_allclose([1.0, 2.0], [1.0, 2.0])).not.toThrow();
  });

  it('should pass within default tolerance', () => {
    expect(() => testing.assert_allclose([1.0, 2.0], [1.0 + 1e-8, 2.0])).not.toThrow();
  });

  it('should fail outside default tolerance', () => {
    expect(() => testing.assert_allclose([1.0, 2.0], [1.1, 2.1])).toThrow(AssertionError);
  });

  it('should respect custom rtol', () => {
    expect(() => testing.assert_allclose([1.0, 2.0], [1.1, 2.2], 0.2)).not.toThrow();
    expect(() => testing.assert_allclose([1.0, 2.0], [1.1, 2.1], 0.05)).toThrow(AssertionError);
  });

  it('should respect custom atol', () => {
    expect(() => testing.assert_allclose([1.0, 2.0], [1.5, 2.5], 0, 0.6)).not.toThrow();
    expect(() => testing.assert_allclose([1.0, 2.0], [1.5, 2.5], 0, 0.4)).toThrow(AssertionError);
  });

  it('should handle NaN with equal_nan=true', () => {
    expect(() => testing.assert_allclose([NaN], [NaN], 1e-7, 0, true)).not.toThrow();
  });

  it('should fail NaN comparison with equal_nan=false', () => {
    expect(() => testing.assert_allclose([NaN], [NaN], 1e-7, 0, false)).toThrow(AssertionError);
  });

  it('should handle Infinity correctly', () => {
    expect(() => testing.assert_allclose([Infinity], [Infinity])).not.toThrow();
    expect(() => testing.assert_allclose([-Infinity], [-Infinity])).not.toThrow();
    expect(() => testing.assert_allclose([Infinity], [-Infinity])).toThrow(AssertionError);
  });

  it('should fail for different shapes', () => {
    expect(() => testing.assert_allclose([1, 2], [1, 2, 3])).toThrow(AssertionError);
  });
});

describe('assert_almost_equal', () => {
  it('should pass for values equal to specified decimal places', () => {
    expect(() => testing.assert_almost_equal(1.23456789, 1.23456780, 7)).not.toThrow();
  });

  it('should fail for values different beyond decimal places', () => {
    expect(() => testing.assert_almost_equal(1.23456789, 1.23456780, 8)).toThrow(AssertionError);
  });

  it('should use default decimal=7', () => {
    expect(() => testing.assert_almost_equal(1.00000001, 1.0)).not.toThrow();
    expect(() => testing.assert_almost_equal(1.0000001, 1.0)).toThrow(AssertionError);
  });

  it('should handle NaN', () => {
    expect(() => testing.assert_almost_equal(NaN, NaN)).not.toThrow();
  });
});

describe('assert_approx_equal', () => {
  it('should pass for values equal to significant figures', () => {
    expect(() => testing.assert_approx_equal(1234567.0, 1234568.0, 6)).not.toThrow();
  });

  it('should fail for values different beyond significant figures', () => {
    expect(() => testing.assert_approx_equal(1234567.0, 1234568.0, 7)).toThrow(AssertionError);
  });

  it('should handle zero correctly', () => {
    expect(() => testing.assert_approx_equal(0, 0)).not.toThrow();
    expect(() => testing.assert_approx_equal(0.0001, 0)).toThrow(AssertionError);
  });

  it('should handle NaN', () => {
    expect(() => testing.assert_approx_equal(NaN, NaN)).not.toThrow();
  });
});

describe('assert_array_almost_equal', () => {
  it('should pass for arrays almost equal', () => {
    expect(() => testing.assert_array_almost_equal([1.0, 2.33333], [1.0, 2.33334], 4)).not.toThrow();
  });

  it('should fail for arrays not almost equal', () => {
    expect(() => testing.assert_array_almost_equal([1.0, 2.33333], [1.0, 2.33400], 4)).toThrow(AssertionError);
  });

  it('should use default decimal=6', () => {
    expect(() => testing.assert_array_almost_equal([1.0000001], [1.0])).not.toThrow();
    // 0.00001 > 10^(-6), so should fail
    expect(() => testing.assert_array_almost_equal([1.00001], [1.0])).toThrow(AssertionError);
  });
});

describe('assert_array_less', () => {
  it('should pass when all elements are less', () => {
    expect(() => testing.assert_array_less([1, 2, 3], [2, 3, 4])).not.toThrow();
  });

  it('should fail when any element is not less', () => {
    expect(() => testing.assert_array_less([1, 2, 3], [1, 3, 4])).toThrow(AssertionError); // 1 >= 1
    expect(() => testing.assert_array_less([1, 2, 3], [2, 2, 4])).toThrow(AssertionError); // 2 >= 2
  });

  it('should handle scalars', () => {
    expect(() => testing.assert_array_less(1, 2)).not.toThrow();
    expect(() => testing.assert_array_less(2, 2)).toThrow(AssertionError);
  });
});

describe('assert_array_max_ulp', () => {
  it('should pass for equal values (0 ULP)', () => {
    const ulp = testing.assert_array_max_ulp([1.0], [1.0]);
    expect(ulp).toBe(0);
  });

  it('should pass for values within maxulp', () => {
    const a = 1.0;
    const b = 1.0 + Number.EPSILON;
    expect(() => testing.assert_array_max_ulp([a], [b], 1)).not.toThrow();
  });

  it('should return max ULP difference found', () => {
    const ulp = testing.assert_array_max_ulp([1.0, 2.0], [1.0, 2.0], 10);
    expect(ulp).toBe(0);
  });

  it('should handle NaN (returns Infinity)', () => {
    expect(() => testing.assert_array_max_ulp([NaN], [NaN], 0)).toThrow(AssertionError);
  });
});

describe('assert_array_compare', () => {
  it('should pass when comparison function returns true for all', () => {
    const greaterOrEqual = (a: number, b: number) => a >= b;
    expect(() => testing.assert_array_compare(greaterOrEqual, [1, 2, 3], [0, 1, 2])).not.toThrow();
  });

  it('should fail when comparison returns false', () => {
    const lessThan = (a: number, b: number) => a < b;
    expect(() => testing.assert_array_compare(lessThan, [1, 2, 3], [1, 3, 4])).toThrow(AssertionError);
  });

  it('should include custom header in error', () => {
    const equal = (a: number, b: number) => a === b;
    expect(() => testing.assert_array_compare(equal, [1], [2], 'Values should be equal'))
      .toThrow('Values should be equal');
  });
});

describe('assert_raises', () => {
  it('should pass when expected exception is thrown', () => {
    expect(() => testing.assert_raises(Error, () => { throw new Error('test'); })).not.toThrow();
  });

  it('should fail when no exception is thrown', () => {
    expect(() => testing.assert_raises(Error, () => { /* no throw */ })).toThrow(AssertionError);
  });

  it('should fail when wrong exception type is thrown', () => {
    expect(() => testing.assert_raises(TypeError, () => { throw new Error('test'); })).toThrow(AssertionError);
  });

  it('should pass arguments to callable', () => {
    const fn = (a: number, b: number) => {
      if (a > b) throw new Error('a > b');
    };
    expect(() => testing.assert_raises(Error, fn, 5, 3)).not.toThrow();
    expect(() => testing.assert_raises(Error, fn, 1, 3)).toThrow(AssertionError);
  });
});

describe('assert_raises_regex', () => {
  it('should pass when exception and message match', () => {
    expect(() => testing.assert_raises_regex(
      Error,
      /dimension/i,
      () => { throw new Error('Dimension mismatch'); }
    )).not.toThrow();
  });

  it('should fail when message does not match', () => {
    expect(() => testing.assert_raises_regex(
      Error,
      /shape/i,
      () => { throw new Error('Dimension mismatch'); }
    )).toThrow(AssertionError);
  });

  it('should work with string pattern', () => {
    expect(() => testing.assert_raises_regex(
      Error,
      'mismatch',
      () => { throw new Error('Dimension mismatch'); }
    )).not.toThrow();
  });

  it('should fail when no exception is thrown', () => {
    expect(() => testing.assert_raises_regex(Error, /test/, () => {})).toThrow(AssertionError);
  });
});

describe('assert_string_equal', () => {
  it('should pass for equal strings', () => {
    expect(() => testing.assert_string_equal('hello', 'hello')).not.toThrow();
  });

  it('should fail for different strings', () => {
    expect(() => testing.assert_string_equal('hello', 'world')).toThrow(AssertionError);
  });

  it('should show diff position in error', () => {
    expect(() => testing.assert_string_equal('hello', 'hallo'))
      .toThrow(/position 1/);
  });

  it('should handle empty strings', () => {
    expect(() => testing.assert_string_equal('', '')).not.toThrow();
    expect(() => testing.assert_string_equal('', 'a')).toThrow(AssertionError);
  });
});

describe('measure', () => {
  it('should return average time in milliseconds', () => {
    const time = testing.measure(() => {
      // Simple computation
      let sum = 0;
      for (let i = 0; i < 1000; i++) sum += i;
    }, 10);

    expect(typeof time).toBe('number');
    expect(time).toBeGreaterThanOrEqual(0);
  });

  it('should use default times=1', () => {
    const time = testing.measure(() => {});
    expect(time).toBeGreaterThanOrEqual(0);
  });
});

describe('print_assert_equal', () => {
  it('should not throw', () => {
    // Just verify it doesn't throw
    expect(() => testing.print_assert_equal('test', 1, 1)).not.toThrow();
    expect(() => testing.print_assert_equal('test', [1, 2], [1, 2])).not.toThrow();
    expect(() => testing.print_assert_equal('test', 'a', 'b')).not.toThrow();
  });
});

describe('Integration with NDArray', () => {
  it('should compare NDArrays with assert_array_equal', async () => {
    const arr1 = await NDArray.fromArray([[1, 2], [3, 4]]);
    const arr2 = await NDArray.fromArray([[1, 2], [3, 4]]);
    const arr3 = await NDArray.fromArray([[1, 2], [3, 5]]);

    expect(() => testing.assert_array_equal(arr1, arr2)).not.toThrow();
    expect(() => testing.assert_array_equal(arr1, arr3)).toThrow(AssertionError);

    arr1.dispose();
    arr2.dispose();
    arr3.dispose();
  });

  it('should compare NDArrays with assert_allclose', async () => {
    const arr1 = await NDArray.fromArray([1.0, 2.0, 3.0]);
    const arr2 = await NDArray.fromArray([1.0001, 2.0001, 3.0001]);

    expect(() => testing.assert_allclose(arr1, arr2, 1e-3)).not.toThrow();
    expect(() => testing.assert_allclose(arr1, arr2, 1e-5)).toThrow(AssertionError);

    arr1.dispose();
    arr2.dispose();
  });

  it('should compare mixed NDArray and plain array', async () => {
    const arr = await NDArray.fromArray([1, 2, 3]);

    expect(() => testing.assert_array_equal(arr, [1, 2, 3])).not.toThrow();
    expect(() => testing.assert_array_equal([1, 2, 3], arr)).not.toThrow();

    arr.dispose();
  });
});
