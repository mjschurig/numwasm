Test support (numpy.testing)
Common test support for all numpy test scripts.

This single module should provide all the common functionality for numpy tests in a single location, so that test scripts can just import it and work right away. For background, see the Testing guidelines

Asserts
assert_allclose(actual, desired[, rtol, ...])

Raises an AssertionError if two objects are not equal up to desired tolerance.

assert_array_almost_equal_nulp(x, y[, nulp])

Compare two arrays relatively to their spacing.

assert_array_max_ulp(a, b[, maxulp, dtype])

Check that all items of arrays differ in at most N Units in the Last Place.

assert_array_equal(actual, desired[, ...])

Raises an AssertionError if two array_like objects are not equal.

assert_array_less(x, y[, err_msg, verbose, ...])

Raises an AssertionError if two array_like objects are not ordered by less than.

assert_equal(actual, desired[, err_msg, ...])

Raises an AssertionError if two objects are not equal.

assert_raises(assert_raises)

Fail unless an exception of class exception_class is thrown by callable when invoked with arguments args and keyword arguments kwargs.

assert_raises_regex(exception_class, ...)

Fail unless an exception of class exception_class and with message that matches expected_regexp is thrown by callable when invoked with arguments args and keyword arguments kwargs.

assert_warns(warning_class, \*args, \*\*kwargs)

Fail unless the given callable throws the specified warning.

assert_no_warnings(\*args, \*\*kwargs)

Fail if the given callable produces any warnings.

assert_no_gc_cycles(\*args, \*\*kwargs)

Fail if the given callable produces any reference cycles.

assert_string_equal(actual, desired)

Test if two strings are equal.

Asserts (not recommended)
It is recommended to use one of assert_allclose, assert_array_almost_equal_nulp or assert_array_max_ulp instead of these functions for more consistent floating point comparisons.

assert\_(val[, msg])

Assert that works in release mode.

assert_almost_equal(actual, desired[, ...])

Raises an AssertionError if two items are not equal up to desired precision.

assert_approx_equal(actual, desired[, ...])

Raises an AssertionError if two items are not equal up to significant digits.

assert_array_almost_equal(actual, desired[, ...])

Raises an AssertionError if two objects are not equal up to desired precision.

print_assert_equal(test_string, actual, desired)

Test if two objects are equal, and print an error message if test fails.

Decorators
decorate_methods(cls, decorator[, testmatch])

Apply a decorator to all methods in a class matching a regular expression.

Test running
clear_and_catch_warnings([record, modules])

Context manager that resets warning registry for catching warnings

measure(code_str[, times, label])

Return elapsed time for executing code in the namespace of the caller.

rundocs([filename, raise_on_error])

Run doctests found in the given file.

suppress_warnings([forwarding_rule])

Context manager and decorator doing much the same as warnings.catch_warnings.

Testing custom array containers (numpy.testing.overrides)
These functions can be useful when testing custom array container implementations which make use of **array_ufunc**/**array_function**.

allows_array_function_override(func)

Determine if a Numpy function can be overridden via **array_function**

allows_array_ufunc_override(func)

Determine if a function can be overridden via **array_ufunc**

get_overridable_numpy_ufuncs()

List all numpy ufuncs overridable via **array_ufunc**

get_overridable_numpy_array_functions()

List all numpy functions overridable via **array_function**

Guidelines
Testing guidelines
Introduction
Testing NumPy
Running tests from inside Python
Running tests from the command line
Running doctests
Other methods of running tests
Writing your own tests
Using C code in tests
build_and_import_extension
Labeling tests
Easier setup and teardown functions / methods
Parametric tests
Doctests
tests/
**init**.py and setup.py
Tips & Tricks
Known failures & skipping tests
Tests on random data
Documentation for numpy.test
test
