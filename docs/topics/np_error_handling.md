Floating point error handling
Error handling settings are stored in contextvars allowing different threads or async tasks to have independent configurations. For more information, see Thread Safety.

How numpy handles numerical exceptions
The default is to 'warn' for invalid, divide, and overflow and 'ignore' for underflow. But this can be changed, and it can be set individually for different kinds of exceptions. The different behaviors are:

'ignore' : Take no action when the exception occurs.

'warn' : Print a RuntimeWarning (via the Python warnings module).

'raise' : Raise a FloatingPointError.

'call' : Call a specified function.

'print' : Print a warning directly to stdout.

'log' : Record error in a Log object.

These behaviors can be set for all kinds of errors or specific ones:

all : apply to all numeric exceptions

invalid : when NaNs are generated

divide : divide by zero (for integers as well!)

overflow : floating point overflows

underflow : floating point underflows

Note that integer divide-by-zero is handled by the same machinery.

The error handling mode can be configured numpy.errstate context manager.

Examples
with np.errstate(all='warn'):
np.zeros(5, dtype=np.float32) / 0.0
<python-input-1>:2: RuntimeWarning: invalid value encountered in divide
array([nan, nan, nan, nan, nan], dtype=float32)
with np.errstate(under='ignore'):
np.array([1.e-100])\*\*10
array([0.])
with np.errstate(invalid='raise'):
np.sqrt(np.array([-1.]))

Traceback (most recent call last):
File "<python-input-1>", line 2, in <module>
np.sqrt(np.array([-1.]))
~~~~~~~^^^^^^^^^^^^^^^^^
FloatingPointError: invalid value encountered in sqrt
def errorhandler(errstr, errflag):
print("saw stupid error!")
with np.errstate(call=errorhandler, all='call'):
np.zeros(5, dtype=np.int32) / 0
saw stupid error!
array([nan, nan, nan, nan, nan])
Setting and getting error handling
seterr([all, divide, over, under, invalid])

Set how floating-point errors are handled.

geterr()

Get the current way of handling floating-point errors.

seterrcall(func)

Set the floating-point error callback function or log object.

geterrcall()

Return the current callback function used on floating-point errors.

errstate(\*\*kwargs)

Context manager for floating-point error handling.
