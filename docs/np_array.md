Array objects
NumPy provides an N-dimensional array type, the ndarray, which describes a collection of “items” of the same type. The items can be indexed using for example N integers.

All ndarrays are homogeneous: every item takes up the same size block of memory, and all blocks are interpreted in exactly the same way. How each item in the array is to be interpreted is specified by a separate data-type object, one of which is associated with every array. In addition to basic types (integers, floats, etc.), the data type objects can also represent data structures.

An item extracted from an array, e.g., by indexing, is represented by a Python object whose type is one of the array scalar types built in NumPy. The array scalars allow easy manipulation of also more complicated arrangements of data.

../\_images/threefundamental.png
Figure Conceptual diagram showing the relationship between the three fundamental objects used to describe the data in an array: 1) the ndarray itself, 2) the data-type object that describes the layout of a single fixed-size element of the array, 3) the array-scalar Python object that is returned when a single element of the array is accessed.
The N-dimensional array (ndarray)
Constructing arrays
Indexing arrays
Internal memory layout of an ndarray
Array attributes
Array methods
Arithmetic, matrix multiplication, and comparison operations
Special methods
Scalars
Built-in scalar types
Attributes
Indexing
Methods
Defining new types
Data type objects (dtype)
Specifying and constructing data types
Checking the data type
dtype
Data type promotion in NumPy
Detailed behavior of Python scalars
Numerical promotion
Exceptions to the general promotion rules
Promotion of non-numerical datatypes
Details of promoted dtype instances
Iterating over arrays
Single array iteration
Broadcasting array iteration
Putting the inner loop in Cython
Standard array subclasses
Special attributes and methods
Matrix objects
Memory-mapped file arrays
Character arrays (numpy.char)
Record arrays
Masked arrays (numpy.ma)
Standard container class
Array iterators
Masked arrays
The numpy.ma module
Using numpy.ma
Examples
Constants of the numpy.ma module
The MaskedArray class
MaskedArray methods
Masked array operations
The array interface protocol
Python side
C-struct access
Type description examples
Differences with array interface (version 2)
Datetimes and timedeltas
Datetime64 conventions and assumptions
Basic datetimes
Datetime and timedelta arithmetic
Datetime units
Business day functionality
Datetime64 shortcomings
