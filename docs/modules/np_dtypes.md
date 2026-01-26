Data type classes (numpy.dtypes)
This module is home to specific dtypes related functionality and their classes. For more general information about dtypes, also see numpy.dtype and Data type objects (dtype).

Similar to the builtin types module, this submodule defines types (classes) that are not widely used directly.

New in version NumPy: 1.25

The dtypes module is new in NumPy 1.25. Previously DType classes were only accessible indirectly.

DType classes
The following are the classes of the corresponding NumPy dtype instances and NumPy scalar types. The classes can be used in isinstance checks and can also be instantiated or used directly. Direct use of these classes is not typical, since their scalar counterparts (e.g. np.float64) or strings like "float64" can be used.

Boolean
numpy.dtypes.BoolDType[source]
Bit-sized integers
numpy.dtypes.Int8DType[source]
numpy.dtypes.UInt8DType
numpy.dtypes.Int16DType
numpy.dtypes.UInt16DType
numpy.dtypes.Int32DType
numpy.dtypes.UInt32DType
numpy.dtypes.Int64DType
numpy.dtypes.UInt64DType
C-named integers (may be aliases)
numpy.dtypes.ByteDType[source]
numpy.dtypes.UByteDType
numpy.dtypes.ShortDType
numpy.dtypes.UShortDType
numpy.dtypes.IntDType
numpy.dtypes.UIntDType
numpy.dtypes.LongDType
numpy.dtypes.ULongDType
numpy.dtypes.LongLongDType
numpy.dtypes.ULongLongDType
Floating point
numpy.dtypes.Float16DType[source]
numpy.dtypes.Float32DType
numpy.dtypes.Float64DType
numpy.dtypes.LongDoubleDType
Complex
numpy.dtypes.Complex64DType[source]
numpy.dtypes.Complex128DType
numpy.dtypes.CLongDoubleDType
Strings and Bytestrings
numpy.dtypes.StrDType[source]
numpy.dtypes.BytesDType
numpy.dtypes.StringDType
Times
numpy.dtypes.DateTime64DType[source]
numpy.dtypes.TimeDelta64DType
Others
numpy.dtypes.ObjectDType[source]
numpy.dtypes.VoidDType
