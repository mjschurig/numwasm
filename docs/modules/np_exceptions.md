Exceptions and Warnings (numpy.exceptions)
General exceptions used by NumPy. Note that some exceptions may be module specific, such as linear algebra errors.

New in version NumPy: 1.25

The exceptions module is new in NumPy 1.25. Older exceptions remain available through the main NumPy namespace for compatibility.

Warnings
ComplexWarning

The warning raised when casting a complex dtype to a real dtype.

VisibleDeprecationWarning

Visible deprecation warning.

RankWarning

Matrix rank warning.

Exceptions
AxisError(axis[, ndim, msg_prefix])

Axis supplied was invalid.

DTypePromotionError

Multiple DTypes could not be converted to a common one.

TooHardError

max_work was exceeded.
