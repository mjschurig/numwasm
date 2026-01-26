Miscellaneous routines
Performance tuning
setbufsize(size)

Set the size of the buffer used in ufuncs.

getbufsize()

Return the size of the buffer used in ufuncs.

Memory ranges
shares_memory(a, b, /[, max_work])

Determine if two arrays share memory.

may_share_memory(a, b, /[, max_work])

Determine if two arrays might share memory

Utility
get_include()

Return the directory that contains the NumPy \*.h header files.

show_config([mode])

Show libraries and system information on which NumPy was built and is being used

show_runtime()

Print information about various resources in the system including available intrinsic support and BLAS/LAPACK library in use

broadcast_shapes(\*args)

Broadcast the input shapes into a single shape.

NumPy-specific help function
info([object, maxwidth, output, toplevel])

Get help information for an array, function, class, or module.
