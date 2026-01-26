Record Arrays (numpy.rec)
Record arrays expose the fields of structured arrays as properties.

Most commonly, ndarrays contain elements of a single type, e.g. floats, integers, bools etc. However, it is possible for elements to be combinations of these using structured types, such as:

import numpy as np
a = np.array([(1, 2.0), (1, 2.0)],
dtype=[('x', np.int64), ('y', np.float64)])
a
array([(1, 2.), (1, 2.)], dtype=[('x', '<i8'), ('y', '<f8')])
Here, each element consists of two fields: x (and int), and y (a float). This is known as a structured array. The different fields are analogous to columns in a spread-sheet. The different fields can be accessed as one would a dictionary:

a['x']
array([1, 1])
a['y']
array([2., 2.])
Record arrays allow us to access fields as properties:

ar = np.rec.array(a)
ar.x
array([1, 1])
ar.y
array([2., 2.])
Functions
array(obj[, dtype, shape, offset, strides, ...])

Construct a record array from a wide-variety of objects.

find_duplicate(list)

Find duplication in a list, return a list of duplicated elements

format_parser(formats, names, titles[, ...])

Class to convert formats, names, titles description to a dtype.

fromarrays(arrayList[, dtype, shape, ...])

Create a record array from a (flat) list of arrays

fromfile(fd[, dtype, shape, offset, ...])

Create an array from binary file data

fromrecords(recList[, dtype, shape, ...])

Create a recarray from a list of records in text form.

fromstring(datastring[, dtype, shape, ...])

Create a record array from binary data

Also, the numpy.recarray class and the numpy.record scalar dtype are present in this
