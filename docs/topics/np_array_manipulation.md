Array manipulation routines
Basic operations
copyto(dst, src[, casting, where])

Copies values from one array to another, broadcasting as necessary.

ndim(a)

Return the number of dimensions of an array.

shape(a)

Return the shape of an array.

size(a[, axis])

Return the number of elements along a given axis.

Changing array shape
reshape(a, /[, shape, newshape, order, copy])

Gives a new shape to an array without changing its data.

ravel(a[, order])

Return a contiguous flattened array.

ndarray.flat

A 1-D iterator over the array.

ndarray.flatten([order])

Return a copy of the array collapsed into one dimension.

Transpose-like operations
moveaxis(a, source, destination)

Move axes of an array to new positions.

rollaxis(a, axis[, start])

Roll the specified axis backwards, until it lies in a given position.

swapaxes(a, axis1, axis2)

Interchange two axes of an array.

ndarray.T

View of the transposed array.

transpose(a[, axes])

Returns an array with axes transposed.

permute_dims(a[, axes])

Returns an array with axes transposed.

matrix_transpose(x, /)

Transposes a matrix (or a stack of matrices) x.

Changing number of dimensions
atleast_1d(\*arys)

Convert inputs to arrays with at least one dimension.

atleast_2d(\*arys)

View inputs as arrays with at least two dimensions.

atleast_3d(\*arys)

View inputs as arrays with at least three dimensions.

broadcast

Produce an object that mimics broadcasting.

broadcast_to(array, shape[, subok])

Broadcast an array to a new shape.

broadcast_arrays(\*args[, subok])

Broadcast any number of arrays against each other.

expand_dims(a, axis)

Expand the shape of an array.

squeeze(a[, axis])

Remove axes of length one from a.

Changing kind of array
asarray(a[, dtype, order, device, copy, like])

Convert the input to an array.

asanyarray(a[, dtype, order, like])

Convert the input to an ndarray, but pass ndarray subclasses through.

asmatrix(data[, dtype])

Interpret the input as a matrix.

asfortranarray(a[, dtype, like])

Return an array (ndim >= 1) laid out in Fortran order in memory.

ascontiguousarray(a[, dtype, like])

Return a contiguous array (ndim >= 1) in memory (C order).

asarray_chkfinite(a[, dtype, order])

Convert the input to an array, checking for NaNs or Infs.

require(a[, dtype, requirements, like])

Return an ndarray of the provided type that satisfies requirements.

Joining arrays
concatenate([axis, out, dtype, casting])

Join a sequence of arrays along an existing axis.

concat([axis, out, dtype, casting])

Join a sequence of arrays along an existing axis.

stack(arrays[, axis, out, dtype, casting])

Join a sequence of arrays along a new axis.

block(arrays)

Assemble an nd-array from nested lists of blocks.

vstack(tup, \*[, dtype, casting])

Stack arrays in sequence vertically (row wise).

hstack(tup, \*[, dtype, casting])

Stack arrays in sequence horizontally (column wise).

dstack(tup)

Stack arrays in sequence depth wise (along third axis).

column_stack(tup)

Stack 1-D arrays as columns into a 2-D array.

Splitting arrays
split(ary, indices_or_sections[, axis])

Split an array into multiple sub-arrays as views into ary.

array_split(ary, indices_or_sections[, axis])

Split an array into multiple sub-arrays.

dsplit(ary, indices_or_sections)

Split array into multiple sub-arrays along the 3rd axis (depth).

hsplit(ary, indices_or_sections)

Split an array into multiple sub-arrays horizontally (column-wise).

vsplit(ary, indices_or_sections)

Split an array into multiple sub-arrays vertically (row-wise).

unstack(x, /, \*[, axis])

Split an array into a sequence of arrays along the given axis.

Tiling arrays
tile(A, reps)

Construct an array by repeating A the number of times given by reps.

repeat(a, repeats[, axis])

Repeat each element of an array after themselves

Adding and removing elements
delete(arr, obj[, axis])

Return a new array with sub-arrays along an axis deleted.

insert(arr, obj, values[, axis])

Insert values along the given axis before the given indices.

append(arr, values[, axis])

Append values to the end of an array.

resize(a, new_shape)

Return a new array with the specified shape.

trim_zeros(filt[, trim])

Trim the leading and/or trailing zeros from a 1-D array or sequence.

unique(ar[, return_index, return_inverse, ...])

Find the unique elements of an array.

pad(array, pad_width[, mode])

Pad an array.

Rearranging elements
flip(m[, axis])

Reverse the order of elements in an array along the given axis.

fliplr(m)

Reverse the order of elements along axis 1 (left/right).

flipud(m)

Reverse the order of elements along axis 0 (up/down).

reshape(a, /[, shape, newshape, order, copy])

Gives a new shape to an array without changing its data.

roll(a, shift[, axis])

Roll array elements along a given axis.

rot90(m[, k, axes])

Rotate an array by 90 degrees in the plane specified by axes.
