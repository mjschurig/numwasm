Sorting, searching, and counting
Sorting
sort(a[, axis, kind, order, stable])

Return a sorted copy of an array.

lexsort(keys[, axis])

Perform an indirect stable sort using a sequence of keys.

argsort(a[, axis, kind, order, stable])

Returns the indices that would sort an array.

ndarray.sort([axis, kind, order])

Sort an array in-place.

sort_complex(a)

Sort a complex array using the real part first, then the imaginary part.

partition(a, kth[, axis, kind, order])

Return a partitioned copy of an array.

argpartition(a, kth[, axis, kind, order])

Perform an indirect partition along the given axis using the algorithm specified by the kind keyword.

Searching
argmax(a[, axis, out, keepdims])

Returns the indices of the maximum values along an axis.

nanargmax(a[, axis, out, keepdims])

Return the indices of the maximum values in the specified axis ignoring NaNs.

argmin(a[, axis, out, keepdims])

Returns the indices of the minimum values along an axis.

nanargmin(a[, axis, out, keepdims])

Return the indices of the minimum values in the specified axis ignoring NaNs.

argwhere(a)

Find the indices of array elements that are non-zero, grouped by element.

nonzero(a)

Return the indices of the elements that are non-zero.

flatnonzero(a)

Return indices that are non-zero in the flattened version of a.

where(condition, [x, y], /)

Return elements chosen from x or y depending on condition.

searchsorted(a, v[, side, sorter])

Find indices where elements should be inserted to maintain order.

extract(condition, arr)

Return the elements of an array that satisfy some condition.

Counting
count_nonzero(a[, axis, keepdims])

Counts the number of non-zero values in the array a.
