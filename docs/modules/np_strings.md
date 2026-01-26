String functionality
The numpy.strings module provides a set of universal functions operating on arrays of type numpy.str* or numpy.bytes*. For example

np.strings.add(["num", "doc"], ["py", "umentation"])
array(['numpy', 'documentation'], dtype='<U13')
These universal functions are also used in numpy.char, which provides the numpy.char.chararray array subclass, in order for those routines to get the performance benefits as well.

Note

Prior to NumPy 2.0, all string functionality was in numpy.char, which only operated on fixed-width strings. That module will not be getting updates and will be deprecated at some point in the future.

String operations
add(x1, x2, /[, out, where, casting, order, ...])

Add arguments element-wise.

center(a, width[, fillchar])

Return a copy of a with its elements centered in a string of length width.

capitalize(a)

Return a copy of a with only the first character of each element capitalized.

decode(a[, encoding, errors])

Calls bytes.decode element-wise.

encode(a[, encoding, errors])

Calls str.encode element-wise.

expandtabs(a[, tabsize])

Return a copy of each string element where all tab characters are replaced by one or more spaces.

ljust(a, width[, fillchar])

Return an array with the elements of a left-justified in a string of length width.

lower(a)

Return an array with the elements converted to lowercase.

lstrip(a[, chars])

For each element in a, return a copy with the leading characters removed.

mod(a, values)

Return (a % i), that is pre-Python 2.6 string formatting (interpolation), element-wise for a pair of array_likes of str or unicode.

multiply(a, i)

Return (a \* i), that is string multiple concatenation, element-wise.

partition(a, sep)

Partition each element in a around sep.

replace(a, old, new[, count])

For each element in a, return a copy of the string with occurrences of substring old replaced by new.

rjust(a, width[, fillchar])

Return an array with the elements of a right-justified in a string of length width.

rpartition(a, sep)

Partition (split) each element around the right-most separator.

rstrip(a[, chars])

For each element in a, return a copy with the trailing characters removed.

strip(a[, chars])

For each element in a, return a copy with the leading and trailing characters removed.

swapcase(a)

Return element-wise a copy of the string with uppercase characters converted to lowercase and vice versa.

title(a)

Return element-wise title cased version of string or unicode.

translate(a, table[, deletechars])

For each element in a, return a copy of the string where all characters occurring in the optional argument deletechars are removed, and the remaining characters have been mapped through the given translation table.

upper(a)

Return an array with the elements converted to uppercase.

zfill(a, width)

Return the numeric string left-filled with zeros.

Comparison
The numpy.strings module also exports the comparison universal functions that can now operate on string arrays as well.

equal(x1, x2, /[, out, where, casting, ...])

Return (x1 == x2) element-wise.

not_equal(x1, x2, /[, out, where, casting, ...])

Return (x1 != x2) element-wise.

greater_equal(x1, x2, /[, out, where, ...])

Return the truth value of (x1 >= x2) element-wise.

less_equal(x1, x2, /[, out, where, casting, ...])

Return the truth value of (x1 <= x2) element-wise.

greater(x1, x2, /[, out, where, casting, ...])

Return the truth value of (x1 > x2) element-wise.

less(x1, x2, /[, out, where, casting, ...])

Return the truth value of (x1 < x2) element-wise.

String information
count(a, sub[, start, end])

Returns an array with the number of non-overlapping occurrences of substring sub in the range [start, end).

endswith(a, suffix[, start, end])

Returns a boolean array which is True where the string element in a ends with suffix, otherwise False.

find(a, sub[, start, end])

For each element, return the lowest index in the string where substring sub is found, such that sub is contained in the range [start, end).

index(a, sub[, start, end])

Like find, but raises ValueError when the substring is not found.

isalnum(x, /[, out, where, casting, order, ...])

Returns true for each element if all characters in the string are alphanumeric and there is at least one character, false otherwise.

isalpha(x, /[, out, where, casting, order, ...])

Returns true for each element if all characters in the data interpreted as a string are alphabetic and there is at least one character, false otherwise.

isdecimal(x, /[, out, where, casting, ...])

For each element, return True if there are only decimal characters in the element.

isdigit(x, /[, out, where, casting, order, ...])

Returns true for each element if all characters in the string are digits and there is at least one character, false otherwise.

islower(x, /[, out, where, casting, order, ...])

Returns true for each element if all cased characters in the string are lowercase and there is at least one cased character, false otherwise.

isnumeric(x, /[, out, where, casting, ...])

For each element, return True if there are only numeric characters in the element.

isspace(x, /[, out, where, casting, order, ...])

Returns true for each element if there are only whitespace characters in the string and there is at least one character, false otherwise.

istitle(x, /[, out, where, casting, order, ...])

Returns true for each element if the element is a titlecased string and there is at least one character, false otherwise.

isupper(x, /[, out, where, casting, order, ...])

Return true for each element if all cased characters in the string are uppercase and there is at least one character, false otherwise.

rfind(a, sub[, start, end])

For each element, return the highest index in the string where substring sub is found, such that sub is contained in the range [start, end).

rindex(a, sub[, start, end])

Like rfind, but raises ValueError when the substring sub is not found.

startswith(a, prefix[, start, end])

Returns a boolean array which is True where the string element in a starts with prefix, otherwise False.

str_len(x, /[, out, where, casting, order, ...])

Returns the length of each element.
