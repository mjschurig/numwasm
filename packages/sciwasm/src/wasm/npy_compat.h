/**
 * NumPy compatibility shim for Emscripten/WASM compilation.
 *
 * Provides type aliases and macros so that scipy's sparsetools C++ headers
 * can compile without numpy installed. Maps npy_ types to standard C/C++ types.
 */
#ifndef NPY_COMPAT_H
#define NPY_COMPAT_H

#include <cstdint>
#include <cstddef>
#include <climits>

/* Index types */
typedef int32_t  npy_int32;
typedef int64_t  npy_int64;
typedef intptr_t npy_intp;

/* Integer types */
typedef int8_t   npy_byte;
typedef uint8_t  npy_ubyte;
typedef int16_t  npy_short;
typedef uint16_t npy_ushort;
typedef int32_t  npy_int;
typedef uint32_t npy_uint;
typedef long     npy_long;
typedef unsigned long npy_ulong;
typedef long long npy_longlong;
typedef unsigned long long npy_ulonglong;

/* Floating types */
typedef float       npy_float;
typedef double      npy_double;
typedef long double npy_longdouble;

/* Boolean type */
typedef unsigned char npy_bool;

/* Size constants */
#define NPY_MAX_INTP INTPTR_MAX
#define NPY_SIZEOF_LONGDOUBLE sizeof(long double)
#define NPY_SIZEOF_DOUBLE sizeof(double)

/* Complex types - simple struct-based (no numpy dependency) */
typedef struct { float real; float imag; } npy_cfloat;
typedef struct { double real; double imag; } npy_cdouble;
typedef struct { long double real; long double imag; } npy_clongdouble;

/* Complex access macros */
#define npy_crealf(z)  ((z).real)
#define npy_cimagf(z)  ((z).imag)
#define npy_creal(z)   ((z).real)
#define npy_cimag(z)   ((z).imag)
#define npy_creall(z)  ((z).real)
#define npy_cimagl(z)  ((z).imag)

#define NPY_CSETREALF(z, v) ((z)->real = (v))
#define NPY_CSETIMAGF(z, v) ((z)->imag = (v))
#define NPY_CSETREAL(z, v)  ((z)->real = (v))
#define NPY_CSETIMAG(z, v)  ((z)->imag = (v))
#define NPY_CSETREALL(z, v) ((z)->real = (v))
#define NPY_CSETIMAGL(z, v) ((z)->imag = (v))

#endif /* NPY_COMPAT_H */
