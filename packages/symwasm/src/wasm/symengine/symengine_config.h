#ifndef SYMENGINE_CONFIG_HPP
#define SYMENGINE_CONFIG_HPP

// SymEngine version
#define SYMENGINE_MAJOR_VERSION 0
#define SYMENGINE_MINOR_VERSION 11
#define SYMENGINE_PATCH_VERSION 0
#define SYMENGINE_VERSION "0.11.0"

// Core features for WASM build
#ifndef WITH_SYMENGINE_RCP
#define WITH_SYMENGINE_RCP
#endif
#ifndef WITH_SYMENGINE_THREAD_SAFE
#define WITH_SYMENGINE_THREAD_SAFE
#endif
#ifndef HAVE_SYMENGINE_GMP
#define HAVE_SYMENGINE_GMP
#endif

// C++11 features (supported by Emscripten)
#define HAVE_DEFAULT_CONSTRUCTORS
#define HAVE_SYMENGINE_NOEXCEPT
#define HAVE_SYMENGINE_IS_CONSTRUCTIBLE
#define HAVE_SYMENGINE_RESERVE
#define HAVE_SYMENGINE_STD_TO_STRING

// RTTI support (disabled to avoid serialization dependencies)
#define HAVE_SYMENGINE_RTTI 0

// Integer class configuration
#define SYMENGINE_GMPXX 0
#define SYMENGINE_PIRANHA 1
#define SYMENGINE_FLINT 2
#define SYMENGINE_GMP 3
#define SYMENGINE_BOOSTMP 4

// Use GMP as integer class
#define SYMENGINE_INTEGER_CLASS SYMENGINE_GMP

// Size of long double
#define SYMENGINE_SIZEOF_LONG_DOUBLE 16

// Noexcept macro
#ifdef HAVE_SYMENGINE_NOEXCEPT
#  define SYMENGINE_NOEXCEPT noexcept
#else
#  define SYMENGINE_NOEXCEPT
#endif

// Export definitions
#ifndef SYMENGINE_EXPORT_H
#define SYMENGINE_EXPORT_H

#ifdef __EMSCRIPTEN__
#  define SYMENGINE_EXPORT __attribute__((visibility("default")))
#  define SYMENGINE_NO_EXPORT __attribute__((visibility("hidden")))
#else
#  define SYMENGINE_EXPORT
#  define SYMENGINE_NO_EXPORT
#endif

#endif // SYMENGINE_EXPORT_H

// Emscripten specific
#ifdef __EMSCRIPTEN__
// Disable features not compatible with WASM
#endif

#endif // SYMENGINE_CONFIG_HPP
