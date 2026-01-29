#ifndef TEUCHOS_CONFIG_H
#define TEUCHOS_CONFIG_H

// Minimal Teuchos configuration for SymEngine WASM build

// Platform detection
#ifdef __EMSCRIPTEN__
#define __linux__
#endif

// C++ standard features
#define HAVE_TEUCHOS_CXX11
#define HAVE_STD_TYPE_TRAITS

// Threading support
#ifdef WITH_SYMENGINE_THREAD_SAFE
#define HAVE_TEUCHOS_THREAD_SAFE
#define HAVE_TEUCHOS_C99_STDINT_H
#define HAVE_TEUCHOS_STDINT_H
#endif

// Ordinal type (used for array indexing)
#define TEUCHOS_ORDINAL_TYPE int

// Deprecation support (empty for now)
#define TEUCHOS_DEPRECATED

// Disable features not needed for WASM
// #define HAVE_TEUCHOS_DEBUG
// #define HAVE_TEUCHOS_EXTENDED

#endif // TEUCHOS_CONFIG_H
