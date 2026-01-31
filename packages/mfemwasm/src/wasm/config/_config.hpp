// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#ifndef MFEM_CONFIG_HEADER
#define MFEM_CONFIG_HEADER

// ============================================================================
// MFEM WebAssembly Configuration
// Serial-only build with minimal dependencies
// ============================================================================

// MFEM version: integer of the form: (major*100 + minor)*100 + patch.
#define MFEM_VERSION 40901

// MFEM version string of the form "3.3" or "3.3.1".
#define MFEM_VERSION_STRING "4.9.1"

// MFEM version type, see the MFEM_VERSION_TYPE_* constants below.
#define MFEM_VERSION_TYPE ((MFEM_VERSION)%2)

// MFEM version type constants.
#define MFEM_VERSION_TYPE_RELEASE 0
#define MFEM_VERSION_TYPE_DEVELOPMENT 1

// Separate MFEM version numbers for major, minor, and patch.
#define MFEM_VERSION_MAJOR ((MFEM_VERSION)/10000)
#define MFEM_VERSION_MINOR (((MFEM_VERSION)/100)%100)
#define MFEM_VERSION_PATCH ((MFEM_VERSION)%100)

// MFEM source directory.
#define MFEM_SOURCE_DIR ""

// MFEM install directory.
#define MFEM_INSTALL_DIR ""

// Git string for version info (not tracked in WASM build)
#define MFEM_GIT_STRING ""

// ============================================================================
// Core Settings - ENABLED
// ============================================================================

// Use double-precision floating point type
#define MFEM_USE_DOUBLE

// Throw an exception on errors (better for WASM error handling)
#define MFEM_USE_EXCEPTIONS

// Internal MFEM option: enable group/batch allocation for some small objects.
#define MFEM_USE_MEMALLOC

// Which library functions to use in class StopWatch for measuring time.
// 0 = std::clock (portable, low resolution)
// 2 = POSIX clock_gettime
// 3 = QueryPerformanceCounter (Windows)
// 4 = gettimeofday (POSIX)
// 6 = mach_absolute_time (Mac)
// For WASM, use std::clock which is portable
#define MFEM_TIMER_TYPE 0

// ============================================================================
// DISABLED Features - Not supported in WebAssembly
// ============================================================================

// MPI - Not available in browsers
// #define MFEM_USE_MPI

// GPU backends - Not available in WASM
// #define MFEM_USE_CUDA
// #define MFEM_USE_HIP

// OpenMP - Limited support in WASM
// #define MFEM_USE_OPENMP
// #define MFEM_USE_LEGACY_OPENMP

// Thread safety overhead not needed for single-threaded WASM
// #define MFEM_THREAD_SAFE

// SIMD - Emscripten has limited SIMD support
// #define MFEM_USE_SIMD

// Debug mode - disable for release builds
// #define MFEM_DEBUG

// Single precision - we use double
// #define MFEM_USE_SINGLE

// ============================================================================
// DISABLED External Libraries
// ============================================================================

// Linear algebra libraries
// #define MFEM_USE_LAPACK
// #define MFEM_USE_SUITESPARSE
// #define MFEM_USE_SUPERLU
// #define MFEM_USE_MUMPS
// #define MFEM_USE_STRUMPACK
// #define MFEM_USE_GINKGO
// #define MFEM_USE_AMGX
// #define MFEM_USE_MAGMA
// #define MFEM_USE_MKL_CPARDISO
// #define MFEM_USE_MKL_PARDISO

// Mesh partitioning
// #define MFEM_USE_METIS
// #define MFEM_USE_METIS_5

// ODE/DAE solvers
// #define MFEM_USE_SUNDIALS
// #define MFEM_USE_PETSC
// #define MFEM_USE_SLEPC

// I/O libraries
// #define MFEM_USE_ZLIB
// #define MFEM_USE_HDF5
// #define MFEM_USE_NETCDF
// #define MFEM_USE_ADIOS2
// #define MFEM_USE_CONDUIT
// #define MFEM_USE_FMS
// #define MFEM_USE_SIDRE

// Mesh infrastructure
// #define MFEM_USE_PUMI
// #define MFEM_USE_MOONOLITH
// #define MFEM_USE_GSLIB
// #define MFEM_USE_SIMMETRIX

// Network/Security
// #define MFEM_USE_GNUTLS

// Performance libraries
// #define MFEM_USE_RAJA
// #define MFEM_USE_OCCA
// #define MFEM_USE_CEED
// #define MFEM_USE_UMPIRE
// #define MFEM_USE_CALIPER

// Numerical libraries
// #define MFEM_USE_MPFR
// #define MFEM_USE_HIOP
// #define MFEM_USE_ALGOIM

// Automatic differentiation
// #define MFEM_USE_ADFORWARD
// #define MFEM_USE_CODIPACK
// #define MFEM_USE_ENZYME

// Backtraces
// #define MFEM_USE_LIBUNWIND

// Benchmarking
// #define MFEM_USE_BENCHMARK

#endif // MFEM_CONFIG_HEADER
