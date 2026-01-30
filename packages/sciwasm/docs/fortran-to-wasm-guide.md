# Compiling Fortran to WebAssembly: A Complete Guide

This guide explains how to compile Fortran code (including BLAS and LAPACK) to WebAssembly using LLVM Flang. Based on the excellent work documented at https://gws.phd/posts/fortran_wasm/

## Overview of Available Tools

| Tool | Status | Best For |
|------|--------|----------|
| **f2c** | Mature but limited | Fortran 77 only, used by Pyodide |
| **LFortran** | Alpha | Future option, not production-ready |
| **Dragonegg** | Deprecated | Historical interest only (requires GCC 4.8) |
| **Classic Flang** | End of life | No 32-bit (wasm32) support |
| **LLVM Flang** | Active development | **Best current option** (requires patches) |

### Why LLVM Flang?

- Modern Fortran support (F77, F90, F95, F2003, F2008, F2018)
- Active development as part of official LLVM project
- Can target WebAssembly with small patches
- Produces high-quality optimized code

## Prerequisites

1. **Development Environment**
   - Linux or macOS (or WSL on Windows)
   - ~50GB disk space for LLVM build
   - 8+ GB RAM recommended
   - CMake 3.20+
   - Ninja build system
   - Python 3.x

2. **Emscripten SDK**
   ```bash
   git clone https://github.com/emscripten-core/emsdk.git
   cd emsdk
   ./emsdk install latest
   ./emsdk activate latest
   source ./emsdk_env.sh
   ```

3. **Native Fortran Compiler** (for testing)
   ```bash
   # Ubuntu/Debian
   sudo apt-get install gfortran

   # macOS
   brew install gcc
   ```

## Step 1: Clone LLVM Source

```bash
git clone --depth=1 --branch=llvmorg-18.1.1 https://github.com/llvm/llvm-project.git
```

> **Note:** Use a specific tag (like 18.1.1) rather than main branch for stability.

## Step 2: Apply WebAssembly Target Patch

The key issue is that Flang doesn't know how to handle wasm32 target ABI (Application Binary Interface). We need to add support for it.

### Patch 1: Add wasm32 Target Support

Edit `llvm-project/flang/lib/Optimizer/CodeGen/Target.cpp`:

Add the following `TargetWasm32` structure (similar to existing targets):

```cpp
namespace {
// ... existing code ...

/// WebAssembly 32-bit target specifics
struct TargetWasm32 : public GenericTarget<TargetWasm32> {
  using GenericTarget::GenericTarget;

  // Complex numbers are passed as {real, imag} tuples
  CodeGenSpecifics::Marshalling
  complexArgumentType(mlir::Location loc, mlir::Type eleTy) const override {
    CodeGenSpecifics::Marshalling marshal;
    // wasm32 uses 4-byte alignment
    auto structTy = mlir::TupleType::get(eleTy.getContext(), {eleTy, eleTy});
    marshal.emplace_back(fir::ReferenceType::get(structTy),
                         AT{/*align=*/4, /*byval=*/true});
    return marshal;
  }

  // Complex return values use struct return (sret)
  CodeGenSpecifics::Marshalling
  complexReturnType(mlir::Location loc, mlir::Type eleTy) const override {
    CodeGenSpecifics::Marshalling marshal;
    auto structTy = mlir::TupleType::get(eleTy.getContext(), {eleTy, eleTy});
    marshal.emplace_back(fir::ReferenceType::get(structTy),
                         AT{/*align=*/4, /*sret=*/true});
    return marshal;
  }
};

} // namespace
```

Then register it in `CodeGenSpecifics::get()`:

```cpp
// In the switch statement for target triples, add:
case llvm::Triple::wasm32:
  return std::make_unique<TargetWasm32>(
      ctx, std::move(trp), std::move(kindMap), targetCPU, targetFeatures, dl);
```

### Patch 2: Fix Size Mismatch for Cross-Compilation

The problem: When cross-compiling from 64-bit host to 32-bit wasm32, `sizeof(long)` differs:
- Host (x86_64): `sizeof(long) = 8`
- Target (wasm32): `sizeof(long) = 4`

Edit `llvm-project/flang/include/flang/Optimizer/Builder/Runtime/RTBuilder.h`:

```cpp
// Find the getModel<long>() specialization and change it:
template <>
constexpr TypeBuilderFunc getModel<long>() {
  return [](mlir::MLIRContext *context) -> mlir::Type {
    // Force 4-byte (i32) for wasm32 target instead of using sizeof()
    return mlir::IntegerType::get(context, 8 * 4);  // 32 bits
  };
}

template <>
constexpr TypeBuilderFunc getModel<unsigned long>() {
  return [](mlir::MLIRContext *context) -> mlir::Type {
    return mlir::IntegerType::get(context, 8 * 4);  // 32 bits
  };
}
```

Also edit `llvm-project/flang/lib/Optimizer/CodeGen/CodeGen.cpp` to fix malloc size arguments.

## Step 3: Build Patched LLVM Flang

```bash
cmake -G Ninja -S llvm-project/llvm -B build \
  -DCMAKE_INSTALL_PREFIX=$PWD/llvm-install \
  -DCMAKE_BUILD_TYPE=MinSizeRel \
  -DLLVM_DEFAULT_TARGET_TRIPLE="wasm32-unknown-emscripten" \
  -DLLVM_TARGETS_TO_BUILD="WebAssembly" \
  -DLLVM_ENABLE_PROJECTS="clang;flang;mlir"

cmake --build build -j$(nproc)
```

> **Build time:** Expect 1-3 hours depending on hardware.

After building, verify:
```bash
./build/bin/flang-new --version
```

## Step 4: Build Fortran Runtime Library

Flang needs a runtime library for I/O, memory allocation, etc. We need to compile it with Emscripten.

Create a Makefile in `llvm-project/flang/runtime/`:

```makefile
CXX = em++
AR = emar
RANLIB = emranlib

CXXFLAGS = -O2 -I../include -I../../build/include \
           -DFLANG_LITTLE_ENDIAN=1

SOURCES = $(wildcard *.cpp)
OBJECTS = $(SOURCES:.cpp=.o)

libFortranRuntime.a: $(OBJECTS)
	$(AR) rcs $@ $^
	$(RANLIB) $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f *.o libFortranRuntime.a
```

Build it:
```bash
cd llvm-project/flang/runtime
make
```

## Step 5: Build BLAS for WebAssembly

Download BLAS:
```bash
wget https://www.netlib.org/blas/blas-3.12.0.tgz
tar xzf blas-3.12.0.tgz
cd BLAS-3.12.0
```

Edit `make.inc`:
```makefile
FC = /path/to/build/bin/flang-new
FFLAGS = -O2
AR = emar
ARFLAGS = cr
RANLIB = emranlib
```

Build:
```bash
make
# Produces: blas_LINUX.a
```

## Step 6: Build LAPACK for WebAssembly

Download LAPACK:
```bash
wget https://github.com/Reference-LAPACK/lapack/archive/refs/tags/v3.12.0.tar.gz
tar xzf v3.12.0.tar.gz
cd lapack-3.12.0
cp make.inc.example make.inc
```

Edit `make.inc`:
```makefile
FC = /full/path/to/build/bin/flang-new
FFLAGS = -O2 -frecursive
AR = emar
ARFLAGS = cr
RANLIB = emranlib

# Use internal timer (not system timer)
TIMER = INT_CPU_TIME

# Point to your BLAS
BLASLIB = /path/to/BLAS-3.12.0/blas_LINUX.a
```

Build:
```bash
make lib
# Produces: liblapack.a
```

## Step 7: Link Everything Together

Create a simple test program `test.f90`:

```fortran
program test
  implicit none
  double precision :: A(3,3), B(3), X(3)
  integer :: ipiv(3), info

  ! Simple 3x3 system: A*X = B
  A = reshape([1,0,0, 0,2,0, 0,0,3], [3,3])
  B = [1, 4, 9]

  call DGESV(3, 1, A, 3, ipiv, B, 3, info)

  print *, 'Solution:', B
  print *, 'Info:', info
end program
```

Compile and link:
```bash
# Compile Fortran to object file
./build/bin/flang-new -c test.f90 -o test.o

# Link with Emscripten
emcc test.o \
  -L./lapack-3.12.0 -llapack \
  -L./BLAS-3.12.0 -lblas \
  -L./build/flang/runtime -lFortranRuntime \
  -s EXPORTED_FUNCTIONS='["_main","_malloc","_free"]' \
  -s EXPORTED_RUNTIME_METHODS='["ccall","cwrap"]' \
  -o test.js
```

Run:
```bash
node test.js
```

## Alternative: Using Docker

The author provides a pre-built Docker container with all patches applied:

```bash
docker pull ghcr.io/pgerber/flang-wasm:llvm-18.1.1

docker run -it --rm -v $(pwd):/work ghcr.io/pgerber/flang-wasm:llvm-18.1.1 \
  flang-new -c /work/mycode.f90 -o /work/mycode.o
```

## Calling Fortran from JavaScript

### Memory Management

Fortran functions expect pointers. In JavaScript, you must:
1. Allocate memory with `Module._malloc()`
2. Write values to the WebAssembly heap
3. Call the function
4. Read results from the heap
5. Free memory with `Module._free()`

### Example: Calling a Simple Subroutine

Fortran code (`add.f90`):
```fortran
subroutine add(x, y, z)
  integer, intent(in) :: x, y
  integer, intent(out) :: z
  z = x + y
end subroutine
```

JavaScript code:
```javascript
// Allocate memory for three integers (4 bytes each)
const xPtr = Module._malloc(4);
const yPtr = Module._malloc(4);
const zPtr = Module._malloc(4);

// Write input values
Module.HEAP32[xPtr / 4] = 123;
Module.HEAP32[yPtr / 4] = 456;

// Call Fortran subroutine (note: name has trailing underscore)
Module._add_(xPtr, yPtr, zPtr);

// Read result
const result = Module.HEAP32[zPtr / 4];
console.log(result);  // 579

// Free memory
Module._free(xPtr);
Module._free(yPtr);
Module._free(zPtr);
```

### Working with Arrays

```javascript
// Allocate a 10-element double array
const n = 10;
const arrayPtr = Module._malloc(n * 8);  // 8 bytes per double

// Write values
for (let i = 0; i < n; i++) {
  Module.HEAPF64[(arrayPtr / 8) + i] = i * 1.5;
}

// Call Fortran function...

// Read values back
for (let i = 0; i < n; i++) {
  console.log(Module.HEAPF64[(arrayPtr / 8) + i]);
}

Module._free(arrayPtr);
```

## Key Technical Details

### Name Mangling

Fortran compilers add a trailing underscore to symbol names:
- `subroutine foo` → `_foo_`
- `function bar` → `_bar_`

### Fortran Array Layout

Fortran uses **column-major** order (opposite of C/JavaScript):
- A(i,j) in Fortran = A[j*rows + i] in C
- When passing 2D arrays, be aware of this difference

### Complex Numbers

Complex numbers in Fortran are passed as structs:
```c
struct complex_t {
  double real;
  double imag;
};
```

## Limitations

1. **Patched LLVM Required**: The patches are not upstream, so you must maintain a fork
2. **Cross-compilation quirks**: Some platform-specific assumptions may cause issues
3. **No I/O to files**: WebAssembly has no direct file system access
4. **Build complexity**: The full toolchain takes significant time and space to build

## Summary

| Step | Description |
|------|-------------|
| 1 | Clone LLVM 18.1.1 source |
| 2 | Apply wasm32 target patch |
| 3 | Apply size mismatch patch |
| 4 | Build LLVM with Flang |
| 5 | Build Fortran runtime with emcc |
| 6 | Build BLAS with flang-new |
| 7 | Build LAPACK with flang-new |
| 8 | Link everything with emcc |

## References

- Original article: https://gws.phd/posts/fortran_wasm/
- LLVM Flang: https://flang.llvm.org/
- Emscripten: https://emscripten.org/
- BLAS/LAPACK: https://netlib.org/lapack/
