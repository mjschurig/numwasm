#!/bin/bash
set -e

echo "Building GMP for WebAssembly..."

# Configuration
GMP_VERSION="6.3.0"
GMP_DIR="gmp-${GMP_VERSION}"
GMP_TARBALL="${GMP_DIR}.tar.xz"
GMP_URL="https://gmplib.org/download/gmp/${GMP_TARBALL}"
BUILD_DIR=".gmp-build"

# Check if already built
if [ -f "$BUILD_DIR/lib/libgmp.a" ]; then
  echo "GMP already built at $BUILD_DIR"
  echo "To rebuild, run: rm -rf $BUILD_DIR $GMP_DIR"
  exit 0
fi

# Create build directory
mkdir -p "$BUILD_DIR"
BUILD_DIR_ABS="$(cd "$BUILD_DIR" && pwd)"

# Download GMP if not present
if [ ! -f "$GMP_TARBALL" ]; then
  echo "Downloading GMP $GMP_VERSION..."
  curl -L -o "$GMP_TARBALL" "$GMP_URL"
fi

# Extract if not already extracted
if [ ! -d "$GMP_DIR" ]; then
  echo "Extracting GMP..."
  tar xf "$GMP_TARBALL"
fi

# Build GMP with Emscripten
cd "$GMP_DIR"

echo "Configuring GMP for Emscripten..."
# Use emconfigure to set up the build environment
emconfigure ./configure \
  --prefix="$BUILD_DIR_ABS" \
  --host=none \
  --disable-shared \
  --enable-static \
  --disable-assembly \
  --disable-fft \
  CC_FOR_BUILD=gcc \
  CFLAGS="-O2"

echo "Building GMP..."
emmake make -j$(nproc 2>/dev/null || echo 4)

echo "Installing GMP..."
emmake make install

cd ..

echo "GMP build complete!"
echo "Headers: $BUILD_DIR/include/gmp.h"
echo "Library: $BUILD_DIR/lib/libgmp.a"
