#!/bin/bash
set -e

echo "Building GMP for WebAssembly..."

# Configuration
GMP_VERSION="6.3.0"
GMP_DIR="gmp-${GMP_VERSION}"
GMP_TARBALL="${GMP_DIR}.tar.xz"
BUILD_DIR=".gmp-build"

# Mirror URLs (try in order)
GMP_MIRRORS=(
  "https://ftp.gnu.org/gnu/gmp/${GMP_TARBALL}"
  "https://mirrors.kernel.org/gnu/gmp/${GMP_TARBALL}"
  "https://mirror.us-midwest-1.nexcess.net/gnu/gmp/${GMP_TARBALL}"
  "https://ftpmirror.gnu.org/gmp/${GMP_TARBALL}"
  "https://gmplib.org/download/gmp/${GMP_TARBALL}"
)

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
  downloaded=0
  for mirror in "${GMP_MIRRORS[@]}"; do
    echo "  Trying: $mirror"
    if curl -L --connect-timeout 10 --max-time 120 -o "$GMP_TARBALL" "$mirror" 2>/dev/null; then
      # Verify file was actually downloaded (not empty or error page)
      if [ -s "$GMP_TARBALL" ] && file "$GMP_TARBALL" | grep -q "XZ compressed"; then
        echo "  Success!"
        downloaded=1
        break
      else
        echo "  Downloaded file is invalid, trying next mirror..."
        rm -f "$GMP_TARBALL"
      fi
    else
      echo "  Failed, trying next mirror..."
      rm -f "$GMP_TARBALL"
    fi
  done

  if [ $downloaded -eq 0 ]; then
    echo "Error: Failed to download GMP from all mirrors"
    exit 1
  fi
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
