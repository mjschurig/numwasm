#!/bin/bash
set -e

SRC_DIR="reference/symengine/symengine"
DEST_DIR="src/wasm/symengine"

# Check if files already exist
if [ -f "$DEST_DIR/cwrapper.h" ]; then
  echo "SymEngine source files already exist in $DEST_DIR"
  echo "Skipping copy (use 'pnpm run clean:all' to force re-copy)"
  exit 0
fi

echo "Copying SymEngine source files from $SRC_DIR to $DEST_DIR..."

# Create destination directory
mkdir -p "$DEST_DIR"

# Copy all .h and .cpp files from main directory
echo "Copying header and source files..."
cp "$SRC_DIR"/*.h "$DEST_DIR/" 2>/dev/null && echo "  ✓ All .h files" || echo "  ⚠ Some .h files not found"
cp "$SRC_DIR"/*.cpp "$DEST_DIR/" 2>/dev/null && echo "  ✓ All .cpp files" || echo "  ⚠ Some .cpp files not found"
cp "$SRC_DIR"/*.inc "$DEST_DIR/" 2>/dev/null && echo "  ✓ All .inc files" || echo "  ⚠ Some .inc files not found"

# Copy essential subdirectories
echo "Copying subdirectories..."

# Utilities (Teuchos RCP) - Required
if [ -d "$SRC_DIR/utilities" ]; then
  cp -r "$SRC_DIR/utilities" "$DEST_DIR/"
  echo "  ✓ utilities/ directory"
fi

# Printers (string conversion) - Required
if [ -d "$SRC_DIR/printers" ]; then
  cp -r "$SRC_DIR/printers" "$DEST_DIR/"
  echo "  ✓ printers/ directory"
fi

# Parser (string parsing) - Required for basic_parse
if [ -d "$SRC_DIR/parser" ]; then
  cp -r "$SRC_DIR/parser" "$DEST_DIR/"
  echo "  ✓ parser/ directory"
fi

# Polys (polynomials) - Required for visitor.h
if [ -d "$SRC_DIR/polys" ]; then
  cp -r "$SRC_DIR/polys" "$DEST_DIR/"
  echo "  ✓ polys/ directory"
fi

# Matrices - Required for some operations
if [ -d "$SRC_DIR/matrices" ]; then
  cp -r "$SRC_DIR/matrices" "$DEST_DIR/"
  echo "  ✓ matrices/ directory"
fi

# Copy SymEngine license
if [ -f "reference/symengine/LICENSE" ]; then
  cp "reference/symengine/LICENSE" "src/wasm/SYMENGINE_LICENSE"
  echo "  ✓ LICENSE file"
fi

echo "✓ SymEngine source files copied successfully"
echo ""
echo "Files copied to: $DEST_DIR"
echo "Total files: $(find "$DEST_DIR" -type f | wc -l)"
