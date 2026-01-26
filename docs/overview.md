# NumJS-WASM Overview

NumJS-WASM is a TypeScript/WebAssembly implementation of NumPy's n-dimensional array library. It provides the same features and, as far as possible, the same syntax and naming conventions as NumPy, while running entirely in JavaScript/WASM environments.

## Goals

- **NumPy Compatibility**: Mirror NumPy's API, naming conventions, and behavior
- **High Performance**: Leverage WebAssembly for near-native computation speed
- **Numerical Accuracy**: Implement NumPy's algorithms (like pairwise summation) for precise results
- **Type Safety**: Full TypeScript support with comprehensive type definitions

## Architecture

```
┌─────────────────────────────────────────────────────────┐
│                    TypeScript API                        │
│              (NDArray, factory methods)                  │
├─────────────────────────────────────────────────────────┤
│                   WASM Loader                            │
│           (Module initialization, memory)                │
├─────────────────────────────────────────────────────────┤
│                  WebAssembly Module                      │
│         (C implementation via Emscripten)                │
└─────────────────────────────────────────────────────────┘
```

### Components

| Directory | Purpose |
|-----------|---------|
| `src/ts/` | TypeScript source - NDArray class, types, WASM loader |
| `src/wasm/` | C source - array implementation, pairwise summation |
| `tests/` | Vitest test suite with NumPy-generated test vectors |
| `benchmark/` | Performance benchmarks comparing to NumPy |
| `docs/` | Documentation and NumPy reference materials |

## Current Features

### NDArray Class

**Factory Methods:**
- `NDArray.zeros(shape, options)` - Create zero-filled arrays
- `NDArray.ones(shape, options)` - Create one-filled arrays
- `NDArray.fromArray(data, shape, options)` - Create from existing data
- `NDArray.arange(start, end, step)` - Create evenly-spaced sequences

**Properties:**
- `shape` - Array dimensions
- `ndim` - Number of dimensions
- `size` - Total element count
- `dtype` - Data type (Float32, Float64, Int32, Int64)

**Operations:**
- `sum()` - Sum all elements using pairwise summation
- `fill(value)` - Fill with a constant value
- `toArray()` - Export to JavaScript array
- `dispose()` - Free WASM memory

### Supported Data Types

```typescript
enum DType {
  Float32 = 0,
  Float64 = 1,
  Int32 = 2,
  Int64 = 3,
}
```

## Usage Example

```typescript
import { NDArray, DType } from 'numjs-wasm';

// Create arrays
const zeros = NDArray.zeros([3, 4]);           // 3x4 zero matrix
const ones = NDArray.ones([2, 2, 2]);          // 2x2x2 cube of ones
const range = NDArray.arange(0, 10, 1);        // [0, 1, 2, ..., 9]

// Create from data
const data = NDArray.fromArray(
  [1, 2, 3, 4, 5, 6],
  [2, 3],
  { dtype: DType.Float64 }
);

// Operations
console.log(data.sum());    // 21
console.log(data.shape);    // [2, 3]
console.log(data.ndim);     // 2

// Memory management
data.dispose();
```

## Numerical Accuracy

NumJS implements NumPy's **pairwise summation** algorithm for floating-point reductions:

- **O(log n) rounding error** instead of O(n) for naive summation
- 8-way loop unrolling for CPU pipelining
- Tree-structured accumulation to minimize precision loss
- Validated against NumPy with generated test vectors

## Build & Development

```bash
# Install dependencies
npm install

# Build WASM and TypeScript
npm run build

# Run tests
npm test

# Run benchmarks
npm run benchmark
```

### Requirements

- Node.js 18+
- Emscripten SDK (for WASM compilation)
- Python 3 with NumPy (for test vector generation)

## Testing

The test suite validates NumJS against NumPy:

1. **Unit tests** - Core functionality, edge cases, memory safety
2. **Comparison tests** - Results validated against NumPy-generated test vectors
3. **Benchmarks** - Performance comparison across array sizes (10² to 10⁷ elements)

## Roadmap

The architecture supports incremental addition of NumPy features:

- Element-wise operations (add, subtract, multiply, divide)
- Array indexing and slicing
- Reshaping and transposition
- Broadcasting
- Linear algebra operations
- Browser support (currently Node.js only)

## License

See LICENSE file for details.
