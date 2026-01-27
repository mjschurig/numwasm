# numwasm

NumPy-inspired n-dimensional array operations in TypeScript, with performance-critical operations compiled to WebAssembly.

**[Documentation & Demos](https://numwasm.quebi.de/)** | **[npm](https://www.npmjs.com/package/numjs-wasm)**

## Features

- **NumPy-compatible API** — familiar function names and semantics (`zeros`, `ones`, `linspace`, `matmul`, `fft`, `mean`, `reshape`, ...)
- **WebAssembly acceleration** — BLAS, LAPACK, FFT, sorting, statistics, and random number generation implemented in C and compiled to WASM
- **Full TypeScript types** — complete type definitions with documented parameters
- **Dual module support** — ESM and CommonJS, works in Node.js and browsers
- **600+ functions** across linear algebra, FFT, statistics, random, polynomials, string operations, masked arrays, and more

## Installation

```bash
npm install numjs-wasm
```

## Quick Start

```typescript
import {
  loadWasmModule,
  array,
  zeros,
  ones,
  linspace,
  reshape,
  add,
  matmul,
  mean,
} from "numjs-wasm";

// Initialize the WASM module (required once before use)
await loadWasmModule();

// Create arrays
const a = array([1, 2, 3, 4, 5, 6]);
const b = zeros([3, 3]);
const c = ones([2, 3]);
const d = linspace(0, 1, 100);

// Reshape and compute
const matrix = reshape(a, [2, 3]);
const result = add(matrix, c);
const product = await matmul(reshape(a, [2, 3]), reshape(a, [3, 2]));

// Statistics
const avg = mean(d);
```

## Modules

### Linear Algebra (`linalg`)

```typescript
import { linalg } from "numjs-wasm";

const result = await linalg.matmul(a, b);
const { values, vectors } = await linalg.eig(matrix);
const solution = await linalg.solve(coefficients, constants);
const determinant = await linalg.det(matrix);
```

`matmul`, `dot`, `inv`, `det`, `solve`, `eig`, `eigh`, `svd`, `qr`, `cholesky`, `norm`, `cond`, `lstsq`, `matrix_rank`, `matrix_power`, `cross`, `kron`, `tensordot`, ...

### FFT

```typescript
import { fftModule } from "numjs-wasm";

const spectrum = fftModule.fft(signal);
const freqs = fftModule.fftfreq(n, dt);
```

`fft`, `ifft`, `rfft`, `irfft`, `fft2`, `ifft2`, `fftn`, `ifftn`, `fftfreq`, `rfftfreq`, `fftshift`, `ifftshift`, `hfft`, `ihfft`

### Random

```typescript
import { default_rng, Generator } from "numjs-wasm";

const rng = default_rng(42);
const samples = rng.normal(0, 1, [1000]);
const uniform = rng.uniform(0, 1, [100]);
```

Bit generators: `PCG64`, `MT19937`, `Philox`, `SFC64`. Distributions: normal, uniform, exponential, gamma, beta, binomial, poisson, and 20+ more.

### Masked Arrays (`ma`)

```typescript
import { ma } from "numjs-wasm";

const masked = ma.array(data, { mask: [false, false, true, false] });
const avg = ma.average(masked);
```

### Polynomials

```typescript
import { Polynomial, Chebyshev } from "numjs-wasm";

const p = new Polynomial([1, 2, 3]); // 1 + 2x + 3x^2
const roots = p.roots();
```

`Polynomial`, `Chebyshev`, `Legendre`, `Hermite`, `HermiteE`, `Laguerre` with full arithmetic, fitting, roots, and conversions.

### Other Modules

- **Strings** (`strings`) — element-wise string operations on arrays
- **Record Arrays** (`rec`) — structured/tabular data with named fields
- **Testing** (`testing`) — `assert_allclose`, `assert_array_equal`, `assert_raises`, ...

## Core API

| Category           | Functions                                                                                                                            |
| ------------------ | ------------------------------------------------------------------------------------------------------------------------------------ |
| **Array Creation** | `array`, `zeros`, `ones`, `empty`, `full`, `arange`, `linspace`, `logspace`, `geomspace`, `eye`, `identity`, `diag`, `meshgrid`, ... |
| **Manipulation**   | `reshape`, `transpose`, `concatenate`, `stack`, `split`, `flip`, `roll`, `rot90`, `tile`, `repeat`, `pad`, ...                       |
| **Math**           | `add`, `subtract`, `multiply`, `divide`, `power`, `sqrt`, `exp`, `log`, `sin`, `cos`, `tan`, `abs`, `clip`, ...                      |
| **Statistics**     | `mean`, `median`, `std`, `var_`, `min`, `max`, `sum`, `prod`, `histogram`, `percentile`, `quantile`, ...                             |
| **Sorting**        | `sort`, `argsort`, `argmax`, `argmin`, `searchsorted`, `partition`, ...                                                              |
| **Logic**          | `all`, `any`, `where`, `logical_and`, `logical_or`, `logical_not`, ...                                                               |
| **Comparison**     | `equal`, `greater`, `less`, `allclose`, `isclose`, `isnan`, `isinf`, ...                                                             |
| **Set Operations** | `unique`, `union1d`, `intersect1d`, `setdiff1d`, `isin`, ...                                                                         |
| **I/O**            | `save`, `load`, `loadtxt`, `savetxt`, `genfromtxt`, `frombuffer`, ...                                                                |
| **Constants**      | `pi`, `e`, `inf`, `nan`, `euler_gamma`, `newaxis`                                                                                    |

## AI / LLM Access

The documentation site serves machine-readable files following the [llms.txt](https://llmstxt.org/) convention:

- **[llms.txt](https://numwasm.quebi.de/llms.txt)** — project overview with module links
- **[llms-full.txt](https://numwasm.quebi.de/llms-full.txt)** — complete API reference (all 600+ functions with signatures and descriptions)

### MCP Server

The `numwasm-mcp` package provides a [Model Context Protocol](https://modelcontextprotocol.io/) server that gives AI coding assistants searchable access to the full API docs. It ships with a bundled docs index — no network calls at runtime.

**Tools exposed:**

- `search_numwasm_docs` — search by function name, module, category, or keyword
- `list_numwasm_modules` — list all modules and categories

**Claude Desktop** — add to `claude_desktop_config.json`:

```json
{
  "mcpServers": {
    "numwasm-docs": {
      "command": "npx",
      "args": ["-y", "numwasm-mcp"]
    }
  }
}
```

**Claude Code** — add `.mcp.json` to your project root:

```json
{
  "mcpServers": {
    "numwasm-docs": {
      "command": "npx",
      "args": ["-y", "numwasm-mcp"]
    }
  }
}
```

## Development

### Prerequisites

- Node.js >= 18
- [pnpm](https://pnpm.io/)
- [Emscripten](https://emscripten.org/) (for building WASM from C sources)
- Python 3 + NumPy (for generating test fixtures and running comparison benchmarks)

### Build

```bash
# Install dependencies
pnpm install

# Build WASM + TypeScript library
npm run build

# Or step by step:
npm run build:wasm    # Compile C → WebAssembly
npm run build:lib     # Bundle TypeScript with Vite
```

### Test

```bash
# Run all tests
npm test

# Watch mode
npm run test:watch

# Run comparison tests against NumPy reference vectors
npm run test:compare

# Browser tests (Playwright)
npm run test:browser
```

### Benchmark

```bash
# Full benchmark pipeline (build + NumPy + NumJS + report)
npm run benchmark

# Individual steps
npm run benchmark:numpy
npm run benchmark:numjs
npm run benchmark:combine
```

### Documentation Site

```bash
# Generate TypeDoc JSON
npm run docs

# Dev server
npm run dev:docs

# Full static build (SSG + sitemap + llms.txt)
cd docs-site && pnpm run build:ssg
```

## Project Structure

```
src/
  ts/              TypeScript implementation (NDArray, ufuncs, modules)
  wasm/            C source files compiled to WebAssembly
scripts/
  build-wasm.sh    Emscripten build script
tests/
  ts/              Vitest unit & integration tests
  browser/         Playwright browser tests
  python/          NumPy test fixture generators
benchmark/
  ts/              TypeScript benchmark suites
  python/          NumPy comparison benchmarks
docs-site/         React + Vite documentation website
dist/              Build output (library + WASM binary)
```

## License

MIT
