# wasm-sci

Scientific computing for TypeScript. NumPy-style arrays, SciPy-style algorithms, and SymPy-style symbolic math — all with WebAssembly acceleration.

**[Documentation & Demos](https://wasm-sci.dev/)** | **[npm: numwasm](https://www.npmjs.com/package/numwasm)** | **[npm: sciwasm](https://www.npmjs.com/package/sciwasm)** | **[npm: symwasm](https://www.npmjs.com/package/symwasm)**

## Packages

| Package | Description | Status |
|---------|-------------|--------|
| **[numwasm](packages/numwasm/)** | NumPy-compatible n-dimensional arrays. Linear algebra, FFT, random, broadcasting, 600+ functions. | Stable |
| **[sciwasm](packages/sciwasm/)** | SciPy-compatible scientific computing. Optimization, integration, interpolation, signal processing, statistics. | Scaffolded |
| **[symwasm](packages/symwasm/)** | SymPy-compatible symbolic math. Symbolic expressions, simplification, equation solving, calculus, printing. | Scaffolded |

## Installation

```bash
npm install numwasm        # N-dimensional arrays
npm install sciwasm        # Scientific computing (requires numwasm)
npm install symwasm        # Symbolic math
```

## Quick Start

```typescript
import * as nw from 'numwasm';
import * as sci from 'sciwasm';
import * as sym from 'symwasm';

// numwasm: N-dimensional arrays
const a = nw.array([[1, 2], [3, 4]]);
const b = nw.linalg.matmul(a, a);

// sciwasm: Scientific computing
const result = sci.optimize.minimize(f, x0);

// symwasm: Symbolic math
const x = new sym.core.Symbol('x');
const deriv = sym.calculus.diff(expr, x);
```

## numwasm

NumPy-inspired n-dimensional array operations in TypeScript, with performance-critical operations compiled to WebAssembly.

### Features

- **NumPy-compatible API** — familiar function names and semantics
- **WebAssembly acceleration** — BLAS, LAPACK, FFT, sorting, statistics, and random number generation in WASM
- **Full TypeScript types** — complete type definitions with documented parameters
- **Dual module support** — ESM and CommonJS, works in Node.js and browsers
- **600+ functions** across linear algebra, FFT, statistics, random, polynomials, string operations, masked arrays, and more

### Modules

- **Linear Algebra** (`linalg`) — matmul, inv, det, solve, eig, svd, qr, cholesky, norm, ...
- **FFT** (`fft`) — fft, ifft, rfft, fft2, fftfreq, fftshift, ...
- **Random** — Generator, PCG64, MT19937, normal, uniform, exponential, 20+ distributions
- **Masked Arrays** (`ma`) — masked array operations
- **Polynomials** — Polynomial, Chebyshev, Legendre, Hermite, Laguerre with full arithmetic
- **Strings** (`strings`) — element-wise string operations
- **Record Arrays** (`rec`) — structured/tabular data
- **Testing** (`testing`) — assert_allclose, assert_array_equal, ...

## sciwasm

SciPy-compatible scientific computing. Currently scaffolded with module stubs.

### Modules

- **Optimization** — minimize, least_squares, root_scalar, linprog, curve_fit
- **Integration** — quad, dblquad, tplquad, trapezoid, simpson, odeint
- **Interpolation** — interp1d, CubicSpline, PchipInterpolator, griddata
- **Statistics** — describe, distributions (norm, t, f, chi2), tests (pearsonr, ttest, kstest)
- **Signal Processing** — convolve, fftconvolve, butter, sosfilt, welch, spectrogram
- **Spatial** — KDTree, Delaunay, ConvexHull, Voronoi, distance
- **Special Functions** — gamma, beta, erf, bessel, factorial, comb, perm
- **Sparse Matrices** — csr_matrix, csc_matrix, eye, diags
- **Clustering** — kmeans, hierarchical clustering
- **N-D Image** — convolve, gaussian_filter, label, morphology
- **I/O** — loadmat, savemat
- **Constants** — physical constants (c, h, G, k, ...)

## symwasm

SymPy-compatible symbolic math. Currently scaffolded with module stubs.

### Modules

- **Core** — Symbol, Expr, Integer, Rational, Float, Add, Mul, Pow, pi, E, I
- **Simplify** — simplify, expand, factor, collect, cancel, trigsimp
- **Solvers** — solve, solveset, linsolve, nonlinsolve, dsolve
- **Calculus** — diff, integrate, limit, series, summation
- **Matrices** — Matrix, eye, zeros, ones, diag, det, inv, eigenvals
- **Printing** — latex, mathml, pretty, sstr

## AI / LLM Access

### llms.txt

- **[llms.txt](https://wasm-sci.dev/llms.txt)** — project overview
- **[llms-full.txt](https://wasm-sci.dev/llms-full.txt)** — complete API reference

### MCP Servers

Three MCP servers provide AI coding assistants with searchable access to the documentation:

| Server | Tools |
|--------|-------|
| `numwasm-mcp` | `search_numwasm_docs`, `list_numwasm_modules` |
| `sciwasm-mcp` | `search_sciwasm_docs`, `list_sciwasm_modules` |
| `symwasm-mcp` | `search_symwasm_docs`, `list_symwasm_modules` |

**Claude Desktop** — add to `claude_desktop_config.json`:

```json
{
  "mcpServers": {
    "numwasm": { "command": "npx", "args": ["-y", "numwasm-mcp"] },
    "sciwasm": { "command": "npx", "args": ["-y", "sciwasm-mcp"] },
    "symwasm": { "command": "npx", "args": ["-y", "symwasm-mcp"] }
  }
}
```

**Claude Code** — add `.mcp.json` to your project root:

```json
{
  "mcpServers": {
    "numwasm": { "command": "npx", "args": ["-y", "numwasm-mcp"] },
    "sciwasm": { "command": "npx", "args": ["-y", "sciwasm-mcp"] },
    "symwasm": { "command": "npx", "args": ["-y", "symwasm-mcp"] }
  }
}
```

## Development

### Prerequisites

- Node.js >= 18
- [pnpm](https://pnpm.io/)
- [Emscripten](https://emscripten.org/) (for building numwasm WASM from C sources)

### Build

```bash
# Install dependencies
pnpm install

# Build all packages
pnpm run build

# Build individual packages
pnpm run build:numwasm   # WASM + TypeScript
pnpm run build:sciwasm   # TypeScript only
pnpm run build:symwasm   # TypeScript only
```

### Test

```bash
# Run all tests
pnpm run test

# Test individual packages
pnpm --filter numwasm run test
pnpm --filter sciwasm run test
pnpm --filter symwasm run test
```

### Documentation Site

```bash
# Generate TypeDoc JSON for all packages
pnpm run docs

# Dev server
pnpm run dev:docs

# Full static build (SSG + sitemap + llms.txt)
cd docs-site && pnpm run build:ssg
```

## Project Structure

```
packages/
  numwasm/           NumPy-compatible arrays (TypeScript + WebAssembly)
  sciwasm/           SciPy-compatible scientific computing
  symwasm/           SymPy-compatible symbolic math
mcp-servers/
  numwasm-mcp/       MCP server for numwasm docs
  sciwasm-mcp/       MCP server for sciwasm docs
  symwasm-mcp/       MCP server for symwasm docs
docs-site/           React + Vite documentation website
docs/                Markdown documentation
```

## License

MIT
