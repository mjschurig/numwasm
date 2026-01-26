# Phase 25: Advanced Linear Algebra Implementation Plan

Complete implementation roadmap for tensor operations, Einstein summation, and additional matrix operations.

---

## ⚠️ Implementation Guideline

**The original NumPy source code is available in `/numpy`.**

Key reference files:
- `numpy/_core/einsumfunc.py` - Einstein summation (~1,650 lines)
- `numpy/_core/numeric.py` - tensordot, kron (~2,711 lines)
- `numpy/_core/fromnumeric.py` - cross product
- `numpy/linalg/_linalg.py` - multi_dot, tensorsolve, tensorinv

Implementation should follow NumPy's algorithms and error handling for consistency.

---

## Current State (Pre-Phase 25)

```
Already Implemented (numpy.linalg):
├── dot(a, b) - vector/matrix dot product
├── vdot(a, b) - vector dot (flattened)
├── inner(a, b) - inner product
├── outer(a, b) - outer product
├── matmul(a, b) - matrix multiplication
├── solve, inv, det, etc.
└── BLAS/LAPACK infrastructure

Missing:
├── tensordot(a, b, axes)
├── einsum(subscripts, *operands)
├── einsum_path(subscripts, *operands, optimize)
├── kron(a, b)
├── cross(a, b, axis)
├── multi_dot(arrays)
├── tensorsolve(a, b, axes)
├── tensorinv(a, ind)
├── matrix_norm(x, ord, keepdims)
└── vector_norm(x, ord, axis, keepdims)
```

---

## Phase 25 Dependency Tree

```
PHASE 25: ADVANCED LINEAR ALGEBRA
│
├── 25.1 Tensor Contraction (TypeScript)
│   ├── 25.1.1 tensordot(a, b, axes)
│   │   ├── Axis specification handling
│   │   ├── Shape computation
│   │   └── Reduction via matmul
│   │
│   └── 25.1.2 multi_dot(arrays)
│       ├── Optimal order computation (dynamic programming)
│       └── Chain multiplication execution
│
│   Dependencies: matmul, transpose, reshape
│
├── 25.2 Einstein Summation (TypeScript)
│   ├── 25.2.1 einsum(subscripts, *operands)
│   │   ├── Subscript parser
│   │   ├── Implicit output computation
│   │   ├── Contraction execution
│   │   └── Optimization paths
│   │
│   └── 25.2.2 einsum_path(subscripts, *operands, optimize)
│       ├── Greedy path finding
│       ├── Optimal path finding
│       └── Path cost estimation
│
│   Dependencies: 25.1.1 (tensordot)
│
├── 25.3 Special Products (TypeScript + WASM)
│   ├── 25.3.1 kron(a, b)
│   │   └── Kronecker product via block expansion
│   │
│   └── 25.3.2 cross(a, b, axis)
│       ├── 2D cross product (scalar)
│       ├── 3D cross product (vector)
│       └── Batched cross product
│
│   Dependencies: multiply, subtract, reshape
│
├── 25.4 Tensor Solving (TypeScript)
│   ├── 25.4.1 tensorsolve(a, b, axes)
│   │   └── Reshape to matrix equation
│   │
│   └── 25.4.2 tensorinv(a, ind)
│       └── Reshape to matrix inverse
│
│   Dependencies: solve, inv, prod
│
└── 25.5 Norm Functions (TypeScript)
    ├── 25.5.1 matrix_norm(x, ord, keepdims)
    └── 25.5.2 vector_norm(x, ord, axis, keepdims)

    Dependencies: svdvals, sum, max, abs
```

---

## Detailed Implementation Specifications

### 25.1 Tensor Contraction

#### 25.1.1 tensordot

**File:** `src/ts/linalg.ts` (additions)

```typescript
/**
 * Compute tensor dot product along specified axes.
 *
 * For tensors a and b, tensordot(a, b, axes) computes the sum of products
 * over the specified axes.
 *
 * @param a - First tensor
 * @param b - Second tensor
 * @param axes - Axes to sum over:
 *   - number N: last N axes of a with first N axes of b
 *   - [axesA, axesB]: specific axes for each tensor
 * @returns Tensor contraction result
 *
 * @example
 * // Matrix multiplication
 * tensordot([[1, 2], [3, 4]], [[5, 6], [7, 8]], axes=1)
 * // Same as: [[1, 2], [3, 4]] @ [[5, 6], [7, 8]]
 *
 * // Outer product
 * tensordot([1, 2, 3], [4, 5], axes=0)
 * // Shape: (3, 2)
 *
 * // Specified axes
 * tensordot(a, b, [[1, 2], [0, 1]])
 * // Contract a's axes 1,2 with b's axes 0,1
 */
export function tensordot(
  a: NDArray,
  b: NDArray,
  axes: number | [number[], number[]] = 2
): NDArray {
  // Parse axes specification
  let axesA: number[];
  let axesB: number[];

  if (typeof axes === 'number') {
    if (axes < 0) {
      throw new ValueError('axes must be non-negative');
    }

    // Last N axes of a, first N axes of b
    axesA = [];
    axesB = [];
    for (let i = 0; i < axes; i++) {
      axesA.push(a.ndim - axes + i);
      axesB.push(i);
    }
  } else {
    [axesA, axesB] = axes;

    if (axesA.length !== axesB.length) {
      throw new ValueError('axes lists must have same length');
    }
  }

  // Normalize negative axes
  axesA = axesA.map(ax => ax < 0 ? a.ndim + ax : ax);
  axesB = axesB.map(ax => ax < 0 ? b.ndim + ax : ax);

  // Validate axes match in size
  for (let i = 0; i < axesA.length; i++) {
    if (a.shape[axesA[i]] !== b.shape[axesB[i]]) {
      throw new ValueError(
        `shape mismatch: contracted axes have different lengths: ` +
        `${a.shape[axesA[i]]} vs ${b.shape[axesB[i]]}`
      );
    }
  }

  // Compute output shape
  const outShapeA = a.shape.filter((_, i) => !axesA.includes(i));
  const outShapeB = b.shape.filter((_, i) => !axesB.includes(i));
  const outShape = [...outShapeA, ...outShapeB];

  // Rearrange a: move contracted axes to end
  const aFreeAxes = a.shape.map((_, i) => i).filter(i => !axesA.includes(i));
  const aNewOrder = [...aFreeAxes, ...axesA];
  const aT = a.transpose(aNewOrder);

  // Rearrange b: move contracted axes to beginning
  const bFreeAxes = b.shape.map((_, i) => i).filter(i => !axesB.includes(i));
  const bNewOrder = [...axesB, ...bFreeAxes];
  const bT = b.transpose(bNewOrder);

  // Reshape for matrix multiplication
  const aFreeSize = aFreeAxes.reduce((p, i) => p * a.shape[i], 1);
  const bFreeSize = bFreeAxes.reduce((p, i) => p * b.shape[i], 1);
  const contractSize = axesA.reduce((p, i) => p * a.shape[i], 1);

  const aReshaped = aT.reshape([aFreeSize, contractSize]);
  const bReshaped = bT.reshape([contractSize, bFreeSize]);

  // Matrix multiplication
  const result = matmul(aReshaped, bReshaped);

  // Reshape to output shape
  if (outShape.length === 0) {
    return result.reshape([]);  // Scalar
  }

  return result.reshape(outShape);
}
```

#### 25.1.2 multi_dot

**File:** `src/ts/linalg.ts` (additions)

```typescript
/**
 * Compute the dot product of two or more arrays in a single function call,
 * while automatically selecting the fastest evaluation order.
 *
 * Uses dynamic programming to find the optimal parenthesization that
 * minimizes the total number of scalar multiplications.
 *
 * @param arrays - List of arrays to multiply together
 * @param out - Output array (optional)
 * @returns Product of all arrays
 *
 * @example
 * // For arrays A (10x100), B (100x5), C (5x50)
 * // multi_dot([A, B, C]) is much faster than A @ B @ C
 * // because it computes (A @ B) @ C instead of A @ (B @ C)
 *
 * multi_dot([A, B, C, D])
 */
export function multi_dot(
  arrays: NDArray[],
  out: NDArray | null = null
): NDArray {
  const n = arrays.length;

  if (n < 2) {
    throw new ValueError('multi_dot requires at least 2 arrays');
  }

  if (n === 2) {
    return dot(arrays[0], arrays[1]);
  }

  // Get dimensions for cost computation
  const dims: number[] = [];

  // First array contributes its first dimension
  if (arrays[0].ndim === 1) {
    dims.push(1);
    dims.push(arrays[0].shape[0]);
  } else {
    dims.push(arrays[0].shape[0]);
    dims.push(arrays[0].shape[1]);
  }

  // Middle arrays contribute their last dimension
  for (let i = 1; i < n - 1; i++) {
    dims.push(arrays[i].shape[arrays[i].ndim - 1]);
  }

  // Last array
  if (arrays[n - 1].ndim === 1) {
    dims.push(arrays[n - 1].shape[0]);
  } else {
    dims.push(arrays[n - 1].shape[1]);
  }

  // Use dynamic programming to find optimal order
  const order = _optimalOrder(dims);

  // Execute multiplications in optimal order
  return _executeMultiDot(arrays, order, 0, n - 1);
}

/**
 * Find optimal parenthesization using dynamic programming.
 * Returns split points for recursive multiplication.
 */
function _optimalOrder(dims: number[]): number[][] {
  const n = dims.length - 1;  // Number of matrices

  // m[i][j] = minimum cost to multiply matrices i..j
  const m: number[][] = Array(n).fill(null).map(() => Array(n).fill(0));

  // s[i][j] = optimal split point for matrices i..j
  const s: number[][] = Array(n).fill(null).map(() => Array(n).fill(0));

  // Chain length
  for (let len = 2; len <= n; len++) {
    for (let i = 0; i <= n - len; i++) {
      const j = i + len - 1;
      m[i][j] = Infinity;

      // Try all split points
      for (let k = i; k < j; k++) {
        const cost = m[i][k] + m[k + 1][j] +
                     dims[i] * dims[k + 1] * dims[j + 1];

        if (cost < m[i][j]) {
          m[i][j] = cost;
          s[i][j] = k;
        }
      }
    }
  }

  return s;
}

/**
 * Execute matrix chain multiplication using precomputed order.
 */
function _executeMultiDot(
  arrays: NDArray[],
  s: number[][],
  i: number,
  j: number
): NDArray {
  if (i === j) {
    return arrays[i];
  }

  const k = s[i][j];

  const left = _executeMultiDot(arrays, s, i, k);
  const right = _executeMultiDot(arrays, s, k + 1, j);

  return dot(left, right);
}
```

---

### 25.2 Einstein Summation

#### 25.2.1 einsum

**File:** `src/ts/einsum.ts` (new file)

```typescript
import { NDArray } from './NDArray.js';
import { DType } from './types.js';
import { asarray, zeros, ones, transpose, reshape, sum } from './index.js';
import { tensordot } from './linalg.js';

/**
 * Evaluates the Einstein summation convention on the operands.
 *
 * Using the Einstein summation convention, many common multi-dimensional,
 * linear algebraic array operations can be represented in a simple fashion.
 *
 * @param subscripts - Subscripts string in format "ij,jk->ik"
 * @param operands - Input arrays
 * @param optimize - Optimization strategy ('greedy', 'optimal', true, false)
 * @param dtype - Output dtype
 * @returns Result of the einstein summation
 *
 * @example
 * // Matrix multiplication
 * einsum('ij,jk->ik', A, B)
 *
 * // Trace
 * einsum('ii->', A)
 *
 * // Transpose
 * einsum('ij->ji', A)
 *
 * // Sum over axis
 * einsum('ij->i', A)  // sum over j
 *
 * // Outer product
 * einsum('i,j->ij', a, b)
 *
 * // Batch matrix multiply
 * einsum('bij,bjk->bik', A, B)
 */
export function einsum(
  subscripts: string,
  ...operands: NDArray[]
): NDArray {
  // Parse subscripts
  const { inputSubs, outputSubs, dimensions } = _parseEinsum(subscripts, operands);

  // Validate operands match subscripts
  if (operands.length !== inputSubs.length) {
    throw new ValueError(
      `Number of operands (${operands.length}) doesn't match ` +
      `number of subscript groups (${inputSubs.length})`
    );
  }

  // Build dimension map
  const dimMap = new Map<string, number>();
  for (let i = 0; i < inputSubs.length; i++) {
    const subs = inputSubs[i];
    const arr = operands[i];

    if (subs.length !== arr.ndim) {
      throw new ValueError(
        `Operand ${i} has ${arr.ndim} dimensions but subscripts ` +
        `specify ${subs.length}`
      );
    }

    for (let j = 0; j < subs.length; j++) {
      const label = subs[j];
      const size = arr.shape[j];

      if (dimMap.has(label)) {
        if (dimMap.get(label) !== size) {
          throw new ValueError(
            `Dimension mismatch for label '${label}': ` +
            `${dimMap.get(label)} vs ${size}`
          );
        }
      } else {
        dimMap.set(label, size);
      }
    }
  }

  // Determine output subscripts if not specified
  const finalOutputSubs = outputSubs ?? _implicitOutput(inputSubs);

  // Execute contraction
  return _executeEinsum(operands, inputSubs, finalOutputSubs, dimMap);
}

/**
 * Parse einsum subscripts string.
 */
function _parseEinsum(
  subscripts: string,
  operands: NDArray[]
): {
  inputSubs: string[];
  outputSubs: string | null;
  dimensions: Map<string, number>;
} {
  // Remove whitespace
  const clean = subscripts.replace(/\s/g, '');

  // Split on '->'
  const parts = clean.split('->');

  if (parts.length > 2) {
    throw new ValueError(`Invalid subscripts: more than one '->'`);
  }

  // Parse input subscripts
  const inputPart = parts[0];
  const inputSubs = inputPart.split(',');

  // Parse output subscripts
  const outputSubs = parts.length === 2 ? parts[1] : null;

  // Validate characters
  const validChars = /^[a-zA-Z]*$/;
  for (const sub of inputSubs) {
    if (!validChars.test(sub)) {
      throw new ValueError(`Invalid subscript character in '${sub}'`);
    }
  }

  if (outputSubs !== null && !validChars.test(outputSubs)) {
    throw new ValueError(`Invalid output subscript character in '${outputSubs}'`);
  }

  return { inputSubs, outputSubs, dimensions: new Map() };
}

/**
 * Compute implicit output subscripts (alphabetical order of non-repeated labels).
 */
function _implicitOutput(inputSubs: string[]): string {
  const counts = new Map<string, number>();

  for (const subs of inputSubs) {
    for (const label of subs) {
      counts.set(label, (counts.get(label) ?? 0) + 1);
    }
  }

  // Labels that appear exactly once, sorted alphabetically
  return Array.from(counts.entries())
    .filter(([_, count]) => count === 1)
    .map(([label, _]) => label)
    .sort()
    .join('');
}

/**
 * Execute einsum contraction.
 */
function _executeEinsum(
  operands: NDArray[],
  inputSubs: string[],
  outputSubs: string,
  dimMap: Map<string, number>
): NDArray {
  // Simple case: single operand
  if (operands.length === 1) {
    return _einsumSingle(operands[0], inputSubs[0], outputSubs, dimMap);
  }

  // Two operands: use optimized path
  if (operands.length === 2) {
    return _einsumPair(
      operands[0], inputSubs[0],
      operands[1], inputSubs[1],
      outputSubs, dimMap
    );
  }

  // Multiple operands: reduce pairwise
  let result = operands[0];
  let resultSubs = inputSubs[0];

  for (let i = 1; i < operands.length; i++) {
    // Compute intermediate output subscripts
    const intermediateSubs = i < operands.length - 1
      ? _intermediateOutput(resultSubs, inputSubs[i], outputSubs)
      : outputSubs;

    result = _einsumPair(
      result, resultSubs,
      operands[i], inputSubs[i],
      intermediateSubs, dimMap
    );
    resultSubs = intermediateSubs;
  }

  return result;
}

/**
 * Einsum for single operand (trace, transpose, sum, diagonal).
 */
function _einsumSingle(
  arr: NDArray,
  inputSub: string,
  outputSub: string,
  dimMap: Map<string, number>
): NDArray {
  // Find repeated indices (diagonal/trace)
  const labelPositions = new Map<string, number[]>();
  for (let i = 0; i < inputSub.length; i++) {
    const label = inputSub[i];
    if (!labelPositions.has(label)) {
      labelPositions.set(label, []);
    }
    labelPositions.get(label)!.push(i);
  }

  // Check for trace (repeated index not in output)
  const traceLabels = Array.from(labelPositions.entries())
    .filter(([label, positions]) =>
      positions.length > 1 && !outputSub.includes(label)
    )
    .map(([label, _]) => label);

  // Check for diagonal (repeated index in output)
  const diagLabels = Array.from(labelPositions.entries())
    .filter(([label, positions]) =>
      positions.length > 1 && outputSub.includes(label)
    )
    .map(([label, _]) => label);

  // Check for sum (index not in output)
  const sumLabels = Array.from(labelPositions.entries())
    .filter(([label, positions]) =>
      positions.length === 1 && !outputSub.includes(label)
    )
    .map(([label, _]) => label);

  let result = arr;

  // Handle traces
  for (const label of traceLabels) {
    const positions = labelPositions.get(label)!;
    result = _trace(result, positions[0], positions[1]);
    // Update tracking after trace
  }

  // Handle diagonals
  for (const label of diagLabels) {
    const positions = labelPositions.get(label)!;
    result = _diagonal(result, positions[0], positions[1]);
  }

  // Handle sums
  for (const label of sumLabels) {
    const positions = labelPositions.get(label)!;
    result = sum(result, positions[0], null, true) as NDArray;
    result = squeeze(result, positions[0]);
  }

  // Handle transpose
  const currentOrder = inputSub.split('').filter(l =>
    outputSub.includes(l) &&
    labelPositions.get(l)!.length === 1
  );

  if (currentOrder.join('') !== outputSub) {
    const newOrder = outputSub.split('').map(l =>
      currentOrder.indexOf(l)
    );
    result = transpose(result, newOrder);
  }

  return result;
}

/**
 * Einsum for pair of operands.
 */
function _einsumPair(
  a: NDArray, aSub: string,
  b: NDArray, bSub: string,
  outputSub: string,
  dimMap: Map<string, number>
): NDArray {
  // Find contracted indices (in both inputs, not in output)
  const aLabels = new Set(aSub.split(''));
  const bLabels = new Set(bSub.split(''));

  const contracted = Array.from(aLabels)
    .filter(l => bLabels.has(l) && !outputSub.includes(l));

  // If no contraction, compute outer product
  if (contracted.length === 0) {
    return _outerEinsum(a, aSub, b, bSub, outputSub, dimMap);
  }

  // Find contraction axes
  const axesA = contracted.map(l => aSub.indexOf(l));
  const axesB = contracted.map(l => bSub.indexOf(l));

  // Use tensordot
  let result = tensordot(a, b, [axesA, axesB]);

  // Build result subscripts after tensordot
  const remainingA = aSub.split('').filter(l => !contracted.includes(l));
  const remainingB = bSub.split('').filter(l => !contracted.includes(l));
  const resultSub = [...remainingA, ...remainingB].join('');

  // Transpose if needed
  if (resultSub !== outputSub) {
    const order = outputSub.split('').map(l => resultSub.indexOf(l));
    result = transpose(result, order);
  }

  return result;
}

/**
 * Compute intermediate output subscripts during pairwise reduction.
 */
function _intermediateOutput(
  sub1: string,
  sub2: string,
  finalOutput: string
): string {
  // Keep all labels that are needed for final output or next contraction
  const all = new Set([...sub1.split(''), ...sub2.split('')]);
  const final = new Set(finalOutput.split(''));

  return Array.from(all)
    .filter(l => final.has(l) || sub2.includes(l))
    .sort()
    .join('');
}

/**
 * Outer product for einsum.
 */
function _outerEinsum(
  a: NDArray, aSub: string,
  b: NDArray, bSub: string,
  outputSub: string,
  dimMap: Map<string, number>
): NDArray {
  // Expand dimensions and broadcast
  const result = outer(a.ravel(), b.ravel());

  // Reshape to proper dimensions
  const shape = outputSub.split('').map(l => dimMap.get(l)!);
  return result.reshape(shape);
}

/**
 * Compute trace over two axes.
 */
function _trace(arr: NDArray, axis1: number, axis2: number): NDArray {
  return trace(arr, 0, axis1, axis2);
}

/**
 * Extract diagonal over two axes.
 */
function _diagonal(arr: NDArray, axis1: number, axis2: number): NDArray {
  return diagonal(arr, 0, axis1, axis2);
}
```

#### 25.2.2 einsum_path

**File:** `src/ts/einsum.ts` (additions)

```typescript
/**
 * Evaluates the lowest cost contraction order for an einsum expression.
 *
 * @param subscripts - Subscripts string
 * @param operands - Input arrays
 * @param optimize - Optimization method ('greedy', 'optimal')
 * @returns [path, string_representation]
 *
 * @example
 * const [path, info] = einsum_path('ij,jk,kl->il', A, B, C);
 * // path: [[0, 1], [0, 1]]  // First contract A,B then result,C
 */
export function einsum_path(
  subscripts: string,
  ...operands: NDArray[]
): [number[][], string] {
  const { inputSubs, outputSubs, dimensions } = _parseEinsum(subscripts, operands);

  // Build dimension map
  const dimMap = new Map<string, number>();
  for (let i = 0; i < inputSubs.length; i++) {
    for (let j = 0; j < inputSubs[i].length; j++) {
      dimMap.set(inputSubs[i][j], operands[i].shape[j]);
    }
  }

  // Find optimal path
  const path = _greedyPath(inputSubs, outputSubs ?? _implicitOutput(inputSubs), dimMap);

  // Generate info string
  const info = _pathInfo(subscripts, operands, path, dimMap);

  return [path, info];
}

/**
 * Find contraction path using greedy algorithm.
 */
function _greedyPath(
  inputSubs: string[],
  outputSubs: string,
  dimMap: Map<string, number>
): number[][] {
  const path: number[][] = [];
  const subs = [...inputSubs];
  const indices = inputSubs.map((_, i) => i);

  while (subs.length > 1) {
    // Find best pair to contract
    let bestPair = [0, 1];
    let bestCost = Infinity;

    for (let i = 0; i < subs.length; i++) {
      for (let j = i + 1; j < subs.length; j++) {
        const cost = _contractionCost(subs[i], subs[j], dimMap);
        if (cost < bestCost) {
          bestCost = cost;
          bestPair = [i, j];
        }
      }
    }

    // Record the contraction
    path.push([indices[bestPair[0]], indices[bestPair[1]]]);

    // Update subscripts
    const newSub = _resultSubscripts(subs[bestPair[0]], subs[bestPair[1]], outputSubs);
    subs.splice(bestPair[1], 1);
    subs.splice(bestPair[0], 1);
    subs.push(newSub);

    // Update indices
    const newIdx = Math.max(...indices) + 1;
    indices.splice(bestPair[1], 1);
    indices.splice(bestPair[0], 1);
    indices.push(newIdx);
  }

  return path;
}

/**
 * Estimate cost of contracting two subscript strings.
 */
function _contractionCost(
  sub1: string,
  sub2: string,
  dimMap: Map<string, number>
): number {
  // Cost = product of all dimensions involved
  const allLabels = new Set([...sub1.split(''), ...sub2.split('')]);
  let cost = 1;

  for (const label of allLabels) {
    cost *= dimMap.get(label) ?? 1;
  }

  return cost;
}

/**
 * Compute result subscripts after contracting two operands.
 */
function _resultSubscripts(sub1: string, sub2: string, finalOutput: string): string {
  const labels1 = sub1.split('');
  const labels2 = sub2.split('');

  // Labels in both (contracted)
  const contracted = labels1.filter(l => labels2.includes(l));

  // Remaining labels
  const remaining = [...labels1, ...labels2].filter(l =>
    !contracted.includes(l) || finalOutput.includes(l)
  );

  return [...new Set(remaining)].sort().join('');
}

/**
 * Generate path info string.
 */
function _pathInfo(
  subscripts: string,
  operands: NDArray[],
  path: number[][],
  dimMap: Map<string, number>
): string {
  const lines = [
    `  Complete contraction:  ${subscripts}`,
    `         Naive scaling:  ${_naiveScaling(subscripts, operands)}`,
    `  Optimized scaling:     ${path.length}`,
    ``,
    `  Contraction path:`,
  ];

  for (let i = 0; i < path.length; i++) {
    lines.push(`    ${i}: ${path[i]}`);
  }

  return lines.join('\n');
}

/**
 * Compute naive (left-to-right) scaling exponent.
 */
function _naiveScaling(subscripts: string, operands: NDArray[]): number {
  // Count unique indices
  const allLabels = new Set(subscripts.replace(/[^a-zA-Z]/g, '').split(''));
  return allLabels.size;
}
```

---

### 25.3 Special Products

#### 25.3.1 kron

**File:** `src/ts/linalg.ts` (additions)

```typescript
/**
 * Kronecker product of two arrays.
 *
 * Computes the Kronecker product, a composite array made of blocks of the
 * second array scaled by the first.
 *
 * @param a - First array
 * @param b - Second array
 * @returns Kronecker product
 *
 * @example
 * kron([1, 10, 100], [5, 6, 7])
 * // [5, 6, 7, 50, 60, 70, 500, 600, 700]
 *
 * kron([[1, 2], [3, 4]], [[1, 0], [0, 1]])
 * // [[1, 0, 2, 0],
 * //  [0, 1, 0, 2],
 * //  [3, 0, 4, 0],
 * //  [0, 3, 0, 4]]
 */
export function kron(a: NDArray, b: NDArray): NDArray {
  // Ensure at least 1D
  const aArr = atleast_1d(a);
  const bArr = atleast_1d(b);

  // Match dimensions by prepending 1s to smaller array
  const maxNdim = Math.max(aArr.ndim, bArr.ndim);
  const aShape = [...Array(maxNdim - aArr.ndim).fill(1), ...aArr.shape];
  const bShape = [...Array(maxNdim - bArr.ndim).fill(1), ...bArr.shape];

  const aReshaped = aArr.reshape(aShape);
  const bReshaped = bArr.reshape(bShape);

  // Interleave dimensions for Kronecker product
  // New shape: [a0, b0, a1, b1, ...]
  const interleavedShape: number[] = [];
  for (let i = 0; i < maxNdim; i++) {
    interleavedShape.push(aShape[i]);
    interleavedShape.push(bShape[i]);
  }

  // Expand a and b to interleaved shape and multiply
  const aExpanded = _expandForKron(aReshaped, maxNdim, true);
  const bExpanded = _expandForKron(bReshaped, maxNdim, false);

  const product = multiply(aExpanded, bExpanded);

  // Collapse pairs of dimensions
  const resultShape: number[] = [];
  for (let i = 0; i < maxNdim; i++) {
    resultShape.push(aShape[i] * bShape[i]);
  }

  return product.reshape(resultShape);
}

/**
 * Expand array for Kronecker product computation.
 */
function _expandForKron(arr: NDArray, ndim: number, isFirst: boolean): NDArray {
  // Add dimensions for broadcasting
  const newShape: number[] = [];

  for (let i = 0; i < ndim; i++) {
    newShape.push(arr.shape[i]);
    newShape.push(1);
  }

  if (isFirst) {
    // a goes in positions 0, 2, 4, ...
    return arr.reshape(newShape);
  } else {
    // b goes in positions 1, 3, 5, ...
    const shiftedShape = newShape.map((_, i) =>
      i % 2 === 0 ? 1 : arr.shape[Math.floor(i / 2)]
    );
    return arr.reshape(shiftedShape);
  }
}
```

#### 25.3.2 cross

**File:** `src/ts/linalg.ts` (additions)

```typescript
/**
 * Return the cross product of two (arrays of) vectors.
 *
 * The cross product of a and b in R^3 is a vector perpendicular to both
 * a and b. If a and b are arrays of vectors, the vectors are defined by
 * the last axis (by default).
 *
 * @param a - First vector(s)
 * @param b - Second vector(s)
 * @param axisa - Axis of a that defines the vector(s) (default: -1)
 * @param axisb - Axis of b that defines the vector(s) (default: -1)
 * @param axisc - Axis of c that contains the cross product vectors (default: -1)
 * @param axis - If specified, overrides axisa, axisb, and axisc
 * @returns Cross product vector(s)
 *
 * @example
 * // 3D cross product
 * cross([1, 0, 0], [0, 1, 0])  // [0, 0, 1]
 *
 * // 2D cross product (returns scalar)
 * cross([1, 0], [0, 1])  // 1
 *
 * // Batch cross product
 * cross([[1, 0, 0], [0, 1, 0]], [[0, 1, 0], [0, 0, 1]])
 * // [[0, 0, 1], [1, 0, 0]]
 */
export function cross(
  a: NDArray | ArrayLike<number>,
  b: NDArray | ArrayLike<number>,
  axisa: number = -1,
  axisb: number = -1,
  axisc: number = -1,
  axis: number | null = null
): NDArray {
  const aArr = asarray(a);
  const bArr = asarray(b);

  // Override with axis if provided
  if (axis !== null) {
    axisa = axis;
    axisb = axis;
    axisc = axis;
  }

  // Normalize axes
  axisa = axisa < 0 ? aArr.ndim + axisa : axisa;
  axisb = axisb < 0 ? bArr.ndim + axisb : axisb;

  // Get vector dimensions
  const dimA = aArr.shape[axisa];
  const dimB = bArr.shape[axisb];

  // Validate dimensions
  if (dimA !== dimB) {
    throw new ValueError(
      `incompatible dimensions for cross product: ${dimA} and ${dimB}`
    );
  }

  if (dimA < 2 || dimA > 3) {
    throw new ValueError(
      `cross product only defined for 2D and 3D vectors, got ${dimA}D`
    );
  }

  // Move vector axis to last position
  const aMoved = moveaxis(aArr, axisa, -1);
  const bMoved = moveaxis(bArr, axisb, -1);

  // Broadcast shapes (excluding last axis)
  const batchShapeA = aMoved.shape.slice(0, -1);
  const batchShapeB = bMoved.shape.slice(0, -1);
  const batchShape = broadcastShapes(batchShapeA, batchShapeB);

  // Broadcast arrays
  const aBC = broadcastTo(aMoved, [...batchShape, dimA]);
  const bBC = broadcastTo(bMoved, [...batchShape, dimB]);

  let result: NDArray;

  if (dimA === 2) {
    // 2D cross product: scalar result
    // a × b = a[0]*b[1] - a[1]*b[0]
    const a0 = aBC.slice([...Array(aBC.ndim - 1).fill(null), 0]);
    const a1 = aBC.slice([...Array(aBC.ndim - 1).fill(null), 1]);
    const b0 = bBC.slice([...Array(bBC.ndim - 1).fill(null), 0]);
    const b1 = bBC.slice([...Array(bBC.ndim - 1).fill(null), 1]);

    result = subtract(multiply(a0, b1), multiply(a1, b0));
  } else {
    // 3D cross product: vector result
    // c[0] = a[1]*b[2] - a[2]*b[1]
    // c[1] = a[2]*b[0] - a[0]*b[2]
    // c[2] = a[0]*b[1] - a[1]*b[0]
    const a0 = aBC.slice([...Array(aBC.ndim - 1).fill(null), 0]);
    const a1 = aBC.slice([...Array(aBC.ndim - 1).fill(null), 1]);
    const a2 = aBC.slice([...Array(aBC.ndim - 1).fill(null), 2]);
    const b0 = bBC.slice([...Array(bBC.ndim - 1).fill(null), 0]);
    const b1 = bBC.slice([...Array(bBC.ndim - 1).fill(null), 1]);
    const b2 = bBC.slice([...Array(bBC.ndim - 1).fill(null), 2]);

    const c0 = subtract(multiply(a1, b2), multiply(a2, b1));
    const c1 = subtract(multiply(a2, b0), multiply(a0, b2));
    const c2 = subtract(multiply(a0, b1), multiply(a1, b0));

    result = stack([c0, c1, c2], -1);

    // Move result axis if needed
    axisc = axisc < 0 ? result.ndim + axisc : axisc;
    if (axisc !== result.ndim - 1) {
      result = moveaxis(result, -1, axisc);
    }
  }

  return result;
}
```

---

### 25.4 Tensor Solving

**File:** `src/ts/linalg.ts` (additions)

```typescript
/**
 * Solve the tensor equation a x = b for x.
 *
 * It is assumed that all indices of x are summed over in the product,
 * together with the rightmost indices of a, as is done in, for example,
 * tensordot(a, x, axes=x.ndim).
 *
 * @param a - Coefficient tensor
 * @param b - Right-hand side tensor
 * @param axes - Axes in a to reorder to the right, for the solve
 * @returns Solution tensor x
 *
 * @example
 * const a = eye(4).reshape([2, 2, 2, 2]);
 * const b = fromArray([1, 2, 3, 4]).reshape([2, 2]);
 * const x = tensorsolve(a, b);
 * // tensordot(a, x, 2) ≈ b
 */
export function tensorsolve(
  a: NDArray,
  b: NDArray,
  axes: number[] | null = null
): NDArray {
  // Determine the shape of x
  const Q = a.ndim;
  const N = b.size;

  // Reorder axes if specified
  let aWork = a;
  if (axes !== null) {
    const allAxes = Array.from({ length: Q }, (_, i) => i);
    const newOrder = allAxes.filter(i => !axes.includes(i)).concat(axes);
    aWork = transpose(a, newOrder);
  }

  // The last axes of a must multiply to equal N
  let M = 1;
  const xShape: number[] = [];

  for (let i = b.ndim; i < Q; i++) {
    M *= aWork.shape[i];
    xShape.push(aWork.shape[i]);
  }

  if (M !== N) {
    throw new LinAlgError('Incompatible dimensions for tensorsolve');
  }

  // Reshape to 2D and solve
  const aMatrix = aWork.reshape([N, M]);
  const bVector = b.ravel();

  const xVector = solve(aMatrix, bVector);

  return xVector.reshape(xShape);
}

/**
 * Compute the 'inverse' of an N-dimensional array.
 *
 * The result is an inverse for a with respect to the tensordot operation
 * tensordot(a, b, ind), i.e., up to floating-point accuracy,
 * tensordot(tensorinv(a), a, ind) is the "identity" tensor.
 *
 * @param a - Tensor to pseudo-invert
 * @param ind - Number of first indices that are involved in the inverse sum
 * @returns Tensor inverse
 *
 * @example
 * const a = eye(4).reshape([2, 2, 2, 2]);
 * const ainv = tensorinv(a);
 * // tensordot(ainv, a, ind=2) ≈ eye(4).reshape([2, 2, 2, 2])
 */
export function tensorinv(a: NDArray, ind: number = 2): NDArray {
  // Compute dimensions
  const oldShape = a.shape;
  const prod1 = oldShape.slice(0, ind).reduce((p, x) => p * x, 1);
  const prod2 = oldShape.slice(ind).reduce((p, x) => p * x, 1);

  if (prod1 !== prod2) {
    throw new LinAlgError(
      `Product of first ${ind} dimensions (${prod1}) must equal ` +
      `product of remaining dimensions (${prod2})`
    );
  }

  // Reshape to 2D
  const aMatrix = a.reshape([prod1, prod2]);

  // Compute inverse
  const aInvMatrix = inv(aMatrix);

  // Reshape back
  const newShape = [...oldShape.slice(ind), ...oldShape.slice(0, ind)];
  return aInvMatrix.reshape(newShape);
}
```

---

### 25.5 Norm Functions

**File:** `src/ts/linalg.ts` (additions)

```typescript
/**
 * Compute the matrix norm.
 *
 * This function is able to return one of eight different matrix norms,
 * depending on the value of the ord parameter.
 *
 * @param x - Input matrix (..., M, N)
 * @param ord - Order of the norm (see norm() for details)
 * @param keepdims - If True, axes are left with size one
 * @returns Matrix norm
 *
 * @example
 * matrix_norm([[1, 2], [3, 4]])  // Frobenius norm
 * matrix_norm([[1, 2], [3, 4]], ord=2)  // Spectral norm (largest singular value)
 */
export function matrix_norm(
  x: NDArray,
  ord: number | 'fro' | 'nuc' = 'fro',
  keepdims: boolean = false
): NDArray | number {
  _assertStacked2d(x, 'x');

  const axis = [x.ndim - 2, x.ndim - 1] as [number, number];
  return norm(x, ord, axis, keepdims);
}

/**
 * Compute the vector norm.
 *
 * @param x - Input array
 * @param ord - Order of the norm (default: 2, Euclidean)
 * @param axis - Axis along which to compute. None = flatten
 * @param keepdims - If True, axes are left with size one
 * @returns Vector norm
 *
 * @example
 * vector_norm([3, 4])  // 5 (Euclidean)
 * vector_norm([3, 4], ord=1)  // 7 (Manhattan)
 * vector_norm([3, 4], ord=Infinity)  // 4 (max)
 */
export function vector_norm(
  x: NDArray,
  ord: number = 2,
  axis: number | null = null,
  keepdims: boolean = false
): NDArray | number {
  const arr = asarray(x);

  if (axis === null) {
    return _vectorNormFlat(arr.ravel(), ord);
  }

  return _vectorNormAxis(arr, ord, axis, keepdims);
}

/**
 * Compute vector norm on flattened array.
 */
function _vectorNormFlat(arr: NDArray, ord: number): number {
  const data = arr.toArray();

  if (ord === Infinity) {
    return Math.max(...data.map(Math.abs));
  }

  if (ord === -Infinity) {
    return Math.min(...data.map(Math.abs));
  }

  if (ord === 0) {
    return data.filter(x => x !== 0).length;
  }

  if (ord === 1) {
    return data.reduce((sum, x) => sum + Math.abs(x), 0);
  }

  if (ord === 2) {
    // Use stable algorithm
    return Math.sqrt(data.reduce((sum, x) => sum + x * x, 0));
  }

  // General p-norm
  const pSum = data.reduce((sum, x) => sum + Math.pow(Math.abs(x), ord), 0);
  return Math.pow(pSum, 1 / ord);
}

/**
 * Compute vector norm along axis.
 */
function _vectorNormAxis(
  arr: NDArray,
  ord: number,
  axis: number,
  keepdims: boolean
): NDArray | number {
  const normalizedAxis = axis < 0 ? arr.ndim + axis : axis;

  if (ord === 2) {
    // Euclidean norm: sqrt(sum(x^2))
    const sqSum = sum(multiply(arr, arr), normalizedAxis, null, keepdims);
    return sqrt(sqSum);
  }

  if (ord === 1) {
    // Manhattan norm: sum(|x|)
    return sum(abs(arr), normalizedAxis, null, keepdims);
  }

  if (ord === Infinity) {
    // Max norm: max(|x|)
    return max(abs(arr), normalizedAxis, keepdims);
  }

  if (ord === -Infinity) {
    // Min norm: min(|x|)
    return min(abs(arr), normalizedAxis, keepdims);
  }

  if (ord === 0) {
    // Count non-zero
    return countNonzero(arr, normalizedAxis);
  }

  // General p-norm
  const pSum = sum(power(abs(arr), ord), normalizedAxis, null, keepdims);
  return power(pSum, 1 / ord);
}
```

---

## File Changes Summary

### New Files to Create

```
src/ts/
└── einsum.ts            # Einstein summation

tests/ts/
└── advanced-linalg.test.ts  # Test suite
```

### Files to Modify

```
src/ts/linalg.ts
├── Add tensordot
├── Add multi_dot
├── Add kron
├── Add cross
├── Add tensorsolve
├── Add tensorinv
├── Add matrix_norm
└── Add vector_norm

src/ts/index.ts
├── Export tensordot
├── Export einsum, einsum_path
├── Export multi_dot
├── Export kron
├── Export cross
├── Export tensorsolve
├── Export tensorinv
├── Export matrix_norm
└── Export vector_norm
```

---

## Implementation Order

```
Phase 25.1: Tensor Contraction (Day 1-2)
├── Day 1: tensordot
│   ├── Axis handling
│   ├── Shape computation
│   └── Tests
│
└── Day 2: multi_dot
    ├── Optimal order algorithm
    ├── Chain multiplication
    └── Tests

Phase 25.2: Einstein Summation (Day 3-5)
├── Day 3: einsum parser
│   ├── Subscript parsing
│   ├── Dimension validation
│   └── Tests
│
├── Day 4: einsum execution
│   ├── Single operand (trace, transpose)
│   ├── Pair contraction
│   └── Tests
│
└── Day 5: einsum_path
    ├── Greedy path finding
    ├── Cost estimation
    └── Tests

Phase 25.3: Special Products (Day 6)
├── kron
├── cross
└── Tests

Phase 25.4: Tensor Solving (Day 7)
├── tensorsolve
├── tensorinv
└── Tests

Phase 25.5: Norm Functions (Day 8)
├── matrix_norm
├── vector_norm
└── Tests

Phase 25.6: Polish (Day 9-10)
├── Edge cases
├── NumPy comparison
└── Documentation
```

---

## Verification Plan

After Phase 25 completion, verify:

```bash
# Build
npm run build

# Run tests
npm test

# Phase 25 specific tests:

# tensordot
✓ tensordot(a, b, 1) ≈ matmul(a, b)
✓ tensordot(a, b, 0) = outer product
✓ tensordot with custom axes works

# einsum
✓ einsum('ij,jk->ik', A, B) ≈ matmul(A, B)
✓ einsum('ii->', A) = trace(A)
✓ einsum('ij->ji', A) = transpose(A)
✓ einsum('i,j->ij', a, b) = outer(a, b)
✓ einsum('bij,bjk->bik', A, B) (batch matmul)

# multi_dot
✓ multi_dot([A, B, C]) ≈ A @ B @ C
✓ multi_dot finds optimal order for chain

# kron
✓ kron([[1,2],[3,4]], eye(2)) has correct shape
✓ kron with different-dim arrays works

# cross
✓ cross([1,0,0], [0,1,0]) = [0,0,1]
✓ cross([1,0], [0,1]) = 1 (2D scalar)
✓ Batch cross product works

# tensorsolve/tensorinv
✓ tensordot(a, tensorsolve(a, b)) ≈ b
✓ tensordot(tensorinv(a), a) ≈ identity

# matrix_norm/vector_norm
✓ matrix_norm(A, 'fro') matches NumPy
✓ vector_norm([3,4]) = 5
```

Generate NumPy comparison vectors:

```python
import numpy as np
import json

A = np.array([[1, 2], [3, 4]], dtype=np.float64)
B = np.array([[5, 6], [7, 8]], dtype=np.float64)

tests = {
    "tensordot_1": {
        "a": A.tolist(),
        "b": B.tolist(),
        "axes": 1,
        "expected": np.tensordot(A, B, axes=1).tolist()
    },
    "einsum_matmul": {
        "subscripts": "ij,jk->ik",
        "operands": [A.tolist(), B.tolist()],
        "expected": np.einsum('ij,jk->ik', A, B).tolist()
    },
    "einsum_trace": {
        "subscripts": "ii->",
        "operands": [A.tolist()],
        "expected": float(np.einsum('ii->', A))
    },
    "kron": {
        "a": [[1, 2], [3, 4]],
        "b": [[1, 0], [0, 1]],
        "expected": np.kron([[1, 2], [3, 4]], [[1, 0], [0, 1]]).tolist()
    },
    "cross_3d": {
        "a": [1, 0, 0],
        "b": [0, 1, 0],
        "expected": np.cross([1, 0, 0], [0, 1, 0]).tolist()
    },
    "cross_2d": {
        "a": [1, 0],
        "b": [0, 1],
        "expected": float(np.cross([1, 0], [0, 1]))
    }
}

with open("tests/fixtures/advanced_linalg_vectors.json", "w") as f:
    json.dump(tests, f, indent=2)
```

---

## API Compatibility Notes

### NumPy Signature Match

```typescript
tensordot(a, b, axes=2)
einsum(subscripts, *operands, optimize=False, dtype=None)
einsum_path(subscripts, *operands, optimize='greedy')
multi_dot(arrays, out=None)
kron(a, b)
cross(a, b, axisa=-1, axisb=-1, axisc=-1, axis=None)
tensorsolve(a, b, axes=None)
tensorinv(a, ind=2)
matrix_norm(x, ord='fro', keepdims=False)
vector_norm(x, ord=2, axis=None, keepdims=False)
```

### Differences from NumPy

1. **einsum optimize**: Initially only supports basic optimization, not full 'optimal' mode.

2. **out parameter**: Not supported for most functions initially.

3. **Complex numbers**: einsum with complex conjugation ('...') not supported initially.
