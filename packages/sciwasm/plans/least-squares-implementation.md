# Implementation Plan: scipy.optimize.least_squares

## Overview
Implement `least_squares(fun, x0, options)` for nonlinear least-squares optimization in sciwasm, focusing on the **Trust Region Reflective (TRF)** algorithm with full feature support including bounds, finite difference Jacobians, and robust loss functions.

## Scope

### Included Features
- ✅ TRF algorithm (default scipy method, supports bounds)
- ✅ Bounds support (lower and upper bounds per variable)
- ✅ Finite difference Jacobian approximation (2-point, 3-point, complex-step)
- ✅ Analytical Jacobian support
- ✅ Robust loss functions (linear, soft_l1, huber, cauchy, arctan)
- ✅ Comprehensive result object with cost, residuals, Jacobian, gradient, optimality
- ✅ All standard tolerances (ftol, xtol, gtol)
- ✅ Variable scaling (x_scale, f_scale)

### Excluded/Future Features
- ⬜ MINPACK Levenberg-Marquardt (method='lm') - future enhancement
- ⬜ Dogbox algorithm (method='dogbox') - future enhancement
- ⬜ Sparse Jacobian optimization - future enhancement
- ⬜ Custom loss functions - future enhancement

## Implementation Strategy

### Approach: TypeScript-first with C optimization path

**Rationale:**
1. TRF algorithm is 587 lines of readable Python (trf.py) - easier to port to TypeScript than C initially
2. TypeScript allows faster iteration and debugging during development
3. Can optimize critical sections to C/WASM later if needed
4. scipy's TRF relies on heavy linear algebra (SVD/LSMR) which numwasm already provides

**Architecture:**
```
TypeScript Layer (src/ts/optimize/):
├── least_squares.ts       - Main entry point, dispatcher, validation
├── trf.ts                 - TRF algorithm implementation
├── lsq_common.ts          - Common utilities (bounds, scaling, termination)
├── lsq_numdiff.ts         - Finite difference Jacobian computation
├── lsq_loss.ts            - Robust loss functions
└── lsq_types.ts           - Type definitions

Dependencies:
└── Use numwasm for linear algebra (matrix ops, SVD, norms)
```

## Critical Files to Implement

### 1. Type Definitions
**File:** `src/ts/optimize/lsq_types.ts` (~150 lines)

```typescript
export interface LeastSquaresResult extends OptimizeResult {
  x: number[];           // Solution vector
  cost: number;          // 0.5 * sum(residuals^2)
  fun: number[];         // Final residuals (length m)
  jac: number[][];       // Jacobian m x n matrix
  grad: number[];        // Gradient: J^T @ residuals (length n)
  optimality: number;    // Infinity norm of gradient
  active_mask: number[]; // 0=not active, -1=lower bound, 1=upper bound
  nfev: number;          // Function evaluations
  njev: number | null;   // Jacobian evaluations (null if finite diff)
  status: number;        // 0-4 termination code
  message: string;       // Human-readable termination message
  success: boolean;      // True if status > 0
}

export type ResidualFunction = (x: number[]) => number[];
export type JacobianFunction = (x: number[]) => number[][];

export interface LeastSquaresOptions {
  jac?: JacobianFunction | '2-point' | '3-point' | 'cs';
  bounds?: Bounds | Array<[number | null, number | null]>;
  method?: 'trf';  // Only TRF in v1
  ftol?: number;   // Function tolerance (default: 1e-8)
  xtol?: number;   // Variable tolerance (default: 1e-8)
  gtol?: number;   // Gradient tolerance (default: 1e-8)
  x_scale?: number[] | number | 'jac';  // Variable scaling
  loss?: 'linear' | 'soft_l1' | 'huber' | 'cauchy' | 'arctan';
  f_scale?: number;  // Loss function scale (default: 1.0)
  max_nfev?: number; // Max function evaluations
  diff_step?: number | number[]; // Finite difference step size
  tr_solver?: 'exact' | 'lsmr';  // Trust region solver
  verbose?: 0 | 1 | 2;
}

export interface Bounds {
  lb: number[];  // Lower bounds
  ub: number[];  // Upper bounds
}
```

### 2. Main Entry Point
**File:** `src/ts/optimize/least_squares.ts` (~300 lines)

**Responsibilities:**
- Validate inputs (check dimensions, tolerances, bounds)
- Normalize bounds format (array-of-tuples → Bounds object)
- Handle x_scale preprocessing ('jac' vs numeric)
- Set up VectorFunction wrapper for residuals + Jacobian
- Dispatch to TRF algorithm
- Convert result status codes to messages
- Return LeastSquaresResult

**Key validation logic:**
```typescript
// Check tolerance validity
if (ftol < EPS && xtol < EPS && gtol < EPS) {
  throw new Error('At least one tolerance must be > machine epsilon');
}

// Validate bounds dimensions
if (bounds.lb.length !== n || bounds.ub.length !== n) {
  throw new Error('Bounds dimension mismatch');
}

// Check bounds ordering
for (let i = 0; i < n; i++) {
  if (bounds.lb[i] > bounds.ub[i]) {
    throw new Error(`Lower bound > upper bound at index ${i}`);
  }
}
```

**Pseudocode:**
```typescript
export async function least_squares(
  fun: ResidualFunction,
  x0: number[],
  options?: LeastSquaresOptions
): Promise<LeastSquaresResult> {
  // 1. Extract and validate options
  const method = options?.method ?? 'trf';
  const ftol = options?.ftol ?? 1e-8;
  const xtol = options?.xtol ?? 1e-8;
  const gtol = options?.gtol ?? 1e-8;

  // 2. Prepare bounds
  const bounds = prepareBounds(options?.bounds, x0.length);

  // 3. Set up Jacobian function (analytical or finite diff)
  const jacFn = prepareJacobian(options?.jac, fun, options?.diff_step);

  // 4. Set up loss function
  const lossFn = prepareLossFunction(options?.loss, options?.f_scale);

  // 5. Make initial point feasible
  const x = makeStrictlyFeasible(x0, bounds);

  // 6. Call TRF algorithm
  const result = await trf(
    fun, jacFn, x, bounds, lossFn,
    ftol, xtol, gtol, options
  );

  // 7. Convert status to message
  result.message = TERMINATION_MESSAGES[result.status];
  result.success = result.status > 0;

  return result;
}
```

### 3. TRF Algorithm Implementation
**File:** `src/ts/optimize/trf.ts` (~600 lines)

**Based on:** `/packages/sciwasm/reference/scipy/scipy/optimize/_lsq/trf.py`

**Algorithm overview:**
1. Compute residuals f(x) and Jacobian J(x)
2. Apply loss function if non-linear: ρ(f) → transformed cost
3. Compute gradient: g = J^T @ f
4. Scale variables using diagonal matrix D based on bounds
5. Solve trust-region subproblem: min ||J @ p + f||^2 s.t. ||p|| ≤ Δ
6. Evaluate step quality ratio ρ
7. Update trust radius Δ and accept/reject step
8. Check termination criteria
9. Repeat until convergence

**Key functions to implement:**
```typescript
// Main TRF iteration loop
export async function trf(
  fun: ResidualFunction,
  jac: JacobianFunction,
  x0: number[],
  bounds: Bounds,
  lossFn: LossFunction | null,
  ftol: number, xtol: number, gtol: number,
  options: LeastSquaresOptions
): Promise<LeastSquaresResult>

// Compute scaling matrix D based on bounds and gradient
function computeScaling(x: number[], g: number[], lb: number[], ub: number[]): number[]

// Solve trust-region subproblem using SVD (exact) or LSMR
function solveTrustRegion(
  J: number[][], f: number[], Delta: number, solver: 'exact' | 'lsmr'
): { p: number[], hits_boundary: boolean }

// Evaluate predicted vs actual reduction
function evaluateStepQuality(
  f_old: number[], f_new: number[], J: number[][], p: number[]
): number  // Returns ρ = actual_reduction / predicted_reduction

// Update trust radius based on step quality
function updateTrustRadius(Delta: number, rho: number, ratio: number): number

// Check termination conditions
function checkTermination(
  f: number[], J: number[], g: number[], x: number[], x_old: number[],
  ftol: number, xtol: number, gtol: number
): { status: number, message: string }

// Step size to boundary (used for reflection)
function stepSizeToBound(x: number[], p: number[], lb: number[], ub: number[]): number

// Reflected step if trust region step hits boundary
function selectStep(
  x: number[], J_scaled: number[][], diag: number[],
  g_scaled: number[], p: number[], p_h: number[], Delta: number,
  lb: number[], ub: number[], theta: number
): { step: number[], step_h: number[], predicted_reduction: number }
```

**Critical implementation details:**
- **Variable scaling:** Construct diagonal matrix D where D[i] = (ub[i] - x[i]) if g[i] < 0 and ub[i] is finite, else D[i] = (x[i] - lb[i]) if g[i] > 0 and lb[i] is finite, else D[i] = 1
- **Reflective step:** If trust region step hits boundary, try reflected direction
- **Strict feasibility:** Always maintain x strictly inside bounds using parameter theta (typically 0.995)
- **Trust radius control:** Increase if ρ > 0.75, decrease if ρ < 0.25

### 4. Common Utilities
**File:** `src/ts/optimize/lsq_common.ts` (~200 lines)

**Port from:** `/packages/sciwasm/reference/scipy/scipy/optimize/_lsq/common.py`

**Key functions:**
```typescript
// Machine epsilon constant
export const EPS = 2.220446049250313e-16;

// Make point strictly feasible (move away from bounds)
export function makeStrictlyFeasible(
  x: number[], bounds: Bounds, theta: number = 0.99
): number[]

// Check if point is within bounds
export function inBounds(x: number[], lb: number[], ub: number[]): boolean

// Compute step from x to boundary along direction p
export function stepSizeToBound(
  x: number[], p: number[], lb: number[], ub: number[]
): { step: number, hits_bound: boolean, ind: number }

// Evaluate quadratic function: 0.5 * p^T @ J^T @ J @ p + g^T @ p
export function evaluateQuadratic(
  J: number[][], g: number[], p: number[]
): number

// Compute gradient from Jacobian and residuals
export function computeGrad(J: number[][], f: number[]): number[]

// Compute Jacobian scaling (for x_scale='jac')
export function computeJacScale(J: number[][], scale_inv_old?: number[]): number[]

// CL (Catenary-Length) scaling vector for trust region
export function CLScalingVector(
  x: number[], g: number[], lb: number[], ub: number[]
): number[]

// Check termination conditions
export function checkTermination(
  dF: number, F: number, dx_norm: number, x_norm: number, ratio: number,
  ftol: number, xtol: number
): number  // Returns status code

// Build 1D quadratic function along direction
export function buildQuadratic1d(
  J: number[][], g: number[], p: number[]
): { a: number, b: number, c: number }

// Minimize 1D quadratic on interval [0, t]
export function minimizeQuadratic1d(
  a: number, b: number, lb: number, ub: number
): number
```

### 5. Finite Difference Jacobian
**File:** `src/ts/optimize/lsq_numdiff.ts` (~250 lines)

**Port from:** `/packages/sciwasm/reference/scipy/scipy/optimize/_numdiff.py`

**Methods to implement:**
```typescript
// Main entry point for finite difference approximation
export function approximateJacobian(
  fun: ResidualFunction,
  x: number[],
  f0: number[],
  method: '2-point' | '3-point' | 'cs',
  rel_step?: number | number[],
  abs_step?: number | number[]
): { J: number[][], nfev: number }

// 2-point forward differences: J[:, i] ≈ (f(x + h*e_i) - f0) / h
function approxJacobian2Point(
  fun: ResidualFunction, x: number[], f0: number[], h: number[]
): { J: number[][], nfev: number }

// 3-point central differences: J[:, i] ≈ (f(x + h*e_i) - f(x - h*e_i)) / (2*h)
function approxJacobian3Point(
  fun: ResidualFunction, x: number[], h: number[]
): { J: number[][], nfev: number }

// Complex-step derivatives: J[:, i] ≈ Im(f(x + ih*e_i)) / h
function approxJacobianCS(
  fun: ResidualFunction, x: number[], h: number[]
): { J: number[][], nfev: number }

// Compute relative step sizes
export function computeRelativeStep(
  x: number[], method: '2-point' | '3-point' | 'cs',
  rel_step?: number | number[]
): number[]

// Convert relative to absolute step
function computeAbsoluteStep(
  rel_step: number[], x: number[], f0: number[], method: string
): number[]
```

**Default step sizes:**
- 2-point: eps^(1/2) ≈ 1.49e-8
- 3-point: eps^(1/3) ≈ 6.06e-6
- complex-step: eps ≈ 2.22e-16 (machine epsilon)

### 6. Robust Loss Functions
**File:** `src/ts/optimize/lsq_loss.ts` (~200 lines)

**Port from:** `/packages/sciwasm/reference/scipy/scipy/optimize/_lsq/least_squares.py` (lines 183-252)

**Loss function interface:**
```typescript
export interface LossFunction {
  // Transform residuals and compute cost
  apply(f: number[], f_scale: number): {
    rho: number;  // Cost: sum(ρ(z))
    rho_1: number[];  // First derivative: ρ'(z) for each residual
    rho_2: number[];  // Second derivative: ρ''(z) for each residual
  };
}

// Apply loss to transform Jacobian: J_transformed = diag(ρ'(z)^0.5) @ J
export function scaleForRobustLossFunction(
  J: number[][], f: number[], rho_1: number[]
): number[][]
```

**Loss implementations:**
```typescript
// Linear (no transformation): ρ(z) = z
export class LinearLoss implements LossFunction

// Soft L1: ρ(z) = 2((1 + z)^0.5 - 1)
export class SoftL1Loss implements LossFunction

// Huber: ρ(z) = z if z ≤ 1 else 2*z^0.5 - 1
export class HuberLoss implements LossFunction

// Cauchy: ρ(z) = ln(1 + z)
export class CauchyLoss implements LossFunction

// Arctan: ρ(z) = arctan(z)
export class ArctanLoss implements LossFunction

// Factory function
export function createLossFunction(
  loss: string, f_scale: number
): LossFunction | null
```

**Loss transformation math:**
Given residuals f, compute z = (f / f_scale)^2, then:
- Cost = 0.5 * f_scale^2 * sum(ρ(z))
- Modified Jacobian: J_loss = diag(ρ'(z)^0.5) @ J
- Modified gradient: g_loss = J_loss^T @ f_loss
where f_loss = ρ'(z)^0.5 * f

### 7. Trust Region Subproblem Solvers

**SVD-based exact solver** (for small/medium problems, n < ~100):
```typescript
// Solve: min ||J @ p + f||^2 subject to ||p|| ≤ Delta
// Using SVD of J = U @ S @ V^T
function solveLsqTrustRegionExact(
  J: number[][], f: number[], Delta: number
): { p: number[], hits_boundary: boolean }
```

**Implementation:**
1. Compute SVD: J = U @ diag(s) @ V^T using numwasm's SVD
2. Transform: g_svd = V^T @ J^T @ f = V^T @ U^T @ (U @ s @ V^T)^T @ f
3. If unconstrained solution ||p|| ≤ Delta, return it
4. Otherwise, solve secular equation to find lambda: sum((s[i] * g_svd[i] / (s[i]^2 + lambda))^2) = Delta^2
5. Compute p = -V @ diag(s / (s^2 + lambda)) @ U^T @ f

**LSMR-based solver** (for large problems):
```typescript
// Use iterative LSMR from scipy.sparse.linalg
// Solve in 2D subspace spanned by gradient and Gauss-Newton direction
function solveLsqTrustRegionLSMR(
  J: number[][], f: number[], Delta: number
): { p: number[], hits_boundary: boolean }
```

## Dependencies and Integration

### NumWASM Dependencies
- **Matrix operations:** Use NumJS NDArray for matrix-vector products
- **SVD:** Use `numwasm.linalg.svd()` for exact trust region solver
- **Norms:** Use `numwasm.linalg.norm()` for vector norms
- **Linear solve:** May need `numwasm.linalg.lstsq()` for least squares

### Build System Changes
**File:** `scripts/build-wasm.sh`
- No WASM changes needed for TypeScript implementation
- Future: If we optimize to C, add trf.c compilation

**File:** `package.json`
- No new npm dependencies needed

## Testing Strategy

### Test File Structure
**File:** `tests/ts/optimize/least_squares.test.ts` (~400 lines)

**Port from:** `/packages/sciwasm/reference/scipy/scipy/optimize/tests/test_least_squares.py`

### Test Categories

**1. Basic Functionality Tests**
```typescript
describe('least_squares - basic', () => {
  test('trivial quadratic residuals', async () => {
    // f(x) = x^2, should converge to x = 0
    const fun = (x: number[]) => x.map(xi => xi * xi);
    const result = await least_squares(fun, [1.0]);
    expect(result.success).toBe(true);
    expect(result.x[0]).toBeCloseTo(0, 6);
  });

  test('Rosenbrock function', async () => {
    // f(x) = [10*(x[1] - x[0]^2), 1 - x[0]]
    const fun = (x: number[]) => [
      10 * (x[1] - x[0] * x[0]),
      1 - x[0]
    ];
    const result = await least_squares(fun, [-1, -1]);
    expect(result.success).toBe(true);
    expect(result.x[0]).toBeCloseTo(1, 3);
    expect(result.x[1]).toBeCloseTo(1, 3);
  });
});
```

**2. Bounds Tests**
```typescript
describe('least_squares - bounds', () => {
  test('respects lower bounds', async () => {
    const fun = (x: number[]) => x.map(xi => xi * xi);
    const result = await least_squares(fun, [5, 5], {
      bounds: { lb: [2, 2], ub: [10, 10] }
    });
    expect(result.success).toBe(true);
    expect(result.x[0]).toBeCloseTo(2, 5);  // At lower bound
    expect(result.active_mask[0]).toBe(-1);  // Lower bound active
  });

  test('respects upper bounds', async () => {
    const fun = (x: number[]) => x.map(xi => (xi - 5) * (xi - 5));
    const result = await least_squares(fun, [0, 0], {
      bounds: { lb: [-10, -10], ub: [2, 2] }
    });
    expect(result.x[0]).toBeCloseTo(2, 5);  // At upper bound
  });
});
```

**3. Jacobian Tests**
```typescript
describe('least_squares - Jacobian', () => {
  test('analytical Jacobian', async () => {
    const fun = (x: number[]) => [10 * (x[1] - x[0]**2), 1 - x[0]];
    const jac = (x: number[]) => [[-20 * x[0], 10], [-1, 0]];

    const result = await least_squares(fun, [-1, -1], { jac });
    expect(result.success).toBe(true);
    expect(result.njev).toBeGreaterThan(0);  // Analytical jac was used
  });

  test('2-point finite differences', async () => {
    const fun = (x: number[]) => [10 * (x[1] - x[0]**2), 1 - x[0]];

    const result = await least_squares(fun, [-1, -1], { jac: '2-point' });
    expect(result.success).toBe(true);
    expect(result.njev).toBeNull();  // Finite diff doesn't count
    expect(result.nfev).toBeGreaterThan(10);  // Extra evals for FD
  });
});
```

**4. Loss Function Tests**
```typescript
describe('least_squares - loss functions', () => {
  test('soft_l1 loss with outliers', async () => {
    // Data with outliers
    const xdata = [0, 1, 2, 3, 4];
    const ydata = [0, 2, 4, 100, 8];  // y[3] is outlier

    const fun = (p: number[]) => xdata.map((x, i) =>
      p[0] * x + p[1] - ydata[i]
    );

    const result = await least_squares(fun, [0, 0], {
      loss: 'soft_l1'
    });

    // Should fit line ignoring outlier
    expect(result.x[0]).toBeCloseTo(2, 1);  // slope ≈ 2
    expect(result.x[1]).toBeCloseTo(0, 1);  // intercept ≈ 0
  });
});
```

**5. Tolerance Tests**
```typescript
describe('least_squares - tolerances', () => {
  test('ftol termination', async () => {
    const fun = (x: number[]) => [x[0] * x[0]];
    const result = await least_squares(fun, [1], {
      ftol: 1e-4, xtol: 0, gtol: 0
    });
    expect(result.status).toBe(2);  // ftol satisfied
  });

  test('gtol termination', async () => {
    const fun = (x: number[]) => [x[0] * x[0]];
    const result = await least_squares(fun, [1], {
      ftol: 0, xtol: 0, gtol: 1e-4
    });
    expect(result.status).toBe(1);  // gtol satisfied
  });
});
```

**6. Edge Cases and Error Handling**
```typescript
describe('least_squares - edge cases', () => {
  test('max_nfev exceeded', async () => {
    const fun = (x: number[]) => [x[0] * x[0]];
    const result = await least_squares(fun, [1], { max_nfev: 5 });
    expect(result.status).toBe(0);  // Max iterations
    expect(result.success).toBe(false);
  });

  test('invalid bounds throw error', () => {
    const fun = (x: number[]) => [x[0]];
    expect(() =>
      least_squares(fun, [1], { bounds: { lb: [5], ub: [2] } })
    ).toThrow(/lower bound.*upper bound/i);
  });

  test('dimension mismatch throws error', () => {
    const fun = (x: number[]) => [x[0], x[1]];
    expect(() =>
      least_squares(fun, [1], {})
    ).toThrow(/dimension/i);
  });
});
```

### Validation Against SciPy

**Create reference outputs:**
```python
# Generate test fixtures from scipy
import numpy as np
from scipy.optimize import least_squares
import json

def rosenbrock(x):
    return [10 * (x[1] - x[0]**2), 1 - x[0]]

result = least_squares(rosenbrock, [-1, -1])

fixture = {
    'x': result.x.tolist(),
    'cost': result.cost,
    'fun': result.fun.tolist(),
    'jac': result.jac.tolist(),
    'optimality': result.optimality,
    'nfev': result.nfev,
    'status': result.status
}

with open('rosenbrock_reference.json', 'w') as f:
    json.dump(fixture, f)
```

**Compare in tests:**
```typescript
test('matches scipy reference for Rosenbrock', async () => {
  const reference = require('./fixtures/rosenbrock_reference.json');
  const result = await least_squares(rosenbrock, [-1, -1]);

  expect(result.x[0]).toBeCloseTo(reference.x[0], 5);
  expect(result.x[1]).toBeCloseTo(reference.x[1], 5);
  expect(result.cost).toBeCloseTo(reference.cost, 8);
  expect(result.optimality).toBeCloseTo(reference.optimality, 8);
});
```

## Implementation Sequence

### Phase 1: Core Infrastructure (Days 1-2)
1. ✅ Create type definitions (lsq_types.ts)
2. ✅ Implement common utilities (lsq_common.ts)
3. ✅ Write basic tests for utilities
4. ✅ Set up test fixtures from scipy

### Phase 2: Finite Difference Jacobian (Days 3-4)
1. ✅ Implement 2-point finite differences (lsq_numdiff.ts)
2. ✅ Add 3-point and complex-step methods
3. ✅ Test against scipy's approx_derivative
4. ✅ Validate step size selection logic

### Phase 3: Loss Functions (Day 5)
1. ✅ Implement all 5 loss functions (lsq_loss.ts)
2. ✅ Test loss transformations
3. ✅ Verify derivatives (ρ', ρ'')

### Phase 4: TRF Algorithm - Core (Days 6-9)
1. ✅ Implement main TRF iteration loop (trf.ts)
2. ✅ Compute variable scaling matrix D
3. ✅ Implement trust region subproblem solver (SVD method)
4. ✅ Add step selection logic (trust region vs reflected step)
5. ✅ Test on simple problems without bounds

### Phase 5: TRF Algorithm - Bounds (Days 10-11)
1. ✅ Add bounds handling in scaling
2. ✅ Implement strict feasibility enforcement
3. ✅ Add active set detection
4. ✅ Test with various bound configurations

### Phase 6: Main Interface (Day 12)
1. ✅ Implement least_squares.ts dispatcher
2. ✅ Add input validation
3. ✅ Connect all components
4. ✅ Export from optimize/index.ts

### Phase 7: Testing & Validation (Days 13-15)
1. ✅ Port all relevant scipy tests
2. ✅ Create comprehensive test suite
3. ✅ Validate against scipy reference outputs
4. ✅ Fix bugs and edge cases
5. ✅ Performance benchmarking

### Phase 8: Documentation (Day 16)
1. ✅ Write API documentation
2. ✅ Add usage examples
3. ✅ Document limitations (no dogbox/lm yet)
4. ✅ Update package README

## Critical Dependencies on Reference Code

### Must Copy/Port These Files:

1. **`/packages/sciwasm/reference/scipy/scipy/optimize/_lsq/trf.py`** (587 lines)
   - Main TRF algorithm implementation
   - Trust region iteration loop
   - Variable scaling logic
   - Step selection (trust region vs reflected)

2. **`/packages/sciwasm/reference/scipy/scipy/optimize/_lsq/common.py`** (731 lines)
   - Utility functions (step_size_to_bound, in_bounds, etc.)
   - Termination checking
   - Quadratic model evaluation
   - Scaling computations

3. **`/packages/sciwasm/reference/scipy/scipy/optimize/_numdiff.py`** (from main scipy)
   - Finite difference Jacobian approximation
   - Step size selection
   - 2-point, 3-point, complex-step methods

4. **`/packages/sciwasm/reference/scipy/scipy/optimize/_lsq/least_squares.py`** (1044 lines)
   - Loss function implementations (lines 183-252)
   - Parameter validation logic
   - Bounds preprocessing
   - Main dispatcher

5. **`/packages/sciwasm/reference/scipy/scipy/optimize/tests/test_least_squares.py`** (999 lines)
   - Test functions (rosenbrock, bvp, etc.)
   - Comprehensive test cases
   - Expected behaviors and tolerances

### Linear Algebra Dependencies from NumWASM:
- `numwasm.linalg.svd()` - Singular value decomposition
- `numwasm.linalg.norm()` - Vector norms
- `numwasm.linalg.lstsq()` - Least squares solve (for LSMR alternative)
- Matrix-vector products (can use NumJS NDArray)

## Trade-offs and Risks

### Advantages of This Approach
- ✅ **Bounds support from day 1** - Most important scipy feature
- ✅ **Finite differences** - Users don't need analytical Jacobians
- ✅ **Robust loss functions** - Handles outliers, matches scipy capabilities
- ✅ **TypeScript implementation** - Easier debugging, faster iteration
- ✅ **TRF is scipy's default** - Most commonly used algorithm

### Limitations
- ⚠️ **No MINPACK LM** - Slightly less efficient for unconstrained problems
- ⚠️ **TypeScript may be slower than C** - Can optimize later if needed
- ⚠️ **No sparse Jacobian optimization** - May be slow for very large problems
- ⚠️ **SVD for trust region** - Limited to medium-scale problems (n < 1000)

### Risks and Mitigation

**Risk 1: Performance**
- TypeScript implementation may be 2-10x slower than scipy's optimized C/Fortran
- **Mitigation:** Profile and optimize critical sections, consider porting hot paths to WASM later

**Risk 2: Numerical Stability**
- Porting math-heavy algorithms can introduce subtle bugs
- **Mitigation:** Extensive testing against scipy reference outputs, validate on benchmark problems

**Risk 3: Linear Algebra Dependencies**
- TRF requires SVD which must be robust
- **Mitigation:** Validate numwasm SVD implementation, fall back to LSMR if issues arise

**Risk 4: Complexity of TRF**
- 587 lines of non-trivial math, reflective steps, scaling matrices
- **Mitigation:** Port incrementally, test each component, extensive comments in code

**Risk 5: Test Coverage**
- Need to port ~999 lines of Python tests
- **Mitigation:** Focus on high-value tests first, add others iteratively

## Success Criteria

### Minimum Viable Product (MVP)
- ✅ `least_squares()` works for Rosenbrock function
- ✅ Supports bounds (lower, upper, mixed)
- ✅ Finite difference Jacobian (2-point)
- ✅ Linear loss function
- ✅ Converges within 2x function evaluations of scipy
- ✅ 10+ tests passing

### Full Feature Set
- ✅ All 5 loss functions implemented and tested
- ✅ 2-point, 3-point, complex-step finite differences
- ✅ Analytical Jacobian support
- ✅ x_scale and f_scale working
- ✅ All tolerance termination conditions
- ✅ 50+ tests passing from scipy suite
- ✅ Results match scipy within 1e-6 relative error

### Production Ready
- ✅ Full test coverage (>90%)
- ✅ API documentation complete
- ✅ Performance benchmarked vs scipy
- ✅ Edge cases handled gracefully
- ✅ Examples and tutorials written

## Future Enhancements (Post-MVP)

### Priority 1: MINPACK Levenberg-Marquardt
- Add method='lm' for unconstrained problems
- Compile minpack.c to WASM
- ~3-5 days of work

### Priority 2: Dogbox Algorithm
- Add method='dogbox' as alternative to TRF
- Port dogbox.py (345 lines)
- ~2-3 days of work

### Priority 3: Sparse Jacobian Optimization
- Support scipy.sparse matrices for Jacobian
- Use LSMR for large-scale problems
- Column grouping for finite differences
- ~5-7 days of work

### Priority 4: Performance Optimization
- Profile TypeScript TRF
- Port hot paths to C/WASM if needed
- Optimize memory allocations
- ~3-5 days of work

## Verification Plan

All implementation can be verified by running tests:

```bash
cd /workspaces/numwasm/packages/sciwasm
pnpm test tests/ts/optimize/least_squares.test.ts
```

Expected outputs:
- All basic tests pass (Rosenbrock, trivial problems)
- Bounds tests pass (lower, upper, mixed)
- Finite difference Jacobian matches scipy
- Loss functions produce correct costs
- Termination criteria work correctly
- Results match scipy reference within tolerance

Performance benchmark:
```bash
pnpm bench least_squares
```

Compare against scipy timings for Rosenbrock (should be within 3x).
