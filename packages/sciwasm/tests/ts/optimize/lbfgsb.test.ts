/**
 * Tests for L-BFGS-B optimization.
 *
 * Ported from scipy/optimize/tests/test_optimize.py
 */

import { describe, it, expect } from 'vitest';
import { minimize } from '../../../src/ts/optimize/minimize.js';
import { minimizeLBFGSB } from '../../../src/ts/optimize/lbfgsb.js';
import { rosenbrock, booth, quadratic, himmelblau } from './_test_functions.js';

describe('minimize with L-BFGS-B', () => {
  it('should minimize a simple quadratic (unbounded)', async () => {
    const result = await minimize(quadratic, [1, 2, 3], {
      method: 'L-BFGS-B',
    });

    expect(result.success).toBe(true);
    expect(result.fun).toBeLessThan(1e-10);
    for (const xi of result.x) {
      expect(Math.abs(xi)).toBeLessThan(1e-4);
    }
  });

  it('should minimize the Rosenbrock function (2D)', async () => {
    const result = await minimize(rosenbrock, [-1, -1], {
      method: 'L-BFGS-B',
      options: { gtol: 1e-6 },
    });

    expect(result.success).toBe(true);
    expect(result.fun).toBeLessThan(1e-8);
    expect(Math.abs(result.x[0] - 1)).toBeLessThan(1e-3);
    expect(Math.abs(result.x[1] - 1)).toBeLessThan(1e-3);
  });

  it('should minimize the Booth function', async () => {
    const result = await minimize(booth, [0, 0], {
      method: 'L-BFGS-B',
      options: { gtol: 1e-6 },
    });

    expect(result.success).toBe(true);
    expect(result.fun).toBeLessThan(1e-8);
    expect(Math.abs(result.x[0] - 1)).toBeLessThan(1e-3);
    expect(Math.abs(result.x[1] - 3)).toBeLessThan(1e-3);
  });

  it('should minimize Himmelblau function', async () => {
    const result = await minimize(himmelblau, [0, 0], {
      method: 'L-BFGS-B',
      options: { gtol: 1e-6 },
    });

    expect(result.success).toBe(true);
    expect(result.fun).toBeLessThan(1e-6);
  });

  it('should respect bounds', async () => {
    // Minimize quadratic with bounds that exclude the true minimum (origin)
    const result = await minimize(quadratic, [5, 5], {
      method: 'L-BFGS-B',
      bounds: {
        lb: [2, 2],
        ub: [10, 10],
      },
    });

    expect(result.success).toBe(true);
    // Solution should be at the lower bound
    expect(Math.abs(result.x[0] - 2)).toBeLessThan(1e-3);
    expect(Math.abs(result.x[1] - 2)).toBeLessThan(1e-3);
  });

  it('should handle lower-bound only', async () => {
    const result = await minimize(quadratic, [5, 5], {
      method: 'L-BFGS-B',
      bounds: {
        lb: [3, 3],
        ub: [Infinity, Infinity],
      },
    });

    expect(result.success).toBe(true);
    expect(Math.abs(result.x[0] - 3)).toBeLessThan(1e-3);
    expect(Math.abs(result.x[1] - 3)).toBeLessThan(1e-3);
  });

  it('should handle upper-bound only', async () => {
    // f(x) = (x[0]-10)^2 + (x[1]-10)^2, min at (10,10)
    const f = (x: number[]) => (x[0] - 10) ** 2 + (x[1] - 10) ** 2;
    const result = await minimize(f, [0, 0], {
      method: 'L-BFGS-B',
      bounds: {
        lb: [-Infinity, -Infinity],
        ub: [5, 5],
      },
    });

    expect(result.success).toBe(true);
    // Minimum should be at upper bound
    expect(Math.abs(result.x[0] - 5)).toBeLessThan(1e-3);
    expect(Math.abs(result.x[1] - 5)).toBeLessThan(1e-3);
  });

  it('should count function and gradient evaluations', async () => {
    const result = await minimize(quadratic, [1, 2], {
      method: 'L-BFGS-B',
    });

    expect(result.nfev).toBeGreaterThan(0);
    expect(result.njev).toBeGreaterThan(0);
    expect(result.nit).toBeGreaterThan(0);
  });

  it('should return jac at solution', async () => {
    const result = await minimize(quadratic, [1, 2], {
      method: 'L-BFGS-B',
    });

    expect(result.jac).toBeDefined();
    expect(result.jac!.length).toBe(2);
    // Gradient should be near zero at optimum
    for (const gi of result.jac!) {
      expect(Math.abs(gi)).toBeLessThan(1e-3);
    }
  });

  it('should work with user-provided Jacobian', async () => {
    const rosenGrad = (x: number[]): number[] => {
      const [x0, x1] = x;
      return [
        -400 * x0 * (x1 - x0 * x0) + 2 * (x0 - 1),
        200 * (x1 - x0 * x0),
      ];
    };

    const result = await minimize(rosenbrock, [-1, -1], {
      method: 'L-BFGS-B',
      jac: rosenGrad,
      options: { gtol: 1e-6 },
    });

    expect(result.success).toBe(true);
    expect(result.fun).toBeLessThan(1e-8);
    expect(Math.abs(result.x[0] - 1)).toBeLessThan(1e-3);
    expect(Math.abs(result.x[1] - 1)).toBeLessThan(1e-3);
  });

  it('should handle 1D optimization', async () => {
    const f = (x: number[]) => (x[0] - 3) ** 2;
    const result = await minimize(f, [0], {
      method: 'L-BFGS-B',
    });

    expect(result.success).toBe(true);
    expect(Math.abs(result.x[0] - 3)).toBeLessThan(1e-4);
  });

  it('should be the default method when bounds are provided', async () => {
    const result = await minimize(quadratic, [5, 5], {
      bounds: {
        lb: [2, 2],
        ub: [10, 10],
      },
    });

    expect(result.success).toBe(true);
    expect(Math.abs(result.x[0] - 2)).toBeLessThan(1e-3);
    expect(Math.abs(result.x[1] - 2)).toBeLessThan(1e-3);
  });
});

describe('minimizeLBFGSB direct', () => {
  it('should work when called directly', async () => {
    const result = await minimizeLBFGSB(quadratic, [1, 2]);

    expect(result.success).toBe(true);
    expect(result.fun).toBeLessThan(1e-4);
  });

  it('should work with bounds', async () => {
    const bounds = { lb: [1, 1], ub: [10, 10] };
    const result = await minimizeLBFGSB(quadratic, [5, 5], undefined, bounds);

    expect(result.success).toBe(true);
    expect(Math.abs(result.x[0] - 1)).toBeLessThan(1e-3);
    expect(Math.abs(result.x[1] - 1)).toBeLessThan(1e-3);
  });
});
