/**
 * Tests for BFGS optimization.
 *
 * Ported from scipy/optimize/tests/test_optimize.py
 */

import { describe, it, expect } from 'vitest';
import { minimize } from '../../../src/ts/optimize/minimize.js';
import { minimizeBFGS } from '../../../src/ts/optimize/bfgs.js';
import { rosenbrock, booth, quadratic, himmelblau } from './_test_functions.js';

describe('minimize with BFGS', () => {
  it('should minimize a simple quadratic', async () => {
    const result = await minimize(quadratic, [1, 2, 3], {
      method: 'BFGS',
    });

    expect(result.success).toBe(true);
    expect(result.status).toBe(0);
    expect(result.fun).toBeLessThan(1e-10);
    for (const xi of result.x) {
      expect(Math.abs(xi)).toBeLessThan(1e-4);
    }
  });

  it('should minimize the Rosenbrock function (2D)', async () => {
    const result = await minimize(rosenbrock, [-1, -1], {
      method: 'BFGS',
      options: { gtol: 1e-6 },
    });

    expect(result.success).toBe(true);
    expect(result.fun).toBeLessThan(1e-8);
    expect(Math.abs(result.x[0] - 1)).toBeLessThan(1e-3);
    expect(Math.abs(result.x[1] - 1)).toBeLessThan(1e-3);
  });

  it('should minimize the Booth function', async () => {
    const result = await minimize(booth, [0, 0], {
      method: 'BFGS',
      options: { gtol: 1e-6 },
    });

    expect(result.success).toBe(true);
    expect(result.fun).toBeLessThan(1e-8);
    expect(Math.abs(result.x[0] - 1)).toBeLessThan(1e-3);
    expect(Math.abs(result.x[1] - 3)).toBeLessThan(1e-3);
  });

  it('should minimize Himmelblau function', async () => {
    const result = await minimize(himmelblau, [0, 0], {
      method: 'BFGS',
      options: { gtol: 1e-6 },
    });

    expect(result.success).toBe(true);
    expect(result.fun).toBeLessThan(1e-6);
  });

  it('should return jac and hess_inv', async () => {
    const result = await minimize(quadratic, [1, 2], {
      method: 'BFGS',
    });

    expect(result.jac).toBeDefined();
    expect(result.jac!.length).toBe(2);
    expect(result.hess_inv).toBeDefined();
    expect(result.hess_inv!.length).toBe(2);
    expect(result.hess_inv![0].length).toBe(2);
    // Gradient should be near zero at optimum
    for (const gi of result.jac!) {
      expect(Math.abs(gi)).toBeLessThan(1e-3);
    }
  });

  it('should count function and gradient evaluations', async () => {
    const result = await minimize(quadratic, [1, 2], {
      method: 'BFGS',
    });

    expect(result.nfev).toBeGreaterThan(0);
    expect(result.njev).toBeGreaterThan(0);
    expect(result.nit).toBeGreaterThan(0);
  });

  it('should respect maxiter limit', async () => {
    const result = await minimize(rosenbrock, [-1, -1], {
      method: 'BFGS',
      options: { maxiter: 2 },
    });

    expect(result.status).toBe(1); // maxiter exceeded
    expect(result.success).toBe(false);
    expect(result.nit).toBeLessThanOrEqual(2);
  });

  it('should work with user-provided Jacobian', async () => {
    // Rosenbrock gradient
    const rosenGrad = (x: number[]): number[] => {
      const [x0, x1] = x;
      return [
        -400 * x0 * (x1 - x0 * x0) + 2 * (x0 - 1),
        200 * (x1 - x0 * x0),
      ];
    };

    const result = await minimize(rosenbrock, [-1, -1], {
      method: 'BFGS',
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
      method: 'BFGS',
    });

    expect(result.success).toBe(true);
    expect(Math.abs(result.x[0] - 3)).toBeLessThan(1e-4);
  });

  it('should work via dispatcher as default (no bounds)', async () => {
    const result = await minimize(quadratic, [5, -3]);

    expect(result.success).toBe(true);
    expect(result.fun).toBeLessThan(1e-8);
  });
});

describe('minimizeBFGS direct', () => {
  it('should work when called directly', async () => {
    const result = await minimizeBFGS(quadratic, [1, 2]);

    expect(result.success).toBe(true);
    expect(result.fun).toBeLessThan(1e-4);
  });
});
