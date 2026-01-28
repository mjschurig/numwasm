/**
 * Tests for Nelder-Mead optimization.
 *
 * Ported from scipy/optimize/tests/test_optimize.py
 */

import { describe, it, expect } from 'vitest';
import { minimize } from '../../../src/ts/optimize/minimize.js';
import { minimizeNelderMead } from '../../../src/ts/optimize/nelder_mead.js';
import { rosenbrock, booth, quadratic, himmelblau } from './_test_functions.js';

describe('minimize with Nelder-Mead', () => {
  it('should minimize a simple quadratic', async () => {
    const result = await minimize(quadratic, [1, 2, 3], {
      method: 'Nelder-Mead',
    });

    expect(result.success).toBe(true);
    expect(result.status).toBe(0);
    expect(result.fun).toBeLessThan(1e-8);
    for (const xi of result.x) {
      expect(Math.abs(xi)).toBeLessThan(1e-4);
    }
  });

  it('should minimize the Rosenbrock function (2D)', async () => {
    const result = await minimize(rosenbrock, [-1, -1], {
      method: 'Nelder-Mead',
      options: { xatol: 1e-8, fatol: 1e-8 },
    });

    expect(result.success).toBe(true);
    expect(result.fun).toBeLessThan(1e-8);
    expect(Math.abs(result.x[0] - 1)).toBeLessThan(1e-3);
    expect(Math.abs(result.x[1] - 1)).toBeLessThan(1e-3);
  });

  it('should minimize the Booth function', async () => {
    const result = await minimize(booth, [0, 0], {
      method: 'Nelder-Mead',
      options: { xatol: 1e-8, fatol: 1e-8 },
    });

    expect(result.success).toBe(true);
    expect(result.fun).toBeLessThan(1e-8);
    expect(Math.abs(result.x[0] - 1)).toBeLessThan(1e-3);
    expect(Math.abs(result.x[1] - 3)).toBeLessThan(1e-3);
  });

  it('should minimize Himmelblau function', async () => {
    const result = await minimize(himmelblau, [0, 0], {
      method: 'Nelder-Mead',
      options: { xatol: 1e-8, fatol: 1e-8 },
    });

    expect(result.success).toBe(true);
    expect(result.fun).toBeLessThan(1e-6);
  });

  it('should return final_simplex', async () => {
    const result = await minimize(quadratic, [1, 2], {
      method: 'Nelder-Mead',
    });

    expect(result.final_simplex).toBeDefined();
    expect(result.final_simplex!.vertices.length).toBe(3); // n+1 = 3
    expect(result.final_simplex!.values.length).toBe(3);
    // Values should be sorted (best first)
    expect(result.final_simplex!.values[0]).toBeLessThanOrEqual(
      result.final_simplex!.values[1]
    );
  });

  it('should count function evaluations', async () => {
    const result = await minimize(quadratic, [1, 2], {
      method: 'Nelder-Mead',
    });

    expect(result.nfev).toBeGreaterThan(0);
    expect(result.nit).toBeGreaterThan(0);
  });

  it('should respect maxfev limit', async () => {
    const result = await minimize(rosenbrock, [-1, -1], {
      method: 'Nelder-Mead',
      options: { maxfev: 10 },
    });

    expect(result.status).toBe(1); // maxfev exceeded
    expect(result.success).toBe(false);
    expect(result.nfev).toBeLessThanOrEqual(10 + 5); // some slack for initial simplex eval
  });

  it('should respect maxiter limit', async () => {
    const result = await minimize(rosenbrock, [-1, -1], {
      method: 'Nelder-Mead',
      options: { maxiter: 5 },
    });

    expect(result.status).toBe(2); // maxiter exceeded
    expect(result.success).toBe(false);
    expect(result.nit).toBeLessThanOrEqual(6);
  });

  it('should support adaptive parameters', async () => {
    const result = await minimize(quadratic, [1, 2, 3, 4, 5], {
      method: 'Nelder-Mead',
      options: { adaptive: true, xatol: 1e-8, fatol: 1e-8 },
    });

    expect(result.success).toBe(true);
    expect(result.fun).toBeLessThan(1e-8);
  });

  it('should support bounds', async () => {
    // Minimize quadratic with bounds that exclude the true minimum
    const result = await minimize(quadratic, [5, 5], {
      method: 'Nelder-Mead',
      bounds: {
        lb: [2, 2],
        ub: [10, 10],
      },
      options: { xatol: 1e-8, fatol: 1e-8 },
    });

    expect(result.success).toBe(true);
    // Solution should be at the lower bound
    expect(Math.abs(result.x[0] - 2)).toBeLessThan(1e-3);
    expect(Math.abs(result.x[1] - 2)).toBeLessThan(1e-3);
  });

  it('should support initial_simplex', async () => {
    // Use an asymmetric simplex that surrounds the origin
    const result = await minimize(quadratic, [0, 0], {
      method: 'Nelder-Mead',
      options: {
        initial_simplex: [
          [1, 0],
          [0, 1],
          [-1, -1],
        ],
        xatol: 1e-8,
        fatol: 1e-8,
        maxiter: 1000,
      },
    });

    expect(result.success).toBe(true);
    expect(result.fun).toBeLessThan(1e-8);
  });

  it('should handle 1D optimization', async () => {
    const f = (x: number[]) => (x[0] - 3) ** 2;
    const result = await minimize(f, [0], {
      method: 'Nelder-Mead',
      options: { xatol: 1e-8, fatol: 1e-8 },
    });

    expect(result.success).toBe(true);
    expect(Math.abs(result.x[0] - 3)).toBeLessThan(1e-3);
  });

  it('should work via dispatcher with default tol', async () => {
    const result = await minimize(quadratic, [5, -3], {
      method: 'Nelder-Mead',
    });

    expect(result.success).toBe(true);
    expect(result.fun).toBeLessThan(1e-4);
  });
});

describe('minimizeNelderMead direct', () => {
  it('should work when called directly', async () => {
    const result = await minimizeNelderMead(quadratic, [1, 2]);

    expect(result.success).toBe(true);
    expect(result.fun).toBeLessThan(1e-4);
  });
});
