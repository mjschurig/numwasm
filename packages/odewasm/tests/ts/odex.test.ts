/**
 * Tests for ODEX solver (GBS Extrapolation Method)
 *
 * ODEX is based on the explicit midpoint rule with extrapolation.
 * It features automatic order selection and is highly efficient for
 * smooth problems requiring high accuracy.
 */

import { describe, it, expect, beforeAll } from "vitest";
import {
  odex,
  solve_ivp,
  loadODEModule,
  resetODEModule,
  ExplicitMethod,
} from "../../src/index.js";
import type { ODEFunction } from "../../src/index.js";

beforeAll(async () => {
  resetODEModule();
  await loadODEModule();
});

describe("odex", () => {
  describe("basic functionality", () => {
    it("solves exponential decay: y' = -y", async () => {
      // dy/dt = -y, y(0) = 1
      // Analytical solution: y(t) = e^(-t)
      const fun: ODEFunction = (_t, y) => [-y[0]];

      const result = await odex(fun, [0, 5], [1], {
        rtol: 1e-6,
        atol: 1e-8,
      });

      expect(result.success).toBe(true);
      expect(result.status).toBeGreaterThan(0);

      // Check final value
      const tFinal = result.t[result.t.length - 1];
      const yFinal = result.y[0][result.y[0].length - 1];
      const expected = Math.exp(-tFinal);

      expect(Math.abs(yFinal - expected)).toBeLessThan(1e-5);
    });

    it("solves harmonic oscillator: y'' = -y", async () => {
      // [y, v] where v = y'
      // dy/dt = v, dv/dt = -y
      // y(0) = 1, v(0) = 0
      // Analytical: y(t) = cos(t), v(t) = -sin(t)
      const fun: ODEFunction = (_t, [y, v]) => [v, -y];

      const result = await odex(fun, [0, 2 * Math.PI], [1, 0], {
        rtol: 1e-8,
        atol: 1e-10,
      });

      expect(result.success).toBe(true);

      // After one period, should return to initial state
      const yFinal = result.y[0][result.y[0].length - 1];
      const vFinal = result.y[1][result.y[1].length - 1];

      expect(Math.abs(yFinal - 1)).toBeLessThan(1e-6);
      expect(Math.abs(vFinal - 0)).toBeLessThan(1e-6);
    });

    it("solves linear ODE: y' = y (exponential growth)", async () => {
      // dy/dt = y, y(0) = 1
      // Analytical solution: y(t) = e^t
      const fun: ODEFunction = (_t, y) => [y[0]];

      const result = await odex(fun, [0, 3], [1], {
        rtol: 1e-8,
        atol: 1e-10,
      });

      expect(result.success).toBe(true);

      const tFinal = result.t[result.t.length - 1];
      const yFinal = result.y[0][result.y[0].length - 1];
      const expected = Math.exp(tFinal);

      expect(Math.abs(yFinal - expected) / expected).toBeLessThan(1e-7);
    });
  });

  describe("high accuracy", () => {
    it("achieves high accuracy for smooth problems", async () => {
      // ODEX should excel at smooth problems with high accuracy requirements
      const fun: ODEFunction = (_t, y) => [-y[0]];

      const result = await odex(fun, [0, 5], [1], {
        rtol: 1e-12,
        atol: 1e-14,
      });

      expect(result.success).toBe(true);

      const tFinal = result.t[result.t.length - 1];
      const yFinal = result.y[0][result.y[0].length - 1];
      const expected = Math.exp(-tFinal);

      // Should achieve very high accuracy
      expect(Math.abs(yFinal - expected)).toBeLessThan(1e-10);
    });

    it("solves Kepler orbit with high precision", async () => {
      // Two-body problem in 2D (circular orbit)
      // x'' = -x/r^3, y'' = -y/r^3 where r = sqrt(x^2 + y^2)
      // Initial conditions for circular orbit: x=1, y=0, vx=0, vy=1
      const fun: ODEFunction = (_t, [x, y, vx, vy]) => {
        const r3 = Math.pow(x * x + y * y, 1.5);
        return [vx, vy, -x / r3, -y / r3];
      };

      // One full orbit
      const result = await odex(fun, [0, 2 * Math.PI], [1, 0, 0, 1], {
        rtol: 1e-10,
        atol: 1e-12,
      });

      expect(result.success).toBe(true);

      // After one period, should return to initial position
      const xFinal = result.y[0][result.y[0].length - 1];
      const yFinal = result.y[1][result.y[1].length - 1];
      const vxFinal = result.y[2][result.y[2].length - 1];
      const vyFinal = result.y[3][result.y[3].length - 1];

      expect(Math.abs(xFinal - 1)).toBeLessThan(1e-8);
      expect(Math.abs(yFinal - 0)).toBeLessThan(1e-8);
      expect(Math.abs(vxFinal - 0)).toBeLessThan(1e-8);
      expect(Math.abs(vyFinal - 1)).toBeLessThan(1e-8);
    });
  });

  describe("multi-dimensional systems", () => {
    it("solves Lorenz system", async () => {
      // Lorenz attractor (chaotic system)
      const sigma = 10,
        rho = 28,
        beta = 8 / 3;

      const fun: ODEFunction = (_t, [x, y, z]) => [
        sigma * (y - x),
        x * (rho - z) - y,
        x * y - beta * z,
      ];

      const result = await odex(fun, [0, 1], [1, 1, 1], {
        rtol: 1e-6,
        atol: 1e-8,
      });

      expect(result.success).toBe(true);
      expect(result.nfev).toBeGreaterThan(0);

      // Solution should exist and be finite
      for (let i = 0; i < 3; i++) {
        const yLast = result.y[i][result.y[i].length - 1];
        expect(isFinite(yLast)).toBe(true);
      }
    });

    it("solves Lotka-Volterra predator-prey system", async () => {
      // dx/dt = alpha*x - beta*x*y (prey)
      // dy/dt = delta*x*y - gamma*y (predator)
      const alpha = 1.5,
        beta = 1.0,
        delta = 1.0,
        gamma = 3.0;

      const fun: ODEFunction = (_t, [x, y]) => [
        alpha * x - beta * x * y,
        delta * x * y - gamma * y,
      ];

      const result = await odex(fun, [0, 15], [10, 5], {
        rtol: 1e-6,
        atol: 1e-8,
      });

      expect(result.success).toBe(true);
      expect(result.nfev).toBeGreaterThan(0);

      // Populations should remain positive
      for (let i = 0; i < result.t.length; i++) {
        expect(result.y[0][i]).toBeGreaterThan(0);
        expect(result.y[1][i]).toBeGreaterThan(0);
      }
    });

    it("solves 10-dimensional linear system", async () => {
      // y' = A*y where A is diagonal with -1, -2, ..., -10
      const n = 10;
      const fun: ODEFunction = (_t, y) => {
        const dydt = new Array(n);
        for (let i = 0; i < n; i++) {
          dydt[i] = -(i + 1) * y[i];
        }
        return dydt;
      };

      const y0 = new Array(n).fill(1);

      const result = await odex(fun, [0, 1], y0, {
        rtol: 1e-8,
        atol: 1e-10,
      });

      expect(result.success).toBe(true);

      // Check analytical solution: y_i(t) = exp(-(i+1)*t)
      const tFinal = result.t[result.t.length - 1];
      for (let i = 0; i < n; i++) {
        const yFinal = result.y[i][result.y[i].length - 1];
        const expected = Math.exp(-(i + 1) * tFinal);
        expect(Math.abs(yFinal - expected)).toBeLessThan(1e-7);
      }
    });
  });

  describe("dense output and t_eval", () => {
    it("collects solution at multiple time points with dense_output", async () => {
      const fun: ODEFunction = (_t, y) => [-y[0]];

      const result = await odex(fun, [0, 5], [1], {
        rtol: 1e-6,
        atol: 1e-8,
        dense_output: true,
      });

      expect(result.success).toBe(true);
      // Should have multiple time points collected
      expect(result.t.length).toBeGreaterThan(2);

      // Time should be monotonically increasing
      for (let i = 1; i < result.t.length; i++) {
        expect(result.t[i]).toBeGreaterThan(result.t[i - 1]);
      }
    });

    it("returns only start and end points without dense_output", async () => {
      const fun: ODEFunction = (_t, y) => [-y[0]];

      const result = await odex(fun, [0, 5], [1], {
        rtol: 1e-6,
        atol: 1e-8,
        dense_output: false,
      });

      expect(result.success).toBe(true);
      // Without dense output, should have only start and end
      expect(result.t.length).toBe(2);
      expect(result.t[0]).toBe(0);
      expect(Math.abs(result.t[1] - 5)).toBeLessThan(1e-10);
    });
  });

  describe("ODEX-specific options", () => {
    it("respects max_step option", async () => {
      const fun: ODEFunction = (_t, y) => [-y[0]];

      const result = await odex(fun, [0, 10], [1], {
        rtol: 1e-6,
        atol: 1e-8,
        max_step: 0.1,
        dense_output: true,
      });

      expect(result.success).toBe(true);

      // With max_step=0.1 over [0,10], should have at least 100 steps
      expect(result.t.length).toBeGreaterThanOrEqual(100);
    });

    it("respects first_step option", async () => {
      const fun: ODEFunction = (_t, y) => [-y[0]];

      const result = await odex(fun, [0, 1], [1], {
        rtol: 1e-6,
        atol: 1e-8,
        first_step: 0.001,
        dense_output: true,
      });

      expect(result.success).toBe(true);
      expect(result.t.length).toBeGreaterThan(2);
    });

    it("works with custom max_columns setting", async () => {
      const fun: ODEFunction = (_t, y) => [-y[0]];

      // Use higher max_columns for potentially higher order extrapolation
      const result = await odex(fun, [0, 5], [1], {
        rtol: 1e-10,
        atol: 1e-12,
        max_columns: 12,
      });

      expect(result.success).toBe(true);

      const tFinal = result.t[result.t.length - 1];
      const yFinal = result.y[0][result.y[0].length - 1];
      const expected = Math.exp(-tFinal);

      expect(Math.abs(yFinal - expected)).toBeLessThan(1e-9);
    });
  });

  describe("vector tolerances", () => {
    it("supports vector absolute tolerances", async () => {
      // Two-component system with different scales
      const fun: ODEFunction = (_t, [y1, y2]) => [-y1, -0.001 * y2];

      const result = await odex(fun, [0, 100], [1, 1000], {
        rtol: 1e-6,
        atol: [1e-8, 1e-5], // Different tolerances for each component
      });

      expect(result.success).toBe(true);

      const y1Final = result.y[0][result.y[0].length - 1];
      const y2Final = result.y[1][result.y[1].length - 1];

      const tFinal = result.t[result.t.length - 1];
      const expected1 = Math.exp(-tFinal);
      const expected2 = 1000 * Math.exp(-0.001 * tFinal);

      expect(Math.abs(y1Final - expected1)).toBeLessThan(1e-5);
      expect(Math.abs(y2Final - expected2) / expected2).toBeLessThan(1e-3);
    });
  });

  describe("backward integration", () => {
    it("integrates backwards in time", async () => {
      const fun: ODEFunction = (_t, y) => [-y[0]];

      // y(5) = e^(-5), solve backwards to t=0
      const y5 = Math.exp(-5);
      const result = await odex(fun, [5, 0], [y5], {
        rtol: 1e-8,
        atol: 1e-10,
      });

      expect(result.success).toBe(true);

      // Should recover y(0) = 1
      const yFinal = result.y[0][result.y[0].length - 1];
      expect(Math.abs(yFinal - 1)).toBeLessThan(1e-6);
    });
  });

  describe("comparison with other solvers", () => {
    it("matches DOPRI5 solution for smooth problem", async () => {
      const fun: ODEFunction = (_t, y) => [-y[0]];

      const resultOdex = await odex(fun, [0, 5], [1], {
        rtol: 1e-8,
        atol: 1e-10,
      });

      const resultDopri5 = await solve_ivp(fun, [0, 5], [1], {
        method: "RK45",
        rtol: 1e-8,
        atol: 1e-10,
      });

      expect(resultOdex.success).toBe(true);
      expect(resultDopri5.success).toBe(true);

      // Final values should match closely
      const yOdex = resultOdex.y[0][resultOdex.y[0].length - 1];
      const yDopri5 = resultDopri5.y[0][resultDopri5.y[0].length - 1];

      expect(Math.abs(yOdex - yDopri5)).toBeLessThan(1e-8);
    });
  });
});

describe("solve_ivp with ODEX method", () => {
  it("solves via unified API with string method", async () => {
    const fun: ODEFunction = (_t, y) => [-y[0]];

    const result = await solve_ivp(fun, [0, 5], [1], {
      method: "ODEX",
      rtol: 1e-6,
      atol: 1e-8,
    });

    expect(result.success).toBe(true);

    const tFinal = result.t[result.t.length - 1];
    const yFinal = result.y[0][result.y[0].length - 1];
    const expected = Math.exp(-tFinal);

    expect(Math.abs(yFinal - expected)).toBeLessThan(1e-5);
  });

  it("solves via unified API with enum method", async () => {
    const fun: ODEFunction = (_t, y) => [-y[0]];

    const result = await solve_ivp(fun, [0, 5], [1], {
      method: ExplicitMethod.ODEX,
      rtol: 1e-6,
      atol: 1e-8,
    });

    expect(result.success).toBe(true);
    expect(result.status).toBeGreaterThan(0);
  });

  it("works with dense_output via unified API", async () => {
    const fun: ODEFunction = (_t, y) => [-y[0]];

    const result = await solve_ivp(fun, [0, 5], [1], {
      method: ExplicitMethod.ODEX,
      rtol: 1e-6,
      atol: 1e-8,
      dense_output: true,
    });

    expect(result.success).toBe(true);
    expect(result.t.length).toBeGreaterThan(2);
  });
});

describe("ODEX edge cases", () => {
  it("handles zero initial condition", async () => {
    // y' = 1, y(0) = 0 => y(t) = t
    const fun: ODEFunction = (_t, _y) => [1];

    const result = await odex(fun, [0, 5], [0], {
      rtol: 1e-6,
      atol: 1e-8,
    });

    expect(result.success).toBe(true);

    const tFinal = result.t[result.t.length - 1];
    const yFinal = result.y[0][result.y[0].length - 1];

    expect(Math.abs(yFinal - tFinal)).toBeLessThan(1e-5);
  });

  it("handles time-dependent ODE", async () => {
    // y' = t, y(0) = 0 => y(t) = t^2/2
    const fun: ODEFunction = (t, _y) => [t];

    const result = await odex(fun, [0, 4], [0], {
      rtol: 1e-8,
      atol: 1e-10,
    });

    expect(result.success).toBe(true);

    const tFinal = result.t[result.t.length - 1];
    const yFinal = result.y[0][result.y[0].length - 1];
    const expected = (tFinal * tFinal) / 2;

    expect(Math.abs(yFinal - expected)).toBeLessThan(1e-7);
  });

  it("handles rapidly oscillating solution", async () => {
    // y'' + 100*y = 0 => y = cos(10t), high frequency
    const omega = 10;
    const fun: ODEFunction = (_t, [y, v]) => [v, -(omega * omega) * y];

    const result = await odex(fun, [0, 2 * Math.PI], [1, 0], {
      rtol: 1e-8,
      atol: 1e-10,
    });

    expect(result.success).toBe(true);

    // After one period (at omega=10), 2*pi corresponds to 10 oscillations
    const yFinal = result.y[0][result.y[0].length - 1];
    const expected = Math.cos(10 * 2 * Math.PI); // = 1

    expect(Math.abs(yFinal - expected)).toBeLessThan(1e-5);
  });

  it("handles negative times", async () => {
    // Solve from t=-1 to t=1
    const fun: ODEFunction = (t, _y) => [t];

    const result = await odex(fun, [-1, 1], [0.5], {
      // y(-1) = 0.5, y' = t
      // y(t) = t^2/2 + C, y(-1) = 0.5 => C = 0
      // y(1) = 0.5
      rtol: 1e-8,
      atol: 1e-10,
    });

    expect(result.success).toBe(true);

    const yFinal = result.y[0][result.y[0].length - 1];
    expect(Math.abs(yFinal - 0.5)).toBeLessThan(1e-7);
  });
});
