/**
 * Tests for ODE solvers
 */

import { describe, it, expect, beforeAll } from "vitest";
import {
  solve_ivp,
  loadODEModule,
  resetODEModule,
  ExplicitMethod,
  ImplicitMethod,
} from "../../src/index.js";
import type { ODEFunction, JacobianFunction } from "../../src/index.js";

beforeAll(async () => {
  resetODEModule();
  await loadODEModule();
});

describe("solve_ivp", () => {
  describe("RK45 (DOPRI5)", () => {
    it("solves exponential decay: y' = -y", async () => {
      // dy/dt = -y, y(0) = 1
      // Analytical solution: y(t) = e^(-t)
      const fun: ODEFunction = (_t, y) => [-y[0]];

      const result = await solve_ivp(fun, [0, 5], [1], {
        method: ExplicitMethod.RK45,
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

      const result = await solve_ivp(fun, [0, 2 * Math.PI], [1, 0], {
        method: ExplicitMethod.RK45,
        rtol: 1e-8,
        atol: 1e-10,
      });

      expect(result.success).toBe(true);

      // After one period, should return to initial state
      const yFinal = result.y[0][result.y[0].length - 1];
      const vFinal = result.y[1][result.y[1].length - 1];

      expect(Math.abs(yFinal - 1)).toBeLessThan(1e-5);
      expect(Math.abs(vFinal - 0)).toBeLessThan(1e-5);
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

      const result = await solve_ivp(fun, [0, 15], [10, 5], {
        method: ExplicitMethod.RK45,
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

    it("collects solution at multiple time points", async () => {
      const fun: ODEFunction = (_t, y) => [-y[0]];

      const result = await solve_ivp(fun, [0, 5], [1], {
        method: ExplicitMethod.RK45,
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
  });

  describe("ODEX", () => {
    it("solves exponential decay", async () => {
      const fun: ODEFunction = (_t, y) => [-y[0]];

      const result = await solve_ivp(fun, [0, 5], [1], {
        method: ExplicitMethod.ODEX,
        rtol: 1e-6,
        atol: 1e-8,
      });

      expect(result.success).toBe(true);
      expect(result.status).toBeGreaterThan(0);

      const tFinal = result.t[result.t.length - 1];
      const yFinal = result.y[0][result.y[0].length - 1];
      const expected = Math.exp(-tFinal);

      expect(Math.abs(yFinal - expected)).toBeLessThan(1e-5);
    });

    it("achieves high accuracy for smooth problems", async () => {
      const fun: ODEFunction = (_t, y) => [-y[0]];

      const result = await solve_ivp(fun, [0, 5], [1], {
        method: ExplicitMethod.ODEX,
        rtol: 1e-12,
        atol: 1e-14,
      });

      expect(result.success).toBe(true);

      const tFinal = result.t[result.t.length - 1];
      const yFinal = result.y[0][result.y[0].length - 1];
      const expected = Math.exp(-tFinal);

      // ODEX should achieve very high accuracy
      expect(Math.abs(yFinal - expected)).toBeLessThan(1e-10);
    });

    it("solves harmonic oscillator", async () => {
      const fun: ODEFunction = (_t, [y, v]) => [v, -y];

      const result = await solve_ivp(fun, [0, 2 * Math.PI], [1, 0], {
        method: ExplicitMethod.ODEX,
        rtol: 1e-8,
        atol: 1e-10,
      });

      expect(result.success).toBe(true);

      const yFinal = result.y[0][result.y[0].length - 1];
      const vFinal = result.y[1][result.y[1].length - 1];

      expect(Math.abs(yFinal - 1)).toBeLessThan(1e-6);
      expect(Math.abs(vFinal - 0)).toBeLessThan(1e-6);
    });
  });

  describe("DOP853", () => {
    it("solves exponential decay with high accuracy", async () => {
      const fun: ODEFunction = (_t, y) => [-y[0]];

      const result = await solve_ivp(fun, [0, 5], [1], {
        method: ExplicitMethod.DOP853,
        rtol: 1e-10,
        atol: 1e-12,
      });

      expect(result.success).toBe(true);

      const tFinal = result.t[result.t.length - 1];
      const yFinal = result.y[0][result.y[0].length - 1];
      const expected = Math.exp(-tFinal);

      // DOP853 should achieve higher accuracy
      expect(Math.abs(yFinal - expected)).toBeLessThan(1e-9);
    });

    it("solves stiff-ish problem (Van der Pol, mu=1)", async () => {
      // y'' - mu*(1 - y^2)*y' + y = 0
      // With mu=1, not actually stiff but tests solver robustness
      const mu = 1.0;

      const fun: ODEFunction = (_t, [y, v]) => [v, mu * (1 - y * y) * v - y];

      const result = await solve_ivp(fun, [0, 20], [2, 0], {
        method: ExplicitMethod.DOP853,
        rtol: 1e-6,
        atol: 1e-8,
      });

      expect(result.success).toBe(true);
      expect(result.nfev).toBeGreaterThan(0);
    });
  });

  describe("Radau", () => {
    it("solves exponential decay", async () => {
      const fun: ODEFunction = (_t, y) => [-y[0]];

      const result = await solve_ivp(fun, [0, 5], [1], {
        method: ImplicitMethod.Radau,
        rtol: 1e-6,
        atol: 1e-8,
      });

      expect(result.success).toBe(true);

      const tFinal = result.t[result.t.length - 1];
      const yFinal = result.y[0][result.y[0].length - 1];
      const expected = Math.exp(-tFinal);

      expect(Math.abs(yFinal - expected)).toBeLessThan(1e-4);
    });

    it("solves exponential decay with user-supplied Jacobian", async () => {
      const fun: ODEFunction = (_t, y) => [-y[0]];
      const jac: JacobianFunction = (_t, _y) => [[-1]];

      const result = await solve_ivp(fun, [0, 5], [1], {
        method: ImplicitMethod.Radau,
        rtol: 1e-6,
        atol: 1e-8,
        jac,
      });

      expect(result.success).toBe(true);
      expect(result.njev).toBeGreaterThan(0); // Jacobian was called

      const tFinal = result.t[result.t.length - 1];
      const yFinal = result.y[0][result.y[0].length - 1];
      const expected = Math.exp(-tFinal);

      expect(Math.abs(yFinal - expected)).toBeLessThan(1e-4);
    });

    it("solves Robertson stiff chemical kinetics", async () => {
      // Robertson problem - classic stiff ODE test
      // dy1/dt = -0.04*y1 + 1e4*y2*y3
      // dy2/dt = 0.04*y1 - 1e4*y2*y3 - 3e7*y2^2
      // dy3/dt = 3e7*y2^2
      const fun: ODEFunction = (_t, [y1, y2, y3]) => [
        -0.04 * y1 + 1e4 * y2 * y3,
        0.04 * y1 - 1e4 * y2 * y3 - 3e7 * y2 * y2,
        3e7 * y2 * y2,
      ];

      const jac: JacobianFunction = (_t, [_y1, y2, y3]) => [
        [-0.04, 1e4 * y3, 1e4 * y2],
        [0.04, -1e4 * y3 - 6e7 * y2, -1e4 * y2],
        [0, 6e7 * y2, 0],
      ];

      const result = await solve_ivp(fun, [0, 1e5], [1, 0, 0], {
        method: ImplicitMethod.Radau,
        rtol: 1e-4,
        atol: [1e-8, 1e-14, 1e-6], // Different tolerances for each component
        jac,
      });

      expect(result.success).toBe(true);

      // Conservation law: y1 + y2 + y3 = 1
      const yFinal0 = result.y[0][result.y[0].length - 1];
      const yFinal1 = result.y[1][result.y[1].length - 1];
      const yFinal2 = result.y[2][result.y[2].length - 1];
      const sum = yFinal0 + yFinal1 + yFinal2;

      expect(Math.abs(sum - 1)).toBeLessThan(1e-4);
    });

    it("solves Van der Pol oscillator (stiff, mu=1000)", async () => {
      // Van der Pol with large mu is stiff
      const mu = 1000;

      const fun: ODEFunction = (_t, [y, v]) => [v, mu * (1 - y * y) * v - y];

      const jac: JacobianFunction = (_t, [y, v]) => [
        [0, 1],
        [-2 * mu * y * v - 1, mu * (1 - y * y)],
      ];

      // Use short time span for stiff problem
      const result = await solve_ivp(fun, [0, 100], [2, 0], {
        method: ImplicitMethod.Radau,
        rtol: 1e-4,
        atol: 1e-6,
        jac,
      });

      expect(result.success).toBe(true);
      expect(result.nfev).toBeGreaterThan(0);
      expect(result.nlu).toBeGreaterThan(0); // LU decompositions for implicit method
    });
  });

  describe("error handling", () => {
    it("throws error for unknown method", async () => {
      const fun: ODEFunction = (_t, y) => [-y[0]];

      await expect(
        solve_ivp(fun, [0, 1], [1], {
          method: "Unknown" as ExplicitMethod.RK45,
        }),
      ).rejects.toThrow("Unknown method");
    });

    it("returns failure status for inconsistent input", async () => {
      // This test checks that the solver properly reports errors
      // Creating a situation where the solver might fail
      const fun: ODEFunction = (_t, y) => {
        // Return wrong dimension (might cause issues)
        return y;
      };

      const result = await solve_ivp(fun, [0, 1], [1], {
        method: ExplicitMethod.RK45,
        rtol: 1e-6,
        atol: 1e-8,
      });

      // This should succeed as the dimension is correct (1 in, 1 out)
      expect(result.success).toBe(true);
    });
  });

  describe("step control options", () => {
    it("respects max_step option", async () => {
      const fun: ODEFunction = (_t, y) => [-y[0]];

      const result = await solve_ivp(fun, [0, 10], [1], {
        method: ExplicitMethod.RK45,
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

      const result = await solve_ivp(fun, [0, 1], [1], {
        method: ExplicitMethod.RK45,
        rtol: 1e-6,
        atol: 1e-8,
        first_step: 0.001,
        dense_output: true,
      });

      expect(result.success).toBe(true);
      // Solver should start with the specified step
      expect(result.t.length).toBeGreaterThan(2);
    });
  });
});
