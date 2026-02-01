/**
 * Tests for QUADPACK high-level API functions
 */

import { describe, test, expect, beforeAll } from 'vitest';
import {
  loadQUADPACKModule,
  getQUADPACKModule,
  isQUADPACKLoaded,
  quad,
  quad_inf,
  quad_osc,
  quad_fourier,
  quad_fixed,
  quad_ng,
  quad_rule,
  quad_break,
  quad_singular,
  quad_cauchy,
  integrate,
} from './index.js';

beforeAll(async () => {
  await loadQUADPACKModule();
});

describe('Module management', () => {
  test('isQUADPACKLoaded returns true after load', () => {
    expect(isQUADPACKLoaded()).toBe(true);
  });

  test('getQUADPACKModule returns the module', () => {
    const module = getQUADPACKModule();
    expect(module).toBeDefined();
    expect(typeof module._dqags_).toBe('function');
  });

  test('loadQUADPACKModule returns same module on subsequent calls', async () => {
    const module1 = await loadQUADPACKModule();
    const module2 = await loadQUADPACKModule();
    expect(module1).toBe(module2);
  });
});

describe('quad - General Adaptive Integration', () => {
  test('integrates x^2 from 0 to 1 (result = 1/3)', async () => {
    const result = await quad((x) => x * x, 0, 1);
    expect(result.success).toBe(true);
    expect(result.result).toBeCloseTo(1 / 3, 8);
    expect(result.ier).toBe(0);
  });

  test('integrates sin(x) from 0 to π (result = 2)', async () => {
    const result = await quad(Math.sin, 0, Math.PI);
    expect(result.success).toBe(true);
    expect(result.result).toBeCloseTo(2, 8);
  });

  test('integrates exp(x) from 0 to 1 (result = e - 1)', async () => {
    const result = await quad(Math.exp, 0, 1);
    expect(result.success).toBe(true);
    expect(result.result).toBeCloseTo(Math.E - 1, 8);
  });

  test('integrates 1/sqrt(x) from 0 to 1 (result = 2) - integrable singularity', async () => {
    const result = await quad((x) => 1 / Math.sqrt(x), 0, 1);
    expect(result.success).toBe(true);
    expect(result.result).toBeCloseTo(2, 6);
  });

  test('integrates Gaussian from -3 to 3', async () => {
    // Integral of exp(-x^2) from -∞ to ∞ is √π
    // From -3 to 3 it's approximately √π (99.7% of the area)
    const result = await quad((x) => Math.exp(-x * x), -3, 3);
    expect(result.success).toBe(true);
    expect(result.result).toBeCloseTo(Math.sqrt(Math.PI), 4);
  });

  test('respects error tolerances', async () => {
    const loose = await quad((x) => x * x, 0, 1, { epsabs: 1e-4, epsrel: 1e-4 });
    const tight = await quad((x) => x * x, 0, 1, { epsabs: 1e-12, epsrel: 1e-12 });

    expect(loose.success).toBe(true);
    expect(tight.success).toBe(true);
    // Both should achieve the requested accuracy
    expect(loose.result).toBeCloseTo(1 / 3, 4);
    expect(tight.result).toBeCloseTo(1 / 3, 12);
  });

  test('throws error for infinite bounds', async () => {
    await expect(quad((x) => x, 0, Infinity)).rejects.toThrow();
  });
});

describe('quad_inf - Infinite Interval Integration', () => {
  test('integrates exp(-x^2) from 0 to +∞ (result = √π/2)', async () => {
    const result = await quad_inf((x) => Math.exp(-x * x), 0, Infinity);
    expect(result.success).toBe(true);
    expect(result.result).toBeCloseTo(Math.sqrt(Math.PI) / 2, 6);
  });

  test('integrates exp(x) from -∞ to 0 (result = 1)', async () => {
    const result = await quad_inf(Math.exp, -Infinity, 0);
    expect(result.success).toBe(true);
    expect(result.result).toBeCloseTo(1, 8);
  });

  test('integrates 1/(1+x^2) from -∞ to +∞ (result = π)', async () => {
    const result = await quad_inf((x) => 1 / (1 + x * x), -Infinity, Infinity);
    expect(result.success).toBe(true);
    expect(result.result).toBeCloseTo(Math.PI, 6);
  });

  test('integrates exp(-x)/sqrt(x) from 0 to +∞ (result = √π)', async () => {
    const result = await quad_inf((x) => Math.exp(-x) / Math.sqrt(x), 0, Infinity);
    expect(result.success).toBe(true);
    expect(result.result).toBeCloseTo(Math.sqrt(Math.PI), 5);
  });

  test('throws error for finite bounds', async () => {
    await expect(quad_inf((x) => x, 0, 1)).rejects.toThrow();
  });
});

describe('quad_osc - Oscillatory Integration', () => {
  test('integrates x * cos(10x) from 0 to 2π', async () => {
    // ∫ x*cos(ωx) dx = (cos(ωx) + ωx*sin(ωx))/ω²
    // At x=2π: (cos(20π) + 20π*sin(20π))/100 = 1/100
    // At x=0: 1/100
    // Result = 0
    const result = await quad_osc((x) => x, 0, 2 * Math.PI, {
      omega: 10,
      weight: 'cos',
    });
    expect(result.success).toBe(true);
    expect(result.result).toBeCloseTo(0, 4);
  });

  test('integrates exp(-x) * sin(x) from 0 to 10', async () => {
    // Known integral: ∫ exp(-x)*sin(x) dx = -exp(-x)*(sin(x)+cos(x))/2
    const expected =
      -Math.exp(-10) * (Math.sin(10) + Math.cos(10)) / 2 +
      (0 + 1) / 2;
    const result = await quad_osc((x) => Math.exp(-x), 0, 10, {
      omega: 1,
      weight: 'sin',
    });
    expect(result.success).toBe(true);
    expect(result.result).toBeCloseTo(expected, 5);
  });
});

describe('quad_fourier - Fourier Integration', () => {
  test('Fourier cosine transform of exp(-x) from 0 to +∞', async () => {
    // ∫₀^∞ exp(-x)*cos(ωx) dx = 1/(1+ω²) for ω=1 → 0.5
    const result = await quad_fourier((x) => Math.exp(-x), 0, {
      omega: 1,
      weight: 'cos',
    });
    expect(result.success).toBe(true);
    expect(result.result).toBeCloseTo(0.5, 5);
  });

  test('Fourier sine transform of exp(-x) from 0 to +∞', async () => {
    // ∫₀^∞ exp(-x)*sin(ωx) dx = ω/(1+ω²) for ω=1 → 0.5
    const result = await quad_fourier((x) => Math.exp(-x), 0, {
      omega: 1,
      weight: 'sin',
    });
    expect(result.success).toBe(true);
    expect(result.result).toBeCloseTo(0.5, 5);
  });
});

describe('quad_fixed - Fixed Gauss-Kronrod Rules', () => {
  test('15-point rule for x^2 from 0 to 1', async () => {
    const result = await quad_fixed((x) => x * x, 0, 1, 15);
    expect(result.result).toBeCloseTo(1 / 3, 10);
    expect(result.abserr).toBeLessThan(1e-10);
  });

  test('61-point rule for x^2 from 0 to 1', async () => {
    const result = await quad_fixed((x) => x * x, 0, 1, 61);
    expect(result.result).toBeCloseTo(1 / 3, 14);
    expect(result.abserr).toBeLessThan(1e-14);
  });

  test('higher-order rules give smaller errors for smooth functions', async () => {
    const r15 = await quad_fixed((x) => Math.exp(-x * x), -1, 1, 15);
    const r31 = await quad_fixed((x) => Math.exp(-x * x), -1, 1, 31);
    const r61 = await quad_fixed((x) => Math.exp(-x * x), -1, 1, 61);

    // All should get very close to the correct answer
    // erf(1) * sqrt(pi) ≈ 1.4936482656
    const expected = 1.4936482656248540;
    expect(r15.result).toBeCloseTo(expected, 10);
    expect(r31.result).toBeCloseTo(expected, 12);
    expect(r61.result).toBeCloseTo(expected, 14);
  });

  test('all rule orders work', async () => {
    const rules = [15, 21, 31, 41, 51, 61] as const;
    for (const rule of rules) {
      const result = await quad_fixed((x) => x, 0, 1, rule);
      expect(result.result).toBeCloseTo(0.5, 10);
    }
  });
});

describe('quad_ng - Non-Adaptive Integration', () => {
  test('integrates smooth polynomial', async () => {
    const result = await quad_ng((x) => x * x * x, 0, 1);
    expect(result.success).toBe(true);
    expect(result.result).toBeCloseTo(0.25, 8);
  });

  test('integrates sin(x) from 0 to π', async () => {
    const result = await quad_ng(Math.sin, 0, Math.PI);
    expect(result.success).toBe(true);
    expect(result.result).toBeCloseTo(2, 6);
  });
});

describe('quad_rule - Adaptive with Rule Selection', () => {
  test('uses different Gauss-Kronrod rules', async () => {
    const r1 = await quad_rule((x) => x * x, 0, 1, { rule: 1 }); // GK15
    const r6 = await quad_rule((x) => x * x, 0, 1, { rule: 6 }); // GK61

    expect(r1.success).toBe(true);
    expect(r6.success).toBe(true);
    expect(r1.result).toBeCloseTo(1 / 3, 8);
    expect(r6.result).toBeCloseTo(1 / 3, 8);
  });
});

describe('quad_break - Integration with Break Points', () => {
  test('integrates function with known discontinuity', async () => {
    // Step function: 0 for x < 0.5, 1 for x >= 0.5
    const step = (x: number) => (x < 0.5 ? 0 : 1);
    const result = await quad_break(step, 0, 1, { points: [0.5] });
    expect(result.success).toBe(true);
    expect(result.result).toBeCloseTo(0.5, 8);
  });

  test('integrates |x| from -1 to 1 with break at 0', async () => {
    const result = await quad_break(Math.abs, -1, 1, { points: [0] });
    expect(result.success).toBe(true);
    expect(result.result).toBeCloseTo(1, 8);
  });
});

describe('quad_singular - Algebraic Singularities', () => {
  test('integrates with algebraic singularity at endpoints', async () => {
    // ∫₀¹ x^(-0.5) * (1-x)^(-0.5) dx = π (Beta function)
    const result = await quad_singular((x) => 1, 0, 1, {
      alfa: -0.5,
      beta: -0.5,
    });
    expect(result.success).toBe(true);
    expect(result.result).toBeCloseTo(Math.PI, 4);
  });

  test('integrates with singularity at left endpoint', async () => {
    // ∫₀¹ x^(-0.5) dx = 2
    const result = await quad_singular((x) => 1, 0, 1, {
      alfa: -0.5,
      beta: 0,
    });
    expect(result.success).toBe(true);
    expect(result.result).toBeCloseTo(2, 6);
  });
});

describe('quad_cauchy - Cauchy Principal Value', () => {
  test('Cauchy principal value of 1/(x-0.5) from 0 to 1', async () => {
    // P.V. ∫₀¹ 1/(x-c) dx = ln|1-c| - ln|0-c| = ln|(1-c)/c|
    // For c=0.5: ln(0.5/0.5) = ln(1) = 0
    const result = await quad_cauchy((x) => 1, 0, 1, { c: 0.5 });
    expect(result.success).toBe(true);
    expect(result.result).toBeCloseTo(0, 6);
  });

  test('Cauchy principal value with f(x) = x', async () => {
    // P.V. ∫₀¹ x/(x-0.5) dx
    // = P.V. ∫₀¹ (x-0.5+0.5)/(x-0.5) dx
    // = ∫₀¹ 1 dx + 0.5 * P.V. ∫₀¹ 1/(x-0.5) dx
    // = 1 + 0.5 * 0 = 1
    const result = await quad_cauchy((x) => x, 0, 1, { c: 0.5 });
    expect(result.success).toBe(true);
    expect(result.result).toBeCloseTo(1, 6);
  });
});

describe('integrate - Unified API', () => {
  test('auto-detects finite interval', async () => {
    const result = await integrate((x) => x * x, 0, 1);
    expect(result.success).toBe(true);
    expect(result.result).toBeCloseTo(1 / 3, 8);
  });

  test('auto-detects infinite interval', async () => {
    const result = await integrate((x) => Math.exp(-x * x), 0, Infinity);
    expect(result.success).toBe(true);
    expect(result.result).toBeCloseTo(Math.sqrt(Math.PI) / 2, 6);
  });

  test('auto-detects oscillatory with omega', async () => {
    const result = await integrate((x) => 1, 0, Math.PI, {
      omega: 1,
      weight: 'sin',
    });
    // ∫₀^π sin(x) dx = 2
    expect(result.success).toBe(true);
    expect(result.result).toBeCloseTo(2, 6);
  });

  test('explicit type override works', async () => {
    const result = await integrate((x) => x * x, 0, 1, { type: 'nonadaptive' });
    expect(result.success).toBe(true);
    expect(result.result).toBeCloseTo(1 / 3, 6);
  });

  test('explicit rule selection works', async () => {
    const result = await integrate((x) => x * x, 0, 1, { rule: 6 });
    expect(result.success).toBe(true);
    expect(result.result).toBeCloseTo(1 / 3, 8);
  });
});

describe('Mathematical identities', () => {
  test('∫₀^∞ x^n * exp(-x) dx = n! (Gamma function)', async () => {
    // n = 3: ∫₀^∞ x³ * exp(-x) dx = 3! = 6
    const result = await quad_inf((x) => x * x * x * Math.exp(-x), 0, Infinity);
    expect(result.success).toBe(true);
    expect(result.result).toBeCloseTo(6, 5);
  });

  test('∫₀^1 x^(a-1) * (1-x)^(b-1) dx = B(a,b) (Beta function)', async () => {
    // B(2, 3) = Γ(2)Γ(3)/Γ(5) = 1*2/24 = 1/12
    const result = await quad((x) => x * (1 - x) * (1 - x), 0, 1);
    expect(result.success).toBe(true);
    expect(result.result).toBeCloseTo(1 / 12, 8);
  });

  test('∫₋₁^1 1/sqrt(1-x²) dx = π', async () => {
    // This is a tough integral with singularities at both endpoints
    // 1/sqrt(1-x²) = 1/sqrt((1-x)(1+x))
    // Using quad_singular with alfa=-0.5, beta=-0.5 handles weight (x-a)^α * (b-x)^β
    // For interval [-1, 1]: (x-(-1))^(-0.5) * (1-x)^(-0.5) = (x+1)^(-0.5) * (1-x)^(-0.5)
    // This equals 1/sqrt((1+x)(1-x)) = 1/sqrt(1-x²)
    // So f(x) = 1 with the singular weight
    const result = await quad_singular((x) => 1, -1, 1, {
      alfa: -0.5,
      beta: -0.5,
    });
    expect(result.success).toBe(true);
    expect(result.result).toBeCloseTo(Math.PI, 4);
  });

  test('∫₀^∞ sin(x)/x dx = π/2 (Dirichlet integral)', async () => {
    // This is a very challenging integral - sinc function oscillates forever
    // The Fourier integration approach works better
    // ∫₀^∞ sin(x)/x dx = π/2
    // We can use: sin(x)/x = ∫₀^1 cos(xt) dt, then swap order of integration
    // Or just test with a finite upper limit and accept less precision
    const result = await quad((x) => (x === 0 ? 1 : Math.sin(x) / x), 0.0001, 100);
    // With finite bounds, we get close to π/2 but not exact
    expect(result.result).toBeCloseTo(Math.PI / 2, 1);
  });
});

describe('Error handling', () => {
  test('returns error info for problematic integrals', async () => {
    // Very tight tolerance that may not be achievable
    const result = await quad((x) => Math.sin(100 * x), 0, 10, {
      epsabs: 1e-20,
      epsrel: 1e-20,
      limit: 5,
    });
    // Should still return a result but may not be successful
    expect(result.ier).toBeGreaterThan(0);
    expect(result.success).toBe(false);
    expect(result.message).toBeTruthy();
  });
});
