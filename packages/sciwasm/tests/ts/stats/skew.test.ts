import { describe, it, expect } from 'vitest';
import { NDArray } from 'numwasm';
import { stats } from '../../../dist/sciwasm.mjs';

describe('skew', () => {
  it('should return 0 for symmetric data', async () => {
    const a = await NDArray.fromArray([1, 2, 3, 4, 5]);
    const result = await stats.skew(a);
    expect(result).toBeCloseTo(0.0, 10);
    a.dispose();
  });

  it('should compute skewness for asymmetric data', async () => {
    // Data: 3 ones and 2 twos along axis 0 => positively skewed
    const a = await NDArray.fromArray([
      [1, 1, 1, 1],
      [1, 1, 1, 1],
      [1, 1, 1, 1],
      [2, 2, 2, 2],
      [2, 2, 2, 2],
    ]);
    const result = await stats.skew(a) as NDArray;
    // scipy gives ~0.4082 for this distribution
    for (let i = 0; i < 4; i++) {
      expect(result.getFlat(i)).toBeCloseTo(0.40824829046386357, 10);
    }
    result.dispose();
    a.dispose();
  });

  it('should return NaN for constant data', async () => {
    const a = await NDArray.fromArray([[5, 5], [5, 5], [5, 5]]);
    const result = await stats.skew(a) as NDArray;
    expect(result.getFlat(0)).toBeNaN();
    expect(result.getFlat(1)).toBeNaN();
    result.dispose();
    a.dispose();
  });

  it('should apply bias correction when bias=false', async () => {
    const a = await NDArray.fromArray([1, 2, 2, 3, 3, 3, 4, 4, 4, 4]);
    const biased = await stats.skew(a, 0, true) as number;
    const unbiased = await stats.skew(a, 0, false) as number;
    // Unbiased correction: sqrt(n*(n-1))/(n-2) * biased
    const n = 10;
    const expected = Math.sqrt(n * (n - 1)) / (n - 2) * biased;
    expect(unbiased).toBeCloseTo(expected, 10);
    a.dispose();
  });
});
