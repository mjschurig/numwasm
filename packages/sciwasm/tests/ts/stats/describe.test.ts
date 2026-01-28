import { describe, it, expect } from 'vitest';
import { NDArray } from 'numwasm';
import { stats } from '../../../dist/sciwasm.mjs';

describe('describe', () => {
  it('should compute descriptive statistics for a scalar', async () => {
    const a = await NDArray.fromArray([4.0]);
    const result = await stats.describe(a);

    expect(result.nobs).toBe(1);
    expect(result.minmax[0]).toBeCloseTo(4.0);
    expect(result.minmax[1]).toBeCloseTo(4.0);
    expect(result.mean).toBeCloseTo(4.0);
    // Variance with ddof=1 for single element is NaN
    expect(result.variance).toBeNaN();
    expect(result.skewness).toBeNaN();
    expect(result.kurtosis).toBeNaN();

    a.dispose();
  });

  it('should compute statistics for 2D array along axis=0', async () => {
    // 3 rows of ones, 2 rows of twos, 4 columns
    const ones = [
      [1, 1, 1, 1],
      [1, 1, 1, 1],
      [1, 1, 1, 1],
    ];
    const twos = [
      [2, 2, 2, 2],
      [2, 2, 2, 2],
    ];
    const a = await NDArray.fromArray([...ones, ...twos]);

    const result = await stats.describe(a);

    expect(result.nobs).toBe(5);

    // minmax
    const minArr = result.minmax[0] as NDArray;
    const maxArr = result.minmax[1] as NDArray;
    for (let i = 0; i < 4; i++) {
      expect(minArr.getFlat(i)).toBeCloseTo(1.0);
      expect(maxArr.getFlat(i)).toBeCloseTo(2.0);
    }

    // mean = 1.4
    const meanArr = result.mean as NDArray;
    for (let i = 0; i < 4; i++) {
      expect(meanArr.getFlat(i)).toBeCloseTo(1.4, 10);
    }

    // variance with ddof=1 = 0.3
    const varArr = result.variance as NDArray;
    for (let i = 0; i < 4; i++) {
      expect(varArr.getFlat(i)).toBeCloseTo(0.3, 10);
    }

    // skewness ≈ 0.4082
    const skArr = result.skewness as NDArray;
    for (let i = 0; i < 4; i++) {
      expect(skArr.getFlat(i)).toBeCloseTo(0.40824829046386357, 6);
    }

    // kurtosis ≈ -1.8333
    const kurtArr = result.kurtosis as NDArray;
    for (let i = 0; i < 4; i++) {
      expect(kurtArr.getFlat(i)).toBeCloseTo(-1.833333333333333, 6);
    }

    a.dispose();
  });

  it('should compute statistics along axis=1', async () => {
    const ones = [
      [1, 1, 1, 1],
      [1, 1, 1, 1],
      [1, 1, 1, 1],
    ];
    const twos = [
      [2, 2, 2, 2],
      [2, 2, 2, 2],
    ];
    const a = await NDArray.fromArray([...ones, ...twos]);
    // Transpose: 4 rows x 5 cols, then axis=1 => same result
    const aT = a.T;

    const result = await stats.describe(aT, { axis: 1 });

    expect(result.nobs).toBe(5);

    const meanArr = result.mean as NDArray;
    for (let i = 0; i < 4; i++) {
      expect(meanArr.getFlat(i)).toBeCloseTo(1.4, 10);
    }

    a.dispose();
  });

  it('should handle ddof=0 (population variance)', async () => {
    const ones = [
      [1, 1, 1, 1],
      [1, 1, 1, 1],
      [1, 1, 1, 1],
    ];
    const twos = [
      [2, 2, 2, 2],
      [2, 2, 2, 2],
    ];
    const a = await NDArray.fromArray([...ones, ...twos]);

    const result = await stats.describe(a, { ddof: 0 });

    expect(result.nobs).toBe(5);

    // Population variance = 0.24 (vs sample variance 0.3)
    const varArr = result.variance as NDArray;
    for (let i = 0; i < 4; i++) {
      expect(varArr.getFlat(i)).toBeCloseTo(0.24, 10);
    }

    a.dispose();
  });

  it('should compute over all elements when axis=null', async () => {
    const ones = [
      [1, 1, 1, 1],
      [1, 1, 1, 1],
      [1, 1, 1, 1],
    ];
    const twos = [
      [2, 2, 2, 2],
      [2, 2, 2, 2],
    ];
    const a = await NDArray.fromArray([...ones, ...twos]);

    const result = await stats.describe(a, { axis: null });

    expect(result.nobs).toBe(20);
    expect(result.minmax[0] as number).toBeCloseTo(1.0);
    expect(result.minmax[1] as number).toBeCloseTo(2.0);
    expect(result.mean as number).toBeCloseTo(1.4, 10);
    expect(result.variance as number).toBeCloseTo(0.25263157894736848, 10);
    expect(result.skewness as number).toBeCloseTo(0.4082482904638634, 6);
    expect(result.kurtosis as number).toBeCloseTo(-1.8333333333333333, 6);

    a.dispose();
  });

  it('should throw on empty input', async () => {
    const a = await NDArray.fromArray([]);
    await expect(stats.describe(a)).rejects.toThrow('The input must not be empty.');
    a.dispose();
  });

  it('should propagate NaN by default', async () => {
    const data = [0, 1, 2, 3, 4, 5, 6, 7, 8, NaN];
    const a = await NDArray.fromArray(data);
    const result = await stats.describe(a);

    expect(result.nobs).toBe(10);
    // With NaN propagation, min, max, mean should be NaN
    expect(result.minmax[0] as number).toBeNaN();
    expect(result.minmax[1] as number).toBeNaN();
    expect(result.mean as number).toBeNaN();
    expect(result.variance as number).toBeNaN();
    expect(result.skewness as number).toBeNaN();
    expect(result.kurtosis as number).toBeNaN();

    a.dispose();
  });

  it('should raise on NaN when nanPolicy=raise', async () => {
    const data = [0, 1, 2, 3, 4, 5, 6, 7, 8, NaN];
    const a = await NDArray.fromArray(data);
    await expect(
      stats.describe(a, { nanPolicy: 'raise' }),
    ).rejects.toThrow('The input contains nan values');
    a.dispose();
  });

  it('should throw for invalid nanPolicy', async () => {
    const a = await NDArray.fromArray([1, 2, 3]);
    await expect(
      stats.describe(a, { nanPolicy: 'foobar' as any }),
    ).rejects.toThrow("nan_policy must be one of");
    a.dispose();
  });

  it('should have all result attributes', async () => {
    const a = await NDArray.fromArray([0, 1, 2, 3, 4]);
    const result = await stats.describe(a);

    expect(result).toHaveProperty('nobs');
    expect(result).toHaveProperty('minmax');
    expect(result).toHaveProperty('mean');
    expect(result).toHaveProperty('variance');
    expect(result).toHaveProperty('skewness');
    expect(result).toHaveProperty('kurtosis');
    expect(Array.isArray(result.minmax)).toBe(true);
    expect(result.minmax.length).toBe(2);

    a.dispose();
  });
});
