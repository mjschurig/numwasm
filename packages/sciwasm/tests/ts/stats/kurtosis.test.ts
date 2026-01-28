import { describe, it, expect } from 'vitest';
import { NDArray } from 'numwasm';
import { stats } from '../../../dist/sciwasm.mjs';

describe('kurtosis', () => {
  it('should compute Fisher kurtosis for known data', async () => {
    // 3 ones and 2 twos along axis 0
    const a = await NDArray.fromArray([
      [1, 1, 1, 1],
      [1, 1, 1, 1],
      [1, 1, 1, 1],
      [2, 2, 2, 2],
      [2, 2, 2, 2],
    ]);
    const result = await stats.kurtosis(a) as NDArray;
    // scipy gives -1.833333... for this distribution
    for (let i = 0; i < 4; i++) {
      expect(result.getFlat(i)).toBeCloseTo(-1.833333333333333, 10);
    }
    result.dispose();
    a.dispose();
  });

  it('should return NaN for constant data', async () => {
    const a = await NDArray.fromArray([[3, 3], [3, 3], [3, 3]]);
    const result = await stats.kurtosis(a) as NDArray;
    expect(result.getFlat(0)).toBeNaN();
    expect(result.getFlat(1)).toBeNaN();
    result.dispose();
    a.dispose();
  });

  it('should compute Pearson kurtosis when fisher=false', async () => {
    const a = await NDArray.fromArray([1, 2, 3, 4, 5]);
    const fisher = await stats.kurtosis(a, 0, true) as number;
    const pearson = await stats.kurtosis(a, 0, false) as number;
    expect(pearson).toBeCloseTo(fisher + 3.0, 10);
    a.dispose();
  });

  it('should compute kurtosis for a range', async () => {
    // np.arange(10) has known kurtosis
    const a = await NDArray.fromArray([0, 1, 2, 3, 4, 5, 6, 7, 8, 9]);
    const result = await stats.kurtosis(a) as number;
    // scipy.stats.kurtosis(np.arange(10)) ≈ -1.2242424...
    expect(result).toBeCloseTo(-1.2242424242424244, 10);
    a.dispose();
  });
});
