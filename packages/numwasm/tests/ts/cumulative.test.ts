/**
 * Cumulative Operations Unit Tests (Phase 22)
 *
 * Tests for cumsum, cumprod, nancumsum, nancumprod functions.
 */

import { describe, it, expect, beforeAll, afterEach } from 'vitest';
import {
  NDArray,
  DType,
  loadWasmModule,
  cumsum,
  cumprod,
  nancumsum,
  nancumprod,
} from 'numwasm';

describe('Cumulative Operations', () => {
  // Load WASM module before all tests
  beforeAll(async () => {
    await loadWasmModule();
  });

  // Track arrays for cleanup
  const arrays: NDArray[] = [];
  afterEach(() => {
    for (const arr of arrays) {
      if (!arr.isDisposed) {
        arr.dispose();
      }
    }
    arrays.length = 0;
  });

  const track = (arr: NDArray) => {
    arrays.push(arr);
    return arr;
  };

  describe('cumsum()', () => {
    describe('1D arrays', () => {
      it('should compute cumulative sum of simple array', async () => {
        const arr = track(await NDArray.fromArray([1, 2, 3, 4]));
        const result = track(await cumsum(arr));
        expect(result.shape).toEqual([4]);
        expect(result.toArray()).toEqual([1, 3, 6, 10]);
      });

      it('should compute cumulative sum of array starting with zero', async () => {
        const arr = track(await NDArray.fromArray([0, 1, 2, 3]));
        const result = track(await cumsum(arr));
        expect(result.toArray()).toEqual([0, 1, 3, 6]);
      });

      it('should handle negative values', async () => {
        const arr = track(await NDArray.fromArray([1, -2, 3, -4]));
        const result = track(await cumsum(arr));
        expect(result.toArray()).toEqual([1, -1, 2, -2]);
      });

      it('should handle floating point values', async () => {
        const arr = track(await NDArray.fromArray([0.5, 1.5, 2.5]));
        const result = track(await cumsum(arr));
        expect(result.toArray()[0]).toBeCloseTo(0.5, 10);
        expect(result.toArray()[1]).toBeCloseTo(2.0, 10);
        expect(result.toArray()[2]).toBeCloseTo(4.5, 10);
      });

      it('should handle single element array', async () => {
        const arr = track(await NDArray.fromArray([42]));
        const result = track(await cumsum(arr));
        expect(result.shape).toEqual([1]);
        expect(result.toArray()).toEqual([42]);
      });
    });

    describe('2D arrays', () => {
      it('should flatten and cumsum when axis=null', async () => {
        const arr = track(await NDArray.fromArray([[1, 2], [3, 4]]));
        const result = track(await cumsum(arr));
        expect(result.shape).toEqual([4]);
        expect(result.toArray()).toEqual([1, 3, 6, 10]);
      });

      it('should cumsum along axis=0 (down columns)', async () => {
        const arr = track(await NDArray.fromArray([[1, 2], [3, 4]]));
        const result = track(await cumsum(arr, 0));
        expect(result.shape).toEqual([2, 2]);
        expect(result.toArray()).toEqual([1, 2, 4, 6]);
      });

      it('should cumsum along axis=1 (across rows)', async () => {
        const arr = track(await NDArray.fromArray([[1, 2], [3, 4]]));
        const result = track(await cumsum(arr, 1));
        expect(result.shape).toEqual([2, 2]);
        expect(result.toArray()).toEqual([1, 3, 3, 7]);
      });

      it('should handle negative axis (-1 = last axis)', async () => {
        const arr = track(await NDArray.fromArray([[1, 2], [3, 4]]));
        const result = track(await cumsum(arr, -1));
        expect(result.shape).toEqual([2, 2]);
        expect(result.toArray()).toEqual([1, 3, 3, 7]);
      });

      it('should handle negative axis (-2 = first axis for 2D)', async () => {
        const arr = track(await NDArray.fromArray([[1, 2], [3, 4]]));
        const result = track(await cumsum(arr, -2));
        expect(result.shape).toEqual([2, 2]);
        expect(result.toArray()).toEqual([1, 2, 4, 6]);
      });
    });

    describe('3D arrays', () => {
      it('should cumsum along axis=0', async () => {
        const arr = track(
          await NDArray.fromArray([
            [[1, 2], [3, 4]],
            [[5, 6], [7, 8]],
          ])
        );
        const result = track(await cumsum(arr, 0));
        expect(result.shape).toEqual([2, 2, 2]);
        // First slice unchanged, second slice = first + second
        expect(result.toArray()).toEqual([1, 2, 3, 4, 6, 8, 10, 12]);
      });

      it('should cumsum along axis=1', async () => {
        const arr = track(
          await NDArray.fromArray([
            [[1, 2], [3, 4]],
            [[5, 6], [7, 8]],
          ])
        );
        const result = track(await cumsum(arr, 1));
        expect(result.shape).toEqual([2, 2, 2]);
        expect(result.toArray()).toEqual([1, 2, 4, 6, 5, 6, 12, 14]);
      });

      it('should cumsum along axis=2 (innermost)', async () => {
        const arr = track(
          await NDArray.fromArray([
            [[1, 2], [3, 4]],
            [[5, 6], [7, 8]],
          ])
        );
        const result = track(await cumsum(arr, 2));
        expect(result.shape).toEqual([2, 2, 2]);
        expect(result.toArray()).toEqual([1, 3, 3, 7, 5, 11, 7, 15]);
      });
    });

    describe('edge cases', () => {
      it('should handle larger arrays', async () => {
        const data = Array.from({ length: 100 }, (_, i) => i + 1);
        const arr = track(await NDArray.fromArray(data));
        const result = track(await cumsum(arr));
        // Sum 1 to n = n*(n+1)/2
        expect(result.toArray()[99]).toBe((100 * 101) / 2);
      });
    });
  });

  describe('cumprod()', () => {
    describe('1D arrays', () => {
      it('should compute cumulative product of simple array', async () => {
        const arr = track(await NDArray.fromArray([1, 2, 3, 4]));
        const result = track(await cumprod(arr));
        expect(result.shape).toEqual([4]);
        expect(result.toArray()).toEqual([1, 2, 6, 24]);
      });

      it('should handle array with zero', async () => {
        const arr = track(await NDArray.fromArray([1, 2, 0, 4]));
        const result = track(await cumprod(arr));
        expect(result.toArray()).toEqual([1, 2, 0, 0]);
      });

      it('should handle negative values', async () => {
        const arr = track(await NDArray.fromArray([1, -2, 3]));
        const result = track(await cumprod(arr));
        expect(result.toArray()).toEqual([1, -2, -6]);
      });

      it('should handle floating point values', async () => {
        const arr = track(await NDArray.fromArray([0.5, 2.0, 3.0]));
        const result = track(await cumprod(arr));
        expect(result.toArray()[0]).toBeCloseTo(0.5, 10);
        expect(result.toArray()[1]).toBeCloseTo(1.0, 10);
        expect(result.toArray()[2]).toBeCloseTo(3.0, 10);
      });

      it('should handle single element array', async () => {
        const arr = track(await NDArray.fromArray([5]));
        const result = track(await cumprod(arr));
        expect(result.shape).toEqual([1]);
        expect(result.toArray()).toEqual([5]);
      });
    });

    describe('2D arrays', () => {
      it('should flatten and cumprod when axis=null', async () => {
        const arr = track(await NDArray.fromArray([[1, 2], [3, 4]]));
        const result = track(await cumprod(arr));
        expect(result.shape).toEqual([4]);
        expect(result.toArray()).toEqual([1, 2, 6, 24]);
      });

      it('should cumprod along axis=0', async () => {
        const arr = track(await NDArray.fromArray([[1, 2], [3, 4]]));
        const result = track(await cumprod(arr, 0));
        expect(result.shape).toEqual([2, 2]);
        expect(result.toArray()).toEqual([1, 2, 3, 8]);
      });

      it('should cumprod along axis=1', async () => {
        const arr = track(await NDArray.fromArray([[1, 2], [3, 4]]));
        const result = track(await cumprod(arr, 1));
        expect(result.shape).toEqual([2, 2]);
        expect(result.toArray()).toEqual([1, 2, 3, 12]);
      });
    });
  });

  describe('nancumsum()', () => {
    describe('1D arrays with NaN', () => {
      it('should treat NaN as zero in cumsum', async () => {
        const arr = track(await NDArray.fromArray([1, NaN, 3, 4]));
        const result = track(await nancumsum(arr));
        expect(result.shape).toEqual([4]);
        expect(result.toArray()[0]).toBe(1);
        expect(result.toArray()[1]).toBe(1); // NaN treated as 0
        expect(result.toArray()[2]).toBe(4);
        expect(result.toArray()[3]).toBe(8);
      });

      it('should handle leading NaN', async () => {
        const arr = track(await NDArray.fromArray([NaN, 2, 3]));
        const result = track(await nancumsum(arr));
        expect(result.toArray()[0]).toBe(0); // Leading NaN -> 0
        expect(result.toArray()[1]).toBe(2);
        expect(result.toArray()[2]).toBe(5);
      });

      it('should handle trailing NaN', async () => {
        const arr = track(await NDArray.fromArray([1, 2, NaN]));
        const result = track(await nancumsum(arr));
        expect(result.toArray()[0]).toBe(1);
        expect(result.toArray()[1]).toBe(3);
        expect(result.toArray()[2]).toBe(3); // NaN doesn't change sum
      });

      it('should handle all NaN', async () => {
        const arr = track(await NDArray.fromArray([NaN, NaN, NaN]));
        const result = track(await nancumsum(arr));
        expect(result.toArray()).toEqual([0, 0, 0]);
      });

      it('should handle consecutive NaN', async () => {
        const arr = track(await NDArray.fromArray([1, NaN, NaN, 4]));
        const result = track(await nancumsum(arr));
        expect(result.toArray()[0]).toBe(1);
        expect(result.toArray()[1]).toBe(1);
        expect(result.toArray()[2]).toBe(1);
        expect(result.toArray()[3]).toBe(5);
      });
    });

    describe('2D arrays with NaN', () => {
      it('should handle NaN along axis=0', async () => {
        const arr = track(await NDArray.fromArray([[1, 2], [NaN, 4]]));
        const result = track(await nancumsum(arr, 0));
        expect(result.shape).toEqual([2, 2]);
        expect(result.toArray()).toEqual([1, 2, 1, 6]);
      });

      it('should handle NaN along axis=1', async () => {
        const arr = track(await NDArray.fromArray([[1, NaN], [3, 4]]));
        const result = track(await nancumsum(arr, 1));
        expect(result.shape).toEqual([2, 2]);
        expect(result.toArray()).toEqual([1, 1, 3, 7]);
      });
    });

    describe('arrays without NaN', () => {
      it('should behave like regular cumsum', async () => {
        const arr = track(await NDArray.fromArray([1, 2, 3, 4]));
        const result = track(await nancumsum(arr));
        expect(result.toArray()).toEqual([1, 3, 6, 10]);
      });
    });
  });

  describe('nancumprod()', () => {
    describe('1D arrays with NaN', () => {
      it('should treat NaN as one in cumprod', async () => {
        const arr = track(await NDArray.fromArray([1, NaN, 3, 4]));
        const result = track(await nancumprod(arr));
        expect(result.shape).toEqual([4]);
        expect(result.toArray()[0]).toBe(1);
        expect(result.toArray()[1]).toBe(1); // NaN treated as 1
        expect(result.toArray()[2]).toBe(3);
        expect(result.toArray()[3]).toBe(12);
      });

      it('should handle leading NaN', async () => {
        const arr = track(await NDArray.fromArray([NaN, 2, 3]));
        const result = track(await nancumprod(arr));
        expect(result.toArray()[0]).toBe(1); // Leading NaN -> 1
        expect(result.toArray()[1]).toBe(2);
        expect(result.toArray()[2]).toBe(6);
      });

      it('should handle trailing NaN', async () => {
        const arr = track(await NDArray.fromArray([2, 3, NaN]));
        const result = track(await nancumprod(arr));
        expect(result.toArray()[0]).toBe(2);
        expect(result.toArray()[1]).toBe(6);
        expect(result.toArray()[2]).toBe(6); // NaN doesn't change product
      });

      it('should handle all NaN', async () => {
        const arr = track(await NDArray.fromArray([NaN, NaN, NaN]));
        const result = track(await nancumprod(arr));
        expect(result.toArray()).toEqual([1, 1, 1]);
      });

      it('should handle consecutive NaN', async () => {
        const arr = track(await NDArray.fromArray([2, NaN, NaN, 3]));
        const result = track(await nancumprod(arr));
        expect(result.toArray()[0]).toBe(2);
        expect(result.toArray()[1]).toBe(2);
        expect(result.toArray()[2]).toBe(2);
        expect(result.toArray()[3]).toBe(6);
      });
    });

    describe('2D arrays with NaN', () => {
      it('should handle NaN along axis=0', async () => {
        const arr = track(await NDArray.fromArray([[1, 2], [NaN, 4]]));
        const result = track(await nancumprod(arr, 0));
        expect(result.shape).toEqual([2, 2]);
        expect(result.toArray()).toEqual([1, 2, 1, 8]);
      });

      it('should handle NaN along axis=1', async () => {
        const arr = track(await NDArray.fromArray([[2, NaN], [3, 4]]));
        const result = track(await nancumprod(arr, 1));
        expect(result.shape).toEqual([2, 2]);
        expect(result.toArray()).toEqual([2, 2, 3, 12]);
      });
    });

    describe('arrays without NaN', () => {
      it('should behave like regular cumprod', async () => {
        const arr = track(await NDArray.fromArray([1, 2, 3, 4]));
        const result = track(await nancumprod(arr));
        expect(result.toArray()).toEqual([1, 2, 6, 24]);
      });
    });
  });

  describe('NumPy compatibility', () => {
    it('cumsum matches NumPy examples', async () => {
      // np.cumsum([1, 2, 3, 4, 5]) = [1, 3, 6, 10, 15]
      const arr = track(await NDArray.fromArray([1, 2, 3, 4, 5]));
      const result = track(await cumsum(arr));
      expect(result.toArray()).toEqual([1, 3, 6, 10, 15]);
    });

    it('cumprod matches NumPy examples', async () => {
      // np.cumprod([1, 2, 3, 4]) = [1, 2, 6, 24]
      const arr = track(await NDArray.fromArray([1, 2, 3, 4]));
      const result = track(await cumprod(arr));
      expect(result.toArray()).toEqual([1, 2, 6, 24]);
    });

    it('cumsum 2D axis=0 matches NumPy', async () => {
      // np.cumsum([[1, 2, 3], [4, 5, 6]], axis=0) = [[1, 2, 3], [5, 7, 9]]
      const arr = track(await NDArray.fromArray([[1, 2, 3], [4, 5, 6]]));
      const result = track(await cumsum(arr, 0));
      expect(result.toArray()).toEqual([1, 2, 3, 5, 7, 9]);
    });

    it('cumsum 2D axis=1 matches NumPy', async () => {
      // np.cumsum([[1, 2, 3], [4, 5, 6]], axis=1) = [[1, 3, 6], [4, 9, 15]]
      const arr = track(await NDArray.fromArray([[1, 2, 3], [4, 5, 6]]));
      const result = track(await cumsum(arr, 1));
      expect(result.toArray()).toEqual([1, 3, 6, 4, 9, 15]);
    });

    it('nancumsum matches NumPy examples', async () => {
      // np.nancumsum([1, np.nan, 3, 4]) = [1, 1, 4, 8]
      const arr = track(await NDArray.fromArray([1, NaN, 3, 4]));
      const result = track(await nancumsum(arr));
      expect(result.toArray()).toEqual([1, 1, 4, 8]);
    });

    it('nancumprod matches NumPy examples', async () => {
      // np.nancumprod([1, np.nan, 3, 4]) = [1, 1, 3, 12]
      const arr = track(await NDArray.fromArray([1, NaN, 3, 4]));
      const result = track(await nancumprod(arr));
      expect(result.toArray()).toEqual([1, 1, 3, 12]);
    });
  });
});
