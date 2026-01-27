import { test, expect } from '@playwright/test';

// Helper to get the base URL for imports - use full URL for dynamic imports in browser
const NUMJS_URL = 'http://localhost:3000/dist/numjs.mjs';

test.describe('NumJS Browser Support', () => {
  test.beforeEach(async ({ page }) => {
    // Navigate to the test page
    await page.goto('/tests/browser/index.html');
  });

  test('loads WASM module successfully', async ({ page }) => {
    const result = await page.evaluate(async (url) => {
      const { loadWasmModule, isWasmLoaded } = await import(url);

      const wasLoadedBefore = isWasmLoaded();
      await loadWasmModule();
      const wasLoadedAfter = isWasmLoaded();

      return { wasLoadedBefore, wasLoadedAfter };
    }, NUMJS_URL);

    expect(result.wasLoadedBefore).toBe(false);
    expect(result.wasLoadedAfter).toBe(true);
  });

  test('creates zeros array', async ({ page }) => {
    const result = await page.evaluate(async (url) => {
      const { loadWasmModule, zeros } = await import(url);
      await loadWasmModule();

      const arr = await zeros([3, 4]);
      const shape = arr.shape;
      const size = arr.size;
      const sum = arr.sum();
      arr.dispose();

      return { shape, size, sum };
    }, NUMJS_URL);

    expect(result.shape).toEqual([3, 4]);
    expect(result.size).toBe(12);
    expect(result.sum).toBe(0);
  });

  test('creates ones array', async ({ page }) => {
    const result = await page.evaluate(async (url) => {
      const { loadWasmModule, ones } = await import(url);
      await loadWasmModule();

      const arr = await ones([2, 3]);
      const shape = arr.shape;
      const sum = arr.sum();
      arr.dispose();

      return { shape, sum };
    }, NUMJS_URL);

    expect(result.shape).toEqual([2, 3]);
    expect(result.sum).toBe(6);
  });

  test('creates array from data', async ({ page }) => {
    const result = await page.evaluate(async (url) => {
      const { loadWasmModule, array } = await import(url);
      await loadWasmModule();

      const arr = await array([1, 2, 3, 4, 5]);
      const shape = arr.shape;
      const sum = arr.sum();
      const data = arr.toArray();
      arr.dispose();

      return { shape, sum, data };
    }, NUMJS_URL);

    expect(result.shape).toEqual([5]);
    expect(result.sum).toBe(15);
    expect(result.data).toEqual([1, 2, 3, 4, 5]);
  });

  test('performs arithmetic operations', async ({ page }) => {
    const result = await page.evaluate(async (url) => {
      const { loadWasmModule, array, add, multiply } = await import(url);
      await loadWasmModule();

      const a = await array([1, 2, 3]);
      const b = await array([4, 5, 6]);

      const sum = await add(a, b);
      const sumData = sum.toArray();

      const prod = await multiply(a, b);
      const prodData = prod.toArray();

      a.dispose();
      b.dispose();
      sum.dispose();
      prod.dispose();

      return { sumData, prodData };
    }, NUMJS_URL);

    expect(result.sumData).toEqual([5, 7, 9]);
    expect(result.prodData).toEqual([4, 10, 18]);
  });

  test('performs matrix multiplication', async ({ page }) => {
    const result = await page.evaluate(async (url) => {
      const { loadWasmModule, array, matmul } = await import(url);
      await loadWasmModule();

      // 2x3 matrix
      const a = await array([
        [1, 2, 3],
        [4, 5, 6],
      ]);
      // 3x2 matrix
      const b = await array([
        [7, 8],
        [9, 10],
        [11, 12],
      ]);

      const c = await matmul(a, b);
      const shape = c.shape;
      const data = c.toArray();

      a.dispose();
      b.dispose();
      c.dispose();

      return { shape, data };
    }, NUMJS_URL);

    expect(result.shape).toEqual([2, 2]);
    // [[1*7+2*9+3*11, 1*8+2*10+3*12], [4*7+5*9+6*11, 4*8+5*10+6*12]]
    // [[7+18+33, 8+20+36], [28+45+66, 32+50+72]]
    // [[58, 64], [139, 154]]
    expect(result.data).toEqual([58, 64, 139, 154]);
  });

  test('reshapes arrays', async ({ page }) => {
    const result = await page.evaluate(async (url) => {
      const { loadWasmModule, arange } = await import(url);
      await loadWasmModule();

      const arr = await arange(12);
      const reshaped = arr.reshape([3, 4]);

      const originalShape = arr.shape;
      const newShape = reshaped.shape;
      const data = reshaped.toArray();

      arr.dispose();

      return { originalShape, newShape, data };
    }, NUMJS_URL);

    expect(result.originalShape).toEqual([12]);
    expect(result.newShape).toEqual([3, 4]);
    expect(result.data).toEqual([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]);
  });

  test('calculates statistics', async ({ page }) => {
    const result = await page.evaluate(async (url) => {
      const { loadWasmModule, array, mean, std, min, max } = await import(url);
      await loadWasmModule();

      const arr = await array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10]);
      const sumVal = arr.sum();
      const meanVal = await mean(arr);
      const minVal = await min(arr);
      const maxVal = await max(arr);
      const stdVal = await std(arr);

      arr.dispose();

      return { sum: sumVal, mean: meanVal, min: minVal, max: maxVal, std: stdVal };
    }, NUMJS_URL);

    expect(result.sum).toBe(55);
    expect(result.mean).toBe(5.5);
    expect(result.min).toBe(1);
    expect(result.max).toBe(10);
    expect(result.std).toBeCloseTo(2.8722813232690143, 5);
  });

  test('supports different dtypes', async ({ page }) => {
    const result = await page.evaluate(async (url) => {
      const { loadWasmModule, zeros, array, DType } = await import(url);
      await loadWasmModule();

      // Create a float64 array (default)
      const f64 = await zeros([2, 2]);

      // Create a float32 array with explicit dtype
      const f32 = await zeros([2, 2], DType.Float32);

      // Create a regular array from integer data (will be float64 by default)
      const intArr = await array([1, 2, 3, 4]);

      const results = {
        f64Shape: f64.shape,
        f64Sum: f64.sum(),
        f32Shape: f32.shape,
        f32Sum: f32.sum(),
        intShape: intArr.shape,
        intSum: intArr.sum(),
        intData: intArr.toArray(),
      };

      f64.dispose();
      f32.dispose();
      intArr.dispose();

      return results;
    }, NUMJS_URL);

    expect(result.f64Shape).toEqual([2, 2]);
    expect(result.f64Sum).toBe(0);
    expect(result.f32Shape).toEqual([2, 2]);
    expect(result.f32Sum).toBe(0);
    expect(result.intShape).toEqual([4]);
    expect(result.intSum).toBe(10);
    expect(result.intData).toEqual([1, 2, 3, 4]);
  });

  test('handles slicing operations', async ({ page }) => {
    const result = await page.evaluate(async (url) => {
      const { loadWasmModule, arange, slice } = await import(url);
      await loadWasmModule();

      const arr = await arange(10);
      // Get elements 2-7 using slice helper
      const sliced = arr.slice([slice(2, 7)]);
      const data = sliced.toArray();

      arr.dispose();

      return { data };
    }, NUMJS_URL);

    expect(result.data).toEqual([2, 3, 4, 5, 6]);
  });

  test('performs element-wise math functions', async ({ page }) => {
    const result = await page.evaluate(async (url) => {
      const { loadWasmModule, array, sqrt } = await import(url);
      await loadWasmModule();

      const arr = await array([1, 4, 9, 16]);

      const sqrtResult = await sqrt(arr);
      const sqrtData = sqrtResult.toArray();

      arr.dispose();
      sqrtResult.dispose();

      return { sqrtData };
    }, NUMJS_URL);

    expect(result.sqrtData).toEqual([1, 2, 3, 4]);
  });

  test('handles 2D array creation and access', async ({ page }) => {
    const result = await page.evaluate(async (url) => {
      const { loadWasmModule, array } = await import(url);
      await loadWasmModule();

      const arr = await array([
        [1, 2, 3],
        [4, 5, 6],
      ]);

      const shape = arr.shape;
      const ndim = arr.ndim;
      const val = arr.get(1, 2); // Should be 6

      arr.dispose();

      return { shape, ndim, val };
    }, NUMJS_URL);

    expect(result.shape).toEqual([2, 3]);
    expect(result.ndim).toBe(2);
    expect(result.val).toBe(6);
  });

  test('random number generation works', async ({ page }) => {
    const result = await page.evaluate(async (url) => {
      const { loadWasmModule, default_rng } = await import(url);
      await loadWasmModule();

      const rng = default_rng(42); // Seeded for reproducibility
      const arr = rng.random(10);

      const shape = arr.shape;
      const size = arr.size;
      const data = arr.toArray();
      const allInRange = data.every((v: number) => v >= 0 && v < 1);

      arr.dispose();

      return { shape, size, allInRange };
    }, NUMJS_URL);

    expect(result.shape).toEqual([10]);
    expect(result.size).toBe(10);
    expect(result.allInRange).toBe(true);
  });
});

test.describe('NumJS Browser Error Handling', () => {
  test('throws error for file I/O operations in browser', async ({ page }) => {
    await page.goto('/tests/browser/index.html');

    const result = await page.evaluate(async (url) => {
      const { loadWasmModule, zeros, save } = await import(url);
      await loadWasmModule();

      const arr = await zeros([2, 2]);

      try {
        await save('/tmp/test.npy', arr);
        return { threw: false, error: null };
      } catch (e) {
        return { threw: true, error: (e as Error).message };
      } finally {
        arr.dispose();
      }
    }, NUMJS_URL);

    expect(result.threw).toBe(true);
    // Should mention that file operations aren't available in browser
    expect(result.error).toBeTruthy();
  });
});
