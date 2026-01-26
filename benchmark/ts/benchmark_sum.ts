/**
 * Benchmark NumJS sum() performance across various array sizes.
 *
 * Measures execution time for array sizes from 10^2 to 10^7 elements.
 * Results are saved as JSON for comparison with NumPy.
 */

import { writeFileSync, mkdirSync } from 'fs';
import { dirname, join } from 'path';
import { fileURLToPath } from 'url';
import { performance } from 'perf_hooks';

// Import from the built package
// The path is relative to the compiled output in dist/benchmark/ts/
// @ts-expect-error - Types are in dist/index.d.ts but TS can't resolve mjs
import { NDArray, loadWasmModule } from '../../numjs.mjs';

const __filename = fileURLToPath(import.meta.url);
const __dirname = dirname(__filename);

// Configuration (match Python exactly)
const ARRAY_SIZES = [1e2, 1e3, 1e4, 1e5, 1e6, 1e7];
const WARMUP_ITERATIONS = 3;
const MIN_ITERATIONS = 10;
const MIN_DURATION_MS = 1000; // 1 second

/**
 * Seeded PRNG (Mulberry32) for reproducibility.
 * This ensures consistent random data across runs.
 */
function createSeededRandom(seed: number): () => number {
  return function () {
    let t = (seed += 0x6d2b79f5);
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

/**
 * Box-Muller transform for normal distribution.
 * Generates numbers similar to numpy.random.randn().
 */
function randn(random: () => number): number {
  const u1 = random();
  const u2 = random();
  return Math.sqrt(-2 * Math.log(u1)) * Math.cos(2 * Math.PI * u2);
}

interface BenchmarkResult {
  size: number;
  iterations: number;
  times_ms: number[];
  mean_ms: number;
  min_ms: number;
  max_ms: number;
  std_ms: number;
  result: number;
}

/**
 * Benchmark sum() for a given array size.
 */
async function benchmarkSum(size: number): Promise<BenchmarkResult> {
  // Create array with reproducible random data
  const random = createSeededRandom(42);
  const data = Array.from({ length: size }, () => randn(random));

  const arr = await NDArray.fromArray(data);

  try {
    // Warmup runs
    for (let i = 0; i < WARMUP_ITERATIONS; i++) {
      arr.sum();
    }

    // Timed runs - adaptive iteration count
    const times: number[] = [];
    const startTotal = performance.now();
    let iterations = 0;
    let result = 0;

    while (
      iterations < MIN_ITERATIONS ||
      performance.now() - startTotal < MIN_DURATION_MS
    ) {
      const start = performance.now();
      result = arr.sum();
      const end = performance.now();
      times.push(end - start);
      iterations++;
    }

    const meanMs = times.reduce((a, b) => a + b, 0) / times.length;
    const variance =
      times.reduce((sum, t) => sum + (t - meanMs) ** 2, 0) / times.length;
    const stdMs = Math.sqrt(variance);

    // Use reduce instead of spread to avoid stack overflow on large arrays
    const minMs = times.reduce((a, b) => Math.min(a, b), Infinity);
    const maxMs = times.reduce((a, b) => Math.max(a, b), -Infinity);

    return {
      size,
      iterations,
      times_ms: times,
      mean_ms: meanMs,
      min_ms: minMs,
      max_ms: maxMs,
      std_ms: stdMs,
      result,
    };
  } finally {
    arr.dispose();
  }
}

async function main() {
  console.log('NumJS sum() Benchmark');
  console.log('Loading WASM module...');
  await loadWasmModule();

  console.log(
    `Array sizes: ${ARRAY_SIZES.map((s) => `10^${Math.log10(s)}`).join(', ')}`
  );
  console.log('-'.repeat(50));

  const results = {
    library: 'numjs',
    version: '0.1.0',
    config: {
      warmup_iterations: WARMUP_ITERATIONS,
      min_iterations: MIN_ITERATIONS,
      min_duration_ms: MIN_DURATION_MS,
    },
    benchmarks: [] as BenchmarkResult[],
  };

  for (const size of ARRAY_SIZES) {
    const sizeStr = size.toLocaleString().padStart(12);
    process.stdout.write(`Benchmarking size ${sizeStr}... `);

    const result = await benchmarkSum(size);
    results.benchmarks.push(result);

    console.log(
      `Mean: ${result.mean_ms.toFixed(4).padStart(10)} ms  (${result.iterations} iterations)`
    );
  }

  // Output to JSON (use project root benchmark/results, not dist)
  const projectRoot = join(__dirname, '..', '..', '..');
  const outputDir = join(projectRoot, 'benchmark', 'results');
  mkdirSync(outputDir, { recursive: true });
  const outputFile = join(outputDir, 'numjs_results.json');

  writeFileSync(outputFile, JSON.stringify(results, null, 2));

  console.log('-'.repeat(50));
  console.log(`Results saved to: ${outputFile}`);
}

main().catch(console.error);
