/**
 * Random number generation benchmarks: uniform, normal, randint
 */

import { BenchmarkRunner } from '../lib/BenchmarkRunner.js';
import { BENCHMARK_CONFIG, ARRAY_SIZES, RANDOM_SEED } from '../lib/config.js';
import { formatSize } from '../lib/utils.js';
import type { OperationResults, BenchmarkDataPoint } from '../lib/types.js';

// @ts-expect-error - Types are in dist/index.d.ts but TS can't resolve mjs
import { loadWasmModule } from '../../../numjs.mjs';
// @ts-expect-error
import * as random from '../../../numjs.mjs';

async function benchmarkUniform(runner: BenchmarkRunner): Promise<OperationResults> {
  const results: BenchmarkDataPoint[] = [];

  for (const size of ARRAY_SIZES) {
    process.stdout.write(`    uniform size=${formatSize(size)}... `);

    // Seed the generator before each benchmark for reproducibility
    random.seed(RANDOM_SEED);

    const result = runner.runSync(
      () => null,
      () => {
        const arr = random.uniform(0, 1, [size]);
        const val = arr.sum();
        arr.dispose();
        return val;
      }
    );

    results.push({
      params: { size },
      meanMs: result.meanMs,
      minMs: result.minMs,
      maxMs: result.maxMs,
      stdMs: result.stdMs,
      iterations: result.iterations,
      result: result.result,
    });

    console.log(`${result.meanMs.toFixed(4)} ms (${result.iterations} iters)`);
  }

  return { id: 'uniform', name: 'Uniform', results };
}

async function benchmarkNormal(runner: BenchmarkRunner): Promise<OperationResults> {
  const results: BenchmarkDataPoint[] = [];

  for (const size of ARRAY_SIZES) {
    process.stdout.write(`    normal size=${formatSize(size)}... `);

    random.seed(RANDOM_SEED);

    const result = runner.runSync(
      () => null,
      () => {
        const arr = random.normal(0, 1, [size]);
        const val = arr.sum();
        arr.dispose();
        return val;
      }
    );

    results.push({
      params: { size },
      meanMs: result.meanMs,
      minMs: result.minMs,
      maxMs: result.maxMs,
      stdMs: result.stdMs,
      iterations: result.iterations,
      result: result.result,
    });

    console.log(`${result.meanMs.toFixed(4)} ms (${result.iterations} iters)`);
  }

  return { id: 'normal', name: 'Normal', results };
}

async function benchmarkRandint(runner: BenchmarkRunner): Promise<OperationResults> {
  const results: BenchmarkDataPoint[] = [];

  for (const size of ARRAY_SIZES) {
    process.stdout.write(`    randint size=${formatSize(size)}... `);

    random.seed(RANDOM_SEED);

    const result = runner.runSync(
      () => null,
      () => {
        const arr = random.randint(0, 100, [size]);
        const val = arr.sum();
        arr.dispose();
        return val;
      }
    );

    results.push({
      params: { size },
      meanMs: result.meanMs,
      minMs: result.minMs,
      maxMs: result.maxMs,
      stdMs: result.stdMs,
      iterations: result.iterations,
      result: result.result,
    });

    console.log(`${result.meanMs.toFixed(4)} ms (${result.iterations} iters)`);
  }

  return { id: 'randint', name: 'Random Integer', results };
}

export async function runRandomBenchmarks(): Promise<OperationResults[]> {
  console.log('  Loading WASM module...');
  await loadWasmModule();

  const runner = new BenchmarkRunner(BENCHMARK_CONFIG);
  const results: OperationResults[] = [];

  console.log('  Benchmarking uniform...');
  results.push(await benchmarkUniform(runner));

  console.log('  Benchmarking normal...');
  results.push(await benchmarkNormal(runner));

  console.log('  Benchmarking randint...');
  results.push(await benchmarkRandint(runner));

  return results;
}
