/**
 * Core array creation benchmarks: zeros, ones, arange, eye, linspace
 */

import { BenchmarkRunner } from '../lib/BenchmarkRunner.js';
import { BENCHMARK_CONFIG, ARRAY_SIZES, MATRIX_SIZES } from '../lib/config.js';
import { formatSize } from '../lib/utils.js';
import type { OperationResults, BenchmarkDataPoint } from '../lib/types.js';

// @ts-expect-error - Types are in dist/index.d.ts but TS can't resolve mjs
import { zeros, ones, arange, eye, linspace, loadWasmModule } from '../../../numjs.mjs';

async function benchmarkZeros(runner: BenchmarkRunner): Promise<OperationResults> {
  const results: BenchmarkDataPoint[] = [];

  for (const size of ARRAY_SIZES) {
    process.stdout.write(`    zeros size=${formatSize(size)}... `);

    const result = await runner.runAsync(
      () => null,
      async () => {
        const arr = await zeros([size]);
        arr.dispose();
        return size;
      }
    );

    results.push({
      params: { size },
      meanMs: result.meanMs,
      minMs: result.minMs,
      maxMs: result.maxMs,
      stdMs: result.stdMs,
      iterations: result.iterations,
    });

    console.log(`${result.meanMs.toFixed(4)} ms (${result.iterations} iters)`);
  }

  return { id: 'zeros', name: 'Zeros', results };
}

async function benchmarkOnes(runner: BenchmarkRunner): Promise<OperationResults> {
  const results: BenchmarkDataPoint[] = [];

  for (const size of ARRAY_SIZES) {
    process.stdout.write(`    ones size=${formatSize(size)}... `);

    const result = await runner.runAsync(
      () => null,
      async () => {
        const arr = await ones([size]);
        arr.dispose();
        return size;
      }
    );

    results.push({
      params: { size },
      meanMs: result.meanMs,
      minMs: result.minMs,
      maxMs: result.maxMs,
      stdMs: result.stdMs,
      iterations: result.iterations,
    });

    console.log(`${result.meanMs.toFixed(4)} ms (${result.iterations} iters)`);
  }

  return { id: 'ones', name: 'Ones', results };
}

async function benchmarkArange(runner: BenchmarkRunner): Promise<OperationResults> {
  const results: BenchmarkDataPoint[] = [];

  for (const size of ARRAY_SIZES) {
    process.stdout.write(`    arange size=${formatSize(size)}... `);

    const result = await runner.runAsync(
      () => null,
      async () => {
        const arr = await arange(size);
        arr.dispose();
        return size;
      }
    );

    results.push({
      params: { size },
      meanMs: result.meanMs,
      minMs: result.minMs,
      maxMs: result.maxMs,
      stdMs: result.stdMs,
      iterations: result.iterations,
    });

    console.log(`${result.meanMs.toFixed(4)} ms (${result.iterations} iters)`);
  }

  return { id: 'arange', name: 'Arange', results };
}

async function benchmarkEye(runner: BenchmarkRunner): Promise<OperationResults> {
  const results: BenchmarkDataPoint[] = [];

  for (const size of MATRIX_SIZES) {
    process.stdout.write(`    eye size=${size}x${size}... `);

    const result = await runner.runAsync(
      () => null,
      async () => {
        const arr = await eye(size);
        arr.dispose();
        return size * size;
      }
    );

    results.push({
      params: { size },
      meanMs: result.meanMs,
      minMs: result.minMs,
      maxMs: result.maxMs,
      stdMs: result.stdMs,
      iterations: result.iterations,
    });

    console.log(`${result.meanMs.toFixed(4)} ms (${result.iterations} iters)`);
  }

  return { id: 'eye', name: 'Eye (Identity)', results };
}

async function benchmarkLinspace(runner: BenchmarkRunner): Promise<OperationResults> {
  const results: BenchmarkDataPoint[] = [];

  for (const size of ARRAY_SIZES) {
    process.stdout.write(`    linspace size=${formatSize(size)}... `);

    const result = await runner.runAsync(
      () => null,
      async () => {
        const arr = await linspace(0, 1, size);
        arr.dispose();
        return size;
      }
    );

    results.push({
      params: { size },
      meanMs: result.meanMs,
      minMs: result.minMs,
      maxMs: result.maxMs,
      stdMs: result.stdMs,
      iterations: result.iterations,
    });

    console.log(`${result.meanMs.toFixed(4)} ms (${result.iterations} iters)`);
  }

  return { id: 'linspace', name: 'Linspace', results };
}

export async function runCoreBenchmarks(): Promise<OperationResults[]> {
  console.log('  Loading WASM module...');
  await loadWasmModule();

  const runner = new BenchmarkRunner(BENCHMARK_CONFIG);
  const results: OperationResults[] = [];

  console.log('  Benchmarking zeros...');
  results.push(await benchmarkZeros(runner));

  console.log('  Benchmarking ones...');
  results.push(await benchmarkOnes(runner));

  console.log('  Benchmarking arange...');
  results.push(await benchmarkArange(runner));

  console.log('  Benchmarking eye...');
  results.push(await benchmarkEye(runner));

  console.log('  Benchmarking linspace...');
  results.push(await benchmarkLinspace(runner));

  return results;
}
