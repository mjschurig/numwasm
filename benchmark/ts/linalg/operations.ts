/**
 * Linear algebra benchmarks: dot, matmul, norm, det, inv
 */

import { BenchmarkRunner } from '../lib/BenchmarkRunner.js';
import { BENCHMARK_CONFIG, MATRIX_SIZES, RANDOM_SEED } from '../lib/config.js';
import { generateRandomData } from '../lib/utils.js';
import type { OperationResults, BenchmarkDataPoint } from '../lib/types.js';

// @ts-expect-error - Types are in dist/index.d.ts but TS can't resolve mjs
import { NDArray, dot, matmul, loadWasmModule } from '../../../numjs.mjs';
// @ts-expect-error
import * as linalg from '../../../numjs.mjs';

async function benchmarkDot(runner: BenchmarkRunner): Promise<OperationResults> {
  const results: BenchmarkDataPoint[] = [];

  for (const size of MATRIX_SIZES) {
    process.stdout.write(`    dot size=${size}... `);

    const data1 = generateRandomData(size, RANDOM_SEED);
    const data2 = generateRandomData(size, RANDOM_SEED + 1);

    const arr1 = await NDArray.fromArray(data1);
    const arr2 = await NDArray.fromArray(data2);

    try {
      const result = await runner.run(
        () => ({ arr1, arr2 }),
        (setup) => dot(setup.arr1, setup.arr2)
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
    } finally {
      arr1.dispose();
      arr2.dispose();
    }
  }

  return { id: 'dot', name: 'Dot Product', results };
}

async function benchmarkMatmul(runner: BenchmarkRunner): Promise<OperationResults> {
  const results: BenchmarkDataPoint[] = [];

  for (const size of MATRIX_SIZES) {
    process.stdout.write(`    matmul size=${size}x${size}... `);

    const data1 = generateRandomData(size * size, RANDOM_SEED);
    const data2 = generateRandomData(size * size, RANDOM_SEED + 1);

    const arr1 = await NDArray.fromArray(data1);
    const arr2 = await NDArray.fromArray(data2);
    const mat1 = arr1.reshape([size, size]);
    const mat2 = arr2.reshape([size, size]);
    arr1.dispose();
    arr2.dispose();

    try {
      const result = await runner.run(
        () => ({ mat1, mat2 }),
        (setup) => {
          const res = matmul(setup.mat1, setup.mat2);
          const val = res.sum();
          res.dispose();
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
    } finally {
      mat1.dispose();
      mat2.dispose();
    }
  }

  return { id: 'matmul', name: 'Matrix Multiply', results };
}

async function benchmarkNorm(runner: BenchmarkRunner): Promise<OperationResults> {
  const results: BenchmarkDataPoint[] = [];

  for (const size of MATRIX_SIZES) {
    process.stdout.write(`    norm size=${size}x${size}... `);

    const data = generateRandomData(size * size, RANDOM_SEED);
    const arr = await NDArray.fromArray(data);
    const mat = arr.reshape([size, size]);
    arr.dispose();

    try {
      const result = await runner.run(
        () => mat,
        (m) => linalg.norm(m)
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
    } finally {
      mat.dispose();
    }
  }

  return { id: 'norm', name: 'Matrix Norm', results };
}

async function benchmarkDet(runner: BenchmarkRunner): Promise<OperationResults> {
  const results: BenchmarkDataPoint[] = [];
  // Use smaller sizes for det (O(n^3) operation)
  const detSizes = MATRIX_SIZES.filter(s => s <= 200);

  for (const size of detSizes) {
    process.stdout.write(`    det size=${size}x${size}... `);

    const data = generateRandomData(size * size, RANDOM_SEED);
    const arr = await NDArray.fromArray(data);
    const mat = arr.reshape([size, size]);
    arr.dispose();

    try {
      const result = await runner.run(
        () => mat,
        (m) => linalg.det(m)
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
    } finally {
      mat.dispose();
    }
  }

  return { id: 'det', name: 'Determinant', results };
}

async function benchmarkInv(runner: BenchmarkRunner): Promise<OperationResults> {
  const results: BenchmarkDataPoint[] = [];
  // Use smaller sizes for inv (O(n^3) operation)
  const invSizes = MATRIX_SIZES.filter(s => s <= 200);

  for (const size of invSizes) {
    process.stdout.write(`    inv size=${size}x${size}... `);

    // Create a well-conditioned matrix (add diagonal dominance)
    const data = generateRandomData(size * size, RANDOM_SEED);
    const arr = await NDArray.fromArray(data);
    const mat = arr.reshape([size, size]);
    arr.dispose();

    // Add identity * size to ensure invertibility
    const eye = linalg.eye(size);
    const scaled = linalg.multiply(eye, size);
    const wellConditioned = linalg.add(mat, scaled);
    mat.dispose();
    eye.dispose();
    scaled.dispose();

    try {
      const result = await runner.run(
        () => wellConditioned,
        (m) => {
          const res = linalg.inv(m);
          const val = res.sum();
          res.dispose();
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
    } finally {
      wellConditioned.dispose();
    }
  }

  return { id: 'inv', name: 'Matrix Inverse', results };
}

export async function runLinalgBenchmarks(): Promise<OperationResults[]> {
  console.log('  Loading WASM module...');
  await loadWasmModule();

  const runner = new BenchmarkRunner(BENCHMARK_CONFIG);
  const results: OperationResults[] = [];

  console.log('  Benchmarking dot...');
  results.push(await benchmarkDot(runner));

  console.log('  Benchmarking matmul...');
  results.push(await benchmarkMatmul(runner));

  console.log('  Benchmarking norm...');
  results.push(await benchmarkNorm(runner));

  console.log('  Benchmarking det...');
  results.push(await benchmarkDet(runner));

  console.log('  Benchmarking inv...');
  results.push(await benchmarkInv(runner));

  return results;
}
