/**
 * Universal function benchmarks: arithmetic and math operations
 */

import { BenchmarkRunner } from '../lib/BenchmarkRunner.js';
import { BENCHMARK_CONFIG, ARRAY_SIZES, RANDOM_SEED } from '../lib/config.js';
import { generateRandomData, generateUniformData, formatSize } from '../lib/utils.js';
import type { OperationResults, BenchmarkDataPoint } from '../lib/types.js';

// @ts-expect-error - Types are in dist/index.d.ts but TS can't resolve mjs
import { NDArray, add, subtract, multiply, divide, sqrt, exp, log, sin, cos, power, loadWasmModule } from '../../../numjs.mjs';

type UfuncName = 'add' | 'subtract' | 'multiply' | 'divide' | 'sqrt' | 'exp' | 'log' | 'sin' | 'cos' | 'power';

const ufuncMap: Record<UfuncName, { fn: Function; binary: boolean; needsPositive?: boolean }> = {
  add: { fn: add, binary: true },
  subtract: { fn: subtract, binary: true },
  multiply: { fn: multiply, binary: true },
  divide: { fn: divide, binary: true },
  sqrt: { fn: sqrt, binary: false, needsPositive: true },
  exp: { fn: exp, binary: false },
  log: { fn: log, binary: false, needsPositive: true },
  sin: { fn: sin, binary: false },
  cos: { fn: cos, binary: false },
  power: { fn: power, binary: true },
};

const ufuncNames: Record<UfuncName, string> = {
  add: 'Add',
  subtract: 'Subtract',
  multiply: 'Multiply',
  divide: 'Divide',
  sqrt: 'Square Root',
  exp: 'Exponential',
  log: 'Logarithm',
  sin: 'Sine',
  cos: 'Cosine',
  power: 'Power',
};

async function benchmarkUfunc(
  name: UfuncName,
  runner: BenchmarkRunner
): Promise<OperationResults> {
  const results: BenchmarkDataPoint[] = [];
  const ufunc = ufuncMap[name];

  for (const size of ARRAY_SIZES) {
    process.stdout.write(`    ${name} size=${formatSize(size)}... `);

    // Generate appropriate data
    let data1: number[];
    if (ufunc.needsPositive) {
      // Use uniform [0.1, 1.1] for positive values (avoid 0 for log/sqrt)
      data1 = generateUniformData(size, RANDOM_SEED).map(x => x + 0.1);
    } else {
      data1 = generateRandomData(size, RANDOM_SEED);
    }

    const arr1 = await NDArray.fromArray(data1);
    let arr2: typeof arr1 | null = null;

    if (ufunc.binary) {
      const data2 = generateRandomData(size, RANDOM_SEED + 1);
      arr2 = await NDArray.fromArray(data2);
    }

    try {
      const result = await runner.run(
        () => ({ arr1, arr2 }),
        (setup) => {
          let res: typeof arr1;
          if (ufunc.binary && setup.arr2) {
            res = ufunc.fn(setup.arr1, setup.arr2);
          } else {
            res = ufunc.fn(setup.arr1);
          }
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
      arr1.dispose();
      if (arr2) arr2.dispose();
    }
  }

  return { id: name, name: ufuncNames[name], results };
}

export async function runUfuncBenchmarks(): Promise<OperationResults[]> {
  console.log('  Loading WASM module...');
  await loadWasmModule();

  const runner = new BenchmarkRunner(BENCHMARK_CONFIG);
  const operations: UfuncName[] = ['add', 'subtract', 'multiply', 'divide', 'sqrt', 'exp', 'log', 'sin', 'cos', 'power'];
  const results: OperationResults[] = [];

  for (const op of operations) {
    console.log(`  Benchmarking ${op}...`);
    results.push(await benchmarkUfunc(op, runner));
  }

  return results;
}
