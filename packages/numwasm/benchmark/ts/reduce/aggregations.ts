/**
 * Reduction benchmarks: sum, mean, std, var, min, max, prod
 */

import { BenchmarkRunner } from '../lib/BenchmarkRunner.js';
import { BENCHMARK_CONFIG, ARRAY_SIZES, RANDOM_SEED } from '../lib/config.js';
import { generateRandomData, formatSize } from '../lib/utils.js';
import type { OperationResults, BenchmarkDataPoint } from '../lib/types.js';

// @ts-expect-error - Types are in dist/index.d.ts but TS can't resolve mjs
import { NDArray, loadWasmModule, sum, mean, std, var_, min, max } from '../../../numjs.mjs';

type ReductionMethod = 'sum' | 'mean' | 'std' | 'var' | 'min' | 'max';

// Map method names to standalone functions
const reductionFunctions: Record<ReductionMethod, (arr: any) => number> = {
  sum: (arr) => sum(arr),
  mean: (arr) => mean(arr),
  std: (arr) => std(arr),
  var: (arr) => var_(arr),
  min: (arr) => min(arr),
  max: (arr) => max(arr),
};

async function benchmarkReduction(
  method: ReductionMethod,
  runner: BenchmarkRunner
): Promise<OperationResults> {
  const results: BenchmarkDataPoint[] = [];
  const reductionFn = reductionFunctions[method];

  for (const size of ARRAY_SIZES) {
    process.stdout.write(`    ${method} size=${formatSize(size)}... `);

    const data = generateRandomData(size, RANDOM_SEED);
    const arr = await NDArray.fromArray(data);

    try {
      const result = await runner.run(
        () => arr,
        (a) => reductionFn(a),
        () => {} // No cleanup during run, we dispose after
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
      arr.dispose();
    }
  }

  const methodNames: Record<ReductionMethod, string> = {
    sum: 'Sum',
    mean: 'Mean',
    std: 'Standard Deviation',
    var: 'Variance',
    min: 'Minimum',
    max: 'Maximum',
  };

  return {
    id: method,
    name: methodNames[method],
    results,
  };
}

export async function runReduceBenchmarks(): Promise<OperationResults[]> {
  console.log('  Loading WASM module...');
  await loadWasmModule();

  const runner = new BenchmarkRunner(BENCHMARK_CONFIG);
  const methods: ReductionMethod[] = ['sum', 'mean', 'std', 'var', 'min', 'max'];
  const results: OperationResults[] = [];

  for (const method of methods) {
    console.log(`  Benchmarking ${method}...`);
    const opResult = await benchmarkReduction(method, runner);
    results.push(opResult);
  }

  return results;
}
