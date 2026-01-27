/**
 * Core benchmark execution framework.
 * Handles warmup, timing, and statistics calculation.
 */

import { performance } from 'perf_hooks';
import type { BenchmarkConfig, BenchmarkResult } from './types.js';

export class BenchmarkRunner {
  private config: BenchmarkConfig;

  constructor(config: BenchmarkConfig) {
    this.config = config;
  }

  /**
   * Run a benchmark with setup, execution, and optional cleanup.
   *
   * @param setupFn - Function to prepare test data (not timed)
   * @param benchFn - Function to benchmark (timed)
   * @param cleanupFn - Optional cleanup function (not timed)
   */
  async run<T>(
    setupFn: () => T | Promise<T>,
    benchFn: (setup: T) => number | void,
    cleanupFn?: (setup: T) => void
  ): Promise<BenchmarkResult> {
    const setup = await setupFn();

    // Warmup runs (not counted)
    for (let i = 0; i < this.config.warmupIterations; i++) {
      benchFn(setup);
    }

    // Timed runs - adaptive iteration count
    const times: number[] = [];
    const startTotal = performance.now();
    let result: number | void = undefined;

    while (
      times.length < this.config.minIterations ||
      performance.now() - startTotal < this.config.minDurationMs
    ) {
      const start = performance.now();
      result = benchFn(setup);
      times.push(performance.now() - start);
    }

    // Cleanup
    if (cleanupFn) {
      cleanupFn(setup);
    }

    return this.computeStats(times, result!);
  }

  /**
   * Run a synchronous benchmark (simpler API for sync functions).
   */
  runSync<T>(
    setupFn: () => T,
    benchFn: (setup: T) => number | void,
    cleanupFn?: (setup: T) => void
  ): BenchmarkResult {
    const setup = setupFn();

    // Warmup runs (not counted)
    for (let i = 0; i < this.config.warmupIterations; i++) {
      benchFn(setup);
    }

    // Timed runs - adaptive iteration count
    const times: number[] = [];
    const startTotal = performance.now();
    let result: number | void = undefined;

    while (
      times.length < this.config.minIterations ||
      performance.now() - startTotal < this.config.minDurationMs
    ) {
      const start = performance.now();
      result = benchFn(setup);
      times.push(performance.now() - start);
    }

    // Cleanup
    if (cleanupFn) {
      cleanupFn(setup);
    }

    return this.computeStats(times, result!);
  }

  /**
   * Run an async benchmark where the benchmark function returns a Promise.
   * Used for operations that are inherently async (like WASM array creation).
   */
  async runAsync<T>(
    setupFn: () => T | Promise<T>,
    benchFn: (setup: T) => Promise<number | void>,
    cleanupFn?: (setup: T) => void | Promise<void>
  ): Promise<BenchmarkResult> {
    const setup = await setupFn();

    // Warmup runs (not counted)
    for (let i = 0; i < this.config.warmupIterations; i++) {
      await benchFn(setup);
    }

    // Timed runs - adaptive iteration count
    const times: number[] = [];
    const startTotal = performance.now();
    let result: number | void = undefined;

    while (
      times.length < this.config.minIterations ||
      performance.now() - startTotal < this.config.minDurationMs
    ) {
      const start = performance.now();
      result = await benchFn(setup);
      times.push(performance.now() - start);
    }

    // Cleanup
    if (cleanupFn) {
      await cleanupFn(setup);
    }

    return this.computeStats(times, result!);
  }

  private computeStats(times: number[], result?: number | void): BenchmarkResult {
    const n = times.length;
    const meanMs = times.reduce((a, b) => a + b, 0) / n;
    const variance = times.reduce((sum, t) => sum + (t - meanMs) ** 2, 0) / n;

    return {
      meanMs,
      minMs: times.reduce((a, b) => Math.min(a, b), Infinity),
      maxMs: times.reduce((a, b) => Math.max(a, b), -Infinity),
      stdMs: Math.sqrt(variance),
      iterations: n,
      result: typeof result === 'number' ? result : undefined,
    };
  }
}

/**
 * Create a default benchmark runner with standard configuration.
 */
export function createRunner(config?: Partial<BenchmarkConfig>): BenchmarkRunner {
  const defaultConfig: BenchmarkConfig = {
    warmupIterations: 3,
    minIterations: 10,
    minDurationMs: 1000,
  };

  return new BenchmarkRunner({ ...defaultConfig, ...config });
}
