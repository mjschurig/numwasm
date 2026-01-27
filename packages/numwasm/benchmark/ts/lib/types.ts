/**
 * TypeScript interfaces for the benchmark system.
 */

// Configuration for benchmark execution
export interface BenchmarkConfig {
  warmupIterations: number;
  minIterations: number;
  minDurationMs: number;
}

// Result from a single benchmark run
export interface BenchmarkResult {
  meanMs: number;
  minMs: number;
  maxMs: number;
  stdMs: number;
  iterations: number;
  result?: number;
}

// A single data point with parameter combination
export interface BenchmarkDataPoint {
  params: Record<string, string | number>;
  meanMs: number;
  minMs: number;
  maxMs: number;
  stdMs: number;
  iterations: number;
  result?: number;
}

// Results for a single operation (e.g., "sum", "mean")
export interface OperationResults {
  id: string;
  name: string;
  results: BenchmarkDataPoint[];
}

// Results for a category (e.g., "reduce", "core")
export interface CategoryResults {
  id: string;
  name: string;
  operations: OperationResults[];
}

// Final combined benchmark results (for benchmarks.json)
export interface BenchmarkResults {
  version: string;
  generatedAt: string;
  environment: {
    numjsVersion: string;
    numpyVersion: string;
    nodeVersion: string;
    pythonVersion: string;
    platform: string;
  };
  categories: BenchmarkCategory[];
}

export interface BenchmarkCategory {
  id: string;
  name: string;
  operations: BenchmarkOperation[];
}

export interface BenchmarkOperation {
  id: string;
  name: string;
  results: CombinedDataPoint[];
}

export interface CombinedDataPoint {
  params: Record<string, string | number>;
  numpy: {
    meanMs: number;
    minMs: number;
    maxMs: number;
    stdMs: number;
    iterations: number;
  };
  numjs: {
    meanMs: number;
    minMs: number;
    maxMs: number;
    stdMs: number;
    iterations: number;
  };
  speedup: number;
}
