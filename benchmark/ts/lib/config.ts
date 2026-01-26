/**
 * Shared configuration for all benchmarks.
 * Based on patterns from NumPy's benchmark suite.
 */

import type { BenchmarkConfig } from './types.js';

// Default benchmark configuration
export const BENCHMARK_CONFIG: BenchmarkConfig = {
  warmupIterations: 3,
  minIterations: 10,
  minDurationMs: 1000, // 1 second
};

// Standard array sizes for benchmarks (powers of 10)
export const ARRAY_SIZES = [100, 1_000, 10_000, 100_000, 1_000_000];

// Smaller sizes for expensive operations (like matrix decompositions)
export const SMALL_SIZES = [10, 50, 100, 500, 1000];

// Matrix sizes for linear algebra benchmarks
export const MATRIX_SIZES = [10, 50, 100, 500];

// Standard dtypes for numeric benchmarks
export const DTYPES = ['float32', 'float64'] as const;

// Extended dtypes (when needed)
export const DTYPES_EXTENDED = [
  'int16',
  'int32',
  'int64',
  'float32',
  'float64',
] as const;

// Random seed for reproducible benchmarks
export const RANDOM_SEED = 42;

// Category names for organization
export const CATEGORY_NAMES: Record<string, string> = {
  core: 'Core Operations',
  reduce: 'Reductions',
  ufunc: 'Universal Functions',
  linalg: 'Linear Algebra',
  random: 'Random',
};
