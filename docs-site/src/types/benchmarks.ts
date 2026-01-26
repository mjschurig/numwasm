/**
 * TypeScript types for benchmark data.
 */

export interface HardwareInfo {
  cpu: string;
  cores: number;
  memory: string;
  os: string;
  arch: string;
}

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
  hardware?: HardwareInfo;
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
  results: BenchmarkDataPoint[];
}

export interface BenchmarkDataPoint {
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
