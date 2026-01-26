"""
Core array creation benchmarks: zeros, ones, arange, eye, linspace
"""

import numpy as np
from ..lib.benchmark_runner import BenchmarkRunner
from ..lib.config import BENCHMARK_CONFIG, ARRAY_SIZES, MATRIX_SIZES


def format_size(size: int) -> str:
    """Format a size number for display."""
    if size >= 1_000_000:
        return f"{size // 1_000_000}M"
    elif size >= 1_000:
        return f"{size // 1_000}K"
    return str(size)


def benchmark_zeros(runner: BenchmarkRunner) -> dict:
    """Benchmark np.zeros across array sizes."""
    results = []

    for size in ARRAY_SIZES:
        print(f"    zeros size={format_size(size)}... ", end='', flush=True)

        result = runner.run(lambda s=size: np.zeros(s))

        results.append({
            'params': {'size': size},
            'meanMs': result.mean_ms,
            'minMs': result.min_ms,
            'maxMs': result.max_ms,
            'stdMs': result.std_ms,
            'iterations': result.iterations,
        })

        print(f"{result.mean_ms:.4f} ms ({result.iterations} iters)")

    return {'id': 'zeros', 'name': 'Zeros', 'results': results}


def benchmark_ones(runner: BenchmarkRunner) -> dict:
    """Benchmark np.ones across array sizes."""
    results = []

    for size in ARRAY_SIZES:
        print(f"    ones size={format_size(size)}... ", end='', flush=True)

        result = runner.run(lambda s=size: np.ones(s))

        results.append({
            'params': {'size': size},
            'meanMs': result.mean_ms,
            'minMs': result.min_ms,
            'maxMs': result.max_ms,
            'stdMs': result.std_ms,
            'iterations': result.iterations,
        })

        print(f"{result.mean_ms:.4f} ms ({result.iterations} iters)")

    return {'id': 'ones', 'name': 'Ones', 'results': results}


def benchmark_arange(runner: BenchmarkRunner) -> dict:
    """Benchmark np.arange across array sizes."""
    results = []

    for size in ARRAY_SIZES:
        print(f"    arange size={format_size(size)}... ", end='', flush=True)

        result = runner.run(lambda s=size: np.arange(s))

        results.append({
            'params': {'size': size},
            'meanMs': result.mean_ms,
            'minMs': result.min_ms,
            'maxMs': result.max_ms,
            'stdMs': result.std_ms,
            'iterations': result.iterations,
        })

        print(f"{result.mean_ms:.4f} ms ({result.iterations} iters)")

    return {'id': 'arange', 'name': 'Arange', 'results': results}


def benchmark_eye(runner: BenchmarkRunner) -> dict:
    """Benchmark np.eye across matrix sizes."""
    results = []

    for size in MATRIX_SIZES:
        print(f"    eye size={size}x{size}... ", end='', flush=True)

        result = runner.run(lambda s=size: np.eye(s))

        results.append({
            'params': {'size': size},
            'meanMs': result.mean_ms,
            'minMs': result.min_ms,
            'maxMs': result.max_ms,
            'stdMs': result.std_ms,
            'iterations': result.iterations,
        })

        print(f"{result.mean_ms:.4f} ms ({result.iterations} iters)")

    return {'id': 'eye', 'name': 'Eye (Identity)', 'results': results}


def benchmark_linspace(runner: BenchmarkRunner) -> dict:
    """Benchmark np.linspace across array sizes."""
    results = []

    for size in ARRAY_SIZES:
        print(f"    linspace size={format_size(size)}... ", end='', flush=True)

        result = runner.run(lambda s=size: np.linspace(0, 1, s))

        results.append({
            'params': {'size': size},
            'meanMs': result.mean_ms,
            'minMs': result.min_ms,
            'maxMs': result.max_ms,
            'stdMs': result.std_ms,
            'iterations': result.iterations,
        })

        print(f"{result.mean_ms:.4f} ms ({result.iterations} iters)")

    return {'id': 'linspace', 'name': 'Linspace', 'results': results}


def run_core_benchmarks() -> list[dict]:
    """Run all core creation benchmarks."""
    runner = BenchmarkRunner(BENCHMARK_CONFIG)
    results = []

    print("  Benchmarking zeros...")
    results.append(benchmark_zeros(runner))

    print("  Benchmarking ones...")
    results.append(benchmark_ones(runner))

    print("  Benchmarking arange...")
    results.append(benchmark_arange(runner))

    print("  Benchmarking eye...")
    results.append(benchmark_eye(runner))

    print("  Benchmarking linspace...")
    results.append(benchmark_linspace(runner))

    return results


if __name__ == '__main__':
    import json
    results = run_core_benchmarks()
    print(json.dumps(results, indent=2))
