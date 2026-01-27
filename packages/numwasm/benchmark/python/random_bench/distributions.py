"""
Random number generation benchmarks: uniform, normal, randint
"""

import numpy as np
from ..lib.benchmark_runner import BenchmarkRunner
from ..lib.config import BENCHMARK_CONFIG, ARRAY_SIZES, RANDOM_SEED


def format_size(size: int) -> str:
    """Format a size number for display."""
    if size >= 1_000_000:
        return f"{size // 1_000_000}M"
    elif size >= 1_000:
        return f"{size // 1_000}K"
    return str(size)


def benchmark_uniform(runner: BenchmarkRunner) -> dict:
    """Benchmark np.random.uniform."""
    results = []

    for size in ARRAY_SIZES:
        print(f"    uniform size={format_size(size)}... ", end='', flush=True)

        np.random.seed(RANDOM_SEED)
        result = runner.run(lambda s=size: np.random.uniform(0, 1, s))

        results.append({
            'params': {'size': size},
            'meanMs': result.mean_ms,
            'minMs': result.min_ms,
            'maxMs': result.max_ms,
            'stdMs': result.std_ms,
            'iterations': result.iterations,
        })

        print(f"{result.mean_ms:.4f} ms ({result.iterations} iters)")

    return {'id': 'uniform', 'name': 'Uniform', 'results': results}


def benchmark_normal(runner: BenchmarkRunner) -> dict:
    """Benchmark np.random.normal."""
    results = []

    for size in ARRAY_SIZES:
        print(f"    normal size={format_size(size)}... ", end='', flush=True)

        np.random.seed(RANDOM_SEED)
        result = runner.run(lambda s=size: np.random.normal(0, 1, s))

        results.append({
            'params': {'size': size},
            'meanMs': result.mean_ms,
            'minMs': result.min_ms,
            'maxMs': result.max_ms,
            'stdMs': result.std_ms,
            'iterations': result.iterations,
        })

        print(f"{result.mean_ms:.4f} ms ({result.iterations} iters)")

    return {'id': 'normal', 'name': 'Normal', 'results': results}


def benchmark_randint(runner: BenchmarkRunner) -> dict:
    """Benchmark np.random.randint."""
    results = []

    for size in ARRAY_SIZES:
        print(f"    randint size={format_size(size)}... ", end='', flush=True)

        np.random.seed(RANDOM_SEED)
        result = runner.run(lambda s=size: np.random.randint(0, 100, s))

        results.append({
            'params': {'size': size},
            'meanMs': result.mean_ms,
            'minMs': result.min_ms,
            'maxMs': result.max_ms,
            'stdMs': result.std_ms,
            'iterations': result.iterations,
        })

        print(f"{result.mean_ms:.4f} ms ({result.iterations} iters)")

    return {'id': 'randint', 'name': 'Random Integer', 'results': results}


def run_random_benchmarks() -> list[dict]:
    """Run all random benchmarks."""
    runner = BenchmarkRunner(BENCHMARK_CONFIG)
    results = []

    print("  Benchmarking uniform...")
    results.append(benchmark_uniform(runner))

    print("  Benchmarking normal...")
    results.append(benchmark_normal(runner))

    print("  Benchmarking randint...")
    results.append(benchmark_randint(runner))

    return results


if __name__ == '__main__':
    import json
    results = run_random_benchmarks()
    print(json.dumps(results, indent=2))
