"""
Reduction benchmarks: sum, mean, std, var, min, max, prod
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


def benchmark_reduction(method: str, runner: BenchmarkRunner) -> dict:
    """Benchmark a reduction method across array sizes."""
    results = []

    method_names = {
        'sum': 'Sum',
        'mean': 'Mean',
        'std': 'Standard Deviation',
        'var': 'Variance',
        'min': 'Minimum',
        'max': 'Maximum',
        'prod': 'Product',
    }

    for size in ARRAY_SIZES:
        print(f"    {method} size={format_size(size)}... ", end='', flush=True)

        # Create array with reproducible random data
        np.random.seed(RANDOM_SEED)
        arr = np.random.randn(size)

        # Get the method to benchmark
        fn = getattr(arr, method)

        result = runner.run(fn)

        results.append({
            'params': {'size': size},
            'meanMs': result.mean_ms,
            'minMs': result.min_ms,
            'maxMs': result.max_ms,
            'stdMs': result.std_ms,
            'iterations': result.iterations,
            'result': result.result,
        })

        print(f"{result.mean_ms:.4f} ms ({result.iterations} iters)")

    return {
        'id': method,
        'name': method_names.get(method, method),
        'results': results,
    }


def run_reduce_benchmarks() -> list[dict]:
    """Run all reduction benchmarks."""
    runner = BenchmarkRunner(BENCHMARK_CONFIG)
    methods = ['sum', 'mean', 'std', 'var', 'min', 'max']
    results = []

    for method in methods:
        print(f"  Benchmarking {method}...")
        op_result = benchmark_reduction(method, runner)
        results.append(op_result)

    return results


if __name__ == '__main__':
    import json
    results = run_reduce_benchmarks()
    print(json.dumps(results, indent=2))
