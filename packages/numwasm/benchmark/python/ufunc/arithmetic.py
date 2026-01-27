"""
Universal function benchmarks: arithmetic and math operations
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


UFUNC_MAP = {
    'add': {'fn': np.add, 'binary': True},
    'subtract': {'fn': np.subtract, 'binary': True},
    'multiply': {'fn': np.multiply, 'binary': True},
    'divide': {'fn': np.divide, 'binary': True},
    'sqrt': {'fn': np.sqrt, 'binary': False, 'needs_positive': True},
    'exp': {'fn': np.exp, 'binary': False},
    'log': {'fn': np.log, 'binary': False, 'needs_positive': True},
    'sin': {'fn': np.sin, 'binary': False},
    'cos': {'fn': np.cos, 'binary': False},
    'power': {'fn': np.power, 'binary': True},
}

UFUNC_NAMES = {
    'add': 'Add',
    'subtract': 'Subtract',
    'multiply': 'Multiply',
    'divide': 'Divide',
    'sqrt': 'Square Root',
    'exp': 'Exponential',
    'log': 'Logarithm',
    'sin': 'Sine',
    'cos': 'Cosine',
    'power': 'Power',
}


def benchmark_ufunc(name: str, runner: BenchmarkRunner) -> dict:
    """Benchmark a ufunc across array sizes."""
    results = []
    ufunc = UFUNC_MAP[name]

    for size in ARRAY_SIZES:
        print(f"    {name} size={format_size(size)}... ", end='', flush=True)

        # Generate appropriate data
        np.random.seed(RANDOM_SEED)
        if ufunc.get('needs_positive'):
            arr1 = np.random.uniform(0.1, 1.1, size)
        else:
            arr1 = np.random.randn(size)

        arr2 = None
        if ufunc['binary']:
            np.random.seed(RANDOM_SEED + 1)
            arr2 = np.random.randn(size)

        if ufunc['binary']:
            result = runner.run(lambda: ufunc['fn'](arr1, arr2))
        else:
            result = runner.run(lambda: ufunc['fn'](arr1))

        results.append({
            'params': {'size': size},
            'meanMs': result.mean_ms,
            'minMs': result.min_ms,
            'maxMs': result.max_ms,
            'stdMs': result.std_ms,
            'iterations': result.iterations,
        })

        print(f"{result.mean_ms:.4f} ms ({result.iterations} iters)")

    return {'id': name, 'name': UFUNC_NAMES[name], 'results': results}


def run_ufunc_benchmarks() -> list[dict]:
    """Run all ufunc benchmarks."""
    runner = BenchmarkRunner(BENCHMARK_CONFIG)
    operations = ['add', 'subtract', 'multiply', 'divide', 'sqrt', 'exp', 'log', 'sin', 'cos', 'power']
    results = []

    for op in operations:
        print(f"  Benchmarking {op}...")
        results.append(benchmark_ufunc(op, runner))

    return results


if __name__ == '__main__':
    import json
    results = run_ufunc_benchmarks()
    print(json.dumps(results, indent=2))
