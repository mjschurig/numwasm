"""
Linear algebra benchmarks: dot, matmul, norm, det, inv
"""

import numpy as np
from ..lib.benchmark_runner import BenchmarkRunner
from ..lib.config import BENCHMARK_CONFIG, MATRIX_SIZES, RANDOM_SEED


def benchmark_dot(runner: BenchmarkRunner) -> dict:
    """Benchmark np.dot for vectors."""
    results = []

    for size in MATRIX_SIZES:
        print(f"    dot size={size}... ", end='', flush=True)

        np.random.seed(RANDOM_SEED)
        arr1 = np.random.randn(size)
        np.random.seed(RANDOM_SEED + 1)
        arr2 = np.random.randn(size)

        result = runner.run(lambda: np.dot(arr1, arr2))

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

    return {'id': 'dot', 'name': 'Dot Product', 'results': results}


def benchmark_matmul(runner: BenchmarkRunner) -> dict:
    """Benchmark np.matmul for matrices."""
    results = []

    for size in MATRIX_SIZES:
        print(f"    matmul size={size}x{size}... ", end='', flush=True)

        np.random.seed(RANDOM_SEED)
        mat1 = np.random.randn(size, size)
        np.random.seed(RANDOM_SEED + 1)
        mat2 = np.random.randn(size, size)

        result = runner.run(lambda: np.matmul(mat1, mat2))

        results.append({
            'params': {'size': size},
            'meanMs': result.mean_ms,
            'minMs': result.min_ms,
            'maxMs': result.max_ms,
            'stdMs': result.std_ms,
            'iterations': result.iterations,
        })

        print(f"{result.mean_ms:.4f} ms ({result.iterations} iters)")

    return {'id': 'matmul', 'name': 'Matrix Multiply', 'results': results}


def benchmark_norm(runner: BenchmarkRunner) -> dict:
    """Benchmark np.linalg.norm for matrices."""
    results = []

    for size in MATRIX_SIZES:
        print(f"    norm size={size}x{size}... ", end='', flush=True)

        np.random.seed(RANDOM_SEED)
        mat = np.random.randn(size, size)

        result = runner.run(lambda: np.linalg.norm(mat))

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

    return {'id': 'norm', 'name': 'Matrix Norm', 'results': results}


def benchmark_det(runner: BenchmarkRunner) -> dict:
    """Benchmark np.linalg.det for matrices."""
    results = []
    # Use smaller sizes for det (O(n^3) operation)
    det_sizes = [s for s in MATRIX_SIZES if s <= 200]

    for size in det_sizes:
        print(f"    det size={size}x{size}... ", end='', flush=True)

        np.random.seed(RANDOM_SEED)
        mat = np.random.randn(size, size)

        result = runner.run(lambda: np.linalg.det(mat))

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

    return {'id': 'det', 'name': 'Determinant', 'results': results}


def benchmark_inv(runner: BenchmarkRunner) -> dict:
    """Benchmark np.linalg.inv for matrices."""
    results = []
    # Use smaller sizes for inv (O(n^3) operation)
    inv_sizes = [s for s in MATRIX_SIZES if s <= 200]

    for size in inv_sizes:
        print(f"    inv size={size}x{size}... ", end='', flush=True)

        np.random.seed(RANDOM_SEED)
        mat = np.random.randn(size, size)
        # Add diagonal dominance for well-conditioned matrix
        mat = mat + np.eye(size) * size

        result = runner.run(lambda: np.linalg.inv(mat))

        results.append({
            'params': {'size': size},
            'meanMs': result.mean_ms,
            'minMs': result.min_ms,
            'maxMs': result.max_ms,
            'stdMs': result.std_ms,
            'iterations': result.iterations,
        })

        print(f"{result.mean_ms:.4f} ms ({result.iterations} iters)")

    return {'id': 'inv', 'name': 'Matrix Inverse', 'results': results}


def run_linalg_benchmarks() -> list[dict]:
    """Run all linear algebra benchmarks."""
    runner = BenchmarkRunner(BENCHMARK_CONFIG)
    results = []

    print("  Benchmarking dot...")
    results.append(benchmark_dot(runner))

    print("  Benchmarking matmul...")
    results.append(benchmark_matmul(runner))

    print("  Benchmarking norm...")
    results.append(benchmark_norm(runner))

    print("  Benchmarking det...")
    results.append(benchmark_det(runner))

    print("  Benchmarking inv...")
    results.append(benchmark_inv(runner))

    return results


if __name__ == '__main__':
    import json
    results = run_linalg_benchmarks()
    print(json.dumps(results, indent=2))
