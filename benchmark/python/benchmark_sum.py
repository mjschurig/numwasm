#!/usr/bin/env python3
"""
Benchmark NumPy sum() performance across various array sizes.

Measures execution time for array sizes from 10^2 to 10^7 elements.
Results are saved as JSON for comparison with NumJS.
"""

import json
import time
import sys
from pathlib import Path

try:
    import numpy as np
except ImportError:
    print("Error: NumPy is not installed. Run: pip install numpy", file=sys.stderr)
    sys.exit(1)

# Configuration
ARRAY_SIZES = [10**2, 10**3, 10**4, 10**5, 10**6, 10**7]
WARMUP_ITERATIONS = 3
MIN_ITERATIONS = 10
MIN_DURATION_SEC = 1.0  # Run for at least 1 second per size


def benchmark_sum(size: int) -> dict:
    """Benchmark sum() for a given array size."""
    # Create array with reproducible random data
    np.random.seed(42)
    arr = np.random.randn(size).astype(np.float64)

    # Warmup runs (primes CPU caches, etc.)
    for _ in range(WARMUP_ITERATIONS):
        _ = arr.sum()

    # Timed runs - adaptive iteration count
    times = []
    start_total = time.perf_counter()
    iterations = 0
    result = 0.0

    while iterations < MIN_ITERATIONS or (time.perf_counter() - start_total) < MIN_DURATION_SEC:
        start = time.perf_counter()
        result = arr.sum()
        end = time.perf_counter()
        times.append(end - start)
        iterations += 1

    return {
        "size": size,
        "iterations": iterations,
        "times_ms": [t * 1000 for t in times],
        "mean_ms": (sum(times) / len(times)) * 1000,
        "min_ms": min(times) * 1000,
        "max_ms": max(times) * 1000,
        "std_ms": (np.std(times) * 1000) if len(times) > 1 else 0,
        "result": float(result),  # For cross-validation
    }


def main():
    print(f"NumPy sum() Benchmark")
    print(f"NumPy version: {np.__version__}")
    print(f"Array sizes: {', '.join(f'10^{len(str(s))-1}' for s in ARRAY_SIZES)}")
    print("-" * 50)

    results = {
        "library": "numpy",
        "version": np.__version__,
        "config": {
            "warmup_iterations": WARMUP_ITERATIONS,
            "min_iterations": MIN_ITERATIONS,
            "min_duration_sec": MIN_DURATION_SEC,
        },
        "benchmarks": []
    }

    for size in ARRAY_SIZES:
        print(f"Benchmarking size {size:>10,}...", end=" ", flush=True)
        result = benchmark_sum(size)
        results["benchmarks"].append(result)
        print(f"Mean: {result['mean_ms']:>10.4f} ms  ({result['iterations']} iterations)")

    # Output to JSON
    output_dir = Path(__file__).parent.parent / "results"
    output_dir.mkdir(exist_ok=True)
    output_file = output_dir / "numpy_results.json"

    with open(output_file, "w") as f:
        json.dump(results, f, indent=2)

    print("-" * 50)
    print(f"Results saved to: {output_file}")


if __name__ == "__main__":
    main()
