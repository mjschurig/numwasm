"""
Core benchmark execution framework for Python/NumPy.
"""

import time
from typing import Callable, TypeVar, Optional, Any
from dataclasses import dataclass

T = TypeVar('T')


@dataclass
class BenchmarkResult:
    """Result from a single benchmark run."""
    mean_ms: float
    min_ms: float
    max_ms: float
    std_ms: float
    iterations: int
    result: Optional[float] = None


class BenchmarkRunner:
    """
    Benchmark runner with warmup, timing, and statistics calculation.
    """

    def __init__(self, config: dict):
        self.warmup_iterations = config.get('warmup_iterations', 3)
        self.min_iterations = config.get('min_iterations', 10)
        self.min_duration_sec = config.get('min_duration_sec', 1.0)

    def run(
        self,
        bench_fn: Callable[[], Any],
        cleanup_fn: Optional[Callable[[], None]] = None
    ) -> BenchmarkResult:
        """
        Run a benchmark function.

        Args:
            bench_fn: Function to benchmark (called repeatedly)
            cleanup_fn: Optional cleanup function (called once after benchmarking)

        Returns:
            BenchmarkResult with timing statistics
        """
        # Warmup runs (not counted)
        for _ in range(self.warmup_iterations):
            bench_fn()

        # Timed runs - adaptive iteration count
        times = []
        start_total = time.perf_counter()
        result = None

        while (len(times) < self.min_iterations or
               time.perf_counter() - start_total < self.min_duration_sec):
            start = time.perf_counter()
            result = bench_fn()
            end = time.perf_counter()
            times.append((end - start) * 1000)  # Convert to milliseconds

        # Cleanup
        if cleanup_fn:
            cleanup_fn()

        return self._compute_stats(times, result)

    def _compute_stats(self, times: list[float], result: Any) -> BenchmarkResult:
        """Compute statistics from timing data."""
        n = len(times)
        mean_ms = sum(times) / n
        variance = sum((t - mean_ms) ** 2 for t in times) / n

        # Convert result to float if it's a scalar
        result_val = None
        if result is not None:
            try:
                result_val = float(result)
            except (TypeError, ValueError):
                pass

        return BenchmarkResult(
            mean_ms=mean_ms,
            min_ms=min(times),
            max_ms=max(times),
            std_ms=variance ** 0.5,
            iterations=n,
            result=result_val,
        )


def create_runner(config: Optional[dict] = None) -> BenchmarkRunner:
    """Create a benchmark runner with default or custom configuration."""
    default_config = {
        'warmup_iterations': 3,
        'min_iterations': 10,
        'min_duration_sec': 1.0,
    }
    if config:
        default_config.update(config)
    return BenchmarkRunner(default_config)
