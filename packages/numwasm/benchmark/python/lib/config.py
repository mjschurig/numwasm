"""
Shared configuration for all Python benchmarks.
"""

# Default benchmark configuration
BENCHMARK_CONFIG = {
    'warmup_iterations': 3,
    'min_iterations': 10,
    'min_duration_sec': 1.0,  # 1 second
}

# Standard array sizes for benchmarks
ARRAY_SIZES = [100, 1_000, 10_000, 100_000, 1_000_000]

# Smaller sizes for expensive operations
SMALL_SIZES = [10, 50, 100, 500, 1000]

# Matrix sizes for linear algebra benchmarks
MATRIX_SIZES = [10, 50, 100, 500]

# Standard dtypes for numeric benchmarks
DTYPES = ['float32', 'float64']

# Extended dtypes
DTYPES_EXTENDED = ['int16', 'int32', 'int64', 'float32', 'float64']

# Random seed for reproducible benchmarks
RANDOM_SEED = 42

# Category names for organization
CATEGORY_NAMES = {
    'core': 'Core Operations',
    'reduce': 'Reductions',
    'ufunc': 'Universal Functions',
    'linalg': 'Linear Algebra',
    'random': 'Random',
}
