#!/usr/bin/env python3
"""
Master benchmark runner for NumPy.
Runs all benchmark categories and outputs results to JSON.
"""

import sys
import os
import json
from pathlib import Path

# Add the benchmark directory to the path
sys.path.insert(0, str(Path(__file__).parent.parent))

from python.reduce.aggregations import run_reduce_benchmarks
from python.core.creation import run_core_benchmarks
from python.ufunc.arithmetic import run_ufunc_benchmarks
from python.linalg.operations import run_linalg_benchmarks
from python.random_bench.distributions import run_random_benchmarks
from python.lib.config import CATEGORY_NAMES


def main():
    # Parse command line arguments
    selected_category = None
    for arg in sys.argv[1:]:
        if arg.startswith('--category='):
            selected_category = arg.split('=')[1]

    print('=' * 60)
    print('NumPy Benchmark Suite')
    print('=' * 60)

    categories = []

    def run_category(cat_id: str, runner):
        if selected_category and selected_category != cat_id:
            return
        print(f"\n[{CATEGORY_NAMES.get(cat_id, cat_id)}]")
        operations = runner()
        categories.append({
            'id': cat_id,
            'name': CATEGORY_NAMES.get(cat_id, cat_id),
            'operations': operations,
        })

    run_category('reduce', run_reduce_benchmarks)
    run_category('core', run_core_benchmarks)
    run_category('ufunc', run_ufunc_benchmarks)
    run_category('linalg', run_linalg_benchmarks)
    run_category('random', run_random_benchmarks)

    # Output results
    project_root = Path(__file__).parent.parent.parent
    output_dir = project_root / 'benchmark' / 'results' / 'numpy'
    output_dir.mkdir(parents=True, exist_ok=True)

    # Write individual category files
    for category in categories:
        output_file = output_dir / f"{category['id']}.json"
        with open(output_file, 'w') as f:
            json.dump(category, f, indent=2)
        print(f"\nResults saved: {output_file}")

    print('\n' + '=' * 60)
    print('NumPy benchmarks complete!')
    print('=' * 60)


if __name__ == '__main__':
    main()
