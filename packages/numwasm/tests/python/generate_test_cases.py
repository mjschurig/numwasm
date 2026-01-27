#!/usr/bin/env python3
"""
Generate test vectors using NumPy for comparison with TypeScript/WASM implementation.

This script creates a JSON file containing test cases with input data and expected
results computed by NumPy. The TypeScript tests will load these vectors and verify
that the WASM implementation produces matching results.

Based on NumPy's own test patterns from numpy/_core/tests/
"""

import json
import sys
from pathlib import Path

try:
    import numpy as np
except ImportError:
    print("Error: NumPy is not installed. Run: pip install numpy", file=sys.stderr)
    sys.exit(1)


def generate_sum_tests():
    """Generate test cases for array sum operations."""
    test_cases = []
    np.random.seed(42)

    # Test case 1: Small array sum
    arr = np.array([1.0, 2.0, 3.0, 4.0, 5.0], dtype=np.float64)
    test_cases.append({
        "name": "small_array_sum",
        "description": "Simple 5-element array",
        "data": arr.tolist(),
        "shape": list(arr.shape),
        "dtype": "float64",
        "expected_sum": float(arr.sum()),
    })

    # Test case 2: Medium array (tests pairwise algorithm with 8-way unrolling)
    arr = np.arange(1000, dtype=np.float64)
    test_cases.append({
        "name": "medium_array_sum",
        "description": "1000 elements (0-999), tests 8-way unrolled loop",
        "data": arr.tolist(),
        "shape": list(arr.shape),
        "dtype": "float64",
        "expected_sum": float(arr.sum()),
    })

    # Test case 3: Large array (tests recursive divide-and-conquer)
    arr = np.random.randn(10000).astype(np.float64)
    test_cases.append({
        "name": "large_random_array_sum",
        "description": "10000 random elements, tests recursive pairwise",
        "data": arr.tolist(),
        "shape": list(arr.shape),
        "dtype": "float64",
        "expected_sum": float(arr.sum()),
    })

    # Test case 4: Precision test - many small numbers
    arr = np.full(10000, 0.1, dtype=np.float64)
    test_cases.append({
        "name": "precision_test",
        "description": "10000 copies of 0.1 - tests precision of pairwise sum",
        "data": arr.tolist(),
        "shape": list(arr.shape),
        "dtype": "float64",
        "expected_sum": float(arr.sum()),
    })

    # Test case 5: Array with negative values
    arr = np.array([-1.5, 2.5, -3.5, 4.5, -5.5], dtype=np.float64)
    test_cases.append({
        "name": "negative_values",
        "description": "Mixed positive and negative values",
        "data": arr.tolist(),
        "shape": list(arr.shape),
        "dtype": "float64",
        "expected_sum": float(arr.sum()),
    })

    # Test case 6: 2D array (flattened for sum)
    arr = np.array([[1, 2, 3], [4, 5, 6]], dtype=np.float64)
    test_cases.append({
        "name": "2d_array_sum",
        "description": "2x3 matrix, sum of all elements",
        "data": arr.flatten().tolist(),
        "shape": list(arr.shape),
        "dtype": "float64",
        "expected_sum": float(arr.sum()),
    })

    # Test case 7: Float32 array
    arr = np.array([1.0, 2.0, 3.0, 4.0, 5.0], dtype=np.float32)
    test_cases.append({
        "name": "float32_array_sum",
        "description": "Float32 array to test single precision",
        "data": [float(x) for x in arr],
        "shape": list(arr.shape),
        "dtype": "float32",
        "expected_sum": float(arr.sum()),
    })

    # Test case 8: Int32 array
    arr = np.array([1, 2, 3, 4, 5], dtype=np.int32)
    test_cases.append({
        "name": "int32_array_sum",
        "description": "Integer array",
        "data": [int(x) for x in arr],
        "shape": list(arr.shape),
        "dtype": "int32",
        "expected_sum": int(arr.sum()),
    })

    # Edge cases
    test_cases.append({
        "name": "single_element",
        "description": "Single element array",
        "data": [42.0],
        "shape": [1],
        "dtype": "float64",
        "expected_sum": 42.0,
    })

    test_cases.append({
        "name": "all_zeros",
        "description": "100 zeros",
        "data": [0.0] * 100,
        "shape": [100],
        "dtype": "float64",
        "expected_sum": 0.0,
    })

    # Boundary cases for pairwise algorithm
    for n in [3, 7, 8, 9, 127, 128, 129, 256]:
        arr = np.arange(n, dtype=np.float64)
        test_cases.append({
            "name": f"boundary_{n}_elements",
            "description": f"{n} elements - pairwise algorithm boundary test",
            "data": arr.tolist(),
            "shape": [n],
            "dtype": "float64",
            "expected_sum": float(arr.sum()),
        })

    return test_cases


def generate_dtype_tests():
    """Generate test cases for dtype operations based on NumPy test_dtype.py patterns."""
    test_cases = []

    # Test dtype sizes (from test_dtype.py)
    dtype_sizes = {
        "bool": 1,
        "int8": 1,
        "int16": 2,
        "int32": 4,
        "int64": 8,
        "uint8": 1,
        "uint16": 2,
        "uint32": 4,
        "uint64": 8,
        "float16": 2,
        "float32": 4,
        "float64": 8,
        "complex64": 8,
        "complex128": 16,
    }

    for dtype_name, expected_size in dtype_sizes.items():
        try:
            dt = np.dtype(dtype_name)
            test_cases.append({
                "name": f"dtype_size_{dtype_name}",
                "dtype": dtype_name,
                "expected_size": dt.itemsize,
                "expected_kind": dt.kind,
            })
        except TypeError:
            pass  # Skip unsupported dtypes

    # Test type promotion (based on test_nep50_promotions.py)
    promotion_tests = [
        ("int32", "float32", "float64"),  # int32 + float32 -> float64
        ("int16", "float32", "float32"),  # int16 + float32 -> float32
        ("int64", "float32", "float64"),  # int64 + float32 -> float64
        ("int32", "int32", "int32"),
        ("float32", "float64", "float64"),
        ("int8", "int16", "int16"),
        ("uint8", "int8", "int16"),
        ("uint16", "int16", "int32"),
        ("uint32", "int32", "int64"),
        ("bool", "int32", "int32"),
        ("bool", "float64", "float64"),
    ]

    for dt1, dt2, expected in promotion_tests:
        try:
            result = np.result_type(np.dtype(dt1), np.dtype(dt2))
            test_cases.append({
                "name": f"promote_{dt1}_{dt2}",
                "dtype1": dt1,
                "dtype2": dt2,
                "expected_result": str(result),
            })
        except TypeError:
            pass

    return test_cases


def generate_creation_tests():
    """Generate test cases for array creation functions based on NumPy test_multiarray.py."""
    test_cases = []

    # zeros
    for shape in [[5], [2, 3], [2, 3, 4]]:
        arr = np.zeros(shape, dtype=np.float64)
        test_cases.append({
            "name": f"zeros_{shape}",
            "function": "zeros",
            "shape": shape,
            "dtype": "float64",
            "expected_data": arr.flatten().tolist(),
            "expected_sum": 0.0,
        })

    # ones
    for shape in [[5], [2, 3], [2, 3, 4]]:
        arr = np.ones(shape, dtype=np.float64)
        test_cases.append({
            "name": f"ones_{shape}",
            "function": "ones",
            "shape": shape,
            "dtype": "float64",
            "expected_data": arr.flatten().tolist(),
            "expected_sum": float(arr.sum()),
        })

    # full
    for fill_value in [0, 1, 42, -3.14, 1e10]:
        arr = np.full([3, 4], fill_value, dtype=np.float64)
        test_cases.append({
            "name": f"full_{fill_value}",
            "function": "full",
            "shape": [3, 4],
            "fill_value": fill_value,
            "dtype": "float64",
            "expected_data": arr.flatten().tolist(),
        })

    # arange
    arange_tests = [
        (0, 10, 1),
        (0, 10, 2),
        (5, 15, 1),
        (0.0, 1.0, 0.1),
        (-5, 5, 1),
    ]
    for start, stop, step in arange_tests:
        arr = np.arange(start, stop, step, dtype=np.float64)
        test_cases.append({
            "name": f"arange_{start}_{stop}_{step}",
            "function": "arange",
            "start": start,
            "stop": stop,
            "step": step,
            "dtype": "float64",
            "expected_data": arr.tolist(),
            "expected_shape": list(arr.shape),
        })

    # linspace
    linspace_tests = [
        (0, 10, 5, True),
        (0, 10, 5, False),
        (0, 1, 50, True),
        (-5, 5, 11, True),
        (0, 10, 1, True),
    ]
    for start, stop, num, endpoint in linspace_tests:
        arr = np.linspace(start, stop, num, endpoint=endpoint, dtype=np.float64)
        test_cases.append({
            "name": f"linspace_{start}_{stop}_{num}_{endpoint}",
            "function": "linspace",
            "start": start,
            "stop": stop,
            "num": num,
            "endpoint": endpoint,
            "dtype": "float64",
            "expected_data": arr.tolist(),
        })

    # logspace
    logspace_tests = [
        (0, 2, 3, True, 10),
        (0, 3, 4, True, 2),
        (1, 3, 5, True, 10),
    ]
    for start, stop, num, endpoint, base in logspace_tests:
        arr = np.logspace(start, stop, num, endpoint=endpoint, base=base, dtype=np.float64)
        test_cases.append({
            "name": f"logspace_{start}_{stop}_{num}_{base}",
            "function": "logspace",
            "start": start,
            "stop": stop,
            "num": num,
            "endpoint": endpoint,
            "base": base,
            "dtype": "float64",
            "expected_data": arr.tolist(),
        })

    # geomspace
    geomspace_tests = [
        (1, 1000, 4, True),
        (1, 8, 4, True),
        (1, 256, 9, True),
    ]
    for start, stop, num, endpoint in geomspace_tests:
        arr = np.geomspace(start, stop, num, endpoint=endpoint, dtype=np.float64)
        test_cases.append({
            "name": f"geomspace_{start}_{stop}_{num}",
            "function": "geomspace",
            "start": start,
            "stop": stop,
            "num": num,
            "endpoint": endpoint,
            "dtype": "float64",
            "expected_data": arr.tolist(),
        })

    # eye
    eye_tests = [
        (3, None, 0),
        (3, 4, 0),
        (3, 4, 1),
        (4, 3, -1),
        (5, 5, 2),
    ]
    for N, M, k in eye_tests:
        M_val = M if M is not None else N
        arr = np.eye(N, M_val, k, dtype=np.float64)
        test_cases.append({
            "name": f"eye_{N}_{M_val}_{k}",
            "function": "eye",
            "N": N,
            "M": M_val,
            "k": k,
            "dtype": "float64",
            "expected_data": arr.flatten().tolist(),
            "expected_shape": list(arr.shape),
        })

    # identity
    for n in [1, 2, 3, 5]:
        arr = np.identity(n, dtype=np.float64)
        test_cases.append({
            "name": f"identity_{n}",
            "function": "identity",
            "n": n,
            "dtype": "float64",
            "expected_data": arr.flatten().tolist(),
        })

    # diag - create diagonal matrix from 1D
    for k in [-1, 0, 1, 2]:
        v = np.array([1, 2, 3], dtype=np.float64)
        arr = np.diag(v, k)
        test_cases.append({
            "name": f"diag_create_{k}",
            "function": "diag_create",
            "input": v.tolist(),
            "k": k,
            "dtype": "float64",
            "expected_data": arr.flatten().tolist(),
            "expected_shape": list(arr.shape),
        })

    # diag - extract diagonal from 2D
    mat = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]], dtype=np.float64)
    for k in [-1, 0, 1]:
        arr = np.diag(mat, k)
        test_cases.append({
            "name": f"diag_extract_{k}",
            "function": "diag_extract",
            "input": mat.flatten().tolist(),
            "input_shape": list(mat.shape),
            "k": k,
            "dtype": "float64",
            "expected_data": arr.tolist(),
        })

    # tri
    for N, M, k in [(3, 3, 0), (3, 4, 0), (4, 3, 1), (3, 3, -1)]:
        arr = np.tri(N, M, k, dtype=np.float64)
        test_cases.append({
            "name": f"tri_{N}_{M}_{k}",
            "function": "tri",
            "N": N,
            "M": M,
            "k": k,
            "dtype": "float64",
            "expected_data": arr.flatten().tolist(),
            "expected_shape": list(arr.shape),
        })

    # tril/triu
    mat = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]], dtype=np.float64)
    for k in [-1, 0, 1]:
        tril = np.tril(mat, k)
        triu = np.triu(mat, k)
        test_cases.append({
            "name": f"tril_{k}",
            "function": "tril",
            "input": mat.flatten().tolist(),
            "input_shape": list(mat.shape),
            "k": k,
            "expected_data": tril.flatten().tolist(),
        })
        test_cases.append({
            "name": f"triu_{k}",
            "function": "triu",
            "input": mat.flatten().tolist(),
            "input_shape": list(mat.shape),
            "k": k,
            "expected_data": triu.flatten().tolist(),
        })

    return test_cases


def generate_element_access_tests():
    """Generate test cases for element access based on NumPy test_indexing.py."""
    test_cases = []

    # 1D indexing
    arr = np.array([10, 20, 30, 40, 50], dtype=np.float64)
    for idx in [0, 2, 4, -1, -3]:
        test_cases.append({
            "name": f"get_1d_{idx}",
            "function": "get",
            "data": arr.tolist(),
            "shape": list(arr.shape),
            "indices": [idx if idx >= 0 else len(arr) + idx],
            "expected_value": float(arr[idx]),
        })

    # 2D indexing
    arr = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]], dtype=np.float64)
    indices_2d = [(0, 0), (0, 2), (1, 1), (2, 0), (2, 2)]
    for i, j in indices_2d:
        test_cases.append({
            "name": f"get_2d_{i}_{j}",
            "function": "get",
            "data": arr.flatten().tolist(),
            "shape": list(arr.shape),
            "indices": [i, j],
            "expected_value": float(arr[i, j]),
        })

    # 3D indexing
    arr = np.arange(24, dtype=np.float64).reshape(2, 3, 4)
    indices_3d = [(0, 0, 0), (0, 1, 2), (1, 2, 3), (1, 0, 0)]
    for i, j, k in indices_3d:
        test_cases.append({
            "name": f"get_3d_{i}_{j}_{k}",
            "function": "get",
            "data": arr.flatten().tolist(),
            "shape": list(arr.shape),
            "indices": [i, j, k],
            "expected_value": float(arr[i, j, k]),
        })

    # Flat indexing
    arr = np.array([[1, 2, 3], [4, 5, 6]], dtype=np.float64)
    for flat_idx in [0, 1, 3, 5]:
        test_cases.append({
            "name": f"get_flat_{flat_idx}",
            "function": "get_flat",
            "data": arr.flatten().tolist(),
            "shape": list(arr.shape),
            "flat_index": flat_idx,
            "expected_value": float(arr.flat[flat_idx]),
        })

    return test_cases


def generate_shape_tests():
    """Generate test cases for shape manipulation based on NumPy test_shape_base.py."""
    test_cases = []

    # reshape
    arr = np.arange(12, dtype=np.float64)
    reshape_tests = [
        ([12], [3, 4]),
        ([12], [2, 6]),
        ([12], [4, 3]),
        ([12], [2, 2, 3]),
        ([12], [1, 12]),
        ([12], [12, 1]),
    ]
    for original_shape, new_shape in reshape_tests:
        test_cases.append({
            "name": f"reshape_{original_shape}_to_{new_shape}",
            "function": "reshape",
            "data": arr.tolist(),
            "original_shape": original_shape,
            "new_shape": new_shape,
            "expected_shape": new_shape,
        })

    # reshape with -1
    test_cases.append({
        "name": "reshape_infer_dim",
        "function": "reshape",
        "data": np.arange(12).tolist(),
        "original_shape": [12],
        "new_shape": [3, -1],
        "expected_shape": [3, 4],
    })

    # transpose
    arr = np.arange(6, dtype=np.float64).reshape(2, 3)
    transposed = arr.T
    test_cases.append({
        "name": "transpose_2d",
        "function": "transpose",
        "data": arr.flatten().tolist(),
        "shape": list(arr.shape),
        "expected_data": transposed.flatten().tolist(),
        "expected_shape": list(transposed.shape),
    })

    # 3D transpose
    arr = np.arange(24, dtype=np.float64).reshape(2, 3, 4)
    transposed = arr.transpose(2, 0, 1)
    test_cases.append({
        "name": "transpose_3d_axes",
        "function": "transpose",
        "data": arr.flatten().tolist(),
        "shape": list(arr.shape),
        "axes": [2, 0, 1],
        "expected_data": transposed.flatten().tolist(),
        "expected_shape": list(transposed.shape),
    })

    # squeeze
    arr = np.array([[[1, 2, 3]]]).astype(np.float64)  # shape (1, 1, 3)
    squeezed = arr.squeeze()
    test_cases.append({
        "name": "squeeze_all",
        "function": "squeeze",
        "data": arr.flatten().tolist(),
        "shape": list(arr.shape),
        "expected_data": squeezed.flatten().tolist(),
        "expected_shape": list(squeezed.shape),
    })

    # expand_dims
    arr = np.array([1, 2, 3], dtype=np.float64)
    for axis in [0, 1, -1]:
        expanded = np.expand_dims(arr, axis)
        test_cases.append({
            "name": f"expand_dims_axis_{axis}",
            "function": "expand_dims",
            "data": arr.tolist(),
            "shape": list(arr.shape),
            "axis": axis,
            "expected_shape": list(expanded.shape),
        })

    # ravel
    arr = np.array([[1, 2, 3], [4, 5, 6]], dtype=np.float64)
    test_cases.append({
        "name": "ravel",
        "function": "ravel",
        "data": arr.flatten().tolist(),
        "shape": list(arr.shape),
        "expected_data": arr.ravel().tolist(),
        "expected_shape": list(arr.ravel().shape),
    })

    # flatten
    test_cases.append({
        "name": "flatten",
        "function": "flatten",
        "data": arr.flatten().tolist(),
        "shape": list(arr.shape),
        "expected_data": arr.flatten().tolist(),
        "expected_shape": list(arr.flatten().shape),
    })

    # swapaxes
    arr = np.arange(24, dtype=np.float64).reshape(2, 3, 4)
    swapped = arr.swapaxes(0, 2)
    test_cases.append({
        "name": "swapaxes_0_2",
        "function": "swapaxes",
        "data": arr.flatten().tolist(),
        "shape": list(arr.shape),
        "axis1": 0,
        "axis2": 2,
        "expected_shape": list(swapped.shape),
    })

    return test_cases


def generate_copy_astype_tests():
    """Generate test cases for copy and type conversion."""
    test_cases = []

    # copy
    arr = np.array([1, 2, 3, 4, 5], dtype=np.float64)
    test_cases.append({
        "name": "copy_1d",
        "function": "copy",
        "data": arr.tolist(),
        "shape": list(arr.shape),
        "dtype": "float64",
        "expected_data": arr.tolist(),
    })

    # astype conversions
    conversions = [
        ("float64", "float32"),
        ("float64", "int32"),
        ("int32", "float64"),
        ("float32", "int32"),
        ("int32", "int16"),
        ("int16", "int32"),
    ]

    for from_dtype, to_dtype in conversions:
        arr = np.array([1.5, 2.7, 3.9, -1.1, -2.9], dtype=getattr(np, from_dtype))
        converted = arr.astype(getattr(np, to_dtype))
        test_cases.append({
            "name": f"astype_{from_dtype}_to_{to_dtype}",
            "function": "astype",
            "data": [float(x) for x in arr],
            "from_dtype": from_dtype,
            "to_dtype": to_dtype,
            "expected_data": [float(x) if 'float' in to_dtype else int(x) for x in converted],
        })

    return test_cases


def main():
    """Generate test vectors and write to JSON file."""
    script_dir = Path(__file__).parent
    output_dir = script_dir.parent / "fixtures"
    output_dir.mkdir(exist_ok=True)

    # Generate all test categories
    all_tests = {
        "sum": generate_sum_tests(),
        "dtype": generate_dtype_tests(),
        "creation": generate_creation_tests(),
        "element_access": generate_element_access_tests(),
        "shape": generate_shape_tests(),
        "copy_astype": generate_copy_astype_tests(),
    }

    total_count = sum(len(tests) for tests in all_tests.values())

    output_file = output_dir / "test_vectors.json"
    with open(output_file, "w") as f:
        json.dump({
            "generated_by": "NumPy",
            "numpy_version": np.__version__,
            "total_test_count": total_count,
            "categories": {
                name: {
                    "count": len(tests),
                    "test_cases": tests
                }
                for name, tests in all_tests.items()
            },
            # Keep backward compatibility with old format
            "test_count": len(all_tests["sum"]),
            "test_cases": all_tests["sum"],
        }, f, indent=2)

    print(f"Generated {total_count} test cases across {len(all_tests)} categories")
    for name, tests in all_tests.items():
        print(f"  {name}: {len(tests)} tests")
    print(f"NumPy version: {np.__version__}")
    print(f"Output file: {output_file}")


if __name__ == "__main__":
    main()
