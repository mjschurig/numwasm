#!/usr/bin/env python3
"""
Extract all exported functions from NumPy reference and NumWasm.
"""
import os
import re
import ast
import json
from pathlib import Path

NUMPY_REF = Path("/workspace/packages/numwasm/reference/numpy/numpy")

def extract_all_from_file(filepath):
    """Extract __all__ list from a Python file."""
    try:
        with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()

        # Try to parse as AST first
        try:
            tree = ast.parse(content)
            all_items = []
            for node in ast.walk(tree):
                if isinstance(node, ast.Assign):
                    for target in node.targets:
                        if isinstance(target, ast.Name) and target.id == '__all__':
                            if isinstance(node.value, ast.List):
                                for elt in node.value.elts:
                                    if isinstance(elt, ast.Constant):
                                        all_items.append(elt.value)
                                    elif isinstance(elt, ast.Str):  # Python < 3.8
                                        all_items.append(elt.s)
                            return all_items
        except:
            pass

        # Fallback to regex
        pattern = r'__all__\s*=\s*\[(.*?)\]'
        match = re.search(pattern, content, re.DOTALL)
        if match:
            items_str = match.group(1)
            items = re.findall(r'["\']([^"\']+)["\']', items_str)
            return items

        return []
    except Exception as e:
        return []

def extract_imports_from_init(filepath):
    """Extract imported names from __init__.py file."""
    try:
        with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()

        imports = set()

        # Parse from ... import ... statements
        pattern = r'from\s+[\w.]+\s+import\s+\((.*?)\)'
        for match in re.finditer(pattern, content, re.DOTALL):
            items = match.group(1)
            names = re.findall(r'\b(\w+)\b', items)
            # Filter out common non-import words
            for name in names:
                if not name.startswith('_') and name not in ['import', 'from', 'as']:
                    imports.add(name)

        # Also handle single-line imports
        pattern2 = r'from\s+[\w.]+\s+import\s+([^(].*?)$'
        for match in re.finditer(pattern2, content, re.MULTILINE):
            items = match.group(1).strip()
            if not items.startswith('*'):
                names = re.findall(r'\b(\w+)\b', items)
                for name in names:
                    if not name.startswith('_') and name not in ['import', 'from', 'as']:
                        imports.add(name)

        return imports
    except:
        return set()

def get_numpy_main_exports():
    """Get main numpy namespace exports."""
    init_file = NUMPY_REF / "__init__.py"

    # Get from imports in __init__.py
    imports = extract_imports_from_init(init_file)

    # Also get __all__
    all_items = extract_all_from_file(init_file)

    return sorted(set(imports) | set(all_items))

def get_submodule_exports(submodule):
    """Get exports from a submodule."""
    init_file = NUMPY_REF / submodule / "__init__.py"
    if not init_file.exists():
        return []

    all_items = extract_all_from_file(init_file)
    if not all_items:
        # Try to get from imports
        imports = extract_imports_from_init(init_file)
        return sorted(imports)
    return all_items

def get_core_exports():
    """Get all exports from numpy._core and related modules."""
    exports = set()

    # Read the core __init__.py
    core_init = NUMPY_REF / "_core" / "__init__.py"
    exports.update(extract_imports_from_init(core_init))
    exports.update(extract_all_from_file(core_init))

    # Read numeric.py __all__
    numeric_file = NUMPY_REF / "_core" / "numeric.py"
    exports.update(extract_all_from_file(numeric_file))

    # Read fromnumeric.py __all__
    fromnumeric_file = NUMPY_REF / "_core" / "fromnumeric.py"
    exports.update(extract_all_from_file(fromnumeric_file))

    # Read function_base.py __all__
    function_base_file = NUMPY_REF / "_core" / "function_base.py"
    exports.update(extract_all_from_file(function_base_file))

    # Read shape_base.py __all__
    shape_base_file = NUMPY_REF / "_core" / "shape_base.py"
    exports.update(extract_all_from_file(shape_base_file))

    # Read getlimits.py __all__
    getlimits_file = NUMPY_REF / "_core" / "getlimits.py"
    exports.update(extract_all_from_file(getlimits_file))

    # Read umath.py __all__
    umath_file = NUMPY_REF / "_core" / "umath.py"
    exports.update(extract_all_from_file(umath_file))

    return sorted(exports)

def get_lib_exports():
    """Get exports from numpy.lib modules."""
    lib_dir = NUMPY_REF / "lib"
    exports = set()

    for pyfile in lib_dir.glob("*_impl.py"):
        exports.update(extract_all_from_file(pyfile))

    return sorted(exports)

def get_linalg_exports():
    """Get numpy.linalg exports."""
    linalg_file = NUMPY_REF / "linalg" / "_linalg.py"
    exports = extract_all_from_file(linalg_file)
    if not exports:
        # Manual list from docstring in __init__.py
        exports = [
            'LinAlgError',
            # Matrix and vector products
            'cross', 'multi_dot', 'matrix_power', 'tensordot', 'matmul', 'outer',
            # Decompositions
            'cholesky', 'qr', 'svd', 'svdvals',
            # Eigenvalues
            'eig', 'eigh', 'eigvals', 'eigvalsh',
            # Norms
            'norm', 'matrix_norm', 'vector_norm', 'cond', 'det', 'matrix_rank', 'slogdet', 'trace',
            # Solving
            'solve', 'tensorsolve', 'lstsq', 'inv', 'pinv', 'tensorinv',
            # Other
            'diagonal', 'matrix_transpose',
        ]
    return exports

def get_fft_exports():
    """Get numpy.fft exports."""
    exports = []
    # From pocketfft
    pocketfft_file = NUMPY_REF / "fft" / "_pocketfft.py"
    exports.extend(extract_all_from_file(pocketfft_file))
    # From helper
    helper_file = NUMPY_REF / "fft" / "_helper.py"
    exports.extend(extract_all_from_file(helper_file))

    if not exports:
        # Manual list
        exports = [
            'fft', 'ifft', 'fft2', 'ifft2', 'fftn', 'ifftn',
            'rfft', 'irfft', 'rfft2', 'irfft2', 'rfftn', 'irfftn',
            'hfft', 'ihfft',
            'fftfreq', 'rfftfreq', 'fftshift', 'ifftshift',
        ]
    return exports

def get_random_exports():
    """Get numpy.random exports."""
    random_init = NUMPY_REF / "random" / "__init__.py"
    return extract_all_from_file(random_init)

def get_ma_exports():
    """Get numpy.ma exports."""
    core_file = NUMPY_REF / "ma" / "core.py"
    extras_file = NUMPY_REF / "ma" / "extras.py"
    exports = []
    exports.extend(extract_all_from_file(core_file))
    exports.extend(extract_all_from_file(extras_file))
    return exports

def get_polynomial_exports():
    """Get numpy.polynomial exports."""
    poly_dir = NUMPY_REF / "polynomial"
    exports = set()

    # Main classes
    exports.update(extract_all_from_file(poly_dir / "__init__.py"))

    # Individual polynomial modules
    for module in ['polynomial', 'chebyshev', 'legendre', 'hermite', 'hermite_e', 'laguerre']:
        mod_file = poly_dir / f"{module}.py"
        if mod_file.exists():
            exports.update(extract_all_from_file(mod_file))

    # Polyutils
    polyutils_file = poly_dir / "polyutils.py"
    if polyutils_file.exists():
        exports.update(extract_all_from_file(polyutils_file))

    return sorted(exports)

def get_strings_exports():
    """Get numpy.strings exports."""
    strings_file = NUMPY_REF / "_core" / "strings.py"
    return extract_all_from_file(strings_file)

def get_testing_exports():
    """Get numpy.testing exports."""
    # Testing utils
    testing_utils = NUMPY_REF / "testing" / "_private" / "utils.py"
    exports = extract_all_from_file(testing_utils)
    if not exports:
        # Manual list
        exports = [
            'assert_equal', 'assert_almost_equal', 'assert_approx_equal',
            'assert_array_equal', 'assert_array_less', 'assert_string_equal',
            'assert_array_almost_equal', 'assert_raises', 'build_err_msg',
            'decorate_methods', 'jiffies', 'memusage', 'print_assert_equal',
            'rundocs', 'runstring', 'verbose', 'measure',
            'assert_', 'assert_array_almost_equal_nulp', 'assert_raises_regex',
            'assert_array_max_ulp', 'assert_warns', 'assert_no_warnings',
            'assert_allclose', 'IgnoreException', 'clear_and_catch_warnings',
            'SkipTest', 'KnownFailureException', 'temppath', 'tempdir', 'IS_PYPY',
            'HAS_REFCOUNT', "IS_WASM", 'suppress_warnings', 'assert_array_compare',
            'assert_no_gc_cycles', 'break_cycles', 'HAS_LAPACK64', 'IS_PYSTON',
            'IS_MUSL', 'check_support_sve', 'NOGIL_BUILD',
            'IS_EDITABLE', 'IS_INSTALLED', 'NUMPY_ROOT', 'run_threaded', 'IS_64BIT',
            'BLAS_SUPPORTS_FPE', 'TestCase', 'overrides',
        ]
    return exports

def main():
    print("Extracting NumPy exports...\n")

    results = {
        "main_namespace": get_numpy_main_exports(),
        "core": get_core_exports(),
        "lib": get_lib_exports(),
        "linalg": get_linalg_exports(),
        "fft": get_fft_exports(),
        "random": get_random_exports(),
        "ma": get_ma_exports(),
        "polynomial": get_polynomial_exports(),
        "strings": get_strings_exports(),
        "testing": get_testing_exports(),
    }

    # Save results
    with open("/workspace/numpy_exports.json", "w") as f:
        json.dump(results, f, indent=2)

    # Print summary
    total = 0
    for category, exports in results.items():
        count = len(exports)
        total += count
        print(f"{category}: {count} exports")

    print(f"\nTotal: {total} exports")
    print("\nResults saved to /workspace/numpy_exports.json")

if __name__ == "__main__":
    main()
