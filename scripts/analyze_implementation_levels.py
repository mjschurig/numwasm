#!/usr/bin/env python3
"""
Analyze NumWasm implementation at three levels:
1. WASM Binary - Low-level functions in WebAssembly
2. TypeScript Implementation - Functions implemented in TS
3. TypeScript API - Public exports from numwasm
"""
import json
import re
from pathlib import Path

def load_json(filepath):
    with open(filepath, 'r') as f:
        return json.load(f)

def extract_wasm_exports():
    """Extract all WASM function exports from wasm-exports.json"""
    wasm_file = Path("/workspace/docs-site/public/wasm-exports.json")
    with open(wasm_file, 'r') as f:
        data = json.load(f)

    wasm_funcs = set()

    def extract_from_dict(d):
        for key, value in d.items():
            if isinstance(value, list):
                for func in value:
                    if isinstance(func, str) and func.startswith('_'):
                        wasm_funcs.add(func)
            elif isinstance(value, dict):
                extract_from_dict(value)

    extract_from_dict(data)
    return wasm_funcs

def categorize_wasm_exports(wasm_funcs):
    """Categorize WASM functions by their purpose"""
    categories = {
        'ndarray_core': [],
        'ufuncs': [],
        'linalg': [],
        'fft': [],
        'random': [],
        'special': [],
        'sparse': [],
        'spatial': [],
        'optimize': [],
        'sorting': [],
        'statistics': [],
        'other': [],
    }

    for func in sorted(wasm_funcs):
        if func.startswith('_ufunc_'):
            categories['ufuncs'].append(func)
        elif func.startswith('_ndarray_'):
            categories['ndarray_core'].append(func)
        elif func.startswith('_linalg_') or func.startswith('_lapack_') or 'matrix' in func.lower():
            categories['linalg'].append(func)
        elif func.startswith('_fft_') or 'fft' in func.lower():
            categories['fft'].append(func)
        elif func.startswith('_random_') or func.startswith('_rng_'):
            categories['random'].append(func)
        elif func.startswith('_sp_') or 'sparse' in func.lower():
            categories['sparse'].append(func)
        elif func.startswith('_kdtree_') or 'spatial' in func.lower():
            categories['spatial'].append(func)
        elif 'sort' in func.lower() or 'partition' in func.lower():
            categories['sorting'].append(func)
        elif 'stat' in func.lower() or func.startswith('_histogram'):
            categories['statistics'].append(func)
        elif any(x in func.lower() for x in ['gamma', 'binom', 'poch', 'comb', 'perm']):
            categories['special'].append(func)
        elif any(x in func.lower() for x in ['minimize', 'bfgs', 'nelder']):
            categories['optimize'].append(func)
        else:
            categories['other'].append(func)

    return categories

def extract_ts_exports():
    """Extract TypeScript exports from index.ts"""
    index_file = Path("/workspace/packages/numwasm/src/ts/index.ts")
    with open(index_file, 'r') as f:
        content = f.read()

    # Remove comments
    content = re.sub(r'//[^\n]*', '', content)

    exports = set()

    # Named exports: export { foo, bar } from './module.js';
    pattern = r'export\s*\{([^}]+)\}\s*from'
    for match in re.finditer(pattern, content, re.DOTALL):
        items = match.group(1)
        for item in items.split(','):
            item = item.strip()
            if not item or item.startswith('type '):
                continue
            if ' as ' in item:
                name = item.split(' as ')[1].strip()
            else:
                name = item.strip()
            name = name.split()[0] if name else ''
            if name and not name.startswith('_'):
                exports.add(name)

    # Namespace exports: export * as testing from './testing';
    pattern2 = r'export\s+\*\s+as\s+(\w+)\s+from'
    for match in re.finditer(pattern2, content):
        exports.add(match.group(1))

    return exports

def map_numpy_to_implementation():
    """Map NumPy functions to their implementation status"""
    numpy_exports = load_json("/workspace/numpy_exports.json")
    ts_exports = extract_ts_exports()
    wasm_funcs = extract_wasm_exports()

    # Normalize for comparison
    ts_exports_lower = {e.lower().replace('_', '') for e in ts_exports}

    # Map NumPy functions to their potential WASM counterparts
    def has_wasm_impl(numpy_name):
        """Check if there's a WASM implementation for this function"""
        patterns = [
            f'_ufunc_{numpy_name}',
            f'_ndarray_{numpy_name}',
            f'_linalg_{numpy_name}',
            f'_fft_{numpy_name}',
            f'_random_{numpy_name}',
        ]
        for pattern in patterns:
            for wasm_func in wasm_funcs:
                if pattern.lower() in wasm_func.lower():
                    return True
        return False

    def has_ts_impl(numpy_name):
        """Check if there's a TypeScript export for this function"""
        normalized = numpy_name.lower().replace('_', '')
        return normalized in ts_exports_lower or numpy_name in ts_exports

    results = {}
    for category, funcs in numpy_exports.items():
        results[category] = []
        for func in funcs:
            wasm = has_wasm_impl(func)
            ts = has_ts_impl(func)
            results[category].append({
                'name': func,
                'wasm': wasm,
                'typescript': ts,
            })

    return results, wasm_funcs, ts_exports

def get_reference_sources():
    """Return reference implementation source locations by category."""
    return {
        'ndarray_core': {
            'description': 'Core ndarray operations (creation, manipulation, indexing)',
            'paths': [
                'numpy/_core/src/multiarray/arrayobject.c',
                'numpy/_core/src/multiarray/ctors.c',
                'numpy/_core/src/multiarray/shape.c',
                'numpy/_core/src/multiarray/item_selection.c',
                'numpy/_core/src/multiarray/iterators.c',
                'numpy/_core/src/multiarray/methods.c',
                'numpy/_core/src/multiarray/calculation.c',
                'numpy/_core/src/multiarray/mapping.c',
                'numpy/_core/src/multiarray/convert.c',
                'numpy/_core/src/multiarray/descriptor.c',
            ],
            'headers': [
                'numpy/_core/include/numpy/arrayobject.h',
                'numpy/_core/include/numpy/ndarrayobject.h',
                'numpy/_core/include/numpy/ndarraytypes.h',
            ],
        },
        'ufuncs': {
            'description': 'Universal functions (math operations)',
            'paths': [
                'numpy/_core/src/umath/umathmodule.c',
                'numpy/_core/src/umath/ufunc_object.c',
                'numpy/_core/src/umath/ufunc_type_resolution.c',
                'numpy/_core/src/umath/loops.c.src',
                'numpy/_core/src/umath/loops_trigonometric.dispatch.cpp',
                'numpy/_core/src/umath/loops_logical.dispatch.cpp',
                'numpy/_core/src/umath/reduction.c',
                'numpy/_core/src/umath/clip.cpp',
            ],
            'headers': [
                'numpy/_core/src/umath/ufunc_object.h',
                'numpy/_core/include/numpy/ufuncobject.h',
            ],
            'math_support': [
                'numpy/_core/src/npymath/npy_math.c',
                'numpy/_core/src/npymath/halffloat.cpp',
                'numpy/_core/src/npymath/ieee754.cpp',
            ],
        },
        'linalg': {
            'description': 'Linear algebra (LAPACK/BLAS)',
            'paths': [
                'numpy/linalg/umath_linalg.cpp',
                'numpy/linalg/lapack_litemodule.c',
            ],
            'lapack_lite': [
                'numpy/linalg/lapack_lite/f2c.c',
                'numpy/linalg/lapack_lite/f2c_blas.c',
                'numpy/linalg/lapack_lite/f2c_lapack.c',
                'numpy/linalg/lapack_lite/f2c_d_lapack.c',
                'numpy/linalg/lapack_lite/f2c_z_lapack.c',
                'numpy/linalg/lapack_lite/python_xerbla.c',
            ],
        },
        'fft': {
            'description': 'Fast Fourier Transform',
            'paths': [
                'numpy/fft/_pocketfft_umath.cpp',
            ],
        },
        'random': {
            'description': 'Random number generation',
            'paths': [
                'numpy/random/src/distributions/distributions.c',
                'numpy/random/src/distributions/logfactorial.c',
                'numpy/random/src/distributions/random_hypergeometric.c',
                'numpy/random/src/legacy/legacy-distributions.c',
            ],
            'generators': {
                'MT19937': [
                    'numpy/random/src/mt19937/mt19937.c',
                    'numpy/random/src/mt19937/randomkit.c',
                ],
                'PCG64': [
                    'numpy/random/src/pcg64/pcg64.c',
                ],
                'Philox': [
                    'numpy/random/src/philox/philox.c',
                ],
                'SFC64': [
                    'numpy/random/src/sfc64/sfc64.c',
                ],
            },
        },
        'sorting': {
            'description': 'Sorting and searching algorithms',
            'paths': [
                'numpy/_core/src/npysort/quicksort.cpp',
                'numpy/_core/src/npysort/mergesort.cpp',
                'numpy/_core/src/npysort/heapsort.cpp',
                'numpy/_core/src/npysort/timsort.cpp',
                'numpy/_core/src/npysort/radixsort.cpp',
                'numpy/_core/src/npysort/binsearch.cpp',
                'numpy/_core/src/npysort/selection.cpp',
            ],
        },
        'strings': {
            'description': 'String operations',
            'paths': [
                'numpy/_core/src/multiarray/strfuncs.c',
                'numpy/_core/src/multiarray/stringdtype/dtype.c',
                'numpy/_core/src/multiarray/stringdtype/casts.cpp',
                'numpy/_core/src/multiarray/stringdtype/static_string.c',
                'numpy/_core/src/multiarray/stringdtype/utf8_utils.c',
                'numpy/_core/src/umath/string_ufuncs.cpp',
            ],
        },
        'sparse': {
            'description': 'Sparse matrix operations (SciPy reference)',
            'note': 'Based on scipy.sparse implementation',
            'paths': [
                'scipy/sparse/sparsetools/*.cpp (external)',
            ],
        },
        'spatial': {
            'description': 'Spatial algorithms (KDTree)',
            'note': 'Based on scipy.spatial implementation',
            'paths': [
                'scipy/spatial/ckdtree/*.cxx (external)',
            ],
        },
        'optimize': {
            'description': 'Optimization algorithms',
            'note': 'Based on scipy.optimize implementation',
            'paths': [
                'scipy/optimize/*.c (external)',
            ],
        },
        'special': {
            'description': 'Special mathematical functions',
            'note': 'Gamma, combinatorics, etc.',
            'paths': [
                'scipy/special/*.c (external)',
            ],
        },
    }


def generate_markdown():
    """Generate comprehensive markdown with both levels"""
    results, wasm_funcs, ts_exports = map_numpy_to_implementation()
    wasm_categories = categorize_wasm_exports(wasm_funcs)
    ref_sources = get_reference_sources()

    md = []
    md.append("# NumWasm Implementation Status")
    md.append("")
    md.append("This document tracks implementation status at two levels:")
    md.append("")
    md.append("| Symbol | Meaning |")
    md.append("|--------|---------|")
    md.append("| ✅ | Has TypeScript API (exported from numwasm) |")
    md.append("| 🔧 | Has WASM implementation (low-level) |")
    md.append("| ⬜ | Not implemented |")
    md.append("")

    # Summary
    total_numpy = sum(len(funcs) for funcs in results.values())
    total_ts = sum(1 for funcs in results.values() for f in funcs if f['typescript'])
    total_wasm = sum(1 for funcs in results.values() for f in funcs if f['wasm'])

    md.append("## Summary")
    md.append("")
    md.append(f"- **Total NumPy functions**: {total_numpy}")
    md.append(f"- **With TypeScript API**: {total_ts} ({total_ts/total_numpy*100:.1f}%)")
    md.append(f"- **With WASM backing**: {total_wasm} ({total_wasm/total_numpy*100:.1f}%)")
    md.append(f"- **Total WASM exports**: {len(wasm_funcs)}")
    md.append(f"- **Total TypeScript exports**: {len(ts_exports)}")
    md.append("")

    # WASM exports summary
    md.append("## WASM Binary Exports by Category")
    md.append("")
    md.append("| Category | Count |")
    md.append("|----------|-------|")
    for cat, funcs in sorted(wasm_categories.items(), key=lambda x: -len(x[1])):
        if funcs:
            md.append(f"| {cat} | {len(funcs)} |")
    md.append("")

    # Detailed WASM exports with reference sources
    md.append("### WASM Categories with Reference Implementation Sources")
    md.append("")
    md.append("Each category shows the C/C++ source files in the NumPy/SciPy reference implementation.")
    md.append("")

    for cat, funcs in sorted(wasm_categories.items()):
        if funcs:
            md.append(f"<details>")
            md.append(f"<summary><b>{cat}</b> ({len(funcs)} WASM functions)</summary>")
            md.append("")

            # Add reference source info if available
            if cat in ref_sources:
                ref = ref_sources[cat]
                md.append(f"**Description:** {ref.get('description', 'N/A')}")
                md.append("")

                if 'note' in ref:
                    md.append(f"*Note: {ref['note']}*")
                    md.append("")

                md.append("**Reference C/C++ Sources:**")
                md.append("```")
                for path in ref.get('paths', []):
                    md.append(path)
                md.append("```")
                md.append("")

                if 'headers' in ref:
                    md.append("**Headers:**")
                    md.append("```")
                    for path in ref['headers']:
                        md.append(path)
                    md.append("```")
                    md.append("")

                if 'lapack_lite' in ref:
                    md.append("**LAPACK Lite:**")
                    md.append("```")
                    for path in ref['lapack_lite']:
                        md.append(path)
                    md.append("```")
                    md.append("")

                if 'generators' in ref:
                    md.append("**RNG Generators:**")
                    for gen_name, gen_paths in ref['generators'].items():
                        md.append(f"- {gen_name}:")
                        for path in gen_paths:
                            md.append(f"  - `{path}`")
                    md.append("")

                if 'math_support' in ref:
                    md.append("**Math Support:**")
                    md.append("```")
                    for path in ref['math_support']:
                        md.append(path)
                    md.append("```")
                    md.append("")

            md.append("**WASM Exports:**")
            md.append("```")
            for func in sorted(funcs):
                md.append(func)
            md.append("```")
            md.append("</details>")
            md.append("")

    # NumPy function status by category
    md.append("---")
    md.append("")
    md.append("## NumPy Function Implementation Status")
    md.append("")

    category_order = ['main_namespace', 'linalg', 'fft', 'random', 'ma', 'polynomial', 'strings', 'testing']
    category_names = {
        'main_namespace': 'numpy (Main Namespace)',
        'linalg': 'numpy.linalg',
        'fft': 'numpy.fft',
        'random': 'numpy.random',
        'ma': 'numpy.ma (Masked Arrays)',
        'polynomial': 'numpy.polynomial',
        'strings': 'numpy.strings',
        'testing': 'numpy.testing',
    }

    for cat in category_order:
        if cat not in results:
            continue
        funcs = results[cat]

        ts_count = sum(1 for f in funcs if f['typescript'])
        wasm_count = sum(1 for f in funcs if f['wasm'])

        md.append(f"### {category_names.get(cat, cat)}")
        md.append("")
        md.append(f"**TypeScript API**: {ts_count}/{len(funcs)} | **WASM backing**: {wasm_count}/{len(funcs)}")
        md.append("")
        md.append("| Function | TS API | WASM |")
        md.append("|----------|--------|------|")

        for f in sorted(funcs, key=lambda x: x['name']):
            ts_mark = "✅" if f['typescript'] else "⬜"
            wasm_mark = "🔧" if f['wasm'] else "⬜"
            md.append(f"| `{f['name']}` | {ts_mark} | {wasm_mark} |")
        md.append("")

    # TypeScript-only exports (not in NumPy)
    numpy_all = set()
    for funcs in results.values():
        numpy_all.update(f['name'].lower() for f in funcs)

    ts_only = [e for e in sorted(ts_exports) if e.lower() not in numpy_all]
    if ts_only:
        md.append("---")
        md.append("")
        md.append("## NumWasm-Specific Exports (not in NumPy)")
        md.append("")
        md.append("These are TypeScript exports unique to NumWasm:")
        md.append("")
        for name in ts_only:
            md.append(f"- `{name}`")
        md.append("")

    return '\n'.join(md)

def main():
    md = generate_markdown()
    output_path = Path("/workspace/packages/numwasm/numwasm-todo.md")
    with open(output_path, 'w') as f:
        f.write(md)
    print(f"Generated: {output_path}")

if __name__ == "__main__":
    main()
