#!/usr/bin/env python3
"""
Generate complete TODO with implementation status.
"""
import json

def load_json(filepath):
    with open(filepath, 'r') as f:
        return json.load(f)

def normalize_name(name):
    """Normalize to lowercase for comparison."""
    return name.lower().replace('_', '')

def main():
    numpy_exports = load_json("/workspace/numpy_exports.json")
    numwasm_exports = load_json("/workspace/numwasm_exports.json")

    # Build NumWasm lookup set (normalized names)
    numwasm_all = set()
    numwasm_original = {}  # Map normalized -> original name

    for category, items in numwasm_exports.items():
        for item in items:
            item = item.replace(' (namespace)', '')
            normalized = normalize_name(item)
            numwasm_all.add(normalized)
            numwasm_original[normalized] = item
            # Also add without underscores
            numwasm_all.add(item.lower())
            numwasm_original[item.lower()] = item

    # Known mappings between NumPy and NumWasm names
    name_mappings = {
        'apply_along_axis': 'applyAlongAxis',
        'apply_over_axes': 'applyOverAxes',
        'array_repr': 'arrayRepr',
        'array_str': 'arrayStr',
        'broadcast_arrays': 'broadcastArrays',
        'broadcast_shapes': 'broadcastShapes',
        'broadcast_to': 'broadcastTo',
        'count_nonzero': 'countNonzero',
        'delete': 'deleteArr',
        'var': 'var_',
        'common_type': 'commonType',
        'promote_types': 'promoteTypes',
        'ravel_multi_index': 'ravelMultiIndex',
        'unravel_index': 'unravelIndex',
        'unique_all': 'uniqueAll',
        'unique_counts': 'uniqueCounts',
        'unique_inverse': 'uniqueInverse',
        'unique_values': 'uniqueValues',
        'set_printoptions': 'setPrintoptions',
        'get_printoptions': 'getPrintoptions',
        'format_float_positional': 'formatFloatPositional',
        'format_float_scientific': 'formatFloatScientific',
        'binary_repr': 'binaryRepr',
        'base_repr': 'baseRepr',
        'amax': 'max',
        'amin': 'min',
        'around': 'round',
        'concatenate': 'concatenate',
        'concat': 'concatenate',
    }

    def is_implemented(name):
        """Check if a NumPy function is implemented in NumWasm."""
        normalized = normalize_name(name)
        if normalized in numwasm_all:
            return True
        if name.lower() in numwasm_all:
            return True
        if name in name_mappings:
            mapped = name_mappings[name]
            if normalize_name(mapped) in numwasm_all:
                return True
        return False

    def get_numwasm_name(numpy_name):
        """Get the NumWasm equivalent name."""
        if numpy_name in name_mappings:
            return name_mappings[numpy_name]
        normalized = normalize_name(numpy_name)
        if normalized in numwasm_original:
            return numwasm_original[normalized]
        if numpy_name.lower() in numwasm_original:
            return numwasm_original[numpy_name.lower()]
        return None

    # Generate markdown
    md = []
    md.append("# NumWasm Implementation Status")
    md.append("")
    md.append("Complete list of NumPy functions with implementation status in NumWasm.")
    md.append("")
    md.append("**Legend:**")
    md.append("- [x] Implemented in NumWasm")
    md.append("- [ ] Not yet implemented")
    md.append("- NumWasm name shown in parentheses if different from NumPy")
    md.append("")

    # Count stats
    total_numpy = 0
    total_implemented = 0

    # ============================================
    # MAIN NAMESPACE
    # ============================================
    md.append("---")
    md.append("")
    md.append("## numpy (Main Namespace)")
    md.append("")

    # Categorize main namespace
    main = numpy_exports['main_namespace']

    # Group by category
    categories = {
        'Array Creation': ['array', 'zeros', 'ones', 'empty', 'full', 'arange', 'linspace', 'logspace',
                          'geomspace', 'eye', 'identity', 'diag', 'diagflat', 'tri', 'tril', 'triu',
                          'vander', 'zeros_like', 'ones_like', 'empty_like', 'full_like', 'fromfunction',
                          'fromiter', 'fromstring', 'frombuffer', 'fromfile', 'fromregex', 'copy',
                          'asarray', 'asanyarray', 'ascontiguousarray', 'asfortranarray', 'asarray_chkfinite',
                          'require', 'from_dlpack'],
        'Array Manipulation': ['reshape', 'ravel', 'flatten', 'squeeze', 'expand_dims', 'transpose',
                              'swapaxes', 'moveaxis', 'rollaxis', 'atleast_1d', 'atleast_2d', 'atleast_3d',
                              'concatenate', 'concat', 'stack', 'vstack', 'hstack', 'dstack', 'column_stack',
                              'row_stack', 'split', 'array_split', 'vsplit', 'hsplit', 'dsplit', 'unstack',
                              'tile', 'repeat', 'pad', 'flip', 'fliplr', 'flipud', 'roll', 'rot90', 'resize',
                              'trim_zeros', 'insert', 'delete', 'append', 'copyto', 'block'],
        'Indexing': ['take', 'put', 'compress', 'extract', 'choose', 'diagonal', 'select', 'where',
                    'nonzero', 'flatnonzero', 'argwhere', 'searchsorted', 'indices', 'ix_', 'meshgrid',
                    'diag_indices', 'diag_indices_from', 'tril_indices', 'triu_indices', 'tril_indices_from',
                    'triu_indices_from', 'mask_indices', 'take_along_axis', 'put_along_axis', 'putmask',
                    'place', 'fill_diagonal', 'c_', 'r_', 's_', 'index_exp', 'mgrid', 'ogrid',
                    'ravel_multi_index', 'unravel_index', 'ndenumerate', 'ndindex'],
        'Math - Basic': ['add', 'subtract', 'multiply', 'divide', 'true_divide', 'floor_divide',
                        'negative', 'positive', 'power', 'pow', 'remainder', 'mod', 'fmod', 'divmod',
                        'absolute', 'abs', 'fabs', 'sign', 'reciprocal', 'sqrt', 'square', 'cbrt'],
        'Math - Exponential': ['exp', 'exp2', 'expm1', 'log', 'log2', 'log10', 'log1p', 'logaddexp', 'logaddexp2'],
        'Math - Trigonometric': ['sin', 'cos', 'tan', 'arcsin', 'arccos', 'arctan', 'arctan2',
                                'asin', 'acos', 'atan', 'atan2', 'hypot', 'degrees', 'radians',
                                'deg2rad', 'rad2deg', 'unwrap'],
        'Math - Hyperbolic': ['sinh', 'cosh', 'tanh', 'arcsinh', 'arccosh', 'arctanh',
                             'asinh', 'acosh', 'atanh'],
        'Math - Rounding': ['around', 'round', 'rint', 'fix', 'floor', 'ceil', 'trunc'],
        'Math - Special': ['clip', 'convolve', 'correlate', 'diff', 'gradient', 'interp', 'trapezoid',
                          'sinc', 'i0', 'heaviside'],
        'Math - Floating': ['copysign', 'frexp', 'ldexp', 'nextafter', 'spacing', 'modf', 'signbit',
                           'isfinite', 'isinf', 'isnan', 'isnat', 'fmax', 'fmin', 'float_power'],
        'Logic': ['all', 'any', 'logical_and', 'logical_or', 'logical_xor', 'logical_not',
                 'allclose', 'isclose', 'array_equal', 'array_equiv', 'greater', 'greater_equal',
                 'less', 'less_equal', 'equal', 'not_equal', 'isreal', 'iscomplex', 'isrealobj',
                 'iscomplexobj', 'isscalar', 'isfortran'],
        'Bitwise': ['bitwise_and', 'bitwise_or', 'bitwise_xor', 'bitwise_not', 'invert',
                   'left_shift', 'right_shift', 'bitwise_left_shift', 'bitwise_right_shift',
                   'bitwise_invert', 'bitwise_count', 'packbits', 'unpackbits'],
        'Statistics': ['sum', 'prod', 'mean', 'std', 'var', 'min', 'max', 'amax', 'amin',
                      'argmin', 'argmax', 'nanmin', 'nanmax', 'nansum', 'nanprod', 'nanmean',
                      'nanstd', 'nanvar', 'nanargmin', 'nanargmax', 'nanmedian', 'nanpercentile',
                      'nanquantile', 'median', 'average', 'percentile', 'quantile', 'ptp',
                      'corrcoef', 'cov', 'histogram', 'histogram2d', 'histogramdd',
                      'histogram_bin_edges', 'bincount', 'digitize', 'count_nonzero',
                      'cumsum', 'cumprod', 'nancumsum', 'nancumprod', 'cumulative_sum', 'cumulative_prod'],
        'Sorting': ['sort', 'argsort', 'partition', 'argpartition', 'lexsort', 'sort_complex'],
        'Linear Algebra': ['dot', 'vdot', 'inner', 'outer', 'matmul', 'tensordot', 'einsum',
                          'einsum_path', 'kron', 'cross', 'trace', 'linalg'],
        'Complex': ['real', 'imag', 'conj', 'conjugate', 'angle', 'real_if_close'],
        'Polynomials': ['poly', 'poly1d', 'polyval', 'polyfit', 'polyadd', 'polysub', 'polymul',
                       'polydiv', 'polyder', 'polyint', 'roots'],
        'Set Operations': ['unique', 'unique_values', 'unique_counts', 'unique_inverse', 'unique_all',
                          'union1d', 'intersect1d', 'setdiff1d', 'setxor1d', 'isin', 'in1d', 'ediff1d'],
        'I/O': ['save', 'load', 'savetxt', 'loadtxt', 'genfromtxt', 'savez', 'savez_compressed',
               'array2string', 'array_repr', 'array_str', 'format_float_positional',
               'format_float_scientific', 'set_printoptions', 'get_printoptions', 'printoptions',
               'binary_repr', 'base_repr'],
        'Window Functions': ['bartlett', 'blackman', 'hamming', 'hanning', 'kaiser'],
        'Type Handling': ['dtype', 'finfo', 'iinfo', 'can_cast', 'result_type', 'min_scalar_type',
                         'promote_types', 'common_type', 'isdtype', 'issubdtype', 'mintypecode',
                         'typename', 'typecodes', 'sctypeDict', 'astype'],
        'Constants': ['e', 'pi', 'euler_gamma', 'inf', 'nan', 'newaxis'],
        'Broadcasting': ['broadcast', 'broadcast_to', 'broadcast_arrays', 'broadcast_shapes'],
        'Memory': ['may_share_memory', 'shares_memory', 'ndim', 'shape', 'size'],
        'Error Handling': ['errstate', 'geterr', 'seterr', 'geterrcall', 'seterrcall',
                          'getbufsize', 'setbufsize'],
        'Iterators': ['nditer', 'ndenumerate', 'ndindex', 'flatiter', 'nested_iters'],
        'Utility': ['iterable', 'vectorize', 'frompyfunc', 'piecewise', 'apply_along_axis',
                   'apply_over_axes', 'nan_to_num', 'isneginf', 'isposinf', 'info',
                   'show_config', 'show_runtime', 'get_include', 'little_endian'],
    }

    # Process each category
    for cat_name, funcs in categories.items():
        cat_implemented = 0
        cat_total = 0
        cat_items = []

        for func in funcs:
            if func in main:
                cat_total += 1
                total_numpy += 1
                impl = is_implemented(func)
                if impl:
                    cat_implemented += 1
                    total_implemented += 1
                    nw_name = get_numwasm_name(func)
                    if nw_name and nw_name.lower() != func.lower():
                        cat_items.append(f"- [x] `{func}` (NumWasm: `{nw_name}`)")
                    else:
                        cat_items.append(f"- [x] `{func}`")
                else:
                    cat_items.append(f"- [ ] `{func}`")

        if cat_items:
            md.append(f"### {cat_name} ({cat_implemented}/{cat_total})")
            md.append("")
            md.extend(cat_items)
            md.append("")

    # Handle remaining items not in categories
    categorized = set()
    for funcs in categories.values():
        categorized.update(funcs)

    remaining = [f for f in main if f not in categorized]
    if remaining:
        md.append("### Other")
        md.append("")
        for func in sorted(remaining):
            total_numpy += 1
            impl = is_implemented(func)
            if impl:
                total_implemented += 1
                nw_name = get_numwasm_name(func)
                if nw_name and nw_name.lower() != func.lower():
                    md.append(f"- [x] `{func}` (NumWasm: `{nw_name}`)")
                else:
                    md.append(f"- [x] `{func}`")
            else:
                md.append(f"- [ ] `{func}`")
        md.append("")

    # ============================================
    # LINALG
    # ============================================
    md.append("---")
    md.append("")
    md.append("## numpy.linalg")
    md.append("")

    linalg_impl = 0
    for func in sorted(numpy_exports['linalg']):
        total_numpy += 1
        impl = is_implemented(func)
        if impl:
            linalg_impl += 1
            total_implemented += 1
            md.append(f"- [x] `{func}`")
        else:
            md.append(f"- [ ] `{func}`")
    md.append("")
    md.append(f"**Status: {linalg_impl}/{len(numpy_exports['linalg'])} implemented**")
    md.append("")

    # ============================================
    # FFT
    # ============================================
    md.append("---")
    md.append("")
    md.append("## numpy.fft")
    md.append("")

    fft_impl = 0
    for func in sorted(numpy_exports['fft']):
        total_numpy += 1
        impl = is_implemented(func)
        if impl:
            fft_impl += 1
            total_implemented += 1
            md.append(f"- [x] `{func}`")
        else:
            md.append(f"- [ ] `{func}`")
    md.append("")
    md.append(f"**Status: {fft_impl}/{len(numpy_exports['fft'])} implemented**")
    md.append("")

    # ============================================
    # RANDOM
    # ============================================
    md.append("---")
    md.append("")
    md.append("## numpy.random")
    md.append("")

    random_impl = 0
    for func in sorted(numpy_exports['random']):
        total_numpy += 1
        impl = is_implemented(func)
        if impl:
            random_impl += 1
            total_implemented += 1
            md.append(f"- [x] `{func}`")
        else:
            md.append(f"- [ ] `{func}`")
    md.append("")
    md.append(f"**Status: {random_impl}/{len(numpy_exports['random'])} implemented**")
    md.append("")

    # ============================================
    # MA (MASKED ARRAYS)
    # ============================================
    md.append("---")
    md.append("")
    md.append("## numpy.ma (Masked Arrays)")
    md.append("")

    ma_impl = 0
    for func in sorted(numpy_exports['ma']):
        total_numpy += 1
        impl = is_implemented(func)
        if impl:
            ma_impl += 1
            total_implemented += 1
            md.append(f"- [x] `{func}`")
        else:
            md.append(f"- [ ] `{func}`")
    md.append("")
    md.append(f"**Status: {ma_impl}/{len(numpy_exports['ma'])} implemented**")
    md.append("")

    # ============================================
    # POLYNOMIAL
    # ============================================
    md.append("---")
    md.append("")
    md.append("## numpy.polynomial")
    md.append("")

    poly_impl = 0
    for func in sorted(numpy_exports['polynomial']):
        total_numpy += 1
        impl = is_implemented(func)
        if impl:
            poly_impl += 1
            total_implemented += 1
            md.append(f"- [x] `{func}`")
        else:
            md.append(f"- [ ] `{func}`")
    md.append("")
    md.append(f"**Status: {poly_impl}/{len(numpy_exports['polynomial'])} implemented**")
    md.append("")

    # ============================================
    # STRINGS
    # ============================================
    md.append("---")
    md.append("")
    md.append("## numpy.strings")
    md.append("")

    strings_impl = 0
    for func in sorted(numpy_exports['strings']):
        total_numpy += 1
        impl = is_implemented(func)
        if impl:
            strings_impl += 1
            total_implemented += 1
            nw_name = get_numwasm_name(func)
            if nw_name and nw_name.lower() != func.lower():
                md.append(f"- [x] `{func}` (NumWasm: `{nw_name}`)")
            else:
                md.append(f"- [x] `{func}`")
        else:
            md.append(f"- [ ] `{func}`")
    md.append("")
    md.append(f"**Status: {strings_impl}/{len(numpy_exports['strings'])} implemented**")
    md.append("")

    # ============================================
    # TESTING
    # ============================================
    md.append("---")
    md.append("")
    md.append("## numpy.testing")
    md.append("")

    testing_impl = 0
    for func in sorted(numpy_exports['testing']):
        total_numpy += 1
        impl = is_implemented(func)
        if impl:
            testing_impl += 1
            total_implemented += 1
            md.append(f"- [x] `{func}`")
        else:
            md.append(f"- [ ] `{func}`")
    md.append("")
    md.append(f"**Status: {testing_impl}/{len(numpy_exports['testing'])} implemented**")
    md.append("")

    # ============================================
    # SUMMARY
    # ============================================
    md.insert(7, "")
    md.insert(8, "## Overall Summary")
    md.insert(9, "")
    md.insert(10, f"**Total NumPy Functions: {total_numpy}**")
    md.insert(11, f"**Implemented in NumWasm: {total_implemented}**")
    md.insert(12, f"**Coverage: {total_implemented/total_numpy*100:.1f}%**")
    md.insert(13, "")

    # Write file
    output_path = "/workspace/packages/numwasm/numwasm-implementation-status.md"
    with open(output_path, 'w') as f:
        f.write('\n'.join(md))

    print(f"Generated: {output_path}")
    print(f"Total NumPy functions: {total_numpy}")
    print(f"Implemented in NumWasm: {total_implemented}")
    print(f"Coverage: {total_implemented/total_numpy*100:.1f}%")

if __name__ == "__main__":
    main()
