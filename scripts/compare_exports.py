#!/usr/bin/env python3
"""
Compare NumPy and NumWasm exports to find missing functionality.
"""
import json

def normalize_name(name):
    """Normalize function names for comparison."""
    # Common TypeScript/JavaScript naming conventions vs Python
    mappings = {
        'applyAlongAxis': 'apply_along_axis',
        'applyOverAxes': 'apply_over_axes',
        'arrayRepr': 'array_repr',
        'arrayStr': 'array_str',
        'broadcastArrays': 'broadcast_arrays',
        'broadcastShapes': 'broadcast_shapes',
        'broadcastTo': 'broadcast_to',
        'countNonzero': 'count_nonzero',
        'deleteArr': 'delete',
        'var_': 'var',
        'commonType': 'common_type',
        'promoteTypes': 'promote_types',
        'ravelMultiIndex': 'ravel_multi_index',
        'unravelIndex': 'unravel_index',
        'uniqueAll': 'unique_all',
        'uniqueCounts': 'unique_counts',
        'uniqueIndex': 'unique_inverse',  # Similar
        'uniqueInverse': 'unique_inverse',
        'uniqueValues': 'unique_values',
        'recarray': 'recarray',
        'setPrintoptions': 'set_printoptions',
        'getPrintoptions': 'get_printoptions',
        'formatFloatPositional': 'format_float_positional',
        'formatFloatScientific': 'format_float_scientific',
        'binaryRepr': 'binary_repr',
        'baseRepr': 'base_repr',
    }
    return mappings.get(name, name)

def load_json(filepath):
    with open(filepath, 'r') as f:
        return json.load(f)

def main():
    numpy_exports = load_json("/workspace/numpy_exports.json")
    numwasm_exports = load_json("/workspace/numwasm_exports.json")

    # Flatten NumWasm exports
    numwasm_all = set()
    for category, items in numwasm_exports.items():
        for item in items:
            # Clean up namespace markers
            item = item.replace(' (namespace)', '')
            numwasm_all.add(item.lower())
            numwasm_all.add(normalize_name(item).lower())

    # Also add the original names
    for category, items in numwasm_exports.items():
        for item in items:
            item = item.replace(' (namespace)', '')
            numwasm_all.add(item)
            numwasm_all.add(normalize_name(item))

    # Categories of missing functionality
    missing = {
        'main_namespace': [],
        'linalg': [],
        'fft': [],
        'random': [],
        'ma': [],
        'polynomial': [],
        'strings': [],
        'testing': [],
    }

    # Compare main namespace
    for item in numpy_exports['main_namespace']:
        if item not in numwasm_all and item.lower() not in numwasm_all:
            missing['main_namespace'].append(item)

    # Compare linalg
    for item in numpy_exports['linalg']:
        if item not in numwasm_all and item.lower() not in numwasm_all:
            missing['linalg'].append(item)

    # Compare fft
    for item in numpy_exports['fft']:
        if item not in numwasm_all and item.lower() not in numwasm_all:
            missing['fft'].append(item)

    # Compare random
    for item in numpy_exports['random']:
        if item not in numwasm_all and item.lower() not in numwasm_all:
            missing['random'].append(item)

    # Compare ma
    for item in numpy_exports['ma']:
        if item not in numwasm_all and item.lower() not in numwasm_all:
            missing['ma'].append(item)

    # Compare polynomial
    for item in numpy_exports['polynomial']:
        if item not in numwasm_all and item.lower() not in numwasm_all:
            missing['polynomial'].append(item)

    # Compare strings
    for item in numpy_exports['strings']:
        if item not in numwasm_all and item.lower() not in numwasm_all:
            missing['strings'].append(item)

    # Compare testing
    for item in numpy_exports['testing']:
        if item not in numwasm_all and item.lower() not in numwasm_all:
            missing['testing'].append(item)

    # Save results
    with open("/workspace/missing_exports.json", "w") as f:
        json.dump(missing, f, indent=2)

    # Print summary
    print("Missing NumPy functionality in NumWasm:")
    print("=" * 50)
    total = 0
    for category, items in missing.items():
        count = len(items)
        total += count
        print(f"{category}: {count} missing")
    print(f"\nTotal missing: {total}")

if __name__ == "__main__":
    main()
