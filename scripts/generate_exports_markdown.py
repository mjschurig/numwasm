#!/usr/bin/env python3
"""
Generate markdown file with all NumPy and NumWasm exports.
"""
import json

def load_json(filepath):
    with open(filepath, 'r') as f:
        return json.load(f)

def main():
    numpy_exports = load_json("/workspace/numpy_exports.json")
    numwasm_exports = load_json("/workspace/numwasm_exports.json")

    md = []
    md.append("# NumPy and NumWasm Exports Reference")
    md.append("")
    md.append("This document contains comprehensive lists of all exports from NumPy (reference) and NumWasm.")
    md.append("")

    # Summary
    md.append("## Summary")
    md.append("")
    md.append("| Category | NumPy Exports | NumWasm Exports |")
    md.append("|----------|---------------|-----------------|")

    numpy_total = sum(len(v) for v in numpy_exports.values())
    numwasm_total = sum(len(v) for v in numwasm_exports.values())

    md.append(f"| **Total** | **{numpy_total}** | **{numwasm_total}** |")
    md.append("")

    # ============================================
    # NUMPY SECTION
    # ============================================
    md.append("---")
    md.append("")
    md.append("# Part 1: NumPy Reference Exports")
    md.append("")

    # Main namespace
    md.append("## numpy (Main Namespace)")
    md.append("")
    md.append(f"**Total: {len(numpy_exports['main_namespace'])} exports**")
    md.append("")

    # Categorize main namespace exports
    main_exports = numpy_exports['main_namespace']

    # Types/Classes (capitalized)
    types_classes = sorted([e for e in main_exports if e[0].isupper() and not e.startswith('__')])
    functions = sorted([e for e in main_exports if not e[0].isupper() and not e.startswith('__')])

    md.append("### Types and Classes")
    md.append("")
    for item in types_classes:
        md.append(f"- `{item}`")
    md.append("")

    md.append("### Functions and Constants")
    md.append("")
    for item in functions:
        md.append(f"- `{item}`")
    md.append("")

    # numpy.linalg
    md.append("## numpy.linalg")
    md.append("")
    md.append(f"**Total: {len(numpy_exports['linalg'])} exports**")
    md.append("")
    for item in sorted(numpy_exports['linalg']):
        md.append(f"- `{item}`")
    md.append("")

    # numpy.fft
    md.append("## numpy.fft")
    md.append("")
    md.append(f"**Total: {len(numpy_exports['fft'])} exports**")
    md.append("")
    for item in sorted(numpy_exports['fft']):
        md.append(f"- `{item}`")
    md.append("")

    # numpy.random
    md.append("## numpy.random")
    md.append("")
    md.append(f"**Total: {len(numpy_exports['random'])} exports**")
    md.append("")
    for item in sorted(numpy_exports['random']):
        md.append(f"- `{item}`")
    md.append("")

    # numpy.ma
    md.append("## numpy.ma (Masked Arrays)")
    md.append("")
    md.append(f"**Total: {len(numpy_exports['ma'])} exports**")
    md.append("")
    for item in sorted(numpy_exports['ma']):
        md.append(f"- `{item}`")
    md.append("")

    # numpy.polynomial
    md.append("## numpy.polynomial")
    md.append("")
    md.append(f"**Total: {len(numpy_exports['polynomial'])} exports**")
    md.append("")
    for item in sorted(numpy_exports['polynomial']):
        md.append(f"- `{item}`")
    md.append("")

    # numpy.strings
    md.append("## numpy.strings")
    md.append("")
    md.append(f"**Total: {len(numpy_exports['strings'])} exports**")
    md.append("")
    for item in sorted(numpy_exports['strings']):
        md.append(f"- `{item}`")
    md.append("")

    # numpy.testing
    md.append("## numpy.testing")
    md.append("")
    md.append(f"**Total: {len(numpy_exports['testing'])} exports**")
    md.append("")
    for item in sorted(numpy_exports['testing']):
        md.append(f"- `{item}`")
    md.append("")

    # numpy.lib
    md.append("## numpy.lib (Internal Library Functions)")
    md.append("")
    md.append(f"**Total: {len(numpy_exports['lib'])} exports**")
    md.append("")
    for item in sorted(numpy_exports['lib']):
        md.append(f"- `{item}`")
    md.append("")

    # ============================================
    # NUMWASM SECTION
    # ============================================
    md.append("---")
    md.append("")
    md.append("# Part 2: NumWasm Exports")
    md.append("")

    # Classes
    md.append("## Classes")
    md.append("")
    md.append(f"**Total: {len(numwasm_exports['classes'])} classes**")
    md.append("")
    for item in numwasm_exports['classes']:
        md.append(f"- `{item}`")
    md.append("")

    # Constants
    md.append("## Constants")
    md.append("")
    md.append(f"**Total: {len(numwasm_exports['constants'])} constants**")
    md.append("")
    for item in numwasm_exports['constants']:
        md.append(f"- `{item}`")
    md.append("")

    # Types
    md.append("## Types")
    md.append("")
    md.append(f"**Total: {len(numwasm_exports['types'])} types**")
    md.append("")
    for item in numwasm_exports['types']:
        md.append(f"- `{item}`")
    md.append("")

    # Functions
    md.append("## Functions")
    md.append("")
    md.append(f"**Total: {len(numwasm_exports['functions'])} functions**")
    md.append("")
    for item in numwasm_exports['functions']:
        md.append(f"- `{item}`")
    md.append("")

    # Write to file
    output_path = "/workspace/docs-site/public/exports-reference.md"
    with open(output_path, 'w') as f:
        f.write('\n'.join(md))

    print(f"Generated: {output_path}")
    print(f"NumPy total exports: {numpy_total}")
    print(f"NumWasm total exports: {numwasm_total}")

if __name__ == "__main__":
    main()
