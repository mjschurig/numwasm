#!/usr/bin/env python3
"""
Generate a dependency tree of all C/C++/Fortran files in NumPy.
Analyzes #include directives to build the dependency graph.
"""
import os
import re
from pathlib import Path
from collections import defaultdict

NUMPY_ROOT = Path("/workspace/packages/numwasm/reference/numpy/numpy")

# File extensions to analyze
SOURCE_EXTS = {'.c', '.cpp', '.cxx', '.cc'}
HEADER_EXTS = {'.h', '.hpp', '.hxx'}
TEMPLATE_EXTS = {'.c.src', '.h.src', '.cpp.src'}
FORTRAN_EXTS = {'.f', '.f90', '.f77', '.for'}

def find_all_files():
    """Find all C/C++/Fortran/Header files."""
    files = {
        'sources': [],
        'headers': [],
        'templates': [],
        'fortran': [],
    }

    for root, dirs, filenames in os.walk(NUMPY_ROOT):
        # Skip test directories
        if 'tests' in root or 'test' in root:
            continue

        for fname in filenames:
            fpath = Path(root) / fname
            rel_path = fpath.relative_to(NUMPY_ROOT)

            if fname.endswith('.src'):
                files['templates'].append(rel_path)
            elif fpath.suffix in SOURCE_EXTS:
                files['sources'].append(rel_path)
            elif fpath.suffix in HEADER_EXTS:
                files['headers'].append(rel_path)
            elif fpath.suffix in FORTRAN_EXTS:
                files['fortran'].append(rel_path)

    return files

def parse_includes(filepath):
    """Parse #include directives from a file."""
    includes = {
        'system': [],  # <...>
        'local': [],   # "..."
    }

    try:
        with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()

        # Match #include <...> and #include "..."
        system_pattern = r'#include\s*<([^>]+)>'
        local_pattern = r'#include\s*"([^"]+)"'

        includes['system'] = re.findall(system_pattern, content)
        includes['local'] = re.findall(local_pattern, content)

    except Exception as e:
        pass

    return includes

def build_dependency_graph(files):
    """Build dependency graph from include directives."""
    graph = defaultdict(lambda: {'includes': [], 'included_by': []})

    all_headers = {str(h) for h in files['headers']}
    header_basenames = {h.name: str(h) for h in files['headers']}

    for src in files['sources'] + files['headers'] + files['templates']:
        full_path = NUMPY_ROOT / src
        includes = parse_includes(full_path)

        for inc in includes['local']:
            # Try to resolve the include
            resolved = None

            # Check if it's a direct match
            if inc in all_headers:
                resolved = inc
            elif inc in header_basenames:
                resolved = header_basenames[inc]
            else:
                # Try to find by basename
                inc_basename = Path(inc).name
                if inc_basename in header_basenames:
                    resolved = header_basenames[inc_basename]

            if resolved:
                graph[str(src)]['includes'].append(resolved)
                graph[resolved]['included_by'].append(str(src))

    return graph

def categorize_files(files):
    """Categorize files by module/component."""
    categories = {
        '_core/src/multiarray': {'name': 'Multiarray (Core)', 'sources': [], 'headers': []},
        '_core/src/umath': {'name': 'Umath (Universal Functions)', 'sources': [], 'headers': []},
        '_core/src/npysort': {'name': 'Sorting Algorithms', 'sources': [], 'headers': []},
        '_core/src/npymath': {'name': 'Math Support', 'sources': [], 'headers': []},
        '_core/src/common': {'name': 'Common Utilities', 'sources': [], 'headers': []},
        '_core/include': {'name': 'Core Headers', 'sources': [], 'headers': []},
        'linalg': {'name': 'Linear Algebra', 'sources': [], 'headers': []},
        'fft': {'name': 'FFT', 'sources': [], 'headers': []},
        'random': {'name': 'Random', 'sources': [], 'headers': []},
        'other': {'name': 'Other', 'sources': [], 'headers': []},
    }

    for src in files['sources']:
        src_str = str(src)
        categorized = False
        for cat_path in categories:
            if cat_path != 'other' and src_str.startswith(cat_path):
                categories[cat_path]['sources'].append(src)
                categorized = True
                break
        if not categorized:
            categories['other']['sources'].append(src)

    for hdr in files['headers']:
        hdr_str = str(hdr)
        categorized = False
        for cat_path in categories:
            if cat_path != 'other' and hdr_str.startswith(cat_path):
                categories[cat_path]['headers'].append(hdr)
                categorized = True
                break
        if not categorized:
            categories['other']['headers'].append(hdr)

    return categories

def generate_markdown(files, graph, categories):
    """Generate markdown documentation."""
    md = []
    md.append("# NumPy C/C++/Fortran Dependency Tree")
    md.append("")
    md.append("This document shows all native code files in NumPy and their dependencies.")
    md.append("")

    # Summary
    md.append("## Summary")
    md.append("")
    md.append(f"- **Source files (.c, .cpp)**: {len(files['sources'])}")
    md.append(f"- **Header files (.h, .hpp)**: {len(files['headers'])}")
    md.append(f"- **Template files (.src)**: {len(files['templates'])}")
    md.append(f"- **Fortran files (.f, .f90)**: {len(files['fortran'])}")
    md.append(f"- **Total**: {sum(len(v) for v in files.values())}")
    md.append("")

    # Category breakdown
    md.append("## Files by Module")
    md.append("")
    md.append("| Module | Sources | Headers |")
    md.append("|--------|---------|---------|")
    for cat_path, cat_info in sorted(categories.items()):
        if cat_info['sources'] or cat_info['headers']:
            md.append(f"| {cat_info['name']} | {len(cat_info['sources'])} | {len(cat_info['headers'])} |")
    md.append("")

    # Detailed file tree by category
    md.append("---")
    md.append("")
    md.append("## Detailed File Tree")
    md.append("")

    for cat_path, cat_info in sorted(categories.items()):
        if not cat_info['sources'] and not cat_info['headers']:
            continue

        md.append(f"### {cat_info['name']}")
        md.append("")

        if cat_info['headers']:
            md.append("**Headers:**")
            md.append("```")
            for h in sorted(cat_info['headers']):
                md.append(f"  {h}")
            md.append("```")
            md.append("")

        if cat_info['sources']:
            md.append("**Sources:**")
            md.append("```")
            for s in sorted(cat_info['sources']):
                md.append(f"  {s}")
            md.append("```")
            md.append("")

    # Templates
    if files['templates']:
        md.append("### Template Files (.src)")
        md.append("")
        md.append("These files are preprocessed to generate actual C/C++ code.")
        md.append("")
        md.append("```")
        for t in sorted(files['templates']):
            md.append(f"  {t}")
        md.append("```")
        md.append("")

    # Fortran
    if files['fortran']:
        md.append("### Fortran Files")
        md.append("")
        md.append("```")
        for f in sorted(files['fortran']):
            md.append(f"  {f}")
        md.append("```")
        md.append("")

    # Dependency graph
    md.append("---")
    md.append("")
    md.append("## Dependency Graph")
    md.append("")
    md.append("Shows which files include which headers.")
    md.append("")

    # Group by category for cleaner output
    for cat_path, cat_info in sorted(categories.items()):
        if not cat_info['sources']:
            continue

        md.append(f"### {cat_info['name']} Dependencies")
        md.append("")

        for src in sorted(cat_info['sources']):
            src_str = str(src)
            if src_str in graph and graph[src_str]['includes']:
                md.append(f"**`{src}`**")
                md.append("```")
                for inc in sorted(graph[src_str]['includes']):
                    md.append(f"  └── {inc}")
                md.append("```")
                md.append("")

    # Core headers and what includes them
    md.append("---")
    md.append("")
    md.append("## Core Header Usage")
    md.append("")
    md.append("Shows which core headers are most widely used.")
    md.append("")

    # Find most-included headers
    header_usage = []
    for hdr in files['headers']:
        hdr_str = str(hdr)
        if hdr_str in graph:
            count = len(graph[hdr_str]['included_by'])
            if count > 0:
                header_usage.append((hdr_str, count, graph[hdr_str]['included_by']))

    header_usage.sort(key=lambda x: -x[1])

    md.append("| Header | Used By (count) |")
    md.append("|--------|-----------------|")
    for hdr, count, _ in header_usage[:30]:  # Top 30
        md.append(f"| `{hdr}` | {count} |")
    md.append("")

    # Detailed usage for key headers
    key_headers = [
        '_core/include/numpy/arrayobject.h',
        '_core/include/numpy/ndarrayobject.h',
        '_core/include/numpy/ndarraytypes.h',
        '_core/include/numpy/ufuncobject.h',
        '_core/include/numpy/npy_math.h',
    ]

    md.append("### Key Header Details")
    md.append("")
    for kh in key_headers:
        if kh in graph and graph[kh]['included_by']:
            md.append(f"<details>")
            md.append(f"<summary><b>{kh}</b> (used by {len(graph[kh]['included_by'])} files)</summary>")
            md.append("")
            md.append("```")
            for user in sorted(graph[kh]['included_by'])[:20]:
                md.append(f"  {user}")
            if len(graph[kh]['included_by']) > 20:
                md.append(f"  ... and {len(graph[kh]['included_by']) - 20} more")
            md.append("```")
            md.append("</details>")
            md.append("")

    return '\n'.join(md)

def main():
    print("Analyzing NumPy C/C++/Fortran files...")

    files = find_all_files()
    print(f"Found {sum(len(v) for v in files.values())} files")

    graph = build_dependency_graph(files)
    print(f"Built dependency graph with {len(graph)} nodes")

    categories = categorize_files(files)

    md = generate_markdown(files, graph, categories)

    output_path = Path("/workspace/packages/numwasm/numpy-dependency-tree.md")
    with open(output_path, 'w') as f:
        f.write(md)

    print(f"Generated: {output_path}")

if __name__ == "__main__":
    main()
