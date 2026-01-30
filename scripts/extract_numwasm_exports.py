#!/usr/bin/env python3
"""
Extract all exported functions from NumWasm TypeScript source.
"""
import os
import re
import json
from pathlib import Path

NUMWASM_SRC = Path("/workspace/packages/numwasm/src/ts")

def extract_exports_from_index(filepath):
    """Extract export statements from TypeScript index file."""
    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()

    # Remove comments first
    # Remove single-line comments
    content = re.sub(r'//[^\n]*', '', content)

    exports = {
        'functions': [],
        'classes': [],
        'types': [],
        'constants': [],
    }

    # Extract named exports like: export { foo, bar, baz } from './module.js';
    pattern = r'export\s*\{([^}]+)\}\s*from'
    for match in re.finditer(pattern, content, re.DOTALL):
        items = match.group(1)
        # Parse individual exports, handling "as" renames
        for item in items.split(','):
            item = item.strip()
            if not item or item.startswith('type '):
                continue
            # Handle "original as alias" pattern
            if ' as ' in item:
                parts = item.split(' as ')
                name = parts[1].strip()
            else:
                name = item.strip()

            # Clean up any trailing comments or whitespace
            name = name.split()[0] if name else ''

            if name and not name.startswith('_'):
                # Categorize
                if name[0].isupper():
                    if name.endswith('Error') or name.endswith('Exception') or name.endswith('Warning'):
                        exports['classes'].append(name)
                    elif any(kw in name for kw in ['Type', 'Options', 'Result', 'Config', 'Mode', 'Kind']):
                        exports['types'].append(name)
                    else:
                        exports['classes'].append(name)
                elif name.isupper() or name in ['NAN', 'PINF', 'NINF', 'PZERO', 'NZERO', 'e', 'pi', 'nan', 'inf']:
                    exports['constants'].append(name)
                else:
                    exports['functions'].append(name)

    # Also extract re-exported modules like: export * as testing from './testing/index.js';
    pattern2 = r'export\s+\*\s+as\s+(\w+)\s+from'
    for match in re.finditer(pattern2, content):
        name = match.group(1)
        if not name.startswith('_'):
            exports['functions'].append(f"{name} (namespace)")

    return exports

def main():
    print("Extracting NumWasm exports...\n")

    index_file = NUMWASM_SRC / "index.ts"
    exports = extract_exports_from_index(index_file)

    # Sort everything and remove duplicates
    for key in exports:
        exports[key] = sorted(set(exports[key]))

    # Save results
    with open("/workspace/numwasm_exports.json", "w") as f:
        json.dump(exports, f, indent=2)

    # Print summary
    total = 0
    for category, items in exports.items():
        count = len(items)
        total += count
        print(f"{category}: {count} exports")

    print(f"\nTotal: {total} exports")
    print("\nResults saved to /workspace/numwasm_exports.json")

if __name__ == "__main__":
    main()
