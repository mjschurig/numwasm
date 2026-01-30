import { useState, useEffect, useMemo } from 'react';
import type {
  WasmExportsData,
  FunctionDocs,
  GraphFilters,
  CyElement,
  CyNode,
  CyEdge,
} from './types';
import { PACKAGE_COLORS, PACKAGE_COLORS_LIGHT, ALL_PACKAGES } from './types';

interface GraphData {
  elements: CyElement[];
  isLoading: boolean;
  error: string | null;
}

/**
 * Find functions that are shared across multiple packages
 */
function findSharedFunctions(data: WasmExportsData): Set<string> {
  const functionPackages: Record<string, Set<string>> = {};

  for (const [packageId, modules] of Object.entries(data.packages)) {
    for (const [, functions] of Object.entries(modules)) {
      for (const fn of functions) {
        const normalized = fn.replace(/^_+/, '').replace(/_+$/, '');
        if (!functionPackages[normalized]) {
          functionPackages[normalized] = new Set();
        }
        functionPackages[normalized].add(packageId);
      }
    }
  }

  for (const [packageId, exports] of Object.entries(data.typescript)) {
    for (const exp of exports) {
      const normalized = exp.name.toLowerCase();
      if (!functionPackages[normalized]) {
        functionPackages[normalized] = new Set();
      }
      functionPackages[normalized].add(packageId);
    }
  }

  const shared = new Set<string>();
  for (const [name, packages] of Object.entries(functionPackages)) {
    if (packages.size > 1) {
      shared.add(name);
    }
  }
  return shared;
}

/**
 * Build Cytoscape elements from WASM exports data
 */
function buildElements(
  data: WasmExportsData,
  filters: GraphFilters,
  docs: FunctionDocs
): CyElement[] {
  const elements: CyElement[] = [];
  const sharedFunctions = findSharedFunctions(data);
  const nodeIds = new Set<string>();

  // Dynamically iterate over all packages present in data
  const availablePackages = ALL_PACKAGES.filter(
    (pkg) => data.packages[pkg] || data.typescript[pkg]
  );

  for (const packageId of availablePackages) {
    if (!filters.packages.has(packageId)) continue;

    const modules = data.packages[packageId] || {};
    const tsExports = data.typescript[packageId] || [];
    const wasmModuleCount = Object.keys(modules).length;
    const tsExportCount = tsExports.length;

    // Package node
    const packageNodeId = `pkg-${packageId}`;
    const packageNode: CyNode = {
      data: {
        id: packageNodeId,
        type: 'package',
        label: `${packageId}\n${wasmModuleCount} modules · ${tsExportCount} exports`,
        packageId,
        color: PACKAGE_COLORS[packageId],
      },
    };
    elements.push(packageNode);
    nodeIds.add(packageNodeId);

    if (filters.showWasmFunctions) {
      for (const [moduleName, functions] of Object.entries(modules)) {
        const moduleNodeId = `mod-${packageId}-${moduleName}`;

        // WASM Module node
        const moduleNode: CyNode = {
          data: {
            id: moduleNodeId,
            type: 'wasmModule',
            label: `${moduleName}.wasm`,
            packageId,
            functionCount: functions.length,
            color: PACKAGE_COLORS_LIGHT[packageId],
          },
        };
        elements.push(moduleNode);
        nodeIds.add(moduleNodeId);

        // Package → Module edge
        const pkgModEdge: CyEdge = {
          data: {
            id: `edge-${packageNodeId}-${moduleNodeId}`,
            source: packageNodeId,
            target: moduleNodeId,
            edgeType: 'containment',
            color: PACKAGE_COLORS[packageId],
            width: 2,
          },
        };
        elements.push(pkgModEdge);

        // Function nodes (no hard filtering - visual indication only)
        for (const fn of functions) {
          const normalizedName = fn.replace(/^_+/, '').replace(/_+$/, '');
          const isShared = sharedFunctions.has(normalizedName);

          // Look up description for WASM function (stored with underscore prefix)
          const description = docs[packageId]?.[fn];

          const fnNodeId = `wasm-${packageId}-${moduleName}-${fn}`;
          const fnNode: CyNode = {
            data: {
              id: fnNodeId,
              type: 'function',
              label: fn,
              packageId,
              moduleName,
              kind: 'wasm',
              isShared,
              hasWasmBinding: false,
              color: PACKAGE_COLORS[packageId],
              description,
            },
          };
          elements.push(fnNode);
          nodeIds.add(fnNodeId);

          // Module → Function edge
          const modFnEdge: CyEdge = {
            data: {
              id: `edge-${moduleNodeId}-${fnNodeId}`,
              source: moduleNodeId,
              target: fnNodeId,
              edgeType: 'containment',
              color: '#10b981',
              width: 1,
            },
          };
          elements.push(modFnEdge);
        }
      }
    }

    if (filters.showTsExports) {
      // Get TS→WASM mappings for this package
      const tsToWasmMappings = data.tsToWasm?.[packageId] || {};

      // No hard filtering - visual indication only for search
      for (const exp of tsExports) {
        const normalizedName = exp.name.toLowerCase();
        const isShared = sharedFunctions.has(normalizedName);

        // Look for WASM bindings
        let wasmBindings: string[] = tsToWasmMappings[exp.name] || [];

        // If it's a class, also collect all ClassName.* method bindings
        if (exp.kind === 'class') {
          const methodPrefix = `${exp.name}.`;
          for (const [key, bindings] of Object.entries(tsToWasmMappings)) {
            if (key.startsWith(methodPrefix)) {
              wasmBindings = [...new Set([...wasmBindings, ...(bindings as string[])])];
            }
          }
        }

        const tsNodeId = `ts-${packageId}-${exp.name}`;
        const description = docs[packageId]?.[exp.name];
        const tsNode: CyNode = {
          data: {
            id: tsNodeId,
            type: 'function',
            label: exp.name,
            packageId,
            kind: 'typescript',
            tsKind: exp.kind,
            isShared,
            hasWasmBinding: wasmBindings.length > 0,
            color: PACKAGE_COLORS[packageId],
            description,
          },
        };
        elements.push(tsNode);
        nodeIds.add(tsNodeId);

        // Package → TS Function edge
        const pkgTsEdge: CyEdge = {
          data: {
            id: `edge-${packageNodeId}-${tsNodeId}`,
            source: packageNodeId,
            target: tsNodeId,
            edgeType: 'containment',
            color: PACKAGE_COLORS[packageId],
            width: 1,
            opacity: 0.6,
          },
        };
        elements.push(pkgTsEdge);

        // TS → WASM binding edges
        if (filters.showWasmFunctions && wasmBindings.length > 0) {
          for (const wasmFn of wasmBindings) {
            // Find the WASM function node
            for (const [moduleName, functions] of Object.entries(modules)) {
              if (functions.includes(wasmFn)) {
                const wasmNodeId = `wasm-${packageId}-${moduleName}-${wasmFn}`;
                // Only add edge if the target node exists
                if (nodeIds.has(wasmNodeId)) {
                  const bindingEdge: CyEdge = {
                    data: {
                      id: `edge-ts-wasm-${tsNodeId}-${wasmNodeId}`,
                      source: tsNodeId,
                      target: wasmNodeId,
                      edgeType: 'binding',
                      color: '#22c55e',
                      width: 2,
                    },
                  };
                  elements.push(bindingEdge);
                }
                break;
              }
            }
          }
        }
      }
    }
  }

  return elements;
}

export function useGraphData(filters: GraphFilters): GraphData {
  const [wasmExports, setWasmExports] = useState<WasmExportsData | null>(null);
  const [functionDocs, setFunctionDocs] = useState<FunctionDocs>({});
  const [isLoading, setIsLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    Promise.all([
      fetch('/wasm-exports.json').then((res) => {
        if (!res.ok) {
          throw new Error(`Failed to load wasm-exports.json: ${res.status}`);
        }
        return res.json();
      }),
      fetch('/function-docs.json')
        .then((res) => (res.ok ? res.json() : {}))
        .catch(() => ({})), // Gracefully handle missing docs
    ])
      .then(([exports, docs]: [WasmExportsData, FunctionDocs]) => {
        setWasmExports(exports);
        setFunctionDocs(docs);
        setIsLoading(false);
      })
      .catch((err) => {
        setError(err.message);
        setIsLoading(false);
      });
  }, []);

  const elements = useMemo(() => {
    if (!wasmExports) {
      return [];
    }
    return buildElements(wasmExports, filters, functionDocs);
  }, [wasmExports, filters, functionDocs]);

  return { elements, isLoading, error };
}
