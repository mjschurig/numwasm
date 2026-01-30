export type PackageId = 'numwasm' | 'sciwasm' | 'symwasm';

export interface WasmExportsData {
  packages: {
    [packageId: string]: {
      [moduleName: string]: string[];
    };
  };
  typescript: {
    [packageId: string]: Array<{
      name: string;
      kind: 'function' | 'class' | 'variable' | 'type';
    }>;
  };
  tsToWasm: {
    [packageId: string]: {
      [tsFunctionName: string]: string[];
    };
  };
}

// Node data types
export interface CyNodeData {
  id: string;
  type: 'package' | 'wasmModule' | 'function';
  label: string;
  packageId: string;
  color: string;
  // For wasmModule nodes
  functionCount?: number;
  // For function nodes
  moduleName?: string;
  kind?: 'wasm' | 'typescript';
  tsKind?: 'function' | 'class' | 'variable' | 'type';
  isShared?: boolean;
  hasWasmBinding?: boolean;
}

// Edge data types
export interface CyEdgeData {
  id: string;
  source: string;
  target: string;
  edgeType: 'containment' | 'binding';
  color: string;
  width: number;
  opacity?: number;
}

// Element types (simplified, no longer depends on Cytoscape)
export interface CyNode {
  data: CyNodeData;
}

export interface CyEdge {
  data: CyEdgeData;
}

export type CyElement = CyNode | CyEdge;

// Package colors
export const PACKAGE_COLORS: Record<PackageId, string> = {
  numwasm: '#2dd4a8', // Teal (primary)
  sciwasm: '#a78bfa', // Purple
  symwasm: '#fb923c', // Orange
};

// Lighter versions for WASM modules
export const PACKAGE_COLORS_LIGHT: Record<PackageId, string> = {
  numwasm: '#5eead4',
  sciwasm: '#c4b5fd',
  symwasm: '#fdba74',
};

// Graph filter state
export interface GraphFilters {
  packages: Set<PackageId>;
  showWasmFunctions: boolean;
  showTsExports: boolean;
  searchQuery: string;
}

export const DEFAULT_FILTERS: GraphFilters = {
  packages: new Set(['numwasm', 'sciwasm', 'symwasm']),
  showWasmFunctions: true,
  showTsExports: true,
  searchQuery: '',
};
