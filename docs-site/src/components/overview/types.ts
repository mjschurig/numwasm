export type PackageId =
  | 'numwasm'
  | 'sciwasm'
  | 'symwasm'
  | 'arwasm'
  | 'lawasm'
  | 'linwasm'
  | 'quadwasm'
  | 'superluwasm'
  | 'xsfwasm'
  | 'odewasm';

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

// Function documentation extracted from TypeDoc
export interface FunctionDocs {
  [packageId: string]: {
    [functionName: string]: string;
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
  description?: string;
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
  numwasm: '#2dd4a8',     // Teal (primary)
  sciwasm: '#a78bfa',     // Purple
  symwasm: '#fb923c',     // Orange
  arwasm: '#f472b6',      // Pink
  lawasm: '#60a5fa',      // Blue
  linwasm: '#34d399',     // Green
  quadwasm: '#fbbf24',    // Yellow
  superluwasm: '#a3e635', // Lime
  xsfwasm: '#e879f9',     // Fuchsia
  odewasm: '#f87171',     // Red
};

// Lighter versions for WASM modules
export const PACKAGE_COLORS_LIGHT: Record<PackageId, string> = {
  numwasm: '#5eead4',
  sciwasm: '#c4b5fd',
  symwasm: '#fdba74',
  arwasm: '#f9a8d4',
  lawasm: '#93c5fd',
  linwasm: '#6ee7b7',
  quadwasm: '#fcd34d',
  superluwasm: '#bef264',
  xsfwasm: '#f0abfc',
  odewasm: '#fca5a5',
};

// Graph filter state
export interface GraphFilters {
  packages: Set<PackageId>;
  showWasmFunctions: boolean;
  showTsExports: boolean;
  searchQuery: string;
}

// All available packages
export const ALL_PACKAGES: PackageId[] = [
  'numwasm',
  'sciwasm',
  'symwasm',
  'arwasm',
  'lawasm',
  'linwasm',
  'quadwasm',
  'superluwasm',
  'xsfwasm',
  'odewasm',
];

export const DEFAULT_FILTERS: GraphFilters = {
  packages: new Set(ALL_PACKAGES),
  showWasmFunctions: true,
  showTsExports: true,
  searchQuery: '',
};
