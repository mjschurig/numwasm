/**
 * WASM Module Loader
 * Handles loading the SymEngine WASM module in both Node.js and browser environments
 */

import type { SymwasmModule, SymwasmModuleFactory } from './wasm-types.js';

// Singleton instance
let wasmModule: SymwasmModule | null = null;

// Detect environment
const isNode =
  typeof process !== 'undefined' &&
  process.versions != null &&
  process.versions.node != null;
const isBrowser = typeof window !== 'undefined';

/**
 * Load the WASM module (singleton pattern)
 * @returns Promise that resolves to the loaded WASM module
 */
export async function loadWasmModule(): Promise<SymwasmModule> {
  if (wasmModule) {
    return wasmModule;
  }

  if (isNode) {
    return loadWasmModuleNode();
  } else if (isBrowser) {
    return loadWasmModuleBrowser();
  } else {
    throw new Error('Unsupported environment for WASM module loading');
  }
}

/**
 * Load WASM module in Node.js environment
 */
async function loadWasmModuleNode(): Promise<SymwasmModule> {
  // Dynamic import to avoid bundling issues
  const path = await import('path');
  const { fileURLToPath } = await import('url');
  const { existsSync } = await import('fs');

  // Find the dist/wasm directory
  // Walk up from current directory to find the package root
  let currentDir: string;
  if (typeof __dirname !== 'undefined') {
    // CommonJS
    currentDir = __dirname;
  } else {
    // ES modules
    currentDir = path.dirname(fileURLToPath(import.meta.url));
  }

  // Walk up to find dist/wasm/symwasm.cjs
  let wasmPath: string | null = null;
  for (let i = 0; i < 5; i++) {
    const candidate = path.join(currentDir, 'wasm', 'symwasm.cjs');
    if (existsSync(candidate)) {
      wasmPath = candidate;
      break;
    }
    currentDir = path.dirname(currentDir);
  }

  if (!wasmPath) {
    throw new Error(
      'Could not find symwasm.cjs. Make sure to run `pnpm run build:wasm` first.'
    );
  }

  // Import the CommonJS module
  const createModule = (await import(wasmPath)) as unknown as {
    default: SymwasmModuleFactory;
  };

  const factory = createModule.default;
  wasmModule = await factory();

  return wasmModule;
}

/**
 * Load WASM module in browser environment
 */
async function loadWasmModuleBrowser(): Promise<SymwasmModule> {
  // In browser, we need to know the path to the WASM files
  // This assumes the WASM files are served alongside the bundle
  const baseUrl = import.meta.url
    ? new URL('.', import.meta.url).href
    : window.location.href;

  // Construct path to symwasm.mjs
  const wasmUrl = new URL('../wasm/symwasm.mjs', baseUrl).href;

  // Dynamic import
  const createModule = (await import(
    /* @vite-ignore */ wasmUrl
  )) as unknown as { default: SymwasmModuleFactory };

  const factory = createModule.default;
  wasmModule = await factory();

  return wasmModule;
}

/**
 * Get the loaded WASM module
 * @throws Error if module is not loaded yet
 * @returns The loaded WASM module
 */
export function getWasmModule(): SymwasmModule {
  if (!wasmModule) {
    throw new Error(
      'WASM module not loaded. Call loadWasmModule() first and await the result.'
    );
  }
  return wasmModule;
}

/**
 * Check if the WASM module is loaded
 * @returns true if module is loaded
 */
export function isWasmModuleLoaded(): boolean {
  return wasmModule !== null;
}

/**
 * Unload the WASM module (for testing/cleanup)
 * WARNING: This will invalidate all existing SymEngine objects
 */
export function unloadWasmModule(): void {
  wasmModule = null;
}
