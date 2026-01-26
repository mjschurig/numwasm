/**
 * WASM Module Loader
 *
 * Handles loading and initialization of the NumJS WebAssembly module.
 * Works in Node.js environments (browser support planned).
 */

import type { WasmModule, WasmModuleFactory } from './types.js';

let wasmModule: WasmModule | null = null;
let modulePromise: Promise<WasmModule> | null = null;

/**
 * Detect if we're running in Node.js
 */
const isNode =
  typeof process !== 'undefined' &&
  process.versions != null &&
  process.versions.node != null;

/**
 * Get the directory path of this module in Node.js.
 * Works for both ESM and CJS contexts.
 */
async function getModuleDir(): Promise<string> {
  const nodePath = await import('path');
  const nodeUrl = await import('url');

  // Try ESM-style import.meta.url first
  try {
    if (import.meta && import.meta.url) {
      return nodePath.dirname(nodeUrl.fileURLToPath(import.meta.url));
    }
  } catch {
    // import.meta not available, fall through to CJS
  }

  // CJS fallback - __dirname should be available
  if (typeof __dirname !== 'undefined') {
    return __dirname;
  }

  // Last resort - use process.cwd()
  return process.cwd();
}

/**
 * Load the WASM glue code using require.
 * Works for both ESM (via createRequire) and CJS (native require).
 */
async function loadWasmGlue(wasmPath: string): Promise<WasmModuleFactory> {
  // Try to use createRequire for ESM context
  try {
    if (import.meta && import.meta.url) {
      const nodeModule = await import('module');
      const req = nodeModule.createRequire(import.meta.url);
      return req(wasmPath);
    }
  } catch {
    // Fall through to CJS require
  }

  // CJS fallback - use global require
  if (typeof require !== 'undefined') {
    return require(wasmPath);
  }

  throw new Error('Unable to load WASM module: no suitable require function available');
}

/**
 * Load the WASM module asynchronously.
 *
 * This function is idempotent - calling it multiple times will return
 * the same module instance. The module is loaded lazily on first call.
 *
 * @returns Promise that resolves to the loaded WASM module
 */
export async function loadWasmModule(): Promise<WasmModule> {
  // Return cached module if already loaded
  if (wasmModule) {
    return wasmModule;
  }

  // Return in-progress promise if currently loading
  if (modulePromise) {
    return modulePromise;
  }

  // Start loading the module
  modulePromise = (async () => {
    try {
      if (!isNode) {
        throw new Error(
          'Browser environment not yet supported. NumJS currently only works in Node.js.'
        );
      }

      const nodePath = await import('path');
      const moduleDir = await getModuleDir();
      const wasmPath = nodePath.join(moduleDir, 'wasm', 'numjs.cjs');

      const createModule = await loadWasmGlue(wasmPath);

      // Initialize the module (this loads the .wasm file)
      wasmModule = await createModule();

      return wasmModule;
    } catch (error) {
      // Reset promise on error to allow retry
      modulePromise = null;
      throw new Error(
        `Failed to load WASM module: ${error instanceof Error ? error.message : String(error)}`
      );
    }
  })();

  return modulePromise;
}

/**
 * Get the loaded WASM module synchronously.
 *
 * @throws Error if the module has not been loaded yet
 * @returns The loaded WASM module
 */
export function getWasmModule(): WasmModule {
  if (!wasmModule) {
    throw new Error(
      'WASM module not loaded. Call loadWasmModule() first and await the result.'
    );
  }
  return wasmModule;
}

/**
 * Check if the WASM module has been loaded.
 *
 * @returns true if the module is loaded and ready to use
 */
export function isWasmLoaded(): boolean {
  return wasmModule !== null;
}

/**
 * Reset the module state (mainly for testing purposes).
 * After calling this, loadWasmModule() will reload the module.
 */
export function resetWasmModule(): void {
  wasmModule = null;
  modulePromise = null;
}
