/**
 * WASM Module Loader
 *
 * Handles loading and initialization of the NumJS WebAssembly module.
 * Works in both Node.js and browser environments.
 */

import type { WasmModule, WasmModuleFactory } from './types.js';

let wasmModule: WasmModule | null = null;
let modulePromise: Promise<WasmModule> | null = null;

/**
 * Configuration options for WASM module loading
 */
export interface WasmLoadConfig {
  /**
   * URL or path to the numjs.wasm file.
   * - In Node.js: defaults to auto-detected path relative to module
   * - In browser: auto-detected via import.meta.url, or can be explicitly set
   */
  wasmUrl?: string;
}

let wasmConfig: WasmLoadConfig = {};

/**
 * Configure WASM loading before initialization.
 * Must be called before loadWasmModule() for custom URLs.
 *
 * @example
 * // Browser with custom WASM location
 * configureWasm({ wasmUrl: '/static/wasm/numjs.wasm' });
 * await loadWasmModule();
 *
 * @example
 * // CDN usage
 * configureWasm({ wasmUrl: 'https://cdn.example.com/numjs/numjs.wasm' });
 * await loadWasmModule();
 */
export function configureWasm(config: WasmLoadConfig): void {
  if (wasmModule || modulePromise) {
    throw new Error('Cannot configure WASM after module has started loading');
  }
  wasmConfig = { ...config };
}

/**
 * Detect if we're running in Node.js
 */
const isNode =
  typeof process !== 'undefined' &&
  process.versions != null &&
  process.versions.node != null;

/**
 * Detect if we're running in a browser or web worker
 */
const isBrowser =
  typeof window !== 'undefined' ||
  (typeof self !== 'undefined' &&
    typeof (self as unknown as { importScripts?: unknown }).importScripts ===
      'function');

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

  throw new Error(
    'Unable to load WASM module: no suitable require function available'
  );
}

/**
 * Load the WASM module in Node.js environment.
 */
async function loadWasmNode(): Promise<WasmModule> {
  const nodePath = await import('path');
  const moduleDir = await getModuleDir();
  const wasmPath = nodePath.join(moduleDir, 'wasm', 'numjs.cjs');

  const createModule = await loadWasmGlue(wasmPath);

  // Initialize the module (this loads the .wasm file)
  return await createModule();
}

/**
 * Get the base URL for loading WASM assets in browser.
 * Uses a dynamic approach to avoid bundler transformation of import.meta.url.
 */
function getBrowserBaseUrl(): string {
  // Try to get the URL of the current script
  // This works because the script is loaded as an ES module
  try {
    // Use indirect eval to prevent bundler from transforming import.meta
    // eslint-disable-next-line no-eval
    const meta = eval('import.meta');
    if (meta && meta.url) {
      return new URL('.', meta.url).href;
    }
  } catch {
    // import.meta not available
  }

  // Fallback: try to find the script URL from document
  if (typeof document !== 'undefined') {
    const scripts = document.querySelectorAll('script[type="module"]');
    for (let i = 0; i < scripts.length; i++) {
      const script = scripts[i];
      const src = script.getAttribute('src');
      if (src && src.includes('numjs')) {
        return new URL('.', new URL(src, location.href)).href;
      }
    }
  }

  // Last resort: use the current page location
  if (typeof location !== 'undefined') {
    // Assume WASM files are at /dist/ relative to origin
    return new URL('/dist/', location.origin).href;
  }

  throw new Error(
    'Could not determine base URL for WASM files. Please use configureWasm({ wasmUrl: "..." })'
  );
}

/**
 * Load the WASM module in browser environment.
 */
async function loadWasmBrowser(): Promise<WasmModule> {
  // Use configured URL or auto-detect
  let wasmUrl = wasmConfig.wasmUrl;
  let glueUrl: string;

  if (wasmUrl) {
    // If wasmUrl is provided, derive glue URL from it
    const wasmUrlObj = new URL(wasmUrl);
    glueUrl = new URL('numjs.mjs', wasmUrlObj.href.replace(/numjs\.wasm$/, '')).href;
  } else {
    // Auto-detect base URL
    const baseUrl = getBrowserBaseUrl();
    wasmUrl = new URL('wasm/numjs.wasm', baseUrl).href;
    glueUrl = new URL('wasm/numjs.mjs', baseUrl).href;
  }

  // Dynamically import the ESM glue code from the resolved URL
  // We use dynamic import with the full URL to load from dist at runtime
  const glueModule = await import(/* @vite-ignore */ glueUrl);
  const createModule: WasmModuleFactory =
    glueModule.default || glueModule.createNumJSModule;

  if (!createModule) {
    throw new Error('Failed to load WASM glue code: no factory function found');
  }

  // Create the module with custom locateFile to resolve the WASM URL
  const module = await createModule({
    locateFile: (path: string) => {
      if (path.endsWith('.wasm')) {
        return wasmUrl!;
      }
      return path;
    },
  });

  return module;
}

/**
 * Load the WASM module asynchronously.
 *
 * This function is idempotent - calling it multiple times will return
 * the same module instance. The module is loaded lazily on first call.
 *
 * @returns Promise that resolves to the loaded WASM module
 *
 * @example
 * // Basic usage (works in both Node.js and browser)
 * await loadWasmModule();
 * const arr = await zeros([3, 3]);
 *
 * @example
 * // Browser with custom WASM location
 * configureWasm({ wasmUrl: '/static/numjs.wasm' });
 * await loadWasmModule();
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
      if (isNode) {
        wasmModule = await loadWasmNode();
      } else if (isBrowser) {
        wasmModule = await loadWasmBrowser();
      } else {
        throw new Error(
          'Unknown environment: neither Node.js nor browser detected'
        );
      }

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
  wasmConfig = {};
}
