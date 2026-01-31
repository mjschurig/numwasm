/**
 * MFEM WASM Module Loader
 *
 * Handles loading and initialization of the MFEM WebAssembly module.
 * Works in both Node.js and browser environments.
 */

import type { MFEMModule, MFEMModuleFactory } from './types.js';

let wasmModule: MFEMModule | null = null;
let modulePromise: Promise<MFEMModule> | null = null;

export interface MFEMLoadConfig {
  wasmUrl?: string;
}

let wasmConfig: MFEMLoadConfig = {};

/**
 * Configure the MFEM WASM loader before initialization.
 * Must be called before loadMFEMModule().
 *
 * @param config Configuration options
 * @throws Error if module has already started loading
 */
export function configureMFEM(config: MFEMLoadConfig): void {
  if (wasmModule || modulePromise) {
    throw new Error(
      'Cannot configure MFEM WASM after module has started loading'
    );
  }
  wasmConfig = { ...config };
}

const isNode =
  typeof process !== 'undefined' &&
  process.versions != null &&
  process.versions.node != null;

const isBrowser =
  typeof window !== 'undefined' ||
  (typeof self !== 'undefined' &&
    typeof (self as unknown as { importScripts?: unknown }).importScripts ===
      'function');

async function getModuleDir(): Promise<string> {
  const nodePath = await import('path');
  const nodeUrl = await import('url');
  try {
    if (import.meta && import.meta.url) {
      return nodePath.dirname(nodeUrl.fileURLToPath(import.meta.url));
    }
  } catch {
    // fall through
  }
  if (typeof __dirname !== 'undefined') {
    return __dirname;
  }
  return process.cwd();
}

async function loadWasmGlue(wasmPath: string): Promise<MFEMModuleFactory> {
  try {
    if (import.meta && import.meta.url) {
      const nodeModule = await import('module');
      const req = nodeModule.createRequire(import.meta.url);
      return req(wasmPath);
    }
  } catch {
    // fall through
  }
  if (typeof require !== 'undefined') {
    return require(wasmPath);
  }
  throw new Error(
    'Unable to load MFEM WASM module: no suitable require function available'
  );
}

async function findWasmCjs(): Promise<string> {
  const nodePath = await import('path');
  const nodeFs = await import('fs');
  const moduleDir = await getModuleDir();

  const candidates = [nodePath.join(moduleDir, 'wasm', 'mfem.cjs')];

  // Walk up to find package root (directory containing package.json)
  let dir = moduleDir;
  for (let i = 0; i < 10; i++) {
    const pkgPath = nodePath.join(dir, 'package.json');
    if (nodeFs.existsSync(pkgPath)) {
      candidates.push(nodePath.join(dir, 'dist', 'wasm', 'mfem.cjs'));
      break;
    }
    const parent = nodePath.dirname(dir);
    if (parent === dir) break;
    dir = parent;
  }

  for (const candidate of candidates) {
    if (nodeFs.existsSync(candidate)) {
      return candidate;
    }
  }

  throw new Error(
    `Cannot find mfem.cjs. Searched:\n${candidates.map((c) => '  - ' + c).join('\n')}\nRun 'pnpm run build:wasm' to compile the MFEM WASM module.`
  );
}

async function loadWasmNode(): Promise<MFEMModule> {
  const wasmPath = await findWasmCjs();
  const createModule = await loadWasmGlue(wasmPath);
  return await createModule();
}

function getBrowserAssetUrls(): { wasmUrl: string; glueUrl: string } {
  const wasmUrl = new URL('./wasm/mfem.wasm', import.meta.url).href;
  const glueUrl = new URL('./wasm/mfem.mjs', import.meta.url).href;
  return { wasmUrl, glueUrl };
}

async function loadWasmBrowser(): Promise<MFEMModule> {
  let wasmUrl: string;
  let glueUrl: string;

  if (wasmConfig.wasmUrl) {
    wasmUrl = wasmConfig.wasmUrl;
    const wasmUrlObj = new URL(wasmUrl);
    glueUrl = new URL(
      'mfem.mjs',
      wasmUrlObj.href.replace(/mfem\.wasm$/, '')
    ).href;
  } else {
    const urls = getBrowserAssetUrls();
    wasmUrl = urls.wasmUrl;
    glueUrl = urls.glueUrl;
  }

  const glueModule = await import(/* @vite-ignore */ glueUrl);
  const createModule: MFEMModuleFactory =
    glueModule.default || glueModule.createMFEMModule;

  if (!createModule) {
    throw new Error(
      'Failed to load MFEM WASM glue code: no factory function found'
    );
  }

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
 * Load and initialize the MFEM WASM module.
 *
 * This function is idempotent - calling it multiple times returns the same
 * module instance. The module is loaded lazily on first call.
 *
 * @returns Promise resolving to the initialized MFEM module
 * @throws Error if loading fails
 *
 * @example
 * ```typescript
 * const mfem = await loadMFEMModule();
 * const meshPtr = mfem._mfem_mesh_create();
 * ```
 */
export async function loadMFEMModule(): Promise<MFEMModule> {
  if (wasmModule) return wasmModule;
  if (modulePromise) return modulePromise;

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
      modulePromise = null;
      throw new Error(
        `Failed to load MFEM WASM module: ${error instanceof Error ? error.message : String(error)}`
      );
    }
  })();

  return modulePromise;
}

/**
 * Get the loaded MFEM module synchronously.
 *
 * @returns The MFEM module if loaded
 * @throws Error if module is not yet loaded
 *
 * @example
 * ```typescript
 * await loadMFEMModule();
 * const mfem = getMFEMModule(); // Now safe to call
 * ```
 */
export function getMFEMModule(): MFEMModule {
  if (!wasmModule) {
    throw new Error(
      'MFEM WASM module not loaded. Call loadMFEMModule() first and await the result.'
    );
  }
  return wasmModule;
}

/**
 * Check if the MFEM module has been loaded.
 *
 * @returns true if the module is loaded and ready
 */
export function isMFEMLoaded(): boolean {
  return wasmModule !== null;
}

/**
 * Reset the module state for testing purposes.
 * After calling this, loadMFEMModule() will load a fresh instance.
 */
export function resetMFEMModule(): void {
  wasmModule = null;
  modulePromise = null;
  wasmConfig = {};
}
