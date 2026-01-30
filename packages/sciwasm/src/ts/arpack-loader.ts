/**
 * ARPACK WASM Module Loader
 *
 * Handles loading and initialization of the ARPACK WebAssembly module.
 * ARPACK is a library for solving large scale eigenvalue problems.
 * Works in both Node.js and browser environments.
 */

import type { ARPACKModule, ARPACKModuleFactory } from './arpack-types.js';

let wasmModule: ARPACKModule | null = null;
let modulePromise: Promise<ARPACKModule> | null = null;

export interface ARPACKLoadConfig {
  wasmUrl?: string;
}

let wasmConfig: ARPACKLoadConfig = {};

export function configureARPACK(config: ARPACKLoadConfig): void {
  if (wasmModule || modulePromise) {
    throw new Error('Cannot configure ARPACK WASM after module has started loading');
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

async function loadWasmGlue(wasmPath: string): Promise<ARPACKModuleFactory> {
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
  throw new Error('Unable to load ARPACK WASM module: no suitable require function available');
}

async function findWasmCjs(): Promise<string> {
  const nodePath = await import('path');
  const nodeFs = await import('fs');
  const moduleDir = await getModuleDir();

  const candidates = [
    nodePath.join(moduleDir, 'wasm', 'arpack.cjs'),
  ];

  // Walk up to find package root (directory containing package.json)
  let dir = moduleDir;
  for (let i = 0; i < 10; i++) {
    const pkgPath = nodePath.join(dir, 'package.json');
    if (nodeFs.existsSync(pkgPath)) {
      candidates.push(nodePath.join(dir, 'dist', 'wasm', 'arpack.cjs'));
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
    `Cannot find arpack.cjs. Searched:\n${candidates.map(c => '  - ' + c).join('\n')}\nRun 'bash scripts/build-arpack.sh' to compile the ARPACK WASM module.`
  );
}

async function loadWasmNode(): Promise<ARPACKModule> {
  const wasmPath = await findWasmCjs();
  const createModule = await loadWasmGlue(wasmPath);
  return await createModule();
}

function getBrowserAssetUrls(): { wasmUrl: string; glueUrl: string } {
  const wasmUrl = new URL('./wasm/arpack.wasm', import.meta.url).href;
  const glueUrl = new URL('./wasm/arpack.mjs', import.meta.url).href;
  return { wasmUrl, glueUrl };
}

async function loadWasmBrowser(): Promise<ARPACKModule> {
  let wasmUrl: string;
  let glueUrl: string;

  if (wasmConfig.wasmUrl) {
    wasmUrl = wasmConfig.wasmUrl;
    const wasmUrlObj = new URL(wasmUrl);
    glueUrl = new URL('arpack.mjs', wasmUrlObj.href.replace(/arpack\.wasm$/, '')).href;
  } else {
    const urls = getBrowserAssetUrls();
    wasmUrl = urls.wasmUrl;
    glueUrl = urls.glueUrl;
  }

  const glueModule = await import(/* @vite-ignore */ glueUrl);
  const createModule: ARPACKModuleFactory =
    glueModule.default || glueModule.createARPACKModule;

  if (!createModule) {
    throw new Error('Failed to load ARPACK WASM glue code: no factory function found');
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

export async function loadARPACKModule(): Promise<ARPACKModule> {
  if (wasmModule) return wasmModule;
  if (modulePromise) return modulePromise;

  modulePromise = (async () => {
    try {
      if (isNode) {
        wasmModule = await loadWasmNode();
      } else if (isBrowser) {
        wasmModule = await loadWasmBrowser();
      } else {
        throw new Error('Unknown environment: neither Node.js nor browser detected');
      }
      return wasmModule;
    } catch (error) {
      modulePromise = null;
      throw new Error(
        `Failed to load ARPACK WASM module: ${error instanceof Error ? error.message : String(error)}`
      );
    }
  })();

  return modulePromise;
}

export function getARPACKModule(): ARPACKModule {
  if (!wasmModule) {
    throw new Error(
      'ARPACK WASM module not loaded. Call loadARPACKModule() first and await the result.'
    );
  }
  return wasmModule;
}

export function isARPACKLoaded(): boolean {
  return wasmModule !== null;
}

export function resetARPACKModule(): void {
  wasmModule = null;
  modulePromise = null;
  wasmConfig = {};
}
